import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
GNOMAD_QUERY = """query($id: String!) {
  variant(variantId: $id, dataset: gnomad_r4) { joint { ac an homozygote_count } }
}"""


def api(url, body=None):
    data = json.dumps(body).encode() if body else None
    req = urllib.request.Request(url, data, {"Content-Type": "application/json",
                                             "Accept": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=300))


def clinvar(chrom, pos, ref, alt):
    term = urllib.parse.quote(f"{chrom}[chr] AND {pos}[chrpos38]")
    ids = api(EUTILS + "esearch.fcgi?db=clinvar&retmode=json&term=" + term)["esearchresult"]["idlist"]
    time.sleep(0.4)
    if not ids:
        return None
    res = api(EUTILS + "esummary.fcgi?db=clinvar&retmode=json&id=" + ",".join(ids))["result"]
    time.sleep(0.4)
    for uid in res["uids"]:
        if res[uid]["variation_set"][0]["canonical_spdi"].endswith(f":{int(pos) - 1}:{ref}:{alt}"):
            return res[uid]
    return None


def gnomad(chrom, pos, ref, alt):
    body = {"query": GNOMAD_QUERY, "variables": {"id": f"{chrom}-{pos}-{ref}-{alt}"}}
    error = ""
    for _ in range(5):
        try:
            r = api("https://gnomad.broadinstitute.org/api", body)
            if r["data"]["variant"]:
                return r["data"]["variant"]
            error = r["errors"][0]["message"]
            if error == "Variant not found":
                break
        except urllib.error.HTTPError as e:
            error = str(e)
        time.sleep(30)
    print(chrom, pos, ref, alt, "gnomAD:", error)
    return {"error": error}


def load_or_fetch(name, fetch):
    path = f"{out}/{name}.json"
    if not os.path.exists(path):
        data = fetch()
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=1, ensure_ascii=False)
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def pick_transcript(rec):
    worst = rec["most_severe_consequence"]
    tr = [t for t in rec.get("transcript_consequences", []) if worst in t["consequence_terms"]]
    tr.sort(key=lambda t: ("mane_select" not in t, "canonical" not in t))
    return tr[0] if tr else {}


def vep_gnomad_af(rec, alt):
    for c in rec.get("colocated_variants", []):
        if c["id"].startswith("rs"):
            return c.get("frequencies", {}).get(alt, {})
    return {}


vcf, out = sys.argv[1], sys.argv[2]
variants = [line.split("\t")[:5] for line in open(vcf) if line.strip() and not line.startswith("#")]

ts = sum(ref + alt in ("AG", "GA", "CT", "TC") for _, _, _, ref, alt in variants)
print(f"записей: {len(variants)}, chr21: {sum(v[0] == '21' for v in variants)}, "
      f"chr22: {sum(v[0] == '22' for v in variants)}, транзиций: {ts}, трансверсий: {len(variants) - ts}")

vep = load_or_fetch("vep", lambda: {r["id"]: r for r in api(
    "https://rest.ensembl.org/vep/homo_sapiens/region",
    {"variants": [" ".join(v) + " . . ." for v in variants],
     "canonical": 1, "mane": 1, "hgvs": 1, "SIFT": "b", "PolyPhen": "b", "AlphaMissense": 1,
     "af_gnomade": 1, "af_gnomadg": 1})})
print("VEP ready")
clinvar_raw = load_or_fetch("clinvar", lambda: {v[2]: clinvar(v[0], v[1], v[3], v[4]) for v in variants})
print("ClinVar reay")
gnomad_raw = load_or_fetch("gnomad", lambda: {v[2]: gnomad(v[0], v[1], v[3], v[4]) for v in variants})
print("gnomAD ready")

rows = []
for chrom, pos, rs, ref, alt in variants:
    t = pick_transcript(vep[rs])
    mane = next((x for x in vep[rs].get("transcript_consequences", [])
                 if "mane_select" in x and x["gene_id"] == t.get("gene_id")), {})
    freq = vep_gnomad_af(vep[rs], alt)
    cv = clinvar_raw[rs]
    cls = cv["germline_classification"] if cv else {}
    g = gnomad_raw[rs].get("joint")
    rows.append({
        "variant": f"{chrom}:{pos} {ref}>{alt}",
        "rsID": rs,
        "gene": t.get("gene_symbol", ""),
        "most_severe_consequence": vep[rs]["most_severe_consequence"],
        "transcript": t.get("transcript_id", ""),
        "is_MANE": "yes" if "mane_select" in t else "no",
        "HGVSc": t.get("hgvsc", ""),
        "HGVSp": t.get("hgvsp", ""),
        "MANE_consequence": ",".join(mane.get("consequence_terms", [])),
        "SIFT": t.get("sift_prediction", ""),
        "PolyPhen": t.get("polyphen_prediction", ""),
        "AlphaMissense": (t.get("alphamissense") or {}).get("am_class", ""),
        "gnomAD_AF": f"{g['ac'] / g['an']:.3g}" if g and g["an"] else "",
        "gnomAD_hom": g["homozygote_count"] if g else "",
        "gnomADe_AF_vep": freq.get("gnomade", ""),
        "gnomADg_AF_vep": freq.get("gnomadg", ""),
        "ClinVar": "yes" if cv else "no",
        "ClinVar_ID": cv["accession"] if cv else "",
        "ClinVar_significance": cls.get("description") or ("no classification" if cv else ""),
        "ClinVar_review": cls.get("review_status", ""),
    })

with open(f"{out}/annotated_variants.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
print("готово:", len(rows), "строк")
