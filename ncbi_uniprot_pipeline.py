#!/usr/bin/env python3
"""
Script Version 9.0: Reconciled Old-Core Framework with Advanced Speed Upgrades
Restores the stable, infallible query-search core of the original working script,
incorporating persistent JSON caching, controlled batch POST chunking, and strict case subdivision.
"""

import sys
import time
import json
import sqlite3
import argparse
import requests
from pathlib import Path
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# ============================================================================
# Global Configuration
# ============================================================================
RATE_LIMIT_DELAY = 0.4
CLUSTERED_NR_DB = Path(__file__).parent / "CLUSTERED_NR_SQLITE" / "cluster_data.sqlite3"
CACHE_DIR = Path(__file__).parent / ".cache"

ipg_cache = {}
uniprot_cache = {}

# ============================================================================
# Cache & Connection Management
# ============================================================================
def ensure_cache_dir():
    CACHE_DIR.mkdir(exist_ok=True)

def load_caches():
    global ipg_cache, uniprot_cache
    ensure_cache_dir()
    try:
        if (CACHE_DIR / 'ipg_cache.json').exists():
            with open(CACHE_DIR / 'ipg_cache.json', 'r') as f:
                ipg_cache.update(json.load(f))
        if (CACHE_DIR / 'uniprot_cache.json').exists():
            with open(CACHE_DIR / 'uniprot_cache.json', 'r') as f:
                uniprot_cache.update(json.load(f))
        print(f"[CACHE] Loaded {len(ipg_cache)} IPG entries and {len(uniprot_cache)} UniProt records.", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Failed to load caches: {e}", file=sys.stderr)

def save_caches():
    ensure_cache_dir()
    try:
        with open(CACHE_DIR / 'ipg_cache.json', 'w') as f:
            json.dump(ipg_cache, f)
        with open(CACHE_DIR / 'uniprot_cache.json', 'w') as f:
            json.dump(uniprot_cache, f)
    except Exception as e:
        print(f"[WARNING] Failed to save caches: {e}", file=sys.stderr)

def get_request_retry(url, params=None, data=None, method="GET"):
    time.sleep(RATE_LIMIT_DELAY)
    s = requests.Session()
    retries = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
    s.mount('https://', HTTPAdapter(max_retries=retries))
    s.mount('http://', HTTPAdapter(max_retries=retries))
    try:
        if method == "POST":
            return s.post(url, data=data, timeout=45)
        return s.get(url, params=params, timeout=30)
    except Exception:
        return None

# ============================================================================
# Stable Infallible UniProt Resolutions (The Original API Engine)
# ============================================================================
def search_direct_uniprot(ncbi_code):
    """Universal cross-reference tracker using the stable search endpoint."""
    if ncbi_code in uniprot_cache:
        return uniprot_cache[ncbi_code]

    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": f"xref:ncbi_protein:{ncbi_code}", "format": "json", "size": 1}
    response = get_request_retry(url, params=params)
    if response and response.status_code == 200:
        data = response.json()
        if data.get("results"):
            up_acc = data["results"][0].get("primaryAccession")
            if up_acc:
                uniprot_cache[ncbi_code] = up_acc
                return up_acc
    return None

def search_uniprot_uniparc(query_code):
    url = "https://rest.uniprot.org/uniparc/search"
    params = {"query": query_code, "format": "json", "size": 1}
    response = get_request_retry(url, params=params)
    if response and response.status_code == 200:
        data = response.json()
        if data.get("results"):
            return data["results"][0].get("uniParcId")
    return None

def check_uniparc_cross_refs(uniparc_id):
    url = f"https://rest.uniprot.org/uniparc/{uniparc_id}"
    response = get_request_retry(url)
    if response and response.status_code == 200:
        data = response.json()
        cross_refs = data.get("uniParcCrossReferences", [])
        for ref in cross_refs:
            db_name = ref.get("database", "")
            if db_name in ("UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL") and ref.get("active", True):
                return ref.get("id")
    return None

def search_uniref90(uniparc_id):
    url = "https://rest.uniprot.org/uniref/search"
    params = {"query": f"uniparc:{uniparc_id} AND identity:0.9", "format": "json", "size": 1}
    response = get_request_retry(url, params=params)
    if response and response.status_code == 200:
        data = response.json()
        if data.get("results"):
            cluster = data["results"][0]
            rep_member = cluster.get("representativeMember", {})
            accessions = rep_member.get("accessions", [])
            if accessions:
                return accessions[0]
    return None

def resolve_via_gateways(ncbi_code):
    """Hierarchical resolution path: Native Search -> UniParc -> UniRef90 -> UniProtKB Extraction."""
    # 1. Native direct attempt
    direct_up = search_direct_uniprot(ncbi_code)
    if direct_up: 
        return direct_up, "-", "-"

    # 2. UniParc gateway transition
    uniparc_id = search_uniprot_uniparc(ncbi_code)
    if uniparc_id:
        # Check if UniParc contains active links pointing to UniProtKB
        up_id = check_uniparc_cross_refs(uniparc_id)
        if up_id: 
            return up_id, uniparc_id, "-"
        
        # Final emergency fallback via UniRef90 cluster
        uniref_id = search_uniref90(uniparc_id)
        if uniref_id: 
            return uniref_id, uniparc_id, "UniRef90_Cluster"

    return None, "-", "-"

# ============================================================================
# Optimized Batch Methods (NCBI Eutils Guidelines Compliance)
# ============================================================================
def get_codes_from_ipg_bulk(list_of_ncbi_codes):
    """
    Extracts IPG twins in bulk using structured HTTP POST requests partitioned 
    into blocks of 200, protecting processing speed with the persistent memory cache.
    """
    clean_codes = sorted(list(set(c.strip() for c in list_of_ncbi_codes if c and c.strip() != "-")))
    if not clean_codes: 
        return {}

    results = {}
    to_fetch = []
    for c in clean_codes:
        if c in ipg_cache: 
            results[c] = set(ipg_cache[c])
        else:
            results[c] = {c}
            to_fetch.append(c)

    if not to_fetch: 
        return results

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    chunk_size = 200
    for i in range(0, len(to_fetch), chunk_size):
        chunk = to_fetch[i:i+chunk_size]
        payload = {"db": "ipg", "id": ",".join(chunk), "rettype": "ipg", "retmode": "text"}
        response = get_request_retry(url, data=payload, method="POST")
        if response and response.status_code == 200:
            lines = response.text.strip().split('\n')
            if len(lines) > 1:
                for line in lines[1:]:
                    columns = line.split('\t')
                    if len(columns) >= 7:
                        q_id = columns[1].strip()
                        twin = columns[6].strip()
                        if q_id in results and twin: 
                            results[q_id].add(twin)
                        for c_orig in chunk:
                            if c_orig in columns: 
                                results[c_orig].add(twin)

    for c in to_fetch: 
        ipg_cache[c] = list(results[c])
    return results

def extract_all_cluster_members_sqlite(list_of_ipgs):
    """Extracts all cluster members sharing a structural representative in bulk."""
    if not list_of_ipgs or not CLUSTERED_NR_DB.exists(): 
        return set()
    clean_ipgs = list(set(str(i).strip() for i in list_of_ipgs if i and str(i).strip() != "-"))
    
    all_members = set()
    try:
        conn = sqlite3.connect(str(CLUSTERED_NR_DB))
        cursor = conn.cursor()
        for i in range(0, len(clean_ipgs), 500):
            chunk = clean_ipgs[i:i+500]
            placeholders = ', '.join(['?'] * len(chunk))
            query = f"""
                SELECT member_accession FROM ClusterInfo 
                WHERE representative IN (
                    SELECT DISTINCT representative FROM ClusterInfo WHERE member_accession IN ({placeholders})
                )
            """
            cursor.execute(query, chunk)
            for row in cursor.fetchall(): 
                all_members.add(row[0])
        conn.close()
    except Exception:
        pass
    return all_members

# ============================================================================
# Main Sequential Cascading Framework (Version 9)
# ============================================================================
def main():
    load_caches()
    parser = argparse.ArgumentParser(description="High-Throughput NCBI to UniProtKB Translation Framework")
    parser.add_argument('-f', '--proteins_file', help="Path to the file containing raw protein accessions")
    parser.add_argument('-p', '--protein_ids', help="Comma-separated list of standalone NCBI accessions")
    args = parser.parse_args()
    
    ids = []
    base = "manual"
    if args.protein_ids: 
        ids = [i.strip() for i in args.protein_ids.split(",") if i.strip()]
    if args.proteins_file:
        raw_lines = open(args.proteins_file).read().splitlines()
        ids = [line.strip() for line in raw_lines if line.strip()]
        base = Path(args.proteins_file).stem

    if not ids: 
        return
    print(f"[INFO] Initializing translation pipeline version 9 for {len(ids)} identifiers.\n", file=sys.stderr)
    
    start_time = time.time()
    all_results = []
    stats = {"Case 0": 0, "Case 1.1": 0, "Case 1.2": 0, "Case 2": 0, "Fail": 0}

    for idx, ncbi_code in enumerate(ids, 1):
        res = {
            "ncbi_code": ncbi_code, "success": False, "case": "Fail",
            "final_ipg_code": "-", "final_uniprot_id": "-", "uniparc_bridge": "-", "uniref_bridge": "-",
            "all_ipg_found": set(), "all_cluster_members": set(), "ipgs_of_members_tried": set()
        }
        prefix = f"[{idx}/{len(ids)}] {ncbi_code:<15}"
        print(f"{prefix} Processing...", end="\r", file=sys.stderr, flush=True)

        # --------------------------------------------------------------------
        # CASE 0: Direct search over original code using stable gateways
        # --------------------------------------------------------------------
        up_id, uparc, uref = resolve_via_gateways(ncbi_code)
        if up_id:
            res.update({"success": True, "case": "Case 0", "final_ipg_code": ncbi_code, "final_uniprot_id": up_id, "uniparc_bridge": uparc, "uniref_bridge": uref})
            stats["Case 0"] += 1
            print(f"{prefix} -> [OK] CASE 0 (Direct Query) -> {up_id}", file=sys.stderr, flush=True)
            all_results.append(res)
            continue

        # --------------------------------------------------------------------
        # CASE 1: Query NCBI IPG and split results into lists for lookup
        # --------------------------------------------------------------------
        ipg_map = get_codes_from_ipg_bulk([ncbi_code])
        twins = ipg_map.get(ncbi_code, {ncbi_code})
        res["all_ipg_found"].update(twins)

        refseq_twins = [str(t) for t in twins if str(t).startswith(("NP_", "WP_", "XP_", "YP_"))]
        genbank_twins = [str(t) for t in twins if not str(t).startswith(("NP_", "WP_", "XP_", "YP_"))]

        # --- Sub-Case 1.1: RefSeq list mapping ---
        case1_found = False
        if refseq_twins:
            for r_twin in refseq_twins:
                up_id, uparc, uref = resolve_via_gateways(r_twin)
                if up_id:
                    res.update({"success": True, "case": "Case 1.1", "final_ipg_code": r_twin, "final_uniprot_id": up_id, "uniparc_bridge": uparc, "uniref_bridge": uref})
                    stats["Case 1.1"] += 1
                    print(f"{prefix} -> [OK] CASE 1.1 (RefSeq Twin: {r_twin}) -> {up_id}", file=sys.stderr, flush=True)
                    case1_found = True
                    break
        if case1_found: 
            all_results.append(res)
            continue

        # --- Sub-Case 1.2: GenBank/Other list mapping ---
        if genbank_twins:
            for g_twin in genbank_twins:
                up_id, uparc, uref = resolve_via_gateways(g_twin)
                if up_id:
                    res.update({"success": True, "case": "Case 1.2", "final_ipg_code": g_twin, "final_uniprot_id": up_id, "uniparc_bridge": uparc, "uniref_bridge": uref})
                    stats["Case 1.2"] += 1
                    print(f"{prefix} -> [OK] CASE 1.2 (GenBank Twin: {g_twin}) -> {up_id}", file=sys.stderr, flush=True)
                    case1_found = True
                    break
        if case1_found: 
            all_results.append(res)
            continue

        # --------------------------------------------------------------------
        # CASE 2: SQLite Clustered Extraction in bulk + chunked POST resolution
        # --------------------------------------------------------------------
        cluster_members = extract_all_cluster_members_sqlite(list(twins))
        if cluster_members:
            res["all_cluster_members"].update(cluster_members)
            
            # Extract IPG twin maps in bulk using mass 200 POST chunking
            members_ipg_map = get_codes_from_ipg_bulk(list(cluster_members))
            
            all_member_twins = set()
            for m_twins in members_ipg_map.values(): 
                all_member_twins.update(m_twins)
            
            if all_member_twins:
                res["ipgs_of_members_tried"].update(all_member_twins)
                
                case2_found = False
                for m_twin in sorted(list(all_member_twins)):
                    up_id, uparc, uref = resolve_via_gateways(m_twin)
                    if up_id:
                        res.update({"success": True, "case": "Case 2", "final_ipg_code": m_twin, "final_uniprot_id": up_id, "uniparc_bridge": uparc, "uniref_bridge": uref})
                        stats["Case 2"] += 1
                        print(f"{prefix} -> [OK] CASE 2 (Cluster Hit via: {m_twin}) -> {up_id}", file=sys.stderr, flush=True)
                        case2_found = True
                        break
                
                if case2_found: 
                    all_results.append(res)
                    continue

        # --------------------------------------------------------------------
        # FAIL: All hierarchical layers exhausted
        # --------------------------------------------------------------------
        stats["Fail"] += 1
        print(f"{prefix} -> [X] FAILED MAPPING", file=sys.stderr, flush=True)
        all_results.append(res)

        if idx % 100 == 0: 
            save_caches()

    print("\n" + "="*60, file=sys.stderr)
    print(" SUMMARY OF CASCADING RESOLUTION (VERSION 9)", file=sys.stderr)
    print("="*60, file=sys.stderr)
    for c_name, count in stats.items(): 
        print(f" * {c_name:<10}: {count} proteins", file=sys.stderr)
    print("="*60 + "\n", file=sys.stderr)

    for r in all_results:
        for k in ["all_ipg_found", "all_cluster_members", "ipgs_of_members_tried"]: 
            r[k] = sorted(list(r[k]))

    total_time = time.time() - start_time
    print_out(all_results, base, total_time)
    save_caches()

# ============================================================================
# Output Functions
# ============================================================================
def clean_visual(data):
    if not data or data == "-": 
        return "-"
    return "; ".join(str(x) for x in data) if isinstance(data, list) else str(data)

def print_out(results, base_filename, total_time):
    out_dir = Path("output") / "version9"
    trace_dir = out_dir / "trace"
    fail_dir = out_dir / "failures"
    for d in [out_dir, trace_dir, fail_dir]: 
        d.mkdir(exist_ok=True, parents=True)

    v = 1
    while (out_dir / f"{base_filename}_extended_mapping_v{v}.tsv").exists(): 
        v += 1
    f_main = out_dir / f"{base_filename}_extended_mapping_v{v}.tsv"
    f_trace = trace_dir / f"{base_filename}_v{v}_trace.tsv"
    f_json = trace_dir / f"{base_filename}_v{v}_trace.json"
    f_fails = fail_dir / f"{base_filename}_v{v}_failures.tsv"

    with open(f_main, 'w') as f:
        f.write("\t".join(["NCBI_Protein_ID", "Status", "Case", "UniProt_ID", "UniParc_Bridge", "UniRef90_Bridge"]) + "\n")
        for r in results:
            f.write("\t".join([r['ncbi_code'], "SUCCESS" if r['success'] else "FAILED", r['case'], r['final_uniprot_id'], r['uniparc_bridge'], r['uniref_bridge']]) + "\n")

    with open(f_trace, 'w') as f:
        f.write("\t".join(["NCBI_Protein_ID", "IPGs_Found", "Members_Extracted", "Deep_IPGs_Tried", "UniParc_Used", "UniRef_Used", "Success"]) + "\n")
        for r in results:
            f.write("\t".join([r['ncbi_code'], clean_visual(r['all_ipg_found']), clean_visual(r['all_cluster_members']), clean_visual(r['ipgs_of_members_tried']), r['uniparc_bridge'], r['uniref_bridge'], "YES" if r['success'] else "NO"]) + "\n")

    with open(f_json, 'w') as f: 
        json.dump(results, f, indent=4)
        
    with open(f_fails, 'w') as f:
        for r in results:
            if not r['success']: 
                f.write(r['ncbi_code'] + "\n")

    print(f"\n" + "="*50)
    print(f"COMPLETED IN: {total_time:.2f} seconds")
    print(f"Main output path (Version 9):    {f_main}")
    print(f"Trace details file:              {f_trace}")
    print(f"="*50)

if __name__ == '__main__':
    main()