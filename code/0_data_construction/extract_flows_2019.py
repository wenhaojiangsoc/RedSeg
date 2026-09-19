# Extract POI x visitor-home-CBG monthly flows from SafeGraph 2019 patterns backfill.
# Usage: python3 extract_flows_2019.py <month_dir> <out_csv_gz>
# Output rows: safegraph_place_id,cbg,visitors   (one row per POI x home CBG for that month)
# US CBGs only (12-digit FIPS keys; Canadian "CA:..." keys dropped).
import csv, gzip, json, sys, glob, os
csv.field_size_limit(sys.maxsize)
month_dir, out_path = sys.argv[1], sys.argv[2]
parts = sorted(glob.glob(os.path.join(month_dir, "core_poi-patterns-part*.csv.gz")))
n_poi, n_rows = 0, 0
with gzip.open(out_path, "wt", compresslevel=6) as out:
    out.write("safegraph_place_id,cbg,visitors\n")
    for p in parts:
        with gzip.open(p, "rt", newline="") as f:
            for row in csv.DictReader(f):
                vh = row.get("visitor_home_cbgs")
                if not vh:
                    continue
                try:
                    d = json.loads(vh)
                except json.JSONDecodeError:
                    continue
                pid = row["safegraph_place_id"]
                wrote = False
                for k, v in d.items():
                    if len(k) == 12 and k.isdigit():
                        out.write(f"{pid},{k},{v}\n")
                        n_rows += 1
                        wrote = True
                if wrote:
                    n_poi += 1
print(f"{month_dir}: {n_poi} POIs, {n_rows} flow rows -> {out_path}", flush=True)
