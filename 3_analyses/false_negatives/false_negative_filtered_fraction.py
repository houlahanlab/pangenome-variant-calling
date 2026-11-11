#!/usr/bin/env python3

# Text statistics only (Results, false negatives): makes no figure or table, and only prints the percentages quoted
# in the manuscript (Strelka2 91.1%, 1,893 / 2,079; Mutect2 54.8%, 401 / 732 of VG pangenome FNs called but filtered).
#
# Of the VG pangenome FN sites, the variants that the caller called but removed with its own filter, versus never
# called, versus passed the caller's filter but excluded by the depth / target filter.
#   Called but filtered: raw record with a non-PASS FILTER (LowEVS for Strelka2; weak_evidence, strand_bias, ... for Mutect2)
#   Never proposed:      no raw record (NO_RECORD)
#   Raw record PASS:     the caller passed it, so it was lost to the depth / target filter
#
# Input:  <data_dir>/False_Negatives_Analysis/Strelka_EVS/strelka_evs_filter_results.tsv (evs_measure_strelka.py)
#         <data_dir>/False_Negatives_Analysis/Mutect2_TLOD/mutect2_tlod_filter_results.tsv (tlod_measure_mutect2.py)
# Output: printed only
#
# Usage:  python3 3_analyses/false_negatives/false_negative_filtered_fraction.py

import argparse
import os

import yaml

FILES = {
    "Strelka2": os.path.join("Strelka_EVS", "strelka_evs_filter_results.tsv"),
    "Mutect2": os.path.join("Mutect2_TLOD", "mutect2_tlod_filter_results.tsv"),
}


def summarize(path):
    filter_counts = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if cols[idx["alignment"]] == "VG_T2T_pangenome" and cols[idx["status"]] == "FN":
                filt = cols[idx["filter"]]
                filter_counts[filt] = filter_counts.get(filt, 0) + 1

    total = sum(filter_counts.values())
    no_record = filter_counts.get("NO_RECORD", 0)
    liftover_failed = filter_counts.get("LIFTOVER_FAILED", 0)
    raw_pass = filter_counts.get("PASS", 0)
    return {"total_fn": total, "called_but_filtered": total - no_record - liftover_failed - raw_pass,
            "no_record": no_record, "liftover_failed": liftover_failed, "raw_pass": raw_pass,
            "filter_counts": filter_counts}


def main():
    parser = argparse.ArgumentParser(description="Why the VG pangenome FN sites were missed.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    fn_dir = os.path.join(config["data_dir"], "False_Negatives_Analysis")
    summaries = {caller: summarize(os.path.join(fn_dir, name)) for caller, name in FILES.items()}

    for caller, s in summaries.items():
        pct = 100 * s["called_but_filtered"] / s["total_fn"] if s["total_fn"] else 0.0
        print(f"\n{'=' * 60}\n{caller}: FN sites missed by VG_T2T_pangenome\n{'=' * 60}")
        print(f"Total FN sites:                {s['total_fn']}")
        print(f"  Called, but filtered out:    {s['called_but_filtered']}  ({pct:.1f}%)")
        print(f"  Never proposed as candidate: {s['no_record']}")
        print(f"  Liftover failed:             {s['liftover_failed']}")
        print(f"  Raw record PASS (excluded for depth/region, not confidence): {s['raw_pass']}")
        print("\n  FILTER tag breakdown:")
        for tag, n in sorted(s["filter_counts"].items(), key=lambda kv: -kv[1]):
            print(f"    {tag:35s} {n}")

    print(f"\n{'=' * 60}")
    print("Sentence-ready (all categories, % of total FN):")
    for caller, s in summaries.items():
        t = s["total_fn"]

        def pct(n):
            return 100 * n / t if t else 0.0

        print(f"\n  {caller} (n={t} FN sites missed by VG_T2T_pangenome):")
        print(f"    Called, but filtered for low confidence: {pct(s['called_but_filtered']):.1f}% ({s['called_but_filtered']}/{t})")
        print(f"    Never proposed as a candidate at all:    {pct(s['no_record']):.1f}% ({s['no_record']}/{t})")
        print(f"    Excluded by depth/target-region criteria despite passing the caller's own filter: {pct(s['raw_pass']):.1f}% ({s['raw_pass']}/{t})")
        print(f"    Liftover mapping failed (untestable):    {pct(s['liftover_failed']):.1f}% ({s['liftover_failed']}/{t})")


if __name__ == "__main__":
    main()