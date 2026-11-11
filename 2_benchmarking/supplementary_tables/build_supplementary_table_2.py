#!/usr/bin/env python3

# Supplementary Table 2: runtime, CPU hours and memory of each alignment step and variant caller, from a rerun of
# pangenome-pipeline on 6 BLCA tumour-normal pairs (2 per ancestry) with the same resources for every alignment.
# BWA GRCh38 BAMs are the TCGA BAMs, so they have no alignment rows.
#
# Input:  <data_dir>/Memory_Test/ (job scripts and logs), sacct accounting saved in
#         <data_dir>/Supplementary_Tables/Supplementary_Table_2_sacct_snapshots/
# Output: <data_dir>/Supplementary_Tables/Supplementary_Table_2/Supplementary_Table_2.xlsx (median across jobs per tool and
#         alignment); prints the number of jobs behind each row
#
# Usage:  python3 2_benchmarking/supplementary_tables/build_supplementary_table_2.py

import argparse
import csv
import glob
import re
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path

import yaml

TABLE_ALIGNMENTS = {"GRCh38_bwa": "BWA_GRCh38", "GRCh38_linear": "VG_GRCh38", "T2T_pangenome": "VG_Pangenome"}

# Slurm job name prefix -> (stage, what the job runs)
STEPS = {
    "run_alignment": ("Alignment", "vg giraffe + vg surject"),
    "run_postprocessing": ("Alignment", "GATK SortSam, MarkDuplicates, BQSR"),
    "run_strelka": ("Variant calling", "Strelka2 (with Manta)"),
    "run_mutect2": ("Variant calling", "GATK Mutect2 + FilterMutectCalls"),
    "run_somaticsniper": ("Variant calling", "SomaticSniper"),
    "run_deepsomatic": ("Variant calling", "DeepSomatic"),
}

MEMWATCH_STEPS = {"run_strelka", "run_deepsomatic"}


def to_hours(text):
    """Slurm [D-]HH:MM:SS -> hours"""
    days, _, clock = text.rpartition("-")
    h, m, s = (["0", "0"] + clock.split(":"))[-3:]
    return int(days or 0) * 24 + int(h) + int(m) / 60 + float(s) / 3600


def median(values):
    values = [v for v in values if v is not None]
    return round(statistics.median(values), 2) if values else ""


def main():
    parser = argparse.ArgumentParser(description="Build Supplementary Table 2.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    data_dir = Path(config["data_dir"])
    memory_test_dir = data_dir / "Memory_Test"
    out_dir = data_dir / "Supplementary_Tables" / "Supplementary_Table_2"
    snapshot_dir = data_dir / "Supplementary_Tables" / "Supplementary_Table_2_sacct_snapshots"
    sacct_snapshot = snapshot_dir / "sacct_snapshot_memory_test.tsv"

    # Jobs whose commands all finished: job id -> (step, alignment, peak memory GB)
    used = {}
    for script in glob.glob(str(memory_test_dir / "*" / "*" / "logs" / "run_*_TCGA-*" / "script.sh")):
        alignment, folder = Path(script).parts[-5], Path(script).parts[-4]
        step = re.match(r"(run_[a-z0-9_]+?)_TCGA", Path(script).parts[-2]).group(1)
        logs = sorted(glob.glob(str(Path(script).parent / "slurm-*.out")))
        if step not in STEPS or not logs or (step in MEMWATCH_STEPS) != folder.endswith("_memwatch"):
            continue
        log = Path(logs[-1]).read_text(errors="replace")
        text = Path(script).read_text()
        if folder.endswith("_memwatch"):
            expected = text.count("memwatch.sh")
            exits = re.findall(r"Memwatch exit status: (\d+)", log)
            peaks = [int(k) for k in re.findall(r"Peak job memory, all programs \(kbytes\): (\d+)", log)]
        else:
            expected = text.count("/usr/bin/time -v") + text.count("somatic-sniper_timed/build")
            exits = re.findall(r"Exit status: (\d+)", log)
            peaks = [int(k) for k in re.findall(r"Maximum resident set size \(kbytes\): (\d+)", log)]
        if expected and len(exits) == expected and set(exits) == {"0"}:
            used[Path(logs[-1]).stem.split("-", 1)[1]] = (step, alignment, max(peaks) / 1024 ** 2)

    # Runtime, CPUs, CPU hours and allocated memory from sacct, saved once
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    if not sacct_snapshot.exists():
        sacct_snapshot.write_text(subprocess.run(
            ["sacct", "-X", "-P", "--delimiter=\t", "-j", ",".join(sorted(used)),
             "--format=JobID,JobName%60,Elapsed,AllocCPUS,CPUTime,ReqMem,State,NodeList"],
            capture_output=True, text=True, check=True).stdout)
    accounting = {}
    for row in csv.DictReader(sacct_snapshot.read_text().splitlines(), delimiter="\t"):
        accounting[row["JobID"]] = row

    out_dir.mkdir(parents=True, exist_ok=True)
    groups = defaultdict(list)
    for job_id, (step, alignment, memory) in used.items():
        row = accounting.get(job_id)
        if row is None or row["State"] != "COMPLETED" or alignment not in TABLE_ALIGNMENTS:
            continue
        # A few jobs have no start time in sacct: memory kept, runtime and CPU hours left out
        started = row["NodeList"] != "None assigned"
        groups[(step, alignment)].append({"cpus": int(row["AllocCPUS"]),
                                          "runtime": to_hours(row["Elapsed"]) if started else None,
                                          "cpu_time": to_hours(row["CPUTime"]) if started else None, "memory": memory,
                                          "allocated": int(re.match(r"(\d+)G", row["ReqMem"]).group(1))})

    # Table: one row per tool and alignment
    columns = ["Stage", "Tool", "Alignment", "CPUs", "Runtime_Hours", "CPU_Hours", "Memory_GB", "Allocated_Memory_GB"]
    rows = []
    for (step, alignment), group in sorted(groups.items(), key=lambda g: (list(STEPS).index(g[0][0]),
                                                                          list(TABLE_ALIGNMENTS).index(g[0][1]))):
        stage, tools = STEPS[step]
        rows.append([stage, tools, TABLE_ALIGNMENTS[alignment], min(j["cpus"] for j in group),
                     median([j["runtime"] for j in group]), median([j["cpu_time"] for j in group]),
                     median([j["memory"] for j in group]), median([j["allocated"] for j in group])])
        print(f"{tools}\t{TABLE_ALIGNMENTS[alignment]}\tmemory test jobs: {len(group)}")

    tsv = out_dir / ".Supplementary_Table_2.tsv"
    with open(tsv, "w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
    subprocess.run([config["tools"]["rscript"], "-e",
                    f"x <- read.delim('{tsv}', check.names = FALSE, stringsAsFactors = FALSE, na.strings = ''); "
                    f"openxlsx::write.xlsx(list('Supplementary Table 2' = x), '{out_dir / 'Supplementary_Table_2.xlsx'}', "
                    f"keepNA = TRUE, na.string = 'NA', overwrite = TRUE)"], check=True)
    tsv.unlink()
    print(f"-> {out_dir / 'Supplementary_Table_2.xlsx'}")


if __name__ == "__main__":
    main()
