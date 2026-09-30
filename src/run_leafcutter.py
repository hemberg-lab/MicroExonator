"""Run the conda-managed Python LeafCutter commands for one comparison.

Kept next to the cluster table (the rest of the working folder is deleted):
  --effect-sizes  leafcutter_ds_effect_sizes.txt: per intron, PSI in each group
                  and deltapsi (B relative to the baseline group A)
  --introns       comparison_perind_numers.counts.gz: the introns of every
                  cluster (chrom:start:end:clu_N_strand) with counts per
                  replicate, which maps cluster IDs to coordinates
With --exons (leafcutter-gtf-to-exons output) clusters are labelled with genes.
"""

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

# leafcutter-ds defaults. Its Python port refuses to run unless both groups
# have at least -i and -g samples, so small designs (3 vs 3) lower them to the
# smaller group size. LeafCutter's p-values are calibrated down to n=4 per group.
MIN_SAMPLES_PER_INTRON = 5
MIN_SAMPLES_PER_GROUP = 3

if __package__:
    from src.umbrella_tool_inputs import leafcutter_files
else:
    from umbrella_tool_inputs import leafcutter_files


def run(preflight_path, joined, output, log, threads=1, effect_sizes=None, introns=None, exons=None):
    output = Path(output)
    log = Path(log)
    kept = [Path(path) for path in (effect_sizes, introns) if path]
    output.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(preflight_path) as stream:
        preflight = json.load(stream)
    if not preflight["inference_supported"]:
        output.write_text("# inference not supported\n")
        for path in kept:
            path.write_text("# inference not supported\n")
        log.write_text("LeafCutter skipped: comparison preflight does not support inference\n")
        return
    with tempfile.TemporaryDirectory(prefix="leafcutter-", dir=output.parent) as directory:
        work = Path(directory)
        leafcutter_files(joined, preflight, work)
        smallest = min(len(preflight["replicates"][side]) for side in ("a", "b"))
        per_intron = min(MIN_SAMPLES_PER_INTRON, smallest)
        per_group = min(MIN_SAMPLES_PER_GROUP, smallest)
        commands = [
            ["leafcutter-cluster", "--juncfiles", "juncfiles.txt",
             "--outprefix", "comparison", "--maxintronlen", "500000",
             "--rundir", "."],
            ["leafcutter-ds", "comparison_perind_numers.counts.gz", "groups.txt",
             "--baseline_group", "a", "--num_threads", str(threads),
             "--min_samples_per_intron", str(per_intron),
             "--min_samples_per_group", str(per_group),
             "--output_prefix", "leafcutter_ds"]
            + (["--exon_file", str(Path(exons).resolve())] if exons else []),
        ]
        with log.open("w") as report:
            if smallest < 4:
                report.write("Note: smallest group has {} replicates; LeafCutter p-values are "
                             "calibrated down to 4 per group\n".format(smallest))
            for command in commands:
                report.write("$ {}\n".format(" ".join(command)))
                report.flush()
                subprocess.run(command, cwd=work, stdout=report,
                               stderr=subprocess.STDOUT, check=True)
        result = work / "leafcutter_ds_cluster_significance.txt"
        if not result.is_file():
            raise RuntimeError("leafcutter-ds did not produce {}".format(result.name))
        shutil.copyfile(result, output)
        if effect_sizes:
            shutil.copyfile(work / "leafcutter_ds_effect_sizes.txt", effect_sizes)
        if introns:
            shutil.copyfile(work / "comparison_perind_numers.counts.gz", introns)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", required=True)
    parser.add_argument("--joined", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--effect-sizes", help="where to keep leafcutter_ds_effect_sizes.txt")
    parser.add_argument("--introns", help="where to keep the cluster intron counts (.counts.gz)")
    parser.add_argument("--exons", help="exon table from leafcutter-gtf-to-exons, to name genes")
    args = parser.parse_args(argv)
    run(args.preflight, args.joined, args.output, args.log, args.threads,
        args.effect_sizes, args.introns, args.exons)


if __name__ == "__main__":
    main()
