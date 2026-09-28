"""Run the conda-managed Python LeafCutter commands for one comparison."""

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

if __package__:
    from src.umbrella_tool_inputs import leafcutter_files
else:
    from umbrella_tool_inputs import leafcutter_files


def run(preflight_path, joined, output, log, threads=1):
    output = Path(output)
    log = Path(log)
    output.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(preflight_path) as stream:
        preflight = json.load(stream)
    if not preflight["inference_supported"]:
        output.write_text("# inference not supported\n")
        log.write_text("LeafCutter skipped: comparison preflight does not support inference\n")
        return
    with tempfile.TemporaryDirectory(prefix="leafcutter-", dir=output.parent) as directory:
        work = Path(directory)
        leafcutter_files(joined, preflight, work)
        commands = [
            ["leafcutter-cluster", "--juncfiles", "juncfiles.txt",
             "--outprefix", "comparison", "--maxintronlen", "500000",
             "--rundir", "."],
            ["leafcutter-ds", "comparison_perind_numers.counts.gz", "groups.txt",
             "--baseline_group", "a", "--num_threads", str(threads),
             "--output_prefix", "leafcutter_ds"],
        ]
        with log.open("w") as report:
            for command in commands:
                report.write("$ {}\n".format(" ".join(command)))
                report.flush()
                subprocess.run(command, cwd=work, stdout=report,
                               stderr=subprocess.STDOUT, check=True)
        result = work / "leafcutter_ds_cluster_significance.txt"
        if not result.is_file():
            raise RuntimeError("leafcutter-ds did not produce {}".format(result.name))
        shutil.copyfile(result, output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", required=True)
    parser.add_argument("--joined", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)
    run(args.preflight, args.joined, args.output, args.log, args.threads)


if __name__ == "__main__":
    main()
