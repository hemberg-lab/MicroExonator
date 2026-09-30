"""Run comparison tools inside their rule-local Conda environments.

Both read the comparison preflight. When the tool cannot run, each declared
output gets one comment line saying why, and no inferential command runs.

  rmats  stage the kept prep files in a fresh --tmp, write b1/b2 with the BAM
         paths recorded at prep time, run --task inte then --task post, keep
         the *.MATS.JC.txt tables
  suppa  per-replicate TPM tables for A and B, psiPerEvent per event type,
         diffSplice with conditions given B then A (so dPSI is A - B)
"""

import argparse
import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

if __package__:
    from src.umbrella_tool_inputs import suppa_tables, write_delta_inputs
else:
    from umbrella_tool_inputs import suppa_tables, write_delta_inputs


UNSUPPORTED = "# inference not supported"


def marker(outputs, text):
    for path in outputs:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(text + "\n")


def run(command, log):
    with open(log, "a") as stream:
        stream.write("$ " + " ".join(command) + "\n")
        stream.flush()
        subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)


def rmats(args, preflight):
    support = preflight["tools"]["rmats"]
    if not support["supported"]:
        marker(args.outputs, UNSUPPORTED)
        return
    work = Path(args.work)
    shutil.rmtree(work, ignore_errors=True)
    (work / "tmp").mkdir(parents=True)
    (work / "od").mkdir()
    records = {}
    for prep in args.prep:
        with open(prep + ".json") as stream:
            record = json.load(stream)
        records[record["run_id"]] = record
        shutil.copyfile(prep, work / "tmp" / Path(prep).name)
    for side, name in (("a", "b1"), ("b", "b2")):
        (work / (name + ".txt")).write_text(",".join(
            records[run_id]["bam"] for run_id in preflight["runs"][side]) + "\n")
    common = ["--b1", str(work / "b1.txt"), "--b2", str(work / "b2.txt"), "--gtf", args.gtf,
              "--od", str(work / "od"), "--tmp", str(work / "tmp"),
              "-t", "paired" if preflight["layouts"] == ["PE"] else "single",
              "--readLength", str(max(r["read_length"] for r in records.values())),
              "--variable-read-length", "--novelSS", "--nthread", str(args.threads)]
    run(["rmats.py"] + common + ["--task", "inte"], args.log)
    run(["rmats.py"] + common + ["--task", "post"], args.log)
    for path in args.outputs:
        shutil.copyfile(work / "od" / Path(path).name, path)
    shutil.rmtree(work)


def suppa(args, preflight):
    if not preflight["inference_supported"]:
        marker(args.outputs, UNSUPPORTED)
        return
    work = Path(args.work)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    tables = suppa_tables(args.joined, preflight)
    for side in ("a", "b"):
        (work / (side + ".tpm")).write_text(tables[side])
    for events, output in zip(args.events, args.outputs):
        kind = Path(output).stem
        for side in ("a", "b"):
            run(["suppa.py", "psiPerEvent", "-i", events, "-e", str(work / (side + ".tpm")),
                 "-o", str(work / "{}_{}".format(kind, side))], args.log)
        run(["suppa.py", "diffSplice", "-m", "empirical", "-gc", "-i", events,
             "-p", str(work / (kind + "_b.psi")), str(work / (kind + "_a.psi")),
             "-e", str(work / "b.tpm"), str(work / "a.tpm"),
             # absolute: SUPPA 2.3 merges its .dpsi.temp.* files from the
             # current directory unless the prefix is absolute, and then writes
             # no .dpsi at all
             "-o", str(Path(output).with_suffix("").resolve())], args.log)
    shutil.rmtree(work)


def microexonator(args, preflight):
    if not preflight["inference_supported"]:
        marker(args.outputs, UNSUPPORTED)
        return
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    paths = {side: [str(work / (run_id + ".tsv.gz")) for run_id in preflight["runs"][side]]
             for side in ("a", "b")}
    try:
        runs = preflight["runs"]["a"] + preflight["runs"]["b"]
        outputs = paths["a"] + paths["b"]
        write_delta_inputs("microexonator", args.joined, runs, outputs)
        command = ["python3", "src/me_delta.py", "-a", ",".join(paths["a"]),
                   "-b", ",".join(paths["b"]), "--microexons", args.microexons]
        command.extend(shlex.split(args.options))
        with open(args.outputs[0], "w") as stream:
            subprocess.run(command, check=True, stdout=stream)
    finally:
        shutil.rmtree(work)


def microexonator_whippet(args, preflight):
    """The legacy MegaSearch route: Whippet delta on Whippet PSI files whose
    microexon nodes carry MicroExonator's corrected PSI and CI.

    Per run, src/Replace_PSI_whippet2.py swaps MicroExonator's PSI into the
    rebuilt Whippet table (the legacy ME_psi_to_quant rule), whippet-delta.jl
    compares the swapped tables (whippet_delta_ME), and
    src/whippet_delta_to_ME.py keeps the microexon nodes
    (delta_ME_from_MicroExonator). Outputs: the full delta, then the
    microexon rows.
    """
    if not preflight["inference_supported"]:
        marker(args.outputs, UNSUPPORTED)
        return
    work = Path(args.work)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    try:
        runs = preflight["runs"]["a"] + preflight["runs"]["b"]
        whippet = [str(work / (run + ".psi.gz")) for run in runs]
        me = [str(work / (run + ".me.tsv.gz")) for run in runs]
        swapped = {run: str(work / (run + ".psi.ME.gz")) for run in runs}
        write_delta_inputs("whippet", args.joined_whippet, runs, whippet)
        write_delta_inputs("microexonator", args.joined, runs, me)
        for run, me_path, whippet_path in zip(runs, me, whippet):
            run_command(["python3", "src/Replace_PSI_whippet2.py", me_path, whippet_path,
                         swapped[run]], args.log)
        prefix = str(work / "delta")
        run_command([args.julia, str(Path(args.whippet_bin) / "whippet-delta.jl"),
                     "-a", ",".join(swapped[run] for run in preflight["runs"]["a"]),
                     "-b", ",".join(swapped[run] for run in preflight["runs"]["b"]),
                     "-o", prefix], args.log)
        full, microexons = args.outputs
        shutil.copyfile(prefix + ".diff.gz", full)
        with open(microexons, "w") as stream:
            subprocess.run(["python3", "src/whippet_delta_to_ME.py", args.whippet_exons, full,
                            args.microexons], check=True, stdout=stream)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def run_command(command, log):
    run(command, log)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tool", choices=("rmats", "suppa", "microexonator", "microexonator_whippet"))
    parser.add_argument("--preflight", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--outputs", nargs="+", required=True)
    parser.add_argument("--prep", nargs="*", default=[])
    parser.add_argument("--gtf")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--joined")
    parser.add_argument("--microexons")
    parser.add_argument("--options", default="")
    parser.add_argument("--events", nargs="*", default=[])
    parser.add_argument("--joined-whippet")
    parser.add_argument("--whippet-exons", help="the Whippet index's .exons.tab.gz")
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--whippet-bin", default="")
    args = parser.parse_args(argv)
    with open(args.preflight) as stream:
        preflight = json.load(stream)
    {"rmats": rmats, "suppa": suppa, "microexonator": microexonator,
     "microexonator_whippet": microexonator_whippet}[args.tool](args, preflight)


if __name__ == "__main__":
    sys.exit(main())
