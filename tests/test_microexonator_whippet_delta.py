"""The legacy MegaSearch delta: Whippet delta with MicroExonator PSI in microexon nodes."""

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

from src import run_comparison_tools

FAKE_DELTA = '''
import gzip, sys
args = dict(zip(sys.argv[1::2], sys.argv[2::2]))
def psi(paths):
    values = {}
    for path in paths.split(","):
        with gzip.open(path, "rt") as stream:
            header = next(stream).rstrip("\\n").split("\\t")
            for line in stream:
                row = dict(zip(header, line.rstrip("\\n").split("\\t")))
                values.setdefault((row["Gene"], row["Node"], row["Coord"], row["Strand"], row["Type"]), []).append(float(row["Psi"]))
    return {key: sum(v) / len(v) for key, v in values.items()}
a, b = psi(args["-a"]), psi(args["-b"])
with gzip.open(args["-o"] + ".diff.gz", "wt") as out:
    out.write("Gene\\tNode\\tCoord\\tStrand\\tType\\tPsi_A\\tPsi_B\\tDeltaPsi\\tProbability\\tComplexity\\tEntropy\\n")
    for key in sorted(a):
        out.write("\\t".join(list(key) + ["%.2f" % a[key], "%.2f" % b[key], "%.2f" % (a[key] - b[key]), "0.99", "K1", "0.0"]) + "\\n")
'''


class MicroexonatorWhippetDeltaTests(unittest.TestCase):
    def test_microexon_nodes_carry_microexonator_psi_and_only_they_are_reported(self):
        runs = {"a": ["a1", "a2"], "b": ["b1", "b2"]}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fields = ("Psi", "CI_Width", "CI_Lo,Hi", "Total_Reads", "Complexity", "Entropy")
            with gzip.open(root / "whippet.tsv.gz", "wt") as stream:
                stream.write("\t".join(["Gene", "Node", "Coord", "Strand", "Type"] +
                                       ["{}.{}".format(r, f) for r in runs["a"] + runs["b"] for f in fields]) + "\n")
                for node, coord in (("1", "chr1:101-110"), ("2", "chr1:201-300")):
                    stream.write("\t".join(["G", node, coord, "+", "CE"] +
                                           ["0.2", "0.1", "0.15,0.25", "20", "K1", "0.0"] * 4) + "\n")
            me_fields = ("ME_coverages", "excluding_covs", "PSI", "CI_Lo", "CI_Hi")
            with gzip.open(root / "microexonator.tsv.gz", "wt") as stream:
                stream.write("\t".join(["ME"] + ["{}.{}".format(r, f) for r in runs["a"] + runs["b"] for f in me_fields]) + "\n")
                stream.write("\t".join(["chr1_+_100_110"] + ["9", "1", "0.9", "0.8", "0.95"] * 2
                                       + ["1", "9", "0.1", "0.05", "0.2"] * 2) + "\n")
            (root / "fixed.tsv").write_text("ME\tTranscript\nchr1_+_100_110\tT1\n")
            with gzip.open(root / "exons.tab.gz", "wt") as stream:
                stream.write("Gene\tPotential_Exon\tIs_Annotated\tWhippet_Nodes\n"
                             "G\tchr1:101-110:+\tY\t1\nG\tchr1:201-300:+\tY\t2\n")
            bin_dir = root / "bin"
            bin_dir.mkdir()
            (bin_dir / "whippet-delta.jl").write_text(FAKE_DELTA)
            preflight = root / "preflight.json"
            preflight.write_text(json.dumps({"inference_supported": True, "runs": runs}))
            full, microexons = root / "me_whippet.diff.gz", root / "me_whippet.microexons.tsv"
            run_comparison_tools.main([
                "microexonator_whippet", "--preflight", str(preflight), "--work", str(root / "work"),
                "--log", str(root / "log"), "--outputs", str(full), str(microexons),
                "--joined", str(root / "microexonator.tsv.gz"), "--joined-whippet", str(root / "whippet.tsv.gz"),
                "--whippet-exons", str(root / "exons.tab.gz"), "--microexons", str(root / "fixed.tsv"),
                "--julia", sys.executable, "--whippet-bin", str(bin_dir)])
            with gzip.open(full, "rt") as stream:
                rows = {line.split("\t")[1]: line.rstrip("\n").split("\t") for line in list(stream)[1:]}
            # node 1 is the microexon: MicroExonator PSI (0.9 vs 0.1); node 2 keeps Whippet's 0.2
            self.assertEqual(rows["1"][5:8], ["0.90", "0.10", "0.80"])
            self.assertEqual(rows["2"][5:7], ["0.20", "0.20"])
            lines = microexons.read_text().splitlines()
            self.assertEqual(lines[0].split("\t")[:2], ["exon_ID", "Gene"])
            self.assertEqual(len(lines), 2)
            self.assertEqual(lines[1].split("\t")[0], "chr1_+_100_110")
            self.assertFalse((root / "work").exists())


if __name__ == "__main__":
    unittest.main()
