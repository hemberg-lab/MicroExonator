import gzip
import json
import os
import re
import subprocess
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

from src.circrna_hisat_delta import boundary_counts, delta, prepare, read_rows


class BoundaryDeltaTests(unittest.TestCase):
    def test_no_cross_boundary_or_strand_pooling(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'junctions.gz'
            with gzip.open(path, 'wt') as out:
                out.write('chrom\tstart\tend\tstrand\tlabel\tmax_anchor\tmulti_total\tshort_anchor_total\tu\n')
                out.write('chr1\t50\t99\t+\tnovel\t20\t0\t0\t7\n')
                out.write('chr1\t201\t250\t+\tnovel\t20\t0\t0\t11\n')
                out.write('chr1\t50\t99\t-\tnovel\t20\t0\t0\t99\n')
            counts = boundary_counts([path], ['u'])
            self.assertEqual(counts['u', 'chr1', 100, '+', 'left'], 7)
            self.assertEqual(counts['u', 'chr1', 200, '+', 'right'], 11)
            self.assertEqual(counts['u', 'chr1', 100, '-', 'left'], 99)
            with self.assertRaisesRegex(ValueError, 'duplicate junction runs'):
                boundary_counts([path, path], ['u'])
            with self.assertRaisesRegex(ValueError, 'missing HISAT2'):
                boundary_counts([path], ['u', 'missing'])

    def test_direction_coverage_and_reproducibility(self):
        a, b = [(80, 20), (85, 15)], [(20, 80), (15, 85)]
        result = delta(a, b)
        self.assertGreater(result['delta'], .5)
        self.assertGreater(result['probability_increase'], .95)
        self.assertEqual(result, delta(a, b))
        self.assertLess(delta(b, a)['delta'], -.5)
        self.assertEqual(delta([(0, 0)]*2, b)['status'], 'insufficient_coverage')
        self.assertEqual(delta(a[:1], b)['status'], 'insufficient_coverage')
        self.assertEqual(delta([(0, 30)]*2, b)['samples_a'], 2)

    def test_identical_groups_not_high_confidence(self):
        result = delta([(10, 20)]*3, [(10, 20)]*3)
        self.assertLess(result['probability'], .55)
        self.assertLess(result['lower'], 0)
        self.assertGreater(result['upper'], 0)

    def test_native_support_to_count_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with gzip.open(root / 'support.tsv.gz', 'wt') as out:
                out.write('circ_id\tchrom\tstart\tend\tstrand\ta1\ta2\tb1\tb2\n')
                out.write('chr1:100|200\tchr1\t100\t200\t+\t1\t2\t0\t4\n')
            preflight = dict(replicates=dict(a=['a1', 'a2'], b=['b1', 'b2']), collapse={})
            prepare('dcc', preflight, root, root / 'diff')
            rows = list(read_rows(root / 'diff/counts.tsv.gz'))
            self.assertEqual(rows[0]['b1'], '0')
            self.assertEqual(list(read_rows(root / 'diff/hisat_delta.tsv.gz')), [])

    def test_minus_strand_and_technical_units_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with gzip.open(root / 'support.tsv.gz', 'wt') as out:
                out.write('circ_id\tchrom\tstart\tend\tstrand\ta1\ta2\tb1\tb2\n')
                out.write('chr1:100|200\tchr1\t100\t200\t-\t10\t10\t1\t1\n')
            with gzip.open(root / 'shard.gz', 'wt') as out:
                out.write('chrom\tstart\tend\tstrand\tlabel\tmax_anchor\tmulti_total\tshort_anchor_total\tu1\tu2\tu3\tu4\tu5\n')
                out.write('chr1\t50\t99\t-\tnovel\t20\t0\t0\t2\t3\t5\t20\t20\n')
                out.write('chr1\t201\t250\t-\tnovel\t20\t0\t0\t7\t7\t14\t20\t20\n')
            preflight = dict(replicates=dict(a=['a1', 'a2'], b=['b1', 'b2']),
                             collapse=dict(u1='a1', u2='a1', u3='a2', u4='b1', u5='b2'))
            prepare('dcc', preflight, root, root / 'diff', [root / 'shard.gz'])
            evidence = {(r['replicate'], r['boundary']): r for r in read_rows(root / 'diff/hisat_evidence.tsv.gz')}
            self.assertEqual(evidence['a1', 'donor']['linear_unique'], '5')
            self.assertEqual(evidence['a1', 'acceptor']['linear_unique'], '14')
            self.assertEqual(evidence['a1', 'donor']['bsj'], '10')
            results = list(read_rows(root / 'diff/hisat_delta.tsv.gz'))
            self.assertEqual(len(results), 2)
            self.assertTrue(all(r['samples_a'] == '2' and r['samples_b'] == '2' for r in results))

    def test_sample_order_does_not_change_draw_assignment(self):
        a, b = [(50, 15), (20, 12)], [(1, 20), (3, 30)]
        self.assertEqual(delta(a, b), delta(a[::-1], b[::-1]))


ENV = os.environ.get("CIRCRNA_CIRIQUANT_ENV")
ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(ENV, "set CIRCRNA_CIRIQUANT_ENV (hisat2 + samtools) for the retained-shard test")
class RetainedShardTests(unittest.TestCase):
    """Real HISAT2 alignments -> the pipeline's own junction tables and shard -> boundary evidence.

    Five technical units (one fixture library split by read) make two biological
    replicates per side, A1 merging two units. Both circles are checked: C1 on +
    (left boundary = acceptor) and C2 on - (left boundary = donor)."""

    @classmethod
    def setUpClass(cls):
        from tests import circrna_fixture as fixture
        cls.tmp = tempfile.TemporaryDirectory()
        work = cls.work = Path(cls.tmp.name)
        env = cls.env = dict(os.environ, PATH=str(Path(ENV) / "bin") + os.pathsep + os.environ["PATH"])
        fx = cls.fx = fixture.build(work / "fx")
        units = cls.units = ["u1", "u2", "u3", "u4", "u5"]
        pairs = []
        with gzip.open(fx["r1"], "rt") as one, gzip.open(fx["r2"], "rt") as two:
            lines1, lines2 = one.read().splitlines(True), two.read().splitlines(True)
        for i in range(0, len(lines1), 4):
            pairs.append((lines1[i:i + 4], lines2[i:i + 4]))
        subprocess.run(["hisat2-build", "-q", str(fx["genome"]), str(work / "idx")], env=env, check=True)
        cls.tables = {}
        for k, unit in enumerate(units):
            mine = pairs[k::len(units)]
            for mate, index in (("R1", 0), ("R2", 1)):
                with gzip.open(work / "{}_{}.fq.gz".format(unit, mate), "wt") as stream:
                    stream.write("".join("".join(pair[index]) for pair in mine))
            bam = work / (unit + ".bam")
            subprocess.run("hisat2 -x {0}/idx -1 {0}/{1}_R1.fq.gz -2 {0}/{1}_R2.fq.gz 2> /dev/null | samtools sort -o {2} - "
                           "&& samtools index {2}".format(work, unit, bam), shell=True, env=env, check=True)
            table = work / (unit + ".junctions.tsv.gz")
            subprocess.run("samtools view {} | python3 src/summarize_junctions.py extract --table {} --summary {}".format(
                bam, table, work / (unit + ".summary.json")), shell=True, cwd=ROOT, env=env, check=True)
            cls.tables[unit] = table
        # intron catalog from the fixture annotation (HISAT2 hints: 0-based last/first exon bases)
        hints = []
        for gene, chrom, strand, exons in (("G1", "chr1", "+", fixture.G1_EXONS), ("G2", "chr2", "-", fixture.G2_EXONS)):
            ordered = sorted(exons)
            for (_, left_end), (right_start, _) in zip(ordered, ordered[1:]):
                hints.append("{}\t{}\t{}\t{}\n".format(chrom, left_end - 1, right_start - 1, strand))
        (work / "hints.txt").write_text("".join(hints))
        subprocess.run(["python3", "src/summarize_junctions.py", "catalog", "--whippet-gtf", str(fx["gtf"]),
                        "--annotation-gtf", str(fx["gtf"]), "--splice-sites", str(work / "hints.txt"),
                        "--output", str(work / "catalog.tsv.gz")], cwd=ROOT, env=env, check=True)
        cls.shard = work / "group.junctions.tsv.gz"
        subprocess.run(["python3", "src/summarize_junctions.py", "shard", "--catalog", str(work / "catalog.tsv.gz"),
                        "--prefix", str(work / "group"), "--reference-id", "ref", "--manifest-sha256", "x" * 64]
                       + ["--table={}={}".format(u, t) for u, t in cls.tables.items()], cwd=ROOT, env=env, check=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def shard_rows(self):
        with gzip.open(self.shard, "rt") as stream:
            return list(csv_reader(stream))

    def test_shard_uses_the_conventions_the_adapter_assumes(self):
        rows = {(r["chrom"], int(r["start"]), int(r["end"])): r for r in self.shard_rows()}
        # introns flanking C1 (+) and C2 (-) are 1-based closed, strand from XS
        self.assertEqual(rows[("chr1", 401, 1000)]["strand"], "+")
        self.assertEqual(rows[("chr1", 1801, 2600)]["strand"], "+")
        self.assertEqual(rows[("chr2", 1401, 2200)]["strand"], "-")
        self.assertEqual(rows[("chr2", 2401, 3200)]["strand"], "-")

    def test_boundary_counts_equal_the_shard_and_the_alignments(self):
        counts = boundary_counts([self.shard], self.units)
        rows = {(r["chrom"], int(r["start"]), int(r["end"])): r for r in self.shard_rows()}
        expected = {("chr1", 1001, "+", "left"): ("chr1", 401, 1000), ("chr1", 1800, "+", "right"): ("chr1", 1801, 2600),
                    ("chr2", 2201, "-", "left"): ("chr2", 1401, 2200), ("chr2", 2400, "-", "right"): ("chr2", 2401, 3200)}
        for site, intron in expected.items():
            for unit in self.units:
                self.assertEqual(counts[(unit,) + site], int(rows[intron][unit]), (site, unit))
            direct = sum(self.spliced_fragments(unit, intron) for unit in self.units)
            total = sum(counts[(unit,) + site] for unit in self.units)
            self.assertGreater(total, 0, site)
            # anchors < 8 nt are counted as short_anchor by the extractor, never as unique
            self.assertLessEqual(total, direct, site)
            self.assertGreaterEqual(total, 0.75 * direct, site)

    def spliced_fragments(self, unit, intron):
        chrom, first, last = intron
        result = subprocess.run(["samtools", "view", "-F", "0x904", str(self.work / (unit + ".bam")), chrom],
                                env=self.env, capture_output=True, text=True, check=True)
        names = set()
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            pos = int(fields[3])
            for length, op in re.findall(r"(\d+)([MIDNSHP=X])", fields[5]):
                if op == "N" and (pos, pos + int(length) - 1) == (first, last):
                    names.add(fields[0])
                if op in "MDN=X":
                    pos += int(length)
        return len(names)

    def test_merged_units_and_both_strands_reach_the_evidence_table(self):
        root = self.work / "dcc"
        root.mkdir(exist_ok=True)
        with gzip.open(root / "support.tsv.gz", "wt") as out:
            out.write("circ_id\tchrom\tstart\tend\tstrand\tA1\tA2\tB1\tB2\n")
            out.write("chr1:1001|1800\tchr1\t1001\t1800\t+\t9\t5\t4\t6\n")
            out.write("chr2:2201|2400\tchr2\t2201\t2400\t-\t7\t3\t2\t8\n")
        preflight = dict(replicates=dict(a=["A1", "A2"], b=["B1", "B2"]),
                         collapse=dict(u1="A1", u2="A1", u3="A2", u4="B1", u5="B2"))
        prepare("dcc", preflight, root, root / "diff", [self.shard], dict(size=500))
        evidence = {(r["circ_id"], r["boundary"], r["replicate"]): r for r in read_rows(root / "diff/hisat_evidence.tsv.gz")}
        counts = boundary_counts([self.shard], self.units)
        # C1 (+): acceptor = left (1001), donor = right (1800); C2 (-): donor = left (2201)
        self.assertEqual(int(evidence["chr1:1001|1800", "acceptor", "A1"]["linear_unique"]),
                         counts["u1", "chr1", 1001, "+", "left"] + counts["u2", "chr1", 1001, "+", "left"])
        self.assertEqual(int(evidence["chr1:1001|1800", "donor", "B2"]["linear_unique"]),
                         counts["u5", "chr1", 1800, "+", "right"])
        self.assertEqual(int(evidence["chr2:2201|2400", "donor", "A2"]["linear_unique"]),
                         counts["u3", "chr2", 2201, "-", "left"])
        self.assertEqual(evidence["chr2:2201|2400", "donor", "A1"]["bsj"], "7")
        results = list(read_rows(root / "diff/hisat_delta.tsv.gz"))
        self.assertEqual(sorted((r["circ_id"], r["boundary"]) for r in results),
                         [("chr1:1001|1800", "acceptor"), ("chr1:1001|1800", "donor"),
                          ("chr2:2201|2400", "acceptor"), ("chr2:2201|2400", "donor")])
        metadata = json.loads((root / "diff/preparation.json").read_text())
        self.assertEqual(metadata["bayesian_status"], "experimental_cross_aligner")


def csv_reader(stream):
    import csv
    return csv.DictReader(stream, delimiter="\t")


if __name__ == '__main__':
    unittest.main()
