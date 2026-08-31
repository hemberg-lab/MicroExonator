import csv
import gzip
import os
import tempfile
import unittest


class DetectionThresholdTests(unittest.TestCase):
    def test_requires_majority_ci_support_and_three_spanning_reads(self):
        try:
            from src.robustness_filter import detect_group_microexons
        except ImportError as exc:
            self.fail("robustness filtering is not implemented: {}".format(exc))

        with tempfile.TemporaryDirectory() as temp_dir:
            psi_1 = os.path.join(temp_dir, "sample_1.PSI.gz")
            psi_2 = os.path.join(temp_dir, "sample_2.PSI.gz")
            reads = os.path.join(temp_dir, "reads.tsv")

            with gzip.open(psi_1, "wt") as handle:
                handle.write("ME\tCI_Lo\n")
                handle.write("me_pass\t0.20\n")
                handle.write("me_half\t0.20\n")
                handle.write("me_low_ci\t0.05\n")
                handle.write("me_low_reads\t0.20\n")
            with gzip.open(psi_2, "wt") as handle:
                handle.write("ME\tCI_Lo\n")
                handle.write("me_pass\t0.10\n")
                handle.write("me_half\t0.05\n")
                handle.write("me_low_ci\t0.05\n")
                handle.write("me_low_reads\t0.10\n")
            with open(reads, "w") as handle:
                handle.write("ME\tSpanning_cov\n")
                handle.write("me_pass\t3\n")
                handle.write("me_half\t10\n")
                handle.write("me_low_ci\t10\n")
                handle.write("me_low_reads\t2\n")

            detected = detect_group_microexons([psi_1, psi_2], [reads])

        self.assertEqual(
            detected,
            {
                "me_pass": {
                    "total_measurements": 2,
                    "detected_samples": 2,
                    "spanning_reads": 3,
                }
            },
        )


class FinalRobustnessFilterTests(unittest.TestCase):
    def test_min_detected_samples_is_applied_within_each_group(self):
        try:
            from src.robustness_filter import collect_robust_microexons
        except ImportError as exc:
            self.fail("final robustness filtering is not implemented: {}".format(exc))

        with tempfile.TemporaryDirectory() as temp_dir:
            detected_file = os.path.join(temp_dir, "control.detected.txt")
            with open(detected_file, "w") as handle:
                handle.write(
                    "ME\ttotal_measurements\tdetected_samples\tspanning_reads\n"
                )
                handle.write("me_pass\t3\t2\t5\n")
                handle.write("me_too_few_samples\t3\t1\t8\n")

            robust = collect_robust_microexons(
                [detected_file], min_detected_samples=2
            )

        self.assertEqual(robust, {"me_pass"})

    def test_robust_output_preserves_all_twelve_me_centric_columns(self):
        try:
            from src.robustness_filter import write_robust_microexons
        except ImportError as exc:
            self.fail("robust output writing is not implemented: {}".format(exc))

        with tempfile.TemporaryDirectory() as temp_dir:
            me_centric = os.path.join(temp_dir, "TOTAL.ME_centric.txt")
            output = os.path.join(temp_dir, "out.robustly_detected.txt")
            with open(me_centric, "w") as handle:
                handle.write(
                    "me_pass\ttx1\t10\tsj1\t5\t4\tACGT\t2\t80\t1.5\t0.01\tmatch_a,match_b\n"
                )
                handle.write(
                    "me_filtered\ttx2\t4\tsj2\t2\t5\tAACGT\t1\t60\t0.5\t0.02\tmatch_c\n"
                )

            write_robust_microexons(me_centric, output, {"me_pass"})

            with open(output) as handle:
                rows = list(csv.reader(handle, delimiter="\t"))

        self.assertEqual(
            rows[0],
            [
                "ME",
                "Transcript",
                "Total_coverage",
                "Total_SJs",
                "ME_coverages",
                "ME_length",
                "ME_seq",
                "ME_matches",
                "U2_score",
                "Mean_conservation",
                "P_MEs",
                "Total_ME",
            ],
        )
        self.assertEqual(
            rows[1],
            [
                "me_pass",
                "tx1",
                "10",
                "sj1",
                "5",
                "4",
                "ACGT",
                "2",
                "80",
                "1.5",
                "0.01",
                "match_a,match_b",
            ],
        )
        self.assertEqual(len(rows), 2)


if __name__ == "__main__":
    unittest.main()
