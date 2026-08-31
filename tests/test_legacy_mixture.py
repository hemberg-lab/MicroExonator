import os
import tempfile
import unittest


class LegacyMixtureValidationTests(unittest.TestCase):
    def test_too_few_scores_recommend_robustness_filter(self):
        try:
            from src.legacy_mixture import validate_legacy_mixture_scores
        except ImportError as exc:
            self.fail("legacy mixture validation is not implemented: {}".format(exc))

        with self.assertRaises(ValueError) as raised:
            validate_legacy_mixture_scores([72.0])

        message = str(raised.exception)
        self.assertIn("too few microexons", message.lower())
        self.assertIn("filter_method: robustness", message)

    def test_identical_scores_recommend_robustness_filter(self):
        from src.legacy_mixture import validate_legacy_mixture_scores

        with self.assertRaisesRegex(ValueError, "filter_method: robustness"):
            validate_legacy_mixture_scores([72.0, 72.0, 72.0])

    def test_preflight_uses_only_microexons_passing_sample_coverage(self):
        try:
            from src.legacy_mixture import load_legacy_mixture_scores
        except ImportError as exc:
            self.fail("legacy mixture preflight is not implemented: {}".format(exc))

        with tempfile.TemporaryDirectory() as temp_dir:
            coverage = os.path.join(temp_dir, "coverage.tsv")
            matches = os.path.join(temp_dir, "matches.tsv")
            with open(coverage, "w") as handle:
                handle.write("ME\tN_samples\n")
                handle.write("me_pass_1\t2\n")
                handle.write("me_pass_2\t3\n")
                handle.write("me_fail\t1\n")
            with open(matches, "w") as handle:
                handle.write(
                    "ME\tU2_score\tVertebrate_conservation\tME_len\tME_max_U2\n"
                )
                handle.write("me_pass_1\t80\t1.0\t4\t80\n")
                handle.write("me_pass_2\t70\t1.0\t4\t70\n")
                handle.write("me_fail\t40\t0.0\t4\t40\n")

            scores = load_legacy_mixture_scores(
                coverage, matches, min_number_files_detected=2
            )

        self.assertEqual(scores, [80.0, 70.0])


if __name__ == "__main__":
    unittest.main()
