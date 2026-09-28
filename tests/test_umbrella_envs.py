import os
import re
import unittest


REPOSITORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class UmbrellaEnvTests(unittest.TestCase):
    def test_suppa_rules_use_the_pinned_statsmodels_env(self):
        rules = open(os.path.join(REPOSITORY, "rules", "umbrella_comparisons.smk")).read()
        for rule in ("umbrella_suppa_events", "umbrella_suppa"):
            body = re.search(r"^rule {}:(.*?)(?=^rule |\Z)".format(rule), rules, re.S | re.M).group(1)
            self.assertIn("envs/umbrella-suppa.yaml", body, rule)
        env = open(os.path.join(REPOSITORY, "envs", "umbrella-suppa.yaml")).read()
        # SUPPA 2.3 needs statsmodels.sandbox.stats.multicomp.multipletests, gone in 0.14
        self.assertRegex(env, r"statsmodels=0\.13")
        self.assertNotIn("suppa", open(os.path.join(REPOSITORY, "envs", "umbrella-comparisons.yaml")).read())


if __name__ == "__main__":
    unittest.main()
