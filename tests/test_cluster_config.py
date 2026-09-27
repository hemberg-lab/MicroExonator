import glob
import json
import os
import re
import unittest


REPOSITORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PBS = os.path.join(REPOSITORY, "Examples", "Cluster_config", "pbs", "cluster.PBS.json")


def rule_threads():
    """rule name -> highest `threads:` among its definitions (1 when unset)."""
    threads = {}
    for path in glob.glob(os.path.join(REPOSITORY, "rules", "*.smk")) + [os.path.join(REPOSITORY, "MicroExonator.smk")]:
        text = open(path).read()
        for match in re.finditer(r"^\s*rule\s+(\w+)\s*:(.*?)(?=^\s*rule\s+\w+\s*:|\Z)", text, re.S | re.M):
            value = re.search(r"^\s*threads:\s*(\d+)", match.group(2), re.M)
            threads[match.group(1)] = max(threads.get(match.group(1), 1), int(value.group(1)) if value else 1)
    return threads


class PbsClusterConfigTests(unittest.TestCase):
    def test_every_key_is_a_rule_with_enough_cores(self):
        config = json.load(open(PBS))
        threads = rule_threads()
        for rule, entry in config.items():
            if rule == "__default__":
                continue
            self.assertIn(rule, threads, "cluster key matches no rule")
            self.assertGreaterEqual(int(entry.get("ppn", 1)), threads[rule], rule)

    def test_every_threaded_umbrella_rule_is_listed(self):
        config = json.load(open(PBS))
        missing = [rule for rule, count in rule_threads().items()
                   if rule.startswith("umbrella_") and count > 1 and rule not in config]
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main()
