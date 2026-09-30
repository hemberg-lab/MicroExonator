import glob
import json
import os
import re
import unittest


REPOSITORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PBS = os.path.join(REPOSITORY, "Examples", "Cluster_config", "pbs", "cluster.PBS.json")
LSF = os.path.join(REPOSITORY, "Examples", "Cluster_config", "lsf", "cluster.json")


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

    def test_memory_is_given_as_mem_for_pbs_submit(self):
        # pbs-submit.py style profiles read cluster["mem"], not mem_mb; without
        # it every job gets the PBS default memory (often 2 GB)
        config = json.load(open(PBS))
        for rule, entry in config.items():
            if "mem_mb" in entry:
                self.assertEqual(entry.get("mem"), "{}mb".format(int(entry["mem_mb"])), rule)


class LsfClusterConfigTests(unittest.TestCase):
    def test_file_is_valid_json_and_every_key_is_a_rule(self):
        config = json.load(open(LSF))
        threads = rule_threads()
        for rule, entry in config.items():
            if rule == "__default__":
                continue
            self.assertIn(rule, threads, "cluster key matches no rule")
            self.assertGreaterEqual(int(entry.get("nCPUs", 1)), threads[rule], rule)


if __name__ == "__main__":
    unittest.main()
