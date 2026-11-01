#!/usr/bin/env python3
"""Snakemake 7 --cluster-status for Slurm: prints running, success or failed for one job id.

Without it Snakemake only sees jobs that end by themselves; a job Slurm kills for time or
memory never writes its marker and Snakemake waits for it forever.
"""
import subprocess
import sys
import time

RUNNING = {"PENDING", "CONFIGURING", "RUNNING", "COMPLETING", "SUSPENDED", "REQUEUED",
           "REQUEUE_HOLD", "REQUEUE_FED", "RESIZING", "SIGNALING", "STAGE_OUT", "RESV_DEL_HOLD"}


def sacct_state(jobid):
    out = subprocess.run(["sacct", "-j", jobid, "-X", "-n", "-P", "-o", "State"],
                         capture_output=True, text=True).stdout.split()
    return out[0] if out else ""


def scontrol_state(jobid):
    out = subprocess.run(["scontrol", "show", "job", jobid], capture_output=True, text=True).stdout
    for field in out.split():
        if field.startswith("JobState="):
            return field.split("=", 1)[1]
    return ""


jobid = sys.argv[1]
state = ""
for attempt in range(5):
    state = sacct_state(jobid) or scontrol_state(jobid)
    if state:
        break
    time.sleep(2)

state = state.split()[0].rstrip("+") if state else ""   # "CANCELLED by 123" / "CANCELLED+"
if not state or state in RUNNING:
    print("running")   # just submitted and not yet in the accounting database
elif state == "COMPLETED":
    print("success")
else:
    print("failed")    # FAILED, TIMEOUT, OUT_OF_MEMORY, CANCELLED, NODE_FAIL, ...
