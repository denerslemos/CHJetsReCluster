#!/usr/bin/env python3
import os
import subprocess
from math import ceil
import argparse

parser = argparse.ArgumentParser(description="Submit ePIC jobs to Condor")
parser.add_argument("--tag", required=True, help="Job tag prefix, e.g. eAu_test")
parser.add_argument("--exec", default="./job.sh", help="Shell script executable (job.sh)")
parser.add_argument("--input-list", default="./input.list", help="Text file containing input ROOT files")
parser.add_argument("--output-dir", default="./results", help="Directory where output ROOT files go")
parser.add_argument("--njobs", type=int, default=1, help="Number of Condor jobs to split input list")
parser.add_argument("--job-args", required=True, help="Arguments after input/output: '0.2,0.4,0.6 1 5'")

args = parser.parse_args()

Exec = args.exec
Tag = args.tag
OutputDir = args.output_dir
JOB_ARGS = args.job_args
NJobs = args.njobs

# Read input file list
with open(args.input_list) as f:
    files = [line.strip() for line in f if line.strip()]

# Prepare directories
os.makedirs(OutputDir, exist_ok=True)
job_folder = os.path.join(OutputDir, "job")
os.makedirs(job_folder, exist_ok=True)

# Split input files across jobs
chunk_size = ceil(len(files) / NJobs)

for i in range(NJobs):
    filenumber = i + 1
    start = i * chunk_size
    end = start + chunk_size
    files_for_job = files[start:end]

    if not files_for_job:
        continue

    # Create per-job input list
    job_input_txt = os.path.join(job_folder, f"{Tag}_{filenumber}_input.txt")
    with open(job_input_txt, "w") as ftxt:
        ftxt.write("\n".join(files_for_job))

    # Output ROOT file
    output_file = os.path.join(OutputDir, f"{Tag}_{filenumber}.root")

    # Condor log files
    job_tag = f"{Tag}_{filenumber}"
    OutFile = os.path.join(job_folder, f"{job_tag}.out")
    ErrFile = os.path.join(job_folder, f"{job_tag}.err")
    LogFile = os.path.join(job_folder, f"{job_tag}.log")
    condor_file = os.path.join(job_folder, f"CondorFile_{job_tag}")

    # Arguments passed to job.sh
    # Format: input.txt output.root RVALS removeelectrons nhitcut
    args_str = f"{job_input_txt} {output_file} {JOB_ARGS}"

    # Write Condor submit file
    with open(condor_file, "w") as f:
        f.write("Universe       = vanilla\n")
        f.write(f"Executable     = {Exec}\n")
        f.write("getenv         = true\n")
        f.write("request_memory = 4G\n")
        f.write(f"Arguments      = {args_str}\n")
        f.write(f"Output         = {OutFile}\n")
        f.write(f"Error          = {ErrFile}\n")
        f.write(f"Log            = {LogFile}\n")
        f.write("Queue\n")

    print(f"Submitting job {filenumber}: {Exec} {args_str}")
    subprocess.run(["condor_submit", condor_file], check=True)

