#!/bin/bash
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no -t 1 script.gh_probe_small.jl 8 50 100 > ../../data/gh_probe/log_L8.txt 2>&1
julia --startup-file=no -t 1 script.gh_probe_small.jl 10 50 100 > ../../data/gh_probe/log_L10.txt 2>&1
powershell -NoProfile -ExecutionPolicy Bypass -File 'C:\Users\Ritam\.jobq\jobq.ps1' done gh_probe_small
