#!/bin/bash

conda activate parallel_kofam_scan

mkdir -p ~/sdb/parallel_kofam_scan_test_outs/main

rm -rf ~/sdb/parallel_kofam_scan_test_outs/main/* && \
    ~/miniconda3/envs/parallel_kofam_scan/bin/python /home/navi/My-Github-Projects/Parallel-kofamscan/parallel_kofam_scan/pkofamscan \
        -i /home/navi/My-Github-Projects/Parallel-kofamscan/test/data/Sample-1.fna \
        -o ~/sdb/parallel_kofam_scan_test_outs/main \
        -db /mnt/kofamscan-db/ \
        -p 56