#!/usr/bin/env python3

import os
import argparse
import sys
from subprocess import run
from glob import glob

# for a given input dir, glob all the "unigenes_hashed.faa" files from sample dirs under the input dir
# concatenate them into a single file, shuffle the lines, and split into multiple files using seqkit


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Parallel-kofamscan',
        description="A sub-script for spliting unigenes into chunks"
    )

    parser.add_argument('-i', '--validate_input_dir', required=True,
                        help="Path to the outdir of the validate_input rule, which contains the hashed faa files stored in sample-id subdirs")
    parser.add_argument('-c', '--chunk_size', required=True,
                        help="Number of sequences in each chunk, passed to seqkit split2")
    parser.add_argument('-t', '--threads', required=False, default=1,
                        help="Number of threads to use -fastx_unique")
    parser.add_argument('-o', '--output_path', required=True,
                        help="Directory for this script")

    return parser.parse_args()


def concat_shuffle_split(validate_input_dir, chunk_size, output_path, threads):
    # glob all the "unigenes_hashed.faa" files from sample dirs under the input dir

    # get fps
    individual_hashed_unigenes_fps = glob(
        os.path.join(validate_input_dir, "*", "unigene_uuid.faa"))

    # concat all inputs
    all_hashed_unigenes_fp = os.path.join(
        output_path, "all_hashed_unigenes.faa")
    with open(all_hashed_unigenes_fp, "wt") as outfile:
        for fp in individual_hashed_unigenes_fps:
            with open(fp, "rt") as infile:
                outfile.write(infile.read())
            os.remove(fp)

    # we need to make sure no duplicate in the concatenated file
    all_hashed_unigenes_no_dup_fp = os.path.join(
        output_path, "all_hashed_unigenes_no_dup.faa")
    remove_dup = run(
        f"conda run -n parallel-kofamscan.dependency.kofamscan usearch -fastx_uniques {all_hashed_unigenes_fp} -fastaout {all_hashed_unigenes_no_dup_fp} -threads {int(threads)}", shell=True)
    if remove_dup.returncode != 0:
        sys.exit(
            f"Error: Failed to remove duplicate sequences using usearch, please check the input files and parameters")
    os.remove(all_hashed_unigenes_fp)

    # run seqkit shuffle and split
    all_hashed_unigenes_no_dup_split_dir = os.path.join(
        output_path, "all_hashed_unigenes_no_dup_split")
    # seqkit will create this dir
    get_seqkit_exec_fp = run(
        f"conda run -n parallel-kofamscan.dependency.kofamscan which seqkit", shell=True, capture_output=True, text=True)
    if get_seqkit_exec_fp.returncode != 0:
        sys.exit(
            f"Error: Failed to get seqkit executable path, please check your conda environment and seqkit installation")
    seqkit_exec_fp = get_seqkit_exec_fp.stdout.strip()
    concat_shuffle_split = run(
        f"{seqkit_exec_fp} shuffle --two-pass --non-deterministic {all_hashed_unigenes_no_dup_fp} | {seqkit_exec_fp} split2 --by-size {chunk_size} --out-dir {all_hashed_unigenes_no_dup_split_dir}", shell=True)
    if concat_shuffle_split.returncode != 0:
        sys.exit(f"Error: Failed to concatenate, shuffle and split unigenes using seqkit, please check the input files and parameters")
    os.remove(all_hashed_unigenes_no_dup_fp)
    os.remove(all_hashed_unigenes_no_dup_fp + ".seqkit.fai")

    # we'll create a "done" file in the output dir to indicate the process is done, and the downstream kofamscan rule can start
    done_fp = os.path.join(output_path, "done")
    with open(done_fp, "w") as f:
        f.write("done")


def main():
    args = parse_args()

    os.makedirs(args.output_path, exist_ok=True)

    # not need to overkill here
    if int(args.threads) > 20:
        threads = 20
    else:
        threads = args.threads

    concat_shuffle_split(args.validate_input_dir,
                         args.chunk_size, args.output_path, threads)


if __name__ == "__main__":
    sys.exit(main())
