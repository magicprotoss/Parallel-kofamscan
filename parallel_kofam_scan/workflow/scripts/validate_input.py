#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import pandas as pd
import hashlib
import datetime
from subprocess import run


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Parallel-kofamscan',
        description="A sub-script for validating the input unigenes and format it"
    )

    parser.add_argument('-i', '--unigenes_fp', required=True,
                        help="Path to unigenes to annotate, filenames will be converted to sample-ids")
    parser.add_argument('-o', '--output_path', required=True,
                        help="Directory to save the faa and header map parquet file, path is parsed from snakemake rule")

    return parser.parse_args()


def validate_input_and_convert_fna_to_faa(unigenes_fp, output_path):
    # translate validate and translate fna into faa
    unigene_faa_fp = os.path.join(output_path, 'unigene.faa')
    check_fna = run(
        f"conda run -n parallel-kofamscan.dependency.kofamscan seqkit seq {unigenes_fp} -v -t dna", shell=True, capture_output=True, text=True)
    if check_fna.returncode == 0:
        # if the input is in fna format, convert it to faa format using seqkit
        print(
            f"Input unigenes {unigenes_fp} is in fna format, converting it to faa format using seqkit...")
        fna_to_faa = run(
            f"conda run -n parallel-kofamscan.dependency.kofamscan seqkit translate -T 11 {unigenes_fp} -o {unigene_faa_fp}", shell=True)
        if fna_to_faa.returncode != 0:
            sys.exit(
                f"Error: Failed to convert {unigenes_fp} from fna to faa format")

    elif check_fna.returncode != 0:
        check_faa = run(
            f"conda run -n parallel-kofamscan.dependency.kofamscan seqkit seq {unigenes_fp} -v -t protein", shell=True, capture_output=True, text=True)
        if check_faa.returncode == 0:
            print(
                f"Input unigenes {unigenes_fp} is in faa format, proceeding...")
            shutil.copy(unigenes_fp, unigene_faa_fp)
        else:
            sys.exit(
                f"Error: Input unigenes {unigenes_fp} is neither in fna format nor in faa format, please check your input file")


def make_header_map_df(output_path):
    sample_id = os.path.basename(output_path)
    unigene_faa_fp = os.path.join(output_path, 'unigene.faa')

    # use seqkit to tabulate unigenes
    # col ids
    # ID, sequence, MD5(sequence)
    # headless input

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # the sub-cmd fx2tab is bugged, sometimes it produce an enpty col before the md5 hash
    # and the sub-cmd fx2tab seemed to be the only way for us to get the md5 hash of the sequence
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # screw it, just use uuid
    # get the hased headers first
    original_headers_fp = os.path.join(output_path, 'ori_headers.txt')
    get_original_headers = run(
        f"conda run -n parallel-kofamscan.dependency.kofamscan usearch -fastx_getlabels {unigene_faa_fp} -output {original_headers_fp}", shell=True)
    if get_original_headers.returncode != 0:
        sys.exit(
            f"Error: Failed to get original headers from {unigene_faa_fp} using usearch")
    # relab using uuid, format: >uuid
    unigene_faa_uuid_fp = os.path.join(output_path, 'unigene_uuid.faa')
    relabel_headers = run(
        "conda run -n parallel-kofamscan.dependency.kofamscan seqkit replace -p .+ -r '{uuid}'" + f" {unigene_faa_fp} > {unigene_faa_uuid_fp}", shell=True)
    if relabel_headers.returncode != 0:
        sys.exit(
            f"Error: Failed to relabel headers in {unigene_faa_fp}")
    os.remove(unigene_faa_fp)
    # get uuid headers
    uuid_headers_fp = os.path.join(output_path, 'uuid_headers.txt')
    get_uuid_headers = run(
        f"conda run -n parallel-kofamscan.dependency.kofamscan usearch -fastx_getlabels {unigene_faa_uuid_fp} -output {uuid_headers_fp}", shell=True)
    if get_uuid_headers.returncode != 0:
        sys.exit(
            f"Error: Failed to get uuid headers from {unigene_faa_uuid_fp} using usearch")
    # perp header map df
    with open(original_headers_fp, "rt") as f_ori, open(uuid_headers_fp, "rt") as f_uuid:
        ori_headers = [line.strip() for line in f_ori]
        uuid_headers = [line.strip() for line in f_uuid]
    if len(ori_headers) != len(uuid_headers):
        sys.exit(
            f"Error: The number of original headers and uuid headers do not match, please check the input file and the relabeling step")
    header_map_df = pd.DataFrame({"id": ori_headers, "uuid": uuid_headers})
    # save the header map dataframe as a parquet file
    header_map_parquet_out = os.path.join(output_path, 'header_map.parquet')
    header_map_df.to_parquet(header_map_parquet_out, index=False)
    os.remove(original_headers_fp)
    os.remove(uuid_headers_fp)


def main():
    args = parse_args()

    os.makedirs(args.output_path, exist_ok=True)

    validate_input_and_convert_fna_to_faa(args.unigenes_fp, args.output_path)

    make_header_map_df(args.output_path)

    with open(os.path.join(args.output_path, "done"), "w") as f:
        f.write("done")


if __name__ == "__main__":
    sys.exit(main())
