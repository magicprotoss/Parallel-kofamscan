#!/usr/bin/env python3

import argparse
import hashlib
import os
from pathlib import Path
from subprocess import DEVNULL, run

import pyarrow as pa
import pyarrow.parquet as pq


BATCH_SIZE = 8192


def parse_args():
    parser = argparse.ArgumentParser(
        prog="Parallel-kofamscan",
        description="Validate inputs and label protein sequences by content hash",
    )
    parser.add_argument("-i", "--unigenes_fp", required=True)
    parser.add_argument("-o", "--output_path", required=True)
    return parser.parse_args()


def validate_input_and_convert_fna_to_faa(unigenes_fp, output_path):
    unigene_faa_fp = Path(output_path) / "unigene.faa"
    seqkit = ["conda", "run", "-n", "parallel-kofamscan.dependency.kofamscan", "seqkit"]
    # An explicit FAA extension prevents nucleotide-alphabet-only proteins from
    # being mistaken for DNA. Other FASTA inputs retain content-based detection.
    input_name = str(unigenes_fp).lower().removesuffix(".gz")
    if not input_name.endswith(".faa"):
        check_fna = run(seqkit + ["seq", str(unigenes_fp), "-v", "-t", "dna", "-o", os.devnull],
                        stdout=DEVNULL, stderr=DEVNULL)
        if check_fna.returncode == 0:
            print(f"Input {unigenes_fp} is FNA; translating with genetic code 11...")
            run(seqkit + ["translate", "-T", "11", str(unigenes_fp), "-o", str(unigene_faa_fp)],
                check=True)
            return "FNA"
    # seqkit writes uncompressed FASTA even when the original FAA is gzipped.
    check_faa = run(seqkit + ["seq", str(unigenes_fp), "-v", "-t", "protein",
                              "-o", str(unigene_faa_fp)], capture_output=True, text=True)
    if check_faa.returncode != 0:
        raise ValueError(f"Input {unigenes_fp} is neither valid FNA nor FAA: {check_faa.stderr}")
    print(f"Input {unigenes_fp} is FAA; proceeding...")
    return "FAA"


def fasta_records(path):
    """Yield the full defline (without '>') and canonical protein sequence."""
    header = None
    fragments = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    if not fragments:
                        raise ValueError(f"{path}: empty sequence for {header!r}")
                    yield header, "".join(fragments).upper()
                header, fragments = line[1:], []
                if not header.strip():
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
            elif line.strip():
                if header is None:
                    raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
                fragments.append("".join(line.split()))
    if header is not None:
        if not fragments:
            raise ValueError(f"{path}: empty sequence for {header!r}")
        yield header, "".join(fragments).upper()


def make_header_map_df(output_path, input_type="FAA"):
    """Stream the map; identical proteins share a hash across all samples."""
    if input_type not in {"FAA", "FNA"}:
        raise ValueError(f"Unknown input type: {input_type}")
    output_path = Path(output_path)
    schema = pa.schema([
        ("sequence_header", pa.string()), ("sequence_hash", pa.string()),
    ], metadata={b"hash_algorithm": b"sha256", b"input_type": input_type.encode("ascii")})
    rows = []
    with pq.ParquetWriter(output_path / "header_map.parquet", schema, compression="zstd") as writer:
        with (output_path / "unigene_hash.faa").open("w", encoding="utf-8") as fasta:
            for header, sequence in fasta_records(output_path / "unigene.faa"):
                sequence_hash = hashlib.sha256(sequence.encode("ascii")).hexdigest()
                fasta.write(f">{sequence_hash}\n{sequence}\n")
                rows.append({"sequence_header": header, "sequence_hash": sequence_hash})
                if len(rows) == BATCH_SIZE:
                    writer.write_table(pa.Table.from_pylist(rows, schema=schema))
                    rows = []
            if rows:
                writer.write_table(pa.Table.from_pylist(rows, schema=schema))
    (output_path / "unigene.faa").unlink()


def main():
    args = parse_args()
    Path(args.output_path).mkdir(parents=True, exist_ok=True)
    input_type = validate_input_and_convert_fna_to_faa(args.unigenes_fp, args.output_path)
    make_header_map_df(args.output_path, input_type)
    (Path(args.output_path) / "done").write_text("done", encoding="utf-8")


if __name__ == "__main__":
    main()
