"""VM-only integration check using the real installed SeqKit commands."""

import gzip
import hashlib
from pathlib import Path
import sys

import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "parallel_kofam_scan/workflow/scripts"))
from validate_input import validate_input_and_convert_fna_to_faa, make_header_map_df
from concat_shuffle_split import concat_shuffle_split


def main(root):
    root.mkdir(parents=True, exist_ok=False)
    dna = root / "equivalent coding sequence.fna"
    protein = root / "equivalent protein.faa.gz"
    dna.write_text(">DNA full header # description\nATGAAAACTGCTGCTTAA\n")
    with gzip.open(protein, "wt") as handle:
        handle.write(">protein full header # description\nmkt\naa*\n")
    expected_hash = hashlib.sha256(b"MKTAA*").hexdigest()
    for number, (source, expected_type, header) in enumerate([
        (dna, "FNA", "DNA full header # description"),
        (protein, "FAA", "protein full header # description"),
    ]):
        destination = root / "validate_input" / f"sample{number}"
        destination.mkdir(parents=True)
        input_type = validate_input_and_convert_fna_to_faa(source, destination)
        assert input_type == expected_type, input_type
        make_header_map_df(destination, input_type)
        table = pq.read_table(destination / "header_map.parquet")
        assert table.to_pylist() == [{"sequence_header": header, "sequence_hash": expected_hash}]
        assert table.schema.metadata[b"input_type"].decode() == expected_type
    split = root / "split"
    split.mkdir()
    concat_shuffle_split(str(root / "validate_input"), 10, str(split), 2)
    parts = list((split / "all_hashed_unigenes_no_dup_split").glob("stdin.part_*.fasta"))
    assert len(parts) == 1, parts
    assert parts[0].read_text() == f">{expected_hash}\nMKTAA*\n"
    assert (split / "done").is_file()
    print("PASS: real FNA translation, gzip FAA normalization, shared hashes, and dedup/shuffle/split")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
