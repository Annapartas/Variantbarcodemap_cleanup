#!/usr/bin/env python3
"""
BC_FC_variant_filter.py

Input:
    A plain text or gzipped file with two columns:
        DNA_code    Variant_sequence

    DNA_code can be either:
        - barcode: 33-bp DNA barcode, not translated
        - flycode: variable-length DNA flycode, translated to peptide

The library type is specified in the JSON parameter file:
    "library_type": "barcoded"
or:
    "library_type": "flycoded"

Output:
    1. *_filtered.txt.gz
       Same two-column format as input:
           DNA_code    Variant_sequence

    2. *_peptide_variantcode.tsv
       DNA_code, Peptide_code, Variant_code, Variant_sequence

    3. *_peptide_variantcode_filteredOUT.tsv
       DNA_code, Peptide_code, Variant_sequence, Reason

Variant filtering:
    - exact WT DNA sequence is accepted and marked as Variant_code = wt
    - exactly one codon substitution is accepted only if:
        * the amino-acid position is in allowed_mutation_positions
        * the mutant codon is in allowed_substitution_codons
        * the substitution does not create a stop codon
    - no additional mutations are allowed
    - the same DNA code may not be assigned to more than one accepted variant

Barcode/flycode filtering:
    - barcoded library:
        barcode_regex = ^([ATGC]{3}[AT][GC]){6}[ATGC]{3}$
        code must be 33 bp
        code is not translated; Peptide_code remains empty

    - flycoded library:
        barcode_regex = ^GGTAGT(GCA|GTT|GAT|CCA|GAA|ACT|GGT|TCT|TAC|CTG|TGG|CAG|TTC|AAC){6,8}(TGGCGG|TGGCTGCGG|TGGCAGTCTCGG|TGGCAGGAAGGAGGTCGG)$
        code is translated to peptide

Requirements:
    Python 3 only.

Run:
    python BC_FC_variant_filter.py --parameters params.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import TextIO


BARCODE_REGEX = r"^([ATGC]{3}[AT][GC]){6}[ATGC]{3}$"

FLYCODE_REGEX = (
    r"^GGTAGT"
    r"(GCA|GTT|GAT|CCA|GAA|ACT|GGT|TCT|TAC|CTG|TGG|CAG|TTC|AAC){6,8}"
    r"(TGGCGG|TGGCTGCGG|TGGCAGTCTCGG|TGGCAGGAAGGAGGTCGG)$"
)

GENETIC_CODE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C",
    "TGA": "*", "TGG": "W",

    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G",
    "GGA": "G", "GGG": "G",
}


def fail(message: str) -> None:
    raise ValueError(message)


def open_text(path: Path, mode: str) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t")
    return path.open(mode, encoding="utf-8")


def load_json(path: Path) -> dict:
    if not path.exists():
        fail(f"Parameter file not found: {path}")

    try:
        with path.open("r", encoding="utf-8") as handle:
            config = json.load(handle)
    except json.JSONDecodeError as error:
        fail(
            f"Invalid JSON in {path}, line {error.lineno}, "
            f"column {error.colno}: {error.msg}"
        )

    if not isinstance(config, dict):
        fail("Top level of parameter file must be a JSON object.")

    return config


def clean_dna(sequence: str, name: str) -> str:
    sequence = str(sequence).strip().upper().replace(" ", "")

    if not sequence:
        fail(f"{name} is empty.")

    invalid = sorted(set(sequence) - set("ACGT"))
    if invalid:
        fail(
            f"{name} contains invalid DNA characters: "
            + ", ".join(invalid)
        )

    return sequence


def compile_regex(pattern: str, name: str) -> re.Pattern[str]:
    try:
        return re.compile(pattern, re.IGNORECASE)
    except re.error as error:
        fail(f"Invalid {name}: {error}")


def translate_dna(
    dna: str,
    sequence_name: str,
    allow_terminal_stop: bool = False,
) -> str:
    dna = clean_dna(dna, sequence_name)

    if len(dna) % 3 != 0:
        fail(
            f"{sequence_name} length is not divisible by 3: "
            f"{len(dna)} nt."
        )

    amino_acids = []

    for index in range(0, len(dna), 3):
        codon = dna[index:index + 3]
        aa = GENETIC_CODE.get(codon)

        if aa is None:
            fail(f"{sequence_name} contains invalid codon {codon}.")

        if aa == "*":
            codon_number = index // 3 + 1

            if allow_terminal_stop and index == len(dna) - 3:
                continue

            fail(
                f"{sequence_name} contains stop codon {codon} "
                f"at codon {codon_number}."
            )

        amino_acids.append(aa)

    return "".join(amino_acids)


def parse_input_line(
    line: str,
    line_number: int,
    delimiter: str | None,
) -> tuple[str, str]:
    line = line.rstrip("\n\r")

    if delimiter is None or delimiter == "auto":
        parts = line.split()
    else:
        parts = line.split(delimiter)

    if len(parts) < 2:
        fail(
            f"Line {line_number} has fewer than two columns. "
            "Expected: DNA_code Variant_sequence"
        )

    return parts[0].strip(), parts[1].strip()


def default_output_prefix(input_file: Path) -> str:
    name = input_file.name

    for suffix in (".txt.gz", ".tsv.gz", ".gz", ".txt", ".tsv"):
        if name.endswith(suffix):
            return name[:-len(suffix)]

    return input_file.stem


def normalize_positions(values) -> set[int]:
    if not isinstance(values, list) or not values:
        fail("allowed_mutation_positions must be a non-empty list of integers.")

    positions = set()

    for value in values:
        position = int(value)

        if position < 1:
            fail("Allowed mutation positions use 1-based numbering and must be >= 1.")

        positions.add(position)

    return positions


def normalize_codons(values) -> set[str]:
    if not isinstance(values, list) or not values:
        fail("allowed_substitution_codons must be a non-empty list of codons.")

    codons = set()

    for value in values:
        codon = clean_dna(str(value), "allowed substitution codon")

        if len(codon) != 3:
            fail(f"Allowed substitution codon must be 3 nt: {codon}")

        if GENETIC_CODE[codon] == "*":
            fail(f"Stop codons are not allowed in allowed_substitution_codons: {codon}")

        codons.add(codon)

    return codons


def get_library_settings(config: dict) -> dict:
    library_type = str(config.get("library_type", "")).strip().lower()

    if library_type not in {"barcoded", "flycoded"}:
        fail('library_type must be either "barcoded" or "flycoded".')

    code_regex_override = config.get("barcode_regex")

    if library_type == "barcoded":
        expected_regex = BARCODE_REGEX
        translate_code_default = False
        required_code_length_default = 33
        code_name = "barcode"
    else:
        expected_regex = FLYCODE_REGEX
        translate_code_default = True
        required_code_length_default = None
        code_name = "flycode"

    if code_regex_override is None:
        barcode_regex = expected_regex
    else:
        barcode_regex = str(code_regex_override).strip()
        if barcode_regex != expected_regex:
            print(
                "WARNING: barcode_regex in the parameter file overrides the "
                f"default regex for library_type={library_type}.",
                file=sys.stderr,
            )

    translate_code = bool(config.get("translate_code", translate_code_default))

    if library_type == "barcoded" and translate_code:
        fail(
            'For library_type="barcoded", translate_code must be false. '
            "Barcodes are checked as DNA but not translated."
        )

    if library_type == "flycoded" and not translate_code:
        fail(
            'For library_type="flycoded", translate_code must be true. '
            "Flycodes are checked as DNA and translated into peptide sequence."
        )

    required_code_length = config.get(
        "required_code_length",
        required_code_length_default,
    )

    if required_code_length is not None:
        required_code_length = int(required_code_length)
        if required_code_length < 1:
            fail("required_code_length must be a positive integer or null.")

    return {
        "library_type": library_type,
        "code_name": code_name,
        "barcode_regex": barcode_regex,
        "compiled_code_regex": compile_regex(barcode_regex, "barcode_regex"),
        "translate_code": translate_code,
        "required_code_length": required_code_length,
    }


def validate_dna_code(
    dna_code: str,
    line_number: int,
    library_settings: dict,
) -> str:
    dna_code = clean_dna(
        dna_code,
        f"DNA {library_settings['code_name']} on line {line_number}",
    )

    required_code_length = library_settings["required_code_length"]

    if required_code_length is not None and len(dna_code) != required_code_length:
        fail(
            f"DNA {library_settings['code_name']} on line {line_number} "
            f"has length {len(dna_code)} bp, expected {required_code_length} bp."
        )

    if library_settings["compiled_code_regex"].fullmatch(dna_code) is None:
        fail(
            f"DNA {library_settings['code_name']} on line {line_number} "
            "does not match barcode_regex."
        )

    return dna_code


def classify_variant_gene(
    variant_gene_dna: str,
    wt_gene_dna: str,
    wt_protein: str,
    allowed_positions: set[int],
    allowed_codons: set[str],
    position_offset: int,
    allow_synonymous_as_wt: bool,
) -> tuple[str | None, str]:
    """
    Return:
        (variant_code, reason)

    Accepted variant_code:
        - "wt"
        - "A512Q"

    Rejected:
        - returns (None, reason)
    """
    variant_gene_dna = clean_dna(variant_gene_dna, "variant gene DNA")

    if len(variant_gene_dna) != len(wt_gene_dna):
        return None, (
            f"gene_length_mismatch:"
            f"wt={len(wt_gene_dna)},variant={len(variant_gene_dna)}"
        )

    if len(variant_gene_dna) % 3 != 0:
        return None, "variant_gene_length_not_divisible_by_3"

    if variant_gene_dna == wt_gene_dna:
        return "wt", "accepted_exact_wt"

    variant_protein = translate_dna(
        variant_gene_dna,
        "variant gene DNA",
        allow_terminal_stop=False,
    )

    if variant_protein == wt_protein:
        if allow_synonymous_as_wt:
            return "wt", "accepted_synonymous_as_wt"
        return None, "synonymous_but_not_exact_wt"

    differing_codons = []

    for index in range(0, len(wt_gene_dna), 3):
        wt_codon = wt_gene_dna[index:index + 3]
        variant_codon = variant_gene_dna[index:index + 3]

        if wt_codon != variant_codon:
            aa_position = index // 3 + 1
            differing_codons.append((aa_position, wt_codon, variant_codon))

    if len(differing_codons) == 0:
        return "wt", "accepted_exact_wt"

    if len(differing_codons) > 1:
        return None, f"multiple_codon_mutations:{len(differing_codons)}"

    aa_position, wt_codon, variant_codon = differing_codons[0]
    reported_position = aa_position + position_offset

    if aa_position not in allowed_positions:
        return None, f"mutation_position_not_allowed:{aa_position}"

    if variant_codon not in allowed_codons:
        return None, f"substitution_codon_not_allowed:{variant_codon}"

    wt_aa = GENETIC_CODE[wt_codon]
    variant_aa = GENETIC_CODE[variant_codon]

    if variant_aa == "*":
        return None, "stop_codon_substitution"

    if wt_aa == variant_aa:
        if allow_synonymous_as_wt:
            return "wt", "accepted_synonymous_as_wt"
        return None, "synonymous_but_not_exact_wt"

    return (
        f"{wt_aa}{reported_position}{variant_aa}",
        "accepted_single_allowed_substitution",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Filter barcode/flycode-variant maps, optionally translate "
            "flycodes, assign variant codes, and remove redundant DNA-code "
            "assignments."
        )
    )

    parser.add_argument(
        "-p",
        "--parameters",
        required=True,
        type=Path,
        help="JSON parameter file",
    )

    args = parser.parse_args()
    config = load_json(args.parameters)

    input_file = Path(config.get("input_file", ""))
    if not input_file.exists():
        fail(f"Input file not found: {input_file}")

    library_settings = get_library_settings(config)

    wt_gene_dna = clean_dna(config.get("wt_gene_dna", ""), "wt_gene_dna")

    if len(wt_gene_dna) % 3 != 0:
        fail(f"wt_gene_dna length must be divisible by 3: {len(wt_gene_dna)} nt")

    allow_terminal_stop_in_wt = bool(config.get("allow_terminal_stop_in_wt", False))
    wt_protein = translate_dna(
        wt_gene_dna,
        "wt_gene_dna",
        allow_terminal_stop=allow_terminal_stop_in_wt,
    )

    delimiter = config.get("delimiter", "auto")
    has_header = bool(config.get("has_header", False))
    allow_synonymous_as_wt = bool(config.get("allow_synonymous_as_wt", False))
    position_offset = int(config.get("position_offset", 0))

    allowed_positions = normalize_positions(
        config.get("allowed_mutation_positions", [])
    )
    allowed_codons = normalize_codons(
        config.get("allowed_substitution_codons", [])
    )

    prefix = config.get("output_prefix", default_output_prefix(input_file))

    filtered_txt_gz = Path(
        config.get("filtered_txt_gz", f"{prefix}_filtered.txt.gz")
    )
    accepted_tsv = Path(
        config.get("peptide_variantcode_tsv", f"{prefix}_peptide_variantcode.tsv")
    )
    filtered_out_tsv = Path(
        config.get(
            "peptide_variantcode_filtered_out_tsv",
            f"{prefix}_peptide_variantcode_filteredOUT.tsv",
        )
    )

    for output_path in (filtered_txt_gz, accepted_tsv, filtered_out_tsv):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    initially_accepted = 0
    initially_filtered = 0
    code_pattern_filtered = 0
    barcode_to_variant_keys: dict[str, set[tuple[str, str]]] = defaultdict(set)
    accepted_records: list[dict] = []
    filtered_records: list[dict] = []

    with open_text(input_file, "r") as source:
        if has_header:
            next(source, None)

        for line_number, line in enumerate(source, start=2 if has_header else 1):
            if not line.strip():
                continue

            total_rows += 1

            try:
                dna_code, variant_sequence = parse_input_line(
                    line=line,
                    line_number=line_number,
                    delimiter=delimiter,
                )

                dna_code = validate_dna_code(
                    dna_code=dna_code,
                    line_number=line_number,
                    library_settings=library_settings,
                )

                variant_sequence = clean_dna(
                    variant_sequence,
                    f"variant sequence on line {line_number}",
                )

                if library_settings["translate_code"]:
                    peptide_code = translate_dna(
                        dna_code,
                        f"DNA {library_settings['code_name']} on line {line_number}",
                        allow_terminal_stop=False,
                    )
                else:
                    peptide_code = ""

                variant_code, reason = classify_variant_gene(
                    variant_gene_dna=variant_sequence,
                    wt_gene_dna=wt_gene_dna,
                    wt_protein=wt_protein,
                    allowed_positions=allowed_positions,
                    allowed_codons=allowed_codons,
                    position_offset=position_offset,
                    allow_synonymous_as_wt=allow_synonymous_as_wt,
                )

                if variant_code is None:
                    initially_filtered += 1
                    filtered_records.append(
                        {
                            "dna_code": dna_code,
                            "peptide_code": peptide_code,
                            "variant_sequence": variant_sequence,
                            "reason": reason,
                        }
                    )
                    continue

                initially_accepted += 1

                record = {
                    "dna_code": dna_code,
                    "peptide_code": peptide_code,
                    "variant_sequence": variant_sequence,
                    "variant_code": variant_code,
                    "reason": reason,
                }
                accepted_records.append(record)
                barcode_to_variant_keys[dna_code].add((variant_code, variant_sequence))

            except ValueError as error:
                initially_filtered += 1

                reason = str(error)
                if "does not match barcode_regex" in reason or "expected 33 bp" in reason:
                    code_pattern_filtered += 1

                try:
                    dna_code, variant_sequence = parse_input_line(
                        line=line,
                        line_number=line_number,
                        delimiter=delimiter,
                    )
                except ValueError:
                    dna_code = "NA"
                    variant_sequence = line.strip()

                filtered_records.append(
                    {
                        "dna_code": dna_code,
                        "peptide_code": "",
                        "variant_sequence": variant_sequence,
                        "reason": reason,
                    }
                )

    final_accepted_records = []
    redundancy_filtered_records = []

    for record in accepted_records:
        dna_code = record["dna_code"]
        variant_keys = barcode_to_variant_keys[dna_code]

        if len(variant_keys) > 1:
            redundancy_filtered_records.append(
                {
                    "dna_code": record["dna_code"],
                    "peptide_code": record["peptide_code"],
                    "variant_sequence": record["variant_sequence"],
                    "reason": (
                        "dna_code_assigned_to_multiple_variants:"
                        + str(len(variant_keys))
                    ),
                }
            )
        else:
            final_accepted_records.append(record)

    filtered_records.extend(redundancy_filtered_records)

    with open_text(filtered_txt_gz, "w") as handle:
        for record in final_accepted_records:
            handle.write(
                f"{record['dna_code']}\t{record['variant_sequence']}\n"
            )

    with accepted_tsv.open("w", encoding="utf-8") as handle:
        handle.write(
            "DNA_code\tPeptide_code\tVariant_code\tVariant_sequence\n"
        )

        for record in final_accepted_records:
            handle.write(
                f"{record['dna_code']}\t"
                f"{record['peptide_code']}\t"
                f"{record['variant_code']}\t"
                f"{record['variant_sequence']}\n"
            )

    with filtered_out_tsv.open("w", encoding="utf-8") as handle:
        handle.write(
            "DNA_code\tPeptide_code\tVariant_sequence\tReason\n"
        )

        for record in filtered_records:
            handle.write(
                f"{record['dna_code']}\t"
                f"{record['peptide_code']}\t"
                f"{record['variant_sequence']}\t"
                f"{record['reason']}\n"
            )

    unique_accepted_codes = len({r["dna_code"] for r in final_accepted_records})
    unique_accepted_variants = len(
        {(r["variant_code"], r["variant_sequence"]) for r in final_accepted_records}
    )

    print("Barcode/flycode variant filtering complete")
    print(f"Input file:                         {input_file}")
    print(f"Library type:                       {library_settings['library_type']}")
    print(f"DNA code regex:                     {library_settings['barcode_regex']}")
    print(f"Translate DNA code:                 {library_settings['translate_code']}")
    print(f"Required DNA code length:           {library_settings['required_code_length']}")
    print(f"WT gene length:                     {len(wt_gene_dna)} nt")
    print(f"WT protein length:                  {len(wt_protein)} aa")
    print(f"Total input rows:                   {total_rows}")
    print(f"Initially accepted rows:            {initially_accepted}")
    print(f"Initially filtered rows:            {initially_filtered}")
    print(f"Code pattern/length filtered rows:  {code_pattern_filtered}")
    print(f"Redundancy-filtered rows:           {len(redundancy_filtered_records)}")
    print(f"Final accepted rows:                {len(final_accepted_records)}")
    print(f"Unique accepted DNA codes:          {unique_accepted_codes}")
    print(f"Unique accepted variants:           {unique_accepted_variants}")
    print(f"Filtered txt.gz:                    {filtered_txt_gz}")
    print(f"Peptide/variant-code TSV:           {accepted_tsv}")
    print(f"Filtered-out TSV:                   {filtered_out_tsv}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
