#!/usr/bin/env python3
"""Filter barcode/flycode variant maps, writing strict two-column Enrich2 inputs.

Usage: python BC_FC_variant_filter.py --parameters barcoded_variant_filter_params.json
"""

from __future__ import annotations

import argparse
import glob
import os
import sqlite3
import tempfile
from collections import Counter
import gzip
import json
import re
import sys
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


def resolve_inputs(config: dict) -> list[Path]:
    """Expand explicit input paths and globs, keeping each file only once."""
    patterns = config.get("input_files", config.get("input_file"))
    if isinstance(patterns, str):
        patterns = [patterns]
    if not isinstance(patterns, list) or not patterns or any(
        not isinstance(p, str) or not p.strip() for p in patterns
    ):
        fail("input_files must be a nonempty list of paths or wildcard patterns.")
    paths = []
    seen = set()
    for pattern in patterns:
        matches = sorted(glob.glob(os.path.expanduser(pattern)))
        if not matches:
            fail(f"No input files match: {pattern}")
        for match in matches:
            path = Path(match).resolve()
            if not path.is_file():
                fail(f"Input is not a regular file: {path}")
            if path not in seen:
                paths.append(path)
                seen.add(path)
    return paths


def output_stem(path: Path) -> str:
    # Preserve the 'unique' part, so X_unique.txt.gz -> X_unique_filtered.txt.gz.
    return default_output_prefix(path)


def process_file(input_file: Path, output_dir: Path, settings: dict) -> dict:
    """Filter one input using a temporary SQLite database, then publish outputs."""
    stem = output_stem(input_file)
    final_map = output_dir / f"{stem}_filtered.txt.gz"
    final_accepted = output_dir / f"{stem}_protein_variantcode.tsv"
    final_rejected = output_dir / f"{stem}_protein_variantcode_filteredOUT.tsv"
    final_stats = output_dir / f"{stem}_filter_summary.json"
    finals = [final_map, final_accepted, final_rejected, final_stats]
    # Never accidentally overwrite existing results (including previous runs).
    if any(path.exists() for path in finals):
        fail(f"Output already exists for {input_file.name}. Move/delete old results before rerunning: "
             + ", ".join(str(p) for p in finals if p.exists()))

    with tempfile.TemporaryDirectory(prefix=f".{stem}_", dir=output_dir) as work:
        work = Path(work)
        db = sqlite3.connect(work / "records.sqlite")
        try:
            db.execute("CREATE TABLE rows (rowid INTEGER PRIMARY KEY, code TEXT, peptide TEXT, "
                       "variant TEXT, sequence TEXT, reason TEXT, accepted INTEGER)")
            db.execute("CREATE TABLE assignments (code TEXT, sequence TEXT, PRIMARY KEY(code, sequence))")
            counts = Counter()
            rejected_reasons = Counter()
            with open_text(input_file, "r") as source:
                if settings["has_header"]:
                    next(source, None)
                for line_no, line in enumerate(source, start=2 if settings["has_header"] else 1):
                    if not line.strip():
                        continue
                    counts["total_rows"] += 1
                    code, variant, peptide = "NA", line.strip(), ""
                    try:
                        code, variant = parse_input_line(line, line_no, settings["delimiter"])
                        code = validate_dna_code(code, line_no, settings["library_settings"])
                        variant = clean_dna(variant, f"variant sequence on line {line_no}")
                        # Important: record *all* well-formed barcode/variant associations,
                        # even those rejected for mutation reasons. One code associated with
                        # multiple DNA sequences must not reach the Enrich2 input.
                        db.execute("INSERT OR IGNORE INTO assignments VALUES (?,?)", (code, variant))
                        if settings["library_settings"]["translate_code"]:
                            peptide = translate_dna(code, f"flycode on line {line_no}")
                        variant_label, reason = classify_variant_gene(
                            variant_gene_dna=variant,
                            wt_gene_dna=settings["wt_gene_dna"],
                            wt_protein=settings["wt_protein"],
                            allowed_positions=settings["allowed_positions"],
                            allowed_codons=settings["allowed_codons"],
                            position_offset=settings["position_offset"],
                            allow_synonymous_as_wt=settings["allow_synonymous_as_wt"],
                        )
                    except ValueError as error:
                        variant_label, reason = None, str(error)
                        if "does not match barcode_regex" in reason or "expected 33 bp" in reason:
                            counts["code_pattern_filtered_rows"] += 1
                    accepted = int(variant_label is not None)
                    counts["initially_accepted_rows" if accepted else "initially_filtered_rows"] += 1
                    if not accepted:
                        rejected_reasons[reason] += 1
                    db.execute("INSERT INTO rows (code,peptide,variant,sequence,reason,accepted) "
                               "VALUES (?,?,?,?,?,?)",
                               (code, peptide, variant_label or "", variant, reason, accepted))
                    if counts["total_rows"] % 10000 == 0:
                        db.commit()
            db.commit()
            # Count distinct variant DNA sequences per valid code; no row is annotated in the
            # Enrich2 map. Exclusions and reasons appear only in the rejection TSV.
            db.execute("CREATE TABLE ambiguous AS SELECT code, COUNT(*) AS n FROM assignments "
                       "GROUP BY code HAVING COUNT(*) > 1")
            db.execute("CREATE UNIQUE INDEX ambiguous_code ON ambiguous(code)")
            counts["ambiguous_dna_codes"] = db.execute("SELECT COUNT(*) FROM ambiguous").fetchone()[0]
            tmp_map, tmp_accepted, tmp_rejected, tmp_stats = [work / f.name for f in finals]
            with gzip.open(tmp_map, "wt", encoding="utf-8") as valid_map, \
                 tmp_accepted.open("w", encoding="utf-8") as accepted_tsv, \
                 tmp_rejected.open("w", encoding="utf-8") as rejected_tsv:
                accepted_tsv.write("DNA_code\tPeptide_code\tVariant_code\tVariant_sequence\n")
                rejected_tsv.write("DNA_code\tPeptide_code\tVariant_sequence\tReason\n")
                for code, peptide, label, sequence, reason, was_accepted, ambiguity in db.execute(
                    "SELECT r.code,r.peptide,r.variant,r.sequence,r.reason,r.accepted,a.n "
                    "FROM rows r LEFT JOIN ambiguous a ON r.code=a.code ORDER BY r.rowid"
                ):
                    if ambiguity is not None:
                        # Rejected initially invalid rows also get an ambiguity reason for transparency.
                        counts["redundancy_filtered_rows"] += int(was_accepted)
                        reason = f"dna_code_assigned_to_multiple_variant_sequences:{ambiguity}" + (
                            f";original_reason:{reason}" if not was_accepted else "")
                        rejected_reasons["ambiguous_dna_code"] += 1
                    if was_accepted and ambiguity is None:
                        valid_map.write(f"{code}\t{sequence}\n")
                        accepted_tsv.write(f"{code}\t{peptide}\t{label}\t{sequence}\n")
                        counts["final_accepted_rows"] += 1
                    else:
                        rejected_tsv.write(f"{code}\t{peptide}\t{sequence}\t{reason}\n")
                        counts["final_rejected_rows"] += 1
            counts["unique_accepted_codes"] = db.execute(
                "SELECT COUNT(DISTINCT r.code) FROM rows r LEFT JOIN ambiguous a ON r.code=a.code "
                "WHERE r.accepted=1 AND a.code IS NULL").fetchone()[0]
            counts["unique_accepted_variants"] = db.execute(
                "SELECT COUNT(*) FROM (SELECT DISTINCT r.variant,r.sequence FROM rows r "
                "LEFT JOIN ambiguous a ON r.code=a.code WHERE r.accepted=1 AND a.code IS NULL)"
            ).fetchone()[0]
            summary = {"input_file": str(input_file), "library_type": settings["library_settings"]["library_type"],
                       **{key: counts[key] for key in (
                           "total_rows", "initially_accepted_rows", "initially_filtered_rows",
                           "code_pattern_filtered_rows", "ambiguous_dna_codes",
                           "redundancy_filtered_rows", "final_accepted_rows", "final_rejected_rows",
                           "unique_accepted_codes", "unique_accepted_variants")},
                       "rejection_reasons": dict(rejected_reasons),
                       "filtered_map": str(final_map), "accepted_tsv": str(final_accepted),
                       "filtered_out_tsv": str(final_rejected)}
            tmp_stats.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
            # Only completed files are placed in final location.
            for temp_path, final_path in zip((tmp_map,tmp_accepted,tmp_rejected,tmp_stats), finals):
                os.replace(temp_path, final_path)
            return summary
        finally:
            db.close()


def main() -> int:
    parser = argparse.ArgumentParser(description="Filter multiple barcode/flycode variant maps.")
    parser.add_argument("-p", "--parameters", required=True, type=Path, help="JSON parameter file")
    args = parser.parse_args()
    config = load_json(args.parameters)
    if "output_prefix" in config or any(k in config for k in (
        "filtered_txt_gz", "protein_variantcode_tsv", "protein_variantcode_filtered_out_tsv"
    )):
        fail("Remove legacy output_prefix / per-file output names from the JSON; "
             "output filenames are now derived from each input filename.")
    inputs = resolve_inputs(config)
    output_dir = Path(config.get("output_directory", "results")).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stems = [output_stem(p) for p in inputs]
    if len(stems) != len(set(stems)):
        fail("Two input files produce the same output basename; rename inputs to prevent collisions.")
    # Fail before any processing if a previous run would be overwritten.
    for stem in stems:
        expected = [
            output_dir / f"{stem}_filtered.txt.gz",
            output_dir / f"{stem}_protein_variantcode.tsv",
            output_dir / f"{stem}_protein_variantcode_filteredOUT.tsv",
            output_dir / f"{stem}_filter_summary.json",
        ]
        if any(p.exists() for p in expected):
            fail(f"Existing outputs for {stem}; move/delete them before rerunning.")
    if (output_dir / "filter_summary.tsv").exists():
        fail("Combined filter_summary.tsv already exists; move/delete it before rerunning.")
    wt = clean_dna(config.get("wt_gene_dna", ""), "wt_gene_dna")
    if len(wt) % 3:
        fail("WT DNA length must be divisible by three.")
    protein = translate_dna(wt, "wt_gene_dna", config.get("allow_terminal_stop_in_wt", False))
    positions = normalize_positions(config.get("allowed_mutation_positions", []))
    if max(positions) > len(protein):
        fail(f"Allowed mutation position {max(positions)} exceeds WT protein length {len(protein)}.")
    settings = {
        "library_settings": get_library_settings(config),
        "wt_gene_dna": wt, "wt_protein": protein,
        "allowed_positions": positions,
        "allowed_codons": normalize_codons(config.get("allowed_substitution_codons", [])),
        "delimiter": config.get("delimiter", "auto"),
        "has_header": bool(config.get("has_header", False)),
        "allow_synonymous_as_wt": bool(config.get("allow_synonymous_as_wt", False)),
        "position_offset": int(config.get("position_offset", 0)),
    }
    summaries = []
    for index, path in enumerate(inputs, 1):
        print(f"[{index}/{len(inputs)}] Filtering {path}", flush=True)
        summary = process_file(path, output_dir, settings)
        summaries.append(summary)
        print(f"  Accepted: {summary['final_accepted_rows']}; rejected: "
              f"{summary['final_rejected_rows']}; ambiguous codes: {summary['ambiguous_dna_codes']}", flush=True)
    overall = output_dir / "filter_summary.tsv"
    if overall.exists():
        fail(f"Combined summary already exists: {overall}. Move/delete it before rerunning.")
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output_dir,
                                     prefix=".summary_", delete=False) as handle:
        tmp_summary = Path(handle.name)
        cols = ["input_file", "library_type", "total_rows", "initially_accepted_rows",
                "initially_filtered_rows", "code_pattern_filtered_rows", "ambiguous_dna_codes",
                "redundancy_filtered_rows", "final_accepted_rows", "final_rejected_rows",
                "unique_accepted_codes", "unique_accepted_variants"]
        handle.write("\t".join(cols) + "\n")
        for summary in summaries:
            handle.write("\t".join(str(summary[c]) for c in cols) + "\n")
    os.replace(tmp_summary, overall)
    print(f"Completed {len(summaries)} files. Summary: {overall}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, sqlite3.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
