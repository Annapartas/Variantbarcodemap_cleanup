# BC/FC variant filter

Filter DNA barcode–variant or flycode–variant maps and assign protein mutation labels.

## What it does

- Validates DNA barcodes (33 bp) or translates flycodes into peptide codes.
- Accepts an exact wild-type gene (`wt`) or a **single permitted codon substitution** at an allowed position.
- Rejects invalid codes, multiple substitutions, disallowed mutations, and ambiguous DNA-code-to-variant assignments.
- Writes a filtered two-column map, an accepted-variant table, and a rejection table with reasons.

## Requirements

Python **3.10+** (standard library only; no additional Python packages).

## Configure

Choose and edit one of the included files:

- `barcoded_variant_filter_params.json` — DNA barcode library.
- `flycoded_variant_filter_params.json` — DNA flycode library.

Set `input_file` to your two-column `.txt`, `.tsv`, or `.txt.gz` input (no header unless `has_header` is `true`), and adjust the `wt_gene_dna`, permitted mutation positions and codons, and output paths as needed. **Variant sequences must include the full gene, in the same reading frame and length as `wt_gene_dna`.**

Absolute input/output paths work on clusters. Relative paths are interpreted from the directory in which you run the pipeline; the supplied shell script changes into the repository directory.

## Run

**Python (local or interactive cluster session)**

```bash
python3 BC_FC_variant_filter.py --parameters barcoded_variant_filter_params.json
# or
python3 BC_FC_variant_filter.py --parameters flycoded_variant_filter_params.json
```

**SLURM cluster** — from the repository directory:

```bash
sbatch run_BC_FC_variant_filter.sh  # barcoded configuration
sbatch --export=ALL,PARAMS=flycoded_variant_filter_params.json run_BC_FC_variant_filter.sh
```

The job log is `BC_FC_filter_JOBID.log` in the directory where you submitted the job. To request notification emails (if supported by your cluster), append `--mail-user=you@university.edu` to `sbatch`. Set `CONDA_BASE` and optionally `BC_FILTER_ENV` if your cluster requires Conda activation; otherwise a working `python3` is sufficient. Adjust SLURM memory and time settings in the script to suit your data.

## Outputs

| File | Contents |
| --- | --- |
| `*_filtered.txt.gz` | Accepted `DNA_code` and `Variant_sequence` |
| `*_peptide_variantcode.tsv` (or configured filename) | Accepted codes, peptides (flycodes only), variant labels, sequences |
| `*_peptide_variantcode_filteredOUT.tsv` (or configured filename) | Rejected entries and reasons |

The actual filenames come from the selected JSON file. The completion summary and counts are printed to the terminal or SLURM log.

**Numbering:** `allowed_mutation_positions` uses 1-based positions in the **supplied WT gene**; `position_offset` changes only the displayed mutation label. If your upstream pipeline omitted the first six coding nucleotides, correct that upstream and use matching full-length reference and variant sequences here.
