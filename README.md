# BC/FC variant filter

Filter barcode–variant or flycode–variant maps and prepare **two-column Enrich2 input** files. Works with plain-text or gzipped input and can process multiple files in one run. Uses Python 3 standard library only.

## What it does

- Checks barcode or flycode sequence design (and translates flycodes).
- Accepts exact wild-type DNA or **one** allowed codon substitution at an allowed amino-acid position. Every other nucleotide must match the WT reference. Synonymous changes are rejected by default.
- Excludes ambiguous DNA codes linked to **multiple distinct variant DNA sequences**, including initially rejected variants. The reason is recorded **only in the filtered-out report**, never in the Enrich2 map.
- Processes each input separately and generates per-file and combined QC statistics.

## Input and configuration

The input is a two-column file **without a header** (unless `has_header` is set):

```text
DNA_code    Variant_sequence
```

Edit one of the supplied configurations:
- `barcoded_variant_filter_params.json`
- `flycoded_variant_filter_params.json`

Specify files or wildcard patterns (absolute paths work on clusters):

```json
"input_files": ["/path/to/data/*_unique.txt.gz"],
"output_directory": "/path/to/results"
```

Set `wt_gene_dna`, `allowed_mutation_positions`, and `allowed_substitution_codons` for your experiment. Mutation positions are **1-based** and relative to the full WT gene supplied. The supplied configurations are examples based on MsbA.

## Run

```bash
python3 BC_FC_variant_filter.py --parameters barcoded_variant_filter_params.json
```

For SLURM, submit **from the repository directory**:

```bash
sbatch run_BC_FC_variant_filter.sh
# or flycoded input:
sbatch --export=ALL,PARAMS=flycoded_variant_filter_params.json run_BC_FC_variant_filter.sh
```

Check `BC_FC_filter_JOBID.log` and `sacct -j JOBID`. Email events (`BEGIN,END,FAIL`) can be enabled with `sbatch --mail-user=you@example.org ...`, if your cluster supports mail. The shell script optionally activates Conda when `CONDA_BASE` is set, but Conda is not required for the pipeline itself.

## Outputs

For each `sample_unique.txt.gz`:

| File | Contents |
| --- | --- |
| `sample_unique_filtered.txt.gz` | **Enrich2 input:** exactly two tab-separated DNA columns, no header, no extra annotations |
| `sample_unique_protein_variantcode.tsv` | Accepted barcode/flycode, peptide (for flycodes), variant name, variant DNA |
| `sample_unique_protein_variantcode_filteredOUT.tsv` | Rejected records and explicit rejection reasons, including ambiguous codes |
| `sample_unique_filter_summary.json` | Per-file counts and rejection reasons |

`filter_summary.tsv` aggregates results for all inputs. Outputs are never silently overwritten: move/delete existing results to rerun. Inputs with identical basenames from different directories must be renamed to avoid output collisions.

## Tests

```bash
python3 -m unittest discover -s tests -v
```
