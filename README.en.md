# Exploratory association networks for small-sample amplicon data

[中文](README.md) | **English**

**Current development version: `0.1.0`.** The planned first open-source tag is `v0.1.0`, marked as a pre-release; no tag or release has been published yet. See [CHANGELOG.md](CHANGELOG.md).

This free, local tool constructs exploratory Spearman association networks from ASV, OTU, or other feature-abundance tables. It focuses on small-sample amplicon studies, but also accepts other data that meet the input and statistical requirements.

> Statistical associations generate hypotheses; they do not establish direct microbial interactions or causality. Small samples, high dimensionality, compositional effects, and study design can affect results. This software is under development and does not guarantee a particular false-positive rate.

## Project structure

```text
network_tool/
├── cli.py / gui.py / network_core.py      # CLI, GUI, and shared analysis core
├── app_version.py / CHANGELOG.md          # Software version and changes
├── requirements.txt / environment.yml    # Python and Conda dependencies
├── README.md / README.en.md               # Chinese and English documentation
├── LICENSE / CONTRIBUTING.md              # License and contribution guide
├── .gitignore / .gitattributes           # Git ignore rules and data line endings
├── .github/workflows/tests.yml            # GitHub Actions configuration
├── docs/                                  # Provenance, sample list, notices, release checklist
│   ├── prjna1137280_t20_samples.tsv
│   ├── DATA_PROVENANCE.md
│   ├── THIRD_PARTY_NOTICES.md
│   └── RELEASE_CHECKLIST.md
├── examples/                              # Synthetic example inputs
├── tests/                                 # Automated tests
└── test_data/                             # Real-data case inputs and provenance
```

`test_data/` contains three derived inputs for the real-data case. Raw FASTQ files and local analysis outputs are not included. Inspect the actual file list before committing.

## Installation and quick start

Python 3.10 or newer, NumPy, and pandas are required. The GUI also requires Tkinter and a graphical display. If an existing environment meets these requirements, use it directly. To create a dedicated Conda environment from the project directory:

```bash
conda env create -f environment.yml
conda activate microbial-network
```

Alternatively install dependencies with `python -m pip install -r requirements.txt` in another Python environment. On Windows, `py -3` may be used if `python` is unavailable. The supplied examples are synthetic and do not represent real microbial communities:

```bash
python cli.py --counts examples/abundance.tsv --metadata examples/metadata.tsv --taxonomy-path examples/taxonomy.tsv --output demo_output
```

By default, samples are analyzed by group without p/q values. The output directory must not already exist; choose a new name for each run. Start the GUI with `python gui.py`. The interface and CLI messages are currently in Chinese. WSL 2 can run the GUI with WSLg; Linux requires a graphical session. On a headless server use the CLI; Conda is optional. See [Microsoft's WSL GUI guide](https://learn.microsoft.com/en-us/windows/wsl/tutorials/gui-apps).

Run `python cli.py --version` to inspect the version. Newly generated `report.json` files also include `software_version`; previous reports are not modified.

## Input

- **Abundance table:** features as rows and samples as columns; the first column contains unique feature IDs and the header contains unique sample IDs. TSV and CSV are supported, including a QIIME 2 BIOM-export `# Constructed from biom file` comment. Values must be finite and nonnegative, with a positive total for each sample.
- **Metadata (optional):** TSV/CSV with `sample-id` and `group` columns by default; column names are configurable. Sample IDs must match the abundance table one-to-one. Without metadata, all samples form one group.
- **Taxonomy (optional):** a QIIME 2 TSV export with unique `Feature ID` and `Taxon` columns and an optional `Confidence` column. Unmatched annotations are blank.

Each sample is divided by its total abundance. Users must assess whether this is appropriate for inputs already rarefied, transformed, or otherwise processed. Each analyzed group must contain at least four samples.

## GUI and CLI

In the GUI, select the abundance table, optional metadata and taxonomy, parameters, and output directory, then start the analysis. A new output directory is the default. An existing directory may be selected, but conflicting report or group paths are never overwritten. A failed run may leave a partial directory; inspect it and choose a new output path before retrying.

```bash
# Analyze groups without p/q filtering
python cli.py --counts abundance.tsv --metadata metadata.tsv --taxonomy-path taxonomy.tsv --output grouped_results

# Pool all samples; p/q are not calculated by default
python cli.py --counts abundance.tsv --metadata metadata.tsv --output pooled_results --use-all-samples yes

# Filter on raw p when permutation assumptions are appropriate; q is also reported
python cli.py --counts abundance.tsv --metadata metadata.tsv --output p_results --independent --p-max 0.05
```

Add `--q-max 0.05` to filter on BH-adjusted q, alone or together with a p threshold. Threshold comparisons are strict and require `--independent`. The user must assess the permutation assumptions from the study design; the software cannot establish them automatically. The GUI pooling option is unchecked by default. CLI pooling accepts `--use-all-samples y` or `yes`; `n` or `no` disables it.

## Parameters and computation

| CLI parameter | Default | Meaning |
|---|---:|---|
| `--sample-column` / `--group-column` | `sample-id` / `group` | Metadata sample and group columns |
| `--use-all-samples` | `no` | Pool all samples with `yes` or `y` |
| `--min-prevalence` | 0.5 | Minimum fraction of nonzero samples within a group |
| `--min-abundance` | 0.001 | Minimum within-group mean relative abundance (0.1%) |
| `--top-n` | 50 | Select the top N filtered features by mean relative abundance; minimum 3 |
| `--max-features` | 200 | Configured upper limit for N |
| `--pseudocount` | 0.5 | CPM pseudocount before CLR; currently configurable via CLI only |
| `--rho-min` | 0.7 | Minimum candidate-edge `|ρ|`, inclusive |
| `--min-stability` | 0.8 | Minimum reselected recovery frequency, inclusive |
| `--p-max` / `--q-max` | disabled | Optional strict raw p / BH q filters, separately or together |
| `--permutations` | 9999 | Random permutations when a group has more than nine samples |
| `--repeats` | 36 | Subsample repetitions; 0 disables subsampling |
| `--seed` | 2026 | Random seed; permutation and subsampling streams are independent |
| `--independent` | disabled | User confirms applicable permutation assumptions and enables p/q calculation |
| `--existing-output` | disabled | Use an existing directory without overwriting conflicting results |

Features are filtered by prevalence and mean relative abundance, then ranked by within-group mean relative abundance. Ties are broken by feature ID. Fewer than N eligible features are all retained; constant features after CLR are removed. CPM uses the full sample total. Selected features receive a pseudocount and CLR transformation before Spearman rank correlation; the CLR reference depends on the selected feature set.

When permutation assumptions are confirmed, groups of 4–9 samples use all permutations for two-sided exact p values. Larger groups use random permutations with `(b+1)/(B+1)`, followed by BH adjustment across all within-group pairs. Small random p values may be sensitive to the permutation budget; interpretation also depends on the null hypothesis, filtering, and dependence among tests. Leaving p/q thresholds blank disables filtering on them. If assumptions are confirmed they are still calculated; otherwise p/q fields are blank.

Stability includes leave-one-out and subsample assessments. Each subset contains `max(3, floor(0.8 × sample count))` samples. All combinations are used when their number does not exceed the requested repetitions; otherwise subsets are randomly sampled and may repeat. Reselected stability repeats filtering, top-N selection, and CLR for each subset. Missing or uncomputable edges count as not recovered in the denominator of all subsets. Frequency is not a probability that an edge is true. With `--repeats 0`, subsample recovery is unavailable and the candidate table is empty.

**Grouping and permutation limitations:** a pooled network measures associations across all samples; between-group differences can drive correlations and do not demonstrate stability within each treatment. Independent samples do not automatically satisfy the required exchangeability. Paired or repeated samples, batches, time, and treatment structures must be assessed. Stratified, paired, and restricted permutations are not implemented; pooled subsampling is not stratified either.

## Output

Each analyzed group has its own directory; pooled output uses `all_samples`. `report.json` records the software version, input hashes, parameters, groups, original group sizes, and NumPy/pandas versions. During development, retain the corresponding source code or commit as well.

| File | Contents |
|---|---|
| `edges.csv` | All pairs, including `source`, `target`, Spearman `rho`, `p_value`, `q_value`, leave-one-out and subsampling metrics |
| `candidate_edges.csv` | Edges satisfying `|ρ|`, reselected stability, and any enabled p/q thresholds; not confirmed interactions |
| `nodes.csv` | Final feature IDs and taxonomy |
| `selected_taxonomy.csv` | Final features, annotations, confidence, and within-group mean relative abundance |

`rho` may be used as a signed edge weight. For visualization software requiring nonnegative weights, use `|rho|` and retain the sign separately. CSV tables can be imported into tools such as Gephi. An empty candidate table is a valid result; thresholds are not automatically relaxed.

## Real-data case: PRJNA1137280

The project used public data for functional and methodological exploration: Hoefle D, Sperber S, Marzban N, Antonietti M, Hoffmann T, Wicaksono WA, Abdelfattah A, Berg G. 2024. *Artificial Humic Acid Diminishes the Effect of Drought on the Soil Microbiome*. Journal of Sustainable Agriculture and Environment 3(4):e70034. [doi:10.1002/sae2.70034](https://doi.org/10.1002/sae2.70034). Public accession: [PRJNA1137280](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1137280).

We selected 36 samples labeled `t20` in the authors' metadata, with nine per group. See the [sample-to-SRR accession manifest](docs/prjna1137280_t20_samples.tsv).

| Treatment | Sample IDs | Count |
|---|---|---:|
| Humic acid + drought | `S_5A1`–`S_5A9` | 9 |
| Humic acid + normal watering | `S_5B1`–`S_5B9` | 9 |
| No humic acid + drought | `S_5C1`–`S_5C9` | 9 |
| No humic acid + normal watering | `S_5D1`–`S_5D9` | 9 |

The publication reports 515F/806R amplification of the 16S rRNA V4 region, Illumina NovaSeq 6000 2 × 150 bp sequencing, and single-end processing in QIIME 2. **This project likewise uses a single-end strategy: R1 was denoised with DADA2 without paired-end merging.** The publication's main text does not specify R1 or R2, so we do not claim that its direction matches ours. The 130 bp truncation and settings below are specific to this project and do not imply an exact reproduction of the publication.

| Step | Settings actually used in this project |
|---|---|
| Primer trimming | QIIME 2 2025.7 / Cutadapt 5.1; anchored 5′ forward `^GTGYCAGCMGCCGCGGTAA` and reverse `^GGACTACNVGGGTWTCTAAT`; error rate 0.1, discard untrimmed reads, 4 cores |
| DADA2 | `denoise-single`, R1; `trim-left=0`, `trunc-len=130`, `max-ee=2`, `trunc-q=2`, independent denoising, `consensus` chimera removal, 16 threads; DADA2 1.30.0 |
| Error model | `learn_min_reads=1,000,000`; the log reports 1,814,611 reads and 235,899,430 bases from two samples actually used |
| Taxonomy | QIIME 2 `classify-sklearn`, SILVA 138.2 SSURef NR99 full-length uniform classifier; confidence 0.7, orientation `same`, 4 processes |
| Network input | Exclude chloroplast and mitochondrial annotations: 1,026 of 45,721 ASVs removed, leaving 44,695; no rarefaction; unclassified taxa and archaea retained |

These degenerate primer sequences are what this project used for trimming, not a claim about the original PCR sequences without further evidence. The publication reports BLAST classification and negative-control-based contaminant identification; this case used `classify-sklearn` and did not complete the same contamination analysis. Taxonomic resolution of 130 bp single-end reads is limited.

`test_data/` contains the derived abundance table, metadata and taxonomy for this case; raw FASTQ files are not included. `examples/` contains synthetic data. See [DATA_PROVENANCE.md](docs/DATA_PROVENANCE.md) and [test_data/README.md](test_data/README.md) for provenance, processing and licensing.

## References and license

- Bolyen E, et al. 2019. QIIME 2. *Nature Biotechnology* 37:852–857. [doi:10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9).
- Callahan BJ, et al. 2016. DADA2. *Nature Methods* 13:581–583. [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869).
- Quast C, et al. 2013. The SILVA rRNA gene database project. *Nucleic Acids Research* 41:D590–D596. [doi:10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219). See the [SILVA license](https://www.arb-silva.de/silva-license-information/).

Developed and maintained by **Xianyu Yin**, with development assistance from **Codex (OpenAI)**. Code is released under the [MIT License](LICENSE). See [THIRD_PARTY_NOTICES.md](docs/THIRD_PARTY_NOTICES.md) for third-party dependencies and data sources, and [CONTRIBUTING.md](CONTRIBUTING.md) for contributions. Run `python -m unittest discover -s tests -v` from the project directory.
