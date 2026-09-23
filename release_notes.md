## RAGE-toolkit/Artic-nf-GUI v2.0.0

Updated to follow the **fieldbioinformatics v1.11.2** workflow, with **Clair3** replacing
Medaka/Longshot for variant calling. Ported from `RAGE-toolkit/rabv-artic-nf` v2.0.0,
adapted to this workflow's typed, `sampleId`-keyed channel/tuple architecture (rather than
rabv-artic-nf's shared-filesystem, string-path style).

### Highlights
- **Clair3 variant calling** per primer pool (pool 1, pool 2), plus calling on unmatched
  reads, then merged, filtered, normalised and turned into a consensus, mirroring
  fieldbioinformatics v1.11.2
- **Single-pass primer trimming**: `align_trim` now runs once per sample (RG-tagged for
  pool 1 / pool 2 / unmatched), replacing the old two-stage `ALIGN_TRIM_1` / `ALIGN_TRIM_2`
  custom scripts
- **Real Nextflow channels throughout**: every new/updated process takes and emits typed
  `tuple(sampleId, path(...))` channels joined on `sampleId`, consistent with the rest of
  this workflow - no reconstructed file paths or shared `currDir` side effects
- **Same multi-architecture Docker images** as rabv-artic-nf: `rage2025/artic-nf-arm64:v2.0`
  and `rage2025/artic-nf-amd64:v2.0`, chosen automatically from the host architecture

### New configurable parameters
| Parameter | Default | Purpose |
|---|---|---|
| `normalise` | 100 | Reads kept per amplicon |
| `primer_match_threshold` | 35 | Max distance (bp) between a read end and its primer site |
| `min_mapq` | 20 | Minimum mapping quality |
| `min_variant_quality` | 10 | Clair3 QUAL cutoff |
| `min_allele_frequency` | 0.6 | AF needed to enter the consensus |
| `min_mask_allele_frequency` | 0.1 | Below this AF, variant discarded |
| `min_frameshift_quality` | 50 | QUAL needed for frameshifting indels |
| `min_minor_allele_count` | 4 | Minimum alt-supporting reads |
| `model_path` | `/opt/miniforge/envs/clair3env/bin/models/r1041_e82_400bps_hac_v520` | Clair3 model directory |

`mask_depth` now also sets the VCF filter's minimum depth.

### Changes
- `main.nf`: the medaka/longshot variant-calling chain (`ALIGN_TRIM_1`/`ALIGN_TRIM_2` →
  `MEDAKA_1`/`MEDAKA_2` → `MEDAKA_SNP_1`/`MEDAKA_SNP_2` → `VCF_MERGE` → `LONGSHOT`) is
  replaced by `ALIGN_TRIM` → `SPLIT_UNMATCHED` → `CLAIR3`/`CLAIR3_2`/`CLAIR3_UNMATCHED` →
  `VCF_MERGE` → `VCF_FILTER` → `COMPRESS_AND_INDEX_VCF` → `MAKE_DEPTH_MASK` → `MASK` →
  `BCFTOOLS_NORM` → `BCFTOOLS_CONSENSUS`
- Per-sample `MUSCLE`/`CONCAT_FOR_MUSCLE` alignment step removed; `CONCAT` now combines the
  `FASTA_HEADER` consensus outputs directly (via `consensus_combiner.py`), matching
  fieldbioinformatics - global alignment across samples still happens once, in `MAFFT`
- `scripts/align_trim.py`, `make_depth_mask.py`, `mask.py`, `vcf_merge.py`, `vcf_filter.py`
  replaced with the fieldbioinformatics-style versions from `rabv-artic-nf`; added
  `consensus_combiner.py`
- `ALIGN_TRIM` calls the `align_trim` console script (bioconda `align_trim>=1.2.0`, already
  present in the shared v2.0 Docker image) instead of a bundled Python script
- `vcf_merge.py` now needs `primalbedtools`, which lives only in the image's `clair3env`
  conda environment - `VCF_MERGE` invokes `/opt/miniforge/envs/clair3env/bin/python3`
  directly rather than the system `python`
- Medaka, Longshot and MUSCLE modules moved to `modules/deprecated/`; the scripts they
  superseded moved to `scripts/deprecated/`
- `PLEX_DIRS` now calls `scripts/plex.py` (the parallelised, deduplicating
  fieldbioinformatics-style `guppyplex` rewrite from rabv-artic-nf) instead of the old
  `directory_plex.py`, and now also applies `--max-length ${params.seq_max_len}` (a filter
  the old script supported but the module never passed); `directory_plex.py` moved to
  `scripts/deprecated/`

### Fixes
- `nextflow.config`'s Docker image selection picked the wrong architecture image
  (`artic-nf-amd64` on Apple Silicon) when Nextflow's JVM itself was an x86_64 build
  running under Rosetta - `System.getProperty('os.arch')` reflects the JVM binary's own
  architecture, not the host's, and shelling out to `uname -m` doesn't help either since
  macOS propagates Rosetta translation to child processes by default. Detection now
  queries `sysctl hw.optional.arm64`, which reads the kernel's view of the physical
  hardware directly and is unaffected by the calling process's translation state
- `BCFTOOLS_CONSENSUS` failed with `could not load index` for the normalised VCF -
  `BCFTOOLS_NORM`'s `.tbi` output was declared but never joined into the consensus
  channel, so the index never got staged into that task's work dir. `bcftools consensus`
  (unlike `samtools`) does not regenerate a missing index itself. Fixed by joining
  `BCFTOOLS_NORM.out.normalised_tbi` in `main.nf` and adding a defensive
  `tabix`-if-missing fallback in `BCFTOOLS_CONSENSUS`

### Known issues
- Not yet run end-to-end against real data in this environment (no local `nextflow`
  binary, and the workflow requires the `rage2025/artic-nf-*:v2.0` container) - review
  before trusting in production
- The `conda` profile is not defined for the Clair3 modules (same limitation as
  rabv-artic-nf). Use `-profile docker`.
