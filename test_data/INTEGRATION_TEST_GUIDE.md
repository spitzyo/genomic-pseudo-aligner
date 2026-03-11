# Integration Test Guide

## Biological Scenario

This test dataset simulates a **metagenomic profiling experiment** for a soil
bacterial community. You have sequenced DNA extracted directly from soil and
want to determine which bacterial species are present and in what proportions.
The three species in the reference database are:

| Genome identifier        | Organism                   | Length   | GC content |
|--------------------------|----------------------------|----------|------------|
| `ecoli_K12`              | *Escherichia coli* K-12    | 10,000 bp | 50.8 %     |
| `salmonella_LT2`         | *Salmonella enterica* LT2  | 10,200 bp | 51.4 %     |
| `bacillus_subtilis_168`  | *Bacillus subtilis* 168    |  9,800 bp | 44.3 %     |

The sequences are realistic in size and GC content (random but seeded, so
every run is fully reproducible).

---

## Read categories in `reads.fastq`

| Read IDs              | Category | Source                                   | Expected result          |
|-----------------------|----------|------------------------------------------|--------------------------|
| `READ_0001–READ_0080` | A        | Exact substrings of `ecoli_K12`          | **unique** → ecoli_K12   |
| `READ_0081–READ_0160` | B        | Exact substrings of `salmonella_LT2`     | **unique** → salmonella_LT2 |
| `READ_0161–READ_0200` | C        | Exact substrings of `bacillus_subtilis_168` | **unique** → bacillus_subtilis_168 |
| `READ_0201–READ_0230` | D        | Completely random sequences              | **unmapped**             |
| `READ_0231–READ_0245` | E        | `ecoli_K12` reads with a G→A SNP at position 2000 | **unique** → ecoli_K12, variant detected |

All reads are 150 bp (Illumina HiSeq 2500 format), with realistic Phred33
quality scores (Q25–Q40 with a 3′-end quality drop). Total: **245 reads**.

---

## Prerequisites

```bash
# From the repository root
pip install pytest   # only dependency
```

All `python main.py` commands below must be run from the **repository root**
(the directory that contains `main.py`).

---

## Step 0 — Verify the input files parse correctly

```bash
python - << 'EOF'
import sys; sys.path.insert(0, '.')
from helper_functions import import_fasta, import_fastq

refs  = import_fasta('test_data/reference_genomes.fa')
reads = import_fastq('test_data/reads.fastq')
for r in refs:
    print(f"  genome  : {r.identifier:30s}  {r.total_bases:,} bp")
print(f"  reads   : {len(reads)} reads parsed, first = {reads[0].identifier}")
EOF
```

**Expected output:**
```
  genome  : ecoli_K12                        10,000 bp
  genome  : salmonella_LT2                   10,200 bp
  genome  : bacillus_subtilis_168             9,800 bp
  reads   : 245 reads parsed, first = READ_0001|ecoli_K12_pos00200
```

---

## Step 1 — Build the k-mer reference database

Index all three genomes with k = 31 (standard for 150 bp bacterial reads)
and save to a compressed `.kdb` file.

```bash
python main.py \
  -t reference \
  -g test_data/reference_genomes.fa \
  -r test_data/bacteria.kdb \
  -k 31
```

No output is printed on success. Verify the file was created:

```bash
ls -lh test_data/bacteria.kdb
```

**Expected:** a `bacteria.kdb` file of roughly **360–380 KB**.

> **Why k = 31?**  A 31-mer has 4^31 ≈ 4.6 × 10^18 possible values. With
> genomes of ~10 kbp each, the probability of a k-mer appearing in more than
> one randomly generated genome is essentially zero. This guarantees that every
> k-mer in our test dataset is *species-specific*, giving clean, unambiguous
> mapping results.

---

## Step 2 — Inspect the reference database

```bash
python main.py \
  -t dumpref \
  -r test_data/bacteria.kdb
```

The full output includes every k-mer sequence and its genomic positions.  
Pipe through `python -c "…"` to see only the summary:

```bash
python main.py -t dumpref -r test_data/bacteria.kdb 2>&1 \
  | python -c "
import sys, json
data = json.load(sys.stdin)
print(json.dumps({'Summary': data['Summary']}, indent=4))
print('Total k-mers:', len(data['Kmers']))
"
```

**Expected output:**
```json
{
    "Summary": {
        "ecoli_K12": {
            "total_bases": 10000,
            "unique_kmers": 9970,
            "multi_mapping_kmers": 0,
            "soft_masked_bases": 0,
            "soft_masked_fraction": 0.0
        },
        "salmonella_LT2": {
            "total_bases": 10200,
            "unique_kmers": 10170,
            "multi_mapping_kmers": 0,
            "soft_masked_bases": 0,
            "soft_masked_fraction": 0.0
        },
        "bacillus_subtilis_168": {
            "total_bases": 9800,
            "unique_kmers": 9770,
            "multi_mapping_kmers": 0,
            "soft_masked_bases": 0,
            "soft_masked_fraction": 0.0
        }
    }
}
Total k-mers: 29910
```

**What to check:**

- `unique_kmers` equals `total_bases − k` (= total_bases − 31) for each genome.
  With 10,000 bp this gives 9,970 — the 30 k-mers "lost" at the edges are the
  ones that would extend beyond the sequence end.
- `multi_mapping_kmers` = 0 for every genome confirms that no 31-mer is shared
  between species — each read will map to exactly one reference, or not at all.
- `soft_masked_bases = 0` and `soft_masked_fraction = 0.0` for all three genomes
  because the synthetic reference sequences are fully uppercase. When working
  with a real genome downloaded from Ensembl or NCBI, soft-masked (repeat)
  regions will be in lowercase, and these fields will be non-zero.

---

## Step 3 — Align reads (save to file)

```bash
python main.py \
  -t align \
  -r test_data/bacteria.kdb \
  --reads test_data/reads.fastq \
  -a test_data/alignment.aln
```

No output is printed. The alignment is stored in `test_data/alignment.aln`
(a gzip-compressed pickle, ~35–40 KB).

---

## Step 4 — Inspect alignment results

```bash
python main.py \
  -t dumpalign \
  -a test_data/alignment.aln
```

**Expected output:**
```json
{
    "Statistics": {
        "unique_mapped_reads": 215,
        "ambiguous_mapped_reads": 0,
        "unmapped_reads": 30,
        "total_read_bases": 36750,
        "soft_masked_read_bases": 0,
        "soft_masked_read_fraction": 0.0
    },
    "Summary": {
        "ecoli_K12": {
            "unique_reads": 95,
            "ambiguous_reads": 0
        },
        "salmonella_LT2": {
            "unique_reads": 80,
            "ambiguous_reads": 0
        },
        "bacillus_subtilis_168": {
            "unique_reads": 40,
            "ambiguous_reads": 0
        }
    }
}
```

**How to verify this is correct:**

| Statistic | Expected | Explanation |
|-----------|----------|-------------|
| `unique_mapped_reads` | **215** | Cat A (80) + Cat B (80) + Cat C (40) + Cat E (15) |
| `ambiguous_mapped_reads` | **0** | No k-mer is shared between species (confirmed in Step 2) |
| `unmapped_reads` | **30** | All Cat D random reads |
| `total_read_bases` | **36,750** | 245 reads × 150 bp |
| `soft_masked_read_bases` | **0** | All reads in the test file are fully uppercase |
| `soft_masked_read_fraction` | **0.0** | No soft-masking in the test reads |
| `ecoli_K12 unique_reads` | **95** | Cat A (80 regular) + Cat E (15 SNP) |
| `salmonella_LT2 unique_reads` | **80** | Exactly Cat B |
| `bacillus_subtilis_168 unique_reads` | **40** | Exactly Cat C |

**Biological interpretation:** 88 % of reads mapped (215 / 245). The 12 %
unmapped reads correspond exactly to the synthetic random reads inserted to
simulate sequencing noise / contamination. The relative read counts
(95 : 80 : 40) reflect the sequencing depth assigned to each species, which
in a real experiment would indicate their relative abundance in the sample.

The `soft_masked_read_fraction` field tells you what fraction of all read bases
originated from soft-masked (lowercase, typically repetitive or low-complexity)
regions in the source FASTQ. A high value could indicate that many reads come
from repeat regions, which are harder to align uniquely.

---

## Step 5 — Soft-masking statistics (real-world reference)

The test genomes are fully uppercase, but real Ensembl/NCBI genomes contain
soft-masked regions (lowercase = repeat elements). To see the masking stats,
create a small example file with soft-masking:

```bash
cat > /tmp/repeat_genome.fa << 'EOF'
>chr1_with_repeats
ACGTacgtacgtACGTACGT
EOF

python main.py -t dumpref -g /tmp/repeat_genome.fa -k 4 2>&1 \
  | python -c "import sys,json; d=json.load(sys.stdin); print(json.dumps({'Summary':d['Summary']},indent=4))"
```

**Expected output:**
```json
{
    "Summary": {
        "chr1_with_repeats": {
            "total_bases": 20,
            "unique_kmers": 0,
            "multi_mapping_kmers": 4,
            "soft_masked_bases": 8,
            "soft_masked_fraction": 0.4
        }
    }
}
```

**Interpretation:**  
- `soft_masked_bases = 8`: the 8 lowercase characters (`acgtacgt`) were
  flagged as soft-masked by an upstream tool such as RepeatMasker.  
- `soft_masked_fraction = 0.4`: 40 % of this genome is in a repeat/low-
  complexity region.  
- Importantly, **all k-mers are still indexed and used for alignment** —
  the masking information is purely statistical. This is consistent with the
  pseudo-aligner convention: sequence content from soft-masked regions is
  retained but flagged so users can assess how much of their mapping comes
  from repetitive sequence.

To see soft-masking stats for reads, pass a FASTQ with lowercase bases:

```bash
cat > /tmp/masked_reads.fq << 'EOF'
@from_repeat_element
acgtACGT
+
IIIIIIII
@normal_read
ACGTACGT
+
IIIIIIII
EOF

python main.py -t dumpalign \
  -g /tmp/repeat_genome.fa -k 4 \
  --reads /tmp/masked_reads.fq 2>&1
```

**Expected output:**
```json
{
    "Statistics": {
        "unique_mapped_reads": 0,
        "ambiguous_mapped_reads": 0,
        "unmapped_reads": 2,
        "total_read_bases": 16,
        "soft_masked_read_bases": 4,
        "soft_masked_read_fraction": 0.25
    },
    "Summary": {}
}
```

`soft_masked_read_fraction = 0.25` (4 soft-masked bases out of 16 total)
shows that 25 % of the sequenced bases came from regions flagged as repetitive
in the read data — useful for downstream quality control.

---

## Step 6 — Coverage analysis

Run alignment and coverage in a single command (no `.aln` file needed):

```bash
python main.py \
  -t dumpalign \
  -r test_data/bacteria.kdb \
  --reads test_data/reads.fastq \
  --coverage
```

**Expected output** (alignment JSON first, then coverage JSON):
```json
{
    "Statistics": { ... },        ← same as Step 4
    "Summary": { ... }
}
{
    "Coverage": {
        "ecoli_K12": {
            "covered_bases_unique": 2536,
            "covered_bases_ambiguous": 0,
            "mean_coverage_unique": 0.3,
            "mean_coverage_ambiguous": 0.0
        },
        "salmonella_LT2": {
            "covered_bases_unique": 2480,
            "covered_bases_ambiguous": 0,
            "mean_coverage_unique": 0.2,
            "mean_coverage_ambiguous": 0.0
        },
        "bacillus_subtilis_168": {
            "covered_bases_unique": 1240,
            "covered_bases_ambiguous": 0,
            "mean_coverage_unique": 0.1,
            "mean_coverage_ambiguous": 0.0
        }
    }
}
```

**Understanding the coverage numbers:**

The pseudo-aligner tracks coverage at the k-mer *anchor* position — for each
read mapped to genome position G, it records k = 31 consecutive base positions
(G … G + 30). This is different from traditional full-read-extent coverage
mappers (like BWA or Bowtie) that would mark all 150 bases.

| Genome | Reads | Covered bases | Calculation |
|--------|-------|---------------|-------------|
| ecoli_K12 | 95 (80 regular + 15 SNP) | 2,536 | 80 non-overlapping × 31 bp + 87 bp from the 15 closely spaced SNP reads |
| salmonella_LT2 | 80 | 2,480 | 80 non-overlapping reads × 31 bp = 2,480 |
| bacillus_subtilis_168 | 40 | 1,240 | 40 non-overlapping reads × 31 bp = 1,240 |

To view per-base coverage for a single genome:

```bash
python main.py \
  -t dumpalign \
  -r test_data/bacteria.kdb \
  --reads test_data/reads.fastq \
  --coverage \
  --full-coverage \
  --genomes ecoli_K12
```

In the `Details.ecoli_K12.unique_cov` array you will see:

- positions 200–230: coverage = **1** (first read's 31-bp anchor window)
- positions 231–274: coverage = **0** (gap between reads; stride 75 > k 31)
- positions 275–305: coverage = **1** (second read's anchor window)
- etc.

---

## Step 6 — Variant detection

Detect the engineered G→A SNP at position 2000 in `ecoli_K12`.
Set `--min-variant-coverage 5` (needs ≥5 supporting reads) and
`--min-variant-quality 25` (Phred ≥ 25 at the variant base).

```bash
python main.py \
  -t dumpalign \
  -r test_data/bacteria.kdb \
  --reads test_data/reads.fastq \
  --detect-variants \
  --min-variant-coverage 5 \
  --min-variant-quality 25
```

**Expected variant output** (appended after the alignment JSON):
```json
{
    "Variants": {
        "ecoli_K12": {
            "2000": {
                "reference": "G",
                "alternate": "A",
                "quality_score": 35,
                "coverage": 15
            }
        }
    },
    "Statistics": {
        "ecoli_K12": {
            "total_variants": 1,
            "filtered_variants": 0
        }
    }
}
```

**What to verify:**

| Field | Value | Explanation |
|-------|-------|-------------|
| position key | `"2000"` | Exact position where the SNP was injected |
| `reference` | `"G"` | The base in the `ecoli_K12` reference at pos 2000 |
| `alternate` | `"A"` | The G→A transition introduced in Cat E reads |
| `quality_score` | `35` | Q35 was set for the SNP base in the generator |
| `coverage` | `15` | All 15 Cat E reads pass both quality and coverage filters |
| `total_variants` | `1` | Exactly one unique position shows a mismatch |
| `filtered_variants` | `0` | All detected variants meet the coverage threshold |

**Why no SNPs for salmonella or bacillus?**  
Cat B and Cat C reads are exact substrings of their respective reference
genomes — zero mismatches — so no variants are reported for those species.
This correctly mirrors real experiments where sequenced strains match
their reference perfectly.

**Increasing sensitivity:** lower `--min-variant-coverage` to `1` to see
every single-read mismatch (useful for very-low-coverage samples):

```bash
python main.py \
  -t dumpalign \
  -r test_data/bacteria.kdb \
  --reads test_data/reads.fastq \
  --detect-variants \
  --min-variant-coverage 1 \
  --min-variant-quality 20
```

---

## Step 7 — Quality filtering

Run with strict quality filters to see how they affect mapped read counts:

```bash
python main.py \
  -t dumpalign \
  -r test_data/bacteria.kdb \
  --reads test_data/reads.fastq \
  --min-read-quality 30 \
  --min-kmer-quality 25
```

With these thresholds a small number of reads whose mean quality falls
below Q30 (or whose individual k-mers fall below Q25) will be filtered out
and reported as `"filtered"` rather than `"unmapped"`. Because our synthetic
reads use Q25–Q40, the mapped counts will be similar to the unfiltered run
but may drop slightly depending on how many 3′-end bases fall below the
k-mer quality threshold.

---

## Complete one-liner (all features combined)

```bash
python main.py \
  -t dumpalign \
  -r test_data/bacteria.kdb \
  --reads test_data/reads.fastq \
  --coverage \
  --detect-variants \
  --min-variant-coverage 5 \
  --min-variant-quality 25 \
  --min-read-quality 25
```

---

## A note on ambiguous reads

The `ambiguous` classification occurs when the pseudo-aligner finds roughly
equal numbers of *species-specific* k-mers pointing to two or more reference
genomes in the same read. In practice this arises when:

1. **Two closely related strains** are in the reference database (e.g. *E. coli*
   K-12 and O157:H7) and a read comes from a conserved region that has slightly
   different flanking sequences in each strain.
2. **Horizontal gene transfer** has placed a genomic island in multiple species
   at different loci, and a read spans the island-to-unique-flank boundary.

In this test dataset every k-mer belongs to exactly one genome
(`multi_mapping_kmers = 0`) so **no ambiguous reads are expected**.  
To artificially trigger the ambiguous path for testing, pass
`--unique-threshold 200`: any genome pair whose specific k-mer counts differ
by fewer than 200 will be classified as ambiguous rather than unique.

```bash
python main.py \
  -t dumpalign \
  -r test_data/bacteria.kdb \
  --reads test_data/reads.fastq \
  --unique-threshold 200
```

With such a high threshold many reads that would normally be unique will
fall below it and be re-classified as ambiguous, demonstrating the threshold's
effect without requiring specially crafted input data.

---

## Summary checklist

After completing all steps, check the following:

- [ ] `bacteria.kdb` created (~360–380 KB)
- [ ] `dumpref` reports 0 `multi_mapping_kmers` for every genome
- [ ] `dumpref` includes `soft_masked_bases` and `soft_masked_fraction` fields (0 for these uppercase test genomes)
- [ ] `dumpalign` from `.aln` file shows **215 unique, 0 ambiguous, 30 unmapped**
- [ ] `dumpalign` Statistics includes `total_read_bases = 36750`, `soft_masked_read_bases = 0`, `soft_masked_read_fraction = 0.0`
- [ ] `ecoli_K12` has **95** unique reads, `salmonella_LT2` **80**, `bacillus_subtilis_168` **40**
- [ ] `--coverage` reports non-zero `covered_bases_unique` for all three genomes
- [ ] `--detect-variants` reports exactly **one variant** at `ecoli_K12` position `2000` (G→A, coverage 15)
- [ ] Increasing `--unique-threshold` shifts reads from unique → ambiguous
- [ ] A FASTA with lowercase bases (Step 5) shows correct non-zero `soft_masked_bases` and `soft_masked_fraction`
