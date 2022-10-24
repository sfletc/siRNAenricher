# siRNAenricher
Python package for identifying genomic regions with siRNA enrichment

### Installation

Clone repository
```pip install -e .``` from root directory

### Usage

```help(sir.extract_enriched_seqs)```

```extract_enriched_seqs(scram_alignment_file, reference_fa, output_fa, window=200, cutoff=30, abund_count=5, padding=30)```
