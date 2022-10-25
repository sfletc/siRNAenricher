# siRNAenricher
Python package for identifying genomic regions with siRNA enrichment. Takes a SCRAM2 alignment file and genome reference (FASTA) as inputs and writes enriched regions to file (FASTA).

### Installation

Clone repository
```pip install -e .``` from root directory

### Usage

```help(sir.extract_enriched_seqs)```

```extract_enriched_seqs(scram_alignment_file, reference_fa, output_fa, window=200, cutoff=30, abund_count=5, strand_ratio=0.2, padding=30)```

```scram_alignment_file```: Output file of a single alignment length from the SCRAM2 aligner.  Generally pick 21 or 22 nt files for siRNAs.

```reference_fa```: Genome reference file in FASTA format (the same reference as used by the SCRAM2 aligner)

```output_fa```: Output file to results identified regions in FASTA format

```window```: window size in nt to scan in.  This will be automatically expanded if enriched windows are adjacent (within len window).  Default = 200.

```cutoff```: Minimum alignment count for an siRNA to be included in identification. Default = 30.

```abund_count```: Minimum number of the siRNAs with an alignment count above cutoff for a window to be identified as an enriched region. Default = 5.

```strand_ration```: Minimum ratio of abund_count of siRNAs on the lower abundance strand.  Prevents identification of degraded single-stranded RNA regions. Default = 0.2.

```padding```: number of nucleotides to added to each end of the window for writing to file.  Default = 30.

