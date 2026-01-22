# STARsolo reference

## 10x droplet (CB_UMI_Simple / Droplet)

Minimal example:
```
STAR --genomeDir /path/to/genomeDir \
  --readFilesIn Read2.fastq.gz Read1.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist /path/to/whitelist.txt
```

Notes:
- Always pass cDNA read first, barcode/UMI read second.
- For multiple lanes, pass comma-separated read lists in the same order.
- v2 chemistry defaults: CB=16, UMI=10.
- v3 chemistry: add `--soloUMIlen 12`.

### Match CellRanger 3.x
```
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR
```

### Match CellRanger 4/5
```
--clipAdapterType CellRanger4 --outFilterScoreMin 30 \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR
```

## Reference build guidance

Use the same FASTA/GTF as CellRanger to match counts.
Example:
```
STAR --runMode genomeGenerate \
  --genomeDir /path/to/genomeDir \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/genes.gtf
```
Optional CellRanger-like index:
```
--genomeSAsparseD 3
```

## Barcode geometry

Simple barcodes:
```
--soloCBstart 1 --soloCBlen 16 \
--soloUMIstart 17 --soloUMIlen 10
```

Barcode + cDNA on same mate (example: 10x 5'):
```
--soloBarcodeMate 1 --clip5pNbases 39 0 \
--soloType CB_UMI_Simple \
--soloCBstart 1 --soloCBlen 16 \
--soloUMIstart 17 --soloUMIlen 10 \
--readFilesIn read1.fq read2.fq
```

Complex barcodes:
```
--soloType CB_UMI_Complex \
--soloCBposition 0_0_2_-1 3_1_3_8 \
--soloUMIposition 3_9_3_14 \
--soloAdapterSequence <adapter>
```

## Cell filtering (cell calling)

Knee filtering (default):
```
--soloCellFilter CellRanger2.2
```

EmptyDrops-like filtering:
```
--soloCellFilter EmptyDrops_CR
```

Filter existing raw matrix without remapping:
```
STAR --runMode soloCellFiltering /path/to/count/dir/raw/ /path/to/output/prefix \
  --soloCellFilter EmptyDrops_CR
```

## Feature quantification

- Genes only (default): `--soloFeatures Gene`
- Pre-mRNA (exon+intron): `--soloFeatures GeneFull`
- Splice junctions: `--soloFeatures SJ`
- Velocyto-like: `--soloFeatures Gene Velocyto`
- All together: `--soloFeatures Gene GeneFull SJ Velocyto`

## Multi-gene reads

Choose one or more strategies with `--soloMultiMappers`:
- `Uniform`: split UMIs uniformly across genes.
- `PropUnique`: weighted by unique UMI counts.
- `EM`: expectation-maximization.
- `Rescue`: weighted by unique + uniform multi-gene.

## BAM input

Use BAM as input and point STARsolo to barcode tags:
```
STAR --readFilesIn input.bam --readFilesType SAM SE \
  --readFilesCommand samtools view -F 0x100 \
  --soloInputSAMattrBarcodeSeq CR UR \
  --soloInputSAMattrBarcodeQual CY UY
```

Avoid duplicate tags in output BAM:
```
--readFilesSAMattrKeep None
```

## BAM tags and output

Add tags to output BAM:
```
--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
```

For corrected CB/UB tags, output coordinate-sorted BAM:
```
--outSAMtype BAM SortedByCoordinate
```

## SmartSeq (plate-based)

Use separate FASTQs per cell and manifest file:
```
--soloType SmartSeq \
--readFilesManifest /path/to/manifest.tsv \
--soloUMIdedup Exact \
--soloStrand Unstranded
```

Manifest format (paired-end):
```
Read1.fastq.gz<TAB>Read2.fastq.gz<TAB>CellID
```

Manifest format (single-end):
```
Read1.fastq.gz<TAB>-<TAB>CellID
```
