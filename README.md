# DE-Preprocessing

## Description

Estimate the expression level of features (genes or transcripts) and perform the associated QC.

## Testing

```
cd data/test
bash get_test_data.sh
cd ../..

nextflow run . --input data/test/samplesheet.csv --kallisto_idx data/test/ggal/kallisto.idx --kraken2_db data/test/k2_minusb_20240112/
```

## Sample sheet 

The sample sheet should define the input samples with their read files, that can either be pair-end or single end (and potentially compressed with spring).

#### Required: Sample Name `sample_name`

Sample name must contains only alphanumeric chars or underscore.

#### Required: Read 1 file `path_r1`

Path for R1 must be provided (either fq fq.gz or spring), cannot contain spaces and must end with .fastq, .fq, .spring optionally followed by .gz

#### Optional: Read 2 file `path_r2`

Path for R2 must be provided (either fq fq.gz), cannot contain spaces and must end with .fastq, .fq optionally followed by .gz
Can be left empty.

#### Optional: Strand information `strand`

Strand information must be provided and be one of '' same as 'unstranded', 'fr-stranded', 'rf-stranded'
Default is "", synonym of unstranded.

#### Optional: Is 3 prime ? `is_3prime`

Whether the sample has been sequenced using 3 prime library.
Default is false.

#### Optional: Mean fragment size `frag_size`

Mean fragment size (default not given)

#### Optional: Standard deviation of fragment size `frag_size_sd`

Standard deviation of fragment size (default not given)

#### Optional: Ignore size info for EM `force_kallisto_single_overhang`

Force `--single-overhang` option for `kallisto`
Default is true.


