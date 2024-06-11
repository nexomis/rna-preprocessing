# DE-Preprocessing

## Description

Estimate the expression level of features (genes or transcripts) and perform the associated QC. Starting from cleaned reads and one pre-built reference database in standard format (cf. Readme.md).

## input sample sheet

| header  | sample_id | path_R1                        | path_R2                                        | strand                                               | kall_single_overhang | mean_FlD                            | sd_FlD                              |
|---------|-----------|--------------------------------|------------------------------------------------|------------------------------------------------------|----------------------|------------------------------------|------------------------------------|
| **content** | string    | path *.{fastq,fq,spring}{.gz,} | same as 'path_R1', but empty if   SE or spring | array ["unstranded",   "fr-stranded', "rf-stranded"] | boolean              | float                              | float                              |
| **status**  | recquired | recquired                      | optional                                       | recquired                                            | recquired            | optional: default value (totally arbitrary!) in case of SE, automatically set in case of PE | optional (arbitrary default value) |
| **description**  | sample nid (must be uniq) | path of fastq (can be springed) file | path of second fastq file in case of PE library and not spring format | library orientation | active or not '--single-overhang' kallisto option (must be actived in case of UTRseq !) | mean of fragment length distribution (preferred to fill in if 'sd_FlD' is filled) | standard deviation of fragment length distribution (preferred to fill in if 'mean_FlD' is filled) |


## Test data sets

### Chicken: Fastq files + transcriptome.fa (subset of chr1)
```
git clone https://github.com/nextflow-io/training.git --branch=master --single-branch --depth=1 tmp
mv tmp/nf-training/data/ggal data
rm -rf tmp/
rm -rf data/liver_1.fq
```

### Kallisto index for the transcriptome

#### Option 1: Building from the downloaded transcriptome (only subset of chr1 chicken)
```
docker run -v $PWD/data/:/data/ -w /data/ -u $(id -u):$(id -g) quay.io/biocontainers/kallisto:0.50.1--h6de1650_2 kallisto index -t 4 --index /data/kallisto.idx /data/transcriptome.fa
```

#### Option 2: Downloading (about full Chichen transriptome)
```
wget -O data/chicken_index_standard.tar.gz https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/gallus_gallus.tar.gz
tar -xzvf data/chicken_index_standard.tar.gz -C tmp/ \
tar -xzvf data/chicken_index_standard.tar.gz --strip-components=1 -C data gallus_gallus/transcriptome.idx \
  && rm data/chicken_index_standard.tar.gz
mv data/transcriptome.idx data/kallisto.idx
```

### Example Usage:
```
nextflow run rna-preprocessing/main.nf --input_dir data/ --kallisto_idx data/kallisto.idx -resume -profile docker --samplesheet data/samplesheet.csv
```
