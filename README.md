# DE-Preprocessing

## Description

Estimate the expression level of features (genes or transcripts) and perform the associated QC. Starting from cleaned reads and one pre-built reference database in standard format (cf. Readme.md).

## input sample sheet

| header  | sample_id | path_R1                        | path_R2                                        | strand                                               | kall_single_overhang | mean_FlD                            | sd_FlD                              |
|---------|-----------|--------------------------------|------------------------------------------------|------------------------------------------------------|----------------------|------------------------------------|------------------------------------|
| **content** | string    | path *.{fastq,fq,spring}{.gz,} | same as 'path_R1', but empty if   SE or spring | array ["unstranded",   "fr-stranded', "rf-stranded"] | boolean              | float                              | float                              |
| **status**  | recquired | recquired                      | optional                                       | recquired                                            | recquired            | optional: default value (totally arbitrary!) in case of SE, automatically set in case of PE | optional (arbitrary default value) |
| **description**  | sample nid (must be uniq) | path of fastq (can be springed) file | path of second fastq file in case of PE library and not spring format | library orientation | active or not '--single-overhang' kallisto option (must be actived in case of UTRseq !) | mean of fragment length distribution (preferred to fill in if 'sd_FlD' is filled) | standard deviation of fragment length distribution (preferred to fill in if 'mean_FlD' is filled) |


# Downloading the dataset

## Fastq files + transcriptome.fa
```
git clone https://github.com/nextflow-io/training.git --branch=master --single-branch --depth=1 tmp
mv tmp/nf-training/data/ggal data
rm -rf tmp/
rm -rf data/liver_1.fq
```

## Kallisto index for the transcriptome (2 options)
### Option 1: Downloading
```
docker run -v data/:/data/ -w /data/ -u $(id -u):$(id -g) quay.io/biocontainers/kallisto:0.50.1--h6de1650_2 kallisto index -t 4 --index /data/kallisto.idx /data/transcriptome.fa
```

### Option 2: Building from the downloaded transcriptome
```
wget -O data/kallisto/human_index_standard.tar.xz https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_standard.tar.xz
tar -xf data/human_index_standard.tar.xz -C data/ && rm data/human_index_standard.tar.xz && mv data/index.idx data/kallisto.idx
```

## Example Usage:
```
nextflow run rna-preprocessing/main.nf --input_dir data/ --kallisto_idx data/kallisto.idx -resume -profile docker --samplesheet data/samplesheet.csv
```
