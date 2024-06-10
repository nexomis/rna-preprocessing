# DE-Preprocessing

## Description

Estimate the expression level of features (genes or transcripts) and perform the associated QC. Starting from cleaned reads and one pre-built reference database in standard format (cf. Readme.md).

## input sample sheet

| header  | sample_id | path_R1                        | path_R2                                        | strand                                               | kall_single_overhang | mean_FlD                            | sd_FlD                              |
|---------|-----------|--------------------------------|------------------------------------------------|------------------------------------------------------|----------------------|------------------------------------|------------------------------------|
| **content** | string    | path *.{fastq,fq,spring}{.gz,} | same as 'path_R1', but empty if   SE or spring | array ["unstranded",   "fr-stranded', "rf-stranded"] | boolean              | float                              | float                              |
| **status**  | recquired | recquired                      | optional                                       | recquired                                            | recquired            | optional: default value (totally arbitrary!) in case of SE, automatically set in case of PE | optional (arbitrary default value) |
| **description**  | sample nid (must be uniq) | path of fastq (can be springed) file | path of second fastq file in case of PE library and not spring format | library orientation | active or not '--single-overhang' kallisto option (must be actived in case of UTRseq !) | mean of fragment length distribution (preferred to fill in if 'sd_FlD' is filled) | standard deviation of fragment length distribution (preferred to fill in if 'mean_FlD' is filled) |
