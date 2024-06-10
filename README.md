# DE-Preprocessing

## Description

Estimate the expression level of features (genes or transcripts) and perform the associated QC. Starting from cleaned reads and one pre-built reference database in standard format (cf. Readme.md).


| header  | sample_id | path_R1                        | path_R2                                        | strand                                               | kall_single_overhang | meanFlD                            | sdFlD                              |
|---------|-----------|--------------------------------|------------------------------------------------|------------------------------------------------------|----------------------|------------------------------------|------------------------------------|
| **content** | string    | path *.{fastq,fq,spring}{.gz,} | same as 'path_R1', but empty if   SE or spring | array ["unstranded",   "fr-stranded', "rf-stranded"] | boolean              | float                              | float                              |
| **status**  | recquired | recquired                      | optional                                       | recquired                                            | recquired            | optional (arbitrary default value) | optional (arbitrary default value) |
