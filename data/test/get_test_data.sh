wget -q https://zenodo.org/record/6457007/files/GSM461177_1_subsampled.fastqsanger
wget -q https://zenodo.org/record/6457007/files/GSM461177_2_subsampled.fastqsanger
wget -q https://zenodo.org/record/6457007/files/GSM461180_1_subsampled.fastqsanger
wget -q https://zenodo.org/record/6457007/files/GSM461180_2_subsampled.fastqsanger

mv GSM461177_1_subsampled.fastqsanger WT_1.fq
mv GSM461177_2_subsampled.fastqsanger WT_2.fq
mv GSM461180_1_subsampled.fastqsanger MUT_1.fq
mv GSM461180_2_subsampled.fastqsanger MUT_2.fq

apptainer run docker://quay.io/biocontainers/slimfastq:2.04--h503566f_5 slimfastq WT_1.fq WT_1.sfq
apptainer run docker://quay.io/biocontainers/slimfastq:2.04--h503566f_5 slimfastq WT_2.fq WT_2.sfq

gzip *.fq

wget https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20240112.tar.gz
mkdir k2_minusb_20240112
cd k2_minusb_20240112/
tar xvzf ../k2_minusb_20240112.tar.gz
cd ../
rm k2_minusb_20240112.tar.gz

#wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz
#wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz

wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_rna.fna.gz
mv GCF_000001215.4_Release_6_plus_ISO1_MT_rna.fna.gz transcriptome.fa.gz
apptainer run docker://quay.io/biocontainers/kallisto:0.50.1--h6de1650_2 kallisto index -i idx transcriptome.fa.gz

