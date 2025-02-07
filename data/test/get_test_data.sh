git clone https://github.com/nextflow-io/training.git --branch=master --single-branch --depth=1 tmp
mv tmp/nf-training/data/ggal .
rm -rf tmp/
gzip ggal/gut_*
apptainer run docker://quay.io/biocontainers/slimfastq:2.04--h503566f_5 slimfastq ggal/liver_1.fq ggal/liver_1.sfq
apptainer run docker://quay.io/biocontainers/slimfastq:2.04--h503566f_5 slimfastq ggal/liver_2.fq ggal/liver_2.sfq
rm ggal/liver*fq
rm ggal/lung_1.fq
apptainer run docker://quay.io/biocontainers/kallisto:0.50.1--h6de1650_2 kallisto index -t 4 --index ggal/kallisto.idx ggal/transcriptome.fa
wget https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20240112.tar.gz
mkdir k2_minusb_20240112
cd k2_minusb_20240112/
tar xvzf ../k2_minusb_20240112.tar.gz
cd ../
rm k2_minusb_20240112.tar.gz
