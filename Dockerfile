
# get base image from biocontainers
FROM biocontainers/biocontainers:latest

# install FASTQC
# based on https://hub.docker.com/r/biocontainers/fastqc/dockerfile
ENV ZIP=fastqc_v0.11.5.zip

RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc//$ZIP -O /tmp/$ZIP && \
  unzip - /tmp/$ZIP -d /tmp && \
  rm /tmp/$ZIP && \
  cd /tmp/FastQC && \
  chmod 755 fastqc && \
  ln -s /tmp/FastQC/fastqc /usr/local/bin/fastqc

# add programs to path
ENV PATH /usr/local/bin:$PATH
