# docker build -t jackgisby/fastqc .
# docker push jackgisby/fastqc

# get base image from biocontainers
FROM biocontainers/biocontainers:v1.1.0_cv2

# install FASTQC
# based on https://hub.docker.com/r/biocontainers/fastqc/dockerfile

USER root

RUN apt-get update && apt-get install -y openjdk-8-jre-headless
RUN mkdir -p /opt/fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip --no-check-certificate -O /opt/fastqc/fastqc_v0.11.9.zip
RUN unzip /opt/fastqc/fastqc_v0.11.9.zip -d /opt/fastqc
RUN rm /opt/fastqc/fastqc_v0.11.9.zip
RUN mv /opt/fastqc/FastQC/* /opt/fastqc
RUN rmdir /opt/fastqc/FastQC && chmod +x /opt/fastqc/fastqc
RUN ln -s /opt/fastqc/fastqc /usr/local/bin/fastqc

# add programs to path
ENV PATH /usr/local/bin:$PATH
