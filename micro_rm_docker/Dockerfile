# Nate Olson Bioinformatic Container for Micro RM Characterization

# Solazyme Bioinformatics Container, base on Timothy Laurent 
# https://github.com/timothyjlaurent/Docker-bioinformatics-image
#
# A portable analysis environment capable of running on any system or virtual machine
# with linux and Docker installed
# version 0.1

FROM phusion/baseimage

MAINTAINER Nate Olson

## Adding VOLUME for micro_rm
VOLUME ["/micro_rm"]

#install required libs
RUN apt-get install -y -q git build-essential autotools-dev automake \
                          pkg-config curl wget bwa unzip zlib1g-dev \
                          libncurses5-dev libyaml-dev cpanminus r-base

# setup Rprofile
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

## install Anaconda
RUN wget http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/Anaconda-1.9.2-Linux-x86_64.sh && bash Anaconda-1.9.2-Linux-x86_64.sh -b && rm -rf Anaconda-1.9.2-Linux-x86_64.sh;\
    mv /root/anaconda/lib/python2.7/platform.py /root/anaconda/lib/python2.7/platform.pybak && cp /usr/lib/python2.7/platform.py /root/anaconda/lib/python2.7/platform.py

ENV PATH $PATH:/root/anaconda/bin/
RUN mkdir /software && cd /software/

# install bowtie2 
RUN wget http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.2/bowtie2-2.2.2-linux-x86_64.zip;\
    unzip bowtie2-2.2.2-linux-x86_64.zip && rm -rf bowtie2-2.2.2-linux-x86_64.zip;\
    ln -s `pwd`/bowtie*/bowtie* /usr/local/bin

# install picard
# this also sets the jars as 'executable' and puts them in the path so they can be found with the 'which ' command
RUN wget http://downloads.sourceforge.net/project/picard/picard-tools/1.111/picard-tools-1.111.zip;\
    unzip picard-tools-1.111.zip && rm -rf picard-tools-1.111.zip;\
    chmod +x picard*/*.jar && ln -s `pwd`/picard*/*.jar /usr/local/bin/

# install htslib samtools and bcftools
RUN wget http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2;\
    tar -xaf samtools-1.1.tar.bz2 && rm -rf samtools-1.1.tar.bz2 ;\
    cd samtools-1.1;\
    make && ln -f -s `pwd`/* /usr/bin && cd /software

RUN wget http://sourceforge.net/projects/samtools/files/samtools/1.1/bcftools-1.1.tar.bz2;\
    tar -xaf bcftools-1.1.tar.bz2 && rm -rf bcftools-1.1.tar.bz2 ;\
    cd bcftools-1.1;\
    make && ln -f -s `pwd`/bcftools /usr/bin && cd /software

#sratoolkit
RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.2-1/sratoolkit.2.4.2-ubuntu64.tar.gz &&\
        tar -xf sratoolkit.2.4.2-ubuntu64.tar.gz &&\
        rm -rf sratoolkit.2.4.2-ubuntu64.tar.gz;\
        ln -f -s sratoolkit.2.4.2-ubuntu64/bin/fastq-dump.2.4.2 /usr/bin/fastq-dump;\
        cd /software

# ### to do setup ipython notebook
