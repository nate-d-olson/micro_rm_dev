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
RUN apt-get update;\
    apt-get install -y -q git build-essential autotools-dev automake \
                          pkg-config curl wget bwa unzip zlib1g-dev \
                          libncurses5-dev libyaml-dev cpanminus r-base
    # mkdir /install && cd install;\
    # apt-get install python-setuptools python-dev -y -q;\
 #    cd /install ;\
 #    git clone https://github.com/timothyjlaurent/cloudbiolinux.git ;\
 #    cd cloudbiolinux ;\
 #    git checkout 1a1fc199f8a3729e0231dd8413e33bbea1a0f1b9 ;\
 #    ls ;\
 #    pwd;\
 #    python setup.py build ;\
 #    sudo python setup.py install;\
 #    cd cloudbiolinux;\
 #    fab -f fabfile.py -H localhost -c config/fabricrc.txt install_biolinux:packages; \
 #    fab -f fabfile.py -H localhost install_biolinux:libraries && apt-get clean ;\
	# fab -f fabfile.py -H localhost -c install_biolinux:cleanup ;\
	# rm -rf /var/lib/apt/lists/* /var/tmp/* /var/tmp/*

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
# install reapr
# install Perl modules pre-reqs
RUN cpanm -v File::Basename File::Copy File::Spec  File::Spec::Link Getopt::Long List::Util
RUN wget ftp://ftp.sanger.ac.uk/pub/resources/software/reapr/Reapr_1.0.17.tar.gz;\
         tar -xf Reapr_1.0.17.tar.gz && rm -rf Reapr_1.0.17.tar.gz;\
    cd Reapr_1.0.17;\
    bash install.sh && ln -f -s `pwd`/src/reapr.pl /usr/bin/reapr && cd /software

# pathoscope2
RUN wget http://sourceforge.net/projects/pathoscope/files/pathoscope_2.0.2_clinical_pathoscope_1.0.3.tar.gz  &&\
        tar -xf pathoscope_2.0.2_clinical_pathoscope_1.0.3.tar.gz &&\
        rm -rf pathoscope_2.0.2_clinical_pathoscope_1.0.3.tar.gz ;\
    chmod -R 755 pathoscope2 && cd /software
## not sure how to change to put in path

#sratoolkit
RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.2-1/sratoolkit.2.4.2-ubuntu64.tar.gz &&\
        tar -xf sratoolkit.2.4.2-ubuntu64.tar.gz &&\
        rm -rf sratoolkit.2.4.2-ubuntu64.tar.gz;\
        ln -f -s sratoolkit.2.4.2-ubuntu64/bin/fastq-dump.2.4.2 /usr/bin/fastq-dump;\
        cd /software


# ### to do setup ipython notebook
# need to install sratoolkit
