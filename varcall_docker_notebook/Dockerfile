FROM ipython/ipython

#MAINTAINER IPython Project <ipython-dev@scipy.org>
MAINTAINER Nate Olson <nolson@nist.gov>
#modified ipython/notebook

# Adding Bioinformatics Tools
RUN apt-get install -y -q git automake \
                          pkg-config curl wget bwa unzip zlib1g-dev \
                          libncurses5-dev libyaml-dev openjdk-7-jdk openjdk-7-jre
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
    make && ln -f -s `pwd`/* /usr/bin && cd ../

RUN wget http://sourceforge.net/projects/samtools/files/samtools/1.1/bcftools-1.1.tar.bz2;\
    tar -xaf bcftools-1.1.tar.bz2 && rm -rf bcftools-1.1.tar.bz2 ;\
    cd bcftools-1.1;\
    make && ln -f -s `pwd`/bcftools /usr/bin && cd ../

# install Plion
RUN wget https://github.com/broadinstitute/pilon/releases/download/v1.10/pilon-1.10.jar;\
    chmod +x pilon-1.10.jar; ln -s `pwd`/pilon-1.10.jar /usr/local/bin/

VOLUME /notebooks
WORKDIR /notebooks

EXPOSE 8888

# You can mount your own SSL certs as necessary here
ENV PEM_FILE /key.pem
# $PASSWORD will get `unset` within notebook.sh, turned into an IPython style hash
ENV PASSWORD micro_rm_pw

ADD notebook.sh /
RUN chmod u+x /notebook.sh

CMD ["/notebook.sh"]

