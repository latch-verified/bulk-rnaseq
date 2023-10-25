FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:dd8f-main

RUN apt-get update --yes &&\
    apt-get install --yes --no-install-recommends wget software-properties-common dirmngr

# >>> Install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN debian_codename=$(lsb_release --codename --short) &&\
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/debian ${debian_codename}-cran40/"

run wget https://github.com/r-lib/rig/releases/download/latest/rig-linux-latest.tar.gz
run tar \
    --extract \
    --gunzip \
    --file rig-linux-latest.tar.gz \
    --directory /usr/local/
run rm rig-linux-latest.tar.gz
run rig add release

# >>> R packages
env R_PKG_SYSREQS2="true"

COPY txImports.R /root/txImports.R
RUN /root/txImports.R
COPY lcImports.R /root/lcImports.R
RUN /root/lcImports.R

RUN apt-get install -y curl vim default-jre-headless zlib1g zlib1g-dev unzip cmake
RUN python3 -m pip install cutadapt
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz &&\
    tar xvzf trim_galore.tar.gz &&\
    mv TrimGalore-0.6.6/trim_galore /bin &&\
    rm -rf TrimGalore-0.6.6

RUN curl -L https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz -o 2.7.10a.tar.gz &&\
    tar xvzf 2.7.10a.tar.gz &&\
    mv STAR-2.7.10a/bin/Linux_x86_64_static/STAR /bin &&\
    rm -rf 2.7.10a 2.7.10a.tar.gz

RUN curl -L https://github.com/COMBINE-lab/salmon/releases/download/v1.8.0/salmon-1.8.0_linux_x86_64.tar.gz -o salmon-1.8.0_linux_x86_64.tar.gz &&\
  tar xvzf salmon-1.8.0_linux_x86_64.tar.gz &&\
  mv salmon-1.8.0_linux_x86_64/bin/salmon /bin &&\
  mv salmon-1.8.0_linux_x86_64/lib/* /lib/x86_64-linux-gnu/ &&\
  rm -rf salmon-1.8.0 salmon-1.8.0_linux_x86_64.tar.gz

RUN curl -L https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.3.tar.gz -o v1.3.3.tar.gz &&\
    tar -xzvf v1.3.3.tar.gz &&\
    cd RSEM-1.3.3 &&\
    make

RUN curl -s https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -O &&\
    unzip *.zip &&\
    chmod u+x ./FastQC/fastqc

RUN curl -L https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 -o samtools-1.16.tar.bz2 &&\
    tar -vxjf samtools-1.16.tar.bz2 &&\
    cd samtools-1.16 &&\
    ./configure --without-curses &&\
    make &&\
    make install

# RUN apt-get upgrade --yes liblzma-dev
RUN curl -L https://tukaani.org/xz/xz-5.2.6.tar.gz -o xz-5.2.6.tar.gz &&\
    tar -xzvf xz-5.2.6.tar.gz &&\
    cd xz-5.2.6 &&\
    ./configure --enable-shared &&\
    make &&\
    make install &&\
    ldconfig
RUN curl -L https://github.com/griffithlab/regtools/archive/refs/tags/0.5.2.tar.gz -o 0.5.2.tar.gz &&\
    tar -vxzf 0.5.2.tar.gz  &&\
    cd regtools-0.5.2 &&\
    mkdir build &&\
    cd build &&\
    cmake ..  &&\
    make &&\
    mv regtools /bin/

RUN wget https://sourceforge.net/projects/libpng/files/zlib/1.2.9/zlib-1.2.9.tar.gz/download &&\
    tar -xzvf download &&\
    cd zlib-1.2.9 &&\
    ./configure && make && make install

COPY gentrome.sh /root/gentrome.sh

RUN python3 -m pip install --upgrade multiqc matplotlib numpy scipy lgenome

RUN pip install latch==2.36.0
COPY wf/ /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
