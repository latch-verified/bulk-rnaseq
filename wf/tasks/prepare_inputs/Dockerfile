FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:dd8f-main

RUN apt-get update --yes &&\
    apt-get install --yes --no-install-recommends wget software-properties-common dirmngr

RUN apt-get install --yes libcurl4-openssl-dev libgsl-dev libssl-dev libxml2-dev

RUN apt-get install -y curl vim default-jre-headless zlib1g zlib1g-dev unzip cmake

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

COPY vendor/gentrome.sh /root/gentrome.sh

RUN python3 -m pip install --upgrade multiqc matplotlib numpy scipy lgenome

RUN python3 -m pip uninstall latch -y
COPY latch /root/latch
RUN cd latch && pip install -e .
COPY wf/ /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
