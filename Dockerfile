FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.9 python3-pip python3-setuptools python3-dev
RUN apt-get -y update; apt-get -y install curl
RUN apt-get install -y gfortran
RUN apt-get install -y xz-utils

RUN apt-get install -y libbz2-dev apt-utils gcc-multilib software-properties-common

RUN apt-get update && apt-get install -y \
	libpq-dev \
	build-essential \
	libcurl4-gnutls-dev \
	libxml2-dev \
	libssl-dev \
	lbzip2 \
	libfftw3-dev \
	libgdal-dev \
	libgeos-dev \
	libgsl0-dev \
	libgl1-mesa-dev \
	libglu1-mesa-dev \
	libhdf4-alt-dev \
	libhdf5-dev \
	libjq-dev \
#	liblwgeom-dev \
	libpq-dev \
	libproj-dev \
	libprotobuf-dev \
	libnetcdf-dev \
	libsqlite3-dev \
	libudunits2-dev \
	netcdf-bin \
	postgis \
	protobuf-compiler \
	sqlite3 \
	tk-dev \
	unixodbc-dev \
	libxml2-dev \
	libcairo2-dev \
	libsqlite3-dev \
#	libmariadbd-dev \
	libpq-dev \
	libssh2-1-dev \
	unixodbc-dev \
#	libcurl4-openssl-dev \
	libssl-dev \
	libnlopt-dev \
	liblapack-dev \
	libblas-dev \
	libxt-dev \
	r-cran-tidyverse

RUN mkdir -p /opt/software/setup/R
ADD installPkgs.R /opt/software/setup/R/
RUN Rscript /opt/software/setup/R/installPkgs.R

RUN mkdir -p /opt/software/setup/python
ADD environment.yml /opt/software/setup/python/
RUN apt-get update && apt-get -y upgrade \
	&& apt-get install -y --no-install-recommends \
		git \
		wget \
		g++ \
		ca-certificates \
		&& rm -rf /var/lib/apt/lists/*
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
	&& mkdir /root/.conda \
	&& bash Miniconda3-latest-Linux-x86_64.sh -b \
	&& rm -f Miniconda3-latest-Linux-x86_64.sh \
	&& echo "Running $(conda --version)" && \
	conda init bash && \
	. /root/.bashrc && \
	conda update conda && \
	conda env create -f /opt/software/setup/python/environment.yml && \
	conda activate dnabert2 && \
	conda install python=3.9 pip && \
	conda install -c bioconda subread && \
	conda install -c bioconda samtools && \
	conda install -c bioconda bedtools && \
	conda install -c bioconda fastqc && \
	conda install -c bioconda star=2.7.9a && \
	conda install -c bioconda bowtie2

ENTRYPOINT [ "/bin/bash", "-l", "-c" ]

ENV DEBIAN_FRONTEND=interactive

RUN apt-get -y update
RUN apt-get -y install git
SHELL [ "/bin/bash", "-c" ]
RUN echo "source activate dnabert2" > /root/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

# CMD ["git clone https://github.com/venkatajonnakuti/PolyAMiner-Bulk.git"]

RUN mkdir -p /root/PolyAMiner-Bulk
ADD Demo /root/PolyAMiner-Bulk/Demo
ADD lib /root/PolyAMiner-Bulk/lib
ADD PolyA-miner.py /root/PolyAMiner-Bulk/
ADD README.md /root/PolyAMiner-Bulk/
ADD ReferenceFiles /root/PolyAMiner-Bulk/ReferenceFiles

RUN . /root/.bashrc && \
	cd /root/PolyAMiner-Bulk/lib/DNABERT && \
	Rscript /root/PolyAMiner-Bulk/lib/DNABERT/installPkgs.R && \
	python3 -m pip install --editable /root/PolyAMiner-Bulk/lib/DNABERT/ && \
	add-apt-repository universe && \
	dpkg --add-architecture i386 && \
	apt-get update && \
	apt-get install -y libncurses5 libncurses5:i386 && \
	cd /root/PolyAMiner-Bulk/

# ENV CONDA_DEFAULT_ENV dnabert2

# CMD ["/bin/sh", "-c", "conda activate dnabert2"]
# SHELL ["conda", "run", "-n", "dnabert2", "/bin/bash", "-c"]
RUN echo "\n" > /root/.bashrc \
&& echo -e "#! /bin/bash\n\n# script to activate the conda environment" > /root/.bashrc \
&& conda init bash \
&& echo -e "conda activate dnabert2\n" >> /root/.bashrc 

RUN apt-get install -y pigz unzip

SHELL ["/bin/bash", "-l", "-c"]
ENV BASH_ENV /root/.bashrc
# ENTRYPOINT ["conda", "run", "-n", "dnabert2", "python3" ,"/root/PolyAMiner-Bulk/PolyA-miner.py"]
# CMD python3 /root/PolyAMiner-Bulk/PolyA-miner.py
# ENTRYPOINT ["/bin/bash", "-l", "-c", "python3", "/root/PolyAMiner-Bulk/PolyA-miner.py"]
# ENV PATH /bin/bash:$PATH

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

