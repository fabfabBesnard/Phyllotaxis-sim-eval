FROM r-base:4.2.0 
#starting from the r-base image
LABEL maintener="Fabrice Besnard <fabrice.besnard@ens-lyon.fr>"
LABEL repository="https://github.com/fabfabBesnard/Phyllotaxis-sim-eval"
LABEL description="Environment to simulate perturbed phyllotaxis sequences and their alignment to a ground truths sequence with sm-dtw"
LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Dirk Eddelbuettel <edd@debian.org>"
#source: https://towardsdatascience.com/creating-sandbox-environments-for-r-with-docker-def54e3491a3

# Make ~/.R 
RUN mkdir -p $HOME/.R

# Change Shell to 'bash', default is 'sh'
SHELL [ "/bin/bash", "-c" ]

# Main Installs
#       - 'apt-get' installs (misc programs + python3)
#       - R + R packages installs (optparse, ggplot2, reshape2, gridExtra, IRkernel for Jupyter notebooks)
#       - pip installs: jupyter notebooks and bash+R kernels
RUN apt-get update -qq && \
    apt install python3 -y && \
    apt-get -y --no-install-recommends install \
        liblzma-dev \
        libbz2-dev \
        clang  \
        ccache \
        g++ \
        gcc \
        default-jdk \
        default-jre \
        git \
        wget \
        python3-pip && \
    R CMD javareconf && \
    install2.r --error \
        ggstance ggrepel ggthemes \
        ###My packages are below this line
        optparse ggplot2 reshape2 gridExtra IRkernel && \
    rm -rf /tmp/downloaded_packages/ /tmp/*.rds && \
    rm -rf /var/lib/apt/lists/*  && \
    python3 -m pip install --upgrade pip && \
    pip install notebook && \
    pip install ipykernel && \
    pip install bash_kernel && \
    python3 -m bash_kernel.install && \
    R -e "IRkernel::installspec()"

# Build a context directory (e.g.: https://newbedev.com/what-is-app-working-directory-for-a-dockerfile)
WORKDIR /myapp

# Conda install for sm-dtw: create a virtual environment called "dtw" with miniconda
# the sm-dtw is installed in this environment using the conda packaging designed for this program
# This virtual environment will be activated automatically at start in the container
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && echo "Running $(conda --version)" && \
    conda init bash && \
    . /root/.bashrc && \
    conda update conda && \
    git clone https://github.com/romi/dtw.git && \
    git clone https://github.com/fabfabBesnard/Phyllotaxis-sim-eval.git && \
    cd dtw && \
    conda env create -f conda/env/dtw.yaml && \
    conda activate dtw && \
    python -m pip install -e .

RUN echo 'conda activate dtw' >> /root/.bashrc