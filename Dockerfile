FROM r-base:4.4.2

# wget, python, samtools and bedtools
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    python3 python3-pip \
    samtools bedtools \
    libcurl4-openssl-dev \
    libssl-dev \
    pkg-config \
    pandoc \
    && apt-get clean

# pysam via pip (force in system instead of venv)
RUN pip install pysam --break-system-packages \
    && apt-get remove --purge -y python3-pip

# R packages
RUN Rscript -e "install.packages(c('ggplot2', 'dplyr', 'gggenomes', 'patchwork', 'cowplot', 'plotly', 'htmlwidgets'))"

# specific R packages (IRanges package via BiocManager)
RUN Rscript -e "install.packages('BiocManager')" \
    && Rscript -e "BiocManager::install('IRanges')" \
    && Rscript -e "remove.packages('BiocManager')"

# copy scripts from specific commit repository
ENV REPO_NAME=plot_pe_mapping_status_from_bam
ENV COMMIT_HASH=546f701b99bc8e266e6c3a8202a3a8a33956bc8a
RUN wget -O ${REPO_NAME}.tar.gz https://github.com/nexomis/${REPO_NAME}/archive/${COMMIT_HASH}.tar.gz
RUN tar -xvzf ${REPO_NAME}.tar.gz \
    && rm ${REPO_NAME}.tar.gz
RUN cp -r ${REPO_NAME}-${COMMIT_HASH}/scripts /usr/src/${REPO_NAME}-${COMMIT_HASH} \
    && ln -s /usr/src/${REPO_NAME}-${COMMIT_HASH}/* /usr/local/bin/ \
    && rm -r ${REPO_NAME}-${COMMIT_HASH}/

# default rocker img run R as entrypoint
CMD ["bash"]
