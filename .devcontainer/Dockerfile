FROM condaforge/miniforge3:26.1.1-3
RUN mamba install -y -c bioconda -c conda-forge \
    bcftools tabix cmake make c-compiler 'nim>=2.2.6' \
    htslib python \
    && mamba clean -afy
