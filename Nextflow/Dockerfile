FROM mambaorg/micromamba:1.5.8

# Set HOME to a writable directory for Nextflow
ENV HOME=/home/mamba
USER root

RUN mkdir -p /home/mamba /pipeline \
    && chmod -R 777 /home/mamba /pipeline

COPY /docker/base_env.yaml /tmp/base_env.yaml

RUN micromamba create -y -n pypeline -f /tmp/base_env.yaml && \
    micromamba install -y -n pypeline -c bioconda -c conda-forge nextflow && \
    micromamba clean --all --yes

SHELL ["micromamba", "run", "-n", "pypeline", "/bin/bash", "-lc"]

ENTRYPOINT ["micromamba", "run", "-n", "pypeline", "bash"]