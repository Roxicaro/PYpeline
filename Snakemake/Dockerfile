FROM mambaorg/micromamba:1.5.8

ENV HOME=/home/mamba
USER root

RUN mkdir -p /home/mamba /pipeline \
    && chmod -R 777 /home/mamba /pipeline

COPY docker/base_env.yaml /tmp/base_env.yaml

# Create micromamba env + install snakemake + awscli
RUN micromamba create -y -n pypeline -f /tmp/base_env.yaml && \
    micromamba install -y -n pypeline -c bioconda -c conda-forge snakemake awscli && \
    micromamba clean --all --yes

SHELL ["micromamba", "run", "-n", "pypeline", "--", "bash", "-c"]

# Default: end inside a shell (useful for local execution)
ENTRYPOINT ["micromamba", "run", "-n", "pypeline"]

# If no command is passed (local use), open shell.
CMD ["bash"]