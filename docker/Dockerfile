FROM mambaorg/micromamba:1.5.8

WORKDIR /workspace

COPY docker/base_env.yaml /tmp/environment.yaml

RUN micromamba create -y -n pypeline -f /tmp/base_env.yaml && \
    micromamba clean --all --yes

SHELL ["micromamba", "run", "-n", "pypeline", "/bin/bash", "-lc"]

ENTRYPOINT ["micromamba", "run", "-n", "pypeline", "snakemake"]
CMD ["--help"]