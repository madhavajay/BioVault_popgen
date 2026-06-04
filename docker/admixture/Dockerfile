FROM condaforge/miniforge3:25.3.1-0

ARG ADMIXTURE_VERSION=1.4.0
ARG ADMIXTURE_URL=https://dalexander.github.io/admixture/binaries/admixture_linux-1.4.0.tar.gz

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates curl tar && \
    rm -rf /var/lib/apt/lists/*

RUN mamba install -y -c conda-forge -c bioconda \
      python=3.11 \
      numpy \
      pandas \
      matplotlib-base \
      plink2 && \
    mamba clean -afy

RUN mkdir -p /opt/admixture && \
    curl -fsSL "${ADMIXTURE_URL}" -o /tmp/admixture.tar.gz && \
    tar -xzf /tmp/admixture.tar.gz -C /opt/admixture --strip-components=1 && \
    find /opt/admixture -maxdepth 2 -type f -name admixture -exec install -m 0755 {} /usr/local/bin/admixture \; && \
    rm -f /tmp/admixture.tar.gz && \
    (admixture 2>&1 || true) | grep -q "ADMIXTURE Version ${ADMIXTURE_VERSION}"

COPY 03_individual_level/sex_biased_admixture_find_k/scripts /opt/biovault/scripts/sex_biased_admixture_find_k

ENV PATH=/opt/conda/bin:${PATH}

ENTRYPOINT ["/usr/local/bin/admixture"]
CMD ["--help"]
