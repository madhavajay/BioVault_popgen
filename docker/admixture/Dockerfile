FROM condaforge/miniforge3:25.3.1-0 AS tools

ARG ADMIXTURE_VERSION=1.4.0
ARG ADMIXTURE_URL=https://dalexander.github.io/admixture/binaries/admixture_linux-1.4.0.tar.gz

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates curl tar && \
    rm -rf /var/lib/apt/lists/*

# plink (1.9, for --bmerge), plink2, bcftools (reference pre-gen), python libs.
RUN mamba install -y -c conda-forge -c bioconda \
      python=3.11 \
      numpy \
      pandas \
      matplotlib-base \
      scikit-learn \
      plink \
      plink2 \
      bcftools && \
    mamba clean -afy

RUN mkdir -p /opt/admixture && \
    curl -fsSL "${ADMIXTURE_URL}" -o /tmp/admixture.tar.gz && \
    tar -xzf /tmp/admixture.tar.gz -C /opt/admixture --strip-components=1 && \
    find /opt/admixture -maxdepth 2 -type f -name admixture -exec install -m 0755 {} /usr/local/bin/admixture \; && \
    rm -f /tmp/admixture.tar.gz && \
    (admixture 2>&1 || true) | grep -q "ADMIXTURE Version ${ADMIXTURE_VERSION}"

# Baked analysis scripts (find_k kept for back-compat; hgp1k is the new flow).
COPY 03_individual_level/sex_biased_admixture_find_k/scripts /opt/biovault/scripts/sex_biased_admixture_find_k
COPY 03_individual_level/sex_biased_admixture_hgp1k/scripts  /opt/biovault/scripts/sex_biased_admixture_hgp1k

ENV PATH=/opt/conda/bin:${PATH}
ENV HGP1K_ADMIXTURE_REF=/opt/biovault/reference/hgp1k_admixture

# The flow runs as the host UID (non-root); main.nf appends a passwd/group
# entry so subprocess pwd-name lookups resolve. Make those files writable —
# same as the main biovault-popgen image.
RUN chmod 666 /etc/passwd /etc/group

ENTRYPOINT ["/usr/local/bin/admixture"]
CMD ["--help"]

# ---------------------------------------------------------------------------
# Final image: bake the pre-generated HGP1K ADMIXTURE reference BED.
# Generate it first with:
#   ./build_admixture_docker.sh        (auto-generates .docker/reference/hgp1k_admixture
#                                        via the `tools` stage if missing)
# ---------------------------------------------------------------------------
FROM tools AS runtime

COPY .docker/reference/hgp1k_admixture/ /opt/biovault/reference/hgp1k_admixture/
