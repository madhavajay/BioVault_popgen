FROM condaforge/miniforge3:24.11.3-0 AS tools

SHELL ["/bin/bash", "-lc"]

ENV CONDA_ENV=biovault_popgen
ENV PATH=/opt/conda/envs/biovault_popgen/bin:/opt/conda/bin:$PATH
ENV MPLCONFIGDIR=/tmp/matplotlib
ENV XDG_CACHE_HOME=/tmp/.cache
ENV LOADINGS_HT=/opt/biovault/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht
ENV LOADINGS_VARIANTS_TSV=/opt/biovault/reference/pca_loadings/loadings_variants.tsv
ENV GCS_CONNECTOR_JAR=/opt/hadoop/gcs-connector-hadoop3.jar

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      bash \
      ca-certificates \
      coreutils \
      findutils \
      gawk \
      gzip \
      less \
      procps \
      sed \
      time \
      tini \
      wget && \
    rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml

RUN conda env create -p "/opt/conda/envs/${CONDA_ENV}" -f /tmp/environment.yml && \
    conda clean -afy && \
    rm -f /tmp/environment.yml

COPY 03_individual_level/gnomad_projection/scripts/extract_loadings_variants.py /tmp/extract_loadings_variants.py

RUN mkdir -p /opt/hadoop && \
    wget -q -O "${GCS_CONNECTOR_JAR}" \
      https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop3-latest.jar

FROM tools AS runtime

RUN mkdir -p /opt/biovault/reference/pca_loadings /opt/biovault/scripts && \
    chmod 666 /etc/passwd /etc/group

COPY .docker/reference/pca_loadings/loadings_variants.tsv ${LOADINGS_VARIANTS_TSV}
COPY .docker/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht ${LOADINGS_HT}
COPY 03_individual_level/gnomad_projection/scripts /opt/biovault/scripts/gnomad_projection
COPY 03_individual_level/gnomad_projection_fast/scripts /opt/biovault/scripts/gnomad_projection_fast
COPY 03_individual_level/pca_qc_fast/scripts /opt/biovault/scripts/pca_qc_fast
RUN chmod +x /opt/biovault/scripts/gnomad_projection/*.sh \
             /opt/biovault/scripts/gnomad_projection_fast/*.sh

WORKDIR /work

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["bash"]
