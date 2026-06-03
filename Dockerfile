FROM condaforge/miniforge3:24.11.3-0 AS tools

SHELL ["/bin/bash", "-lc"]

ENV CONDA_ENV=biovault_popgen
ENV PATH=/opt/conda/envs/biovault_popgen/bin:/opt/conda/bin:$PATH
ENV MPLCONFIGDIR=/tmp/matplotlib
ENV XDG_CACHE_HOME=/tmp/.cache
ENV LOADINGS_HT=/opt/biovault/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht
ENV LOADINGS_VARIANTS_TSV=/opt/biovault/reference/pca_loadings/loadings_variants.tsv
ENV LOADINGS_NPZ=/opt/biovault/reference/pca_loadings/loadings.npz
ENV HGP1K_MATRIX_NPZ=/opt/biovault/reference/hgp1k/hgp1k_dosage.npz
ENV HGP1K_METADATA_TSV=/opt/biovault/reference/hgp1k/20130606_g1k_3202_samples_ped_population.txt
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

FROM pgscatalog/plink2:2.00a5.10 AS plink2-bin

FROM ghcr.io/openmined/biosynth:0.1.31 AS biosynth-bin

FROM tools AS runtime

ARG HGP1K_REFERENCE_SUBDIR=matrix

RUN mkdir -p /opt/biovault/reference/pca_loadings /opt/biovault/reference/aims /opt/biovault/reference/hgp1k /opt/biovault/scripts && \
    chmod 666 /etc/passwd /etc/group

COPY .docker/reference/pca_loadings/loadings_variants.tsv ${LOADINGS_VARIANTS_TSV}
COPY .docker/reference/pca_loadings/loadings.npz ${LOADINGS_NPZ}
COPY .docker/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht.tar.gz /tmp/gnomad.v3.1.pca_loadings.ht.tar.gz
RUN tar -xzf /tmp/gnomad.v3.1.pca_loadings.ht.tar.gz -C /opt/biovault/reference/pca_loadings && \
    rm -f /tmp/gnomad.v3.1.pca_loadings.ht.tar.gz
COPY .docker/reference/aims/gnomad_af_per_locus.tsv /opt/biovault/reference/aims/gnomad_af_per_locus.tsv
COPY .docker/reference/hgp1k/${HGP1K_REFERENCE_SUBDIR}/ /opt/biovault/reference/hgp1k/
COPY .docker/reference/hgp1k/20130606_g1k_3202_samples_ped_population.txt ${HGP1K_METADATA_TSV}
COPY tools /opt/biovault/tools
COPY 00_qc_all_files/scripts /opt/biovault/scripts/qc_all_files
COPY 03_individual_level/gnomad_projection/scripts /opt/biovault/scripts/gnomad_projection
COPY 03_individual_level/gnomad_projection_fast/scripts /opt/biovault/scripts/gnomad_projection_fast
COPY 03_individual_level/hgp1k_projection_fast/scripts /opt/biovault/scripts/hgp1k_projection_fast
COPY 03_individual_level/pca_qc_fast/scripts /opt/biovault/scripts/pca_qc_fast
COPY 03_individual_level/sex_biased_admixture/scripts /opt/biovault/scripts/sex_biased_admixture
COPY 03_individual_level/sex_biased_admixture_fast/scripts /opt/biovault/scripts/sex_biased_admixture_fast
COPY 03_individual_level/sex_biased_admixture_find_k/scripts /opt/biovault/scripts/sex_biased_admixture_find_k
COPY 04_population_level/fst_aims_fast/scripts /opt/biovault/scripts/population_level
# NOTE: flow definitions are intentionally NOT baked. Only analysis script
# trees live in the image. The population flow's split step runs in the
# biosynth container (just `bvs` calls, inlined in main.nf — no baked script
# needed there); its FST/AIMs step consumes the baked population_level/.
RUN chmod +x /opt/biovault/scripts/gnomad_projection/*.sh \
             /opt/biovault/scripts/gnomad_projection_fast/*.sh \
             /opt/biovault/scripts/hgp1k_projection_fast/*.sh

WORKDIR /work

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["bash"]

FROM condaforge/miniforge3:25.3.1-0 AS fast-tools

SHELL ["/bin/bash", "-lc"]

ENV CONDA_ENV=biovault_popgen
ENV PATH=/opt/conda/envs/biovault_popgen/bin:/opt/conda/bin:$PATH
ENV MPLCONFIGDIR=/tmp/matplotlib
ENV XDG_CACHE_HOME=/tmp/.cache
ENV LOADINGS_NPZ=/opt/biovault/reference/pca_loadings/loadings.npz
ENV HGP1K_MATRIX_NPZ=/opt/biovault/reference/hgp1k/hgp1k_dosage.npz
ENV HGP1K_METADATA_TSV=/opt/biovault/reference/hgp1k/20130606_g1k_3202_samples_ped_population.txt

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      bash \
      ca-certificates \
      coreutils \
      findutils \
      gawk \
      gzip \
      less \
      libblas3 \
      libgfortran5 \
      liblapack3 \
      procps \
      sed \
      time \
      tini && \
    rm -rf /var/lib/apt/lists/*

COPY environment.fast.yml /tmp/environment.fast.yml

RUN conda env create -p "/opt/conda/envs/${CONDA_ENV}" -f /tmp/environment.fast.yml && \
    conda clean -afy && \
    rm -f /tmp/environment.fast.yml

FROM fast-tools AS fast-runtime

ARG HGP1K_REFERENCE_SUBDIR=matrix

RUN mkdir -p /opt/biovault/reference/pca_loadings /opt/biovault/reference/aims /opt/biovault/reference/hgp1k /opt/biovault/scripts && \
    chmod 666 /etc/passwd /etc/group

COPY --from=plink2-bin /usr/local/bin/plink2 /usr/local/bin/plink2
COPY --from=biosynth-bin /usr/local/bin/bvs /usr/local/bin/bvs
COPY --from=biosynth-bin /app/data/genostats.sqlite /app/data/genostats.sqlite
ENV BVS_READ_ONLY_DB=1

COPY .docker/reference/pca_loadings/loadings.npz ${LOADINGS_NPZ}
COPY .docker/reference/aims/gnomad_af_per_locus.tsv /opt/biovault/reference/aims/gnomad_af_per_locus.tsv
COPY .docker/reference/hgp1k/${HGP1K_REFERENCE_SUBDIR}/ /opt/biovault/reference/hgp1k/
COPY .docker/reference/hgp1k/20130606_g1k_3202_samples_ped_population.txt ${HGP1K_METADATA_TSV}
COPY tools /opt/biovault/tools
COPY 00_qc_all_files/scripts /opt/biovault/scripts/qc_all_files
COPY 03_individual_level/gnomad_projection_fast/scripts /opt/biovault/scripts/gnomad_projection_fast
COPY 03_individual_level/hgp1k_projection_fast/scripts /opt/biovault/scripts/hgp1k_projection_fast
COPY 03_individual_level/pca_qc_fast/scripts /opt/biovault/scripts/pca_qc_fast
COPY 03_individual_level/sex_biased_admixture/scripts /opt/biovault/scripts/sex_biased_admixture
COPY 03_individual_level/sex_biased_admixture_fast/scripts /opt/biovault/scripts/sex_biased_admixture_fast
COPY 03_individual_level/sex_biased_admixture_find_k/scripts /opt/biovault/scripts/sex_biased_admixture_find_k
COPY 04_population_level/fst_aims_fast/scripts /opt/biovault/scripts/population_level

RUN chmod +x /opt/biovault/scripts/gnomad_projection_fast/*.sh \
             /opt/biovault/scripts/hgp1k_projection_fast/*.sh

WORKDIR /work

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["bash"]
