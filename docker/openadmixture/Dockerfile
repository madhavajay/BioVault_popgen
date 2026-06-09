ARG BASE_IMAGE=ghcr.io/madhavajay/openadmixture.jl:latest
FROM ${BASE_IMAGE}

ENV HGP1K_ADMIXTURE_REF=/opt/biovault/reference/hgp1k_admixture
RUN cp -a /root/.julia /opt/julia-depot && chmod -R a+rX /opt/julia-depot
ENV JULIA_DEPOT_PATH=/tmp/.julia:/opt/julia-depot

COPY flows/08_bv_paper_sex_biased_openadmixture/scripts /opt/biovault/scripts/sex_biased_openadmixture_hgp1k
COPY .docker/reference/openadmixture_hgp1k/ /opt/biovault/reference/hgp1k_admixture/
