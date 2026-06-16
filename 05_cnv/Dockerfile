FROM romanhaa/penncnv:latest

LABEL org.opencontainers.image.title="BioVault CNV"
LABEL org.opencontainers.image.description="IDAT-to-PennCNV workflow wrapper. Mount a licensed IAAP CLI binary at runtime."

SHELL ["/bin/bash", "-lc"]

ENV PENNCNV_HOME=/opt/PennCNV-1.0.5
ENV PATH=/opt/PennCNV-1.0.5:/opt/illumina/bin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV IAAP_BIN=/opt/illumina/bin/iaap-cli

RUN mkdir -p /opt/biovault-cnv/scripts /opt/illumina/bin /input /out

COPY scripts/ /opt/biovault-cnv/scripts/

RUN chmod +x /opt/biovault-cnv/scripts/*.sh

WORKDIR /work

ENTRYPOINT ["/opt/biovault-cnv/scripts/run_pipeline.sh"]
