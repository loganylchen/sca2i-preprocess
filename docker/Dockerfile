# Reproducible analysis / visualization environment for sca2i-preprocess.
#
# Build:
#   docker build -f docker/Dockerfile -t sca2i-analyze:latest .
#
# Run (mount repo + results into /workspace):
#   docker run --rm -it -v "$PWD:/workspace" -w /workspace sca2i-analyze:latest \
#       python tools/analyze_with_scanpy.py
#
# Notes:
# - micromamba base keeps the image small and fast to build.
# - PYTHONDONTWRITEBYTECODE=1 avoids junk .pyc files on the host bind mount.
# - PYTHONNOUSERSITE=1 mirrors the production rule (avoids host ~/.local pandas
#   shadowing inside the container if /home is mounted by accident).

FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="sca2i-analyze"
LABEL org.opencontainers.image.description="scanpy + pyranges + matplotlib env for sca2i-preprocess QC / visualization"
LABEL org.opencontainers.image.source="https://github.com/loganylchen/sca2i-preprocess"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl git \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONNOUSERSITE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /workspace
CMD ["python", "-c", "import scanpy, pyranges, anndata; print(scanpy.__version__, pyranges.__version__, anndata.__version__)"]
