# Prebuilt Binder image for operator_gb.
#
# This is NOT the Dockerfile Binder builds -- it is built and pushed to GHCR by
# .github/workflows/binder-image.yml, and binder/Dockerfile then adds a thin
# layer on top of it. Building the heavy parts (Sage image, CaDiCaL, the Cython
# extensions) here keeps mybinder.org launches fast.
#
# Base image: the official SageMath Binder environment, see
# https://github.com/sagemath/sage-binder-env
# Ubuntu 24.04, Sage in /home/sage/sage, user `sage` with uid 1000,
# jupyterlab + notebook preinstalled, linux/amd64 only.
FROM ghcr.io/sagemath/sage-binder-env:10.9

LABEL org.opencontainers.image.source=https://github.com/ClemensHofstadler/operator_gb
LABEL org.opencontainers.image.description="SageMath 10.9 with operator_gb installed, for mybinder.org"

USER root

# The base image is a *runtime* image: it ships gcc/g++/pkg-config, but neither
# `make` (needed by setup.py to build the bundled CaDiCaL) nor the Python
# headers (needed to compile the Cython extensions -- /usr/include/python3.12
# does not exist in the base image).
RUN apt-get -qq update \
    && apt-get -qq install -y --no-install-recommends make python3-dev \
    && apt-get -qq clean \
    && rm -rf /var/lib/apt/lists/*

# Build and install as uid 1000, so nothing ends up root-owned inside the Sage
# venv that the Binder user later owns.
USER sage

COPY --chown=sage:sage . /tmp/operator_gb

# sage.env.sage_include_directories() does not include $SAGE_LOCAL, but
# rational_linear_algebra.pyx cimports sage.libs.gmp.mpq, which pulls in
# `#include <gmp.h>` and `# distutils: libraries = gmp`. Sage's own GMP lives in
# $SAGE_LOCAL, which is on no default search path, so point the compiler and the
# linker at it.
# --no-build-isolation makes pip build against the installed Sage library, as
# required by setup.py (it imports sage.env).
RUN SAGE_LOCAL="$(sage -sh -c 'echo $SAGE_LOCAL')" \
    && CPATH="$SAGE_LOCAL/include${CPATH:+:$CPATH}" \
       LIBRARY_PATH="$SAGE_LOCAL/lib${LIBRARY_PATH:+:$LIBRARY_PATH}" \
       sage -pip install --no-cache-dir --no-build-isolation /tmp/operator_gb \
    && rm -rf /tmp/operator_gb

# Fail the image build loudly rather than shipping something broken.
RUN sage -c "\
from operator_gb import *; \
import os, operator_gb; \
assert os.path.isfile(os.path.join(os.path.dirname(operator_gb.__file__), 'cadical')), 'CaDiCaL binary missing'; \
F = FreeAlgebra(QQ, ['a', 'b', 'a_adj', 'b_adj']); \
a, b, a_adj, b_adj = F.gens(); \
proof = certify(add_adj(pinv(a, b, a_adj, b_adj)), a*b*a - a); \
print('smoke test OK, proof has %d steps' % len(proof))"

WORKDIR /home/sage
