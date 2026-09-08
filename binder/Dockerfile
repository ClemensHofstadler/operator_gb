# Dockerfile for binder
# Reference: https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html
#
# The heavy lifting (SageMath + operator_gb + CaDiCaL) happens in
# binder/base.Dockerfile, which CI builds and pushes to GHCR. This file only
# adds the notebook and the Binder-specific user setup, so a launch is quick.
#
# The user-handling below is taken from https://github.com/sagemath/sage-binder-env
FROM ghcr.io/clemenshofstadler/operator_gb-binder:latest

USER root

# Create user with uid 1000 if it doesn't exist, or rename sage if it exists
ARG NB_USER=user
ARG NB_UID=1000
ENV NB_USER=${NB_USER}
ENV NB_UID=${NB_UID}
ENV HOME=/home/${NB_USER}
ENV PATH="/home/sage/sage:/home/sage/sage/local/bin:/home/sage/sage/venv/bin:$PATH"

RUN if id -u ${NB_UID} >/dev/null 2>&1; then \
        EXISTING_USER=$(id -nu ${NB_UID}); \
        EXISTING_GROUP=$(id -gn ${NB_UID}); \
        if [ "$EXISTING_GROUP" != "$NB_USER" ]; then \
            groupmod -n ${NB_USER} ${EXISTING_GROUP} || true; \
        fi; \
        if [ "$EXISTING_USER" != "$NB_USER" ]; then \
            usermod -l ${NB_USER} -d ${HOME} ${EXISTING_USER}; \
            mkdir -p ${HOME}; \
            chown ${NB_USER}:${NB_USER} ${HOME}; \
        fi; \
    else \
        groupadd -g ${NB_UID} ${NB_USER} || true; \
        useradd -m -s /bin/bash -u ${NB_UID} -g ${NB_USER} ${NB_USER}; \
    fi

# The tutorial is what Binder opens; the case study is the full reference
COPY --chown=${NB_USER}:${NB_USER} Tutorial.ipynb Case-Study.ipynb ${HOME}/

# Switch to the user
USER ${NB_USER}

# Make Sage's Jupyter kernel visible to the notebook server.
#
# JupyterLab resolves a notebook's kernel by kernelspec *name*, and a notebook
# carries the name of whatever install it was last saved with: Case-Study.ipynb
# asks for `sagemath`, while Tutorial.ipynb (saved with the macOS Sage app)
# asks for `sagemath-10.9`. Sage's kernel here is only called `sagemath`, and an
# unresolved name makes the server fall back to the plain python3 kernel without
# warning, so every Sage cell fails. Expose the kernel under `sagemath` plus
# every `sagemath*` name the shipped notebooks request.
RUN set -eu; \
    kernels_dir="$(/usr/bin/sage -sh -c 'jupyter --data-dir')/kernels"; \
    sage_kernel="$(/usr/bin/sage -sh -c 'echo $SAGE_VENV')/share/jupyter/kernels/sagemath"; \
    test -f "$sage_kernel/kernel.json"; \
    mkdir -p "$kernels_dir"; \
    names="$(python3 -c "import glob, json; print(' '.join(sorted({'sagemath'} | {n for f in glob.glob('${HOME}/*.ipynb') for n in [json.load(open(f)).get('metadata', {}).get('kernelspec', {}).get('name', '')] if n.startswith('sagemath')})))")"; \
    for name in $names; do \
        ln -sfn "$sage_kernel" "$kernels_dir/$name"; \
    done; \
    echo "Sage kernel registered as: $names"

# Start in the home directory of the user
WORKDIR ${HOME}

# Create the jupyter_lab_config.py file with a custom logging filter to
# suppress the perpetual nodejs warning
RUN mkdir -p ${HOME}/.jupyter
RUN echo "\
import logging\n\
\n\
class NoNodeJSWarningFilter(logging.Filter):\n\
    def filter(self, record):\n\
        return 'Could not determine jupyterlab build status without nodejs' not in record.getMessage()\n\
\n\
logging.getLogger('LabApp').addFilter(NoNodeJSWarningFilter())\n\
" > ${HOME}/.jupyter/jupyter_lab_config.py

# Should a notebook still ask for a kernelspec this image does not have, start
# Sage rather than the python3 kernel, so Sage code runs instead of erroring.
# jupyter_server_config.py is read by both jupyter-lab and notebook 7.
RUN echo "c.MultiKernelManager.default_kernel_name = 'sagemath'" \
    > ${HOME}/.jupyter/jupyter_server_config.py
