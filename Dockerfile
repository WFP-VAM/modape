FROM continuumio/miniconda3

LABEL maintainer="valentin.pesendorfer@wfp.org"

ENV CONDA_ENV_PATH /opt/conda/envs/
ENV CONDA_ENV "python36"
ENV CPL_ZIP_ENCODING UTF-8

RUN conda install -y python=3.6

RUN apt-get update && apt-get install -y vim gcc build-essential aria2

RUN useradd -m worker

RUN conda install -y \
	pip h5py gdal

RUN conda clean -y -t

ADD . /home/worker
WORKDIR /home/worker

RUN python setup.py install
RUN python setup.py test

RUN rm -rf *

USER worker

CMD ["/bin/bash"]
