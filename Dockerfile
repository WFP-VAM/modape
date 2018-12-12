FROM continuumio/miniconda3:4.3.11

LABEL maintainer="valentin.pesendorfer@wfp.org"

ENV CONDA_ENV_PATH /opt/conda/envs/
ENV CONDA_ENV "python35"
ENV CPL_ZIP_ENCODING UTF-8

RUN conda install -y python=3.5

RUN apt-get update && apt-get install -y vim gcc build-essential

RUN useradd -m worker

#RUN conda update --quiet --yes conda
RUN conda install -y \
	pip numpy cython pytables h5py scipy cython gdal

RUN conda clean -y -t

ADD . /home/worker
WORKDIR /home/worker

RUN python setup.py install

RUN rm -rf *

USER worker

CMD ["/bin/bash"]
