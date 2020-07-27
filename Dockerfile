FROM ubuntu:bionic

LABEL maintainer="valentin.pesendorfer@wfp.org"

RUN apt-get update && apt-get install -y \
    gcc \
    build-essential \
    aria2 \
    software-properties-common \
    python3 \
    python3-dev \
    python3-pip \
  && rm -rf /var/lib/apt/lists/*

RUN add-apt-repository ppa:ubuntugis/ppa

RUN apt-get update && apt-get install -y \
    gdal-bin \
    python3-gdal \
  && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y libpq-dev \
   && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install cython

RUN useradd -m worker
ADD . /home/worker
WORKDIR /home/worker

RUN pip3 install .
RUN python3 setup.py test

RUN rm -rf *

USER worker

ENV LC_ALL=C.UTF-8 \
  LANG=C.UTF-8

CMD ["modape_version"]
