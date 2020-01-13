FROM ubuntu:bionic
LABEL maintainer="valentin.pesendorfer@wfp.org"

RUN apt-get update && apt-get install -y \
    gcc \
    build-essential \
    aria2 \
    software-properties-common \
    python3.6 \
    python3.6-dev \
    python3-pip \
  && rm -rf /var/lib/apt/lists/*

RUN add-apt-repository ppa:ubuntugis/ppa

RUN apt-get update && apt-get install -y \
    gdal-bin=2.4.2+dfsg-1~bionic0 \
    python3.6-gdal \
  && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y libpq-dev \
   && rm -rf /var/lib/apt/lists/*

RUN useradd -m worker
ADD . /home/worker
WORKDIR /home/worker

RUN python3 setup.py install
RUN python3 setup.py test

RUN rm -rf *

USER worker

CMD ["modape_version"]
