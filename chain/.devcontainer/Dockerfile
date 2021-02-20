FROM ubuntu:20.04
ENV TZ=UTC
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update \
  && apt-get install -y build-essential pkg-config \
  && apt-get install -y libssl-dev libcurl4-openssl-dev nginx \
  && apt-get install -y swig sqlite3 \
  && apt-get install -y python3-dev libsqlite3-dev python3-setuptools \
  && apt-get install -y libtiff-dev \
  && apt-get install -y bison flex unzip \
  && apt-get install -y libnetcdf-dev libzstd-dev \
  && apt-get install -y python3-pip python3-numpy python3-venv \
  && apt-get install -y virtualenvwrapper \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 --no-cache-dir install --upgrade pip \
  && rm -rf /var/lib/apt/lists/*

# RUN cd /usr && chmod a+w src && cd /usr/src && rsync –-recursive /var/storage/mediakit/dependencies/* .
COPY dependencies/* /usr/src/

RUN cd /usr/src && tar -xzf hdf-4.2.15.tar.gz && cd hdf-4.2.15 \
  && export CFLAGS="-fPIC" && export CXXFLAGS="-fPIC" && ./configure --disable-netcdf --disable-fortran \
  && make && make install && cp -R hdf4 /opt \
  && ldconfig

RUN cd /usr/src && tar -xzf proj-7.2.0.tar.gz && cd proj-7.2.0 \
  && ./configure && make && make install \
  && ldconfig

RUN cd /usr/src && tar -jxf geos-3.8.1.tar.bz2 && cd geos-3.8.1 \
  && export PYTHON=/usr/bin/python3 \
  && ./configure --enable-python && make && make install \
  && ldconfig

RUN cd /usr/src && tar -xzf gdal-3.2.0.tar.gz && cd gdal-3.2.0 \
  && ./configure --with-hdf4=/opt/hdf4 --with-python --with-netcdf --with-liblzma --with-libzstd \
  && CPPFLAGS=-I/opt/hdf4/include make && make install \
  && ldconfig

RUN pip install numpy h5py>=2.9 pandas>=0.24
RUN pip install cython gunicorn flask gdal==3.2.0

RUN rm -vf /var/lib/apt/lists/* && apt-get update && apt-get install -y git

ARG USERNAME=worker
ARG USER_UID=1000
ARG USER_GID=$USER_UID
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --no-log-init --uid $USER_UID --gid $USER_GID -m $USERNAME
