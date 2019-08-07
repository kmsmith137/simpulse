# build image
FROM python:2.7-slim AS build-img

# set the $(HOME) environment variable
ENV HOME /usr

# add folder contents
ADD . /simpulse

# set working directory
WORKDIR /simpulse

# install necessary dependencies for making simpulse
# and then make & install simpulse
RUN set -ex \
    && apt-get update \
    && apt-get install -y build-essential git cmake autoconf \
    libtool pkg-config fftw-dev libfftw3-dev \
    && pip install --no-cache-dir -r requirements.txt \
    && make all install


# FROM python2.7-slim
# 
# RUN set -ex \
#     && apt-get update \
#     && apt-get install --no-install-recommends -y curl git \
#     && rm -rf \
#         /var/cache/debconf/*-old \
#         /var/lib/apt/lists/* \
#         /tmp/* \
#         /var/tmp/* \
#         /usr/share/man \
#         /usr/share/doc \
#         /usr/share/doc-base
