ARG BASE_IMAGE
FROM ${BASE_IMAGE}
COPY matcha /usr/local/bin/matcha
# Biocontainers ship libhts.so under /usr/local/lib (conda prefix); not on
# the default loader path. Set LD_LIBRARY_PATH so dlopen() finds it.
ENV LD_LIBRARY_PATH=/usr/local/lib
RUN matcha --version && matcha --help > /dev/null
