ARG BASE_IMAGE
FROM ${BASE_IMAGE}
COPY matcha /usr/local/bin/matcha
RUN matcha --version && matcha --help > /dev/null
