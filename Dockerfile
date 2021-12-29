FROM continuumio/miniconda:latest

# Set working directory and copy tDRnamer software into docker container
COPY . /opt/tdrnamer/
ENV PATH /opt/tdrnamer:$PATH

# Install Core tDRnamer Dependencies
RUN conda env update -n base -f /opt/tdrnamer/tdrnamer_env.yaml

# Add empty folder for RNA database docker volumes
RUN mkdir /rnadb && \
     chmod -R 777 /rnadb && \
     chmod -R 777 /home

# Add user for the container
RUN useradd -ms /bin/bash jerry
USER jerry
WORKDIR /home/jerry
