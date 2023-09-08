# Use nvidia/cuda image
FROM nvidia/cuda:10.2-cudnn7-devel-ubuntu18.04

# set bash as current shell
RUN chsh -s /bin/bash
SHELL ["/bin/bash", "-c"]

# install anaconda
RUN apt-get update
RUN apt-get install -y wget bzip2 ca-certificates git libglib2.0-0 && \
        apt-get clean
	
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh -O ~/anaconda.sh && \
        /bin/bash ~/anaconda.sh -b -p /opt/conda && \
        rm ~/anaconda.sh && \
        ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
        echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
        find /opt/conda/ -follow -type f -name '*.a' -delete && \
        find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
        /opt/conda/bin/conda clean -afy

# set path to conda
ENV PATH /opt/conda/bin:$PATH

# setup conda virtual environment
COPY ./environment.yml /tmp/environment.yml
RUN conda update conda \
    && conda env create -f /tmp/environment.yml \
    && conda clean -afy

RUN echo "conda activate tangram-sc" >> ~/.bashrc
ENV PATH /opt/conda/envs/tangram-sc/bin:$PATH
ENV CONDA_DEFAULT_ENV $tangram-sc
