FROM centos:latest

# Maintainer
MAINTAINER Intel-CCC

# update system
RUN yum -y update && yum -y install \
make \
java \
tar \
wget \
gcc-c++ \
glibc-static \
zlib-devel \
vim

# install STAR from https://github.com/alexdobin/STAR
RUN cd /opt && \
    wget -c -P /opt/star https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz && \
    cd /opt/star && \
    tar -xzf /opt/star/STAR_2.4.2a.tar.gz && \
    make STAR -C /opt/star/STAR-STAR_2.4.2a/source

# Add star to PATH
ENV PATH /opt/star/STAR-STAR_2.4.2a/bin/Linux_x86_64:$PATH
