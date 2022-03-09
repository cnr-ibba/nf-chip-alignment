#
# VERSION 0.1
# DOCKER-VERSION  20.10.12
# AUTHOR:         Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
# DESCRIPTION:    A python container with some dependencies
#

# start from a python base image
FROM python:3.9

LABEL maintainer="Paolo Cozzi <paolo.cozzi@ibba.cnr.it>" \
    description="A python container with some dependencies" \
    version="0.1"

# This environment variable force stdin, stdout and stderr to be totally unbuffered
ENV PYTHONUNBUFFERED 1

# Elegantly activating a virtualenv in a Dockerfile
# https://pythonspeed.com/articles/activate-virtualenv-dockerfile/
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Upgrade pip to latest version
RUN python3 -m pip install --upgrade pip wheel

# Install pypi modules (uwsgi)
COPY requirements.txt /root/
RUN pip install --no-cache-dir -r /root/requirements.txt

# Execute a default command
CMD ["/bin/bash"]
