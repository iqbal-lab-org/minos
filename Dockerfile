FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/bioinf-tools/:$PATH
ENV LANG=C.UTF-8
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

ARG MINOS_DIR=/minos
RUN mkdir -p $MINOS_DIR/.ci/
COPY .ci/install_dependencies.sh $MINOS_DIR/.ci/install_dependencies.sh
RUN $MINOS_DIR/.ci/install_dependencies.sh /bioinf-tools

COPY . $MINOS_DIR
RUN cd $MINOS_DIR \
  && pip3 install tox \
  && cd $MINOS_DIR \
  && tox \
  && pip3 install .

CMD minos
