FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/bioinf-tools/:$PATH
ENV LANG=C.UTF-8

ARG VIR_WF_DIR=/viridian_workflow
RUN mkdir -p $VIR_WF_DIR/.ci/
COPY .ci/install_dependencies.sh $VIR_WF_DIR/.ci/install_dependencies.sh
RUN $VIR_WF_DIR/.ci/install_dependencies.sh /bioinf-tools

COPY . $VIR_WF_DIR
RUN pip3 install tox \
  && cd /viridian_workflow \
  && tox \
  && pip3 install .

CMD viridian_workflow
