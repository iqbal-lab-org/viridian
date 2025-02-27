FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/bioinf-tools/:/bioinf-tools/enaBrowserTools/python3/:$PATH
ENV LANG=C.UTF-8

ARG VIR_WF_DIR=/viridian
RUN mkdir -p $VIR_WF_DIR/.ci/
COPY .ci/install_dependencies.sh $VIR_WF_DIR/.ci/install_dependencies.sh
RUN $VIR_WF_DIR/.ci/install_dependencies.sh /bioinf-tools

COPY . $VIR_WF_DIR
RUN cd /viridian \
  && python3 -m pip install pytest \
  && python3 -m pip install -r requirements.txt \
  && pytest \
  && python3 -m pip install .

#RUN apt-get autoremove -y \
#  && apt-get purge -y --auto-remove \
#  && rm -rf /var/lib/apt/lists/*

CMD viridian
