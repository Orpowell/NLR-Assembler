FROM ubuntu:20.04
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y python3

COPY requirements.txt requirements.txt

RUN apt-get install -y python3-pip
RUN pip3 install -r requirements.txt

ENTRYPOINT ["python3"]