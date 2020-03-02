FROM ubuntu:19.04

WORKDIR /app

# Install GCC and the GNU Scientific Library.
RUN apt-get update -y \
    && apt-get install -y build-essential libgsl-dev

ADD Makefile main.cpp /app/

RUN make compile \
    && make clean

CMD ["/app/main.exec"]
