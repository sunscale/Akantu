FROM alpine:edge
LABEL maintainer "Nicolas Richart <nicolas.richart@epfl.ch>"

WORKDIR /usr/src/app

RUN apk --update add --no-cache --upgrade \
    clang\
    clang-extra-tools \
    musl-dev \
    python3 \
    py3-lxml && \
    rm -rf /usr/share/ri && \
    adduser -u 9000 -D -s /bin/false app

#COPY engine.json /
COPY . ./
RUN chown -R app:app ./
USER app

VOLUME /code
WORKDIR /code

CMD ["/usr/src/app/bin/codeclimate-clang-tidy"]
