FROM python:3

WORKDIR /usr/src/app

RUN pip install --upgrade setuptools==47.1.1

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

WORKDIR /usr/src/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 \
    && tar -jxvf samtools-1.15.1.tar.bz2 \
    && cd samtools-1.15.1 \
    && ./configure --prefix=/usr/src/samtools \
    && make \
    && make install

ENV PATH /usr/src/samtools/bin:$PATH

COPY . .

CMD [ "/bin/bash" ]

