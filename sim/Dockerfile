FROM gcc:7.2

RUN wget http://www.met.reading.ac.uk/clouds/adept/adept-2.0.5.tar.gz && \
    tar xvfz adept-2.0.5.tar.gz && \
    cd adept-2.0.5 && \
    ./configure CXXFLAGS="-g1 -O3 -march=native" && \
    make -j2 && make install

RUN git clone https://github.com/stanford-ssi/balloons-VALBAL balloons-VALBAL
ENV VALBAL_DIR balloons-VALBAL

RUN mkdir sim
WORKDIR sim
ADD . .
RUN make