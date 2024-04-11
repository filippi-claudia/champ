# Start from the latest Ubuntu image
FROM ubuntu:latest

# Set work directory
WORKDIR /app

# Install necessary packages
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    cmake \
    gnupg \
    autoconf \
    libtool \
    python3 \
    emacs \
    gcc \
    g++ \
    gfortran \
    openmpi-bin \
    libopenmpi-dev \
    gawk \
    build-essential \
    libscalapack-openmpi-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Install HDF5
RUN curl -fsSLO 'https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_14_3/src/hdf5-1.14.3.tar' \
    && tar -xvf hdf5-1.14.3.tar \
    && cd hdf5-1.14.3 \
    && ./configure --prefix=/usr/local --enable-fortran --enable-parallel --enable-hl FC=mpif90 CC=mpicc  CXX=mpicxx \
    && make -j $(nproc) \
    && make install

# Install trexio
RUN curl -fsSLO 'https://github.com/TREX-CoE/trexio/archive/refs/tags/v2.4.2.tar.gz' \
    && tar -xzvf v2.4.2.tar.gz \
    && cd trexio-2.4.2 \
    && cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc \
    && cd build \
    && make -j $(nproc) \
    && make install

# Install CHAMP
RUN git clone https://github.com/filippi-claudia/champ.git \
    && cd champ \
    && cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90 \
        -DCMAKE_C_COMPILER=mpicc \
        -DENABLE_TREXIO=ON \
        -DTREXIO_INCLUDE_DIR=/usr/local/include \
        -DTREXIO_LIBRARY=/usr/local/lib/libtrexio.so \
    && cd build \
    && make -j $(nproc) \
    && cp /app/champ/bin/vmc.mov1 /app/champ/bin/dmc.mov1 /usr/local/bin/ \
    && cp /app/champ/bin/libparser.so /app/champ/bin/libpspline.a /usr/local/lib/
