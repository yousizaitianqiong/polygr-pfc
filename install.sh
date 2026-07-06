#!/bin/bash

set -euo pipefail

VERSION=3.3.10
PREFIX=${PREFIX:-"$HOME/.local/fftw-$VERSION-mpi"}

wget "https://www.fftw.org/fftw-$VERSION.tar.gz"
tar -xzf "fftw-$VERSION.tar.gz"
cd "fftw-$VERSION"

./configure \
    --prefix="$PREFIX" \
    --enable-mpi \
    --enable-threads \
    --enable-openmp \
    --enable-shared

make -j "${JOBS:-$(nproc)}"
make install

echo "Installed FFTW MPI/thread libraries in $PREFIX"
