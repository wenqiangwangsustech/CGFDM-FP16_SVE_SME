# Project Dependencies

This project requires several C/C++ libraries to be installed manually before building. Please follow the steps below.

## Step 1. Install libtiff SQLite3 PROJ

```bash
wget http://download.osgeo.org/libtiff/tiff-4.0.7.tar.gz
tar -xvf tiff-4.0.7.tar.gz
cd tiff-4.0.7
./configure --prefix=/your/install/path/tiff
make -j
make install
wget https://www.sqlite.org/2024/sqlite-autoconf-3460000.tar.gz
tar -xvf sqlite-autoconf-3460000.tar.gz
cd sqlite-autoconf-3460000
./configure --prefix=/your/install/path/sqlite3
make -j
make install
wget https://download.osgeo.org/proj/proj-8.1.0.tar.gz
tar -xvf proj-8.1.0.tar.gz
cd proj-8.1.0

export TIFF_CFLAGS="-I/your/install/path/tiff/include"
export TIFF_LIBS="-L/your/install/path/tiff/lib -ltiff"
export SQLITE3_CFLAGS="-I/your/install/path/sqlite3/include"
export SQLITE3_LIBS="-L/your/install/path/sqlite3/lib -lsqlite3"

## Step 1. Install libtiff SQLite3 PROJ
## Step 2. Install Anaconda
conda create -n pyenv python=3.10 -y
conda activate pyenv


## Step 3. Install Python libararies
