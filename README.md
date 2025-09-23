# Project Dependencies

This project requires several dependencies to be installed before building and running. Please follow the steps below.

---


## Step 1. Install libtiff, SQLite3, and PROJ (from source)

```bash
# Install libtiff
wget http://download.osgeo.org/libtiff/tiff-4.0.7.tar.gz
tar -xvf tiff-4.0.7.tar.gz
cd tiff-4.0.7
./configure --prefix=/your/install/path/tiff
make -j
make install

# Install SQLite3
wget https://www.sqlite.org/2024/sqlite-autoconf-3460000.tar.gz
tar -xvf sqlite-autoconf-3460000.tar.gz
cd sqlite-autoconf-3460000
./configure --prefix=/your/install/path/sqlite3
make -j
make install

# Install PROJ
wget https://download.osgeo.org/proj/proj-8.1.0.tar.gz
tar -xvf proj-8.1.0.tar.gz
cd proj-8.1.0
````


Then configure, build, and install PROJ:

```bash
./configure --prefix=/your/install/path/proj
make -j
make install
```

---

## Step 2. Install Anaconda and Create Environment

Download and install Anaconda (if not installed already), then create a dedicated environment:

```bash
conda create -n pyenv python=3.10 -y
conda activate pyenv
```

---

## Step 3. Install Python Libraries

Inside the `pyenv` environment, install the required Python libraries:

```bash
conda install pyproj -y
```

(Additional dependencies can be installed using `pip install -r requirements.txt` if provided.)

---

After completing these steps, all required C/C++ and Python dependencies are installed, and the environment `pyenv` is ready for building and running this project.


## Step 4. Setup Build and Run Environment:

```bash
export TIFF_CFLAGS="-I/your/install/path/tiff/include"
export TIFF_LIBS="-L/your/install/path/tiff/lib -ltiff"
export SQLITE3_CFLAGS="-I/your/install/path/sqlite3/include"
export SQLITE3_LIBS="-L/your/install/path/sqlite3/lib -lsqlite3"
```



