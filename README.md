## Step 1. Install libtiff, sqite3, and proj (from source)

This project requires several dependencies to be installed before building and running. Please follow the steps below.
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
All the Makefile* files:
```
CCHOME   := /your/compiler/direcory
MPIHOME  := /your/mpi/direcory
PROJHOME := /your/install/proj-8.1.0
```
All the run* files:
```
export LD_LIBRARY_PATH=/your/install/proj-8.1.0/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/your/install/sqlite3/lib:${LD_LIBRARY_PATH}
export PROJ_LIB=/your/install/proj-8.1.0/share/proj
```

## Step 5. Run provided examples

There are four examples:

* `HalfSpaceModel`
* `GassianShapedModel`
* `PerformanceTest`
* `TibetEarthquake`

To run an example, go into its directory and execute:

```bash
cd HalfSpaceModel
bash runAll
```
Each example generates output data that can be visualized using the provided Python scripts.
```
```

