Project Dependencies

This project requires several dependencies to be installed before building and running. Follow these steps to set up your environment.

1. Install System Dependencies from Source

libtiff Installation

# Download and extract libtiff
wget http://download.osgeo.org/libtiff/tiff-4.0.7.tar.gz
tar -xvf tiff-4.0.7.tar.gz
cd tiff-4.0.7

# Configure, build and install
./configure --prefix=/your/install/path/tiff
make -j$(nproc)
sudo make install


SQLite3 Installation

# Download and extract SQLite3
wget https://www.sqlite.org/2024/sqlite-autoconf-3460000.tar.gz
tar -xvf sqlite-autoconf-3460000.tar.gz
cd sqlite-autoconf-3460000

# Configure, build and install
./configure --prefix=/your/install/path/sqlite3
make -j$(nproc)
sudo make install


PROJ Installation (requires libtiff and SQLite3)

# Download and extract PROJ
wget https://download.osgeo.org/proj/proj-8.1.0.tar.gz
tar -xvf proj-8.1.0.tar.gz
cd proj-8.1.0

# Set environment variables
export TIFF_CFLAGS="-I/your/install/path/tiff/include"
export TIFF_LIBS="-L/your/install/path/tiff/lib -ltiff"
export SQLITE3_CFLAGS="-I/your/install/path/sqlite3/include"
export SQLITE3_LIBS="-L/your/install/path/sqlite3/lib -lsqlite3"

# Configure, build and install
./configure --prefix=/your/install/path/proj
make -j$(nproc)
sudo make install


2. Set Up Python Environment

Create and activate a Conda environment:
conda create -n pyenv python=3.10 -y
conda activate pyenv


3. Install Python Dependencies

Install required Python packages:
conda install pyproj -y
# Install additional dependencies if requirements.txt exists
pip install -r requirements.txt


Environment Configuration

Add these to your shell profile (~/.bashrc or ~/.zshrc):
export PATH="/your/install/path/tiff/bin:/your/install/path/sqlite3/bin:/your/install/path/proj/bin:$PATH"
export LD_LIBRARY_PATH="/your/install/path/tiff/lib:/your/install/path/sqlite3/lib:/your/install/path/proj/lib:$LD_LIBRARY_PATH"


Verification

Verify installations:
tiffinfo --version
sqlite3 --version
proj --version
python -c "import pyproj; print(pyproj.__version__)"


After completing these steps, all dependencies are installed and the environment is ready for building and running the project.

Note: Replace /your/install/path with your actual installation directory (e.g., /usr/local or $HOME/local). The -j$(nproc) flag utilizes all available CPU cores for faster compilation.
