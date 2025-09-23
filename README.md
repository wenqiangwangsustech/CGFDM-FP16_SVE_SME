# Dependency Installation Guide

This project depends on several C/C++ libraries. Please install them in the following order before building the project.

---

## 1. Install **libtiff**

- Download  
  ```bash
  wget http://download.osgeo.org/libtiff/tiff-4.0.7.tar.gz
tar -xvf tiff-4.0.7.tar.gz
cd tiff-4.0.7
make -j
make install

wget https://www.sqlite.org/2024/sqlite-autoconf-3460000.tar.gz
