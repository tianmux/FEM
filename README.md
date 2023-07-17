# FEM
prerequisite
1. Intel MKL.
   The recommended way of installing the required library:
   https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html
   DO NOT FORGET TO ADD LD_LIBRARY_PATH
   add the following line in ~/.bashrc
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/MKL
3. NetCDF4.
   Installation guide:https://docs.unidata.ucar.edu/nug/current/getting_and_building_netcdf.html
   a. Download the source from unidata:
   https://downloads.unidata.ucar.edu/netcdf
   b. download the needed package hdf5:
   https://www.hdfgroup.org/downloads/hdf5
   c. download the needed package zlib:
   https://zlib.net/
   d. Install zlib first:
    $ # Build and install zlib
    $ ZDIR=/usr/local
    $ ./configure --prefix=${ZDIR}
    $ make check
    $ make install   # or sudo make install, if root permissions required
   e. Install hdf5:
    $ # Build and install HDF5
    $ H5DIR=/usr/local
    $ ./configure --with-zlib=${ZDIR} --prefix=${H5DIR} --enable-hl
    $ make check
    $ make install   # or sudo make install, if root permissions required
   f. Install netCDF
    $ # Build and install netCDF-4
    $ NCDIR=/usr/local
    $ CPPFLAGS='-I${H5DIR}/include -I${ZDIR}/include' LDFLAGS='-L${H5DIR}/lib -L${ZDIR}/lib' ./configure --prefix=${NCDIR}
    $ make check
    $ make install  # or sudo make install
   g. Add the library directory to the LD_LIBRARY_PATH
    add the following line in ~/.bashrc
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/ncdf
3. Add the relevant paths to .bashrc for the use of makegcc
   export NCDIR=/path/to/ncdf
   export MKLDIR=/path/to/MKL
