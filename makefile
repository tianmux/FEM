tmake: ./src/main.cpp 
	dpcpp ./src/main.cpp -o FEA -qopenmp -std=c++14 -O3 -g -xCORE-AVX2  -I /opt/intel/advisor/include -L${MKLROOT}/lib/intel64 -liomp5 -lpthread -lm -ldl  -DMKL_ILP64  -qmkl=parallel -lnetcdf 

