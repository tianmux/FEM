// includes, system
#ifndef MY_HEADER_FOR_FEA_H
#define MY_HEADER_FOR_FEA_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <map>
#include <time.h>
#include <omp.h>
#include <mkl.h>
#ifdef _WIN32
#  define WINDOWS_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#endif
#define _USE_MATH_DEFINES
#include <cmath>

// OpenGL Graphics includes
//#include <GL/glew.h>
#if defined (__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#ifndef glutCloseFunc
#define glutCloseFunc glutWMCloseFunc
#endif
#else
//#include <GL/freeglut.h>
#endif

#endif