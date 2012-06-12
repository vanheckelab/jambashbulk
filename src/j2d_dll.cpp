// g++ -shared -Wl,-soname,j2d_dll.so -o j2d_dll.so j2d_dll.cpp

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


extern "C" {
  #include "jamBashbulk.cpp"
  vector <long double> *p_ptr = &p;

  int vector_as_array(vector <long double> *ptr, long double **out_ptr) {
    *out_ptr = &((*ptr)[0]);
    return ptr->capacity();
  }
}
