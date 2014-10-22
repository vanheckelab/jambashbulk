#include <string>
#include <sstream>
#include <cstdlib>

#include "gitversion.h"

#ifdef _WIN32
#define PLATFORM "win32"
typedef double LDBL;
#endif
#ifdef _WIN64
#define PLATFORM "win64"
typedef double LDBL;
#endif
#ifdef __unix__
#define PLATFORM "unix"
typedef long double LDBL;
#endif

std::string header() {
  std::stringstream sstm;

  sstm << "# " __BASE_FILE__ " (" PLATFORM ") \n"
       << "# commit " GIT_HEAD " ";

  if (GIT_CHANGED_FILES) {
    sstm << "+" << GIT_CHANGED_FILES << " ";
  }

  sstm << "(" GIT_HEAD_DATE ")\n";

  srand(0);
  sstm << "# RAND_MAX = " << RAND_MAX << "\n"
       << "# srand(0); rand() = " << rand() << "\n";

  return sstm.str();
}

std::string FILE_HEADER = header();
const char * FILE_HEADER_C = FILE_HEADER.c_str();
