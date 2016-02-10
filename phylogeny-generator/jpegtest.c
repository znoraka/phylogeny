#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <ctime>
/* #include <Magick++.h> */
extern "C" {
#include <jpeglib.h>
}

int main(int argc, char **argv) {
  struct jpeg_decompress_struct info;
  struct jpeg_error_mgr err;
  info.err = jpeg_std_error(&err);     

}
