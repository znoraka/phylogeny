#include <stdio.h>
#include <jpeglib.h>
#include <iostream>

using namespace std;

JBLOCKARRAY rowPtrs[MAX_COMPONENTS];

void read(jpeg_decompress_struct srcinfo, jvirt_barray_ptr * src_coef_arrays, int coeff) {
  //cout << "Started reading DCT" << endl;
  std::cout << "(" << coeff << " (";
  for (JDIMENSION compNum=0; compNum < srcinfo.num_components; compNum++) {
    size_t blockRowSize = (size_t) sizeof(JCOEF) * DCTSIZE2 * srcinfo.comp_info[compNum].width_in_blocks;
    for (JDIMENSION rowNum=0; rowNum < srcinfo.comp_info[compNum].height_in_blocks; rowNum++) {
      // A pointer to the virtual array of dct values
      rowPtrs[compNum] = ((&srcinfo)->mem->access_virt_barray)((j_common_ptr) &srcinfo, src_coef_arrays[compNum],rowNum, (JDIMENSION) 1, FALSE);
      // Loop through the blocks to get the dct values
      for (JDIMENSION blockNum=0; blockNum < srcinfo.comp_info[compNum].width_in_blocks; blockNum++){
	//...iterate over DCT coefficients
	// for (JDIMENSION i=0; i<DCTSIZE2; i++){
	//and print them to standard out - one per line
	// std::cout <<  rowPtrs[compNum][0][blockNum][0] << " ";
	// cout << ((rowPtrs[compNum][0][blockNum][0] / coeff) - 1) * coeff << " ";
	cout << rowPtrs[compNum][0][blockNum][0] * coeff << " ";
	// cout << ") ";
      }
    }
  }
  std::cout << "))" << std::endl;
}

int main(int argc, char **argv) {
  //TODO: change this to read a different file
  const char* filename = argv[1];
  
  FILE * infile;  
  struct jpeg_decompress_struct srcinfo;
  struct jpeg_error_mgr srcerr;

  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  srcinfo.err = jpeg_std_error(&srcerr);
  jpeg_create_decompress(&srcinfo);
  jpeg_stdio_src(&srcinfo, infile);
  (void) jpeg_read_header(&srcinfo, FALSE);

  int coeff = (int) srcinfo.quant_tbl_ptrs[0]->quantval[0];

  //coefficients
  jvirt_barray_ptr * src_coef_arrays = jpeg_read_coefficients(&srcinfo);
  read(srcinfo, src_coef_arrays, coeff);

    
  jpeg_destroy_decompress(&srcinfo);
  fclose(infile);
  return 0; 
}
