#ifndef DCT_EXTRACTOR
#define DCT_EXTRACTOR

#include <stdio.h>
#include <jpeglib.h>
#include <iostream>
#include <vector>

struct q_dct {
  int q;
  std::vector<std::vector<int> > dct;
  std::vector<int> dctDiffZero;
  int numberOfBlocks;
  q_dct(int q, std::vector<std::vector<int> > dct, std::vector<int> dctDiffZero, int n) {
    this->q = q;
    this->dct = dct;
    this->dctDiffZero = dctDiffZero;
    this->numberOfBlocks = n;
  }
};

int zigzag[64];

std::vector<std::vector<int> > readDCT(jpeg_decompress_struct srcinfo, jvirt_barray_ptr * src_coef_arrays) {
  JBLOCKARRAY rowPtrs[MAX_COMPONENTS];
  
  std::vector<std::vector<int> > vec;
  for (int i = 0; i < DCTSIZE2; i++) {
    vec.push_back(std::vector<int>());
  }

  for (JDIMENSION compNum=0; compNum < srcinfo.num_components; compNum++) {
    size_t blockRowSize = (size_t) sizeof(JCOEF) * DCTSIZE2 * srcinfo.comp_info[compNum].width_in_blocks;
    for (JDIMENSION rowNum=0; rowNum < srcinfo.comp_info[compNum].height_in_blocks; rowNum++) {
      rowPtrs[compNum] = ((&srcinfo)->mem->access_virt_barray)((j_common_ptr) &srcinfo, src_coef_arrays[compNum],rowNum, (JDIMENSION) 1, FALSE);
      for (JDIMENSION blockNum=0; blockNum < srcinfo.comp_info[compNum].width_in_blocks; blockNum++){
	for (JDIMENSION i=0; i<DCTSIZE2; i++){
	  vec[i].push_back(rowPtrs[compNum][0][blockNum][zigzag[i]]);
	}
      }
    }
  }
  return vec;
}


std::vector<std::vector<int> > getDctCoeffs(std::string path) {
  FILE * infile;  
  struct jpeg_decompress_struct srcinfo;
  struct jpeg_error_mgr srcerr;

  if ((infile = fopen(path.c_str(), "rb")) == NULL) {
    std::cerr << "can't open " << path << std::endl;
  }

  srcinfo.err = jpeg_std_error(&srcerr);
  jpeg_create_decompress(&srcinfo);
  jpeg_stdio_src(&srcinfo, infile);
  (void) jpeg_read_header(&srcinfo, FALSE);

  jvirt_barray_ptr * src_coef_arrays = jpeg_read_coefficients(&srcinfo);

  // for (int i = 0; i < 64; i++) {
  //   std::cout << (int) srcinfo.quant_tbl_ptrs[0]->quantval[zigzag[i]] << "\n";
  // }

  return readDCT(srcinfo, src_coef_arrays);
}

#endif
