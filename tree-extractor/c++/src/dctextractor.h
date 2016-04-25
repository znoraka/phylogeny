#ifndef DCT_EXTRACTOR
#define DCT_EXTRACTOR

#include <stdio.h>
#include <jpeglib.h>
#include <iostream>
#include <vector>

struct q_dct {
  int q;
  std::vector<int> dct;
  std::vector<int> dctDiffZero;
  int numberOfBlocks;
  q_dct(int q, std::vector<int> dct, std::vector<int> dctDiffZero, int n) {
    this->q = q;
    this->dct = dct;
    this->dctDiffZero = dctDiffZero;
    this->numberOfBlocks = n;
  }
};

int zigzag[64];

std::vector<int> readDCTDiffZero(jpeg_decompress_struct srcinfo, jvirt_barray_ptr * src_coef_arrays, int &numberOfBlocks, int coeff) {
  JBLOCKARRAY rowPtrs[MAX_COMPONENTS];
  
  std::vector<int> vec;
  vec.resize(64, 0);
  for (JDIMENSION compNum=0; compNum < srcinfo.num_components; compNum++) {
    size_t blockRowSize = (size_t) sizeof(JCOEF) * DCTSIZE2 * srcinfo.comp_info[compNum].width_in_blocks;
    for (JDIMENSION rowNum=0; rowNum < srcinfo.comp_info[compNum].height_in_blocks; rowNum++) {
      rowPtrs[compNum] = ((&srcinfo)->mem->access_virt_barray)((j_common_ptr) &srcinfo, src_coef_arrays[compNum],rowNum, (JDIMENSION) 1, FALSE);
      for (JDIMENSION blockNum=0; blockNum < srcinfo.comp_info[compNum].width_in_blocks; blockNum++){
	for (JDIMENSION i=0; i<DCTSIZE2; i++){
	  if(abs(rowPtrs[compNum][0][blockNum][zigzag[i]]) > 5) {
	    vec[i]++;
	  }
	}
      }
    }
  }
  numberOfBlocks = srcinfo.comp_info[0].height_in_blocks * srcinfo.comp_info[0].width_in_blocks;
  return vec;
}


std::vector<int> readDCT(jpeg_decompress_struct srcinfo, jvirt_barray_ptr * src_coef_arrays, int coeff) {
  JBLOCKARRAY rowPtrs[MAX_COMPONENTS];

  // auto avg = [](JBLOCK v) {
  //   double sum = 0;
  //   for (int i = 0; i < DCTSIZE2; i++) {
  //     sum += v[i];
  //     std::cout << "v[" << i << "] = " << v[i] << std::endl;
  //   }
    
  //   return sum / DCTSIZE2;
  // };
  
  std::vector<int> vec;
  for (JDIMENSION compNum=0; compNum < srcinfo.num_components; compNum++) {
    size_t blockRowSize = (size_t) sizeof(JCOEF) * DCTSIZE2 * srcinfo.comp_info[compNum].width_in_blocks;
    for (JDIMENSION rowNum=0; rowNum < srcinfo.comp_info[compNum].height_in_blocks; rowNum++) {
      rowPtrs[compNum] = ((&srcinfo)->mem->access_virt_barray)((j_common_ptr) &srcinfo, src_coef_arrays[compNum],rowNum, (JDIMENSION) 1, FALSE);
      for (JDIMENSION blockNum=0; blockNum < srcinfo.comp_info[compNum].width_in_blocks; blockNum++){
	for (JDIMENSION i=0; i<DCTSIZE2; i++){
	//and print them to standard out - one per line
	// std::cout <<  rowPtrs[compNum][0][blockNum][0] << " ";
	// cout << ((rowPtrs[compNum][0][blockNum][0] / coeff) - 1) * coeff << " ";
	// cout << rowPtrs[compNum][0][blockNum][0] * coeff << " ";
	// cout << ") ";
	/* vec.push_back(rowPtrs[compNum][0][blockNum][0] * coeff); */
	  // std::cout << avg(rowPtrs[compNum][0][blockNum]) << std::endl;
	  vec.push_back(rowPtrs[compNum][0][blockNum][zigzag[0]]);
	}
      }
    }
  }
  return vec;
}


q_dct getQAndDct(std::string path, int quality = 100){
  FILE * infile;  
  struct jpeg_decompress_struct srcinfo;
  struct jpeg_error_mgr srcerr;

  int n = 0;
  for (int i = 0; i < 8 * 2; i++)
    for (int j = (i < 8) ? 0 : i-8+1; j <= i && j < 8; j++)
      zigzag[n++] = (i&1)? j*(8-1)+i : (i-j)*8+j;


  if ((infile = fopen(path.c_str(), "rb")) == NULL) {
    // fprintf(stderr, "can't open %s\n", path);
    std::cerr << "can't open " << path << std::endl;
    // return q_dct;
  }

  srcinfo.err = jpeg_std_error(&srcerr);
  jpeg_create_decompress(&srcinfo);
  jpeg_stdio_src(&srcinfo, infile);
  (void) jpeg_read_header(&srcinfo, FALSE);

  int coeff = (int) srcinfo.quant_tbl_ptrs[0]->quantval[1];

  //coefficients
  jvirt_barray_ptr * src_coef_arrays = jpeg_read_coefficients(&srcinfo);
  int numberOfBlocks = 0;
  auto dct = readDCT(srcinfo, src_coef_arrays, coeff);
  auto dctDiffZero = readDCTDiffZero(srcinfo, src_coef_arrays, numberOfBlocks, coeff);
  q_dct q(coeff, dct, dctDiffZero, numberOfBlocks);
      
  jpeg_destroy_decompress(&srcinfo);
  fclose(infile);
  return q;
}

#endif
