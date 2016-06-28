#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <sys/stat.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <list>

#include "dctextractor.h"
#include "tools.h"
#include "image_ppm.h"

using namespace tools;

struct Image {
  std::vector<OCTET> data;
  int width;
  int height;

  Image(int width, int height) {
    this->width = width;
    this->height = height;
    data.resize(width * height);
  }

  Image(std::string path) {
    lire_nb_lignes_colonnes_image_pgm(const_cast<char*>(path.c_str()), &height, &width);
    data.resize(height * width);
    lire_image_pgm(const_cast<char*>(path.c_str()), &data[0], data.size());
  }
};

static std::vector<int> baseTable = {16,12,14,14,18,24,49,72,11,12,13,17,22,35,64,92,10,14,16,22,37,55,78,95,16,19,24,29,56,64,87,98,24,26,40,51,68,81,103,112,40,58,57,87,109,104,121,100,51,60,69,80,103,113,120,103,61,55,56,62,77,92,101,99};

static std::vector<std::vector<int> > tables;
double meanError;
int overestimate;
int underestimate;
std::string tmpDirectory = "./";

std::vector<int> table(int quality) {
  auto getQ = [&](int val) -> int {
    if(val < 50) {
      return 5000 / val;
    }
    return (200 - (val * 2));
  };

  int q = getQ(quality);
  
  std::vector<int> v(baseTable);
  for(auto &i : v) {
    i = fmax(fmin((i * q + 50) / 100, 255), 1);
  }

  return v;
}

int estimateQ(std::vector<std::vector<int> > dctCoeffs);
std::vector<int> makeHisto (std::vector<int> v);
bool isParent(Image img1, Image img2);
void compressPGMImage(Image img, int quality, std::string path);

/**
 * Estimation du facteur de qualité d'une image JPEG à partir de ses coefficents DCT
 */
int estimateQ(std::vector<std::vector<int> > dctCoeffs) {
  int peakFactor = 5;
  int peakCount = 3;
  int dctCount = 35; //limiter les coefficients dans l'ordre zigzag permet d'avoir une meilleur précision

  auto isPeak = [](std::vector<double> vec, int width, int index) {
    for (int i = index - width; i < index; i++) {
      if(vec[i] >= vec[i + 1]) return false;
    }

    for (int i = index; i < index + width; i++) {
      if(vec[i] <= vec[i + 1]) return false;
    }
    return true;
  };

  auto computeAutocorrelation = [](std::vector<double> histo) {
    auto metric = Distance::autocorrelation;
    std::vector<double> plot2;

    for (int i = 1; i < fmin(255, histo.size()); i++) {
      double nextDistance = Distance::compute(metric, histo.begin(), histo.end(), histo.begin() + i, histo.end());
      plot2.push_back(nextDistance);
    }
    
    return plot2;
  };

  auto findPeaks = [&](int width, std::vector<double> plot){
    std::vector<int> peaks;
    double plotAvg = tools::mean(plot);
    
    for (int i = width; plot.size() > 1 && i < plot.size() - width && peaks.size() < peakCount; i++) {
      if(plot[i] > plotAvg * peakFactor && isPeak(plot, 1, i)) {
	peaks.push_back(i);
      }
    }
    return peaks;
  };

  auto periodFromPeaks = [](std::vector<int> peaks) {
    for (int i = peaks.size() - 1; i > 0; i--) {
      peaks[i] = peaks[i] - peaks[i - 1];
    }
    peaks[0] = peaks[1];

    double vecStdDev = tools::standardDeviation(peaks);
    
    int sum = 0;
    int t = 0;

    for(auto i : peaks) {
      if(i > vecStdDev) {
	sum += i;
	t++;
      }
    }
    return (int) (sum / t);
  };

  auto primaryEstimation = [](std::vector<int> periods) {
    int index = 0;
    double minValue = std::numeric_limits<double>::max();
  
    for (int i = 1; i < 100; i++) {
      std::vector<double> tabledouble;
      for (int j = 0; j < tables[i].size(); j++) {
  	tabledouble.push_back(tables[i][zigzag[j]]);
      }

      std::vector<double> qsdouble;

      for (int i = 0; i < periods.size(); i++) {
  	if(periods[i] == -1) qsdouble.push_back(1);
  	else if(periods[i] <= 1) qsdouble.push_back(tabledouble[i]);
  	else qsdouble.push_back(periods[i]);
      }

      double d = Distance::compute(Distance::euclid, qsdouble, tabledouble);
      if(d < minValue) {
  	minValue = d;
  	index = i;
      }
    }

    return index;
  };

  auto secondaryEstimation = [&](std::vector<int> periods) {
    int index = primaryEstimation(periods);
    if(index < 50)
      return index;

    auto qFromPeriod = [&](double period, int i) {
      double q = (100.0 * period - 50.0) / (double)baseTable[i];
      return floor((200.0 - q) / 2.0);
    };

    std::vector<double> estimatedQs;

    for (int i = 0; i < periods.size(); i++) {
      estimatedQs.push_back(qFromPeriod(periods[i], zigzag[i]));
    }

    int estimated = tools::mean(estimatedQs);

    if(index < 92 && abs(index - estimated) >= 3) return index;
      
    return estimated;
  };
  
  std::vector<int> qs;
  int n = 0;
  for (int k = 0; k < dctCount; k++) {
    auto dct = dctCoeffs[k];
    auto histo = makeHisto(dct);
    std::vector<double> vec(histo.begin(), histo.end());

    std::vector<double> plot = computeAutocorrelation(vec);
    std::vector<int> peaks = findPeaks(2, plot);

    if(peaks.size() == 0) {
      qs.push_back(1);
      continue; 
    } else if (peaks.size() == 1) {
      qs.push_back(peaks[0] + 1);
      continue;
    }
    
    qs.push_back(periodFromPeaks(peaks));
  }

  int p = secondaryEstimation(qs);

  if(p < 5) p = 97;

  return p;
}

/**
 * Transforme un vecteur d'int en histogramme
 */
std::vector<int> makeHisto (std::vector<int> v) {
  int maxElem = *std::max_element(v.begin(), v.end());
  maxElem = fmax(maxElem, abs(*std::min_element(v.begin(), v.end())));
  std::vector<int> out;
  out.resize(maxElem+1, 0);

  for(auto i : v) {
    out[abs(i)]++;
  }
    
  return out;
}

/**
 * Crée une Image à partir d'une image JPEG
 */
Image readJPEGImage(std::string path) {
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW row_pointer[1];

  FILE *infile = fopen(path.c_str(), "rb");
  unsigned long location = 0;
  int i = 0;

  if ( !infile ) {
    printf("Error opening jpeg file %s\n!", path.c_str() );
    return Image(0, 0);
  }
  cinfo.err = jpeg_std_error( &jerr );
  jpeg_create_decompress( &cinfo );
  jpeg_stdio_src( &cinfo, infile );
  jpeg_read_header( &cinfo, TRUE );
  jpeg_start_decompress( &cinfo );

  row_pointer[0] = (unsigned char *)malloc( cinfo.output_width*cinfo.num_components );
  Image img(cinfo.image_width, cinfo.image_height);
  
  while( cinfo.output_scanline < cinfo.image_height ) {
    jpeg_read_scanlines( &cinfo, row_pointer, 1 );
    for( i=0; i<cinfo.image_width*cinfo.num_components;i++)
      img.data[location++] = row_pointer[0][i];
    }
  jpeg_finish_decompress( &cinfo );
  jpeg_destroy_decompress( &cinfo );
  free( row_pointer[0] );
  fclose( infile );
  return img;
}


/**
 * Compresse une bitmap en jpeg
 */
void compressPGMImage(Image img, int quality, std::string path = "./tmp.jpg") {
  FILE* outfile = fopen(path.c_str(), "wb");
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr       jerr;
 
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, outfile);

  cinfo.image_width      = img.width;
  cinfo.image_height     = img.height;
  cinfo.input_components = 1;
  cinfo.in_color_space   = JCS_GRAYSCALE;

  jpeg_set_defaults(&cinfo);
  jpeg_set_quality (&cinfo, quality, true);
  jpeg_start_compress(&cinfo, true);

  int n = 0;
  JSAMPROW row_pointer;  
  
  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer = (JSAMPROW)&img.data[n];
    n += img.height;
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

}

/**
 * Retourne vrai si img1 est le parent de img2
 */			  
bool isParent(Image img1, Image img2) {
  auto makeDistrib = [](std::vector<int> v) {
    std::vector<double> out;
    long sum = 0;

    for(auto i : v) {
      sum += i;
    }
    
    for(auto i : v) {
      double di = i;
      double d = di / (double) sum;
      out.push_back(d);
    }
    return out;
  };

  auto distancesOk = [&](std::vector<std::vector<int> > d1, std::vector<std::vector<int> > d2) {
    auto metric = Distance::kb;
    double sumd = 0;
    double sumdd = 0;

    for (int i = 0; i < d1.size(); i++) {
      double d = Distance::compute(metric, makeDistrib(makeHisto(d1[i])), makeDistrib(makeHisto(d2[i])));
      double dd = Distance::compute(metric, makeDistrib(makeHisto(d2[i])), makeDistrib(makeHisto(d1[i])));

      sumd += d;
      sumdd += dd;
    }

    // std::cout << sumd << " -- " << sumdd << " -- " << fabs(sumd - sumdd) << std::endl;
    return (sumd < 0.6 && sumdd < 0.6 && fabs(sumd - sumdd) < 0.01);
  };

  compressPGMImage(img2, 100, tmpDirectory + "img2.jpg");
  auto coeffs2 = getDctCoeffs(tmpDirectory + "img2.jpg");
  int estimated = estimateQ(coeffs2);

  // std::cout << "estimated Qf for image 2 = " << estimated << "\n";
  
  int step = 0;
  int range = 100;
  int k = 0;
  while(k < range) {
    int cptj = -1;
    step = step + ((k % 2 == 0)? -1 : 1) * k++;

    if(estimated + step > 100 || estimated + step < 0) continue;

    compressPGMImage(img1, estimated + step, tmpDirectory + "img1.jpg");
    // std::cout << "testing Qf = " << estimated + step << " :: ";
    Image img = readJPEGImage(tmpDirectory + "img1.jpg");

    compressPGMImage(img, 100, tmpDirectory + "img1.jpg");    
    
    auto coeffs1 = getDctCoeffs(tmpDirectory + "img1.jpg");
    bool parent = distancesOk(coeffs1, coeffs2);
    if(parent) return true;
  }  
  return false;
}

void displayHelp() {
  std::cout << "Tells whether an image is the parent of the other" << "\n\n";
  std::cout << "Usage: ancestor-estimator image1.pgm image2.pgm [tmp directory]" << "\n";
  std::cout << "\n  Images must be pgm files" << "\n";
  std::cout << "\n  tmp directory is the directory where temporary files will be written, default is ./, a ramdisk might ease your hard drive and speed up the process" << "\n";
}

int main(int argc, char **argv) {
  int n = 0;
  for (int i = 0; i < 8 * 2; i++)
    for (int j = (i < 8) ? 0 : i-8+1; j <= i && j < 8; j++)
      zigzag[n++] = (i&1)? j*(8-1)+i : (i-j)*8+j;


  for (int i = 1; i < 101; i++) {
    tables.push_back(table(i));
  }

  if(argc < 3) {
    displayHelp();
    return -1;
  }

  if(argc == 4) {
    tmpDirectory = argv[3];
    tmpDirectory += "/";
  }
  
  auto f = [&](std::string img1, std::string img2) {
    Image image1(img1);
    Image image2(img2);

    std::string s = "";
    std::cout << "\n";
    bool b = isParent(image1, image2);
    if(!b) {
      s = "not ";
    }
    std::cout << img1 << " is " << s <<  "the parent of " << img2 << "\n";
    return b;
  };

  if(f(argv[1], argv[2])) return 0;
  f(argv[2], argv[1]);
    
  return 0;
}
