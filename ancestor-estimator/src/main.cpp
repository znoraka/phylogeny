#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <sys/stat.h>
#include <assert.h>
#include <Magick++.h>
#include <Magick++/STL.h>
#include <magick/MagickCore.h>
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

// static std::vector<int> baseTable = {11,12,14,12,10,16,14,13,14,18,17,16,19,24,40,26,24,22,22,24,49,35,37,29,40,58,51,61,60,57,51,56,55,64,72,92,78,64,68,87,69,55,56,80,109,81,87,95,98,103,104,103,62,77,113,121,112,100,120,92,101,103,99};

static std::vector<std::vector<int> > tables;
double meanError;
int overestimate;
int underestimate;

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

std::vector<std::vector<bool> > estimateParents(std::string directory);
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
    // std::cout << "rough = " << index << std::endl;
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
 * Transforme un vecteur d'int en histogramme de ces valeurs
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

void displayTable(int i) {
  int cpt = 0;
  for(auto v : table(i)) {
    // std::cout << v << "&";
    if(cpt++ == 7) {
      cpt = 0;
      // std::cout << "\\\\ \\hline" << std::endl;
    }
  }
}

Image readJPEGImage(std::string path) {
  /* these are standard libjpeg structures for reading(decompression) */
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  /* libjpeg data structure for storing one row, that is, scanline of an image */
  JSAMPROW row_pointer[1];

  FILE *infile = fopen(path.c_str(), "rb");
  unsigned long location = 0;
  int i = 0;

  if ( !infile ) {
    printf("Error opening jpeg file %s\n!", path.c_str() );
    // return -1;
    return Image(0, 0);
  }
  /* here we set up the standard libjpeg error handler */
  cinfo.err = jpeg_std_error( &jerr );
  /* setup decompression process and source, then read JPEG header */
  jpeg_create_decompress( &cinfo );
  /* this makes the library read from infile */
  jpeg_stdio_src( &cinfo, infile );
  /* reading the image header which contains image information */
  jpeg_read_header( &cinfo, TRUE );
  /* Uncomment the following to output image information, if needed. */
  /*--
    printf( "JPEG File Information: \n" );
    printf( "Image width and height: %d pixels and %d pixels.\n", cinfo.image_width, cinfo.image_height );
    printf( "Color components per pixel: %d.\n", cinfo.num_components );
    printf( "Color space: %d.\n", cinfo.jpeg_color_space );
    --*/
  /* Start decompression jpeg here */
  jpeg_start_decompress( &cinfo );

  // /* allocate memory to hold the uncompressed image */
  // raw_image = (unsigned char*)malloc( cinfo.output_width*cinfo.output_height*cinfo.num_components );
  // /* now actually read the jpeg into the raw buffer */
  row_pointer[0] = (unsigned char *)malloc( cinfo.output_width*cinfo.num_components );
  /* read one scan line at a time */
  Image img(cinfo.image_width, cinfo.image_height);
  
  while( cinfo.output_scanline < cinfo.image_height ) {
    jpeg_read_scanlines( &cinfo, row_pointer, 1 );
    for( i=0; i<cinfo.image_width*cinfo.num_components;i++)
      // raw_image[location++] = row_pointer[0][i];
      img.data[location++] = row_pointer[0][i];
    }
  /* wrap up decompression, destroy objects, free pointers and close open files */
  jpeg_finish_decompress( &cinfo );
  jpeg_destroy_decompress( &cinfo );
  free( row_pointer[0] );
  fclose( infile );
  /* yup, we succeeded! */
  return img;
}


/**
 * Compresse une bitmap en jpeg
 */
void compressPGMImage(Image img, int quality, std::string path = "tmp.jpg") {
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

    std::cout << sumd << " -- " << sumdd << " -- " << fabs(sumd - sumdd) << std::endl;
    return (sumd < 0.6 && sumdd < 0.6 && fabs(sumd - sumdd) < 0.01);
    // return sumd == sumdd;
  };

  // ecrire_image_pgm("tmp.pgm", &img2.data[0],  img2.height, img2.width);

  // Magick::Image image;
  // image.read("tmp.pgm");
  // image.quality(100);
  // image.write("img2.jpg");

  compressPGMImage(img2, 100, "img2.jpg");
  auto coeffs2 = getDctCoeffs("img2.jpg");
  int estimated = estimateQ(coeffs2);

  std::cout << "estimated Qf for image 2 = " << estimated << "\n";

  // ecrire_image_pgm("tmp.pgm", &img1.data[0],  img1.height, img1.width);

  // compressPGMImage(img1, 50, "img1.jpg");
  // auto coeffs3 = getDctCoeffs("img1.jpg");
  
  int step = 0;
  int range = 100;
  int k = 0;
  while(k < range) {
    int cptj = -1;
    step = step + ((k % 2 == 0)? -1 : 1) * k++;

    compressPGMImage(img1, estimated + step, "img1.jpg");
    // image.read("tmp.pgm");
    // image.quality(estimated + step);
    std::cout << "testing Qf = " << estimated + step << " :: ";
    // image.write("img1.jpg");
    // image.read("img1.jpg");
    // image.quality(100);
    // image.write("img1.jpg");

    Image img = readJPEGImage("img1.jpg");
    ecrire_image_pgm("tmp1.pgm", &img.data[0],  img.height, img.width);

    compressPGMImage(img, 100, "img1.jpg");    
    
    auto coeffs1 = getDctCoeffs("img1.jpg");
    bool parent = distancesOk(coeffs1, coeffs2);
    if(parent) return true;
  }  
  return false;
}


int main(int argc, char **argv) {
  Magick::InitializeMagick(*argv);

  int n = 0;
  for (int i = 0; i < 8 * 2; i++)
    for (int j = (i < 8) ? 0 : i-8+1; j <= i && j < 8; j++)
      zigzag[n++] = (i&1)? j*(8-1)+i : (i-j)*8+j;


  for (int i = 1; i < 101; i++) {
    tables.push_back(table(i));
  }

  Image image1(argv[1]);
  Image image2(argv[2]);

  std::cout << "parents = " << isParent(image1, image2) << std::endl;
  
  return 0;
}
