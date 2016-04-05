#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <ctime>
#include <Magick++.h>
#include <Magick++/STL.h>
#include <magick/MagickCore.h>
#include <map>

#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>

extern "C" {
#include <jpeglib.h>
}

#include "gnuplot_i.hpp"
#include "CImg.h"

#define cimg_use_jpeg

using namespace Magick;
using namespace cimg_library; 
using namespace cv;

//parcours zigzag
int zigzag[64];

struct Node {
  std::string name;
  float entropy;
  float entropy2;
  int compression;
  std::vector<Node*> children;
};

struct ljpegimg {
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr err;
  jvirt_barray_ptr *coeffs;
  FILE *fp;

  ljpegimg(std::string path) {
    cinfo.err = jpeg_std_error(&err);     
    jpeg_create_decompress(&cinfo);
    if((fp = fopen(path.c_str(), "rb")) == NULL) {
      std::cout << "could not open image : " << path << std::endl;
      return;
    }
    jpeg_stdio_src(&cinfo, fp);
    (void) jpeg_read_header(&cinfo, true);
    coeffs = jpeg_read_coefficients(&cinfo);
  }

  int height() {
    return cinfo.comp_info->height_in_blocks;
  }

  int width() {
    return cinfo.comp_info->width_in_blocks;
  }

  int dctCoeff(int i, int j, int dctIndex) {
    return cinfo.mem->access_virt_barray((j_common_ptr)&cinfo, coeffs[0], i, 1, false)[0][j][zigzag[dctIndex]];
  }
};

void createChild(std::string path, int imageCount, std::vector<Node*> &nodes);
Node *createPhylogeny(std::string path, std::string rootImage, int imageCount);
void applyTransform(Image &image, Node *node);
void exportPhylogeny(std::ostream& stream, Node *root);
void plotDctCoefficients(std::ostream &stream, std::string path);
void plotDctCoefficientsMultipleImages(std::vector<std::string> paths);
std::vector<std::vector<float> > getDctHistograms(std::string path, int count);
float computeEntropy(Image &child, Image &parent);

std::vector<std::vector<float> > getDctHistograms(std::string path, int count) {
  ljpegimg img(path);
  std::vector<std::vector<float> > vec;

  for (int dctIndex = 0; dctIndex < count; dctIndex++) {
      std::map<int, int> histo;
      for (int i = 0; i < img.height(); i++) {
	for (int j = 0; j < img.width(); j++) {
	  int dct = img.dctCoeff(i, j, dctIndex);
	  histo[dct]++;
	}
      }

      int maxVal = -999999;
      int minVal = 999999;
      for(auto i : histo) {
	maxVal = fmax(maxVal, i.first);
	minVal = fmin(minVal, i.first);
      }

      int range = fmax(maxVal, abs(minVal));
      std::vector<float> h(2 * range + 1);

      for(auto i : histo) {
	h[i.first + range] = (float) (i.second) / (float) (img.width() * img.height());
      }
      vec.push_back(h);
  }
  return vec;
}

void plotDctCoefficientsMultipleImages(std::vector<std::string> paths) {
  Gnuplot g1("hello");

  std::vector<std::vector<std::vector<float> > > histos;
  
  for(auto path : paths) {
    histos.push_back(getDctHistograms(path, 16));
  }

  auto f = [](std::vector<float> v1, std::vector<float> v2) {
    int diff = (v1.size() - v2.size()) * 0.5;

    std::vector<float> v(v1.size(), 0);

    for (int i = 0; i < v2.size(); i++) {
      v[i + diff] = v2[i];
    }
    return v;
  };

  auto h1 = histos[0][0];
  auto h2 = histos[2][0];

  if(h1.size() < h2.size()) {
    h1 = f(h2, h1);
  } else {
    h2 = f(h1, h2);
  }
  
  cv::Mat histo1(h1);
  cv::Mat histo2(h2);

  std::cout << compareHist(histo1, histo2, CV_COMP_CORREL) << std::endl;
  
  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < paths.size(); j++) {
      auto v = histos[j][i];
      int range = v.size() * 0.5;
      g1.cmd("set xrange [" + std::to_string(-range) + ":" + std::to_string(range) + "]");
      g1.cmd("set style data histogram");
      g1.cmd("set style histogram clustered");
      g1.cmd("set style fill solid border");
      g1.cmd("set term png");
      std::string s;
      s += std::to_string(i) + ".png";

      g1.set_style("boxes");

      std::vector<int> x;
      std::vector<float> y;

      for (int a = 0; a < v.size(); a++) {
	  x.push_back(a - range);
	  y.push_back(v[a]);
	  // std::cout << v[a] << std::endl;
      }
      
      g1.cmd("set output \"" + s + "\"");
      g1.plot_xy(x, y, "");
      // g1.cmd("replot");
    }
    g1.reset_plot();
  }
}



float computeEntropy(Image &child, Image &parent) {
  double entropy = 0;
  
  MagickCore::ExceptionInfo *exceptionInfo = MagickCore::AcquireExceptionInfo();;
  MagickCore::CompareImages(parent.image(), child.image(), MagickCore::PerceptualHashErrorMetric, &entropy, exceptionInfo);
  
  // CImg<unsigned char> img(path.c_str());
  // // CImg<unsigned char> img("out.jpg");
  // std::vector<float> entropies;

  // for (int i = 0; i < img.height() - 8; i+=8) {
  //   for (int j = 0; j < img.width() - 8; j+=8) {
  //     std::map<int, int> histo;
  //     for (int x = 0; x < 8 + 0; x++) {
  // 	for (int y = 0; y < 8 + 0; y++) {
  // 	  histo[img.atXYZ(x + i, y + j, 0)]++;
  // 	}
  //     }
  //     float e = 0;
  //     for (auto p : histo) {
  // 	double freq = static_cast<float>( p.second ) / 64 ;
  // 	e += freq * log2( freq ) ;
  //     }
  //     e *= -1 ;
  //     entropies.push_back(e);
  //   }
  // }

  // double sum = std::accumulate(entropies.begin(), entropies.end(), 0.0);
  // double mean = sum / entropies.size();
  // entropy = mean;
  
  // std::map<Color,unsigned long> histogram;
  // colorHistogram( &histogram, image );
  // int size = image.rows() * image.columns();
  // for (auto p : histogram) {
  //   double freq = static_cast<float>( p.second ) / size ;
  //   entropy += freq * log2( freq ) ;
  // }
  // entropy *= -1 ;

  // struct jpeg_decompress_struct info;
  // struct jpeg_error_mgr err;
  // jvirt_barray_ptr *coeffs;
  
  // info.err = jpeg_std_error(&err);     
  // jpeg_create_decompress(&info);

  // FILE *fp;
  // if((fp = fopen("out.jpg", "rb")) == NULL) {
  //   std::cout << "could not open image : " << path << std::endl;
  //   return -1;
  // }
  
  // jpeg_stdio_src(&info, fp);
  // (void) jpeg_read_header(&info, true);
  // coeffs = jpeg_read_coefficients(&info);

  // std::map<int, int> histo;

  // for (int channel = 0; channel < 1; channel++) {
  //   for (int dctIndex = 0; dctIndex < 1; dctIndex++) {
  //     for (int i = 0; i < info.comp_info->height_in_blocks; i++) {
  // 	for (int j = 0; j < info.comp_info->width_in_blocks; j++) {

  // 	  //aller voir jpeg_snoop
  // 	  int n = abs(info.mem->access_virt_barray((j_common_ptr)&info, coeffs[0], i, 1, false)[0][j][zigzag[dctIndex]]);
  // 	  histo[n]++;
  // 	}
  //     }
  //   }
  // }

  // for (auto p : histo) {
  //   double freq = static_cast<float>( p.second ) / histo.size() ;
  //   entropy += freq * log2( freq ) ;
  // }
  // entropy *= -1;
  // fclose(fp);

  return entropy;
}

void plotDctCoefficients(std::ostream &stream, std::string path) {
  struct jpeg_decompress_struct info;
  struct jpeg_error_mgr err;
  jvirt_barray_ptr *coeffs;

  Gnuplot g1("lines");
  
  info.err = jpeg_std_error(&err);     
  jpeg_create_decompress(&info);

  FILE *fp;
  if((fp = fopen(path.c_str(), "rb")) == NULL) {
    std::cout << "could not open image : " << path << std::endl;
    return;
  }
  
  jpeg_stdio_src(&info, fp);
  (void) jpeg_read_header(&info, true);
  coeffs = jpeg_read_coefficients(&info);

  for (int channel = 0; channel < 1; channel++) {
    // for (int dctIndex = 0; dctIndex < 1; dctIndex++) {
    for (int dctIndex = 0; dctIndex < 16; dctIndex++) {
      std::map<int, int> histo;
      for (int i = 0; i < info.comp_info->height_in_blocks; i++) {
	for (int j = 0; j < info.comp_info->width_in_blocks; j++) {

	  //aller voir jpeg_snoop
	  int n = info.mem->access_virt_barray((j_common_ptr)&info, coeffs[0], i, 1, false)[0][j][zigzag[dctIndex]];
	  histo[n]++;
	}
      }

      std::vector<int> x;
      std::vector<float> y;

      int range = 0; // to center histogram
      
      for(auto i : histo) {
	x.push_back(i.first);
	if(abs(i.first) > range) range = abs(i.first);
	y.push_back((float)i.second / (info.comp_info->height_in_blocks * info.comp_info->width_in_blocks));
      }

      // montage -geometry 640x480+0+0 {0..15}.png out.png
      
      g1.cmd("set xrange [" + std::to_string(-range) + ":" + std::to_string(range) + "]");
      g1.cmd("set style data histogram");
      g1.cmd("set style histogram clustered");
      g1.cmd("set style fill solid border");
      g1.cmd("set term png");
      std::string s;
      s += std::to_string(dctIndex) + ".png";

      g1.set_style("boxes");
      
      g1.cmd("set output \"" + s + "\"");
      g1.plot_xy(x, y, "");
      g1.cmd("replot");
      g1.reset_plot();
    }
  }
  // for (int ci = 0; ci < 1; ci++)
  //   {
  //     JBLOCKARRAY buffer_one;
  //     JCOEFPTR blockptr_one;
  //     jpeg_component_info* compptr_one;
  //     compptr_one = info.comp_info + ci;

  //     for (int by = 0; by < 1; by++)
  //     // for (int by = 0; by < compptr_one->height_in_blocks; by++)
  //       {
  // 	  buffer_one = (info.mem->access_virt_barray) ((j_common_ptr)&info, src_coef_arrays[ci], by, (JDIMENSION)1, FALSE);
  // 	  for (int bx = 0; bx < 1; bx++)
  // 	  // for (int bx = 0; bx < compptr_one->width_in_blocks; bx++)
  //           {
  // 	      blockptr_one = buffer_one[0][bx];
  // 	      // QVector<int> tmp;
  // 	      for (int bi = 0; bi < 64; bi++)
  //               {
  // 		  std::cout << blockptr_one[bi] << std::endl;
  // 	      	  // tmp.append(blockptr_one[bi]);
  //               }
  // 	      // dct_coeff.push_back(tmp);
  //           }
  //       }
  //   }  
}

void applyTransform(Image &image, Node *node) {
  std::vector<std::function<void(Image &image)> > functions;

  //resize
  // functions.push_back([&](Image &image) {
  //     int cols = image.columns() * (0.9 + (std::rand() % 100) * 0.002);
  //     int rows = image.rows() * (0.9 + (std::rand() % 100) * 0.002);
  //     std::string s = std::to_string(cols) + "x" + std::to_string(rows) + "!";      
  //     image.resize(s);
  //   });

  //recompress
  functions.push_back([&](Image &image) {
      int q = 50 + std::rand() % 50;
      image.quality(q);
      node->compression = q;
    });

  //gamma adjustment
  // functions.push_back([&](Image &image) {
  //     image.gamma(0.3 + (std::rand() % 200) * 0.01);
  //   });

  //brightness adjustment
  // functions.push_back([&](Image &image) {
  //     image.modulate((std::rand() % 100) * 0.4 + 80, 100, 100);
  //   });

  //rotation
  // functions.push_back([&](Image &image) {
  //     image.rotate((std::rand() % 10) - 5);
  //   });


  //contrast adjustment
  // functions.push_back([&](Image &image) {
  //     image.contrast(100);
  //   });


  functions[std::rand() % functions.size()](image);
}

Node *createPhylogeny(std::string path, std::string rootImage, int imageCount) {
  Image image;
  int count = 0;
  std::vector<Node*> nodes;

  std::string rootPath = std::to_string(count);
  Node *root = new Node();
  root->name = rootPath + ".jpg";
  nodes.push_back(root);
  
  rootPath = path + rootPath + ".jpg";
  image.read(rootImage);
  // root->entropy = computeEntropy(image, rootPath);
  root->entropy = 0;
  image.write(rootPath);
  
  while (count < imageCount) {
    createChild(path, count++, nodes);
  }

  return root;
}

void createChild(std::string path, int imageCount, std::vector<Node*> &nodes) {
  Image image, parent;
  int parentIndex = (imageCount > 0) ? (std::rand() % imageCount) : 0;
  Node *n = new Node();

  std::string parentPath = std::to_string(parentIndex);
  parentPath = path + parentPath + ".jpg";

  image.read(parentPath);
  parent.read(parentPath);

  applyTransform(image, n);

  std::string childPath = std::to_string(imageCount + 1);
  n->name = childPath + ".jpg";
  childPath = path + childPath + ".jpg";

  image.write(childPath);

  std::vector<std::vector<std::vector<float> > > histos;

  histos.push_back(getDctHistograms(parentPath, 16));
  histos.push_back(getDctHistograms(childPath, 16));

  auto f = [](std::vector<float> v1, std::vector<float> v2) {
    int diff = (v1.size() - v2.size()) * 0.5;

    std::vector<float> v(v1.size(), 0);

    for (int i = 0; i < v2.size(); i++) {
      v[i + diff] = v2[i];
    }
    return v;
  };

  auto h1 = histos[0][0];
  auto h2 = histos[1][0];

  if(h1.size() < h2.size()) {
    h1 = f(h2, h1);
  } else {
    h2 = f(h1, h2);
  }
  
  cv::Mat histo1(h1);
  cv::Mat histo2(h2);

  image.read(childPath);
  
  // n->entropy = compareHist(histo1, histo2, CV_COMP_BHATTACHARYYA);
  
  n->entropy = computeEntropy(image, parent);
  n->entropy2 = computeEntropy(parent, image);

  Node *parentNode = nodes[parentIndex];
  nodes.push_back(n);
  parentNode->children.push_back(n);
}

void exportPhylogeny(std::ostream& stream, Node *root) {
  std::function<void(Node *n)> f;
  f = [&](Node *n) {
    stream << "[" << std::endl;
    // stream << "\\includegraphics{" << n->name << "}" << std::endl;
    // stream << "\\href{run:" << n->name << "}{"<< n->name << " // " << n->compression << " // " << n->entropy << " // " << n->entropy2 << "}" << std::endl;

    stream << "\\href{run:" << n->name << "}{"<< n->name << "}" << std::endl;
    // stream << "\\href{run:" << n->name << "}{\\includegraphics{" << n->name << "}}" << std::endl;
    // stream << n->name << std::endl;
    for(auto i : n->children) {
      f(i);
    }
    stream << "]" << std::endl;
  };
  
  stream << "\\documentclass[tikz,border=10pt]{standalone}" << std::endl;
  stream << "\\usepackage{forest}" << std::endl;
  stream << "\\usepackage{hyperref}" << std::endl;
  stream << "\\usepackage{graphicx}" << std::endl;
  stream << "\\begin{document}" << std::endl;
  stream << "\\begin{forest}" << std::endl;

  f(root);
  
  stream << "\\end{forest}" << std::endl;
  stream << "\\end{document}" << std::endl;
}

int main(int argc, char **argv){
  std::srand(std::time(NULL));
  InitializeMagick(*argv);

  int n = 0;
  for (int i = 0; i < 8 * 2; i++)
    for (int j = (i < 8) ? 0 : i-8+1; j <= i && j < 8; j++)
      zigzag[n++] = (i&1)? j*(8-1)+i : (i-j)*8+j;

  std::vector<std::string> paths {"../data/orig_100.jpg", "../data/orig_50.jpg", "../data/orig_50_100.jpg"};
  plotDctCoefficientsMultipleImages(paths);

  // plotDctCoefficients(std::cout, argv[2]);
  Node *root = createPhylogeny(argv[1], argv[2], atoi(argv[3]));
  // std::cout << "phylogeny created" << std::endl;
  std::string treePath = std::string(argv[1]) + "tree.tex";
  std::ofstream stream(treePath);
  exportPhylogeny(stream, root);
  return 0; 
}
