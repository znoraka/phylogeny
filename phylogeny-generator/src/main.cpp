#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <ctime>
#include <Magick++.h>
#include <map>
extern "C" {
#include <jpeglib.h>
}

#include "gnuplot_i.hpp"

using namespace Magick;

struct Node {
  std::string name;
  std::vector<Node*> children;
};

void createChild(std::string path, int imageCount, std::vector<Node*> &nodes);
Node *createPhylogeny(std::string path, std::string rootImage, int imageCount);
void applyTransform(Image& image);
void exportPhylogeny(std::ostream& stream, Node *root);
void plotDctCoefficients(std::ostream &stream, std::string path);

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


  //parcours zigzag
  int zigzag[64];
  int n = 0;
  for (int i = 0; i < 8 * 2; i++)
    for (int j = (i < 8) ? 0 : i-8+1; j <= i && j < 8; j++)
      zigzag[n++] = (i&1)? j*(8-1)+i : (i-j)*8+j;

  for (int channel = 0; channel < 1; channel++) {
    // for (int dctIndex = 0; dctIndex < 1; dctIndex++) {
    for (int dctIndex = 0; dctIndex < 16; dctIndex++) {
      std::map<int, int> histo;
      for (int i = 0; i < info.comp_info->height_in_blocks; i++) {
	for (int j = 0; j < info.comp_info->width_in_blocks; j++) {

	  
	  int n = info.mem->access_virt_barray((j_common_ptr)&info, coeffs[0], i, 1, false)[0][j][zigzag[dctIndex]];

	  auto s = histo.find(n);
	  
	  if(s != histo.end()) {
	    s->second += 1;
	  } else {
	    histo.insert(std::pair<int, int>(n, 1));
	  }
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

void applyTransform(Image& image) {
  std::vector<std::function<void(Image &image)> > functions;

  //resize
  functions.push_back([&](Image &image) {
      int cols = image.columns() * (0.9 + (std::rand() % 100) * 0.002);
      int rows = image.rows() * (0.9 + (std::rand() % 100) * 0.002);
      std::string s = std::to_string(cols) + "x" + std::to_string(rows) + "!";      
      image.resize(s);
    });

  //recompress
  functions.push_back([&](Image &image) {
      image.quality(25 + std::rand() % 75);
    });

  //gamma adjustment
  functions.push_back([&](Image &image) {
      image.gamma(0.3 + (std::rand() % 200) * 0.01);
    });

  //brightness adjustment
  functions.push_back([&](Image &image) {
      image.modulate((std::rand() % 100) * 0.4 + 80, 100, 100);
    });

  //rotation
  functions.push_back([&](Image &image) {
      image.rotate((std::rand() % 10) - 5);
    });


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
  image.write(rootPath);
  
  while (count < imageCount) {
    createChild(path, count++, nodes);
  }

  return root;
}

void createChild(std::string path, int imageCount, std::vector<Node*> &nodes) {
  Image image;
  int parentIndex = (imageCount > 0) ? (std::rand() % imageCount) : 0;
  Node *n = new Node();

  std::string parentPath = std::to_string(parentIndex);
  parentPath = path + parentPath + ".jpg";

  image.read(parentPath);

  applyTransform(image);

  std::string childPath = std::to_string(imageCount + 1);
  n->name = childPath + ".jpg";
  childPath = path + childPath + ".jpg";
  image.write(childPath);

  Node *parent = nodes[parentIndex];
  nodes.push_back(n);
  parent->children.push_back(n);
}

void exportPhylogeny(std::ostream& stream, Node *root) {
  std::function<void(Node *n)> f;
  f = [&](Node *n) {
    stream << "[" << std::endl;
    stream << "\\href{run:" << n->name << "}{" << n->name << "}" << std::endl;
    // stream << n->name << std::endl;
    for(auto i : n->children) {
      f(i);
    }
    stream << "]" << std::endl;
  };
  
  stream << "\\documentclass[tikz,border=10pt]{standalone}" << std::endl;
  stream << "\\usepackage{forest}" << std::endl;
  stream << "\\usepackage{hyperref}" << std::endl;
  stream << "\\begin{document}" << std::endl;
  stream << "\\begin{forest}" << std::endl;

  f(root);
  
  stream << "\\end{forest}" << std::endl;
  stream << "\\end{document}" << std::endl;
}

int main(int argc, char **argv){
  std::srand(std::time(NULL));
  InitializeMagick(*argv);

  plotDctCoefficients(std::cout, argv[2]);


  Node *root = createPhylogeny(argv[1], argv[2], atoi(argv[3]));
  std::string treePath = std::string(argv[1]) + "tree.tex";
  std::ofstream stream(treePath);
  exportPhylogeny(stream, root);
  return 0; 
}
