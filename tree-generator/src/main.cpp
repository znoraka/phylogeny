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

using namespace Magick;
using namespace cv;

struct Node {
  std::string name;
  int compression;
  std::vector<Node*> children;
};

void createChild(std::string path, int imageCount, std::vector<Node*> &nodes);
Node *createPhylogeny(std::string path, std::string rootImage, int imageCount);
void applyTransform(Image &image, Node *node);
void exportPhylogeny(std::ostream& stream, Node *root);
void recompress(Image &parent, Node *node, int parentQ);

void recompress(Image &image, Node *node, int parentQ) {
  int q = fmax(30, parentQ - (1 + std::rand() % 4) * 5);
  image.quality(q);
  node->compression = q;
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
  recompress(image, root, 100);
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
  Node *parentNode = nodes[parentIndex];


  std::string parentPath = std::to_string(parentIndex);
  parentPath = path + parentPath + ".jpg";

  image.read(parentPath);

  recompress(image, n, parentNode->compression);

  std::string childPath = std::to_string(imageCount + 1);
  n->name = childPath + ".jpg";
  childPath = path + childPath + ".jpg";

  image.write(childPath);
  
  nodes.push_back(n);
  parentNode->children.push_back(n);
}

void exportPhylogeny(std::ostream& stream, Node *root) {
  std::function<void(Node *n)> f;
  f = [&](Node *n) {
    stream << "[" << std::endl;
    stream << "\\href{run:" << n->name << "}{"<< n->name << " // " << n->compression << "}" << std::endl;
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

  Node *root = createPhylogeny(std::string(argv[1]) + "/", argv[2], atoi(argv[3]));
  // std::cout << "phylogeny created" << std::endl;
  std::string treePath = std::string(argv[1]) + "/tree.tex";
  std::ofstream stream(treePath);
  exportPhylogeny(stream, root);
  return 0; 
}
