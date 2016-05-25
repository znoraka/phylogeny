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
#include <algorithm>
#include <random>

#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>

using namespace Magick;
using namespace cv;

struct Node {
  std::string name;
  int compression;
  std::vector<Node*> children;
};

void createChild(std::string path, int imageCount, std::vector<Node*> &nodes, bool randomQ, std::vector<int> randomNames);
Node *createPhylogeny(std::string path, std::string rootImage, int imageCount, bool randomQ);
void applyTransform(Image &image, Node *node);
void exportPhylogeny(std::ostream& stream, Node *root);
void recompress(Image &parent, Node *node, int parentQ);

void recompress(Image &image, Node *node, int parentQ, bool randomQ = false) {
  int q;
  if(randomQ) {
    do {
      q = 40 + std::rand() % 50;
    } while(q == parentQ);
  } else {
    q = fmax(30, parentQ - (1 + std::rand() % 15));
  }
  // int q = fmax(30, parentQ - (1 + std::rand() % 3) * 5);
  // int q;
  
  image.quality(q);
  node->compression = q;
}


Node *createPhylogeny(std::string path, std::string rootImage, int imageCount, bool randomQ) {
  std::vector<int> randomNames;
  for (int i = 0; i < imageCount; i++) {
    randomNames.push_back(i);
  }
  auto engine = std::default_random_engine(std::time(NULL));
  std::shuffle(randomNames.begin(), randomNames.end(), engine);

  Image image;
  int count = 0;
  std::vector<Node*> nodes;
  nodes.resize(imageCount);

  std::string rootPath = std::to_string(randomNames[count]);
  Node *root = new Node();
  root->name = rootPath + ".jpg";
  nodes[randomNames[count]] = root;
  
  rootPath = path + rootPath + ".jpg";
  image.read(rootImage);
  recompress(image, root, 100, randomQ);
  image.write(rootPath);
  
  while (count < imageCount - 1) {
    createChild(path, count++, nodes, randomQ, randomNames);
  }

  return root;
}

void createChild(std::string path, int imageCount, std::vector<Node*> &nodes, bool randomQ, std::vector<int> randomNames) {
  Image image, parent;
  int parentIndex = randomNames[(imageCount > 0) ? (std::rand() % imageCount) : 0];
  Node *n = new Node();
  Node *parentNode = nodes[parentIndex];


  std::string parentPath = std::to_string(parentIndex);
  parentPath = path + parentPath + ".jpg";

  image.read(parentPath);

  recompress(image, n, parentNode->compression, randomQ);

  std::string childPath = std::to_string(randomNames[imageCount + 1]);
  n->name = childPath + ".jpg";
  childPath = path + childPath + ".jpg";

  image.write(childPath);
  
  nodes[randomNames[imageCount + 1]] = n;
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

  system("rm -rf /media/ramdisk/data && mkdir /media/ramdisk/data");
  Node *root = createPhylogeny(std::string(argv[1]) + "/", argv[2], atoi(argv[3]), atoi(argv[4]));
  std::string treePath = std::string(argv[1]) + "/tree.tex";
  std::ofstream stream(treePath);
  exportPhylogeny(stream, root);
  return 0; 
}
