#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <ctime>
#include <Magick++.h>

using namespace Magick;

struct Node {
  std::string name;
  std::vector<Node*> children;
};

void createChild(std::string path, int imageCount, std::vector<Node*> &nodes);
Node *createPhylogeny(std::string path, std::string rootImage, int imageCount);
void applyTransform(Image& image);
void exportPhylogeny(std::ostream& stream, Node *root);

void applyTransform(Image& image) {
  // image.addNoise(GaussianNoise); 
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
  std::srand(std::time(0));
  InitializeMagick(*argv);

  Node *root = createPhylogeny(argv[1], argv[2], atoi(argv[3]));
  std::string treePath = std::string(argv[1]) + "tree.tex";
  std::ofstream stream(treePath);
  exportPhylogeny(stream, root);
  return 0; 
}
