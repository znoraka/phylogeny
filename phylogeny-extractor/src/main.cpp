#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <ctime>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <sstream>

struct Node {
  std::string name;
  int n;
  std::vector<Node *> children;
};

Node *buildTreeFromMatrix(std::vector<std::vector<bool> > matrix);
std::vector<std::vector<bool> > matrixFromTree(Node *root);
void exportPhylogeny(std::ostream& stream, Node *root);
bool compareTrees(Node *r1, Node *r2);
std::vector<std::vector<bool> > parseMatrixFromFile(std::string path);

Node *buildTreeFromMatrix(std::vector<std::vector<bool> > matrix) {
  struct ref {
    int index;
    int sum;
    Node *node;
    bool done;
    ref(int index, int sum, Node *node) {
      this->index = index;
      this->sum = sum;
      this->node = node;
      this->done = false;
    }
  };

  std::vector<ref*> refs;

  Node *root = nullptr;
  ref *nextRoot = nullptr;

  int cpt = 0;
  for(auto i : matrix) {
    Node *n = new Node();
    n->n = cpt;
    ref *r = new ref(cpt++, std::accumulate(i.begin(), i.end(), 0), n);
    refs.push_back(r);

    if(nextRoot == nullptr || nextRoot->sum > r->sum) {
      nextRoot = r;
    }
  }

  for(auto kkk : matrix) {
    int nextSum = matrix.size() + 1;
    ref *r = nextRoot;
    r->done = true;
   
    if(root == nullptr) {
      root = r->node;
      r->sum -= 1;
    }
    
    for (int i = 0; i < matrix.size(); i++) {
      ref *e = refs[i];

      e->sum -= matrix[i][r->index];
      matrix[i][r->index] = 0;

      if(e->sum == 0) {
	e->sum -= 1;
	r->node->children.push_back(e->node);
      }

      if(!e->done && e->sum < nextSum) {
	nextRoot = e;
	nextSum = e->sum;
      }
    }
  }

  for(auto i : refs) {
    delete i;
  }
  
  return root;
}

std::vector<std::vector<bool> > matrixFromTree(Node *root) {
  int n = 0;
  std::function<void(Node *root)> f;
  f = [&](Node *root) {
    n++;
    for(auto i : root->children) {
      f(i);
    }
  };
  f(root);

  std::vector<std::vector<bool> > matrix;

  for (int i = 0; i < n; i++) {
    std::vector<bool> b(n, false);
    matrix.push_back(b);
  }

  std::function<void(Node *root, std::vector<int> parents)> f1;
  f1 = [&](Node *root, std::vector<int> parents) {
    for(auto i : parents) {
      matrix[root->n][i] = true;

      //random pour tester la robustesse de l'algo
      // if(std::rand() % 100 < 20) {
      // 	matrix[root->n][i] = !matrix[root->n][i];
      // }
    }
    
    parents.push_back(root->n);

    // if(parents.size() > 3) {
    //   // parents.pop_back();
    //   parents = std::vector<int>(parents.begin() + 1, parents.end());
    // }

    for(auto i : root->children) {
      f1(i, parents);
    }
  };

  std::vector<int> v;
  f1(root, v);

  return matrix;
}

void exportPhylogeny(std::ostream& stream, Node *root, Node *root2) {
  std::function<void(Node *n)> f;
  f = [&](Node *n) {
    stream << "[" << std::endl;
    stream << "\\href{run:" << std::to_string(n->n) << "}{" << std::to_string(n->n) << ".jpg" << "}" << std::endl;
    // stream << n->n << std::endl;
    for(auto i : n->children) {
      // stream << n->n << "--" << i->n << ";";
      f(i);
    }
    stream << "]" << std::endl;
  };

  // stream << "graph model {";
  // f(root);
  // stream << "}" << std::endl;

  // stream << "graph computed {";
  // f(root2);
  // stream << "}";

  
  stream << "\\documentclass[tikz,border=10pt]{standalone}" << std::endl;
  stream << "\\usepackage{forest}" << std::endl;
  stream << "\\usepackage{hyperref}" << std::endl;
  stream << "\\begin{document}" << std::endl;

  stream << "\\begin{forest}" << std::endl;
  f(root);
  stream << "\\end{forest}" << std::endl;
  
  // stream << "\\begin{forest}" << std::endl;
  // f(root2);
  // stream << "\\end{forest}" << std::endl;

  stream << "\\end{document}" << std::endl;
}

bool compareTrees(Node *r1, Node *r2) {
  if(r1->n != r2->n) {
    std::cout << "error on nodes " << r1->n << " and " << r2->n << std::endl;
    return false;
  }
  if(r1->children.size() != r2->children.size()) return false;
  
  for (int i = 0; i < r1->children.size(); i++) {
    return compareTrees(r1->children[i], r2->children[i]);
  }
  return true;
}

std::vector<std::vector<bool> > parseMatrixFromFile(std::string path) {
  std::vector<std::vector<bool> > matrix;
  std::ifstream infile(path);

  std::string str;

  do {
    std::getline(infile, str);
    str.erase(std::remove(str.begin(), str.end(), '\''), str.end());
    str.erase(std::remove(str.begin(), str.end(), '#'), str.end());
    str.erase(std::remove(str.begin(), str.end(), '('), str.end());
    str.erase(std::remove(str.begin(), str.end(), ')'), str.end());

    std::string s;
    std::stringstream ss;
    ss << str;
    std::vector<bool> v;
    while(std::getline(ss, s, ' ')) {
      if(s.compare("1") == 0) {
	v.push_back(true);
      } else if(s.compare("0") == 0) {
	v.push_back(false);
      }
    }
    if(v.size() > 0)
      matrix.push_back(v);
  } while (str.size() > 1);

  return matrix;
}

int main(int argc, char **argv) {
  std::srand(std::time(NULL));
  std::vector<Node*> nodes;

  // parseMatrixFromFile(argv[1]);

  // auto f = [&](Node *root) {
  //   Node *n = new Node();
  //   n->n = nodes.size();
  //   nodes.push_back(n);
  //   root->children.push_back(n);
  // };

  // Node *root = new Node();
  // root->n = 0;
  // nodes.push_back(root);

  // for (int i = 0; i < atoi(argv[1]); i++) {
  //   f(nodes[std::rand() % (nodes.size())]);
  // }
  
  // Node *out = buildTreeFromMatrix(matrixFromTree(root));
  Node *out = buildTreeFromMatrix(parseMatrixFromFile(argv[1]));

  // // if(compareTrees(out, root)) {
  // //   std::cout << "trees are the same!" << std::endl;
  // // } else {
  // //   std::cout << "trees are not the same" << std::endl;
  // // }
  
  std::string treePath = "tree.tex";
  std::ofstream stream(treePath);
  std::cout << "exporting" << std::endl;
  exportPhylogeny(stream, out, out);
  // exportPhylogeny(std::cout, out, root);
  
  // std::cout << "done" << std::endl;

  return 0;
}
