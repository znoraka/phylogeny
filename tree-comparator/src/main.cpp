#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <functional>
#include <fstream>
#include <ctime>
#include <Magick++.h>
#include <Magick++/STL.h>
#include <magick/MagickCore.h>
#include <map>
#include <algorithm>
#include <random>

using namespace Magick;

struct Node {
  std::string name;
  int n;
  std::vector<Node *> children;
  Node *parent = nullptr;
};

std::vector<Node *> buildTreeFromMatrix(std::vector<std::vector<bool> > matrix);
std::vector<std::vector<bool> > parseMatrixFromFile(std::string path);
bool metric_root(std::vector<Node *> tree1, std::vector<Node *> tree2);
double metric_edges(std::vector<Node *> tree1, std::vector<Node *> tree2);
double metric_leaves(std::vector<Node *> tree1, std::vector<Node *> tree2);
double metric_ancestry(std::vector<Node *> tree1, std::vector<Node *> tree2);

std::vector<Node *> buildTreeFromMatrix(std::vector<std::vector<bool> > matrix) {
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
  std::vector<Node *> nodes;

  Node *root = nullptr;
  ref *nextRoot = nullptr;

  int cpt = 0;
  for(auto i : matrix) {
    Node *n = new Node();
    nodes.push_back(n);
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
	e->node->parent = r->node;
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

  for(auto i : nodes) {
    std::sort(i->children.begin(), i->children.end(), [](Node *n1, Node *n2) {
	return n1->n < n2->n;
      });
  }

  return nodes;
}

bool metric_root(std::vector<Node *> tree1, std::vector<Node *> tree2) {
  // int root1 = -1;
  // int root2 = -2;

  // int n = 0;
  // for(auto i : tree1) {
  //   if(i->parent == nullptr) {
  //     if(root1 != -1) {
  // 	return false;	
  //     }
  //     root1 = i->n;
  //   }

  //   n++;
  // }

  // for(auto i : tree2) {
  //   if(i->parent == nullptr) {
  //     if(root2 != -2) return false;
  //     root2 = i->n;
  //   }
  // }

  // return root1 == root2;

  for(auto i : tree1) {
    for(auto j : tree2) {
      if(i->parent == nullptr &&
	 j->parent == nullptr) {
	if(i->n == j->n) return true;
  	break;
      }
    }
  }
  
  return false;
}
		    
double metric_edges(std::vector<Node *> tree1, std::vector<Node *> tree2) {
  int n = 0;
  
  for(auto i : tree1) {
    for(auto j : tree2) {
      if(i->parent != nullptr &&
	 j->parent != nullptr &&
	 i->n == j->n) {
	n += (i->parent->n == j->parent->n);
	break;
      }
    }
  }

  return (double)n / (double)(tree1.size() - 1);
}

double metric_leaves(std::vector<Node *> tree1, std::vector<Node *> tree2) {
  int l_union = 0;
  int l_inter = 0;

  for(auto i : tree1) {
    for(auto j : tree2) {
      if(i->n == j->n) {
  	l_inter += (i->children.size() == 0 && j->children.size() == 0);
  	l_union += (i->children.size() == 0 || j->children.size() == 0);
  	break;
      }
    }
  }
  
  return (double) l_inter / (double) l_union;
}

double metric_ancestry(std::vector<Node *> tree1, std::vector<Node *> tree2) {
  std::vector<std::string> set_tree1;
  std::vector<std::string> set_tree2;

  for(auto i : tree1) {
    if(i->parent != nullptr) {
      set_tree1.push_back(std::to_string(i->n) + "/" + std::to_string(i->parent->n));
    }
  }

  for(auto i : tree2) {
    if(i->parent != nullptr) {
      set_tree2.push_back(std::to_string(i->n) + "/" + std::to_string(i->parent->n));
    }
  }

  std::vector<std::string> v_union(set_tree1.size() + set_tree2.size());
  std::vector<std::string> v_inter(set_tree1.size() + set_tree2.size());

  std::vector<std::string>::iterator it;

  it = std::set_intersection (set_tree1.begin(), set_tree1.end(), set_tree2.begin(), set_tree2.end(), v_inter.begin());
  v_inter.resize(it-v_inter.begin());

  it = std::set_union (set_tree1.begin(), set_tree1.end(), set_tree2.begin(), set_tree2.end(), v_union.begin());
  v_union.resize(it-v_union.begin());

  return (double)v_inter.size() / (double)v_union.size();
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


int main(int argc, char **argv){
  std::srand(std::time(NULL));

  auto t1 = buildTreeFromMatrix(parseMatrixFromFile(argv[1]));
  auto t2 = buildTreeFromMatrix(parseMatrixFromFile(argv[2]));

  std::cout << "roots = " << metric_root(t1, t2) << std::endl;
  std::cout << "edges = " << metric_edges(t1, t2) << std::endl;
  std::cout << "leaves = " << metric_leaves(t1, t2) << std::endl;
  std::cout << "ancestry = " << metric_ancestry(t1, t2) << std::endl;
  
  return 0; 
}
