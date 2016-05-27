#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <random>
#include <sstream>
#include <dirent.h>
#include <stdio.h>

std::vector<double> parseFile(std::string path) {
  std::ifstream file(path);
  std::string line;
  std::vector<double> v;
  
  while (std::getline(file, line)) {
    // std::cout << line << std::endl;
    std::istringstream iss(line);
      
    std::string s;
    double d;
      
    if (!(iss >> s >> s >> d)) { break; } // error

    v.push_back(d);
  }
  return v;
}

int main(int argc, char **argv){
  std::cout << argv[1] << std::endl;

  /**
   * 0 = mean error
   * 1 = overestimate
   * 2 = underestimate
   * 3 = roots
   * 4 = edges
   * 5 = leaves
   * 6 = ancestry
   */
  std::vector<double> results(7);

  int n = 0;

  const char* PATH = argv[1];

  DIR *dir = opendir(PATH);

  struct dirent *entry = readdir(dir);

  while (entry != NULL) {
    if (entry->d_type == DT_DIR) {
      std::string path(entry->d_name);
	
      if(path[0] != '.') {
	auto v = parseFile(std::string(argv[1]) + "/" + path + "/results.txt");
	n += (v.size() > 0);
	
	for (int i = 0; i < v.size(); i++) {
	  results[i] += v[i];
	}

      }
    }
      
    entry = readdir(dir);
  }

  closedir(dir);


  for(auto &i : results) {
    i /= n;
  }
  
  std::cout << "mean error = " << results[0] << std::endl;
  std::cout << "overestimate = " << results[1] << std::endl;
  std::cout << "underestimate = " << results[2] << std::endl;
  std::cout << "roots = " << results[3] << std::endl;
  std::cout << "edges = " << results[4] << std::endl;
  std::cout << "leaves = " << results[5] << std::endl;
  std::cout << "ancestry = "     << results[6] << std::endl;
    
  return 0; 
}
