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
#include "gnuplot_i.hpp"
#include "fft.h"
#include "tools.h"

using namespace Magick;
using namespace tools;

static std::vector<int> baseTable = {16,12,14,14,18,24,49,72,11,12,13,17,22,35,64,92,10,14,16,22,37,55,78,95,16,19,24,29,56,64,87,98,24,26,40,51,68,81,103,112,40,58,57,87,109,104,121,100,51,60,69,80,103,113,120,103,61,55,56,62,77,92,101,99};

static std::vector<std::vector<int> > tables;

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


struct Node {
  std::string name;
  int n;
  std::vector<Node *> children;
};

std::vector<std::vector<bool> > estimateParents(std::string directory);
int estimateQ(std::vector<int> vec);
void exportPhylogeny(std::ostream& stream, Node *root);
Node *buildTreeFromMatrix(std::vector<std::vector<bool> > matrix);
void plotHistograms(std::string outFile, std::string name1, std::string name2, std::vector<double> d1, std::vector<double> d2);
std::vector<int> makeHisto (std::vector<int> v);

void plotHistograms(std::string outFile, std::string name1, std::string name2, std::vector<double> d1, std::vector<double> d2) {
    Gnuplot g1("hello");

    g1.cmd("set xrange [0:" + std::to_string(fmax(d1.size(), d2.size())) + "]");
    g1.cmd("set style data histogram");
    g1.cmd("set style data lines");
    g1.cmd("set style histogram clustered");
    g1.cmd("set style fill solid border");
    // g1.cmd("set term png");
    g1.cmd("set terminal pngcairo size 960,720");
    g1.set_style("boxes");
      
    outFile += ".png";
    g1.cmd("set output \"/media/ramdisk/" + outFile + "\"");

    std::ofstream of("/media/ramdisk/out.txt", std::ofstream::out);
    
    for (int i = 0; i < fmax(d1.size(), d2.size()); i++) {
      of << ((i < d1.size())?d1[i]:0) << " ";
      of << ((i < d2.size())?d2[i]:0) << std::endl;
    }

    of.close();
    g1.cmd("plot '/media/ramdisk/out.txt' using 1 title \"" + name1 + "\", '/media/ramdisk/out.txt' using 2 title \"" + name2 + "\"");
}

int estimateQ(std::vector<std::vector<int> > dctCoeffs) {
  int peakFactor = 5;
  int peakCount = 3;
  int dctCount = 40; //les 4 premiers coefficients dans l'ordre zigzag sont ceux qui permettent d'avoir la plus grande précision.
  
  auto isPeak = [](std::vector<double> vec, int width, int index) {
    for (int i = index - width; i < index; i++) {
      if(vec[i] >= vec[i + 1]) return false;
    }

    for (int i = index; i < index + width; i++) {
      if(vec[i] <= vec[i + 1]) return false;
    }
    return true;
  };

  auto smooth = [](std::vector<double> &vec, int width, int index) {
      double avg = 0;
      for (int i = -width; i <= width; i++) {
	avg += vec[i + index];
      }
      avg /= (width * 2 + 1);
      return avg;
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

  auto qFromPeriods = [](std::vector<int> periods) {
    int index = 0;
    int minValue = std::numeric_limits<int>::max();
  
    for (int i = 1; i < 100; i++) {
      // std::vector<double> tabledouble(tables[i].begin(), tables[i].end());
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

  auto plotPeaks = [&](std::vector<double> plot, std::vector<int> peaks, int n) {
    double plotAvg = tools::mean(plot);
    std::vector<double> tmp(plot.size(), plotAvg * peakFactor);

    std::string s;
    for (int i = 0; i < peaks.size() && i < 10; i++) {
      s += std::to_string(peaks[i]) + " - ";
    }

    plotHistograms("autocorrelation" + std::to_string(n), "auto correlation", s, plot, tmp);
  };

  
  std::vector<int> qs;
  int n = 0;
  for (int k = 0; k < dctCount; k++) {
  // for (int k = 0; k < dctCoeffs.size(); k++) {
    auto dct = dctCoeffs[k];
    auto histo = makeHisto(dct);
    std::vector<double> vec(histo.begin(), histo.end());

    std::vector<double> plot = computeAutocorrelation(vec);
    std::vector<int> peaks = findPeaks(2, plot);

    // plotPeaks(plot, peaks, n++);

    if(peaks.size() == 0) {
      qs.push_back(1);
      // qs.push_back(-1);
      continue; 
    } else if (peaks.size() == 1) {
      qs.push_back(peaks[0]);

      continue;
    }
    
    qs.push_back(periodFromPeaks(peaks));
  }

  return qFromPeriods(qs);
}

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

void exportPhylogeny(std::ostream& stream, Node *root) {
  std::function<void(Node *n)> f;
  f = [&](Node *n) {
    stream << "[" << std::endl;
    stream << "\\href{run:" << std::to_string(n->n) << "}{" << std::to_string(n->n) << ".jpg" << "}" << std::endl;
    for(auto i : n->children) {
      if(i != n)
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

std::vector<std::vector<bool> > estimateParents(std::string directory) {
  std::vector<std::vector<bool> > matrix;
  

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

  auto toDoubleVector = [](std::vector<int> v, int n) {
    std::vector<double> out;
    for(auto i : v) {
      out.push_back((double) i / (double) n);
    }
    return out;
  };

  auto reajust = [](std::vector<int> v, int q) {
    std::vector<double> out;
    int maxElem = *std::max_element(v.begin(), v.end());
    maxElem = fmax(maxElem, abs(*std::min_element(v.begin(), v.end())));
    out.resize(maxElem+1, 0);
    
    for(auto i : v) {
      int n = (int) (i / q) * q;
      out[abs(n)]++;
    }
    return out;
  };

  auto makeArea = [](std::vector<double> v) {
    std::vector<double> out;
    double sum = 0;
    for(auto i : v) {
      sum += i;
      out.push_back(sum);
    }
    return out;
  };

  auto areaUnderTheCurve = [](std::vector<double> v) {
    double d = 0;

    for(auto i : v) {
      d += i;
    }
   
    return d;
  };
  
  auto exists = [](std::string path) {
    struct stat buffer;   
    return (stat (path.c_str(), &buffer) == 0); 
  };

  auto createPath = [](std::string directory, int imageIndex) {
    std::string s = directory + std::to_string(imageIndex) + ".jpg";
    return s;
  };

  auto getImagePathes = [&](std::string directory) {
    int imageIndex = 0;
    std::string path = createPath(directory, imageIndex);
    std::vector<std::string> pathes;

    while(exists(path)) {
      pathes.push_back(path);
      path = createPath(directory, ++imageIndex);
    }
    return pathes;
  };

  auto distancesOk = [&](std::vector<double> d1, std::vector<double> d2) {
    auto metric = Distance::kb;

    double d = Distance::compute(metric, d1, d2);
    double dd = Distance::compute(metric, d2, d1);

    //les deux distances sont égales => les images sont identiques
    return d == dd;
  };

  auto areaUnderTheCurveOk = [&](std::vector<double> d1, std::vector<double> d2) {
    return !(areaUnderTheCurve(d1) > areaUnderTheCurve(d2));
  };

  auto missingValuesOk = [&](std::vector<double> d1, std::vector<double> d2) {
    //TODO vérifier si les deux vecteurs ont la même taille
    double sum = 0;
    
    for(auto i : d1) {
      sum += i;
    }
    // sum /= d1.size();
    for(auto &i : d1) {
      i /= sum;
    }

    sum = 0;
    for(auto i : d2) {
      sum += i;
    }
    // sum /= d1.size();
    for(auto &i : d2) {
      i /= sum;
    }

    int missingValues = 0;
    int sum2 = 0;
    
    // for (int i = 0; i < fmin(d1.size(), d2.size()); i++) {
    //   sum2 += d1[i];
    //   // std::cout << sum2 << " < " << (double) sum * 0.90 << std::endl;
    //   if(d1[i] > sum * 0.1) {
    //   // if((double) sum2 < (sum * 0.99)) {
    // 	//grand d1 et petit d2, au moins 2x
    // 	if(d1[i] > 5 * d2[i]) {
    // 	  missingValues++;
    // 	}
    //   // } else {
    // 	// break;
    //   }
    // }

    // std::cout << (double) missingValues << " => ";
    // // return (((double) missingValues / (double) fmin(d1.size(), d2.size())) < 0.1);
    // return missingValues == 0;

    double d = Distance::compute(Distance::kb, d1, d2);
    std::cout << d << " => ";
    
    return d == 0;
  };

  auto imageStdDev = [&](Image &image) {
    std::cout << "here" << std::endl;

    double avg = 0;
    for (int i = 0; i < image.columns(); i++) {
      for (int j = 0; j < image.rows(); j++) {
	// std::cout << pixels->red << std::endl;
	avg += ((ColorGray) image.pixelColor(i, j)).shade();
	// std::cout << ((ColorGray) image.pixelColor(i, j)).shade() * 255 << std::endl;
      }
    }

    avg /= (image.columns() * image.rows());
    
    std::cout << "avg = " << avg << std::endl;

  };


  std::vector<std::string> pathes = getImagePathes(directory);
  
  int cpti = 0;
  int error = 0;
  for(auto image_i : pathes) {
    std::vector<bool> vec;
    Image image2, imagepng;
    image2.read(image_i);
    q_dct base = getQAndDct(image_i);
    int quality = image2.quality();
    /*
    image2.write("/media/ramdisk/i.png");
    imagepng.read("/media/ramdisk/i.png");
    imagepng.quality(100);
    imagepng.write("/media/ramdisk/i.jpg");
    //*/
    
    image2.quality(100);
    image2.write("/media/ramdisk/i.jpg");
    //*/
    
    // q_dct i = getQAndDct(image_i);
    q_dct i = getQAndDct("/media/ramdisk/i.jpg");
    auto di = makeDistrib(makeHisto(i.dct[0]));
    // auto di = toDoubleVector(i.dctDiffZero, i.numberOfBlocks);

    // if(cpti == 0) {
    //   // fft(toDoubleVector(makeHisto(i.dct[0]), 1));
    //   plotHistograms("fft", "fft", "histo", fft(toDoubleVector(makeHisto(i.dct[0]), 1)), toDoubleVector(makeHisto(i.dct[0]), 1));
    // }

    // if(cpti == 19) {
    std::cout << cpti  << " : " << std::endl;
    int estimated = estimateQ(i.dct);
    std::cout << "estimated = " << estimated << std::endl;
    std::cout << "real      = " << quality << std::endl;
    error += (abs(estimated - quality));
    std::cout << std::endl;
    // }
      // int n = 0;
      // for(auto dctVec : i.dct) {
	// plotHistograms(std::to_string(n), std::to_string(cpti), std::to_string(n++), toDoubleVector(makeHisto(dctVec), 1), toDoubleVector(makeHisto(dctVec), 1));
      // // plotHistograms(std::to_string(n), std::to_string(cpti), std::to_string(n++), fft(dctVec), fft(dctVec)); // 

      // 	plotHistograms(std::to_string(n), std::to_string(cpti), std::to_string(n++), toDoubleVector(makeHisto(dctVec), 1), toDoubleVector(makeHisto(dctVec), 1));
      // }
      
    // }
      
    // if(cpti == 0) {
      // int n = 0;
      // for(auto dctVec : i.dct) {
      // plotHistonnnngrams(std::to_string(n), std::to_string(cpti), std::to_string(n++), toDoubleVector(makeHisto(dctVec), 1), toDoubleVector(makeHisto(dctVec), 1));
    // 	// plotHistograms(std::to_string(n), std::to_string(cpti), std::to_string(n++), toDoubleVector(makeHisto(dctVec), 1), toDoubleVector(makeHisto(dctVec), 1));
      // }
    // }
      // std::cout << base.q << " // " << estimateQ(i.dct) << std::endl;
    // }

    //********************************//
    /*
    int cptj = 0;
    for(auto image_j : pathes) {
      Image image;
      image.read(image_j);

      if(cpti == cptj) {
    	vec.push_back(false);
      // } else if(i.q >= getQAndDct(image_j).q) {
      // } else if(image2.quality() < image.quality()) {
      } else if (1) {
    	image.quality(quality);
    	image.write("/media/ramdisk/out.jpg");
    	Image image3;
    	image3.read("/media/ramdisk/out.jpg");
    	image3.quality(100);
    	image3.write("/media/ramdisk/j.jpg");
    	// q_dct j = getQAndDct("out.jpg");
    	q_dct j = getQAndDct("/media/ramdisk/j.jpg");
    	q_dct jUncompressed = getQAndDct(image_j);

    	// q_dct j = getQAndDct(image_j);
	
    	auto dj = makeDistrib(makeHisto(j.dct[0]));
    	// auto dj = toDoubleVector(j.dctDiffZero, j.numberOfBlocks);

    	// if(cpti == 8) {
    	// plotHistograms(std::to_string(cpti) + " - " + std::to_string(cptj), std::to_string(cpti), std::to_string(cptj), fft(makeHisto(i.dct[2])), fft(makeHisto(j.dct[2])));
    	  // plotHistograms(std::to_string(cptj), std::to_string(cpti), std::to_string(cptj), toDoubleVector(makeHisto(base.dct[0]), 1), toDoubleVector(makeHisto(i.dct[0]), 1));
    	  // plotHistograms(std::to_string(cptj), std::to_string(cpti), std::to_string(cptj), reajust(i.dct[0], base.q), reajust(j.dct[0], base.q));
    	  // plotHistograms(std::to_string(cptj), std::to_string(cpti), std::to_string(cptj), toDoubleVector(base.dctDiffZero, i.numberOfBlocks), toDoubleVector(jUncompressed.dctDiffZero, j.numberOfBlocks));
    	  // plotHistograms(std::to_string(cptj), toDoubleVector(i.dctDiffZero, i.numberOfBlocks), toDoubleVector(jUncompressed.dctDiffZero, jUncompressed.numberOfBlocks));
    	// }
	
    	std::cout << cpti << " and " << cptj << " = ";
    	bool b1 = distancesOk(di, dj);
    	bool b2 = areaUnderTheCurveOk(toDoubleVector(base.dctDiffZero, i.numberOfBlocks), toDoubleVector(jUncompressed.dctDiffZero, j.numberOfBlocks));
    	bool b3 = missingValuesOk(reajust(i.dct[0], base.q), reajust(j.dct[0], base.q));
    	bool b4 = distancesOk(fft(makeHisto(i.dct[0])), fft(makeHisto(j.dct[0])));
    	bool b5 = missingValuesOk(fft(makeHisto(i.dct[1])), fft(makeHisto(j.dct[1])));

    	vec.push_back(b1 || (b2 && b3 && b4 && b5));
    	// vec.push_back(b1);
    	std::cout << b2 << " " << b3 << " " << b4 << " " << b5 << std::endl;
    	// std::cout << cpti << " and " << cptj << " = ";
    	// canBeChild(toDoubleVector(i.dctDiffZero, i.numberOfBlocks), toDoubleVector(j.dctDiffZero, j.numberOfBlocks), i.q);
    	// vec.push_back(!canBeChild(toDoubleVector(i.dctDiffZero, i.numberOfBlocks), toDoubleVector(j.dctDiffZero, j.numberOfBlocks), i.q));
      } else {
    	vec.push_back(false);
      }
      cptj++;
    }
    //*/
    matrix.push_back(vec);
    cpti++;
  }

  for(auto i : matrix) {
    for(auto j : i) {
      std::cout << j << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "mean error = " << (double)error / (double)cpti << std::endl;

  return matrix;
}

int main(int argc, char **argv) {
  InitializeMagick(*argv);

  for (int i = 1; i < 101; i++) {
    tables.push_back(table(i));
  }


  auto matrix = estimateParents(argv[1]);
  Node *root = buildTreeFromMatrix(matrix);

  // std::vector<int> v = {600, 470, 170, 300, 430};
  
  std::string treePath = "tree.tex";
  std::ofstream stream(treePath);
  std::cout << "exporting" << std::endl;
  exportPhylogeny(stream, root);
  
  return 0;
}
