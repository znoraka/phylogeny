namespace tools {

  template<typename T>
  double mean(T container) {
    double d = 0;
    int n = 0;

    for(auto i : container) {
      d += i;
      n++;
    }

    return d / n;
  }

  template<typename T>
  double standardDeviation(T container, double meanValue = std::numeric_limits<double>::max()) {
    double d = 0;
    int n = 0;

    if(meanValue == std::numeric_limits<double>::max()) {
      meanValue = mean(container);
    }

    for(auto i : container) {
      d += pow(i - meanValue, 2);
      n++;
    }

    return sqrt(d / n);
  }

  #define D Distance::distanceFormula Distance

  class Distance {
  private:
    struct distanceFormula {
      std::function<double(double)> global;
      std::function<double(double, double)> perItem;
    
      distanceFormula(std::function<double(double)> global,
		      std::function<double(double, double)> perItem) {
	this->global = global;
	this->perItem = perItem;
      }
    };
  
  public:
    static double compute(distanceFormula d,
			  std::vector<double>::iterator begin1,
			  std::vector<double>::iterator end1,
			  std::vector<double>::iterator begin2,
			  std::vector<double>::iterator end2) {
      double sum = 0;
      for (int i = 0; begin1 + i != end1 && begin2 + i != end2; i++) {
	double tmp = d.perItem(*(begin1 + i), *(begin2 + i));
	if(tmp >= 0)
	  sum += tmp;
      }
      return d.global(sum);
    }

    static double compute(distanceFormula d, std::vector<double> d1, std::vector<double> d2) {
      return compute(d, d1.begin(), d1.end(), d2.begin(), d2.end());
    }

    static distanceFormula euclid, autocorrelation, bhattacharyya, kb, hellinger, hellinger_b, s_hellinger, jeffreys, kdiv;

  };

  D::euclid = distanceFormula([](double d){return sqrt(d);},
			      [](double i, double j) {
				return pow(i - j, 2);}
			      );

  D::autocorrelation = distanceFormula([](double d){return d;},
				       [](double i, double j) {
					 return i * j;}
				       );

  D::bhattacharyya = distanceFormula([](double d){return log(d);},
				     [](double i, double j) {
				       return sqrt(i * j);}
				     );

  D::kb = distanceFormula([](double d){return d;},
			  [](double i, double j) {
			    if(i == 0 || j == 0) return 0.;
			    return i * log(i / j);}
			  );

  D::hellinger = distanceFormula([](double d){return sqrt(2 * d);},
				 [](double i, double j) {
				   return pow(sqrt(i) - sqrt(j) , 2);}
				 );

  D::hellinger_b = distanceFormula([](double d){return 2 * sqrt(1 - d);},
				   [](double i, double j) {
				     return sqrt(i * j);}
				   );

  D::s_hellinger = distanceFormula([](double d){return 2 * d;},
				   [](double i, double j) {
				     return pow(sqrt(i) - sqrt(j) , 2);}
				   );

  D::jeffreys = distanceFormula([](double d){return d;},
				[](double i, double j) {
				  return (i - j) * log(i / j);}
				);

  D::kdiv = distanceFormula([](double d){return d;},
			    [](double i, double j) {
			      return i * log((2 * i) / (i + j));}
			    );

}
