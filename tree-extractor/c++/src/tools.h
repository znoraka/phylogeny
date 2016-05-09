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
}
