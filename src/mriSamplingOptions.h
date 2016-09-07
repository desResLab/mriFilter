#ifndef MRISAMPLINGOPTIONS_H
#define MRISAMPLINGOPTIONS_H

#include <string>

class MRISamplingOptions
{
  public:
    int totalSamples;
    std::string outputFile;
    int numberOfBins;
    bool useBox;
    double limitBox[6];
    // Constructor and Destructor
    MRISamplingOptions();
    ~MRISamplingOptions();
};

#endif // MRISAMPLINGOPTIONS_H
