#ifndef MRISAMPLINGOPTIONS_H
#define MRISAMPLINGOPTIONS_H

# include <string>
# include <mriTypes.h>

class mriSamplingOptions{
  public:
    int totalSamples;
    std::string outputFile;
    int numberOfBins;
    bool useBox;
    mriDoubleVec limitBox;
    // Constructor and Destructor
    mriSamplingOptions();
    ~mriSamplingOptions();
};

#endif // MRISAMPLINGOPTIONS_H
