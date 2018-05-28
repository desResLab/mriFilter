#ifndef MRISAMPLINGOPTIONS_H
#define MRISAMPLINGOPTIONS_H

# include <string>

# include <mriTypes.h>

class MRISamplingOptions{
  public:
    int totalSamples;
    std::string outputFile;
    int numberOfBins;
    bool useBox;
    MRIDoubleVec limitBox;
    // Constructor and Destructor
    MRISamplingOptions();
    ~MRISamplingOptions();
};

#endif // MRISAMPLINGOPTIONS_H
