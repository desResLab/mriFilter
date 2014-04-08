#ifndef MRIOPTIONS_H
#define MRIOPTIONS_H

#include "mriThresholdCriteria.h"

class MRIOptions{
  public:
    // DATA MEMBERS
  	double tolerance;
    int maxIterations;		
    bool useBCFilter;
    bool useConstantPatterns;
    MRIThresholdCriteria* thresholdCriteria;
    // CONSTRUCTOR AND DESTRUCTOR
    MRIOptions(double tol, int maxIt);
    ~MRIOptions(){}
};
	
#endif // MRIOPTIONS_H
