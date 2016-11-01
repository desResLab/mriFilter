#ifndef MRITHRESHOLDCRITERIA_H
#define MRITHRESHOLDCRITERIA_H

# include <math.h>

# include "mriConstants.h"

class MRIThresholdCriteria{
public:
  // Data Members
  int thresholdType;
  int thresholdQty;
  double thresholdValue;
  // Constructor and Destructor
  MRIThresholdCriteria(int type, int qty, double value);
  ~MRIThresholdCriteria();
  // Data Members
  bool meetsCriteria(double currentValue);
};

#endif // MRITHRESHOLDCRITERIA_H
