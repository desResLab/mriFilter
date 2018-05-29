#ifndef MRITHRESHOLDCRITERIA_H
#define MRITHRESHOLDCRITERIA_H

# include <math.h>

# include "mriConstants.h"

class mriThresholdCriteria{
public:
  // Data Members
  int thresholdType;
  int thresholdQty;
  double thresholdValue;
  // Constructor and Destructor
  mriThresholdCriteria(int type, int qty, double value);
  ~mriThresholdCriteria();
  // Data Members
  bool meetsCriteria(double currentValue);
};

#endif // MRITHRESHOLDCRITERIA_H
