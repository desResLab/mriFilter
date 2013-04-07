#ifndef MRITHRESHOLDCRITERIA_H
#define MRITHRESHOLDCRITERIA_H

class MRIThresholdCriteria
{
public:
  // Data Members
  int thresholdType;
  int thresholdQty;
  double thresholdValue;
  // Constructor and Destructor
  MRIThresholdCriteria(int type, int qty, double value);
  ~MRIThresholdCriteria();
  // Data Members
  bool MeetsCriteria(double currentValue);
};

#endif // MRITHRESHOLDCRITERIA_H
