#include <math.h>
#include "mriThresholdCriteria.h"
#include "mriConstants.h"

MRIThresholdCriteria::MRIThresholdCriteria(int qty, int type, double value){
  thresholdQty = qty;
  thresholdType = type;
  thresholdValue = value;
}

MRIThresholdCriteria::~MRIThresholdCriteria()
{
}

// Check If a Value Meets the Predefined Criteria
bool MRIThresholdCriteria::MeetsCriteria(double currentValue){
  switch(thresholdType){
    case kCriterionLessThen:
      return (currentValue<thresholdValue);
      break;
    case kCriterionGreaterThen:
      return (currentValue>thresholdValue);
      break;
    case kCriterionABSLessThen:
      return (fabs(currentValue)<thresholdValue);
      break;
    case kCriterionABSGreaterThen:
      return (fabs(currentValue)>thresholdValue);
      break;
  }
  return false;
}
  

