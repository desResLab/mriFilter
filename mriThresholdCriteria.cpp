#include <math.h>
#include "mriThresholdCriteria.h"
#include "mriConstants.h"

MRIThresholdCriteria::MRIThresholdCriteria(int type, int qty, double value){
  thresholdType = type;
  thresholdQty = qty;
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
  

