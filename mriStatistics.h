#ifndef MRISTATISTICS_H
#define MRISTATISTICS_H

#include "mriException.h"
#include "mriScan.h"
#include "mriSequence.h"

namespace MRIStatistics
{

void NormalizeBinArray(int size, double* binArray,double currInterval){
  double sum = 0.0;
  // Compute Summation
  for(int loopA=0;loopA<size;loopA++){
    sum += binArray[loopA];
  }
  // Exit if zero sum
  if (fabs(sum)<kMathZero){
    return;
  }
  // Normalize
  for(int loopA=0;loopA<size;loopA++){
    binArray[loopA] /= (sum*currInterval);
  }  
}

// Assign to BIN
void AssignToBin(double currValue, int numberOfBins, double* binMin, double* binMax, double* binArray){
  bool found = false;
  int count = 0;
  bool isMoreThanMin = false;
  bool isLessThanMax = false;
  while ((!found)&&(count<numberOfBins)){
    if (fabs(currValue-binMin[0])<kMathZero){
      isMoreThanMin = (currValue >= binMin[count] - kMathZero);
    }else{
      isMoreThanMin = (currValue > binMin[count]);
    }
    if (fabs(currValue-binMin[numberOfBins-1])<kMathZero){
      isLessThanMax = (currValue <= binMax[count]);
    }else{
      isLessThanMax = (currValue <= binMax[count] + kMathZero);
    }
    found = (isMoreThanMin)&&(isLessThanMax);
    // Update
    if (!found){
      count++;
    }
  }
  if (found){
    // Increase Bin Count
    binArray[count] = binArray[count] + 1.0;
  }else{
    throw MRIStatisticsException("Error: Value Cannot fit in Bin.\n");
  } 
}

// FORM BIN LIMITS FOR SINGLE SCAN
void FormBinLimits(MRIScan* scan, int pdfQuantity, double &currInterval, double* limitBox, int numberOfBins, double* binMin, double* binMax, double* binCenter){
    // Initialize Limits
  double  minRange = std::numeric_limits<double>::max();
  double  maxRange = -std::numeric_limits<double>::max();
  double* cellCoord = NULL;
  double  otherQuantity = 0.0;
  double  refQuantity = 0.0;
  double  currValue = 0.0;
  for(int loopA=0;loopA<scan->totalCellPoints;loopA++){
    // Get Cell Coords
    cellCoord = scan->cellPoints[loopA].position;
    // Get quantity
    currValue = scan->cellPoints[loopA].getQuantity(pdfQuantity);
    // Check If Within the Bin 
    if (MRIUtils::IsPointInsideBox(cellCoord[0],cellCoord[1],cellCoord[2],limitBox)){
      // Assign Values
      if(currValue>maxRange) maxRange = currValue;
      if(currValue<minRange) minRange = currValue;
    }
  }
  // If minRange and maxRange are the same than add something
  if (fabs(maxRange - minRange)<kMathZero){
    minRange = minRange - 1.0;
    maxRange = maxRange + 1.0;
  }
  // Fill the bin arrays
  double currPtr = minRange;
  currInterval = ((maxRange - minRange)/(double)numberOfBins);
  for(int loopA=0;loopA<numberOfBins;loopA++){
    binMin[loopA] = currPtr;
    binMax[loopA] = currPtr + currInterval;
    binCenter[loopA] = 0.5*(binMin[loopA] + binMax[loopA]);
    // Update
    currPtr += currInterval;
  }
}

// FORM BIN LIMITS
void FormDifferenceBinLimits(MRIScan* scanOther, MRIScan* scanRef, int pdfQuantity, double &currInterval, double* limitBox, int numberOfBins, double* binMin, double* binMax, double* binCenter){
  // Initialize Limits
  double minRange = std::numeric_limits<double>::max();
  double maxRange = -std::numeric_limits<double>::max();
  double* cellCoord = NULL;
  double otherQuantity = 0.0;
  double refQuantity = 0.0;
  double currValue = 0.0;
  for(int loopA=0;loopA<scanRef->totalCellPoints;loopA++){
    // Get Cell Coords
    cellCoord = scanRef->cellPoints[loopA].position;
    // Get quantity
    otherQuantity = scanOther->cellPoints[loopA].getQuantity(pdfQuantity);
    refQuantity = scanRef->cellPoints[loopA].getQuantity(pdfQuantity);
    // Get Value
    currValue = (otherQuantity - refQuantity);
    // Check If Within the Bin 
    if (MRIUtils::IsPointInsideBox(cellCoord[0],cellCoord[1],cellCoord[2],limitBox)){
      // Assign Values
      if(currValue>maxRange) maxRange = currValue;
      if(currValue<minRange) minRange = currValue;
    }
  }
  // If minRange and maxRange are the same than add something
  if (fabs(maxRange - minRange)<kMathZero){
    minRange = minRange - 1.0;
    maxRange = maxRange + 1.0;
  }
  // Fill the bin arrays
  currInterval = ((maxRange - minRange)/(double)numberOfBins);
  double currPtr = minRange;
  for(int loopA=0;loopA<numberOfBins;loopA++){
    binMin[loopA] = currPtr;
    binMax[loopA] = currPtr + currInterval;
    binCenter[loopA] = 0.5*(binMin[loopA] + binMax[loopA]);
    // Update
    currPtr += currInterval;
  }
}

// Eval The PDF of a Single Scan
void EvalScanPDF(MRIScan* scan, int pdfQuantity, int numberOfBins, bool useBox, double* limitBox,double* binCenter, double* binArray){
  // Allocate Quantities
  double* binMin = new double[numberOfBins];
  double* binMax = new double[numberOfBins];
  // Form Bin 
  double currInterval = 0.0;
  FormBinLimits(scan,pdfQuantity,currInterval,limitBox,numberOfBins,binMin,binMax,binCenter);
  // Initialize binArray
  for(int loopA=0;loopA<numberOfBins;loopA++){
    binArray[loopA] = 0.0;
  }
  // Loop through the points
  double* cellCoord = NULL;
  double currValue = 0.0;
  for(int loopA=0;loopA<scan->totalCellPoints;loopA++){
    // Get Cell Coords
    cellCoord = scan->cellPoints[loopA].position;
    // Get quantity
    currValue = scan->cellPoints[loopA].getQuantity(pdfQuantity);
    // Assign Value to Bin
    if (MRIUtils::IsPointInsideBox(cellCoord[0],cellCoord[1],cellCoord[2],limitBox)&&fabs(currValue)>1.2e-2){
      // COMPLETE
      AssignToBin(currValue,numberOfBins,binMin,binMax,binArray);
    }
  }
  // Normalize
  NormalizeBinArray(numberOfBins,binArray,currInterval);
  // DeAllocate
  delete [] binMin;
  delete [] binMax;
}

// Eval The Difference PDF of Scans
void EvalScanDifferencePDF(MRISequence* sequence, int otherScan, int refScan, const int pdfQuantity, int numberOfBins, bool useBox, double* limitBox, double* binCenters, double* binArray){
  // Get The Scans out of the sequence
  MRIScan* scanOther = sequence->GetScan(otherScan);
  MRIScan* scanRef = sequence->GetScan(refScan);
  // Allocate Quantities
  double* binMin = new double[numberOfBins];
  double* binMax = new double[numberOfBins];
  // Form Bin 
  double currInterval = 0.0;
  FormDifferenceBinLimits(scanOther,scanRef,pdfQuantity,currInterval,limitBox,numberOfBins,binMin,binMax,binCenters);
  // Intialize bins
  for(int loopA=0;loopA<numberOfBins;loopA++){
    binArray[loopA] = 0.0;
  }
  // Loop through the points
  double sourceValue = 0.0;
  double refValue = 0.0;
  double* cellCoord = NULL;
  double otherQuantity = 0.0;
  double refQuantity = 0.0;
  double currValue = 0.0;
  for(int loopA=0;loopA<scanRef->totalCellPoints;loopA++){
    // Get Cell Coords
    cellCoord = scanRef->cellPoints[loopA].position;
    // Get quantity
    otherQuantity = scanOther->cellPoints[loopA].getQuantity(pdfQuantity);
    refQuantity = scanRef->cellPoints[loopA].getQuantity(pdfQuantity);
    // Get Value
    currValue = (otherQuantity - refQuantity);
    // Assign Value to Bin
    if (MRIUtils::IsPointInsideBox(cellCoord[0],cellCoord[1],cellCoord[2],limitBox)){
      // COMPLETE
      AssignToBin(currValue,numberOfBins,binMin,binMax,binArray);
    }
  }
  // Normalize
  NormalizeBinArray(numberOfBins,binArray,currInterval);
  // DeAllocate
  delete [] binMin;
  delete [] binMax;
}

} // MRIStatistics

#endif // MRISTATISTICS_H
