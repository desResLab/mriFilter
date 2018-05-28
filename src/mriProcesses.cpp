#include <math.h>
#include <boost/lexical_cast.hpp>
#include "mriScan.h"
#include "mriThresholdCriteria.h"
#include "mriUtils.h"

using namespace std;

// ========================================
// PERFORM CONVOLUTION WITH GAUSSIAN KERNEL
// ========================================
double convoluteWithGaussianKernel(MRIIntVec cellNeighbors,MRIDoubleVec neighValues,MRIDoubleVec gaussKernelVector){
  // Check Compatibility
  //printf("cellNeig: %d, neighVals: %d, gaussKernel: %d\n",cellNeighbors.size(),neighValues.size(),gaussKernelVector.size());
  //getchar();
  if((cellNeighbors.size() != neighValues.size())||(gaussKernelVector.size() != neighValues.size())){
    throw MRIException("ERROR: Incompatible size in ConvoluteWithGaussianKernel.\n");
  }
  double res = 0.0;
  double sumValid = 0.0;
  for(size_t loopA=0;loopA<cellNeighbors.size();loopA++){
    if(cellNeighbors[loopA]>-1){
      res += neighValues[loopA] * gaussKernelVector[loopA];
      sumValid += gaussKernelVector[loopA];
    }
  }
  return res/sumValid;
}

// ================================
// GET NEIGHBORS IN STRUCTURED GRID
// ================================
void MRIScan::getStructuredNeighbourCells(int centreCell,int order,MRIThresholdCriteria* threshold,MRIIntVec& cellNeighbors){
 MRIIntVec centreCellCoords(3);
 int nextCell = 0;
 double cellQty = 0.0;
 cellNeighbors.clear();
 if (order == 0){
   cellNeighbors.push_back(centreCell);
   return;
 }else{

   // Get Coords in current cell
   topology->mapIndexToCoords(centreCell,centreCellCoords);
   for(int i=-order;i<=order;i++){
     for(int j=-order;j<=order;j++){
       for(int k=-order;k<=order;k++){
         if((((centreCellCoords[0]+i) > -1)&&((centreCellCoords[0]+i) < topology->cellTotals[0])) &&
            (((centreCellCoords[1]+j) > -1)&&((centreCellCoords[1]+j) < topology->cellTotals[1])) &&
            (((centreCellCoords[2]+k) > -1)&&((centreCellCoords[2]+k) < topology->cellTotals[2]))){
           nextCell = topology->mapCoordsToIndex(centreCellCoords[0] + i, centreCellCoords[1] + j, centreCellCoords[2] + k);
         }else{
           nextCell = -1;
         }
         if(nextCell>-1){
           cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
           if(!threshold->meetsCriteria(cellQty)){
             cellNeighbors.push_back(nextCell);
           }else{
             cellNeighbors.push_back(-1);
           }
         }else{
           cellNeighbors.push_back(-1);
         }
       }
     }
   }
 }
}


// CREATE GAUSSIAN CONVOLUTION KERNEL
void createGaussianKernel(int order,MRIDoubleVec& kernel){
  kernel.clear();
  double sigma = order;
  double currentValue = 0.0;
  double sum = 0.0;
  for(int i=-order;i<=order;i++){
    for(int j=-order;j<=order;j++){
      for(int k=-order;k<=order;k++){
        currentValue = exp(-(i*i + j*j + k*k)/(2*sigma*sigma));
        kernel.push_back(currentValue);
        sum += currentValue;
      }
    }
  }
  // Rescale to 1.0
  for(size_t loopA=0;loopA<kernel.size();loopA++){
    kernel[loopA] = kernel[loopA]/sum;
  }
}

// ===================
// APPLY MEDIAN FILTER
// ===================
void MRIScan::applyMedianFilter(int qtyID,int maxIt,int order,int filterType,MRIThresholdCriteria* threshold){
  MRIDoubleVec gaussKernelVector;
  if(filterType == kMedianFilter){
    writeSchMessage(std::string("Applying Median Filter...\n"));
  }else if(filterType == kMeanFilter){
    writeSchMessage(std::string("Applying Mean Filter...\n"));
  }else if(filterType == kGaussianFilter){
    writeSchMessage(std::string("Applying Gaussian Filter...\n"));    
    createGaussianKernel(order,gaussKernelVector);
  }else{
    throw MRIException("Error in ApplyMedianFilter: Invalid Filter Type.\n");
  }
  MRIDoubleVec tempVec(topology->totalCells);
  double currValue = 0.0;
  double currMedian = 0.0;
  double currError = 0.0;
  double centerCellValue = 0.0;
  int currCell = 0;
  double cellQty = 0.0;
  MRIIntVec neighbours;
  MRIDoubleVec neighValues;
  // PERFORM ITERATIONS
  for(int loop0=0;loop0<maxIt;loop0++){
    // LOOP ON ALL CELLS
    double maxError = 0.0;
    for(int loopA=0;loopA<topology->totalCells;loopA++){

      // Get Value in Current Cell
      centerCellValue = cells[loopA].getQuantity(qtyID);

      // GET NEIGHBOURS
      getStructuredNeighbourCells(loopA,order,threshold,neighbours);      

      // CHECK IF THE CURRENT CELL IS TO BE PROCESSED
      cellQty = cells[loopA].getQuantity(threshold->thresholdQty);
      if(!threshold->meetsCriteria(cellQty)){

        // GET THE VALUES ON NEIGHBOR CELLS
        neighValues.clear();
        for(size_t loopB=0;loopB<neighbours.size();loopB++){
          currCell = neighbours[loopB];
          if(currCell>-1){
            // Get Quantity for Neighbor Cell
            currValue = cells[currCell].getQuantity(qtyID);
            // Store value
            neighValues.push_back(currValue);
          }else if(filterType == kGaussianFilter){
            neighValues.push_back(0.0);
          }
        }

        // FIND MEDIAN VALUE
        if(filterType == kMedianFilter){
          if(neighValues.size() > 0){
            currMedian = MRIUtils::getMedian(neighValues);
          }else{
            currMedian = centerCellValue;
          }
        }else if(filterType == kMeanFilter){
          if(neighValues.size() > 0){
            currMedian = MRIUtils::getMean(neighValues);
          }else{
            currMedian = centerCellValue;
          }
        }else if(filterType == kGaussianFilter){
          currMedian = convoluteWithGaussianKernel(neighbours,neighValues,gaussKernelVector);
        }
      }else{
        currMedian = centerCellValue;
      }
      // EVAL CHANGE
      currError = fabs(((currMedian-centerCellValue)/(double)maxVelModule)*100.0);
      if(currError>maxError){
        maxError = currError;
      }
      // ASSIGN MEDIAN
      tempVec[loopA] = currMedian;
    }
    // ASSIGN VALUES
    for(int loopA=0;loopA<topology->totalCells;loopA++){
      cells[loopA].setQuantity(qtyID,tempVec[loopA]);
    }
    // END OF ITERATION PRINT MAX ERROR
    writeSchMessage(string("Iteration "+MRIUtils::intToStr(loop0+1)+"; Max Error: "+MRIUtils::floatToStr(maxError)+"\n"));
  }
}

// =================================
// EVALUATE AVERAGE VELOCITY MODULUS
// =================================
double MRIScan::evalAverageVelocityMod(){
  // Init Result
  double avVel = 0.0;
  double currentMod = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    currentMod = sqrt((cells[loopA].velocity[0])*(cells[loopA].velocity[0])+
                      (cells[loopA].velocity[1])*(cells[loopA].velocity[1])+
                      (cells[loopA].velocity[2])*(cells[loopA].velocity[2]));
    avVel += currentMod;
  }
  return ((double)currentMod/(double)topology->totalCells);
}

// ==================
// APPLY THRESHOLDING
// ==================
void MRIScan::applyThresholding(MRIThresholdCriteria* thresholdCriteria){
  writeSchMessage(std::string("\n"));
  writeSchMessage(std::string("Applying Thresholding...\n"));
  // Init Number Of Filtered
  int numberOfFiltered = 0;
  // Apply Threshold
  double currentValue = 0.0;
  // If No Quantity then return
  if(thresholdCriteria->thresholdQty == kNoQuantity){
    return;
  }
  // Loop through the cells
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    switch(thresholdCriteria->thresholdQty){
      case kQtyConcentration:
        currentValue = cells[loopA].concentration;
        break;
      case kQtyVelocityX:
        currentValue = cells[loopA].velocity[0];
        break;
      case kQtyVelocityY:
        currentValue = cells[loopA].velocity[1];
        break;
      case kQtyVelocityZ:
        currentValue = cells[loopA].velocity[2];
        break;
    }
    if(thresholdCriteria->meetsCriteria(currentValue)){
      numberOfFiltered++;
      cells[loopA].velocity[0] = 0.0;
      cells[loopA].velocity[1] = 0.0;
      cells[loopA].velocity[2] = 0.0;
    }
  }
  writeSchMessage(string("Cells Modified: "+MRIUtils::intToStr(numberOfFiltered)+"\n"));
  writeSchMessage(string("------------------------------------------------------------------\n"));
}

// =====================
// APPLY LAVISION FILTER
// =====================
void MRIScan::applySmoothingFilter(){
  bool converged = false;
  int itCount = 0;
  double maxDivergence = 0.0;
  // Get Velocities from Neighbors
  double velXPlus = 0.0;
  double velXMinus = 0.0;
  double velYPlus = 0.0;
  double velYMinus = 0.0;
  double velZPlus = 0.0;
  double velZMinus = 0.0;
  // Local Divergence 
  double localDivergence = 0.0;
  MRIIntVec otherCells;
  // LOOP UNTIL CONVERGED
  while(!converged){
    // Increment Iteration Count
    itCount++;
    // Reset Max Divergence
    maxDivergence = 0.0;
    for(int loopA=0;loopA<topology->totalCells;loopA++){
      if(isInnerCell(loopA)){
        // Eval Neighbours
        getCartesianNeighbourCells(loopA,otherCells,false);
        // Get Velocities from Neighbors
        velXPlus =  cells[otherCells[0]].velocity[0];
        velXMinus = cells[otherCells[1]].velocity[0];
        velYPlus =  cells[otherCells[2]].velocity[1];
        velYMinus = cells[otherCells[3]].velocity[1];
        velZPlus =  cells[otherCells[4]].velocity[2];
        velZMinus = cells[otherCells[5]].velocity[2];
        // Compute Local Divergence
        localDivergence = (velXPlus - velXMinus) + (velYPlus - velYMinus) + (velZPlus - velZMinus);
        // Store Max Value
        if(fabs(localDivergence)>maxDivergence) maxDivergence = fabs(localDivergence);
        // Spread The Value Of Divergence
        cells[otherCells[0]].velocity[0] -= (1.0/6.0) * localDivergence;
        cells[otherCells[1]].velocity[0] += (1.0/6.0) * localDivergence;
        cells[otherCells[2]].velocity[1] -= (1.0/6.0) * localDivergence;
        cells[otherCells[3]].velocity[1] += (1.0/6.0) * localDivergence;
        cells[otherCells[4]].velocity[2] -= (1.0/6.0) * localDivergence;
        cells[otherCells[5]].velocity[2] += (1.0/6.0) * localDivergence;
      }
    }
    writeSchMessage("It: "+MRIUtils::intToStr(itCount)+";Max Div: "+MRIUtils::floatToStr(maxDivergence)+"\n");
    // Check Convergence
    converged = (maxDivergence<1.0e-4);
  }
}

// ====================
// APPLY GAUSSIAN NOISE
// ====================
void MRIScan::ApplyGaussianNoise(double stDev, double seed){
  // Multiply by the velocity module
  stDev = stDev * maxVelModule / 100.0;
  MRIUtils::SetSeed(seed);
  // Apply Gaussian Noise
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    // Loop On Cells
    for(int loopB=0;loopB<topology->totalCells;loopB++){
      cells[loopB].velocity[loopA] = cells[loopB].velocity[loopA] + MRIUtils::generateStandardGaussian(stDev);
    }
  }
}
