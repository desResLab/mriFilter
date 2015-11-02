#include <math.h>
#include <boost/lexical_cast.hpp>

#include "mriScan.h"
#include "mriStructuredScan.h"
#include "mriThresholdCriteria.h"
#include "schMessages.h"
#include "mriUtils.h"

// ========================================
// PERFORM CONVOLUTION WITH GAUSSIAN KERNEL
// ========================================
double ConvoluteWithGaussianKernel(MRIIntVec cellNeighbors,MRIDoubleVec neighValues,MRIDoubleVec gaussKernelVector){
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
void MRIStructuredScan::GetStructuredNeighbourCells(int centreCell,int order,MRIThresholdCriteria* threshold,MRIIntVec& cellNeighbors){
 int centreCellCoords[3];
 int nextCell = 0;
 double cellQty = 0.0;
 cellNeighbors.clear();
 if (order == 0){
   cellNeighbors.push_back(centreCell);
   return;
 }else{
   // Get Coords in current cell
   MapIndexToCoords(centreCell,centreCellCoords);
   for(int i=-order;i<=order;i++){
     for(int j=-order;j<=order;j++){
       for(int k=-order;k<=order;k++){
         if((((centreCellCoords[0]+i) > -1)&&((centreCellCoords[0]+i) < cellTotals[0])) &&
            (((centreCellCoords[1]+j) > -1)&&((centreCellCoords[1]+j) < cellTotals[1])) &&
            (((centreCellCoords[2]+k) > -1)&&((centreCellCoords[2]+k) < cellTotals[2]))){
           nextCell = MapCoordsToIndex(centreCellCoords[0]+i,centreCellCoords[1]+j,centreCellCoords[2]+k);
         }else{
           nextCell = -1;
         }
         cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
         if(!threshold->MeetsCriteria(cellQty)){
           cellNeighbors.push_back(nextCell);
         }else{
           cellNeighbors.push_back(-1);
         }
       }
     }
   }
 }
}

// =================================
// GET CARTESIAN NEIGHBORS OF A CELL
// =================================
void MRIStructuredScan::GetCartesianNeighbourCells(int CurrentCell,std::vector<int> &cellNeighbors, bool addself){
  int* coords = new int[3];
  cellNeighbors.clear();
  if(addself){
    cellNeighbors.push_back(CurrentCell);
  }
  //Get The Coordinates of the Current Cell
  MapIndexToCoords(CurrentCell,coords);
  // Get Neighbor
  // coords[0]
  if ((coords[0]-1)>=0){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0]-1,coords[1],coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[1]
  if((coords[0]+1)<cellTotals[0]){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0]+1,coords[1],coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[2]
  if((coords[1]-1)>=0){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0],coords[1]-1,coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[3]
  if((coords[1]+1)<cellTotals[1]){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0],coords[1]+1,coords[2]));
  }else{
    cellNeighbors.push_back(-1);
  }
  // coords[4]
  if((coords[2]-1)>=0){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0],coords[1],coords[2]-1));
  }else{
    cellNeighbors.push_back(-1);
  }
    // coords[5]
  if((coords[2]+1)<cellTotals[2]){
    cellNeighbors.push_back(MapCoordsToIndex(coords[0],coords[1],coords[2]+1));
  }else{
    cellNeighbors.push_back(-1);
  }
  // SCRAMBLE VECTOR !!!
  //std::random_shuffle(&cellNeighbors[0], &cellNeighbors[5]);
  // DEALLOCATE
  delete [] coords;
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
void MRIStructuredScan::ApplyMedianFilter(int qtyID,int maxIt,int order,int filterType,MRIThresholdCriteria* threshold){
  MRIDoubleVec gaussKernelVector;
  if(filterType == kMedianFilter){
    WriteSchMessage(std::string("Applying Median Filter...\n"));
  }else if(filterType == kMeanFilter){
    WriteSchMessage(std::string("Applying Mean Filter...\n"));
  }else if(filterType == kGaussianFilter){
    WriteSchMessage(std::string("Applying Gaussian Filter...\n"));
    createGaussianKernel(order,gaussKernelVector);
    //for(int loopA=0;loopA<gaussKernelVector.size();loopA++){
    //  printf("%f\n",gaussKernelVector[loopA]);
    //}
  }else{
    throw MRIException("Error in ApplyMedianFilter: Invalid Filter Type.\n");
  }
  double* tempVec = new double[totalCellPoints];
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
    for(int loopA=0;loopA<totalCellPoints;loopA++){
      // Get Value in Current Cell
      centerCellValue = cellPoints[loopA].getQuantity(qtyID);
      // GET NEIGHBOURS
      GetStructuredNeighbourCells(loopA,order,threshold,neighbours);

      // CHECK IF THE CURRENT CELL IS TO BE PROCESSED
      cellQty = cellPoints[loopA].getQuantity(threshold->thresholdQty);
      if(!threshold->MeetsCriteria(cellQty)){

        // GET THE VALUES ON NEIGHBOR CELLS
        neighValues.clear();
        for(size_t loopB=0;loopB<neighbours.size();loopB++){
          currCell = neighbours[loopB];
          if(currCell>-1){
            // Get Quantity for Neighbor Cell
            currValue = cellPoints[currCell].getQuantity(qtyID);
            // Store value
            neighValues.push_back(currValue);
          }else if(filterType == kGaussianFilter){
            neighValues.push_back(0.0);
          }
        }

        // FIND MEDIAN VALUE
        if(filterType == kMedianFilter){
          if(neighValues.size() > 0){
            currMedian = MRIUtils::GetMedian(neighValues);
          }else{
            currMedian = centerCellValue;
          }
        }else if(filterType == kMeanFilter){
          if(neighValues.size() > 0){
            currMedian = MRIUtils::GetMean(neighValues);
          }else{
            currMedian = centerCellValue;
          }
        }else if(filterType == kGaussianFilter){
          currMedian = ConvoluteWithGaussianKernel(neighbours,neighValues,gaussKernelVector);
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
    for(int loopA=0;loopA<totalCellPoints;loopA++){
      cellPoints[loopA].setQuantity(qtyID,tempVec[loopA]);
    }
    // END OF ITERATION PRINT MAX ERROR
    WriteSchMessage(std::string("Iteration "+MRIUtils::IntToStr(loop0+1)+"; Max Error: "+MRIUtils::FloatToStr(maxError)+"\n"));
  }
  delete [] tempVec;
}

// =================================
// EVALUATE AVERAGE VELOCITY MODULUS
// =================================
double MRIScan::EvalAverageVelocityMod(){
  // Init Result
  double avVel = 0.0;
  double currentMod = 0.0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    currentMod = sqrt((cellPoints[loopA].velocity[0])*(cellPoints[loopA].velocity[0])+
                      (cellPoints[loopA].velocity[1])*(cellPoints[loopA].velocity[1])+
                      (cellPoints[loopA].velocity[2])*(cellPoints[loopA].velocity[2]));
    avVel += currentMod;
  }
  return ((double)currentMod/(double)totalCellPoints);
}

// ==================
// APPLY THRESHOLDING
// ==================
void MRIScan::ApplyThresholding(MRIThresholdCriteria* thresholdCriteria){
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("Applying Thresholding...\n"));
  // Init Number Of Filtered
  int numberOfFiltered = 0;
  // Apply Threshold
  double currentValue = 0.0;
  // If No Quantity then return
  if(thresholdCriteria->thresholdQty == kNoQuantity){
    return;
  }
  // Loop through the cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    switch(thresholdCriteria->thresholdQty){
      case kQtyPositionX:
        currentValue = cellPoints[loopA].position[0];
        break;
      case kQtyPositionY:
        currentValue = cellPoints[loopA].position[1];
        break;
      case kQtyPositionZ:
        currentValue = cellPoints[loopA].position[2];
        break;
      case kQtyConcentration:
        currentValue = cellPoints[loopA].concentration;
        break;
      case kQtyVelocityX:
        currentValue = cellPoints[loopA].velocity[0];
        break;
      case kQtyVelocityY:
        currentValue = cellPoints[loopA].velocity[1];
        break;
      case kQtyVelocityZ:
        currentValue = cellPoints[loopA].velocity[2];
        break;
    }
    if(thresholdCriteria->MeetsCriteria(currentValue)){
      numberOfFiltered++;
      cellPoints[loopA].velocity[0] = 0.0;
      cellPoints[loopA].velocity[1] = 0.0;
      cellPoints[loopA].velocity[2] = 0.0;
    }
  }
  WriteSchMessage(std::string("Cells Modified: "+MRIUtils::IntToStr(numberOfFiltered)+"\n"));
  WriteSchMessage(std::string("------------------------------------------------------------------\n"));
}

// =====================
// APPLY LAVISION FILTER
// =====================
void MRIScan::ApplySmoothingFilter(){
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
  std::vector<int> otherCells;
  // LOOP UNTIL CONVERGED
  while(!converged){
    // Increment Iteration Count
    itCount++;
    // Reset Max Divergence
    maxDivergence = 0.0;
    for(int loopA=0;loopA<totalCellPoints;loopA++){
      if(IsInnerCell(loopA)){
        // Eval Neighbours
        GetCartesianNeighbourCells(loopA,otherCells,false);
        // Get Velocities from Neighbors
        velXPlus =  cellPoints[otherCells[0]].velocity[0];
        velXMinus = cellPoints[otherCells[1]].velocity[0];
        velYPlus =  cellPoints[otherCells[2]].velocity[1];
        velYMinus = cellPoints[otherCells[3]].velocity[1];
        velZPlus =  cellPoints[otherCells[4]].velocity[2];
        velZMinus = cellPoints[otherCells[5]].velocity[2];
        // Compute Local Divergence
        localDivergence = (velXPlus - velXMinus) + (velYPlus - velYMinus) + (velZPlus - velZMinus);
        // Store Max Value
        if(fabs(localDivergence)>maxDivergence) maxDivergence = fabs(localDivergence);
        // Spread The Value Of Divergence
        cellPoints[otherCells[0]].velocity[0] -= (1.0/6.0) * localDivergence;
        cellPoints[otherCells[1]].velocity[0] += (1.0/6.0) * localDivergence;
        cellPoints[otherCells[2]].velocity[1] -= (1.0/6.0) * localDivergence;
        cellPoints[otherCells[3]].velocity[1] += (1.0/6.0) * localDivergence;
        cellPoints[otherCells[4]].velocity[2] -= (1.0/6.0) * localDivergence;
        cellPoints[otherCells[5]].velocity[2] += (1.0/6.0) * localDivergence;
      }
    }
    WriteSchMessage("It: "+MRIUtils::IntToStr(itCount)+";Max Div: "+MRIUtils::FloatToStr(maxDivergence)+"\n");
    // Check Convergence
    converged = (maxDivergence<1.0e-4);
  }
}

// ====================
// APPLY GAUSSIAN NOISE
// ====================
void MRIScan::ApplyGaussianNoise(double stDev){
  // Multiply by the velocity module
  stDev = stDev * maxVelModule / 100.0;
  // Apply Gaussian Noise
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    // Loop On Cells
    for(int loopB=0;loopB<totalCellPoints;loopB++){
      cellPoints[loopB].velocity[loopA] = cellPoints[loopB].velocity[loopA] + MRIUtils::GenerateStandardGaussian(stDev);
    }
  }
}
