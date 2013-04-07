#include <math.h>
#include <boost/lexical_cast.hpp>

#include "mriScan.h"
#include "mriThresholdCriteria.h"
#include "schMessages.h"
#include "mriUtils.h"

// EVALUATE AVERAGE VELOCITY MODULUS
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

// APPLY THRESHOLDING
void MRIScan::ApplyThresholding(MRIThresholdCriteria thresholdCriteria){
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("APPLYING THRESHOLD - Step : "+MRIUtils::FloatToStr(scanTime)+" ---\n"));
  // Init Number Of Filtered
  int numberOfFiltered = 0;
  // Apply Threshold
  double currentValue = 0.0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    switch(thresholdCriteria.thresholdQty){
      case 0: 
        currentValue = cellPoints[loopA].position[0];
        break;
      case 1: 
        currentValue = cellPoints[loopA].position[1];
        break;
      case 2: 
        currentValue = cellPoints[loopA].position[2];
        break;
      case 3: 
        currentValue = cellPoints[loopA].concentration;
        break;
      case 4: 
        currentValue = cellPoints[loopA].velocity[0];
        break;
      case 5: 
        currentValue = cellPoints[loopA].velocity[1];
        break;
      case 6: 
        currentValue = cellPoints[loopA].velocity[2];
        break;
    }
    if(thresholdCriteria.MeetsCriteria(currentValue)){
      numberOfFiltered++;
      cellPoints[loopA].velocity[0] = 0.0;
      cellPoints[loopA].velocity[1] = 0.0;
      cellPoints[loopA].velocity[2] = 0.0;
    }
  }
  WriteSchMessage(std::string("Cells Modified: "+MRIUtils::IntToStr(numberOfFiltered)+"\n"));
  WriteSchMessage(std::string("------------------------------------------------------------------\n"));
}

// APPLY LAVISION FILTER
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
  int* otherCells = new int[k3DNeighbors];
  // LOOP UNTIL CONVERGED
  while(!converged){
    // Increment Iteration Count
    itCount++;
    // Reset Max Divergence
    maxDivergence = 0.0;
    for(int loopA=0;loopA<totalCellPoints;loopA++){
      if(IsInnerCell(loopA)){
        // Eval Neighbours
        GetNeighbourCells(loopA,otherCells);
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
  // Deallocate
  delete [] otherCells;
}

// APPLY GAUSSIAN NOISE
void MRIScan::ApplyGaussianNoise(double stDev){
  // Apply Gaussian Noise
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    // Loop On Cells
    for(int loopB=0;loopB<totalCellPoints;loopB++){
      cellPoints[loopB].velocity[loopA] = cellPoints[loopB].velocity[loopA] + MRIUtils::GenerateStandardGaussian(stDev);
    }
  }
}
