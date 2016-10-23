#include <math.h>

#include "mriScan.h"
#include "mriConstants.h"
#include "mriSequence.h"
#include "mriUtils.h"
#include "mriException.h"

// ===========================================
// COMPUTE PRESSURE GRADIENTS FOR GENERIC CELL
// ===========================================
void MRIScan::evalCellPressureGradients(int currentCell,
                                        const MRIDoubleVec& timeDeriv, 
                                        const MRIDoubleMat& firstDerivs, 
                                        const MRIDoubleMat& secondDerivs,
                                        const MRIDoubleMat& ReynoldsStressGrad,
                                        MRIDoubleVec& pressureGrad){

  // SET PARAMETER
  int pressureGradientType = 3;

  // INIT
  double gravityTerm = 0.0;
  double viscousTerm = 0.0;
  double convectiveTerm = 0.0;
  double ReynoldsStressTerm = 0.0;
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    // Eval Viscous Term
    viscousTerm = viscosity * (secondDerivs[0][loopA] + secondDerivs[1][loopA] + secondDerivs[2][loopA]);
    // Eval Convective Term
    convectiveTerm = density * (timeDeriv[loopA] +
                                cells[currentCell].velocity[0] * firstDerivs[0][loopA]+
                                cells[currentCell].velocity[1] * firstDerivs[1][loopA]+
                                cells[currentCell].velocity[2] * firstDerivs[2][loopA]);
    // Eval Reynolds Stress Gradients
    if(ReynoldsStressGrad.size() > 0){
      ReynoldsStressTerm = density * (ReynoldsStressGrad[0][loopA] + ReynoldsStressGrad[1][loopA] + ReynoldsStressGrad[2][loopA]);
    }

    // Eval Pressure Gradient
    switch(pressureGradientType){
      case 0:
        pressureGrad[loopA] = viscousTerm;
        break;
      case 1:
        pressureGrad[loopA] = - convectiveTerm;
        break;
      case 2:
        pressureGrad[loopA] = - ReynoldsStressTerm;
        break;
      case 3:
        pressureGrad[loopA] = viscousTerm + gravityTerm - convectiveTerm; // - ReynoldsStressTerm;
        //pressureGrad[loopA] = - convectiveTerm;
        break;
    }
  }
}

// Find Derivatives using a Parabola
double findParabolicDeriv(double y1, double y2, double y3, double h1, double h2, bool isFirstPoint){
  double bValue = ((y3-y2)-(y1-y2)*(h2/h1)*(h2/h1))/((h2*h2)/(h1)+h2);
  double aValue = (y1-y2+bValue*h1)/(h1*h1);
  if (isFirstPoint){
    return (bValue - 2.0*aValue*h1);
  }else{
    return (bValue + 2.0*aValue*h2);
  }
}

// Eval Time Derivatives: TODO: Complete for Transient Case
void MRISequence::evalTimeDerivs(int currentScan, int currentCell, MRIDoubleVec& timeDeriv){
  // Eval Difference Formulae
  if (totalScans>1){
    // multiple Scan
    bool isFirstScan = (currentScan == 0);
    bool isLastScan = (currentScan == (totalScans-1));

    MRIDoubleVec currVel(3);
    MRIDoubleVec prevVel(3);
    MRIDoubleVec nextVel(3);
    MRIDoubleVec nextNextVel(3);

    if (isFirstScan){
      if (isCyclic){
        double deltaT2 = sequence[currentScan+1]->scanTime - sequence[currentScan]->scanTime;
        // Assume Same Time Step
        double deltaT1 = deltaT2;
        // Get The Velocities
        MRIDoubleVec firstPart(3,0.0);
        MRIDoubleVec secondPart(3,0.0);
        //
        currVel[0] = sequence[currentScan]->cells[currentCell].velocity[0];
        currVel[1] = sequence[currentScan]->cells[currentCell].velocity[1];
        currVel[2] = sequence[currentScan]->cells[currentCell].velocity[2];
        // 
        prevVel[0] = sequence[totalScans-1]->cells[currentCell].velocity[0];
        prevVel[1] = sequence[totalScans-1]->cells[currentCell].velocity[1];
        prevVel[2] = sequence[totalScans-1]->cells[currentCell].velocity[2];
        //
        nextVel[0] = sequence[currentScan+1]->cells[currentCell].velocity[0];
        nextVel[1] = sequence[currentScan+1]->cells[currentCell].velocity[1];
        nextVel[2] = sequence[currentScan+1]->cells[currentCell].velocity[2];
        // Use Central Differencing
        // First Part
        firstPart[0] = (currVel[0] - prevVel[0])/(2.0*deltaT1);
        firstPart[1] = (currVel[1] - prevVel[1])/(2.0*deltaT1);
        firstPart[2] = (currVel[2] - prevVel[2])/(2.0*deltaT1);      
        // Second Part
        secondPart[0] = (nextVel[0] - currVel[0])/(2.0*deltaT2);
        secondPart[1] = (nextVel[1] - currVel[1])/(2.0*deltaT2);
        secondPart[2] = (nextVel[2] - currVel[2])/(2.0*deltaT2);
        // Final Time Derivative
        timeDeriv[0] = firstPart[0] + secondPart[0];
        timeDeriv[1] = firstPart[1] + secondPart[1];
        timeDeriv[2] = firstPart[2] + secondPart[2];
      }else if(totalScans>2){
        // Parabola Derivatives
        double deltaT1 = sequence[currentScan+1]->scanTime - sequence[currentScan]->scanTime;
        double deltaT2 = sequence[currentScan+2]->scanTime - sequence[currentScan+1]->scanTime;
        // Get The Velocities
        currVel[0] = sequence[currentScan]->cells[currentCell].velocity[0];
        currVel[1] = sequence[currentScan]->cells[currentCell].velocity[1];
        currVel[2] = sequence[currentScan]->cells[currentCell].velocity[2];
        // 
        nextVel[0] = sequence[currentScan+1]->cells[currentCell].velocity[0];
        nextVel[1] = sequence[currentScan+1]->cells[currentCell].velocity[1];
        nextVel[2] = sequence[currentScan+1]->cells[currentCell].velocity[2];
        //
        nextNextVel[0] = sequence[currentScan+2]->cells[currentCell].velocity[0];
        nextNextVel[1] = sequence[currentScan+2]->cells[currentCell].velocity[1];
        nextNextVel[2] = sequence[currentScan+2]->cells[currentCell].velocity[2];

        // Eval Time Derivatives
        timeDeriv[0] = findParabolicDeriv(currVel[0],nextVel[0],nextNextVel[0],deltaT1,deltaT2,true);
        timeDeriv[1] = findParabolicDeriv(currVel[1],nextVel[1],nextNextVel[1],deltaT1,deltaT2,true);
        timeDeriv[2] = findParabolicDeriv(currVel[2],nextVel[2],nextNextVel[2],deltaT1,deltaT2,true);
      }else{
        // Euler Formula
        double deltaT = sequence[currentScan+1]->scanTime - sequence[currentScan]->scanTime;
        // Use Simple Difference Scheme
        timeDeriv[0] = (sequence[currentScan+1]->cells[currentCell].velocity[0] - sequence[currentScan]->cells[currentCell].velocity[0])/(deltaT);
        timeDeriv[1] = (sequence[currentScan+1]->cells[currentCell].velocity[1] - sequence[currentScan]->cells[currentCell].velocity[1])/(deltaT);
        timeDeriv[2] = (sequence[currentScan+1]->cells[currentCell].velocity[2] - sequence[currentScan]->cells[currentCell].velocity[2])/(deltaT);        
      }
    }else if(isLastScan){
      if (isCyclic){
        double deltaT1 = sequence[currentScan]->scanTime - sequence[currentScan-1]->scanTime;
        // Assume Same Time Step
        double deltaT2 = deltaT1;
        // Get The Velocities
        MRIDoubleVec firstPart(3,0.0);
        MRIDoubleVec secondPart(3,0.0);
        //
        currVel[0] = sequence[currentScan]->cells[currentCell].velocity[0];
        currVel[1] = sequence[currentScan]->cells[currentCell].velocity[1];
        currVel[2] = sequence[currentScan]->cells[currentCell].velocity[2];
        //
        prevVel[0] = sequence[currentScan-1]->cells[currentCell].velocity[0];
        prevVel[1] = sequence[currentScan-1]->cells[currentCell].velocity[1];
        prevVel[2] = sequence[currentScan-1]->cells[currentCell].velocity[2];
        //
        nextVel[0] = sequence[0]->cells[currentCell].velocity[0];      
        nextVel[1] = sequence[0]->cells[currentCell].velocity[1];      
        nextVel[2] = sequence[0]->cells[currentCell].velocity[2];      
        // Use Central Differencing
        // First Part
        firstPart[0] = (currVel[0] - prevVel[0])/(2.0*deltaT1);
        firstPart[1] = (currVel[1] - prevVel[1])/(2.0*deltaT1);
        firstPart[2] = (currVel[2] - prevVel[2])/(2.0*deltaT1);      
        // Second Part
        secondPart[0] = (nextVel[0] - currVel[0])/(2.0*deltaT2);
        secondPart[1] = (nextVel[1] - currVel[1])/(2.0*deltaT2);
        secondPart[2] = (nextVel[2] - currVel[2])/(2.0*deltaT2);
        // Final Time Derivative
        timeDeriv[0] = firstPart[0] + secondPart[0];
        timeDeriv[1] = firstPart[1] + secondPart[1];
        timeDeriv[2] = firstPart[2] + secondPart[2];        
      }else if(totalScans>2){
        // Parabola Derivatives
        double deltaT2 = sequence[currentScan]->scanTime - sequence[currentScan-1]->scanTime;
        double deltaT1 = sequence[currentScan-1]->scanTime - sequence[currentScan-2]->scanTime;
        // Get The Velocities
        currVel[0] = sequence[currentScan-2]->cells[currentCell].velocity[0];
        currVel[1] = sequence[currentScan-2]->cells[currentCell].velocity[1];
        currVel[2] = sequence[currentScan-2]->cells[currentCell].velocity[2];
        //
        nextVel[0] = sequence[currentScan-1]->cells[currentCell].velocity[0];
        nextVel[1] = sequence[currentScan-1]->cells[currentCell].velocity[1];
        nextVel[2] = sequence[currentScan-1]->cells[currentCell].velocity[2];
        //
        nextNextVel[0] = sequence[currentScan]->cells[currentCell].velocity[0];
        nextNextVel[1] = sequence[currentScan]->cells[currentCell].velocity[1];
        nextNextVel[2] = sequence[currentScan]->cells[currentCell].velocity[2];
        // Eval Time Derivatives
        timeDeriv[0] = findParabolicDeriv(currVel[0],nextVel[0],nextNextVel[0],deltaT1,deltaT2,false);
        timeDeriv[1] = findParabolicDeriv(currVel[1],nextVel[1],nextNextVel[1],deltaT1,deltaT2,false);
        timeDeriv[2] = findParabolicDeriv(currVel[2],nextVel[2],nextNextVel[2],deltaT1,deltaT2,false);        
      }else{
        // Euler Formula
        double deltaT = sequence[currentScan]->scanTime - sequence[currentScan-1]->scanTime;
        timeDeriv[0] = (sequence[currentScan]->cells[currentCell].velocity[0] - sequence[currentScan-1]->cells[currentCell].velocity[0])/(deltaT);
        timeDeriv[1] = (sequence[currentScan]->cells[currentCell].velocity[1] - sequence[currentScan-1]->cells[currentCell].velocity[1])/(deltaT);
        timeDeriv[2] = (sequence[currentScan]->cells[currentCell].velocity[2] - sequence[currentScan-1]->cells[currentCell].velocity[2])/(deltaT);                
      }
    }else{
      MRIDoubleVec firstPart(3,0.0);
      MRIDoubleVec secondPart(3,0.0);
      double deltaT1 = sequence[currentScan]->scanTime - sequence[currentScan-1]->scanTime;
      double deltaT2 = sequence[currentScan+1]->scanTime - sequence[currentScan]->scanTime;
      // Get The Velocities
      currVel[0] = sequence[currentScan]->cells[currentCell].velocity[0];
      currVel[1] = sequence[currentScan]->cells[currentCell].velocity[1];
      currVel[2] = sequence[currentScan]->cells[currentCell].velocity[2];
      //
      prevVel[0] = sequence[currentScan-1]->cells[currentCell].velocity[0];
      prevVel[1] = sequence[currentScan-1]->cells[currentCell].velocity[1];
      prevVel[2] = sequence[currentScan-1]->cells[currentCell].velocity[2];
      //
      nextVel[0] = sequence[currentScan+1]->cells[currentCell].velocity[0];
      nextVel[1] = sequence[currentScan+1]->cells[currentCell].velocity[1];
      nextVel[2] = sequence[currentScan+1]->cells[currentCell].velocity[2];
      // Use Central Differencing
      // First Part
      firstPart[0] = (currVel[0] - prevVel[0])/(2.0*deltaT1);
      firstPart[1] = (currVel[1] - prevVel[1])/(2.0*deltaT1);
      firstPart[2] = (currVel[2] - prevVel[2])/(2.0*deltaT1);      
      // Second Part
      secondPart[0] = (nextVel[0] - currVel[0])/(2.0*deltaT2);
      secondPart[1] = (nextVel[1] - currVel[1])/(2.0*deltaT2);
      secondPart[2] = (nextVel[2] - currVel[2])/(2.0*deltaT2);
      // Final Time Derivative
      timeDeriv[0] = firstPart[0] + secondPart[0];
      timeDeriv[1] = firstPart[1] + secondPart[1];
      timeDeriv[2] = firstPart[2] + secondPart[2];
    }
  }else{
    // Only One Scan
    // DVX/DT
    timeDeriv[0] = 0.0;
    // DVY/DT
    timeDeriv[1] = 0.0;
    // DVZ/DT
    timeDeriv[2] = 0.0;
  }
}

// ==================================
// EVAL TIME DERIVATIVES FOR THE SCAN
// ==================================
void MRISequence::evalScanTimeDerivs(int currentScan,MRIDoubleMat& timeDeriv){
  timeDeriv.clear();
  MRIDoubleVec temp;
  MRIDoubleVec cellTimeDeriv(3);
  for(int loopA=0;loopA<sequence[currentScan]->topology->totalCells;loopA++){
    temp.clear();
    evalTimeDerivs(currentScan,loopA,cellTimeDeriv);
    temp.push_back(cellTimeDeriv[0]);
    temp.push_back(cellTimeDeriv[1]);
    temp.push_back(cellTimeDeriv[2]);
    timeDeriv.push_back(temp);
  }
}

// ================================
// EVAL REYNOLDS STRESS DERIVATIVES
// ================================
void MRISequence::evalScanReynoldsStressDerivs(int currentScan,MRIDoubleMat& reynoldsDeriv){
  reynoldsDeriv.clear();
  MRIDoubleVec temp;
  MRIDoubleMat rsg;
  rsg.resize(3);
  for(int loopA=0;loopA<3;loopA++){
    rsg[loopA].resize(3);
  }
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    temp.clear();
    sequence[currentScan]->evalReynoldsStressGradient(loopA, rsg);
    temp.push_back(rsg[0][0] + rsg[1][0] + rsg[2][0]);
    temp.push_back(rsg[0][1] + rsg[1][1] + rsg[2][1]);
    temp.push_back(rsg[0][2] + rsg[1][2] + rsg[2][2]);
    reynoldsDeriv.push_back(temp);
  }
}

// ==========================================
// EVAL FIRST AND SECOND DERIVATIVES IN SPACE
// ==========================================
void MRIScan::evalSpaceDerivs(int currentCell, MRIThresholdCriteria* threshold, MRIDoubleMat& firstDerivs, MRIDoubleMat& secondDerivs){
  // FirstDerivs
  // DVX/DX DVY/DX DVZ/DX
  // DVX/DY DVY/DY DVZ/DY
  // DVX/DZ DVY/DZ DVZ/DZ

  // Get quantity for threshold evaluation
  double cellQty = 0.0;

  // Map Index To Coords
  MRIIntVec currentCellCoords(kNumberOfDimensions);
  topology->mapIndexToCoords(currentCell,currentCellCoords);
  int firstCell,secondCell,thirdCell,fourthCell,nextCell;

  // Assemble Terms
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    switch(loopA){
      case 0:
        // X Deriv
        if((currentCellCoords[0]-1)<0){
          firstCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0]-1,currentCellCoords[1],currentCellCoords[2]);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            firstCell = -1;
          }else{
            firstCell = nextCell;
          }
        }
        if((currentCellCoords[0]+1)>(topology->cellTotals[0]-1)){
          secondCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0]+1,currentCellCoords[1],currentCellCoords[2]);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            secondCell = -1;
          }else{
            secondCell = nextCell;
          }
        }
        if((currentCellCoords[0]-2)<0){
          thirdCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0]-2,currentCellCoords[1],currentCellCoords[2]);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            thirdCell = -1;
          }else{
            thirdCell = nextCell;
          }
        }
        if((currentCellCoords[0]+2)>(topology->cellTotals[0]-1)){
          fourthCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0]+2,currentCellCoords[1],currentCellCoords[2]);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            fourthCell = -1;
          }else{
            fourthCell = nextCell;
          }
        }
        break;
      case 1:
        // Y Deriv
        if((currentCellCoords[1]-1)<0){
          firstCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]-1,currentCellCoords[2]);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            firstCell = -1;
          }else{
            firstCell = nextCell;
          }
        }
        if((currentCellCoords[1]+1)>(topology->cellTotals[1]-1)){
          secondCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]+1,currentCellCoords[2]);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            secondCell = -1;
          }else{
            secondCell = nextCell;
          }
        }
        if((currentCellCoords[1]-2)<0){
          thirdCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]-2,currentCellCoords[2]);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            thirdCell = -1;
          }else{
            thirdCell = nextCell;
          }
        }
        if((currentCellCoords[1]+2)>(topology->cellTotals[1]-1)){
          fourthCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]+2,currentCellCoords[2]);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            fourthCell = -1;
          }else{
            fourthCell = nextCell;
          }
        }
        break;
      case 2:
        // Z Deriv
        if((currentCellCoords[2]-1)<0){
          firstCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]-1);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            firstCell = -1;
          }else{
            firstCell = nextCell;
          }
        }
        if((currentCellCoords[2]+1)>(topology->cellTotals[2]-1)){
          secondCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]+1);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            secondCell = -1;
          }else{
            secondCell = nextCell;
          }
        }
        if((currentCellCoords[2]-2)<0){
          thirdCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]-2);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            thirdCell = -1;
          }else{
            thirdCell = nextCell;
          }
        }
        if((currentCellCoords[2]+2)>(topology->cellTotals[2]-1)){
          fourthCell = -1;
        }else{
          nextCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]+2);
          cellQty = cells[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->meetsCriteria(cellQty)){
            fourthCell = -1;
          }else{
            fourthCell = nextCell;
          }
        }
        break;
    }

    // Get Deltas
    double deltaMinus = 0.0;
    double deltaPlus = 0.0;
    double deltaPlusPlus = 0.0;
    double deltaMinusMinus = 0.0;
    if(firstCell>-1){
      deltaMinus = 0.5*(topology->cellLengths[loopA][currentCellCoords[loopA]] + topology->cellLengths[loopA][currentCellCoords[loopA]-1]);
    }else{
      deltaMinus = 0.0;
    }
    if(secondCell>-1){
      deltaPlus = 0.5*(topology->cellLengths[loopA][currentCellCoords[loopA]] + topology->cellLengths[loopA][currentCellCoords[loopA]+1]);
    }else{
      deltaPlus = 0.0;
    }
    if(thirdCell>-1){
      deltaMinusMinus = 0.5*(topology->cellLengths[loopA][currentCellCoords[loopA]-1] + topology->cellLengths[loopA][currentCellCoords[loopA]-2]);
    }else{
      deltaMinusMinus = 0.0;
    }
    if(fourthCell>-1){
      deltaPlusPlus = 0.5*(topology->cellLengths[loopA][currentCellCoords[loopA]+1] + topology->cellLengths[loopA][currentCellCoords[loopA]+2]);
    }else{
      deltaPlusPlus = 0.0;
    }

    // Find Components
    double currentVComponent = 0.0;
    double firstVComponent = 0.0;
    double secondVComponent = 0.0;
    double thirdVComponent = 0.0;
    double fourthVComponent = 0.0;
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      // First Component
      if(firstCell>-1){
        firstVComponent = cells[firstCell].velocity[loopB];
      }else{
        firstVComponent = 0.0;
      }
      // Second Component
      if(secondCell>-1){
        secondVComponent = cells[secondCell].velocity[loopB];
      }else{
        secondVComponent = 0.0;
      }
      // Third Component
      if(thirdCell>-1){
        thirdVComponent = cells[thirdCell].velocity[loopB];
      }else{
        thirdVComponent = 0.0;
      }
      // Fourth Component
      if(fourthCell>-1){
        fourthVComponent = cells[fourthCell].velocity[loopB];
      }else{
        fourthVComponent = 0.0;
      }

      // Check if the current Cell has a
      cellQty = cells[currentCell].getQuantity(threshold->thresholdQty);
      if(threshold->meetsCriteria(cellQty)){
        firstDerivs[loopA][loopB] = 0.0;
        secondDerivs[loopA][loopB] = 0.0;
      }else{
        // Current Component
        currentVComponent = cells[currentCell].velocity[loopB];

        // FIRST DERIVS
        if(firstCell<0){
          // Simple Euler Formula
          if(fabs(deltaPlus) > kMathZero){
            firstDerivs[loopA][loopB] = (secondVComponent-currentVComponent)/(deltaPlus);
          }else{
            firstDerivs[loopA][loopB] = 0.0;
          }
        }else if (secondCell<0){
          // Simple Euler Formula
          if(fabs(deltaMinus) > kMathZero){
            firstDerivs[loopA][loopB] = (currentVComponent-firstVComponent)/(deltaMinus);
          }else{
            firstDerivs[loopA][loopB] = 0.0;
          }
        }else if((firstCell>-1)&&(secondCell>-1)){
          // Central Difference Formula: CAREFULL: ONLY FIRST ORDER IF GRID SPACING VARIES SIGNIFICANTLY
          if((deltaPlus + deltaMinus) > kMathZero){
            firstDerivs[loopA][loopB] = (secondVComponent-firstVComponent)/(deltaPlus + deltaMinus);
          }else{
            firstDerivs[loopA][loopB] = 0.0;
          }
        }else{
          // Show Error Message
          throw MRIException("Error: Both First and Second Cells are Zero in EvalFirstSpaceDerivs");
        }

        // SECOND DERIVS
        if(firstCell<0){
          if((secondCell>-1)&&(fourthCell>-1)){
            secondDerivs[loopA][loopB] = (currentVComponent-2.0*secondVComponent+fourthVComponent)/(deltaPlus*deltaPlusPlus);
          }else{
            secondDerivs[loopA][loopB] = 0.0;
          }
        }else if (secondCell<0){
          if((firstCell>-1)&&(thirdCell>-1)){
            secondDerivs[loopA][loopB] = (thirdVComponent-2.0*firstVComponent+currentVComponent)/(deltaMinus*deltaMinusMinus);
          }else{
            secondDerivs[loopA][loopB] = 0.0;
          }
        }else if((firstCell>-1)&&(secondCell>-1)){
          // Central Difference Formula
          secondDerivs[loopA][loopB] = (secondVComponent-2.0*currentVComponent+firstVComponent)/(deltaPlus*deltaMinus);
        }else{
          // Show Error Message
          throw MRIException("Error: Both First and Second Cells are Zero in EvalFirstSpaceDerivs");
        }
      }
    }
  }
}

// ==========================================
// EVAL FIRST AND SECOND DERIVATIVES IN SPACE
// ==========================================
void MRIScan::evalSpaceGradient(int currentCell,int qtyID, MRIDoubleVec& gradient){
  // Map Index To Coords
  MRIIntVec currentCellCoords(kNumberOfDimensions);
  topology->mapIndexToCoords(currentCell,currentCellCoords);
  int firstCell,secondCell;
  // Assemble Terms
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    switch(loopA){
      case 0:
        // X Deriv
        if((currentCellCoords[0]-1)<0){
          firstCell = -1;
        }else{
          firstCell = topology->mapCoordsToIndex(currentCellCoords[0]-1,currentCellCoords[1],currentCellCoords[2]);
        }
        if((currentCellCoords[0]+1)>(topology->cellTotals[0]-1)){
          secondCell = -1;
        }else{
          secondCell = topology->mapCoordsToIndex(currentCellCoords[0]+1,currentCellCoords[1],currentCellCoords[2]);
        }
        break;
      case 1:
        // Y Deriv
        if((currentCellCoords[1]-1)<0){
          firstCell = -1;
        }else{
          firstCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]-1,currentCellCoords[2]);
        }
        if((currentCellCoords[1]+1)>(topology->cellTotals[1]-1)){
          secondCell = -1;
        }else{
          secondCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]+1,currentCellCoords[2]);
        }
        break;
      case 2:
        // Z Deriv
        if((currentCellCoords[2]-1)<0){
          firstCell = -1;
        }else{
          firstCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]-1);
        }
        if((currentCellCoords[2]+1)>(topology->cellTotals[2]-1)){
          secondCell = -1;
        }else{
          secondCell = topology->mapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]+1);
        }
        break;
    }

    // Get Deltas
    double deltaMinus = 0.0;
    double deltaPlus = 0.0;
    if(firstCell>-1){
      deltaMinus = 0.5*(topology->cellLengths[loopA][currentCellCoords[loopA]] + topology->cellLengths[loopA][currentCellCoords[loopA]-1]);
    }else{
      deltaMinus = 0.0;
    }
    if(secondCell>-1){
      deltaPlus = 0.5*(topology->cellLengths[loopA][currentCellCoords[loopA]] + topology->cellLengths[loopA][currentCellCoords[loopA]+1]);
    }else{
      deltaPlus = 0.0;
    }

    // Eval Cell Distance for Coordinate LoopA
    double firstVComponent = 0.0;
    double secondVComponent = 0.0;
    double currentVComponent = 0.0;
    // First Component
    if(firstCell>-1){
      firstVComponent = cells[firstCell].getQuantity(qtyID);
    }else{
      firstVComponent = 0.0;
    }
    // Second Component
    if(secondCell>-1){
      secondVComponent = cells[secondCell].getQuantity(qtyID);
    }else{
      secondVComponent = 0.0;
    }
    // Current Component
    currentVComponent = cells[currentCell].getQuantity(qtyID);

    // FIRST DERIVS
    if(firstCell<0){
      // Simple Euler Formula
      gradient[loopA] = (secondVComponent-currentVComponent)/(deltaPlus);
    }else if (secondCell<0){
      // Simple Euler Formula
      gradient[loopA] = (currentVComponent-firstVComponent)/(deltaMinus);
    }else if((firstCell>-1)&&(secondCell>-1)){
      // Central Difference Formula
      gradient[loopA] = (secondVComponent-firstVComponent)/(deltaPlus + deltaMinus);
    }else{
      // Show Error Message
      throw MRIException("Error: Both First and Second Cells are Zero in EvalFirstSpaceDerivs");
    }
  }
}

// Print The Derivatives
void printDerivatives(double** firstDerivs, double** secondDerivs){
  printf("First Derivatives\n");
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    printf("%e %e %e\n",firstDerivs[loopA][0],firstDerivs[loopA][1],firstDerivs[loopA][2]);
  }
  printf("\n");
  printf("Second Derivatives\n");
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    printf("%e %e %e\n",secondDerivs[loopA][0],secondDerivs[loopA][1],secondDerivs[loopA][2]);
  }
  
}

// EVAL PRESSURE GRADIENTS
void MRISequence::computePressureGradients(MRIThresholdCriteria* threshold){
  
  // Allocate Local Velocity Gradients
  MRIDoubleVec timeDerivs(kNumberOfDimensions);
  MRIDoubleVec tmp;
  // First and Second Derivatives
  MRIDoubleMat firstDerivs;
  MRIDoubleMat secondDerivs;
  MRIDoubleMat ReynoldsStressGrad;

  firstDerivs.resize(kNumberOfDimensions);
  secondDerivs.resize(kNumberOfDimensions);
  ReynoldsStressGrad.resize(kNumberOfDimensions);

  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA].resize(kNumberOfDimensions);
    secondDerivs[loopA].resize(kNumberOfDimensions);
    ReynoldsStressGrad[loopA].resize(kNumberOfDimensions);
  }
  MRIDoubleVec currentGradient(kNumberOfDimensions);
    
  // Write Message
  writeSchMessage(string("\n"));
  writeSchMessage(string("PRESS GRADIENT COMPUTATION ------------------------------------\n"));
  
  // Loop through the Sequences
  for(int loopA=0;loopA<totalScans;loopA++){
    
    // Write Message
    writeSchMessage(string("Computing Pressure Gradient: Step "+MRIUtils::intToStr(loopA+1)+"/"+MRIUtils::intToStr(totalScans)+"..."));

    // EVALUATE REYNOLDS STRESSES
    sequence[loopA]->evalReynoldsStress(threshold);

    //Loop Through the Cells
    for(int loopB=0;loopB<topology->totalCells;loopB++){
      
      // Eval Time Derivatives
      evalTimeDerivs(loopA/*Scan*/,loopB/*Cell*/,timeDerivs);

      // Eval First Derivatives in space
      sequence[loopA]->evalSpaceDerivs(loopB/*Cell*/,threshold,firstDerivs,secondDerivs);

      // EVAL REYNOLDS STRESS GRADIENTS
      if(sequence[loopA]->reynoldsStress.size() > 0){
        sequence[loopA]->evalReynoldsStressGradient(loopB/*Cell*/,ReynoldsStressGrad);
      }

      // Eval First Derivatives in space
      sequence[loopA]->evalCellPressureGradients(loopB,timeDerivs,firstDerivs,secondDerivs,ReynoldsStressGrad,currentGradient);
       
      // Store Gradients
      tmp.clear();
      tmp.push_back(currentGradient[0]);
      tmp.push_back(currentGradient[1]);
      tmp.push_back(currentGradient[2]);
      sequence[loopA]->qtyGradient.push_back(tmp);
    }
    writeSchMessage(string("Done.\n"));
  }
}

// CHECK IF NEIGHBOR ARE NOT VISITED
bool MRIScan::areThereNotVisitedNeighbor(int cell, bool* visitedCell){
  std::vector<int> otherCells;
  // Should I consider also the diagonal neighbors!!!
  getCartesianNeighbourCells(cell,otherCells,false);
  bool areThereVisited = false;
  for(size_t loopA=0;loopA<otherCells.size();loopA++){
    if((otherCells[loopA]>-1)&&(isInnerCell(otherCells[loopA]))){
      areThereVisited = ((areThereVisited)||(!visitedCell[otherCells[loopA]]));
    }
  }
  // Return
  return areThereVisited;
}

// CHECK IF NEIGHBOR ARE VISITED
bool MRIScan::areThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor){
  std::vector<int> otherCells;
  // Should I consider also the diagonal neighbors!!!
  getCartesianNeighbourCells(cell,otherCells,false);
  bool areThereVisited = false;
  for(size_t loopA=0;loopA<otherCells.size();loopA++){
    if((otherCells[loopA]>-1)&&(isInnerCell(otherCells[loopA]))){
      areThereVisited = ((areThereVisited)||((visitedCell[otherCells[loopA]])&&(!isBoundaryCell[loopA])));
      if ((visitedCell[otherCells[loopA]])&&(!isBoundaryCell[loopA])){
        visitedNeighbor = otherCells[loopA];
      }
    }
  }
  // Return
  return areThereVisited;
}

// GET CELL FROM STACK
int MRIScan::getCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished, bool& secondStage){
  // Check if Stack is Empty
  if (cellStack.size() ==  0){
    finished = true;
    return -1;
  }
  // Look for next Cell
  bool found = false;
  unsigned int count = 0;
  while ((!found)&&(count<cellStack.size())){
    // CAREFUL: CONSIDER NOISY GRADIENT POINTS!!!
    // FILIPPO JET
    //found = (visitedCell[cellStack[count]])&&
    //        (AreThereNotVisitedNeighbor(cellStack[count],visitedCell)&&
    //        (cellPoints[cellStack[count]].filteredVel[0]<0.2));
    // JET JULIEN
    //if(secondStage){
    //  printf("OK\n");
    found = (visitedCell[cellStack[count]])&&(areThereNotVisitedNeighbor(cellStack[count],visitedCell));
    //}else{
    //  found = (visitedCell[cellStack[count]])&&
    //          (AreThereNotVisitedNeighbor(cellStack[count],visitedCell)&&
    //          //(cellPoints[cellStack[count]].filteredVel[0]<0.001));
    //          (cellPoints[cellStack[count]].filteredVel[0]<1500.0));
    //          //(cellPoints[cellStack[count]].filteredVel[0]<0.35));
    //}
    // Update Count
    if (!found){
      count++;
    }
  }
  if (found){
    // Return
    int newCell = cellStack[count];
    // Erase From Stack
    cellStack.erase(cellStack.begin(),cellStack.begin()+count+1);
    // Return 
    return newCell;
  }else{
    finished = true;
    return -1;
  }
}

// Check if there are still unvisited Cells
int findFirstNotVisited(int cellTotal, bool* visitedCell, std::vector<int> cellStack){
  int newCell = -1;
  int count = 0;
  // Check if the stack size is full
  /*if(cellStack.size()>0){
    printf("Entrato!!! Stack Pieno\n");
    bool found = false;
    while ((!found)&&(count<cellStack.size())){
      //found = !visitedCell[count];
      //found = AreThereNotVisitedNeighbor(count,visitedCell);
      found = ((visitedCell[cellStack[count]])&&(AreThereNotVisitedNeighbor(cellStack[count],visitedCell)));
      // Update
      if (!found){
        count++;
      }else{
        newCell = count;
      }
    }
  }else{*/
    //printf("Entrato!!! Stack Vuoto\n");
    bool found = false;
    while ((!found)&&(count<cellTotal)){
      found = !visitedCell[count];
      //found = AreThereNotVisitedNeighbor(count,visitedCell);
      //found = ((visitedCell[count])&&(notVisitedList[count]));
      // Update
      if (!found){
        count++;
      }else{
        newCell = count;
      }
    }
  /*}*/
  return newCell;
}


