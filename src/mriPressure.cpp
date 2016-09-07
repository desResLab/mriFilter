#include <math.h>

#include "mriScan.h"
#include "mriStructuredScan.h"
#include "mriConstants.h"
#include "mriSequence.h"

#include "mriCellMaterial.h"
#include "schMessages.h"
#include "mriUtils.h"
#include "mriException.h"

// ===========================================
// COMPUTE PRESSURE GRADIENTS FOR GENERIC CELL
// ===========================================
void MRIScan::EvalCellPressureGradients(int currentCell, MRICellMaterial material,
                                        double* timeDeriv, double** firstDerivs, double** secondDerivs,
                                        double** ReynoldsStressGrad,
                                        double* pressureGrad){

  // SET PARAMETER
  int pressureGradientType = 3;

  // INIT
  double gravityTerm = 0.0;
  double viscousTerm = 0.0;
  double convectiveTerm = 0.0;
  double ReynoldsStressTerm = 0.0;
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    // Eval Viscous Term
    viscousTerm = material.viscosity * (secondDerivs[0][loopA] + secondDerivs[1][loopA] + secondDerivs[2][loopA]);
    // Eval Convective Term
    convectiveTerm = material.density * (timeDeriv[loopA] +
                                         cellPoints[currentCell].velocity[0] * firstDerivs[0][loopA]+
                                         cellPoints[currentCell].velocity[1] * firstDerivs[1][loopA]+
                                         cellPoints[currentCell].velocity[2] * firstDerivs[2][loopA]);
    // Eval Reynolds Stress Gradients
    if(hasReynoldsStress){
      ReynoldsStressTerm = material.density * (ReynoldsStressGrad[0][loopA] + ReynoldsStressGrad[1][loopA] + ReynoldsStressGrad[2][loopA]);
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
void MRISequence::EvalTimeDerivs(int currentScan, int currentCell,double* timeDeriv){
  // Eval Difference Formulae
  if (totalScans>1){
    // multiple Scan
    bool isFirstScan = (currentScan == 0);
    bool isLastScan = (currentScan == (totalScans-1));
    if (isFirstScan){
      if (isCyclic){
        double deltaT2 = sequence[currentScan+1]->scanTime - sequence[currentScan]->scanTime;
        // Assume Same Time Step
        double deltaT1 = deltaT2;
        // Get The Velocities
        double firstPart[3] = {0.0};
        double secondPart[3] = {0.0};        
        double* currVel = sequence[currentScan]->cellPoints[currentCell].velocity;
        double* prevVel = sequence[totalScans-1]->cellPoints[currentCell].velocity;
        double* nextVel = sequence[currentScan+1]->cellPoints[currentCell].velocity;      
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
        double* currVel = sequence[currentScan]->cellPoints[currentCell].velocity;
        double* nextVel = sequence[currentScan+1]->cellPoints[currentCell].velocity;
        double* nextNextVel = sequence[currentScan+2]->cellPoints[currentCell].velocity;
        // Eval Time Derivatives
        timeDeriv[0] = findParabolicDeriv(currVel[0],nextVel[0],nextNextVel[0],deltaT1,deltaT2,true);
        timeDeriv[1] = findParabolicDeriv(currVel[1],nextVel[1],nextNextVel[1],deltaT1,deltaT2,true);
        timeDeriv[2] = findParabolicDeriv(currVel[2],nextVel[2],nextNextVel[2],deltaT1,deltaT2,true);
      }else{
        // Euler Formula
        double deltaT = sequence[currentScan+1]->scanTime - sequence[currentScan]->scanTime;
        // Use Simple Difference Scheme
        timeDeriv[0] = (sequence[currentScan+1]->cellPoints[currentCell].velocity[0] - sequence[currentScan]->cellPoints[currentCell].velocity[0])/(deltaT);
        timeDeriv[1] = (sequence[currentScan+1]->cellPoints[currentCell].velocity[1] - sequence[currentScan]->cellPoints[currentCell].velocity[1])/(deltaT);
        timeDeriv[2] = (sequence[currentScan+1]->cellPoints[currentCell].velocity[2] - sequence[currentScan]->cellPoints[currentCell].velocity[2])/(deltaT);        
      }
    }else if(isLastScan){
      if (isCyclic){
        double deltaT1 = sequence[currentScan]->scanTime - sequence[currentScan-1]->scanTime;
        // Assume Same Time Step
        double deltaT2 = deltaT1;
        // Get The Velocities
        double firstPart[3] = {0.0};
        double secondPart[3] = {0.0};
        double* currVel = sequence[currentScan]->cellPoints[currentCell].velocity;
        double* prevVel = sequence[currentScan-1]->cellPoints[currentCell].velocity;
        double* nextVel = sequence[0]->cellPoints[currentCell].velocity;      
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
        double* currVel = sequence[currentScan-2]->cellPoints[currentCell].velocity;
        double* nextVel = sequence[currentScan-1]->cellPoints[currentCell].velocity;
        double* nextNextVel = sequence[currentScan]->cellPoints[currentCell].velocity;
        // Eval Time Derivatives
        timeDeriv[0] = findParabolicDeriv(currVel[0],nextVel[0],nextNextVel[0],deltaT1,deltaT2,false);
        timeDeriv[1] = findParabolicDeriv(currVel[1],nextVel[1],nextNextVel[1],deltaT1,deltaT2,false);
        timeDeriv[2] = findParabolicDeriv(currVel[2],nextVel[2],nextNextVel[2],deltaT1,deltaT2,false);        
      }else{
        // Euler Formula
        double deltaT = sequence[currentScan]->scanTime - sequence[currentScan-1]->scanTime;
        timeDeriv[0] = (sequence[currentScan]->cellPoints[currentCell].velocity[0] - sequence[currentScan-1]->cellPoints[currentCell].velocity[0])/(deltaT);
        timeDeriv[1] = (sequence[currentScan]->cellPoints[currentCell].velocity[1] - sequence[currentScan-1]->cellPoints[currentCell].velocity[1])/(deltaT);
        timeDeriv[2] = (sequence[currentScan]->cellPoints[currentCell].velocity[2] - sequence[currentScan-1]->cellPoints[currentCell].velocity[2])/(deltaT);                
      }
    }else{
      double firstPart[3] = {0.0};
      double secondPart[3] = {0.0};
      double deltaT1 = sequence[currentScan]->scanTime - sequence[currentScan-1]->scanTime;
      double deltaT2 = sequence[currentScan+1]->scanTime - sequence[currentScan]->scanTime;
      // Get The Velocities
      double* currVel = sequence[currentScan]->cellPoints[currentCell].velocity;
      double* prevVel = sequence[currentScan-1]->cellPoints[currentCell].velocity;
      double* nextVel = sequence[currentScan+1]->cellPoints[currentCell].velocity;      
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
void MRISequence::EvalScanTimeDerivs(int currentScan,MRIDoubleMat& timeDeriv){
  timeDeriv.clear();
  MRIDoubleVec temp;
  double cellTimeDeriv[3];
  for(int loopA=0;loopA<sequence[currentScan]->totalCellPoints;loopA++){
    temp.clear();
    EvalTimeDerivs(currentScan,loopA,cellTimeDeriv);
    temp.push_back(cellTimeDeriv[0]);
    temp.push_back(cellTimeDeriv[1]);
    temp.push_back(cellTimeDeriv[2]);
    timeDeriv.push_back(temp);
  }
}

// ================================
// EVAL REYNOLDS STRESS DERIVATIVES
// ================================
void MRISequence::EvalScanReynoldsStressDerivs(int currentScan,MRIDoubleMat& reynoldsDeriv){
  reynoldsDeriv.clear();
  MRIDoubleVec temp;
  double** rsg;
  rsg = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    rsg[loopA] = new double[3];
  }
  for(int loopA=0;loopA<sequence[currentScan]->totalCellPoints;loopA++){
    temp.clear();
    sequence[currentScan]->EvalReynoldsStressGradient(loopA, rsg);
    temp.push_back(rsg[0][0] + rsg[1][0] + rsg[2][0]);
    temp.push_back(rsg[0][1] + rsg[1][1] + rsg[2][1]);
    temp.push_back(rsg[0][2] + rsg[1][2] + rsg[2][2]);
    reynoldsDeriv.push_back(temp);
  }
  // Free Matrix
  for(int loopA=0;loopA<3;loopA++){
    delete [] rsg[loopA];
  }
  delete [] rsg;
}

// ==========================================
// EVAL FIRST AND SECOND DERIVATIVES IN SPACE
// ==========================================
void MRIStructuredScan::EvalSpaceDerivs(int currentCell, MRIThresholdCriteria* threshold, double** firstDerivs, double** secondDerivs){
  // FirstDerivs
  // DVX/DX DVY/DX DVZ/DX
  // DVX/DY DVY/DY DVZ/DY
  // DVX/DZ DVY/DZ DVZ/DZ

  // Get quantity for threshold evaluation
  double cellQty = 0.0;

  // Map Index To Coords
  int* currentCellCoords = new int[kNumberOfDimensions];
  MapIndexToCoords(currentCell,currentCellCoords);
  int firstCell,secondCell,thirdCell,fourthCell,nextCell;

  // Assemble Terms
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    switch(loopA){
      case 0:
        // X Deriv
        if((currentCellCoords[0]-1)<0){
          firstCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0]-1,currentCellCoords[1],currentCellCoords[2]);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
            firstCell = -1;
          }else{
            firstCell = nextCell;
          }
        }
        if((currentCellCoords[0]+1)>(cellTotals[0]-1)){
          secondCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0]+1,currentCellCoords[1],currentCellCoords[2]);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
            secondCell = -1;
          }else{
            secondCell = nextCell;
          }
        }
        if((currentCellCoords[0]-2)<0){
          thirdCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0]-2,currentCellCoords[1],currentCellCoords[2]);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
            thirdCell = -1;
          }else{
            thirdCell = nextCell;
          }
        }
        if((currentCellCoords[0]+2)>(cellTotals[0]-1)){
          fourthCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0]+2,currentCellCoords[1],currentCellCoords[2]);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
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
          nextCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]-1,currentCellCoords[2]);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
            firstCell = -1;
          }else{
            firstCell = nextCell;
          }
        }
        if((currentCellCoords[1]+1)>(cellTotals[1]-1)){
          secondCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]+1,currentCellCoords[2]);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
            secondCell = -1;
          }else{
            secondCell = nextCell;
          }
        }
        if((currentCellCoords[1]-2)<0){
          thirdCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]-2,currentCellCoords[2]);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
            thirdCell = -1;
          }else{
            thirdCell = nextCell;
          }
        }
        if((currentCellCoords[1]+2)>(cellTotals[1]-1)){
          fourthCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]+2,currentCellCoords[2]);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
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
          nextCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]-1);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
            firstCell = -1;
          }else{
            firstCell = nextCell;
          }
        }
        if((currentCellCoords[2]+1)>(cellTotals[2]-1)){
          secondCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]+1);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
            secondCell = -1;
          }else{
            secondCell = nextCell;
          }
        }
        if((currentCellCoords[2]-2)<0){
          thirdCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]-2);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
            thirdCell = -1;
          }else{
            thirdCell = nextCell;
          }
        }
        if((currentCellCoords[2]+2)>(cellTotals[2]-1)){
          fourthCell = -1;
        }else{
          nextCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]+2);
          cellQty = cellPoints[nextCell].getQuantity(threshold->thresholdQty);
          if(threshold->MeetsCriteria(cellQty)){
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
      deltaMinus = 0.5*(cellLengths[loopA][currentCellCoords[loopA]] + cellLengths[loopA][currentCellCoords[loopA]-1]);
    }else{
      deltaMinus = 0.0;
    }
    if(secondCell>-1){
      deltaPlus = 0.5*(cellLengths[loopA][currentCellCoords[loopA]] + cellLengths[loopA][currentCellCoords[loopA]+1]);
    }else{
      deltaPlus = 0.0;
    }
    if(thirdCell>-1){
      deltaMinusMinus = 0.5*(cellLengths[loopA][currentCellCoords[loopA]-1] + cellLengths[loopA][currentCellCoords[loopA]-2]);
    }else{
      deltaMinusMinus = 0.0;
    }
    if(fourthCell>-1){
      deltaPlusPlus = 0.5*(cellLengths[loopA][currentCellCoords[loopA]+1] + cellLengths[loopA][currentCellCoords[loopA]+2]);
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
        firstVComponent = cellPoints[firstCell].velocity[loopB];
      }else{
        firstVComponent = 0.0;
      }
      // Second Component
      if(secondCell>-1){
        secondVComponent = cellPoints[secondCell].velocity[loopB];
      }else{
        secondVComponent = 0.0;
      }
      // Third Component
      if(thirdCell>-1){
        thirdVComponent = cellPoints[thirdCell].velocity[loopB];
      }else{
        thirdVComponent = 0.0;
      }
      // Fourth Component
      if(fourthCell>-1){
        fourthVComponent = cellPoints[fourthCell].velocity[loopB];
      }else{
        fourthVComponent = 0.0;
      }

      // Check if the current Cell has a
      cellQty = cellPoints[currentCell].getQuantity(threshold->thresholdQty);
      if(threshold->MeetsCriteria(cellQty)){
        firstDerivs[loopA][loopB] = 0.0;
        secondDerivs[loopA][loopB] = 0.0;
      }else{
        // Current Component
        currentVComponent = cellPoints[currentCell].velocity[loopB];

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
          throw MRIPressureComputationException("Error: Both First and Second Cells are Zero in EvalFirstSpaceDerivs");
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
          throw MRIPressureComputationException("Error: Both First and Second Cells are Zero in EvalFirstSpaceDerivs");
        }
      }
    }
  }
  // Deallocate
  delete [] currentCellCoords;
}

// ==========================================
// EVAL FIRST AND SECOND DERIVATIVES IN SPACE
// ==========================================
void MRIStructuredScan::EvalSpaceGradient(int currentCell,int qtyID, double* gradient){
  // Map Index To Coords
  int* currentCellCoords = new int[kNumberOfDimensions];
  MapIndexToCoords(currentCell,currentCellCoords);
  int firstCell,secondCell;
  // Assemble Terms
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    switch(loopA){
      case 0:
        // X Deriv
        if((currentCellCoords[0]-1)<0){
          firstCell = -1;
        }else{
          firstCell = MapCoordsToIndex(currentCellCoords[0]-1,currentCellCoords[1],currentCellCoords[2]);
        }
        if((currentCellCoords[0]+1)>(cellTotals[0]-1)){
          secondCell = -1;
        }else{
          secondCell = MapCoordsToIndex(currentCellCoords[0]+1,currentCellCoords[1],currentCellCoords[2]);
        }
        break;
      case 1:
        // Y Deriv
        if((currentCellCoords[1]-1)<0){
          firstCell = -1;
        }else{
          firstCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]-1,currentCellCoords[2]);
        }
        if((currentCellCoords[1]+1)>(cellTotals[1]-1)){
          secondCell = -1;
        }else{
          secondCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]+1,currentCellCoords[2]);
        }
        break;
      case 2:
        // Z Deriv
        if((currentCellCoords[2]-1)<0){
          firstCell = -1;
        }else{
          firstCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]-1);
        }
        if((currentCellCoords[2]+1)>(cellTotals[2]-1)){
          secondCell = -1;
        }else{
          secondCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]+1);
        }
        break;
    }

    // Get Deltas
    double deltaMinus = 0.0;
    double deltaPlus = 0.0;
    if(firstCell>-1){
      deltaMinus = 0.5*(cellLengths[loopA][currentCellCoords[loopA]] + cellLengths[loopA][currentCellCoords[loopA]-1]);
    }else{
      deltaMinus = 0.0;
    }
    if(secondCell>-1){
      deltaPlus = 0.5*(cellLengths[loopA][currentCellCoords[loopA]] + cellLengths[loopA][currentCellCoords[loopA]+1]);
    }else{
      deltaPlus = 0.0;
    }

    // Eval Cell Distance for Coordinate LoopA
    double firstVComponent = 0.0;
    double secondVComponent = 0.0;
    double currentVComponent = 0.0;
    // First Component
    if(firstCell>-1){
      firstVComponent = cellPoints[firstCell].getQuantity(qtyID);
    }else{
      firstVComponent = 0.0;
    }
    // Second Component
    if(secondCell>-1){
      secondVComponent = cellPoints[secondCell].getQuantity(qtyID);
    }else{
      secondVComponent = 0.0;
    }
    // Current Component
    currentVComponent = cellPoints[currentCell].getQuantity(qtyID);

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
      throw MRIPressureComputationException("Error: Both First and Second Cells are Zero in EvalFirstSpaceDerivs");
    }
  }
  // Deallocate
  delete [] currentCellCoords;
}


// GET THE INDEX OF THE REYNOLDS STRESS COMPONENT
int getReynoldsStressIndex(int loopA,int loopB){
  // FirstDerivs
  // DRXX/DX DRYX/DX DRZX/DX
  // DRXY/DY DRYY/DY DRZY/DY
  // DRXZ/DZ DRYZ/DZ DRZZ/DZ
  // REYNOLDS STRESS INDEXES
  // 0-RXX 1-RXY 2-RXZ 3-RYY 4-RYZ 5-RZZ
  int result = 0;
  switch(loopA){
    case 0:
      switch(loopB){
        case 0:
          result = 0;
          break;
        case 1:
          result = 1;
          break;
        case 2:
          result = 2;
          break;
      }
      break;
    case 1:
      switch(loopB){
        case 0:
          result = 1;
          break;
        case 1:
          result = 3;
          break;
        case 2:
          result = 4;
          break;
      }
      break;
    case 2:
      switch(loopB){
        case 0:
          result = 2;
          break;
        case 1:
          result = 4;
          break;
        case 2:
          result = 5;
          break;
      }
      break;
  }
  return result;
}

// ==============================
// EVAL REYNOLDS STRESS GRADIENTS
// ==============================
void MRIStructuredScan::EvalReynoldsStressGradient(int currentCell, double** ReynoldsStressGradient){
  // FirstDerivs
  // DRXX/DX DRYX/DX DRZX/DX
  // DRXY/DY DRYY/DY DRZY/DY
  // DRXZ/DZ DRYZ/DZ DRZZ/DZ
  // Map Index To Coords
  int* currentCellCoords = new int[kNumberOfDimensions];
  MapIndexToCoords(currentCell,currentCellCoords);
  int firstCell,secondCell;
  // Assemble Terms
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    switch(loopA){
      case 0:
        // X Deriv
        if((currentCellCoords[0]-1)<0){
          firstCell = -1;
        }else{
          firstCell = MapCoordsToIndex(currentCellCoords[0]-1,currentCellCoords[1],currentCellCoords[2]);
        }
        if((currentCellCoords[0]+1)>(cellTotals[0]-1)){
          secondCell = -1;
        }else{
          secondCell = MapCoordsToIndex(currentCellCoords[0]+1,currentCellCoords[1],currentCellCoords[2]);
        }
        break;
      case 1:
        // Y Deriv
        if((currentCellCoords[1]-1)<0){
          firstCell = -1;
        }else{
          firstCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]-1,currentCellCoords[2]);
        }
        if((currentCellCoords[1]+1)>(cellTotals[1]-1)){
          secondCell = -1;
        }else{
          secondCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1]+1,currentCellCoords[2]);
        }
        break;
      case 2:
        // Z Deriv
        if((currentCellCoords[2]-1)<0){
          firstCell = -1;
        }else{
          firstCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]-1);
        }
        if((currentCellCoords[2]+1)>(cellTotals[2]-1)){
          secondCell = -1;
        }else{
          secondCell = MapCoordsToIndex(currentCellCoords[0],currentCellCoords[1],currentCellCoords[2]+1);
        }
        break;
    }

    // Get Deltas
    double deltaMinus = 0.0;
    double deltaPlus = 0.0;
    if(firstCell>-1){
      deltaMinus = 0.5*(cellLengths[loopA][currentCellCoords[loopA]] + cellLengths[loopA][currentCellCoords[loopA]-1]);
    }else{
      deltaMinus = 0.0;
    }
    if(secondCell>-1){
      deltaPlus = 0.5*(cellLengths[loopA][currentCellCoords[loopA]] + cellLengths[loopA][currentCellCoords[loopA]+1]);
    }else{
      deltaPlus = 0.0;
    }

    // Eval Cell Distance for Coordinate LoopA
    double firstVComponent = 0.0;
    double secondVComponent = 0.0;
    double currentVComponent = 0.0;
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      int ReStressIndex = getReynoldsStressIndex(loopA,loopB);
      // First Component
      if(firstCell>-1){
        firstVComponent = cellPoints[firstCell].ReStress[ReStressIndex];
      }else{
        firstVComponent = 0.0;
      }
      // Second Component
      if(secondCell>-1){
        secondVComponent = cellPoints[secondCell].ReStress[ReStressIndex];
      }else{
        secondVComponent = 0.0;
      }
      // Current Component
      currentVComponent = cellPoints[currentCell].ReStress[ReStressIndex];

      // FIRST DERIVS
      if(firstCell<0){
        // Simple Euler Formula
        ReynoldsStressGradient[loopA][loopB] = (secondVComponent-currentVComponent)/(deltaPlus);
      }else if (secondCell<0){
        // Simple Euler Formula
        ReynoldsStressGradient[loopA][loopB] = (currentVComponent-firstVComponent)/(deltaMinus);
      }else if((firstCell>-1)&&(secondCell>-1)){
        // Central Difference Formula
        ReynoldsStressGradient[loopA][loopB] = (secondVComponent-firstVComponent)/(deltaPlus + deltaMinus);
      }else{
        // Show Error Message
        throw MRIPressureComputationException("Error: Both First and Second Cells are Zero in EvalFirstSpaceDerivs");
      }

    }
  }
  // Deallocate
  delete [] currentCellCoords;
}


// Print The Derivatives
void PrintDerivatives(double** firstDerivs, double** secondDerivs){
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
void MRISequence::ComputePressureGradients(MRIThresholdCriteria* threshold){
  
  // Allocate Local Velocity Gradients
  double* timeDerivs = new double[kNumberOfDimensions];
  // First and Second Derivatives
  double** firstDerivs = new double*[kNumberOfDimensions];
  double** secondDerivs = new double*[kNumberOfDimensions];
  double** ReynoldsStressGrad = new double*[kNumberOfDimensions];
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA] = new double[kNumberOfDimensions];
    secondDerivs[loopA] = new double[kNumberOfDimensions];
    ReynoldsStressGrad[loopA] = new double[kNumberOfDimensions];
  }
  double* currentGradient = new double[kNumberOfDimensions];
  
  // FILIPPO JET
  MRICellMaterial* material = new MRICellMaterial(1.002e-3,998.2);
  // Create an Unique Material: USE BLOOD
  //MRICellMaterial* material = new MRICellMaterial(3.5e-3,1060.0);
  // JULIEN JET
  //MRICellMaterial* material = new MRICellMaterial(0.021,0.78733);
  
  // Write Message
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("PRESS GRADIENT COMPUTATION ------------------------------------\n"));
  
  // Loop through the Sequences
  for(int loopA=0;loopA<totalScans;loopA++){
    
    // Write Message
    WriteSchMessage(std::string("Computing Pressure Gradient: Step "+MRIUtils::IntToStr(loopA+1)+"/"+MRIUtils::IntToStr(totalScans)+"..."));

    // EVALUATE REYNOLDS STRESSES
    sequence[loopA]->EvalReynoldsStressComponent(threshold);

    //Loop Through the Cells}
    for(int loopB=0;loopB<sequence[loopA]->totalCellPoints;loopB++){
      
      // Eval Time Derivatives
      EvalTimeDerivs(loopA/*Scan*/,loopB/*Cell*/,timeDerivs);

      // Eval First Derivatives in space
      sequence[loopA]->EvalSpaceDerivs(loopB/*Cell*/,threshold,firstDerivs,secondDerivs);

      // EVAL REYNOLDS STRESS GRADIENTS
      if(sequence[loopA]->hasReynoldsStress){
        sequence[loopA]->EvalReynoldsStressGradient(loopB/*Cell*/,ReynoldsStressGrad);
      }

      //if (loopB==100){
      //  PrintDerivatives(firstDerivs,secondDerivs);
      //}

      // Eval First Derivatives in space
      sequence[loopA]->EvalCellPressureGradients(loopB,*material,timeDerivs,firstDerivs,secondDerivs,ReynoldsStressGrad,currentGradient);
       
       // Store Gradients
       for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
         sequence[loopA]->cellPoints[loopB].pressGrad[loopC] = currentGradient[loopC];
       }
    }
    // Set Pressure Gradient Properties
    sequence[loopA]->hasPressureGradient = true;
    WriteSchMessage(std::string("Done.\n"));
  }
  // Deallocate
  delete [] timeDerivs;
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    delete [] firstDerivs[loopA];
    delete [] secondDerivs[loopA];
  }  
  delete [] firstDerivs;
  delete [] secondDerivs;
  delete [] currentGradient;
}

// ITERATIVE EVALUATION OF PRESSURE
void MRIScan::EvalPressureIterative(int currentCell, double currentValue, bool* visitedCell, std::vector<int> &cellStack,int& cellCount){
  std::vector<int> otherCells;
  visitedCell[currentCell] = true;
  // Eval Neighbours
  bool useCount = false;
  int cell = 0;
  // Should I visit also diagonal neighbor!!!
  GetCartesianNeighbourCells(currentCell,otherCells,true);
  double diff[3] = {0.0};
  double avGradient[3] = {0.0};
  // Loop through Neighbours
  for(size_t loopA=0;loopA<otherCells.size();loopA++){
    cell = otherCells[loopA];
    if((cell>-1)&&(cell<totalCellPoints)&&(!visitedCell[cell])&&(IsInnerCell(cell))){
      visitedCell[cell] = true;
      // Get Position
      diff[0] = cellPoints[cell].position[0] - cellPoints[currentCell].position[0];
      diff[1] = cellPoints[cell].position[1] - cellPoints[currentCell].position[1];
      diff[2] = cellPoints[cell].position[2] - cellPoints[currentCell].position[2];
      // Eval Average Gradient
      for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
        avGradient[loopB] = 0.5*(cellPoints[currentCell].pressGrad[loopB] + cellPoints[cell].pressGrad[loopB]);
        //avGradient[loopB] = cellPoints[currentCell].pressGrad[loopB];
        //avGradient[loopB] = cellPoints[cell].pressGrad[loopB];
      }
      // Eval Pressure
      if(useCount){
        cellCount++;
        cellPoints[currentCell].relPressure = cellCount;
      }else{
        cellPoints[cell].relPressure = currentValue;
        for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
          cellPoints[cell].relPressure = cellPoints[cell].relPressure + avGradient[loopB] * diff[loopB];
        }
      }
      // Add it to the Stack
      cellStack.push_back(cell);
    }
  }
}

// CHECK IF NEIGHBOR ARE NOT VISITED
bool MRIScan::AreThereNotVisitedNeighbor(int cell, bool* visitedCell){
  std::vector<int> otherCells;
  // Should I consider also the diagonal neighbors!!!
  GetCartesianNeighbourCells(cell,otherCells,false);
  bool areThereVisited = false;
  for(size_t loopA=0;loopA<otherCells.size();loopA++){
    if((otherCells[loopA]>-1)&&(IsInnerCell(otherCells[loopA]))){
      areThereVisited = ((areThereVisited)||(!visitedCell[otherCells[loopA]]));
    }
  }
  // Return
  return areThereVisited;
}

// CHECK IF NEIGHBOR ARE VISITED
bool MRIScan::AreThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor){
  std::vector<int> otherCells;
  // Should I consider also the diagonal neighbors!!!
  GetCartesianNeighbourCells(cell,otherCells,false);
  bool areThereVisited = false;
  for(size_t loopA=0;loopA<otherCells.size();loopA++){
    if((otherCells[loopA]>-1)&&(IsInnerCell(otherCells[loopA]))){
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
int MRIScan::GetCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished, bool& secondStage){
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
    found = (visitedCell[cellStack[count]])&&(AreThereNotVisitedNeighbor(cellStack[count],visitedCell));
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


// EVAL RELATIVE PRESSURES USING FLOOD-FILL
void MRIScan::EvalRelativePressure(int startingCell, double refPressure){
  double* currVel = NULL;
  double currModulus = 0.0;
  int toScanCellsCount = 0;
  // Allocate 
  bool* visitedCell = new bool[totalCellPoints];
  bool* isBoundaryCell = new bool[totalCellPoints];
  // Stack
  std::vector<int> cellStack;
  
  // Initialize Visited Cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Eval Velocity Magnitude
    currVel = cellPoints[loopA].velocity;
    currModulus = sqrt(currVel[0]*currVel[0] + currVel[1]*currVel[1] + currVel[2]*currVel[2]);        
    // If Inner than OK
    // FILIPPO JET
    //if(IsInnerCell(loopA)&&(currModulus>kMathZero)&&(cellPoints[loopA].concentration>1.0e-3)){
    // LUNG
    if(IsInnerCell(loopA)&&(currModulus>kMathZero)&&(cellPoints[loopA].concentration<0.2)){
    // JET JULIEN
    //if(IsInnerCell(loopA)&&(currModulus>kMathZero)&&(cellPoints[loopA].filteredVel[0]<3500.0)){
    //if(IsInnerCell(loopA)&&(currModulus>kMathZero)){
      toScanCellsCount++;
      visitedCell[loopA] = false;
      isBoundaryCell[loopA] = false;
    }else{
      visitedCell[loopA] = true;
      isBoundaryCell[loopA] = true;
    }
  }
  
  // Start Iterations
  bool finished = false;
  bool secondStage = false;
  int nextCell = 0;
  int cellCount = 0;
  int cellNumber = 0;
  int bookmark = 0;
  int progress = 0;
  int remainingCell = -1;
  int percentCounted = 0;
  printf("Started.\n");
  while(!finished){
    // Update cell count
    cellCount++;
    // Check Progress
    progress = MRIUtils::round((((double)cellCount/(double)toScanCellsCount)*100));
    if (((progress % 10) == 0)&&((progress / 10) != percentCounted)){
      percentCounted = (progress / 10);
      printf("Flood Fill Status: %d%%\n",progress);
    }    
    // Eval Pressures Recursively
    EvalPressureIterative(startingCell,refPressure,visitedCell,cellStack,cellNumber);
    
    // Get Next Starting Cell
    nextCell = GetCellFromStack(cellStack,visitedCell,isBoundaryCell,finished,secondStage);
    
    // Check if all the cells have been visited
    if(finished){
      //secondStage = true;
      remainingCell = findFirstNotVisited(totalCellPoints,visitedCell,cellStack);
      //remainingCell = -1;
      if (remainingCell != -1){
        finished = false;
        nextCell = remainingCell;
        cellPoints[nextCell].relPressure = 0.0;
      }
    }
    // Update
    startingCell = nextCell;
    refPressure = cellPoints[nextCell].relPressure;
  }
  
  // Deallocate
  delete [] visitedCell;
  delete [] isBoundaryCell;
  // Set Relative Pressure Availability
  hasRelativePressure = true;
}

// PERFORM ITERATIVE AVERAGE FOR RELATIVE PRESSURES
void MRIScan::PerformPressureIterations(){
  // Pressure Iterations
  bool converged = false;
  int itCount = 0;
  int avCount = 0;
  std::vector<int> neighbours;
  double maxPressureDiff = 0.0;
  double avPressure = 0.0;
  double newPressure = 0.0;
  double currentPressure = 0.0;
  double neighborPressure = 0.0;
  // Vectors
  double grad[3] = {0.0};
  double diff[3] = {0.0};
  // Message
  WriteSchMessage(std::string("Applying Pressure Smoothing...\n"));
  while(!converged){
    itCount++;
    // Loop On All Cells
    maxPressureDiff = 0.0;
    for(int loopA=0;loopA<totalCellPoints;loopA++){
      if(IsInnerCell(loopA)){
        // Eval Current Pressure
        currentPressure = cellPoints[loopA].relPressure;
        // Eval The Avergage Of Neighbor Pressures
        // Should I consider also diagonal neighbor!!!
        GetCartesianNeighbourCells(loopA,neighbours,false);
        avPressure = 0.0;
        avCount = 0;
        for(int loopB=0;loopB<k3DNeighbors;loopB++){
          if((neighbours[loopB]>-1)&&(neighbours[loopB]<=totalCellPoints)&&(IsInnerCell(neighbours[loopB]))){
            avCount++;
            // Eval Gradient at Neighbor
            for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
              grad[loopC] = cellPoints[loopA].pressGrad[loopC];
            }
            // Eval Distance
            for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
              diff[loopC] = cellPoints[neighbours[loopB]].position[loopC]-cellPoints[loopA].position[loopC];
            }
            // Eval Pressure From Neighbor
            neighborPressure = cellPoints[neighbours[loopB]].relPressure;
            for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
              neighborPressure += diff[loopC]*grad[loopC];
            }
            avPressure += neighborPressure;
          }
        }
        // Eval Average
        avPressure = (avPressure/avCount);
        newPressure = 0.5 * currentPressure + 0.5 * avPressure;
        if(fabs(currentPressure)>kMathZero){
          if (fabs(newPressure-currentPressure)/fabs(currentPressure)>maxPressureDiff){
            maxPressureDiff = fabs(newPressure-currentPressure)/fabs(currentPressure);
          }
        }
        // Update Pressure
        cellPoints[loopA].relPressure = newPressure;
      }
    }
    // Show Message
    //MessageDlg('Current Pressure Norm: '+FloatToStr(MaxPressureDiff),mtInformation,[mbOK],0);
    // Check Convergence
    //converged = (maxPressureDiff < 5.0);
    WriteSchMessage(std::string("Max Pressure Diff: "+MRIUtils::FloatToStr(maxPressureDiff)+"\n"));
    converged = (itCount > 10);
  }
}

