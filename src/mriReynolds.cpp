# include "mriScan.h"

// ========================
// BOUSSINESQ APPROXIMATION
// ========================
// EVAL TURBULENT VISCOSITY
void mriScan::evalEddyViscosity_simple(mriDoubleVec& nuT){
  // Loop through all cells
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // JET Filippo
    nuT[loopA] = 0.01*0.5*0.0058;
    //JET Julien
    //nuT[loopA] = 0.01*86.1*1.0;
  }
}
// EVAL TURBULENT KINETIC ENERGY
void mriScan::evalTurbulentKineticEnergy(mriThresholdCriteria* threshold, mriDoubleVec& turbK){
  // First and Second Cell Derivatives

  mriDoubleMat firstDerivs;
  mriDoubleMat secondDerivs;
  firstDerivs.resize(3);
  secondDerivs.resize(3);
  for(int loopA=0;loopA<3;loopA++){
    firstDerivs[loopA].resize(3);
    secondDerivs[loopA].resize(3);
  }
  // Loop through all cells
  double tensorProduct = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Eval First and Derivative Tensor
    evalSpaceDerivs(loopA,threshold,firstDerivs,secondDerivs);
    // Eval Tensor Product
    tensorProduct = 0.0;
    for(int loopB=0;loopB<3;loopB++){
      for(int loopC=0;loopC<3;loopC++){
        tensorProduct += firstDerivs[loopB][loopC]*firstDerivs[loopB][loopC];
      }
    }
    // Careful SI Units
    turbK[loopA] = 2.0*tensorProduct*0.001;
  }
}

// ==============================================
// GET THE INDEX OF THE REYNOLDS STRESS COMPONENT
// ==============================================
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

// ========================
// EVAL TURBULENT VISCOSITY
// ========================
void mriScan::evalReynoldsStress(mriThresholdCriteria* threshold){
  // Allocate
  mriDoubleVec turbNu(topology->totalCells,0.0);
  mriDoubleVec turbK(topology->totalCells,0.0);
  // First and Second Cell Derivatives
  mriDoubleMat firstDerivs;
  mriDoubleMat secondDerivs;

  firstDerivs.resize(3);
  secondDerivs.resize(3);

  for(int loopA=0;loopA<3;loopA++){
    firstDerivs[loopA].resize(3);
    secondDerivs[loopA].resize(3);
    for(int loopB=0;loopB<3;loopB++){
      firstDerivs[loopA][loopB] = 0.0;
      secondDerivs[loopA][loopB] = 0.0;
    }
  }

  // CHOOSE WHICH REAYNOLDS TO USE
  int ReynoldsCriterion = 0;
  //int ReynoldsCriterion = 1;

  // WRITE MESSAGE
  writeSchMessage("\n");
  writeSchMessage("Computing Reynolds Stresses...");

  // Eval Turbulent Viscosity in current Scan
  evalEddyViscosity_simple(turbNu);

  // Eval Turbulent Kinetic Energy in current Scan
  evalTurbulentKineticEnergy(threshold,turbK);

  // Allocate Reynolds Stress matrix
  reynoldsStress.resize(topology->totalCells);
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    reynoldsStress[loopA].resize(6);
  }

  // Cycle through all cells
  for(int loopA=0;loopA<topology->totalCells;loopA++){

    // Eval First Derivative Tensor
    evalSpaceDerivs(loopA,threshold,firstDerivs,secondDerivs);

    if(ReynoldsCriterion == 0){
      // COMPLETE STRESSES
      // RXXs
      reynoldsStress[loopA][0] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][0]+firstDerivs[0][0]))-((2.0)/(3.0))*turbK[loopA];
      //cellPoints[loopA].ReStress[0] = turbK[loopA];
      // RXY
      reynoldsStress[loopA][1] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][1]+firstDerivs[1][0]));
      // RXZ
      reynoldsStress[loopA][2] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][2]+firstDerivs[2][0]));
      // RYY
      reynoldsStress[loopA][3] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[1][1]+firstDerivs[1][1]))-((2.0)/(3.0))*turbK[loopA];
      // RYZ
      reynoldsStress[loopA][4] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[1][2]+firstDerivs[2][1]));
      // RZZ
      reynoldsStress[loopA][5] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[2][2]+firstDerivs[2][2]))-((2.0)/(3.0))*turbK[loopA];
    }else if(ReynoldsCriterion == 1){
      // ONLY VELOCITY GRADIENTS
      // RXXs
      reynoldsStress[loopA][0] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][0]+firstDerivs[0][0]));
      // RXY
      reynoldsStress[loopA][1] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][1]+firstDerivs[1][0]));
      // RXZ
      reynoldsStress[loopA][2] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][2]+firstDerivs[2][0]));
      // RYY
      reynoldsStress[loopA][3] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[1][1]+firstDerivs[1][1]));
      // RYZ
      reynoldsStress[loopA][4] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[1][2]+firstDerivs[2][1]));
      // RZZ
      reynoldsStress[loopA][5] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[2][2]+firstDerivs[2][2]));
    }else if(ReynoldsCriterion == 2){
      // ONLY HYDROSTATIC PART
      // RXXs
      reynoldsStress[loopA][0] = -((2.0)/(3.0))*turbK[loopA];
      // RXY
      reynoldsStress[loopA][1] = 0.0;
      // RXZ
      reynoldsStress[loopA][2] = 0.0;
      // RYY
      reynoldsStress[loopA][3] = -((2.0)/(3.0))*turbK[loopA];
      // RYZ
      reynoldsStress[loopA][4] = 0.0;
      // RZZ
      reynoldsStress[loopA][5] = -((2.0)/(3.0))*turbK[loopA];
    }
  }

  // DONE
  writeSchMessage("Done.\n");
}

// ================================================
// EVAL SMAGORINSKY LILLY TURBULENT VISCOSITY MODEL
// ================================================
void mriScan::evalSmagorinskyLillyTurbViscosity(double density, double smagorinskyCoeff, mriThresholdCriteria* threshold, mriDoubleMat& turbViscosity){
  double modS = 0.0;
  double currCharDist = 0.0;
  double sTerm = 0.0;
  mriDoubleVec tmp;
  // First and Second Derivatives
  mriDoubleMat firstDerivs;
  mriDoubleMat secondDerivs;
  firstDerivs.resize(kNumberOfDimensions);
  secondDerivs.resize(kNumberOfDimensions);
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA].resize(kNumberOfDimensions);
    secondDerivs[loopA].resize(kNumberOfDimensions);
  }

  turbViscosity.clear();
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Eva Spatial Derivatives
    evalSpaceDerivs(loopA, threshold, firstDerivs, secondDerivs);
    // Evaluate the Module of the strain rate tensor
    sTerm = 0.0;
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
        // Pradtl
        sTerm += (0.5 * (firstDerivs[loopB][loopC] + firstDerivs[loopC][loopB])) * (0.5 * (firstDerivs[loopB][loopC] + firstDerivs[loopC][loopB]));
        // Baldwin and Lomax
        //sTerm += (0.5 * (firstDerivs[loopB][loopC] - firstDerivs[loopC][loopB])) * (0.5 * (firstDerivs[loopB][loopC] - firstDerivs[loopC][loopB]));
      }
    }
    currCharDist = pow(evalCellVolume(loopA),1.0/3.0);
    modS = sqrt(2.0 * sTerm);
    tmp.clear();
    tmp.push_back(density * (smagorinskyCoeff * currCharDist) * (smagorinskyCoeff * currCharDist) * modS);
    turbViscosity.push_back(tmp);
  }
}

// ==============================
// EVAL REYNOLDS STRESS GRADIENTS
// ==============================
void mriScan::evalReynoldsStressGradient(int currentCell, mriDoubleMat& ReynoldsStressGradient){
  // FirstDerivs
  // DRXX/DX DRYX/DX DRZX/DX
  // DRXY/DY DRYY/DY DRZY/DY
  // DRXZ/DZ DRYZ/DZ DRZZ/DZ
  // Map Index To Coords
  mriIntVec currentCellCoords(kNumberOfDimensions);
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
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      int ReStressIndex = getReynoldsStressIndex(loopA,loopB);
      // First Component
      if(firstCell>-1){
        firstVComponent = reynoldsStress[firstCell][ReStressIndex];
      }else{
        firstVComponent = 0.0;
      }
      // Second Component
      if(secondCell>-1){
        secondVComponent = reynoldsStress[secondCell][ReStressIndex];
      }else{
        secondVComponent = 0.0;
      }
      // Current Component
      currentVComponent = reynoldsStress[currentCell][ReStressIndex];

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
        throw mriException("ERROR: Both First and Second Cells are Zero in EvalFirstSpaceDerivs");
      }

    }
  }
}
