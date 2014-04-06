#include "mriUnstructuredScan.h"
#include "mriTurbModel.h"
#include "schMessages.h"

// ========================
// EVAL TURBULENT VISCOSITY
// ========================
void MRIScan::EvalReynoldsStressComponent(){
  // Allocate
  double* turbNu = new double[totalCellPoints];
  double* turbK = new double[totalCellPoints];
  // First and Second Cell Derivatives
  double** firstDerivs = new double*[3];
  double** secondDerivs = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    firstDerivs[loopA] = new double[3];
    secondDerivs[loopA] = new double[3];
  }

  // CHOOSE WHICH REAYNOLDS TO USE
  int ReynoldsCriterion = 0;
  //int ReynoldsCriterion = 1;

  // WRITE MESSAGE
  WriteSchMessage("\n");
  WriteSchMessage("Computing Reynolds Stresses...");

  // Eval Turbulent Viscosity in current Scan
  MRITurbBoussinesq* turbModel = new MRITurbBoussinesq();
  turbModel->EvalTurbulentViscosity(this,turbNu);

  // Eval Turbulent Kinetic Energy in current Scan
  turbModel->EvalTurbulentKineticEnergy(this,turbK);

  // Cycle through all cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){

    // Eval First Derivative Tensor
    EvalSpaceDerivs(loopA,firstDerivs,secondDerivs);

    if(ReynoldsCriterion == 0){
      // COMPLETE STRESSES
      // RXXs
      cellPoints[loopA].ReStress[0] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][0]+firstDerivs[0][0]))-((2.0)/(3.0))*turbK[loopA];
      //cellPoints[loopA].ReStress[0] = turbK[loopA];
      // RXY
      cellPoints[loopA].ReStress[1] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][1]+firstDerivs[1][0]));
      // RXZ
      cellPoints[loopA].ReStress[2] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][2]+firstDerivs[2][0]));
      // RYY
      cellPoints[loopA].ReStress[3] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[1][1]+firstDerivs[1][1]))-((2.0)/(3.0))*turbK[loopA];
      // RYZ
      cellPoints[loopA].ReStress[4] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[1][2]+firstDerivs[2][1]));
      // RZZ
      cellPoints[loopA].ReStress[5] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[2][2]+firstDerivs[2][2]))-((2.0)/(3.0))*turbK[loopA];
    }else if(ReynoldsCriterion == 1){
      // ONLY VELOCITY GRADIENTS
      // RXXs
      cellPoints[loopA].ReStress[0] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][0]+firstDerivs[0][0]));
      // RXY
      cellPoints[loopA].ReStress[1] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][1]+firstDerivs[1][0]));
      // RXZ
      cellPoints[loopA].ReStress[2] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[0][2]+firstDerivs[2][0]));
      // RYY
      cellPoints[loopA].ReStress[3] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[1][1]+firstDerivs[1][1]));
      // RYZ
      cellPoints[loopA].ReStress[4] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[1][2]+firstDerivs[2][1]));
      // RZZ
      cellPoints[loopA].ReStress[5] = 2.0*turbNu[loopA]*(0.5*(firstDerivs[2][2]+firstDerivs[2][2]));
    }else if(ReynoldsCriterion == 2){
      // ONLY HYDROSTATIC PART
      // RXXs
      cellPoints[loopA].ReStress[0] = -((2.0)/(3.0))*turbK[loopA];
      // RXY
      cellPoints[loopA].ReStress[1] = 0.0;
      // RXZ
      cellPoints[loopA].ReStress[2] = 0.0;
      // RYY
      cellPoints[loopA].ReStress[3] = -((2.0)/(3.0))*turbK[loopA];
      // RYZ
      cellPoints[loopA].ReStress[4] = 0.0;
      // RZZ
      cellPoints[loopA].ReStress[5] = -((2.0)/(3.0))*turbK[loopA];
    }
  }

  // SET THE FLAG FOR AVAILABILITY OF REYNOLDS STRESSES
  hasReynoldsStress = true;

  // Delete first and second derivatives
  for(int loopA=0;loopA<3;loopA++){
    delete [] firstDerivs[loopA];
    delete [] secondDerivs[loopA];
  }
  delete [] firstDerivs;
  delete [] secondDerivs;
  delete [] turbNu;
  delete [] turbK;

  // DONE
  WriteSchMessage("Done.\n");
}



