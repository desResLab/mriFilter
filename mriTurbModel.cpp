#include "mriTurbModel.h"

// CONSTRUCTOR FOR BASE CLASS
MRITurbModel::MRITurbModel(){
}

// DISTRUCTOR FOR BASE CLASS
MRITurbModel::~MRITurbModel(){
}

// ========================
// BOUSSINESQ APPROXIMATION
// ========================
// EVAL TURBULENT VISCOSITY
void MRITurbBoussinesq::EvalTurbulentViscosity(MRIScan* currScan, double* nuT){
  // Loop through all cells
  for(int loopA=0;loopA<currScan->totalCellPoints;loopA++){
    // Careful SI Units
    nuT[loopA] = 0.01*0.5*0.0058;
  }
}
// EVAL TURBULENT KINETIC ENERGY
void MRITurbBoussinesq::EvalTurbulentKineticEnergy(MRIScan* currScan, double* kTurb){
  // First and Second Cell Derivatives
  double** firstDerivs = new double*[3];
  double** secondDerivs = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    firstDerivs[loopA] = new double[3];
    secondDerivs[loopA] = new double[3];
  }
  // Loop through all cells
  double tensorProduct = 0.0;
  for(int loopA=0;loopA<currScan->totalCellPoints;loopA++){
    // Eval First and Derivative Tensor
    currScan->EvalSpaceDerivs(loopA,firstDerivs,secondDerivs);
    // Eval Tensor Product
    tensorProduct = 0.0;
    for(int loopB=0;loopB<3;loopB++){
      for(int loopC=0;loopC<3;loopC++){
        tensorProduct += firstDerivs[loopB][loopC]*firstDerivs[loopB][loopC];
      }
    }
    // Careful SI Units
    kTurb[loopA] = 2.0*tensorProduct*0.01;
  }
  // Delete first and second derivatives
  for(int loopA=0;loopA<3;loopA++){
    delete [] firstDerivs[loopA];
    delete [] secondDerivs[loopA];
  }
  delete [] firstDerivs;
  delete [] secondDerivs;
}



