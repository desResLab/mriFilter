#include <stdio.h>
#include <math.h>
#include "mriExpansion.h"

// ===========
// CONSTRUCTOR
// ===========
MRIExpansion::MRIExpansion(int totVortex){
  totalVortices = totVortex;
  constantFluxCoeff = new double[3];
  vortexCoeff = new double[totVortex];
  // Initialize Constant Flux
  for(int loopA=0;loopA<3;loopA++){
    constantFluxCoeff[loopA] = 0.0;
  }
  // Initialize Vortex Array
  for(int loopA=0;loopA<totVortex;loopA++){
    vortexCoeff[loopA] = 0.0;
  }
}

// ================
// COPY CONSTRUCTOR
// ================
MRIExpansion::MRIExpansion(MRIExpansion* otherExp){
  totalVortices = otherExp->totalVortices;
  constantFluxCoeff = new double[3];
  vortexCoeff = new double[totalVortices];
  // Initialize Constant Flux
  for(int loopA=0;loopA<3;loopA++){
    constantFluxCoeff[loopA] = otherExp->constantFluxCoeff[loopA];
  }
  // Initialize Vortex Array
  for(int loopA=0;loopA<totalVortices;loopA++){
    vortexCoeff[loopA] = otherExp->vortexCoeff[loopA];
  }
}

// ==========================
// CONSTRUCT FROM STD::VECTOR
// ==========================
void MRIExpansion::FillFromVector(std::vector<double> Expansion){
  // INITIALIZE VALUES
  for(int loopA=0;loopA<3;loopA++){
    constantFluxCoeff[loopA] = Expansion[loopA];
  }
  for(unsigned int loopA=3;loopA<Expansion.size();loopA++){
    vortexCoeff[loopA-3] = Expansion[loopA];
  }
}

// Distructor
MRIExpansion::~MRIExpansion(){
  delete[] constantFluxCoeff;
  delete[] vortexCoeff;
}

// ===========================
// Apply ratio-based Threshold
// ===========================
void MRIExpansion::ApplyVortexThreshold(double ratio){
  // Init
  double maxCoeff = 0.0;
  double currValue = 0.0;
  // Find Max Value
  for(int loopA=0;loopA<totalVortices;loopA++){
    currValue = fabs(vortexCoeff[loopA]);
    if(currValue>maxCoeff){
      maxCoeff = currValue;
    }
  }
  // Set Threshold
  double currThreshold = maxCoeff * ratio;
  // Apply Threshold
  for(int loopA=0;loopA<totalVortices;loopA++){
    currValue = fabs(vortexCoeff[loopA]);
    //if(currValue<=currThreshold){
    if(currValue<currThreshold){
      vortexCoeff[loopA] = 0.0;
    }
  }
}

// =============
// PRINT TO FILE
// =============
void MRIExpansion::WriteToFile(std::string outFile){
  // Open Output File
  FILE* fid;
  fid = fopen(outFile.c_str(),"w");
  // Write Constant Flux Components
  for(int loopA=0;loopA<3;loopA++){
    fprintf(fid,"%d,%15.6e\n",loopA,constantFluxCoeff[loopA]);
  }
  // Write Vortex Component
  for(int loopA=0;loopA<totalVortices;loopA++){
    fprintf(fid,"%d,%15.6e\n",loopA,vortexCoeff[loopA]);
  }
  // Close Output file
  fclose(fid);
}

// ==========
// GET 2 NORM
// ==========
double MRIExpansion::Get2Norm(bool onlyVortex){
  // Init
  double currNorm = 0.0;
  // Consider Constant Flux Coefficient
  if(!onlyVortex){
    for(int loopA=0;loopA<3;loopA++){
      currNorm += constantFluxCoeff[loopA] * constantFluxCoeff[loopA];
    }
  }
  // Consider Vortices
  for(int loopA=0;loopA<totalVortices;loopA++){
    currNorm += vortexCoeff[loopA] * vortexCoeff[loopA];
  }
  return sqrt(currNorm);
}



