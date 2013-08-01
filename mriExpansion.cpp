#include <math.h>
#include "mriExpansion.h"

// Constructor
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
    if(currValue<=currThreshold){
      vortexCoeff[loopA] = 0.0;
    }
  }
}

