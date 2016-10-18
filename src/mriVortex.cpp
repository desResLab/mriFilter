#include "mriStructuredScan.h"
#include "mriUtils.h"
#include "mriException.h"

// ===================================================
// CRITERIAL FOR IDENTIFICATION OF VORTICAL STRUCTURES
// ===================================================

// ===========================================
// EVAL DECOMPOSITION OF THE VELOCITY GRADIENT
// ===========================================
void MRIScan::evalCellVelocityGradientDecomposition(int currentCell, const MRIDoubleMat& velGrad, MRIDoubleMat& deformation, MRIDoubleMat& rotation){
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      deformation[loopA][loopB] = 0.5*(velGrad[loopA][loopB] + velGrad[loopB][loopA]);
      rotation[loopA][loopB] = 0.5*(velGrad[loopA][loopB] - velGrad[loopB][loopA]);
    }
  }
}

// ================
// EVAL Q CRITERION
// ================
double MRIScan::evalCellQCriterion(int currentCell, const MRIDoubleMat& deformation, const MRIDoubleMat& rotation){
  double traceDEF = 0.0;
  double traceROT = 0.0;
  double opResultDEF = 0.0;
  double opResultROT = 0.0;
  // MULTIPLY
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      opResultDEF = 0.0;
      opResultROT = 0.0;
      for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
        opResultDEF += deformation[loopA][loopC] * deformation[loopB][loopC];
        opResultROT += rotation[loopA][loopC] * rotation[loopB][loopC];
      }
      if(loopA == loopB){
        traceDEF += opResultDEF;
        traceROT += opResultROT;
      }
    }
  }
  // COMPUTE CRITERION
  return 0.5*(traceROT-traceDEF);
}

// =======================
// EVAL LAMBDA-2 CRITERION
// =======================
double MRIScan::evalCellL2Criterion(int currentCell, const MRIDoubleMat& deformation, const MRIDoubleMat& rotation){
  double mat[3][3];
  double eig[3];
  double opResultDEF = 0.0;
  double opResultROT = 0.0;

  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      opResultDEF = 0.0;
      opResultROT = 0.0;
      for(int loopC=0;loopC<kNumberOfDimensions;loopC++){
        opResultDEF += deformation[loopA][loopC] * deformation[loopC][loopB];
        opResultROT += rotation[loopA][loopC] * rotation[loopC][loopB];
      }
      mat[loopA][loopB] = opResultDEF + opResultROT;
    }
  }
  // Compute the second eigenvalue
  MRIUtils::compute3x3MatrixEigenvals(mat,eig);
  return eig[1];
}

// ====================
// EVAL DELTA CRITERION
// ====================
double MRIScan::evalCellDeltaCriterion(int currentCell, const MRIDoubleMat& deformation, const MRIDoubleMat& rotation, const MRIDoubleMat& velGradient){
  double qVal = EvalCellQCriterion(currentCell,deformation,rotation);
  double detVelGrad = (velGradient[0][0] * (velGradient[1][1] * velGradient[2][2] - velGradient[2][1] * velGradient[1][2])
                      -velGradient[1][0] * (velGradient[0][1] * velGradient[2][2] - velGradient[2][1] * velGradient[0][2])
                      +velGradient[2][0] * (velGradient[0][1] * velGradient[1][2] - velGradient[1][1] * velGradient[0][2]));
  return 0.5*detVelGrad*detVelGrad + (1.0/27.0)*qVal*qVal*qVal;
}

// ====================
// EVAL DELTA CRITERION
// ====================
double MRIScan::evalCellVortexCriteria(int currentCell,int criteriaType, const MRIDoubleMat& velGradient, const MRIDoubleMat& deformation, const MRIDoubleMat& rotation){
  double critRes = 0.0;
  // Vortex Criterion
  switch(criteriaType){
    case kVortexQ:
      critRes = EvalCellQCriterion(currentCell,deformation,rotation);
      break;
    case kVortexL2:
      critRes = EvalCellL2Criterion(currentCell,deformation,rotation);
      break;
    case kVortexDelta:
      critRes = EvalCellDeltaCriterion(currentCell,velGradient,deformation,rotation);
      break;
    default:
      throw MRIVortexException("Error: Invalid Vortex Criterion");
      critRes = 0.0;
      break;
  }
  return critRes;
}

// =============================
// EVALUATION OF VORTEX CRITERIA
// =============================
void MRIScan::evalVortexCriteria(MRIThresholdCriteria* threshold){
  // Velocity Gradient and Hessian
  // First and Second Derivatives
  double** deformation = new double*[kNumberOfDimensions];
  double** rotation = new double*[kNumberOfDimensions];
  double** firstDerivs = new double*[kNumberOfDimensions];
  double** secondDerivs = new double*[kNumberOfDimensions];
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA] = new double[kNumberOfDimensions];
    secondDerivs[loopA] = new double[kNumberOfDimensions];
    deformation[loopA] = new double[kNumberOfDimensions];
    rotation[loopA] = new double[kNumberOfDimensions];
  }
  // Create three new outputs
  MRIOutput out1("QCriterion",1);
  MRIOutput out2("L2Criterion",1);
  MRIOutput out3("DeltaCriterion",1);
  // Loop on cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    EvalSpaceDerivs(loopA,threshold,firstDerivs,secondDerivs);
    EvalCellVelocityGradientDecomposition(loopA,deformation,rotation,firstDerivs);
    // Store Criteria
    out1.values.push_back(EvalCellVortexCriteria(loopA,kVortexQ,deformation,rotation,firstDerivs));
    out2.values.push_back(EvalCellVortexCriteria(loopA,kVortexL2,deformation,rotation,firstDerivs));
    out3.values.push_back(EvalCellVortexCriteria(loopA,kVortexDelta,deformation,rotation,firstDerivs));
  }
  // Add output quantities
  outputs.push_back(out1);
  outputs.push_back(out2);
  outputs.push_back(out3);
  // Deallocate
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    delete [] firstDerivs[loopA];
    delete [] secondDerivs[loopA];
    delete [] deformation[loopA];
    delete [] rotation[loopA];
  }
  delete [] firstDerivs;
  delete [] secondDerivs;
  delete [] deformation;
  delete [] rotation;
}

// COMPUTATION OF VORTICITY
void computeVorticity(double** firstDerivs,double* auxVector){
  auxVector[0] = firstDerivs[1][2] - firstDerivs[2][1];
  auxVector[1] = firstDerivs[2][0] - firstDerivs[0][2];
  auxVector[2] = firstDerivs[0][1] - firstDerivs[1][0];
}

// ==============
// EVAL VORTICITY
// ==============
void MRIScan::evalVorticity(MRIThresholdCriteria* threshold){
  // Allocate derivatives
  double** firstDerivs = new double*[kNumberOfDimensions];
  double** secondDerivs = new double*[kNumberOfDimensions];
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA] = new double[kNumberOfDimensions];
    secondDerivs[loopA] = new double[kNumberOfDimensions];
  }
  double vort[3];
  MRIOutput out1("Vorticity",3);
  // Loop on cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    EvalSpaceDerivs(loopA,threshold,firstDerivs,secondDerivs);
    // Store Criteria
    ComputeVorticity(firstDerivs,vort);
    out1.values.push_back(vort[0]);
    out1.values.push_back(vort[1]);
    out1.values.push_back(vort[2]);
  }
  // Add output quantities
  outputs.push_back(out1);
  // Deallocate
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    delete [] firstDerivs[loopA];
    delete [] secondDerivs[loopA];
  }
  delete [] firstDerivs;
  delete [] secondDerivs;
}

// ==============
// EVAL ENSTROPHY
// ==============
void MRIScan::evalEnstrophy(MRIThresholdCriteria* threshold){
  // Allocate derivatives  
  double** firstDerivs = new double*[kNumberOfDimensions];
  double** secondDerivs = new double*[kNumberOfDimensions];
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    firstDerivs[loopA] = new double[kNumberOfDimensions];
    secondDerivs[loopA] = new double[kNumberOfDimensions];
  }
  double vec[3];
  double value = 0.0;
  MRIOutput out1("Enstrophy",1);
  // Loop on cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    EvalSpaceDerivs(loopA,threshold,firstDerivs,secondDerivs);
    // Store Criteria
    ComputeVorticity(firstDerivs,vec);
    // Square Modulus
    out1.values.push_back(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    // Compute Sum
    value += vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
  }
  // Add output quantities
  outputs.push_back(out1);
  // Print
  printf("Total Enstrophy: %e\n",value);
  // Deallocate
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    delete [] firstDerivs[loopA];
    delete [] secondDerivs[loopA];
  }
  delete [] firstDerivs;
  delete [] secondDerivs;
}



