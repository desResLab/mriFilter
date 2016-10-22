# include "mriOperation.h"

// DATA MEMBER
MRIOpApplySmoothing::MRIOpApplySmoothing(int filterNumIterations, int filterType, int filterOrder){
  this->numIterations = filterNumIterations;
  this->filterType = filterType;
  this->filterOrder = filterOrder;
}

// SCALE VELOCITIES
void MRIOpScaleVelocities::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->scaleVelocities(factor);
}

// SCALE CELLS
void MRIOpScaleCells::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->scalePositions(origin,factor);
}

// SAVE SEQUENCE AT CURRENT STATE
void MRIOpSaveState::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->saveVelocity();
}

// APPLY NOISE
void MRIOpApplyNoise::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->applyNoise(noiseIntensity);
}

// APPLY SMOOTHING FILTER
void MRIOpApplySmoothing::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->applyMedianFilter(kQtyVelocityX,numIterations,filterOrder,filterType,thresholdCriteria);
  seq->applyMedianFilter(kQtyVelocityY,numIterations,filterOrder,filterType,thresholdCriteria);
  seq->applyMedianFilter(kQtyVelocityZ,numIterations,filterOrder,filterType,thresholdCriteria);
}

// CLEAN NORMAL COMPONENT ON BOUNDARY
void MRIOpCleanNormalComponentOnBoundary::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->cleanNormalComponentOnBoundary(thresholdCriteria);
}

// INTERPOLATE BOUNDARY VELOCITIES
void MRIOpInterpolateVelocityOnBoundary::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->interpolateBoundaryVelocities(thresholdCriteria);
}

// INTERPOLATE BOUNDARY VELOCITIES
void MRIOpApplySolenoidalFilter::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->applySMPFilter(comm, applyBCFilter, 
                      thresholdCriteria,
                      itTol,maxIt,
                      useConstantPatterns);
}

// APPLY THRESHOLD
void MRIOpApplyThreshold::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->applyThresholding(thresholdCriteria);
}

// COMPUTE VORTEX CRITERIA
void MRIOpComputeVortexCriteria::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  // Eval Popular Vortex Criteria
  seq->evalVortexCriteria(thresholdCriteria);
  // Compute Vorticity
  if(computeVorticity){
    seq->evalVorticity(thresholdCriteria);
  }
  // Compute Enstrophy
  if(computeEnstrophy){
    seq->evalEnstrophy(thresholdCriteria);
  }
}

// COMPUTE SMP VORTEX CRITERION
void MRIOpComputeSMPVortexCriteria::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->evalSMPVortexCriteria();
}

// WRITE EXPANSION COEFFICIENTS
void MRIOpWriteExpansionCoefficients::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->writeExpansionFile(string(outputFileName + "_expCoeff"));
}

// EXPORT FOR FINITE ELEMENT POISSON SOLVER
void MRIOpExportForPoissonSolver::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->exportForPoisson(fileName,density,viscosity,thresholdCriteria,
                        PPE_IncludeAccelerationTerm,
                        PPE_IncludeAdvectionTerm,
                        PPE_IncludeDiffusionTerm,
                        PPE_IncludeReynoldsTerm,
                        readMuTFromFile,
                        muTFile,
                        smagorinskyCoeff);
}

// EXPORT FOR DISTANCE COMPUTATION
void MRIOpExportForDistanceSolver::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  seq->exportForDistancing(distanceFileName,thresholdCriteria);
}

// DATA MEMBER
void MRIOpEvalStatistics::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  // Apply Factors to limitBox
  double xFactor = 1.0;
  double yFactor = 1.0;
  double zFactor = 1.0;
  MRIUtils::applyLimitBoxFactors(xFactor,yFactor,zFactor,limitBox);
  
  // Allocate Bin Arrays
  MRIDoubleVec binCenters(numberOfBins);
  MRIDoubleVec binValues(numberOfBins);
  // Eval Single PDFs
  // FIRST
  seq->getScan(0)->evalScanPDF(kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
  MRIUtils::printBinArrayToFile(statFileNameFirst,numberOfBins,binCenters,binValues);

  // SECOND
  seq->getScan(0)->evalScanPDF(kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
  MRIUtils::printBinArrayToFile(statFileNameSecond,numberOfBins,binCenters,binValues);

  // DIFFERENCE
  seq->evalScanDifferencePDF(1,0,kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
  MRIUtils::printBinArrayToFile(statFileNameDiff,numberOfBins,binCenters,binValues);
}

// COMPUTE SCAN MATRICES
void MRIOpComputeScanMatrices::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){
  // VAR
  int totalERows = 0;
  int totalECols = 0;
  MRIDoubleMat EMat;
  int totalDRows = 0;
  int totalDCols = 0;
  MRIDoubleMat DMat;
  int totalStarRows = 0;
  int totalStarCols = 0;
  MRIDoubleMat StarMatrix;

  // SET PARAMETERS
  bool isIsotropic = true;
  
  // Set Template Parameters
  MRIDoubleVec params(8);
  params[0] = 5;
  params[1] = 5;
  params[2] = 5;
  params[3] = 1.0;
  params[4] = 1.0;
  params[5] = 1.0;
  params[6] = 0.0;
  params[7] = 1.0;
  
  // SET ISOTROPIC OR ANISOTROPIC CASE
  if(isIsotropic){
    // ISOTROPIC
    seq->createSampleCase(kConstantFlow,params);
  }else{
    // ANISOTROPIC
    params[3] = 1.0;
    params[4] = 2.0;
    params[5] = 3.0;
    seq->createSampleCase(kConstantFlow,params);
  }

  // RETRIEVE OPERATORS IN MATRIX FORM
  // ENCODING
  seq->getScan(0)->assembleEncodingMatrix(totalERows,totalECols,EMat);
  MRIUtils::printMatrixToFile("EncodingMat.dat",EMat);
  // DECODING
  seq->getScan(0)->assembleDecodingMatrix(totalDRows,totalDCols,DMat);
  MRIUtils::printMatrixToFile("DecodingMat.dat",DMat);
  // VORTEX FRAME MATRIX
  seq->getScan(0)->assembleStarMatrix(totalStarRows,totalStarCols,StarMatrix);
  MRIUtils::printMatrixToFile("StarMat.dat",StarMatrix);
}

// =======================
// SHOW FACE FLUX PATTERNS
// =======================
void MRIOpShowFaceFluxPatterns::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){

  bool isIsotropic = true;

  MRIDoubleVec params(7);
  params[0] = 5.0;
  params[1] = 5.0;
  params[2] = 5.0;
  params[3] = 1.0;
  params[4] = 1.0;
  params[5] = 1.0;
  params[6] = 1.0;
  if(isIsotropic){
    // ISOTROPIC
    seq->createSampleCase(kConstantFlow,params);
  }else{
   // ANISOTROPIC
   params[3] = 1.0;
   params[4] = 2.0;
   params[5] = 3.0;
   seq->createSampleCase(kConstantFlow,params);
  }

  // READ FACE FLUXES FROM FILE
  int totalRows = 0;
  int totalCols = 0;
  MRIDoubleMat faceFluxMat;
  MRIUtils::readMatrixFromFile(faceFluxFileName,totalRows,totalCols,faceFluxMat);

  // COPY THE INTERESTING COLUMN
  MRIDoubleVec faceFluxVec(totalRows);
  for(int loop0=0;loop0<100;loop0++){
    int selectedCol = loop0;
    for(int loopA=0;loopA<totalRows;loopA++){
      faceFluxVec[loopA] = faceFluxMat[loopA][selectedCol];
    }

    // TRASFORM FACE FLUXES IN VELOCITIES
    seq->getScan(0)->recoverCellVelocitiesRT0(false,faceFluxVec);
    
    // UPDATE VELOCITIES
    seq->getScan(0)->updateVelocities();
  }
}

// EVALUATE CONCENTRATION GRADIENT
void MRIOpEvalConcentrationGradient::processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq){

  bool doPressureSmoothing = false;

  // EVAL REYNOLDS STRESSES AND PRESSURE GRADIENTS
  seq->getScan(0)->computeQuantityGradient(kQtyConcentration);

  // EVAL RELATIVE PRESSURE
  seq->computeRelativePressure(doPressureSmoothing);

}