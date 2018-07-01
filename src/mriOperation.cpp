# include "mriOperation.h"

// CONSTRUCTOR
mriOperation::mriOperation(){

}

// DISTRUCTOR
mriOperation::~mriOperation(){

}

// CONSTRUCTOR
mriOpExportForPoissonSolver::mriOpExportForPoissonSolver(string fileName,
                                                         double density,
                                                         double viscosity,
                                                         bool PPE_IncludeAccelerationTerm,
                                                         bool PPE_IncludeAdvectionTerm,
                                                         bool PPE_IncludeDiffusionTerm,
                                                         bool PPE_IncludeReynoldsTerm,
                                                         bool readMuTFromFile,
                                                         string muTFile,
                                                         double smagorinskyCoeff){
  this->fileName = fileName;
  this->density = density;
  this->viscosity = viscosity;
  this->PPE_IncludeAccelerationTerm = PPE_IncludeAccelerationTerm;
  this->PPE_IncludeAdvectionTerm = PPE_IncludeAdvectionTerm;
  this->PPE_IncludeDiffusionTerm = PPE_IncludeDiffusionTerm;
  this->PPE_IncludeReynoldsTerm = PPE_IncludeReynoldsTerm;
  this->readMuTFromFile = readMuTFromFile;
  this->muTFile = muTFile;
  this->smagorinskyCoeff = smagorinskyCoeff;
}

// DATA MEMBER
mriOpApplySmoothing::mriOpApplySmoothing(int filterNumIterations, int filterType, int filterOrder){
  this->numIterations = filterNumIterations;
  this->filterType = filterType;
  this->filterOrder = filterOrder;
}

// INITIALIZE APPLY NOISE OPERATION
mriOpApplyNoise::mriOpApplyNoise(double noiseIntensity, double seed){
  this->noiseIntensity = noiseIntensity;
  this->seed = seed;
}

// CONSTRUCTOR FOR SOLENOIDAL FILTER OPERATION
mriOpApplySolenoidalFilter::mriOpApplySolenoidalFilter(bool applyBCFilter,bool useConstantPatterns,double itTol,int maxIt){
  this->applyBCFilter = applyBCFilter;
  this->useConstantPatterns = useConstantPatterns;
  this->itTol = itTol;
  this->maxIt = maxIt;
}


// SCALE VELOCITIES
void mriOpScaleVelocities::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->scaleVelocities(factor);
}

// SCALE CELLS
void mriOpScaleCells::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->scalePositions(origin,factor);
}

// SAVE SEQUENCE AT CURRENT STATE
void mriOpSaveState::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->saveVelocity();
}

// APPLY NOISE
void mriOpApplyNoise::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->applyNoise(noiseIntensity,seed);
}

// APPLY SMOOTHING FILTER
void mriOpApplySmoothing::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->applyMedianFilter(kQtyVelocityX,numIterations,filterOrder,filterType,thresholdCriteria);
  seq->applyMedianFilter(kQtyVelocityY,numIterations,filterOrder,filterType,thresholdCriteria);
  seq->applyMedianFilter(kQtyVelocityZ,numIterations,filterOrder,filterType,thresholdCriteria);
}

// CLEAN NORMAL COMPONENT ON BOUNDARY
void mriOpCleanNormalComponentOnBoundary::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->cleanNormalComponentOnBoundary(thresholdCriteria);
}

// INTERPOLATE BOUNDARY VELOCITIES
void mriOpInterpolateVelocityOnBoundary::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->interpolateBoundaryVelocities(thresholdCriteria);
}

// INTERPOLATE BOUNDARY VELOCITIES
void mriOpApplySolenoidalFilter::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->applySMPFilter(comm, false, 
                      thresholdCriteria,
                      itTol,maxIt,
                      useConstantPatterns);
  if(applyBCFilter){
    seq->applySMPFilter(comm, true, 
                        thresholdCriteria,
                        itTol,maxIt,
                        useConstantPatterns);
  }
}

// APPLY THRESHOLD
void mriOpApplyThreshold::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->applyThresholding(thresholdCriteria);
}

// COMPUTE VORTEX CRITERIA
void mriOpComputeVortexCriteria::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
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
void mriOpComputeSMPVortexCriteria::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->evalSMPVortexCriteria();
}

// WRITE EXPANSION COEFFICIENTS
void mriOpWriteExpansionCoefficients::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->writeExpansionFile(string(outputFileName + "_expCoeff"));
}

// EXPORT FOR FINITE ELEMENT POISSON SOLVER
void mriOpExportForPoissonSolver::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
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
void mriOpExportForDistanceSolver::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  seq->exportForDistancing(distanceFileName,thresholdCriteria);
}

// DATA MEMBER
void mriOpEvalStatistics::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  // Apply Factors to limitBox
  double xFactor = 1.0;
  double yFactor = 1.0;
  double zFactor = 1.0;
  mriUtils::applyLimitBoxFactors(xFactor,yFactor,zFactor,limitBox);
  
  // Allocate Bin Arrays
  mriDoubleVec binCenters(numberOfBins);
  mriDoubleVec binValues(numberOfBins);
  // Eval Single PDFs
  // FIRST
  seq->getScan(0)->evalScanPDF(kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
  mriUtils::printBinArrayToFile(statFileNameFirst,numberOfBins,binCenters,binValues);

  // SECOND
  seq->getScan(0)->evalScanPDF(kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
  mriUtils::printBinArrayToFile(statFileNameSecond,numberOfBins,binCenters,binValues);

  // DIFFERENCE
  seq->evalScanDifferencePDF(1,0,kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
  mriUtils::printBinArrayToFile(statFileNameDiff,numberOfBins,binCenters,binValues);
}

// COMPUTE SCAN MATRICES
void mriOpComputeScanMatrices::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){
  // VAR
  int totalERows = 0;
  int totalECols = 0;
  mriDoubleMat EMat;
  int totalDRows = 0;
  int totalDCols = 0;
  mriDoubleMat DMat;
  int totalStarRows = 0;
  int totalStarCols = 0;
  mriDoubleMat StarMatrix;

  // SET PARAMETERS
  bool isIsotropic = true;
  
  // Set Template Parameters
  mriDoubleVec params(8);
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
  mriUtils::printMatrixToFile("EncodingMat.dat",EMat);
  // DECODING
  seq->getScan(0)->assembleDecodingMatrix(totalDRows,totalDCols,DMat);
  mriUtils::printMatrixToFile("DecodingMat.dat",DMat);
  // VORTEX FRAME MATRIX
  seq->getScan(0)->assembleStarMatrix(totalStarRows,totalStarCols,StarMatrix);
  mriUtils::printMatrixToFile("StarMat.dat",StarMatrix);
}

// =======================
// SHOW FACE FLUX PATTERNS
// =======================
void mriOpShowFaceFluxPatterns::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){

  bool isIsotropic = true;

  mriDoubleVec params(7);
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
  mriDoubleMat faceFluxMat;
  mriUtils::readMatrixFromFile(faceFluxFileName,totalRows,totalCols,faceFluxMat);

  // COPY THE INTERESTING COLUMN
  mriDoubleVec faceFluxVec(totalRows);
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
void mriOpEvalConcentrationGradient::processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq){

  bool doPressureSmoothing = false;

  // EVAL REYNOLDS STRESSES AND PRESSURE GRADIENTS
  seq->getScan(0)->computeQuantityGradient(kQtyConcentration);

}
