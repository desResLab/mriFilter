# include "mriOperation.h"

// SCALE VELOCITIES
void MRIOpScaleVelocities::processSequence(MRISequence* seq){
  seq->scaleVelocities(factor);
}

// SCALE CELLS
void MRIOpScaleCells::processSequence(MRISequence* seq){
  seq->scalePositions(origin,factor);
}

// SAVE SEQUENCE AT CURRENT STATE
void MRIOpSaveState::processSequence(MRISequence* seq){
  seq->saveVelocity();
}

// APPLY NOISE
void MRIOpApplyNoise::processSequence(MRISequence* seq){
  seq->applyNoise(noiseIntensity);
}

// APPLY SMOOTHING FILTER
void MRIOpApplySmoothing::processSequence(MRISequence* seq){
  seq->ApplyMedianFilter(kQtyVelocityX,opts->filterNumIterations,opts->medianFilterOrder,opts->medianFilterType,opts->thresholdCriteria);
  seq->ApplyMedianFilter(kQtyVelocityY,opts->filterNumIterations,opts->medianFilterOrder,opts->medianFilterType,opts->thresholdCriteria);
  seq->ApplyMedianFilter(kQtyVelocityZ,opts->filterNumIterations,opts->medianFilterOrder,opts->medianFilterType,opts->thresholdCriteria);
}

// CLEAN NORMAL COMPONENT ON BOUNDARY
void MRIOpCleanNormalComponentOnBoundary::processSequence(MRISequence* seq){
  seq->cleanNormalComponentOnBoundary();
}

// INTERPOLATE BOUNDARY VELOCITIES
class MRIOpInterpolateVelocityOnBoundary::processSequence(MRISequence* seq){
  seq->InterpolateBoundaryVelocities();
}

// INTERPOLATE BOUNDARY VELOCITIES
class MRIOpApplySolenoidalFilter::processSequence(MRISequence* seq){
  seq->ApplySMPFilter(opts,applyBCFilter,comm);
}

// APPLY THRESHOLD
class MRIOpApplyThreshold::processSequence(MRISequence* seq){
  seq->ApplyThresholding(opts->thresholdCriteria);
}

// COMPUTE VORTEX CRITERIA
class MRIOpComputeVortexCriteria::processSequence(MRISequence* seq){
  // Eval Popular Vortex Criteria
  seq->EvalVortexCriteria(opts->thresholdCriteria);
  // Compute Vorticity
  if(computeVorticity){
    seq->EvalVorticity(opts->thresholdCriteria);
  }
  // Compute Enstrophy
  if(computeEnstrophy){
    seq->EvalEnstrophy(opts->thresholdCriteria);
  }
}

// COMPUTE SMP VORTEX CRITERION
class MRIOpComputeSMPVortexCriteria::processSequence(MRISequence* seq){
  // if(opts->applySMPFilter && opts->evalSMPVortexCriterion)
  seq->EvalSMPVortexCriteria();
}

// WRITE EXPANSION COEFFICIENTS
class MRIOpWriteExpansionCoefficients::processSequence(MRISequence* seq){
  if(opts->applySMPFilter && opts->saveExpansionCoeffs){
    seq->eriteExpansionFile(std::string(opts->outputFileName + "_expCoeff"));
  }
}

// EXPORT FOR FINITE ELEMENT POISSON SOLVER
class MRIOpExportForPoissonSolver::processSequence(MRISequence* seq){
  seq->ExportForPoisson(poissonFileName,density,viscosity,thresholdCriteria,
                        PPE_IncludeAccelerationTerm,
                        PPE_IncludeAdvectionTerm,
                        PPE_IncludeDiffusionTerm,
                        PPE_IncludeReynoldsTerm,
                        readMuTFromFile,
                        muTFile,
                        smagorinskyCoeff);
}

// EXPORT FOR DISTANCE COMPUTATION
class MRIOpExportForDistanceSolver::processSequence(MRISequence* seq){
  seq->ExportForDistancing(distanceFileName,thresholdCriteria);
}
