#ifndef MRIOPERATION_H
#define MRIOPERATION_H

# include "mriSequence.h"

class MRISequence;

// GENERIC MRI OPERATION
class MRIOperation{
  public:
    MRIOperation();
    ~MRIOperation();

    // DATA MEMBER
    virtual void processSequence(MRISequence* seq) = 0;
};

// SCALE VELOCITIES
class MRIOpScaleVelocities: public MRIOperation{
  public:
    double factor;

    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// SCALE CELLS
class MRIOpScaleCells: public MRIOperation{
  public:
    double origin[3];
    double factor;

    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// SAVE SEQUENCE AT CURRENT STATE
class MRIOpSaveState: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// APPLY NOISE
class MRIOpApplyNoise: public MRIOperation{
  public:
    double noiseIntensity;

    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// APPLY SMOOTHING FILTER
class MRIOpApplySmoothing: public MRIOperation{
  public:
    // Filtering Options
    int numIterations;
    int filterType;
    int filterOrder;

    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// CLEAN NORMAL COMPONENT ON BOUNDARY
class MRIOpCleanNormalComponentOnBoundary: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// INTERPOLATE BOUNDARY VELOCITIES
class MRIOpInterpolateVelocityOnBoundary: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// INTERPOLATE BOUNDARY VELOCITIES
class MRIOpApplySolenoidalFilter: public MRIOperation{
  public:
    bool applyBCFilter;
    bool useConstantPatterns;
    double itTol;
    int maxIt;

    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// APPLY THRESHOLD
class MRIOpApplyThreshold: public MRIOperation{
  public:
    int thresholdQty;
    int thresholdType;
    double thresholdValue; 

    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// COMPUTE VORTEX CRITERIA
class MRIOpComputeVortexCriteria: public MRIOperation{
  public:
    bool computeVorticity;
    bool computeEnstrophy;

    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// COMPUTE SMP VORTEX CRITERION
class MRIOpComputeSMPVortexCriteria: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// WRITE EXPANSION COEFFICIENTS
class MRIOpWriteExpansionCoefficients: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// EXPORT FOR FINITE ELEMENT POISSON SOLVER
class MRIOpExportForPoissonSolver: public MRIOperation{
  public:
    // Pressure Gradient Components to include
    bool PPE_IncludeAccelerationTerm;
    bool PPE_IncludeAdvectionTerm;
    bool PPE_IncludeDiffusionTerm;
    bool PPE_IncludeReynoldsTerm;
    // Turbulent Viscosity
    bool readMuTFromFile;
    string muTFile;
    double smagorinskyCoeff;

    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

// EXPORT FOR DISTANCE COMPUTATION
class MRIOpExportForDistanceSolver: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRISequence* seq);
};

#endif // MRIOPERATION_H