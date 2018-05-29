#ifndef MRIOPERATION_H
#define MRIOPERATION_H

# include "mriSequence.h"
# include "mriCommunicator.h"
# include "mriThresholdCriteria.h"

class mriSequence;

// GENERIC MRI OPERATION
class mriOperation{
  public:
    mriOperation();
    virtual ~mriOperation();

    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq) = 0;
};

// SCALE VELOCITIES
class mriOpScaleVelocities: public mriOperation{
  public:
    double factor;

    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// SCALE CELLS
class mriOpScaleCells: public mriOperation{
  public:
    mriDoubleVec origin;
    double factor;

    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// SAVE SEQUENCE AT CURRENT STATE
class mriOpSaveState: public mriOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// APPLY NOISE
class mriOpApplyNoise: public mriOperation{
  public:
    double noiseIntensity;
    double seed;
    // CONSTRUCTOR
    mriOpApplyNoise(double noiseIntensity, double seed);
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// APPLY SMOOTHING FILTER
class mriOpApplySmoothing: public mriOperation{
  public:
    // Filtering Options
    int numIterations;
    int filterType;
    int filterOrder;

    // DATA MEMBER
    mriOpApplySmoothing(int filterNumIterations, int filterType, int filterOrder);
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// CLEAN NORMAL COMPONENT ON BOUNDARY
class mriOpCleanNormalComponentOnBoundary: public mriOperation{
  public:
        
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// INTERPOLATE BOUNDARY VELOCITIES
class mriOpInterpolateVelocityOnBoundary: public mriOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// INTERPOLATE BOUNDARY VELOCITIES
class mriOpApplySolenoidalFilter: public mriOperation{
  public:
    bool applyBCFilter;
    bool useConstantPatterns;
    double itTol;
    int maxIt;

    // CONSTRUCTOR
    mriOpApplySolenoidalFilter(bool applyBCFilter,bool useConstantPatterns,double itTol,int maxIt);
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// APPLY THRESHOLD
class mriOpApplyThreshold: public mriOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// COMPUTE VORTEX CRITERIA
class mriOpComputeVortexCriteria: public mriOperation{
  public:
    bool computeVorticity;
    bool computeEnstrophy;
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// COMPUTE SMP VORTEX CRITERION
class mriOpComputeSMPVortexCriteria: public mriOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// WRITE EXPANSION COEFFICIENTS
class mriOpWriteExpansionCoefficients: public mriOperation{
  public:
    string outputFileName;
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// EXPORT FOR FINITE ELEMENT POISSON SOLVER
class mriOpExportForPoissonSolver: public mriOperation{
  public:
    string fileName;
    double density;
    double viscosity;
    
    // Pressure Gradient Components to include
    bool PPE_IncludeAccelerationTerm;
    bool PPE_IncludeAdvectionTerm;
    bool PPE_IncludeDiffusionTerm;
    bool PPE_IncludeReynoldsTerm;
    
    // Turbulent Viscosity
    bool readMuTFromFile;
    string muTFile;
    double smagorinskyCoeff;

    // CONSTRUCTOR
    mriOpExportForPoissonSolver(string fileName,
                                double density,
                                double viscosity,
                                bool PPE_IncludeAccelerationTerm,
                                bool PPE_IncludeAdvectionTerm,
                                bool PPE_IncludeDiffusionTerm,
                                bool PPE_IncludeReynoldsTerm,
                                bool readMuTFromFile,
                                string muTFile,
                                double smagorinskyCoeff);

    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// EXPORT FOR DISTANCE COMPUTATION
class mriOpExportForDistanceSolver: public mriOperation{
  public:
    string distanceFileName;
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// EVALUATE STATISTICS OF TWO SCANS
class mriOpEvalStatistics: public mriOperation{
  public:
    int numberOfBins;
    bool useBox;
    mriDoubleVec limitBox;
    string statFileNameFirst;
    string statFileNameSecond;
    string statFileNameDiff;

    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// COMPUTE SCAN MATRICES
class mriOpComputeScanMatrices: public mriOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// SHOW FACE FLUX PATTERNS
class mriOpShowFaceFluxPatterns: public mriOperation{
  public:
    string faceFluxFileName;
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

// EVALUATE CONCENTRATION GRADIENT
class mriOpEvalConcentrationGradient: public mriOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(mriCommunicator* comm, mriThresholdCriteria* thresholdCriteria, mriSequence* seq);
};

#endif // mriOPERATION_H