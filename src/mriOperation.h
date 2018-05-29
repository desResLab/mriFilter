#ifndef MRIOPERATION_H
#define MRIOPERATION_H

# include "mriSequence.h"
# include "mriCommunicator.h"
# include "mriThresholdCriteria.h"

class MRISequence;

// GENERIC MRI OPERATION
class MRIOperation{
  public:
    MRIOperation();
    virtual ~MRIOperation();

    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq) = 0;
};

// SCALE VELOCITIES
class MRIOpScaleVelocities: public MRIOperation{
  public:
    double factor;

    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// SCALE CELLS
class MRIOpScaleCells: public MRIOperation{
  public:
    MRIDoubleVec origin;
    double factor;

    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// SAVE SEQUENCE AT CURRENT STATE
class MRIOpSaveState: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// APPLY NOISE
class MRIOpApplyNoise: public MRIOperation{
  public:
    double noiseIntensity;
    double seed;
    // CONSTRUCTOR
    MRIOpApplyNoise(double noiseIntensity, double seed);
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// APPLY SMOOTHING FILTER
class MRIOpApplySmoothing: public MRIOperation{
  public:
    // Filtering Options
    int numIterations;
    int filterType;
    int filterOrder;

    // DATA MEMBER
    MRIOpApplySmoothing(int filterNumIterations, int filterType, int filterOrder);
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// CLEAN NORMAL COMPONENT ON BOUNDARY
class MRIOpCleanNormalComponentOnBoundary: public MRIOperation{
  public:
        
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// INTERPOLATE BOUNDARY VELOCITIES
class MRIOpInterpolateVelocityOnBoundary: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// INTERPOLATE BOUNDARY VELOCITIES
class MRIOpApplySolenoidalFilter: public MRIOperation{
  public:
    bool applyBCFilter;
    bool useConstantPatterns;
    double itTol;
    int maxIt;

    // CONSTRUCTOR
    MRIOpApplySolenoidalFilter(bool applyBCFilter,bool useConstantPatterns,double itTol,int maxIt);
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// APPLY THRESHOLD
class MRIOpApplyThreshold: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// COMPUTE VORTEX CRITERIA
class MRIOpComputeVortexCriteria: public MRIOperation{
  public:
    bool computeVorticity;
    bool computeEnstrophy;
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// COMPUTE SMP VORTEX CRITERION
class MRIOpComputeSMPVortexCriteria: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// WRITE EXPANSION COEFFICIENTS
class MRIOpWriteExpansionCoefficients: public MRIOperation{
  public:
    string outputFileName;
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// EXPORT FOR FINITE ELEMENT POISSON SOLVER
class MRIOpExportForPoissonSolver: public MRIOperation{
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
    MRIOpExportForPoissonSolver(string fileName,
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
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// EXPORT FOR DISTANCE COMPUTATION
class MRIOpExportForDistanceSolver: public MRIOperation{
  public:
    string distanceFileName;
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// EVALUATE STATISTICS OF TWO SCANS
class MRIOpEvalStatistics: public MRIOperation{
  public:
    int numberOfBins;
    bool useBox;
    MRIDoubleVec limitBox;
    string statFileNameFirst;
    string statFileNameSecond;
    string statFileNameDiff;

    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// COMPUTE SCAN MATRICES
class MRIOpComputeScanMatrices: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// SHOW FACE FLUX PATTERNS
class MRIOpShowFaceFluxPatterns: public MRIOperation{
  public:
    string faceFluxFileName;
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

// EVALUATE CONCENTRATION GRADIENT
class MRIOpEvalConcentrationGradient: public MRIOperation{
  public:
    // DATA MEMBER
    virtual void processSequence(MRICommunicator* comm, MRIThresholdCriteria* thresholdCriteria, MRISequence* seq);
};

#endif // MRIOPERATION_H