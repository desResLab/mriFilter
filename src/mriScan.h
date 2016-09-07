#ifndef MRISCAN_H
#define MRISCAN_H

#include <stdio.h>
#include <string>
#include <vector>

#include "mriCell.h"
#include "mriVolData.h"
#include "mriCellMaterial.h"
#include "mriSamplingOptions.h"
#include "mriStreamline.h"
#include "mriStreamlineOptions.h"
#include "mriThresholdCriteria.h"
#include "mriCommunicator.h"
#include "mriConstants.h"
#include "mriImagedata.h"
#include "mriExpansion.h"
#include "mriException.h"
#include "mriOptions.h"
#include "mriOutput.h"

using namespace std;

// =====================================================
// MAIN CLASS FOR BOTH UNSTRUCTURED AND STRUCTURED SCANS
// =====================================================
class MRIScan{
  public:
    // CONSTRUCTOR
    MRIScan(double currentTime);
    // COPY CONSTRUCTOR
    MRIScan(const MRIScan &copyScan);
    // DISTRUCTOR
    ~MRIScan(){}
    // ============
    // DATA MEMBERS
    // ============
    // Domain Dimension
    double domainSizeMin[3];
    double domainSizeMax[3];
    // Envelope Velocities
    double maxVelModule;
    // Velocities And Concentrations for all Measure Points
    int totalCellPoints;
    vector<MRICell> cellPoints;
    MRIIntVec mriCellTags;
    // Output Quantities
    vector<MRIOutput> outputs;
    // Utility Functions
    bool hasPressureGradient;
    bool hasRelativePressure;
    bool hasReynoldsStress;
    double scanTime;
    // MRI Expansion
    MRIExpansion* expansion;

    // ================
    // COMMON FUNCTIONS
    // ================

    // EXPORT DATA
    void ExportNodesToFile(std::string fileName);

    // DATA MANIPULATION
    void   ComputeQuantityGradient(int qtyID);
    void   UpdateVelocities();
    double GetDiffNorm(MRIScan* otherScan);

    // TOPOLOGY
    bool IsInnerCell(int Cell);

    // INFO
    virtual int getTotalFaces(){return 0;}
    double EvalAverageVelocityMod();
    virtual std::string WriteStatistics();

    // SCAN TRAVERSAL ALGORITHMS
    bool AreThereNotVisitedNeighbor(int cell, bool* visitedCell);
    bool AreThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor);
    int  findFirstNotVisite(int cellTotal, bool* visitedCell, std::vector<int> cellStack);
    void GetUnitVector(int CurrentCell, double* GlobalFaceCoords, double* &myVect);

    // SAVE QUANTITIES TO OUTPUTS
    void saveVelocity();

    // SAMPLING
    void SampleVelocities(MRISamplingOptions SamplingOptions);

    // PRESSURE
    void EvalRelativePressure(int startingCell, double refPressure);
    void PerformPressureIterations();
    void EvalCellPressureGradients(int currentCell, MRICellMaterial material,double* timeDeriv, double** firstDerivs, double** secondDerivs,double** ReynoldsStressGrad,double* pressureGrad);
    int  GetCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished, bool& secondStage);
    void EvalNoisyPressureGradientPoints();

    // TURBULENCE MODELLING
    void EvalReynoldsStressComponent(MRIThresholdCriteria* threshold);
    void EvalPressureIterative(int currentCell, double currentValue, bool* visitedCell, std::vector<int> &cellStack,int& cellCount);

    // FILTERING
    virtual void ApplyMedianFilter(int qtyID,int maxIt,int order,int filterType,MRIThresholdCriteria* threshold){throw MRIException("Error: Not Implemented!");}
    void ApplyThresholding(MRIThresholdCriteria* thresholdCriteria);
    void ApplySmoothingFilter();
    void ApplyGaussianNoise(double stDev);

    // SMP FILTER
    virtual void applySMPFilter(MRIOptions* options, bool isBC, MRICommunicator* comm){}
    virtual void RecoverGlobalErrorEstimates(double& AvNormError, double& AvAngleError){}

    // RECONSTRUCT FROM EXPANSION
    void ReconstructFromExpansion();

    // TENSORS AND VORTEX CRITERIA
    virtual void   EvalCellVelocityGradientDecomposition(int currentCell, double** deformation, double** rotation, double** firstDerivs){throw MRIException("Error: Not Implemented!");}
    virtual double EvalCellQCriterion(int currentCell, double** deformation, double** rotation){throw MRIException("Error: Not Implemented!");}
    virtual double EvalCellL2Criterion(int currentCell, double** deformation, double** rotation){throw MRIException("Error: Not Implemented!");}
    virtual double EvalCellDeltaCriterion(int currentCell, double** deformation, double** rotation, double** velGradient){throw MRIException("Error: Not Implemented!");}
    virtual double EvalCellVortexCriteria(int currentCell,int criteriaType, double** deformation, double** rotation, double** velGradient){throw MRIException("Error: Not Implemented!");}
    virtual void   EvalVortexCriteria(MRIThresholdCriteria* threshold){throw MRIException("Error: Not Implemented!");}
    virtual void   EvalVorticity(MRIThresholdCriteria* threshold){throw MRIException("Error: Not Implemented!");}
    virtual void   EvalEnstrophy(MRIThresholdCriteria* threshold){throw MRIException("Error: Not Implemented!");}

    // =================
    // VIRTUAL FUNCTIONS
    // =================

    // I/O
    virtual void ReadScanFromVOLFiles(std::string fileNameAn, std::string fileNameX, std::string fileNameY, std::string fileNameZ){throw MRIException("Error: Not Implemented!");}
    virtual void ExportToVOL(std::string FileName){throw MRIException("Error: Not Implemented!");}
    virtual void ExportToVTK(std::string fileName, MRIThresholdCriteria* threshold){throw MRIException("Error: Not Implemented!");}
    virtual void ExportToTECPLOT(std::string FileName, bool isFirstFile){throw MRIException("Error: Not Implemented!");}
    virtual void ExportForPoisson(string inputFileName,double density,double viscosity,MRIThresholdCriteria* threshold, const MRIDoubleMat& timeDerivs,
                                  bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                                  bool readMuTFromFile, string muTFile, double smagorinskyCoeff){throw MRIException("Error: ExportForPOISSONPartial Not Implemented!");}
    virtual void ExportForDistancing(string inputFileName, MRIThresholdCriteria* threshold){throw MRIException("Error: Not Implemented!");}
    virtual void WriteExpansionFile(std::string fileName){throw MRIException("Error: Not Implemented!");}

    // MODEL MANIPULATION
    virtual void Crop(double* limitBox){throw MRIException("Error: Not Implemented!");}
    virtual void ScaleVelocities(double factor){throw MRIException("Error: Not Implemented!");}
    virtual void ScalePositions(double factor){throw MRIException("Error: Not Implemented!");}

    // INFO
    virtual int EvalCentralCell(){throw MRIException("Error: Not Implemented!");}

    // TOPOLOGY
    virtual void GetCartesianNeighbourCells(int CurrentCell, std::vector<int> &coords, bool addself){throw MRIException("Error: Not Implemented!");}

    // COEFFICIENTS
    virtual void RebuildFromExpansion(MRIExpansion* expansion,bool useConstantFlux){throw MRIException("Error: Not Implemented!");}

    // EVAL DERIVATIVES
    virtual void EvalSpaceGradient(int currentCell,int qtyID, double* gradient){throw MRIException("Error: Not Implemented!");}
    virtual void EvalSpaceDerivs(int currentCell, MRIThresholdCriteria* threshold, double** firstDerivs, double** secondDerivs){throw MRIException("Error: Not Implemented!");}

    // EVAL REYNOLDS STRESSES
    virtual void EvalReynoldsStressGradient(int currentCell, double** ReynoldsStressGradient){throw MRIException("Error: Not Implemented!");}

    // COMPATIBILITY OF DIFFERENT SCANS
    bool isCompatibleWith(MRIScan* secondScan){throw MRIException("Error: Not Implemented!");}

    // SMP FILTER
    virtual int    EvalTotalVortex(){return 0;}
    virtual void   AssembleResidualVector(bool useBCFilter, MRIThresholdCriteria* thresholdCriteria, int &totalFaces, double* &ResVec, double* &filteredVec, double &resNorm){throw MRIException("Error: Not Implemented!");}
    virtual void   AssembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs){throw MRIException("Error: Not Implemented!");}
    virtual void   AssembleConstantPatternMPI(int currentDim, int &totalConstantFacesOnProc,
                                              std::vector<int> &facesIDOnProc, std::vector<double> &facesCoeffsOnProc,
                                              int minFaceOnProc, int maxFaceOnProc){throw MRIException("Error: Not Implemented!");}
    virtual void   AssembleStarShape(int vortexNumber, int &totalFaces,std::vector<int> &facesID,std::vector<double> &facesCoeffs){throw MRIException("Error: Not Implemented!");}
    virtual double EvalMaxDivergence(double* filteredVec){return 1000.0;}
    virtual void   RecoverCellVelocitiesRT0(bool useBCFilter, double* filteredVec){throw MRIException("Error: Not Implemented!");}

    // TENSOR AND VORTEX
    virtual void EvalSMPVortexCriteria(MRIExpansion* exp){throw MRIException("Error: Not Implemented!");}

    // MESSAGE PASSING
    virtual void DistributeScanData(MRICommunicator* comm);

    // CLEAN VELOCITIES ON BOUNDARY
    virtual void cleanNormalComponentOnBoundary(){throw MRIException("Error: Not Implemented!");}
    virtual void InterpolateBoundaryVelocities(){throw MRIException("Error: Not Implemented!");}
};

#endif // MRISCAN_H
