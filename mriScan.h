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
#include "mriConstants.h"
#include "mriImagedata.h"
#include "mriExpansion.h"
#include "mriException.h"
#include "mriOptions.h"

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
    std::vector<MRICell> cellPoints;
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

    // SCAN TRAVERSAL ALGORITHMS
    bool AreThereNotVisitedNeighbor(int cell, bool* visitedCell);
    bool AreThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor);
    int  findFirstNotVisite(int cellTotal, bool* visitedCell, std::vector<int> cellStack);
    void GetUnitVector(int CurrentCell, double* GlobalFaceCoords, double* &myVect);

    // SAMPLING
    void SampleVelocities(MRISamplingOptions SamplingOptions);

    // PRESSURE
    void EvalRelativePressure(int startingCell, double refPressure);
    void PerformPressureIterations();
    void EvalCellPressureGradients(int currentCell, MRICellMaterial material,double* timeDeriv, double** firstDerivs, double** secondDerivs,double** ReynoldsStressGrad,double* pressureGrad);
    int  GetCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished, bool& secondStage);
    void EvalNoisyPressureGradientPoints();

    // TURBULENCE MODELLING
    void EvalReynoldsStressComponent();
    void EvalPressureIterative(int currentCell, double currentValue, bool* visitedCell, std::vector<int> &cellStack,int& cellCount);

    // FILTERING
    void ApplyMedianFilter(int qtyID,int maxIt);
    void ApplyThresholding(MRIThresholdCriteria thresholdCriteria);
    void ApplySmoothingFilter();
    void ApplyGaussianNoise(double stDev);

    // SMP FILTER
    void applySMPFilter(MRIOptions* options);
    void RecoverGlobalErrorEstimates(double& AvNormError, double& AvAngleError);

    // RECONSTRUCT FROM EXPANSION
    void ReconstructFromExpansion();

    // TENSORS AND VORTEX CRITERIA
    void   EvalCellVelocityGradientDecomposition(int currentCell, double** deformation, double** rotation, double** firstDerivs);
    double EvalCellQCriterion(int currentCell, double** deformation, double** rotation);
    double EvalCellL2Criterion(int currentCell, double** deformation, double** rotation);
    double EvalCellDeltaCriterion(int currentCell, double** deformation, double** rotation, double** velGradient);
    double EvalCellVortexCriteria(int currentCell,int criteriaType, double** deformation, double** rotation, double** velGradient);
    void   EvalVortexCriteria();
    void   EvalVorticity();
    void   EvalEnstrophy();

    // =================
    // VIRTUAL FUNCTIONS
    // =================

    // I/O
    virtual void ReadScanFromVOLFiles(std::string fileNameAn, std::string fileNameX, std::string fileNameY, std::string fileNameZ){throw MRIException("Error: Not Implemented!");}
    virtual void ExportToVOL(std::string FileName){throw MRIException("Error: Not Implemented!");}
    virtual void ExportToVTK(std::string fileName){throw MRIException("Error: Not Implemented!");}
    virtual void ExportToTECPLOT(std::string FileName, bool isFirstFile){throw MRIException("Error: Not Implemented!");}
    virtual void WriteExpansionFile(std::string fileName){throw MRIException("Error: Not Implemented!");}

    // MODEL MANIPULATION
    virtual void ApplyThresholding(MRIThresholdCriteria* thresholdCriteria){throw MRIException("Error: Not Implemented!");}
    virtual void Crop(double* limitBox){throw MRIException("Error: Not Implemented!");}
    virtual void ScaleVelocities(double factor){throw MRIException("Error: Not Implemented!");}
    virtual void ScalePositions(double factor){throw MRIException("Error: Not Implemented!");}

    // INFO
    virtual int EvalCentralCell(){throw MRIException("Error: Not Implemented!");}

    // TOPOLOGY
    virtual void GetNeighbourCells(int CurrentCell, std::vector<int> &coords){throw MRIException("Error: Not Implemented!");}

    // COEFFICIENTS
    virtual void RebuildFromExpansion(MRIExpansion* expansion,bool useConstantFlux){throw MRIException("Error: Not Implemented!");}

    // EVAL DERIVATIVES
    virtual void EvalSpaceGradient(int currentCell,int qtyID, double* gradient){throw MRIException("Error: Not Implemented!");}
    virtual void EvalSpaceDerivs(int currentCell, double** firstDerivs, double** secondDerivs){throw MRIException("Error: Not Implemented!");}

    // EVAL REYNOLDS STRESSES
    virtual void EvalReynoldsStressGradient(int currentCell, double** ReynoldsStressGradient){throw MRIException("Error: Not Implemented!");}

    // COMPATIBILITY OF DIFFERENT SCANS
    bool isCompatibleWith(MRIScan* secondScan){throw MRIException("Error: Not Implemented!");}

    // SMP FILTER
    virtual int    EvalTotalVortex(){return 0;}
    virtual void   AssembleResidualVector(bool useBCFilter, MRIThresholdCriteria* thresholdCriteria, int &totalFaces, double* &ResVec, double* &filteredVec, double &resNorm){throw MRIException("Error: Not Implemented!");}
    virtual void   AssembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs){throw MRIException("Error: Not Implemented!");}
    virtual void   AssembleStarShape(int vortexNumber, int &totalFaces,std::vector<int> &facesID,std::vector<double> &facesCoeffs){throw MRIException("Error: Not Implemented!");}
    virtual double EvalMaxDivergence(double* filteredVec){return 1000.0;}
    virtual void   RecoverCellVelocitiesRT0(bool useBCFilter, double* filteredVec){throw MRIException("Error: Not Implemented!");}

    // TENSOR AND VORTEX
    void EvalSMPVortexCriteria(MRIExpansion* exp){throw MRIException("Error: Not Implemented!");}
};

#endif // MRISCAN_H
