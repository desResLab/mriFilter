#ifndef MRISCAN_H
#define MRISCAN_H

#include <math.h>
#include <stdio.h>
#include <string>
#include <limits>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "mriCell.h"
#include "mriTypes.h"
#include "mriExpansion.h"
#include "mriThresholdCriteria.h"
#include "mriUtils.h"
#include "mriImagedata.h"
#include "mriConstants.h"
#include "mriVolData.h"
#include "mriException.h"
#include "schMessages.h"
#include "mriCellMaterial.h"
#include "mriSamplingOptions.h"
#include "mriStreamline.h"
#include "mriStreamlineOptions.h"
#include "mriOutput.h"

#include "mriOptions.h"
#include "mriCommunicator.h"

class MRIScan;
class MRICommunicator;

// ========================
// UNSTRUCTURED GRID LAYOUT
// ========================
class MRIScan{
  public:
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
    // MEMBER FUNCTIONS
    // ================
    // Constructor
    // 1 - Empty
    MRIScan(double currentTime);
    // 2 - Create a Zero Scan from another Scan
    MRIScan(const MRIScan& copyScan);
    // Destructor
    ~MRIScan(){}
    // INFO FUNCTIONS
    std::string WriteStatistics();
    int GetTotalFaces();

    void CreateGridFromVTKStructuredPoints(vtkStructuredPointsOptionRecord opts);
  
    // ==============
    // READ FUNCTIONS
    // ==============
    void ReadVTKStructuredPoints(std::string vtkFileName, bool DoReorderCells);
    void ReadPltFile(std::string PltFileName, bool DoReorderCells);
    void ReadScanFromVOLFiles(std::string fileNameAn, std::string fileNameX, std::string fileNameY, std::string fileNameZ);
    void ReadScanFromSingleVOLFile(std::string fileName);
    void ReadFromExpansionFile(std::string fileName,bool applyThreshold,int thresholdType,double thresholdValue);

    // ==================
    // READ FROM RAW DATA
    // ==================
    void ReadRAWFileSequence(std::string fileListName);
    int  ReadRawImage(std::string FileName, MRIImageData &data);
	
    // ===============
    // WRITE FUNCTIONS
    // ===============
    void FillPLTHeader(std::vector<std::string> &pltHeader, bool isFirstFile);
    void ExportToLSDYNA(std::string LSFileName);
    void ExportToCSV(std::string FileName);
    void ExportNodesToFile(std::string FileName);
    void FlushToFile(std::string FileName);
    void ExportVelocitiesToFile(std::string fileName, bool append);
    // VIRTUAL
    void ExportToVOL(std::string FileName);
    void ExportToTECPLOT(std::string FileName, bool isFirstFile);
    void ExportToVTK(std::string fileName, MRIThresholdCriteria* threshold);
    void WriteExpansionFile(std::string fileName);
    // Export to Poisson Solver Only element with significant concentration
    void ExportForDistancing(string inputFileName, MRIThresholdCriteria* threshold);
    void ExportForPoisson(string inputFileName,double density,double viscosity,MRIThresholdCriteria* threshold, const MRIDoubleMat& timeDerivs,
                                  bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                                  bool readMuTFromFile, string muTFile, double smagorinskyCoeff);

    // ========
    // VOL DATA
    // ========
    int  ReadBinVolFile(std::string FileName,MRIVolData & VolData);
    bool ValidateVOLBinData(MRIVolData &VolDataAn, MRIVolData &VolDataX, MRIVolData &VolDataY, MRIVolData &VolDataZ);
    void FormGlobadDataFromVOL(MRIVolData &VolDataAn, MRIVolData &VolDataX, MRIVolData &VolDataY, MRIVolData &VolDataZ);
    void CreateVolDataRecord(int volDataType, MRIVolData &VolData);

    // ========
    // TOPOLOGY
    // ========
    // VORTEX COEFFICIENTS
    double getEdgeFaceVortexCoeff(int edgeID, int faceID);
    // AUXILIARY GEOMETRY
    void   getEdgeDirection(int edgeID, double* edgeDirVector);
    void   getEdgeCenter(int edgeID, double* ec);
    void   getFaceCenter(int faceID, double* fc);
    void   getEdgeToFaceDirection(int edgeID, int faceID, double* edgeFaceVector);
    // VOLUME
    double evalCellVolume(int cellNumber);
    // OTHER
    int GetCellFaceID(int CellId,int FaceId);
    bool hasUniformSpacing();

    // =============================
    // TRANSFORMATIONS AND THRESHOLD
    // =============================
    void Crop(double* limitBox);
    void ScaleVelocities(double factor);
    void ScalePositions(double factor);

    // ==========================================
    // RECONSTRUCTION FROM EXPANSION COEFFICIENTS
    // ==========================================
    void RebuildFromExpansion(MRIExpansion* expansion,bool useConstantFlux);

    // ===============================
    // RECONSTRUCTION FROM FACE FLUXES
    // ==============================
    void RebuildFromFaceFluxes(double* faceFluxes);

    // Reorder Cells    
    void ReorderCells(std::vector<int> Perm);
    // Get Global Permutation
    void GetGlobalPermutation(std::vector<int> &GlobalPerm);
    // REORDER GLOBAL SCAN
    void ReorderScan();

    // MAPPING FUNCTIONS
    // Get Cell Number From Coords
    int  GetCellNumber(double* coords);
    void GetCartesianNeighbourCells(int CurrentCell, std::vector<int> &coords, bool addself);
    void GetStructuredNeighbourCells(int centreCell,int order,MRIThresholdCriteria* threshold,MRIIntVec& cellNeighbors);
    bool isCompatibleWith(MRIScan* secondScan);
    // Get Face from Cell Vector
    int  GetFacewithCellVector(int CurrentCell, double *UnitVector);
    // Get Adjacency Face
    int  GetAdjacentFace(int GlobalNodeNumber, int AdjType);
    // Get Unit Vector From Current Cell To Face Centre
    void GetGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, double* &globalCoords);
    int  FaceLocaltoGlobal(int LocalFace,int DimNumber,int SliceNumber);

    // Sequential Index to Integer Coords
    void MapIndexToCoords(int index, int* intCoords);
    void MapIndexToAuxNodeCoords(int index, int* intCoords);
    int  MapCoordsToIndex(int i, int j, int k);

    // Map Integer Coords to Position
    void MapCoordsToPosition(int* coords, bool addMeshMinima, double* pos);
    void MapAuxCoordsToPosition(int* auxCoords, double* pos);
    void GetLocalStarFaces(int StarNum, int CellsX, int CellsY, int &BottomFace, int &TopFace, int &LeftFace, int &RightFace);
    int  findFirstNotVisited(int cellTotal, bool* visitedCell, std::vector<int> cellStack);
    void formNotVisitedList(int cellTotal, bool* visitedCell,std::vector<bool>& notVisitedList);

    // Map cell vector to face vector
    void cellToFace(bool deleteWalls, MRIThresholdCriteria* thresholdCriteria,MRIDoubleMat cellVec, MRIDoubleVec &faceVec);
    void cellToFacePartial(MRIIntVec elUsageMap, MRIThresholdCriteria* thresholdCriteria,
                           MRIDoubleMat cellVec, MRIDoubleVec &faceVec);

    // =========
    // MP FILTER
    // =========
    int    GetTotalBasisNumber();
    void   applySMPFilter(MRIOptions* options, bool isBC, MRICommunicator* comm);
    void   AssembleResidualVector(bool useBCFilter, MRIThresholdCriteria* thresholdCriteria, int &totalFaces, double* &ResVec, double* &filteredVec, double &resNorm);
    void   AssembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs);
    void   AssembleConstantPatternMPI(int currentDim, int &totalConstantFacesOnProc,
                                              std::vector<int> &facesIDOnProc, std::vector<double> &facesCoeffsOnProc,
                                              int minFaceOnProc, int maxFaceOnProc,MRICommunicator* comm);
    void   AssembleStarShape(int vortexNumber, int &totalFaces,std::vector<int> &facesID,std::vector<double> &facesCoeffs);
    double EvalMaxDivergence(double* filteredVec);
    void   RecoverGlobalErrorEstimates(double& AvNormError, double& AvAngleError);
    void   ExpandStarShape(int totalStarFaces, int* facesID, double* facesCoeffs, double* &fullStarVector);
    void   RecoverCellVelocitiesRT0(bool useBCFilter, double* filteredVec);
    void   ReconstructFromExpansion();
    void   getDimensionSliceStarFromVortex(int vortexNumber,int &dimNumber,int &sliceNumber,int &starNumber);
    
    // FILTER MATRICES
    void AssembleEncodingMatrix(int &totalFaces, int &totalBasis, double** &Mat);
    void AssembleDecodingMatrix(int &totalFaces, int &totalBasis, double** &Mat);
    void AssembleStarMatrix(int &totalFaces, int &totalBasis, double** &Matrix);
      
    // ==============
    // INFO FUNCTIONS
    // ==============
    double EvalAverageVelocityMod();
    void   evalCellAreas(int cellNumber,double* Areas);
    int    getTotalAuxNodes();
    int    EvalTotalVortex();
    int    getTotalFaces();

    // CELL SAMPLING
    void SampleVelocities(MRISamplingOptions SamplingOptions);

    // GRADIENTS AND DERIVATIVES
    void EvalSpaceDerivs(int currentCell, MRIThresholdCriteria* threshold, double** firstDerivs, double** secondDerivs);
    void EvalSpaceGradient(int currentCell,int qtyID, double* gradient);
    void ComputeQuantityGradient(int qtyID);

    // DIVERGENCE
    MRIDoubleVec evalCellDivergences(MRIDoubleVec faceVec);
  
    // PRESSURE COMPUTATION
    void EvalRelativePressure(int startingCell, double refPressure);
    void PerformPressureIterations();
    // Others
    void EvalPressureIterative(int currentCell, double currentValue, bool* visitedCell,int* otherCells, std::vector<int> &cellStack,int& cellCount);
    bool AreThereNotVisitedNeighbor(int cell, bool* visitedCell);
    bool AreThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor);
    int  GetCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished, bool& secondStage);
    int  GetNextStartingCell(int currentCell, bool* visitedCell, bool* isBoundaryCell, bool &finished, int &bookmark);
    int EvalCentralCell();
    int SolvePoissonEquation(MRICommunicator* comm);

    // REYNOLDS STRESS COMPUTATION    
    void EvalReynoldsStressGradient(int currentCell, double** ReynoldsStressGradient);
    void evalPradtlTurbViscosity(MRIDoubleMat cellDistance, MRIThresholdCriteria* threshold, double density, MRIDoubleMat& turbViscosity);
    void evalSmagorinskyLillyTurbViscosity(double density, double smagorinskyCoeff, MRIThresholdCriteria* threshold, MRIDoubleMat& turbViscosity);
    
    // APPLY SMOOTHING FILTER - LAVISION
    void ApplySmoothingFilter();
    void ApplyMedianFilter(int qtyID,int maxIt,int order,int filterType,MRIThresholdCriteria* threshold);

    // THRESHOLD
    void ThresholdQuantity(int qtyID,double threshold);
    void EvalNoisyPressureGradientPoints();
  
    // SAMPLE FLOWS
    void CreateSampleCase(MRISamples sampleType, vector<double> params);
    void AssignVelocitySignature(MRIDirection dir, MRISamples sample, double currTime);
    void AssignConstantSignature(MRIDirection dir);
    void AssignStagnationFlowSignature(MRIDirection dir);
    void AssignPoiseilleSignature(MRIDirection dir);
    void AssignCylindricalFlowSignature(MRIDirection dir);
    void AssignSphericalFlowSignature(MRIDirection dir);
    void AssignToroidalVortexFlowSignature();
    void AssignTimeDependentPoiseilleSignature(double omega, double radius, double viscosity, double currtime, double maxVel);
    void AssignConstantFlowWithStep();
    void AssignTaylorVortexSignature(MRIDirection dir);
    void AssignRandomStandardGaussianFlow();
    void AssignRandomComponent(const int kdirX,stdRndGenerator &generator);
    void AssignZeroVelocities();

    // VORTEX IDENTIFICATION
    void   EvalCellVelocityGradientDecomposition(int currentCell, double** deformation, double** rotation, double** firstDerivs);
    double EvalCellQCriterion(int currentCell, double** deformation, double** rotation);
    double EvalCellL2Criterion(int currentCell, double** deformation, double** rotation);
    double EvalCellDeltaCriterion(int currentCell, double** deformation, double** rotation, double** velGradient);
    double EvalCellVortexCriteria(int currentCell,int criteriaType, double** deformation, double** rotation, double** velGradient);
    void   EvalVortexCriteria(MRIThresholdCriteria* threshold);
    void   EvalVorticity(MRIThresholdCriteria* threshold);
    void   EvalEnstrophy(MRIThresholdCriteria* threshold);
    void   EvalSMPVortexCriteria(MRIExpansion* exp);

    // SPATIAL REPRESENTATION OF VORTEX COEFFICIENTS
    double EvalVortexCriteria(MRIExpansion* exp);
    void   getNeighborVortexes(int cellNumber,int dim, MRIIntVec& idx);
    
    // ADD GAUSSIAN NOISE
    void ApplyGaussianNoise(double stDev);
    
    // TESTING FUNCTIONALITIES
    void TestScanAdjacency(std::string fileName);
    
   // STREAMLINES UTILITIES 
   void ComputeStreamlines(std::string outName, MRIStreamlineOptions &options, std::vector<MRIStreamline*> &streamlines);
   void PerformVelocityLinearInterpolation(double* coords, int CurrentCell, int xCell, int yCell, int zCell, double* &velocity);
   void GetPointVelocity(double xCoord, double yCoord, double zCoord, double* &pointVel);
   void EvalSingleStreamLine(double* start, MRIStreamlineOptions &options, MRIStreamline* &sL);
   void EvalSLTransverseDiffusionWithSpace(int totalSL, std::vector<MRIStreamline*> &streamlines, MRIDirection dir, double minCoord, double maxCoord, int totalSteps, std::vector<double> &space, std::vector<double> &crossDeviations);
   void EvalStreamLineStatistics(std::string outName, MRIDirection dir, MRIStreamlineOptions &options, std::vector<MRIStreamline*> &streamlines);
   void EvalSLArrivalPointDistribution(int totalSL, std::vector<MRIStreamline*> &streamlines, MRIDirection dir, double minCoord, double maxCoord, int totalSlices, std::vector<double> &sliceCenter, std::vector<double> &sliceNormArrivals);

   // COMPARISON BETWEEN SCANS
   double GetDiffNorm(MRIScan* otherScan);

   void buildMetisConnectivities(int *eptr,int *eind);

   // MESSAGE PASSING
   void formVortexList(int totVortex,int* minFace,int* maxFace,MRIIntVec& innerVortexList,MRIIntVec& boundaryVortexList,MRICommunicator* comm);
   void passScanData(MRICommunicator* comm);
   void DistributeScanData(MRICommunicator* comm);

   // BOUNDARY CLEANING
   void cleanNormalComponentOnBoundary(MRIThresholdCriteria* threshold);
   void InterpolateBoundaryVelocities(MRIThresholdCriteria* threshold);
   void projectCellVelocity(int cell,double* normal);
   int getOppositeCell(int cell, double* normal);
   void tagByNeighbour(int tag,int* cellTags, bool* isTaggable,int startingCell);
   bool hasUntaggedNeighbours(int cell,int* cellTags, bool* isTaggable);
   void setWallFluxesToZero(bool* isFaceOnWalls, MRIDoubleVec& poissonSourceFaceVec);
};

#endif // MRISTRUCTUREDSCAN_H
