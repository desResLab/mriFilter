#ifndef MRIUNSTRUCTUREDSCAN_H
#define MRIUNSTRUCTUREDSCAN_H

#include <stdio.h>
#include <string>
#include <vector>

#include "mriScan.h"
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
#include "mriOptions.h"

// ========================
// UNSTRUCTURED GRID LAYOUT
// ========================
class MRIStructuredScan: public MRIScan{
  public:
    // Cells Totals
    int cellTotals[3];
    std::vector<std::vector<double>> cellLengths;
    // Aux Vector to Assemble Vortex Faces
    std::vector<std::vector<int>> cellConnections;
    std::vector<std::vector<int>> cellFaces;
    std::vector<std::vector<int>> faceConnections;
    std::vector<std::vector<int>> egdeConnections;
    // ================
    // MEMBER FUNCTIONS
    // ================
    // Constructor
    // 1 - Empty
    MRIStructuredScan(double currentTime):MRIScan(currentTime){}
    // 2 - Create a Zero Scan from another Scan
    MRIStructuredScan(MRIStructuredScan &copyScan);
    // Destructor
    virtual ~MRIStructuredScan(){}
    // INFO FUNCTIONS
    std::string WriteStatistics();
    int GetTotalFaces();
  
    // ==============
    // READ FUNCTIONS
    // ==============
    void ReadPltFile(std::string PltFileName, bool DoReorderCells);
    virtual void ReadScanFromVOLFiles(std::string fileNameAn, std::string fileNameX, std::string fileNameY, std::string fileNameZ);
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
    virtual void ExportToTECPLOT(std::string FileName, bool isFirstFile);
    void ExportNodesToFile(std::string FileName);
    virtual void ExportToVOL(std::string FileName);
    void FlushToFile(std::string FileName);
    void ExportVelocitiesToFile(std::string fileName, bool append);
    virtual void ExportToVTK(std::string fileName);
    virtual void WriteExpansionFile(std::string fileName);

    // ========
    // VOL DATA
    // ========
    int  ReadBinVolFile(std::string FileName,MRIVolData & VolData);
    bool ValidateVOLBinData(MRIVolData &VolDataAn, MRIVolData &VolDataX, MRIVolData &VolDataY, MRIVolData &VolDataZ);
    void FormGlobadDataFromVOL(MRIVolData &VolDataAn, MRIVolData &VolDataX, MRIVolData &VolDataY, MRIVolData &VolDataZ);
    void CreateVolDataRecord(int volDataType, MRIVolData &VolData);
  
    // =============================
    // TRANSFORMATIONS AND THRESHOLD
    // =============================
    virtual void Crop(double* limitBox);
    virtual void ScaleVelocities(double factor);
    virtual void ScalePositions(double factor);

    // ==========================================
    // RECONSTRUCTION FROM EXPANSION COEFFICIENTS
    // ==========================================
    virtual void RebuildFromExpansion(MRIExpansion* expansion,bool useConstantFlux);

    // ===============================
    // RECONSTRUCTION FROM FACE FLUXES
    // ==============================
    void RebuildFromFaceFluxes(double* faceFluxes);

    // Reorder Cells    
    void ReorderCells(int* Perm);
    // Get Global Permutation
    void GetGlobalPermutation(int* &GlobalPerm);
    // REORDER GLOBAL SCAN
    void ReorderScan();
    // MAPPING FUNCTIONS
    // Get Cell Number From Coords
    int  GetCellNumber(MRIReal* coords);
    void GetNeighbourCells(int CurrentCell, std::vector<int> &coords);
    bool isCompatibleWith(MRIStructuredScan* secondScan);
    // Get Face from Cell Vector
    int  GetFacewithCellVector(int CurrentCell, double *UnitVector);
    // Get Adjacency Face
    int  GetAdjacentFace(int GlobalNodeNumber, int AdjType);
    // Get Unit Vector From Current Cell To Face Centre
    void GetGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, double* &globalCoords);
    int  FaceLocaltoGlobal(int LocalFace,int DimNumber,int SliceNumber);
    // Sequential Index to Integer Coords
    void MapIndexToCoords(int index, int* intCoords);
    int  MapCoordsToIndex(int i, int j, int k);
    // Map Integer Coords to Position
    void MapCoordsToPosition(int* coords,double* pos);

    void GetLocalStarFaces(int StarNum, int CellsX, int CellsY, int &BottomFace, int &TopFace, int &LeftFace, int &RightFace);
    int  findFirstNotVisited(int cellTotal, bool* visitedCell, std::vector<int> cellStack);
    void formNotVisitedList(int cellTotal, bool* visitedCell,std::vector<bool>& notVisitedList);
  
    // =========
    // MP FILTER
    // =========
    int    GetTotalBasisNumber();
    virtual void   AssembleResidualVector(bool useBCFilter, MRIThresholdCriteria thresholdCriteria, int &totalFaces, double* &ResVec, double* &filteredVec, double &resNorm);
    virtual void   AssembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs);
    virtual void   AssembleStarShape(int vortexNumber, int &totalFaces,std::vector<int> &facesID,std::vector<double> &facesCoeffs);
    virtual double EvalMaxDivergence(double* filteredVec);
    void   ExpandStarShape(int totalStarFaces, int* facesID, double* facesCoeffs, double* &fullStarVector);
    void   RecoverGlobalErrorEstimates(double& AvNormError,double& AvAngleError);
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
    bool   isUniform();
    virtual int EvalTotalVortex();
    virtual int getTotalFaces();

    // CELL SAMPLING
    void SampleVelocities(MRISamplingOptions SamplingOptions);

    // GRADIENTS AND DERIVATIVES
    virtual void EvalSpaceDerivs(int currentCell, double** firstDerivs, double** secondDerivs);
    virtual void EvalSpaceGradient(int currentCell,int qtyID, double* gradient);
    void ComputeQuantityGradient(int qtyID);
  
    // PRESSURE COMPUTATION
    void EvalRelativePressure(int startingCell, double refPressure);
    void PerformPressureIterations();
    // Others
    void EvalPressureIterative(int currentCell, double currentValue, bool* visitedCell,int* otherCells, std::vector<int> &cellStack,int& cellCount);
    bool AreThereNotVisitedNeighbor(int cell, bool* visitedCell);
    bool AreThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor);
    int  GetCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished, bool& secondStage);
    int  GetNextStartingCell(int currentCell, bool* visitedCell, bool* isBoundaryCell, bool &finished, int &bookmark);
    virtual int EvalCentralCell();

    // REYNOLDS STRESS COMPUTATION    
    void EvalReynoldsStressGradient(int currentCell, double** ReynoldsStressGradient);
    
    // APPLY SMOOTHING FILTER - LAVISION
    void ApplySmoothingFilter();

    // APPLY MEDIAN FILTER
    void ApplyMedianFilter(int qtyID,int maxIt);
    void ThresholdQuantity(int qtyID,double threshold);
    void EvalNoisyPressureGradientPoints();
  
    // SAMPLE FLOWS
    void CreateSampleCase(MRISamples sample, int sizeX, int sizeY, int sizeZ, double distX, double distY, double distZ, double currTime, MRIDirection dir);    
    void AssignVelocitySignature(MRIDirection dir, MRISamples sample, double currTime);
    void AssignConstantSignature(MRIDirection dir);
    void AssignStagnationFlowSignature(MRIDirection dir);
    void AssignPoiseilleSignature(MRIDirection dir);
    void AssignCylindricalFlowSignature(MRIDirection dir);
    void AssignSphericalFlowSignature(MRIDirection dir);
    void AssignToroidalVortexFlowSignature();
    void AssignTimeDependentPoiseilleSignature(double omega, double radius, double viscosity, double currtime, double maxVel);
    void AssignConstantFlowWithStep();
    void AssignRandomStandardGaussianFlow();
    void AssignRandomComponent(const int kdirX,stdRndGenerator &generator);
    void AssignZeroVelocities();

    // VORTEX IDENTIFICATION
    // With SMP Expansion Coefficients
    virtual void EvalSMPVortexCriteria(MRIExpansion* exp);
    // Used Q,L2,Delta Criteria
    void   EvalVortexCriteria();
    void   EvalCellVelocityGradientDecomposition(int currentCell, double** deformation, double** rotation, double** firstDerivs);
    double EvalCellQCriterion(int currentCell, double** deformation, double** rotation);
    double EvalCellL2Criterion(int currentCell, double** deformation, double** rotation);
    double EvalCellDeltaCriterion(int currentCell, double** deformation, double** rotation, double** velGradient);
    double EvalCellVortexCriteria(int currentCell,int criteriaType, double** deformation, double** rotation, double** velGradient);
    void   EvalVorticity();
    void   EvalEnstrophy();

    // SPATIAL REPRESENTATION OF VORTEX COEFFICIENTS
    double EvalVortexCriteria(MRIExpansion* exp);
    void getNeighborVortexes(int cellNumber,int dim,int& idx1,int& idx2,int& idx3,int& idx4);
    
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
   double GetDiffNorm(MRIStructuredScan* otherScan);
};

#endif // MRIUNSTRUCTUREDSCAN_H