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
#include "mriOptions.h"

class MRIScan{
  public:
    // Cells Totals
    int cellTotals[3];
    double cellLength[3];
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
    double scanTime;
    // Member Functions
    // Constructor
    MRIScan(double currentTime);
    // Destructor
    ~MRIScan(){};
    // INFO FUNCTIONS
    std::string WriteStatistics();
    int GetTotalFaces();
  
    // READ FUNCTIONS
    void ReadPltFile(std::string PltFileName, bool DoReorderCells);
    void ReadScanFromVOLFiles(std::string fileNameAn, std::string fileNameX, std::string fileNameY, std::string fileNameZ);
	
    // WRITE FUNCTIONS
    void FillPLTHeader(bool hasPressureGradient, bool hasRelativePressure,std::vector<std::string> &pltHeader, bool isFirstFile);
    void ExportToLSDYNA(std::string LSFileName);
    void ExportToCSV(std::string FileName);
    void ExportToTECPLOT(std::string FileName, bool isFirstFile);
    void ExportNodesToFile(std::string FileName);
    void ExportToVOL(std::string FileName);
    void FlushToFile(std::string FileName);
    void ExportVelocitiesToFile(std::string fileName, bool append);
    void ExportToVTK(std::string fileName);

    // VOL DATA
    int  ReadBinVolFile(std::string FileName,MRIVolData & VolData);
    bool ValidateVOLBinData(MRIVolData &VolDataAn, MRIVolData &VolDataX, MRIVolData &VolDataY, MRIVolData &VolDataZ);
    void FormGlobadDataFromVOL(MRIVolData &VolDataAn, MRIVolData &VolDataX, MRIVolData &VolDataY, MRIVolData &VolDataZ);
    void CreateVolDataRecord(int volDataType, MRIVolData &VolData);
  
    // TRANSFORMATIONS AND THRESHOLD
    void ApplyThresholding(MRIThresholdCriteria thresholdCriteria);
    void Crop(double* limitBox);
    void ScaleVelocities(double factor);
    void ScalePositions(double factor);
	
    // Reorder Cells
    void ReorderCells(int* Perm);
    // Get Global Permutation
    void GetGlobalPermutation(int* &GlobalPerm);
    // MAPPING FUNCTIONS
    // Get Cell Number From Coords
    int GetCellNumber(MRIReal* coords);
    void GetNeighbourCells(int CurrentCell, int* &coords);
    // Get Face from Cell Vector
    int GetFacewithCellVector(int CurrentCell, double *UnitVector);
    // Get Adjacency Face
    int GetAdjacentFace(int GlobalNodeNumber, int AdjType);
    // Get Unit Vector From Current Cell To Face Centre
    void GetUnitVector(int CurrentCell, double* GlobalFaceCoords, double* &myVect);
    void GetGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, double* &globalCoords);
    int  FaceLocaltoGlobal(int LocalFace,int DimNumber,int SliceNumber);
    void MapIndexToCoords(int index, int* &intCoords);
    int  MapCoordsToIndex(int i, int j, int k);
    void GetLocalStarFaces(int StarNum, int CellsX, int CellsY, int &BottomFace, int &TopFace, int &LeftFace, int &RightFace);
    bool IsInnerCell(int Cell);
  
    // MP FILTER
    void   PerformPhysicsFiltering(MRIOptions Options, bool useBCFilter, bool useConstantPatterns, MRIThresholdCriteria thresholdCriteria);
    int    GetTotalBasisNumber();
    void   AssembleResidualVector(bool useBCFilter, MRIThresholdCriteria thresholdCriteria, int &totalFaces, double* &ResVec, double* &filteredVec, double &resNorm);
    void   AssembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs);
    void   AssembleStarShape(int dimNumber, int sliceNumber, int starNumber, int &totalFaces,std::vector<int> &facesID,std::vector<double> &facesCoeffs);
    void   ExpandStarShape(int totalStarFaces, int* facesID, double* facesCoeffs, double* &fullStarVector);
    void   UpdateVelocities();
    double EvalMaxDivergence(double* filteredVec);
    void   RecoverGlobalErrorEstimates(double& AvNormError,double& AvAngleError);
    void   RecoverCellVelocitiesRT0(bool useBCFilter, double* filteredVec);
    
    // FILTER MATRICES
    void AssembleEncodingMatrix(int &totalFaces, int &totalBasis, double** &Mat);
    void AssembleDecodingMatrix(int &totalFaces, int &totalBasis, double** &Mat);
    void AssembleStarMatrix(int &totalFaces, int &totalBasis, double** &Matrix);
      
    // INFO Functions
    double EvalAverageVelocityMod();
  
    // CELL SAMPLING
    void SampleVelocities(MRISamplingOptions SamplingOptions);
  
    // PRESSURE COMPUTATION
    // MAIN
    void EvalCellPressureGradients(int currentCell, MRICellMaterial material, double* timeDeriv, double** firstDerivs, double** secondDerivs, double* pressureGrad);
    void EvalRelativePressure(int startingCell, double refPressure);
    void PerformPressureIterations();
    // Others
    void EvalSpaceDerivs(int currentCell, double** firstDerivs, double** secondDerivs);
    void EvalPressureIterative(int currentCell, double currentValue, bool* visitedCell,int* otherCells, std::vector<int> &cellStack);
    bool AreThereNotVisitedNeighbor(int cell, bool* visitedCell);
    bool AreThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor);
    int  GetCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished);
    int  GetNextStartingCell(int currentCell, bool* visitedCell, bool* isBoundaryCell, bool &finished, int &bookmark);
    int  EvalCentralCell();
    
    // APPLY SMOOTHING FILTER - LAVISION
    void ApplySmoothingFilter();
  
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
    
    // ADD GAUSSIAN NOISE
    void ApplyGaussianNoise(double stDev);
    
    // TESTING FUNCTIONALITIES
    void TestScanAdjacency(std::string fileName);
    
   // STREAMLINES UTILITIES 
   void ComputeStreamlines(MRIStreamlineOptions &options, std::vector<MRIStreamline*> &streamlines);
   void PerformVelocityLinearInterpolation(double* coords, int CurrentCell, int xCell, int yCell, int zCell, double* &velocity);
   void GetPointVelocity(double xCoord, double yCoord, double zCoord, double* &pointVel);
   void EvalSingleStreamLine(double* start, MRIStreamlineOptions &options, MRIStreamline* &sL);
   void RecoverPlaneMaxMinCoords(MRIPlane plane, double* minCoords, double* maxCoords);
   void EvalSLTransverseDiffusionWithSpace(int totalSL, std::vector<MRIStreamline*> &streamlines, MRIDirection dir, double minCoord, double maxCoord, int totalSteps, std::vector<double> &space, std::vector<double> &crossDeviations);
   void EvalStreamLineStatistics(MRIDirection dir, MRIStreamlineOptions &options, std::vector<MRIStreamline*> &streamlines);
   void EvalSLArrivalPointDistribution(int totalSL, std::vector<MRIStreamline*> &streamlines, MRIDirection dir, double minCoord, double maxCoord, int totalSlices, std::vector<double> &sliceCenter, std::vector<double> &sliceNormArrivals);
};

#endif // MRISCAN_H
