#ifndef MRISTRUCTUREDSCAN_H
#define MRISTRUCTUREDSCAN_H

#include <stdio.h>
#include <string>
#include <vector>

#include "mriScan.h"
#include "mriTypes.h"
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
#include "mriCommunicator.h"

// RELATIVE POSITION BETWEEN EDGE AND FACE
enum EdgeFacePositionType{ptTop,ptBottom,ptLeft,ptRight};

// AUXILIARY FACE
struct mriFace{
  int number;
  std::vector<int> connections;
};

// AUXILIARY EDGE
struct mriEdge{
  int number;
  std::vector<int> connections;
};

// ===================
// TYPES FOR PLT FILES
// ===================
enum pltFileTypes{
  pltUNIFORM,
  pltSTRUCTURED
};

struct PLTOptionRecord{
  int i;
  int j;
  int k;
  int N;
  int E;
  pltFileTypes type;
};

// ===================
// TYPES FOR PLT FILES
// ===================
struct vtkStructuredPointsOptionRecord{
  bool isASCII;
  bool isValidDataset;
  int dimensions[3];
  double origin[3];
  double spacing[3];
  int numDefined;
  bool isDefined[5];
  MRIIntVec dataBlockStart;
  MRIIntVec dataBlockType;
  MRIBoolVec dataBlockRead;
};

class MRIScan;
class MRICommunicator;

// ========================
// UNSTRUCTURED GRID LAYOUT
// ========================
class MRIStructuredScan: public MRIScan{
  protected:
    void CreateGridFromVTKStructuredPoints(vtkStructuredPointsOptionRecord opts);
  public:
    // Cells Totals
    MRIIntVec cellTotals;
    MRIDoubleMat cellLengths;
    // Cells Topology
    MRIIntMat cellConnections;
    MRIIntMat cellFaces;
    // Face Topology
    MRIIntMat faceCells;
    MRIIntMat faceConnections;
    MRIIntMat faceEdges;
    MRIDoubleVec faceArea;
    MRIDoubleMat faceNormal;
    // Edge Topology
    MRIIntMat edgeConnections;
    MRIIntMat edgeFaces;

    // ================
    // MEMBER FUNCTIONS
    // ================
    // Constructor
    // 1 - Empty
    MRIStructuredScan(double currentTime);
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
    void ReadVTKStructuredPoints(std::string vtkFileName, bool DoReorderCells);
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
    void ExportNodesToFile(std::string FileName);
    void FlushToFile(std::string FileName);
    void ExportVelocitiesToFile(std::string fileName, bool append);
    // VIRTUAL
    virtual void ExportToVOL(std::string FileName);
    virtual void ExportToTECPLOT(std::string FileName, bool isFirstFile);
    virtual void ExportToVTK(std::string fileName);
    virtual void WriteExpansionFile(std::string fileName);
    void         ExportForPOISSON(string inputFileName);

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
    // MAIN
    void CreateTopology();
    // CELLS
    void buildCellConnections();
    int addToFaceConnections(std::vector<std::vector<mriFace* >> &AuxFirstNodeFaceList, std::vector<int> faceIds);
    // FACES
    void buildFaceConnections();
    void buildFaceCells();
    void buildFaceAreasAndNormals();
    void getExternalFaceNormal(int cellID, int localFaceID, double* extNormal);
    // EDGES
    void buildEdgeConnections();
    int  addToEdgeConnections(std::vector<std::vector<mriEdge*>> &AuxFirstNodeEdgeList, std::vector<int> edgeIds);
    void assembleEdgeToFaceConnectivity(std::vector<int> &vortexBottomFaces,std::vector<int> &vortexTopFaces,std::vector<int> &vortexLeftFaces,std::vector<int> &vortexRightFaces);
    // VORTEX COEFFICIENTS
    double getEdgeFaceVortexCoeff(int edgeID, int faceID);
    // AUXILIARY GEOMETRY
    void   getAuxNodeCoordinates(int nodeNum, double* pos);
    void   getEdgeDirection(int edgeID, std::vector<double> &edgeDirVector);
    void   getEdgeCenter(int edgeID, double* ec);
    void   getFaceCenter(int faceID, double* fc);
    void   getEdgeToFaceDirection(int edgeID, int faceID, std::vector<double> &edgeFaceVector);
    // VOLUME
    double evalCellVolume(int cellNumber);
    // OTHER
    int GetCellFaceID(int CellId,int FaceId);

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
    void ReorderCells(std::vector<int> Perm);
    // Get Global Permutation
    void GetGlobalPermutation(std::vector<int> &GlobalPerm);
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

    // =========
    // MP FILTER
    // =========
    int    GetTotalBasisNumber();
    virtual void   applySMPFilter(MRIOptions* options, bool isBC, MRICommunicator* comm);
    virtual void   AssembleResidualVector(bool useBCFilter, MRIThresholdCriteria* thresholdCriteria, int &totalFaces, double* &ResVec, double* &filteredVec, double &resNorm);
    virtual void   AssembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs);
    virtual void   AssembleConstantPatternMPI(int currentDim, int &totalConstantFacesOnProc,
                                              std::vector<int> &facesIDOnProc, std::vector<double> &facesCoeffsOnProc,
                                              int minFaceOnProc, int maxFaceOnProc,MRICommunicator* comm);
    virtual void   AssembleStarShape(int vortexNumber, int &totalFaces,std::vector<int> &facesID,std::vector<double> &facesCoeffs);
    virtual double EvalMaxDivergence(double* filteredVec);
    virtual void   RecoverGlobalErrorEstimates(double& AvNormError, double& AvAngleError);
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
    bool   isUniform();
    virtual int  EvalTotalVortex();
    virtual int  getTotalFaces();

    // CELL SAMPLING
    void SampleVelocities(MRISamplingOptions SamplingOptions);

    // GRADIENTS AND DERIVATIVES
    virtual void EvalSpaceDerivs(int currentCell, double** firstDerivs, double** secondDerivs);
    virtual void EvalSpaceGradient(int currentCell,int qtyID, double* gradient);
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
    virtual int EvalCentralCell();
    int SolvePoissonEquation(MRICommunicator* comm);

    // REYNOLDS STRESS COMPUTATION    
    void EvalReynoldsStressGradient(int currentCell, double** ReynoldsStressGradient);
    
    // APPLY SMOOTHING FILTER - LAVISION
    void ApplySmoothingFilter();

    // APPLY MEDIAN FILTER
    void ApplyMedianFilter(int qtyID,int maxIt);
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
    virtual void   EvalCellVelocityGradientDecomposition(int currentCell, double** deformation, double** rotation, double** firstDerivs);
    virtual double EvalCellQCriterion(int currentCell, double** deformation, double** rotation);
    virtual double EvalCellL2Criterion(int currentCell, double** deformation, double** rotation);
    virtual double EvalCellDeltaCriterion(int currentCell, double** deformation, double** rotation, double** velGradient);
    virtual double EvalCellVortexCriteria(int currentCell,int criteriaType, double** deformation, double** rotation, double** velGradient);
    virtual void   EvalVortexCriteria();
    virtual void   EvalVorticity();
    virtual void   EvalEnstrophy();
    virtual void   EvalSMPVortexCriteria(MRIExpansion* exp);

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
   double GetDiffNorm(MRIStructuredScan* otherScan);

   void buildMetisConnectivities(int *eptr,int *eind);

   // MESSAGE PASSING
   void formVortexList(int totVortex,int* minFace,int* maxFace,MRIIntVec& innerVortexList,MRIIntVec& boundaryVortexList,MRICommunicator* comm);
   void passScanData(MRICommunicator* comm);
   virtual void DistributeScanData(MRICommunicator* comm);
};

#endif // MRISTRUCTUREDSCAN_H
