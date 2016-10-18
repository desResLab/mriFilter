#ifndef MRISCAN_H
#define MRISCAN_H

# include <math.h>
# include <stdio.h>
# include <string>
# include <limits>
# include <vector>
# include <iostream>
# include <fstream>
# include <boost/algorithm/string.hpp>

# include "mriUtils.h"
# include "mriCell.h"
# include "mriTypes.h"
# include "mriExpansion.h"
# include "mriThresholdCriteria.h"
# include "mriImagedata.h"
# include "mriConstants.h"
# include "mriException.h"
# include "mriSamplingOptions.h"
# include "mriOutput.h"
# include "mriTopology.h"

# include "mriOptions.h"
# include "mriCommunicator.h"

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
    vector<MRICell> cells;
    MRIIntVec       cellTags;
    // Output Quantities
    vector<MRIOutput> outputs;
    // Utility Functions
    bool hasPressureGradient;
    bool hasRelativePressure;
    bool hasReynoldsStress;
    double scanTime;
    double maxVelModule;
    // MRI Expansion
    MRIExpansion* expansion;

    // POINTER TO A TOPOLOGY 
    MRITopology* topology;

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

    // ==============
    // READ FUNCTIONS
    // ==============
    void createGridFromVTKStructuredPoints(vtkStructuredPointsOptionRecord opts);
    void readVTKStructuredPoints(std::string vtkFileName, bool DoReorderCells);
    void readPLTFile(std::string PltFileName, bool DoReorderCells);
    void readFromExpansionFile(std::string fileName,bool applyThreshold,int thresholdType,double thresholdValue);

    // ==================
    // READ FROM RAW DATA
    // ==================
    void readRAWFileSequence(std::string fileListName);
    int  readRawImage(std::string FileName, MRIImageData &data);
	
    // ===============
    // WRITE FUNCTIONS
    // ===============
    void fillPLTHeader(std::vector<std::string> &pltHeader, bool isFirstFile);
    void exportToLSDYNA(std::string LSFileName, double scale);
    void exportToCSV(std::string FileName);
    void exportNodesToFile(std::string FileName);
    void flushToFile(std::string FileName);
    void exportVelocitiesToFile(std::string fileName, bool append);
    // VIRTUAL
    void exportToVOL(std::string FileName);
    void exportToTECPLOT(std::string FileName, bool isFirstFile);
    void exportToVTK(std::string fileName, MRIThresholdCriteria* threshold);
    void writeExpansionFile(std::string fileName);
    // Export to Poisson Solver Only element with significant concentration
    void exportForDistancing(string inputFileName, MRIThresholdCriteria* threshold);
    void exportForPoisson(string inputFileName,double density,double viscosity,MRIThresholdCriteria* threshold, const MRIDoubleMat& timeDerivs,
                                  bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                                  bool readMuTFromFile, string muTFile, double smagorinskyCoeff);

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
    int  getCellFaceID(int CellId,int FaceId);
    bool hasUniformSpacing();
    // NORMALS
    void getExternalFaceNormal(int cellID, int localFaceID, double* extNormal);

    // =============================
    // TRANSFORMATIONS AND THRESHOLD
    // =============================
    void crop(const MRIDoubleVec& limitBox);
    void scaleVelocities(double factor);
    void scalePositions(double factor);

    // ==========================================
    // RECONSTRUCTION FROM EXPANSION COEFFICIENTS
    // ==========================================
    void rebuildFromExpansion(MRIExpansion* expansion,bool useConstantFlux);

    // ===============================
    // RECONSTRUCTION FROM FACE FLUXES
    // ==============================
    void rebuildFromFaceFluxes(double* faceFluxes);

    // Reorder Cells    
    void reorderCells(std::vector<int> Perm);
    // Get Global Permutation
    void getGlobalPermutation(std::vector<int> &GlobalPerm);
    // REORDER GLOBAL SCAN
    void reorderScan();

    // MAPPING FUNCTIONS
    // Get Cell Number From Coords
    int  getCellNumber(double* coords);
    void getCartesianNeighbourCells(int CurrentCell, std::vector<int> &coords, bool addself);
    void getStructuredNeighbourCells(int centreCell,int order, MRIThresholdCriteria* threshold, MRIIntVec& cellNeighbors);
    bool isCompatibleWith(MRIScan* secondScan);
    // Get Face from Cell Vector
    int  getFacewithCellVector(int CurrentCell, double *UnitVector);
    // Get Adjacency Face
    int  getAdjacentFace(int GlobalNodeNumber, int AdjType);
    // Get Unit Vector From Current Cell To Face Centre
    void getGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, MRIDoubleVec& globalCoords);
    int  faceLocaltoGlobal(int LocalFace, int DimNumber, int SliceNumber);

    // Sequential Index to Integer Coords
    void mapIndexToCoords(int index, MRIIntVec& intCoords);
    void mapIndexToAuxNodeCoords(int index, MRIIntVec& intCoords);
    int  mapCoordsToIndex(int i, int j, int k);

    // Map Integer Coords to Position
    void mapCoordsToPosition(int* coords, bool addMeshMinima, double* pos);
    void mapAuxCoordsToPosition(int* auxCoords, double* pos);
    void getLocalStarFaces(int StarNum, int CellsX, int CellsY, int &BottomFace, int &TopFace, int &LeftFace, int &RightFace);
    int  findFirstNotVisited(int cellTotal, bool* visitedCell, MRIIntVec cellStack);
    void formNotVisitedList(int cellTotal, const MRIBoolVec& visitedCell, MRIBoolVec& notVisitedList);

    // Map cell vector to face vector
    void cellToFace(bool deleteWalls, MRIThresholdCriteria* thresholdCriteria,MRIDoubleMat cellVec, MRIDoubleVec &faceVec);
    void cellToFacePartial(MRIIntVec elUsageMap, MRIThresholdCriteria* thresholdCriteria,
                           MRIDoubleMat cellVec, MRIDoubleVec &faceVec);

    // =========
    // MP FILTER
    // =========
    int    getTotalBasisNumber();
    void   applySMPFilter(MRIOptions* options, bool isBC, MRICommunicator* comm);
    void   assembleResidualVector(bool useBCFilter, MRIThresholdCriteria* thresholdCriteria, int &totalFaces, double* &ResVec, double* &filteredVec, double &resNorm);
    void   assembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs);
    void   assembleConstantPatternMPI(int currentDim, int &totalConstantFacesOnProc,
                                              std::vector<int> &facesIDOnProc, std::vector<double> &facesCoeffsOnProc,
                                              int minFaceOnProc, int maxFaceOnProc,MRICommunicator* comm);
    void   assembleStarShape(int vortexNumber, int &totalFaces,std::vector<int> &facesID,std::vector<double> &facesCoeffs);
    double evalMaxDivergence(double* filteredVec);
    void   recoverGlobalErrorEstimates(double& AvNormError, double& AvAngleError);
    void   expandStarShape(int totalStarFaces, int* facesID, double* facesCoeffs, double* &fullStarVector);
    void   recoverCellVelocitiesRT0(bool useBCFilter, double* filteredVec);
    void   reconstructFromExpansion();
    void   getDimensionSliceStarFromVortex(int vortexNumber,int &dimNumber,int &sliceNumber,int &starNumber);
    
    // FILTER MATRICES
    void   assembleEncodingMatrix(int &totalFaces, int &totalBasis, MRIDoubleMat& Mat);
    void   assembleDecodingMatrix(int &totalFaces, int &totalBasis, MRIDoubleMat& Mat);
    void   assembleStarMatrix(int &totalFaces, int &totalBasis, MRIDoubleMat& Matrix);
      
    // ==============
    // INFO FUNCTIONS
    // ==============
    double evalAverageVelocityMod();
    void   evalCellAreas(int cellNumber,double* Areas);
    int    getTotalAuxNodes();
    int    evalTotalVortex();
    int    getTotalFaces();

    // CELL SAMPLING
    void sampleVelocities(MRISamplingOptions SamplingOptions, MRIIntVec& bins);

    // GRADIENTS AND DERIVATIVES
    void evalSpaceDerivs(int currentCell, MRIThresholdCriteria* threshold, MRIDoubleMat& firstDerivs, MRIDoubleMat& secondDerivs);
    void evalSpaceGradient(int currentCell,int qtyID, double* gradient);
    void computeQuantityGradient(int qtyID);

    // DIVERGENCE
    MRIDoubleVec evalCellDivergences(MRIDoubleVec faceVec);
  
    // PRESSURE COMPUTATION
    void evalRelativePressure(int startingCell, double refPressure);
    void performPressureIterations();

    // OTHERS
    void evalPressureIterative(int currentCell, double currentValue, bool* visitedCell,int* otherCells, std::vector<int> &cellStack,int& cellCount);
    bool areThereNotVisitedNeighbor(int cell, bool* visitedCell);
    bool areThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor);
    int  getCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished, bool& secondStage);
    int  getNextStartingCell(int currentCell, bool* visitedCell, bool* isBoundaryCell, bool &finished, int &bookmark);
    int  evalCentralCell();
    void updateVelocities();

    // REYNOLDS STRESS COMPUTATION    
    void evalReynoldsStress(MRIThresholdCriteria* threshold);
    void evalReynoldsStressGradient(int currentCell, MRIDoubleMat& ReynoldsStressGradient);
    void evalPradtlTurbViscosity(MRIDoubleMat cellDistance, MRIThresholdCriteria* threshold, double density, MRIDoubleMat& turbViscosity);
    void evalSmagorinskyLillyTurbViscosity(double density, double smagorinskyCoeff, MRIThresholdCriteria* threshold, MRIDoubleMat& turbViscosity);
    void evalEddyViscosity_simple(MRIDoubleVec& nuT);
    void evalTurbulentKineticEnergy(MRIThresholdCriteria* threshold, MRIDoubleVec& turbK);
    
    // APPLY SMOOTHING FILTER - LAVISION
    void applySmoothingFilter();
    void applyMedianFilter(int qtyID,int maxIt,int order,int filterType,MRIThresholdCriteria* threshold);

    // THRESHOLD
    void thresholdQuantity(int qtyID,double threshold);
    void evalNoisyPressureGradientPoints();
  
    // ==============
    // TEMPLATE FLOWS
    // ==============
    void createSampleCase(MRISamples sampleType, const MRIDoubleVec& params);
    void assignVelocitySignature(MRIDirection dir, MRISamples sample, double currTime);
    void assignConstantSignature(MRIDirection dir);
    void assignStagnationFlowSignature(MRIDirection dir);
    void assignPoiseilleSignature(MRIDirection dir);
    void assignCylindricalFlowSignature(MRIDirection dir);
    void assignSphericalFlowSignature(MRIDirection dir);
    void assignToroidalVortexFlowSignature();
    void assignTimeDependentPoiseilleSignature(double omega, double radius, double viscosity, double currtime, double maxVel);
    void assignConstantFlowWithStep();
    void assignTaylorVortexSignature(MRIDirection dir);
    void assignRandomStandardGaussianFlow();
    void assignRandomComponent(const int kdirX,stdRndGenerator &generator);
    void assignZeroVelocities();

    // VORTEX IDENTIFICATION
    void   evalCellVelocityGradientDecomposition(int currentCell, double** deformation, double** rotation, double** firstDerivs);
    double evalCellQCriterion(int currentCell, double** deformation, double** rotation);
    double evalCellL2Criterion(int currentCell, double** deformation, double** rotation);
    double evalCellDeltaCriterion(int currentCell, double** deformation, double** rotation, double** velGradient);
    double evalCellVortexCriteria(int currentCell,int criteriaType, double** deformation, double** rotation, double** velGradient);
    void   evalVortexCriteria(MRIThresholdCriteria* threshold);
    void   evalVorticity(MRIThresholdCriteria* threshold);
    void   evalEnstrophy(MRIThresholdCriteria* threshold);
    void   evalSMPVortexCriteria(MRIExpansion* exp);

    // SPATIAL REPRESENTATION OF VORTEX COEFFICIENTS
    double evalVortexCriteria(MRIExpansion* exp);
    void   getNeighborVortexes(int cellNumber,int dim, MRIIntVec& idx);
    
    // ADD GAUSSIAN NOISE
    void applyGaussianNoise(double stDev);
    
    // TESTING FUNCTIONALITIES
    void testScanAdjacency(string fileName);
    
   // COMPARISON BETWEEN SCANS
   double getDiffNorm(MRIScan* otherScan);

   void buildMetisConnectivities(int *eptr,int *eind);

   // MESSAGE PASSING
   void formVortexList(int totVortex,int* minFace,int* maxFace,MRIIntVec& innerVortexList,MRIIntVec& boundaryVortexList,MRICommunicator* comm);
   void passScanData(MRICommunicator* comm);
   void distributeScanData(MRICommunicator* comm);

   // BOUNDARY CLEANING
   void cleanNormalComponentOnBoundary(MRIThresholdCriteria* threshold);
   void interpolateBoundaryVelocities(MRIThresholdCriteria* threshold);
   void projectCellVelocity(int cell,double* normal);
   int  getOppositeCell(int cell, double* normal);
   void tagByNeighbour(int tag,int* cellTags, bool* isTaggable,int startingCell);
   bool hasUntaggedNeighbours(int cell,int* cellTags, bool* isTaggable);
   void setWallFluxesToZero(bool* isFaceOnWalls, MRIDoubleVec& poissonSourceFaceVec);

   // STATISTICS
   void evalScanPDF(int pdfQuantity, int numberOfBins, bool useBox, MRIDoubleVec& limitBox,MRIDoubleVec& binCenter, MRIDoubleVec& binArray);
   void formBinLimits(int pdfQuantity, double& currInterval, const MRIDoubleVec& limitBox, int numberOfBins, MRIDoubleVec& binMin, MRIDoubleVec& binMax, MRIDoubleVec& binCenter);
};

#endif // MRISTRUCTUREDSCAN_H
