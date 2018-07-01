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
# include "mriIO.h"

# include "mriOptions.h"
# include "mriCommunicator.h"

class mriScan;
class mriCommunicator;

// ========================
// UNSTRUCTURED GRID LAYOUT
// ========================
class mriScan{
  public:
    // ============
    // DATA MEMBERS
    // ============

    // Cell Data
    vector<mriCell> cells;
    mriIntVec       cellTags;
    
    // Output Quantities
    vector<mriOutput> outputs;
    
    // Utility Functions
    mriDoubleMat reynoldsStress;
    mriDoubleMat qtyGradient;
    double scanTime;
    double maxVelModule;
    
    // mri Expansion
    mriExpansion* expansion;

    // Pointer to the Common Topology of the Sequence
    mriTopology* topology;

    // Material Properties
    double density;
    double viscosity;

    // ================
    // MEMBER FUNCTIONS
    // ================
    // Constructor
    // 1 - Empty
    mriScan(double currentTime);
    // 2 - Create a Zero Scan from another Scan
    mriScan(const mriScan& copyScan);
    // Destructor
    ~mriScan(){}
    // INFO FUNCTIONS
    string writeStatistics();

    // ==============
    // READ FUNCTIONS
    // ==============
    void readFromVTK_ASCII(string vtkFileName, const vtkStructuredPointsOptionRecord& vtkOptions);
    void readFromPLT_ASCII(std::string PltFileName, const pltOptionRecord& pltOptions);
    void readFromExpansionFile(std::string fileName,bool applyThreshold,int thresholdType,double thresholdValue);

    // ==================
    // READ FROM RAW DATA
    // ==================
    void readRAWFileSequence(std::string fileListName);
    int  readRawImage(std::string FileName, mriImageData &data);
	
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
    void exportToVTK(std::string fileName, mriThresholdCriteria* threshold);
    void writeExpansionFile(std::string fileName);
    // Export to Poisson Solver Only element with significant concentration
    void exportForDistancing(string inputFileName, mriThresholdCriteria* threshold);
    void exportForPoisson(string inputFileName,double density,double viscosity,mriThresholdCriteria* threshold, const mriDoubleMat& timeDerivs,
                          bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                          bool readMuTFromFile, string muTFile, double smagorinskyCoeff);

    // ========
    // TOPOLOGY
    // ========
    // VORTEX COEFFICIENTS
    double getEdgeFaceVortexCoeff(int edgeID, int faceID);
    // AUXILIARY GEOMETRY
    void   getEdgeCenter(int edgeID, double* ec);
    void   getFaceCenter(int faceID, double* fc);
    void   getEdgeToFaceDirection(int edgeID, int faceID, double* edgeFaceVector);
    // VOLUME
    double evalCellVolume(int cellNumber);
    // OTHER
    int  getCellFaceID(int CellId,int FaceId);
    bool hasUniformSpacing();
    // NORMALS
    bool isInnerCell(int cell);

    // =============================
    // TRANSFORMATIONS AND THRESHOLD
    // =============================
    void crop(const mriDoubleVec& limitBox, const mriBoolVec& indexesToCrop);
    void scaleVelocities(double factor);
    void scalePositions(double factor);

    // ==========================================
    // RECONSTRUCTION FROM EXPANSION COEFFICIENTS
    // ==========================================
    void rebuildFromExpansion(mriExpansion* expansion,bool useConstantFlux);

    // ===============================
    // RECONSTRUCTION FROM FACE FLUXES
    // ==============================
    void rebuildFromFaceFluxes(double* faceFluxes);

    // ============
    // THRESHOLDING
    // ============
    void applyThresholding(mriThresholdCriteria* thresholdCriteria);

    // SAVE VELOCITIES
    void saveVelocity();

    // MAPPING FUNCTIONS
    // Get Cell Number From Coords
    int  getCellNumber(double* coords);
    void getCartesianNeighbourCells(int CurrentCell, std::vector<int> &coords, bool addself);
    void getStructuredNeighbourCells(int centreCell,int order, mriThresholdCriteria* threshold, mriIntVec& cellNeighbors);
    // Get Face from Cell Vector
    int  getFacewithCellVector(int CurrentCell, double *UnitVector);
    // Get Unit Vector From Current Cell To Face Centre
    void getGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, mriDoubleVec& globalCoords);
    int  faceLocaltoGlobal(int LocalFace, int DimNumber, int SliceNumber);

    // Map Integer Coords to Position
    void getLocalStarFaces(int StarNum, int CellsX, int CellsY, int &BottomFace, int &TopFace, int &LeftFace, int &RightFace);
    int  findFirstNotVisited(int cellTotal, bool* visitedCell, mriIntVec cellStack);
    void formNotVisitedList(int cellTotal, const mriBoolVec& visitedCell, mriBoolVec& notVisitedList);

    // Map cell vector to face vector
    void cellToFace(bool deleteWalls, mriThresholdCriteria* thresholdCriteria,mriDoubleMat cellVec, mriDoubleVec &faceVec);
    void cellToFacePartial(mriIntVec elUsageMap, mriThresholdCriteria* thresholdCriteria,
                           mriDoubleMat cellVec, mriDoubleVec &faceVec);

    // =========
    // MP FILTER
    // =========
    int    getTotalBasisNumber();
    void   applySMPFilter(mriCommunicator* comm, bool isBC, 
                          mriThresholdCriteria* thresholdCriteria,
                          double itTol,
                          int maxIt,
                          bool useConstantPatterns);
    void   assembleResidualVector(bool useBCFilter, 
                                  mriThresholdCriteria* thresholdCriteria,
                                  int& totalFaces, 
                                  mriDoubleVec& resVec, 
                                  mriDoubleVec& filteredVec, 
                                  double& resNorm);
    void   assembleConstantPattern(int currentDim, int &totalConstantFaces, std::vector<int> &facesID, std::vector<double> &facesCoeffs);
    void   assembleConstantPatternMPI(int currentDim, int &totalConstantFacesOnProc,
                                              std::vector<int> &facesIDOnProc, std::vector<double> &facesCoeffsOnProc,
                                              int minFaceOnProc, int maxFaceOnProc,mriCommunicator* comm);
    void   assembleStarShape(int vortexNumber, int &totalFaces,std::vector<int> &facesID,std::vector<double> &facesCoeffs);
    double evalMaxDivergence(const mriDoubleVec& filteredVec);
    void   recoverGlobalErrorEstimates(double& AvNormError, double& AvAngleError);
    void   expandStarShape(int totalStarFaces, int* facesID, double* facesCoeffs, double* &fullStarVector);
    void   recoverCellVelocitiesRT0(bool useBCFilter, mriDoubleVec& filteredVec);
    void   reconstructFromExpansion();
    void   getDimensionSliceStarFromVortex(int vortexNumber,int &dimNumber,int &sliceNumber,int &starNumber);
    
    // FILTER MATRICES
    void   assembleEncodingMatrix(int &totalFaces, int &totalBasis, mriDoubleMat& Mat);
    void   assembleDecodingMatrix(int &totalFaces, int &totalBasis, mriDoubleMat& Mat);
    void   assembleStarMatrix(int &totalFaces, int &totalBasis, mriDoubleMat& Matrix);
      
    // ==============
    // INFO FUNCTIONS
    // ==============
    double evalAverageVelocityMod();
    void   evalCellAreas(int cellNumber,mriDoubleVec& Areas);
    int    evalTotalVortex();

    // CELL SAMPLING
    void sampleVelocities(mriSamplingOptions SamplingOptions, mriIntVec& bins);

    // GRADIENTS AND DERIVATIVES
    void evalSpaceDerivs(int currentCell, mriThresholdCriteria* threshold, mriDoubleMat& firstDerivs, mriDoubleMat& secondDerivs);
    void evalSpaceGradient(int currentCell,int qtyID, mriDoubleVec& gradient);
    void computeQuantityGradient(int qtyID);

    // DIVERGENCE
    void evalCellDivergences(const mriDoubleVec& faceVec,mriDoubleVec& cellDivs);
  
    // OTHERS
    void evalPressureIterative(int currentCell, double currentValue, bool* visitedCell,int* otherCells, std::vector<int> &cellStack,int& cellCount);
    bool areThereNotVisitedNeighbor(int cell, bool* visitedCell);
    bool areThereVisitedNeighbor(int cell, bool* visitedCell, bool* isBoundaryCell, int &visitedNeighbor);
    int  getCellFromStack(std::vector<int> &cellStack, bool* visitedCell, bool* isBoundaryCell, bool &finished, bool& secondStage);
    int  getNextStartingCell(int currentCell, bool* visitedCell, bool* isBoundaryCell, bool &finished, int &bookmark);
    int  evalCentralCell();
    void updateVelocities();

    // REYNOLDS STRESS COMPUTATION    
    void evalReynoldsStress(mriThresholdCriteria* threshold);
    void evalReynoldsStressGradient(int currentCell, mriDoubleMat& ReynoldsStressGradient);
    void evalPradtlTurbViscosity(mriDoubleMat cellDistance, mriThresholdCriteria* threshold, double density, mriDoubleMat& turbViscosity);
    void evalSmagorinskyLillyTurbViscosity(double density, double smagorinskyCoeff, mriThresholdCriteria* threshold, mriDoubleMat& turbViscosity);
    void evalEddyViscosity_simple(mriDoubleVec& nuT);
    void evalTurbulentKineticEnergy(mriThresholdCriteria* threshold, mriDoubleVec& turbK);
    
    // APPLY SMOOTHING FILTER - LAVISION
    void applySmoothingFilter();
    void applyMedianFilter(int qtyID,int maxIt,int order,int filterType,mriThresholdCriteria* threshold);

    // THRESHOLD
    void thresholdQuantity(int qtyID,double threshold);
    void evalNoisyPressureGradientPoints();
  
    // ==============
    // TEMPLATE FLOWS
    // ==============
    void createFromTemplate(mriTemplateType sampleType,const mriDoubleVec& params);
    void assignVelocitySignature(mriDirection dir, mriTemplateType sample, double currTime, const mriDoubleVec& auxParams);
    void assignConstantSignature(mriDirection dir);
    void assignStagnationFlowSignature(mriDirection dir);
    void assignPoiseuilleSignature(mriDirection dir, const mriDoubleVec& auxParams);
    void assignCylindricalFlowSignature(mriDirection dir);
    void assignSphericalFlowSignature(mriDirection dir);
    void assignToroidalVortexFlowSignature();
    void assignTimeDependentPoiseilleSignature(double omega, double radius, double viscosity, double currtime, double maxVel);
    void assignConstantFlowWithStep();
    void assignTaylorVortexSignature(mriDirection dir);
    void assignRandomStandardGaussianFlow();
    void assignRandomComponent(const int kdirX,stdRndGenerator &generator);
    void assignZeroVelocities();

    // VORTEX IDENTIFICATION
    double evalCellQCriterion(int currentCell, const mriDoubleMat& deformation, const mriDoubleMat& rotation);
    double evalCellL2Criterion(int currentCell, const mriDoubleMat& deformation, const mriDoubleMat& rotation);
    double evalCellDeltaCriterion(int currentCell, const mriDoubleMat& deformation, const mriDoubleMat& rotation, const mriDoubleMat& velGradient);
    double evalCellVortexCriteria(int currentCell,int criteriaType, const mriDoubleMat& velGradient, const mriDoubleMat& deformation, const mriDoubleMat& rotation);
    void   evalVortexCriteria(mriThresholdCriteria* threshold);
    void   evalVorticity(mriThresholdCriteria* threshold);
    void   evalEnstrophy(mriThresholdCriteria* threshold);
    void   evalSMPVortexCriteria(mriExpansion* exp);

    // SPATIAL REPRESENTATION OF VORTEX COEFFICIENTS
    double evalVortexCriteria(mriExpansion* exp);
    
    // ADD GAUSSIAN NOISE
    void applyGaussianNoise(double stDev, double seed);
    
    // TESTING FUNCTIONALITIES
    void testScanAdjacency(string fileName);
    
   // COMPARISON BETWEEN SCANS
   double getDiffNorm(mriScan* otherScan);

   void buildMetisConnectivities(int *eptr,int *eind);

   // MESSAGE PASSING
   void formVortexList(mriCommunicator* comm,
                       int totVortex,
                       const mriIntVec& minFace,
                       const mriIntVec& maxFace,
                       mriIntVec& innerVortexList,
                       mriIntVec& boundaryVortexList);

   void passScanData(mriCommunicator* comm);
   void distributeScanData(mriCommunicator* comm);

   // BOUNDARY CLEANING
   void cleanNormalComponentOnBoundary(mriThresholdCriteria* threshold);
   void interpolateBoundaryVelocities(mriThresholdCriteria* threshold);
   void projectCellVelocity(int cell,double* normal);
   int  getOppositeCell(int cell, double* normal);
   void tagByNeighbour(int tag,int* cellTags, bool* isTaggable,int startingCell);
   bool hasUntaggedNeighbours(int cell,int* cellTags, bool* isTaggable);
   void setWallFluxesToZero(bool* isFaceOnWalls, mriDoubleVec& poissonSourceFaceVec);

   // STATISTICS
   void evalScanPDF(int pdfQuantity, int numberOfBins, bool useBox, mriDoubleVec& limitBox,mriDoubleVec& binCenter, mriDoubleVec& binArray);
   void formBinLimits(int pdfQuantity, double& currInterval, const mriDoubleVec& limitBox, int numberOfBins, mriDoubleVec& binMin, mriDoubleVec& binMax, mriDoubleVec& binCenter);

   // PRESSURE
   void evalCellPressureGradients(int currentCell,
                                  const mriDoubleVec& timeDeriv, 
                                  const mriDoubleMat& firstDerivs, 
                                  const mriDoubleMat& secondDerivs,
                                  const mriDoubleMat& ReynoldsStressGrad,
                                  mriDoubleVec& pressureGrad);
};

#endif // mriSTRUCTUREDSCAN_H
