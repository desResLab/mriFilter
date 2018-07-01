#ifndef MRISEQUENCE_H
#define MRISEQUENCE_H

#include "mriScan.h"
#include "mriCommunicator.h"
#include "mriThresholdCriteria.h"
#include "mriTopology.h"
#include "mriIO.h"

using namespace std;

class mriScan;
class mriTopology;
class mriCommunicator;

// Generic Sequence Containing General 
// Unstructured Data Sets

class mriSequence{
  public:
    vector<mriScan*> sequence;
    string*          fileNames;
    bool             isCyclic;

    // TOPOLOGY COMMON TO EVERY SCAN
    mriTopology* topology;

    // Constructor and
    mriSequence(bool cyclic);
    // Copy Constructor
    mriSequence(mriSequence* copySequence);
    // Destructor
    ~mriSequence();

    // WRITE STATS
    string writeStatistics();
    
    // ADD AND GET FROM SEQUENCE
    void addScan(mriScan* scan);
    int  getTotalScans(){return sequence.size();}
    bool getIsCyclic(){return isCyclic;}
    mriScan* getScan(int scanNumber);

    // READ SEQUENCE FROM FILE
    void readFromASCIISequence(int asciiInputType, 
                               const mriStringVec& vtkFileNames, 
                               const mriDoubleVec& times);    
    void readFromExpansionFiles(const mriStringVec& fileNames, 
                                const mriDoubleVec& Times, 
                                bool applyThreshold, 
                                int thresholdType,
                                double thresholdRatio);
    
    // EXPORT SEQUENCE TO FILE
    void exportToTECPLOT(string outfileName);
    void exportToVOL(string outfileName);
    void exportToVTK(string outfileName,mriThresholdCriteria* thresholdCriteria);    
    void exportForDistancing(string inputFileName, mriThresholdCriteria* threshold);
    void exportForPoisson(string inputFileName,double density,double viscosity,mriThresholdCriteria* threshold,
                          bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                          bool readMuTFromFile, string muTFile, double smagorinskyCoeff);
    void writeExpansionFile(string fileName);    

    // SAVE QUANTITIES TO OUTPUTS
    void   saveVelocity();
    double getVelocityNormAtCell(int cell);

    // TOPOLOGY AND MAPPING
    void createTopology();
    void getUnitVector(int CurrentCell, const mriDoubleVec& GlobalFaceCoords, mriDoubleVec& myVect);
    void getGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, mriDoubleVec& globalCoords);
    int  getCellNumber(const mriDoubleVec& coords);
    void getGlobalPermutation(mriIntVec& GlobalPerm);
    int  getAdjacentFace(int globalNodeNumber, int AdjType);
    void getNeighborVortexes(int cellNumber,int dim,mriIntVec& idx);
    void getEdgeDirection(int edgeID, double* edgeDirVector);
    void getAuxNodeCoordinates(int nodeNum, mriDoubleVec& pos);
    void buildCellConnections();
    int  addToFaceConnections(const mriIntVec& faceIds, vector<vector<mriFace* > >& AuxFirstNodeFaceList);
    int  addToEdgeConnections(const mriIntVec& edgeIds,vector<vector<mriEdge*> >& AuxFirstNodeEdgeList);
    void buildFaceConnections();
    void buildFaceCells();
    void buildEdgeConnections();
    
    // MAPPING     
    int  mapCoordsToIndex(int i, int j, int k);
    void mapCoordsToPosition(const mriIntVec& coords, bool addMeshMinima, mriDoubleVec& pos);
    void mapIndexToAuxNodeCoords(int index, mriIntVec& intCoords);
    void mapAuxCoordsToPosition(const mriIntVec& auxCoords, mriDoubleVec& pos);    
    
    // DIV FREE Filtering
    void applySMPFilter(mriCommunicator* comm, bool isBC, 
                        mriThresholdCriteria* thresholdCriteria,
                        double itTol,
                        int maxIt,
                        bool useConstantPatterns);
    
    // APPLY THRESHOLDING 
    void applyThresholding(mriThresholdCriteria* thresholdCriteria);

    // EVAL VORTEX CRITERIA
    void evalVortexCriteria(mriThresholdCriteria* thresholdCriteria);
    void evalVorticity(mriThresholdCriteria* thresholdCriteria);
    void evalEnstrophy(mriThresholdCriteria* thresholdCriteria);
    void evalSMPVortexCriteria();
    
    // PRESSURE COMPUTATION
    void computePressureGradients(mriThresholdCriteria* threshold);

    // COMPUTING REYNOLDS STRESSES
    void evalReynoldsStresses(mriThresholdCriteria* threshold);
    void evalScanReynoldsStressDerivs(int currentScan,mriDoubleMat& reynoldsDeriv);

    // ADD NOISE
    void applyNoise(double noiseIntensity, double seed);

    // FILTER DATA
    void applyMedianFilter(int qtyID,int maxIt,int order,int filterType,mriThresholdCriteria* threshold);
    
    // File List Printing
    void printSequenceFiles(std::string outFIleName);
    // Operations on single Scans
    void makeScanDifference(int firstScanID, int secondScanID);
    void makeScanAverage(int numberOfMeasures, int firstScanID, int secondScanID);
    // Eval Time Derivatives
    void evalTimeDerivs(int currentScan, int currentCell, mriDoubleVec& timeDeriv);
    void evalScanTimeDerivs(int currentScan,mriDoubleMat& timeDeriv);
  
    // STATISTICS
    void evalScanDifferencePDF(int otherScan, int refScan, const int pdfQuantity, int numberOfBins, bool useBox, mriDoubleVec& limitBox, mriDoubleVec& binCenters, mriDoubleVec& binArray);
    void extractSinglePointTimeCurve(int cellNumber, int exportQty, string fileName);
    void formDifferenceBinLimits(int otherScan, int refScan, 
                                 int pdfQuantity, double& currInterval,
                                 const mriDoubleVec& limitBox, 
                                 int numberOfBins, 
                                 mriDoubleVec& binMin, 
                                 mriDoubleVec& binMax, 
                                 mriDoubleVec& binCenter);
    
    // TRANFORMATION
    void crop(const mriDoubleVec& limitBox);
    void scaleVelocities(double factor);
    void scalePositions(const mriDoubleVec& origin, double factor);

    // MESSAGE PASSING
    void distributeSequenceData(mriCommunicator* comm);

    // CLEAN COMPONENTS ON BOUNDARY
    void cleanNormalComponentOnBoundary(mriThresholdCriteria* threshold);
    void interpolateBoundaryVelocities(mriThresholdCriteria* threshold);

    // TEMPLATE FLOW SEQUENCE
    void createSampleCase(mriTemplateType sampleType,const mriDoubleVec& params);
};

#endif // MRISEQUENCE_H
