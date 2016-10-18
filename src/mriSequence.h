#ifndef MRISEQUENCE_H
#define MRISEQUENCE_H

#include "mriScan.h"
#include "mriCommunicator.h"

using namespace std;

class MRIScan;
class MRICommunicator;

// Generic Sequence Containing General 
// Unstructured Data Sets

class MRISequence{
  public:
    int                   totalScans;
    std::vector<MRIScan*> sequence;
    std::string*          fileNames;
    bool                  isCyclic;

    // TOPOLOGY COMMON TO EVERY SCAN
    MRITopology* topology;

    // Constructor and
    MRISequence(bool cyclic);
    // Copy Constructor
    MRISequence(MRISequence* copySequence);
    // Destructor
    ~MRISequence();
    
    // ADD AND GET FROM SEQUENCE
    void addScan(MRIScan* scan);
    int  getTotalScans(){return totalScans;}
    bool getIsCyclic(){return isCyclic;}
    MRIScan* getScan(int scanNumber);

    // READ SEQUENCE FROM FILE
    void readPLTFile(string PltFileName, bool DoReorderCells);
    void readFromVolSequence(string outfileName);
    void readFromExpansionFile(string fileName,bool applyThreshold, int thresholdType,double thresholdRatio);
    
    // EXPORT SEQUENCE TO FILE
    void exportToTECPLOT(string outfileName);
    void exportToVOL(string outfileName);
    void exportToVTK(string outfileName,MRIThresholdCriteria* thresholdCriteria);
    void writeExpansionFile(string fileName);
    void exportForDistancing(string inputFileName, MRIThresholdCriteria* threshold);
    void exportForPoisson(string inputFileName,double density,double viscosity,MRIThresholdCriteria* threshold,
                          bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                          bool readMuTFromFile, string muTFile, double smagorinskyCoeff);

    // SAVE QUANTITIES TO OUTPUTS
    void   saveVelocity();
    double getVelocityNormAtCell(int cell);

    // TOPOLOGY AND MAPPING
    void createTopology();
    bool isInnerCell(int cell);
    void getCartesianNeighbourCells(int CurrentCell,MRIIntVec& cellNeighbors, bool addself);
    void getUnitVector(int CurrentCell, const MRIDoubleVec& GlobalFaceCoords, MRIDoubleVec& myVect);
    void getGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, MRIDoubleVec& globalCoords);
    int  getCellNumber(const MRIDoubleVec& coords);
    void getGlobalPermutation(MRIIntVec& GlobalPerm);
    int  getAdjacentFace(int globalNodeNumber, int AdjType);
    void getNeighborVortexes(int cellNumber,int dim,MRIIntVec& idx);
    void getEdgeDirection(int edgeID, double* edgeDirVector);
    void getAuxNodeCoordinates(int nodeNum, MRIDoubleVec& pos);
    
    // MAPPING 
    void mapIndexToCoords(int index, MRIIntVec& intCoords);
    int  mapCoordsToIndex(int i, int j, int k);
    void mapCoordsToPosition(const MRIIntVec& coords, bool addMeshMinima, MRIDoubleVec& pos);
    void mapIndexToAuxNodeCoords(int index, MRIIntVec& intCoords);
    void mapAuxCoordsToPosition(const MRIIntVec& auxCoords, MRIDoubleVec& pos);
    
    // DIV FREE Filtering
    void applySMPFilter(MRIOptions* options, bool isBC, MRICommunicator* comm);
    
    // APPLY THRESHOLDING 
    void applyThresholding(MRIThresholdCriteria* thresholdCriteria);

    // EVAL VORTEX CRITERIA
    void evalVortexCriteria(MRIThresholdCriteria* thresholdCriteria);
    void evalVorticity(MRIThresholdCriteria* thresholdCriteria);
    void evalEnstrophy(MRIThresholdCriteria* thresholdCriteria);
    void evalSMPVortexCriteria();
    
    // PRESSURE COMPUTATION
    void computePressureGradients(MRIThresholdCriteria* threshold);
    void computeRelativePressure(bool doPressureSmoothing);

    // COMPUTING REYNOLDS STRESSES
    void evalReynoldsStresses(MRIThresholdCriteria* threshold);
    void evalScanReynoldsStressDerivs(int currentScan,MRIDoubleMat& reynoldsDeriv);

    // ADD NOISE
    void applyNoise(double noiseIntensity);

    // FILTER DATA
    void applyMedianFilter(int qtyID,int maxIt,int order,int filterType,MRIThresholdCriteria* threshold);
    
    // File List Printing
    void printSequenceFiles(std::string outFIleName);
    // Operations on single Scans
    void makeScanDifference(int firstScanID, int secondScanID);
    void makeScanAverage(int numberOfMeasures, int firstScanID, int secondScanID);
    // Eval Time Derivatives
    void evalTimeDerivs(int currentScan, int currentCell,double* timeDeriv);
    void evalScanTimeDerivs(int currentScan,MRIDoubleMat& timeDeriv);
  
    // STATISTICS
    void evalScanDifferencePDF(int otherScan, int refScan, const int pdfQuantity, int numberOfBins, bool useBox, MRIDoubleVec& limitBox, MRIDoubleVec& binCenters, MRIDoubleVec& binArray);
    void extractSinglePointTimeCurve(int cellNumber, int exportQty, string fileName);

    // TRANFORMATION
    void crop(double* limitBox);
    void scaleVelocities(double factor);
    void scalePositions(double factor);

    // MESSAGE PASSING
    void distributeSequenceData(MRICommunicator* comm);

    // CLEAN COMPONENTS ON BOUNDARY
    void cleanNormalComponentOnBoundary();
    void interpolateBoundaryVelocities();

    // TEMPLATE FLOW SEQUENCE
    void createSampleCase(MRISamples sampleType,const MRIDoubleVec& params);
};

#endif // MRISEQUENCE_H
