#ifndef MRISEQUENCE_H
#define MRISEQUENCE_H

#include "mriScan.h"
#include "mriCommunicator.h"
#include "mriThresholdCriteria.h"
#include "mriTopology.h"
#include "mriIO.h"

using namespace std;

class MRIScan;
class MRITopology;
class MRICommunicator;

// Generic Sequence Containing General 
// Unstructured Data Sets

class MRISequence{
  public:
    int              totalScans;
    vector<MRIScan*> sequence;
    string*          fileNames;
    bool             isCyclic;

    // TOPOLOGY COMMON TO EVERY SCAN
    MRITopology* topology;

    // Constructor and
    MRISequence(bool cyclic);
    // Copy Constructor
    MRISequence(MRISequence* copySequence);
    // Destructor
    ~MRISequence();

    // WRITE STATS
    string writeStatistics();
    
    // ADD AND GET FROM SEQUENCE
    void addScan(MRIScan* scan);
    int  getTotalScans(){return totalScans;}
    bool getIsCyclic(){return isCyclic;}
    MRIScan* getScan(int scanNumber);

    // READ SEQUENCE FROM FILE
    void readFromASCIISequence(int asciiInputType, 
                               const MRIStringVec& vtkFileNames, 
                               const MRIDoubleVec& times);    
    void readFromExpansionFiles(const MRIStringVec& fileNames, 
                                const MRIDoubleVec& Times, 
                                bool applyThreshold, 
                                int thresholdType,
                                double thresholdRatio);
    
    // EXPORT SEQUENCE TO FILE
    void exportToTECPLOT(string outfileName);
    void exportToVOL(string outfileName);
    void exportToVTK(string outfileName,MRIThresholdCriteria* thresholdCriteria);    
    void exportForDistancing(string inputFileName, MRIThresholdCriteria* threshold);
    void exportForPoisson(string inputFileName,double density,double viscosity,MRIThresholdCriteria* threshold,
                          bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                          bool readMuTFromFile, string muTFile, double smagorinskyCoeff);
    void writeExpansionFile(string fileName);    

    // SAVE QUANTITIES TO OUTPUTS
    void   saveVelocity();
    double getVelocityNormAtCell(int cell);

    // TOPOLOGY AND MAPPING
    void createTopology();
    void getUnitVector(int CurrentCell, const MRIDoubleVec& GlobalFaceCoords, MRIDoubleVec& myVect);
    void getGlobalCoords(int DimNumber, int SliceNumber, double FaceCoord1, double FaceCoord2, MRIDoubleVec& globalCoords);
    int  getCellNumber(const MRIDoubleVec& coords);
    void getGlobalPermutation(MRIIntVec& GlobalPerm);
    int  getAdjacentFace(int globalNodeNumber, int AdjType);
    void getNeighborVortexes(int cellNumber,int dim,MRIIntVec& idx);
    void getEdgeDirection(int edgeID, double* edgeDirVector);
    void getAuxNodeCoordinates(int nodeNum, MRIDoubleVec& pos);
    void buildCellConnections();
    int  addToFaceConnections(const MRIIntVec& faceIds, vector<vector<mriFace* > >& AuxFirstNodeFaceList);
    int  addToEdgeConnections(const MRIIntVec& edgeIds,vector<vector<mriEdge*> >& AuxFirstNodeEdgeList);
    void buildFaceConnections();
    void buildFaceCells();
    void buildEdgeConnections();
    
    // MAPPING     
    int  mapCoordsToIndex(int i, int j, int k);
    void mapCoordsToPosition(const MRIIntVec& coords, bool addMeshMinima, MRIDoubleVec& pos);
    void mapIndexToAuxNodeCoords(int index, MRIIntVec& intCoords);
    void mapAuxCoordsToPosition(const MRIIntVec& auxCoords, MRIDoubleVec& pos);    
    
    // DIV FREE Filtering
    void applySMPFilter(MRICommunicator* comm, bool isBC, 
                        MRIThresholdCriteria* thresholdCriteria,
                        double itTol,
                        int maxIt,
                        bool useConstantPatterns);
    
    // APPLY THRESHOLDING 
    void applyThresholding(MRIThresholdCriteria* thresholdCriteria);

    // EVAL VORTEX CRITERIA
    void evalVortexCriteria(MRIThresholdCriteria* thresholdCriteria);
    void evalVorticity(MRIThresholdCriteria* thresholdCriteria);
    void evalEnstrophy(MRIThresholdCriteria* thresholdCriteria);
    void evalSMPVortexCriteria();
    
    // PRESSURE COMPUTATION
    void computePressureGradients(MRIThresholdCriteria* threshold);

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
    void evalTimeDerivs(int currentScan, int currentCell, MRIDoubleVec& timeDeriv);
    void evalScanTimeDerivs(int currentScan,MRIDoubleMat& timeDeriv);
  
    // STATISTICS
    void evalScanDifferencePDF(int otherScan, int refScan, const int pdfQuantity, int numberOfBins, bool useBox, MRIDoubleVec& limitBox, MRIDoubleVec& binCenters, MRIDoubleVec& binArray);
    void extractSinglePointTimeCurve(int cellNumber, int exportQty, string fileName);
    void formDifferenceBinLimits(int otherScan, int refScan, 
                                 int pdfQuantity, double& currInterval,
                                 const MRIDoubleVec& limitBox, 
                                 int numberOfBins, 
                                 MRIDoubleVec& binMin, 
                                 MRIDoubleVec& binMax, 
                                 MRIDoubleVec& binCenter);
    
    // TRANFORMATION
    void crop(const MRIDoubleVec& limitBox);
    void scaleVelocities(double factor);
    void scalePositions(const MRIDoubleVec& origin, double factor);

    // MESSAGE PASSING
    void distributeSequenceData(MRICommunicator* comm);

    // CLEAN COMPONENTS ON BOUNDARY
    void cleanNormalComponentOnBoundary(MRIThresholdCriteria* threshold);
    void interpolateBoundaryVelocities(MRIThresholdCriteria* threshold);

    // TEMPLATE FLOW SEQUENCE
    void createSampleCase(MRISamples sampleType,const MRIDoubleVec& params);
};

#endif // MRISEQUENCE_H
