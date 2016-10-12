#ifndef MRISEQUENCE_H
#define MRISEQUENCE_H

#include "mriScan.h"
#include "mriCommunicator.h"

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
    
    // ================
    // MEMBER FUNCTIONS
    // ================
    // ADD AND GET FROM SEQUENCE
    void AddScan(MRIScan* scan);
    int GetTotalScans(){return totalScans;}
    bool GetIsCyclic(){return isCyclic;}
    MRIScan* GetScan(int scanNumber);

    // READ SEQUENCE FROM FILE
    void ReadFromVolSequence(std::string outfileName);
    
    // EXPORT SEQUENCE TO FILE
    void ExportToTECPLOT(std::string outfileName);
    void ExportToVOL(std::string outfileName);
    void ExportToVTK(std::string outfileName,MRIThresholdCriteria* thresholdCriteria);
    void WriteExpansionFile(string fileName);
    void ExportForDistancing(string inputFileName, MRIThresholdCriteria* threshold);
    void ExportForPoisson(string inputFileName,double density,double viscosity,MRIThresholdCriteria* threshold,
                          bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                          bool readMuTFromFile, string muTFile, double smagorinskyCoeff);

    // SAVE QUANTITIES TO OUTPUTS
    void   saveVelocity();
    double getVelocityNormAtCell(int cell);

    // TOPOLOGY AND MAPPING
    void CreateTopology();
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
    void ApplySMPFilter(MRIOptions* options, bool isBC, MRICommunicator* comm);
    
    // APPLY THRESHOLDING 
    void ApplyThresholding(MRIThresholdCriteria* thresholdCriteria);

    // EVAL VORTEX CRITERIA
    void EvalVortexCriteria(MRIThresholdCriteria* thresholdCriteria);
    void EvalVorticity(MRIThresholdCriteria* thresholdCriteria);
    void EvalEnstrophy(MRIThresholdCriteria* thresholdCriteria);
    void EvalSMPVortexCriteria();
    
    // PRESSURE COMPUTATION
    void ComputePressureGradients(MRIThresholdCriteria* threshold);
    void ComputeRelativePressure(bool doPressureSmoothing);

    // COMPUTING REYNOLDS STRESSES
    void EvalReynoldsStresses(MRIThresholdCriteria* threshold);
    void EvalScanReynoldsStressDerivs(int currentScan,MRIDoubleMat& reynoldsDeriv);

    // ADD NOISE
    void applyNoise(double noiseIntensity);

    // FILTER DATA
    void ApplyMedianFilter(int qtyID,int maxIt,int order,int filterType,MRIThresholdCriteria* threshold);
    
    // File List Printing
    void PrintSequenceFiles(std::string outFIleName);
    // Operations on single Scans
    void MakeScanDifference(int firstScanID, int secondScanID);
    void MakeScanAverage(int numberOfMeasures, int firstScanID, int secondScanID);
    // Eval Time Derivatives
    void EvalTimeDerivs(int currentScan, int currentCell,double* timeDeriv);
    void EvalScanTimeDerivs(int currentScan,MRIDoubleMat& timeDeriv);
  
    // STATISTICS
    void ExtractSinglePointTimeCurve(int cellNumber, int exportQty, std::string fileName);

    // TRANFORMATION
    void Crop(double* limitBox);
    void ScaleVelocities(double factor);
    void ScalePositions(double factor);

    // MESSAGE PASSING
    void DistributeSequenceData(MRICommunicator* comm);

    // CLEAN COMPONENTS ON BOUNDARY
    void cleanNormalComponentOnBoundary();
    void InterpolateBoundaryVelocities();
};

#endif // MRISEQUENCE_H
