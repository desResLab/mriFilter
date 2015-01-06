#ifndef MRISEQUENCE_H
#define MRISEQUENCE_H

#include "mriScan.h"
#include "mriCommunicator.h"

class MRISequence
{
  protected:
    int totalScans;
    std::vector<MRIScan*> sequence;
    std::string* fileNames;
    bool isCyclic;
  public:
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
    MRIScan* GetScan(int scanNumber);

    // READ SEQUENCE FROM FILE
    void ReadFromVolSequence(std::string outfileName);
    
    // EXPORT SEQUENCE TO FILE
    void ExportToTECPLOT(std::string outfileName);
    void ExportToVOL(std::string outfileName);
    void ExportToVTK(std::string outfileName);
    void ExportForPOISSON(std::string outfileName);
    void WriteExpansionFile(string fileName);

    // SAVE QUANTITIES TO OUTPUTS
    void saveVelocity();
    
    // DIV FREE Filtering
    void ApplySMPFilter(MRIOptions* options, bool isBC, MRICommunicator* comm);
    
    // APPLY THRESHOLDING 
    void ApplyThresholding(MRIThresholdCriteria* thresholdCriteria);

    // EVAL VORTEX CRITERIA
    void EvalVortexCriteria();
    void EvalVorticity();
    void EvalEnstrophy();
    void EvalSMPVortexCriteria();
    
    // PRESSURE COMPUTATION
    void ComputePressureGradients();
    void ComputeRelativePressure(bool doPressureSmoothing);
    
    // File List Printing
    void PrintSequenceFiles(std::string outFIleName);
    // Operations on single Scans
    void MakeScanDifference(int firstScanID, int secondScanID);
    void MakeScanAverage(int numberOfMeasures, int firstScanID, int secondScanID);
    // Eval Time Derivatives
    void EvalTimeDerivs(int currentScan, int currentCell,double* &timeDeriv);
  
    // STATISTICS
    void ExtractSinglePointTimeCurve(int cellNumber, int exportQty, std::string fileName);

    // TRANFORMATION
    void Crop(double* limitBox);
    void ScaleVelocities(double factor);
    void ScalePositions(double factor);
};

#endif // MRISEQUENCE_H
