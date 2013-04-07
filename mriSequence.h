#ifndef MRISEQUENCE_H
#define MRISEQUENCE_H

#include "mriScan.h"

class MRISequence
{
  protected:
    int totalScans;
    std::vector<MRIScan*> sequence;
    std::string* fileNames;
    bool isCyclic;
  public:
    // Constructor and Destructor
	  MRISequence(bool cyclic);
	  ~MRISequence();
    // ----------------
    // MEMBER FUNCTIONS
    //-----------------
    
    // ADD AND GET FROM SEQUENCE
    void AddScan(MRIScan* scan);
    MRIScan* GetScan(int scanNumber);
    
    // READ SEQUENCE FROM FILE
    void ReadFromVolSequence(std::string outfileName);
    
    // EXPORT SEQUENCE TO FILE
    void ExportToTECPLOT(std::string outfileName);
    void ExportToVOL(std::string outfileName);
    
    // DIV FREE Filtering
    void ApplyMPFilter(MRIOptions Options, bool useBCFilter, bool useConstantPatterns, MRIThresholdCriteria thresholdCriteria);
    
    // APPLY THRESHOLDING 
    void ApplyThresholding(MRIThresholdCriteria thresholdCriteria);
    
    // PRESSURE COMPUTATION
    void ComputePressureGradients();
    void ComputeRelativePressure(bool doPressureSmoothing);
    
    // File List Printing
    void PrintSequenceFiles(std::string outFIleName);
    // Operations on single Scans
    void MakeScanDifference(MRIScan &firstScan, const MRIScan &secondScan);
    void MakeScanAverage(int numberOfMeasures, MRIScan &firstScan, const MRIScan &secondScan);
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
