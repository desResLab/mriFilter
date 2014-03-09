#include <iostream>

#include "mriSequence.h"
#include "mriScan.h"
#include "mriUtils.h"
#include "mriStatistics.h"
#include "mriConstants.h"
#include "schMessages.h"

#include <boost/random.hpp>

// ======================================
// Read Scan Sequence from Raw Data Files
// ======================================
void ConvertDICOMToVTK(std::string inFileName,std::string outfileName){

  // Add File to Sequence
  MRIScan* MyMRIScan = new MRIScan(0.0);

  // Read from Raw Binary File
  MyMRIScan->ReadRAWFileSequence(inFileName);

  // Export to VTK
  MyMRIScan->ExportToVTK("testVTK.vtk");

}

// ============================
// CONVERT TECPLOT ASCII TO VTK
// ============================
void ConvertTECPLOToVTK(std::string inFileName,std::string outfileName){

  // Add File to Sequence
  MRIScan* MyMRIScan = new MRIScan(0.0);

  // Read from Raw Binary File
  MyMRIScan->ReadPltFile(inFileName,true);

  // Export to VTK
  MyMRIScan->ExportToVTK(outfileName.c_str());
  //MyMRIScan->ExportToTECPLOT(outfileName.c_str(),true);

}


// =====================================
// PROCESS SEQUENCE TO PRODUCE PRESSURES
// =====================================
void EvalSequencePressure(std::string inFileName, std::string outfileName){

  // Create New Sequence
  MRISequence* MyMRISequence = new MRISequence(true/*Cyclic Sequence*/);

  // Read From VOL Files
  MyMRISequence->ReadFromVolSequence(inFileName);

  // Apply MP Filter to all Scan Separately
  // SET OPTIONS AND THRESHOLD
  MRIOptions Options(1.0e-4,2000);
  bool useBCFilter = false;
  bool useConstantPatterns = true;
  MRIThresholdCriteria criteria(kCriterionLessThen,kQtyConcentration,500.0);

  // PRELIMINARY THRESHOLDING
  MyMRISequence->ApplyThresholding(criteria);

  // APPLY FULL FILTER
  MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

  // Apply Boundary Condition Filter to all scans
  useBCFilter = true;
  // APPLY FILTER
  MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

  // Apply Final Threshold
  MyMRISequence->ApplyThresholding(criteria);

  // Compute Pressure Gradient
  MyMRISequence->ComputePressureGradients();

  // Compute Relative Pressure
  MyMRISequence->ComputeRelativePressure(false);

  // Export Sequence To VOL File
  MyMRISequence->ExportToVOL(outfileName);
}

// ====================================
// PROCESS SIGNATURE FLOW FOR PRESSURES
// ====================================
void EvalPressureFromSignatureFlow(std::string inFileName,std::string outfileName){

  // Set Parameters
  int totalSlides = 1;
  double currentTime = 0.0;
  double timeIncrement = 0.1;

  // Create New Sequence
  MRISequence* MyMRISequence = new MRISequence(true/*Cyclic Sequence*/);

  // Fill Seguence
  currentTime = 0.0;
  for(int loopA=0;loopA<totalSlides;loopA++){
    MRIScan* MyMRIScan = new MRIScan(currentTime);
    MyMRIScan->CreateSampleCase(kStagnationFlow,20,80,80,0.1,0.1,0.1,currentTime,kdirX);
    //MyMRIScan->CreateSampleCase(kTransientFlow,40,60,80,0.01,0.01,0.01,currentTime,kdirX);
    MyMRISequence->AddScan(MyMRIScan);
    currentTime += timeIncrement;
  }

  // Read Plt File
  //MyMRIScan->ReadPltFile(inFileName,true);

  // Generate Sample Case
  // kPoiseilleFlow, kStagnationFlow, kCylindricalVortex,kSphericalVortex
  //MyMRIScan->CreateSampleCase(kCylindricalVortex,20,50,50,1.0,1.0,1.0,kdirX);
  //MyMRIScan->CreateSampleCase(kSphericalVortex,60,60,60,0.5,0.5,0.5,kdirX);
  //MyMRIScan->CreateSampleCase(kToroidalVortex,40,60,80,1.0,1.0,1.0,kdirX);
  //MyMRIScan->CreateSampleCase(kTransientFlow,40,60,80,0.01,0.01,0.01,0.0,kdirX);

  // TEST: USE VOL FILES
  //MyMRIScan->ReadScanFromVOLFiles(std::string("volDataAn.vol"),std::string("volDataAn.vol"),std::string("volDataAn.vol"),std::string("volDataAn.vol"));

  // USE DIVERGENCE SMOOTHING FILTER
  //MyMRIScan->ApplySmoothingFilter();

  // SET OPTIONS AND THRESHOLD
  MRIOptions Options(1.0e-4,2000);
  bool useBCFilter = false;
  bool useConstantPatterns = true;
  MRIThresholdCriteria criteria(kCriterionLessThen,kQtyConcentration,1000.0);

  // APPLY FILTER
  //MyMRIScan->PerformPhysicsFiltering(Options,useBCFilter,useConstantPatterns,criteria);

  // Update Velocities
  //MyMRIScan->UpdateVelocities();

  // =========================
  // Compute Pressure Gradient
  // =========================
  MyMRISequence->ComputePressureGradients();

  // =========================
  // Compute Relative Pressure
  // =========================
  MyMRISequence->ComputeRelativePressure(false);

  // =========================
  // Select the resulting File
  // =========================
  MyMRISequence->ExportToTECPLOT(outfileName);
}

// ===================
// PROCESS SINGLE SCAN
// ===================
void ProcessSingleScan(std::string inFileName,std::string outfileName,double itTol, int maxIt, std::string thresholdTypeString,double thresholdValue){

  // Create New Sequence
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // Add File to Sequence
  MRIScan* MyMRIScan = new MRIScan(0.0);
  MyMRIScan->ReadPltFile(inFileName, true);
  MyMRISequence->AddScan(MyMRIScan);

  // Echo Inserted Parameters
  WriteSchMessage(std::string("--------------------------------------------\n"));
  WriteSchMessage(std::string("MP Iteration tolerance Value: ")+MRIUtils::FloatToStr(itTol)+"\n");
  WriteSchMessage(std::string("Threshold Value             : ")+MRIUtils::FloatToStr(thresholdValue)+"\n");
  WriteSchMessage(std::string("--------------------------------------------\n"));

  // SET OPTIONS AND THRESHOLD
  MRIOptions Options(itTol,maxIt);
  bool useBCFilter = false;
  bool useConstantPatterns = true;
  int thresholdType = 0;
  if (thresholdTypeString == "conc"){
    thresholdType = kQtyConcentration;
  }else{
    thresholdType = kQtyVelModule;
  }
  MRIThresholdCriteria criteria(kCriterionABSLessThen,thresholdType,thresholdValue);

  // APPLY FULL FILTER
  MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

  // APPLY BOUNDARY CONDITION FILTER
  useBCFilter = true;
  MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

  // SAVE EXPANSION TO FILE
  MyMRISequence->GetScan(0)->WriteExpansionFile(std::string("ExpansionFile.dat"));

  // APPLY FINAL THRESHOLD
  MyMRISequence->ApplyThresholding(criteria);

  // Evaluate Statistics
  //MyMRISequence->EvalStatistics();

  // Compute Pressure Gradient
  //MyMRISequence->ComputePressureGradients();

  // Compute Relative Pressure
  //MyMRISequence->ComputeRelativePressure(false);

  // EXPORT FILE
  MyMRISequence->ExportToTECPLOT(outfileName);
}

// =============================
// EVAL STATISTICS BETWEEN SCANS
// =============================
void ComputeScanStatistics(std::string firstFileName,std::string secondFileName,std::string statFileName){

    // Init File Names
    std::string statFileNameFirst = statFileName+"_First.dat";
    std::string statFileNameSecond = statFileName+"_Second.dat";
    std::string statFileNameDiff = statFileName+"_Diff.dat";

    // Read the File in a Sequence
    // Create New Sequence
    MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

    // Add First File to Sequence
    MRIScan* MyMRIScan = new MRIScan(0.0);
    MyMRIScan->ReadPltFile(firstFileName, true);
    MyMRISequence->AddScan(MyMRIScan);

    // Add Second File to Sequence
    MyMRIScan = new MRIScan(1.0);
    MyMRIScan->ReadPltFile(secondFileName, true);
    MyMRISequence->AddScan(MyMRIScan);

    // Create Statistics
    bool useBox = false;
    int numberOfBins = 301;
    double* limitBox = new double[6];
    // Set Limits
    limitBox[0] = MyMRISequence->GetScan(0)->domainSizeMin[0];
    limitBox[1] = MyMRISequence->GetScan(0)->domainSizeMax[0];
    limitBox[2] = MyMRISequence->GetScan(0)->domainSizeMin[1];
    limitBox[3] = MyMRISequence->GetScan(0)->domainSizeMax[1];
    limitBox[4] = MyMRISequence->GetScan(0)->domainSizeMin[2];
    limitBox[5] = MyMRISequence->GetScan(0)->domainSizeMax[2];
    // Apply Factors to limitBox
    double xFactor = 1.0;
    double yFactor = 1.0;
    double zFactor = 1.0;
    MRIUtils::ApplylimitBoxFactors(xFactor,yFactor,zFactor,limitBox);
    // Allocate Bin Arrays
    double* binCenters = new double[numberOfBins];
    double* binValues =  new double[numberOfBins];
    // Eval Single PDFs
    // FIRST
    MRIStatistics::EvalScanPDF(MyMRISequence->GetScan(0),kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
    MRIUtils::PrintBinArrayToFile(statFileNameFirst,numberOfBins,binCenters,binValues);

    // SECOND
    MRIStatistics::EvalScanPDF(MyMRISequence->GetScan(1),kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
    MRIUtils::PrintBinArrayToFile(statFileNameSecond,numberOfBins,binCenters,binValues);

    // DIFFERENCE
    MRIStatistics::EvalScanDifferencePDF(MyMRISequence,1,0,kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
    MRIUtils::PrintBinArrayToFile(statFileNameDiff,numberOfBins,binCenters,binValues);

    // DEALLOCATE
    delete [] binCenters;
    delete [] binValues;
}

// =============================
// EXPLICITLY EVAL SCAN MATRICES
// =============================
void ComputeScanMatrices(){
  // VAR
  int totalERows = 0;
  int totalECols = 0;
  double** EMat = NULL;
  int totalDRows = 0;
  int totalDCols = 0;
  double** DMat = NULL;
  int totalStarRows = 0;
  int totalStarCols = 0;
  double** StarMatrix = NULL;

  // Print Matrices for Matlab Analysis Of Variance
  MRIScan* MyMRIScan = new MRIScan(0.0);
  // ISOTROPIC
  MyMRIScan->CreateSampleCase(kConstantFlow,5,5,5,1.0,1.0,1.0,0.0,kdirX);
  // ANISOTROPIC
  //MyMRIScan->CreateSampleCase(kConstantFlow,5,5,5,1.0,2.0,3.0,0.0,kdirX);
  // Assemble All Matrices
  // Assemble Encoding
  MyMRIScan->AssembleEncodingMatrix(totalERows,totalECols,EMat);
  MRIUtils::PrintMatrixToFile("EncodingMat.dat",totalERows,totalECols,EMat);
  // Assemble Decoding
  MyMRIScan->AssembleDecodingMatrix(totalDRows,totalDCols,DMat);
  MRIUtils::PrintMatrixToFile("DecodingMat.dat",totalDRows,totalDCols,DMat);
  // Assemble Star Matrix
  MyMRIScan->AssembleStarMatrix(totalStarRows,totalStarCols,StarMatrix);
  MRIUtils::PrintMatrixToFile("StarMat.dat",totalStarRows,totalStarCols,StarMatrix);

  // Deallocate Matrices
  for(int loopA=0;loopA<totalERows;loopA++){
    delete [] EMat[loopA];
  }
  for(int loopA=0;loopA<totalDRows;loopA++){
    delete [] DMat[loopA];
  }
  for(int loopA=0;loopA<totalStarRows;loopA++){
    delete [] StarMatrix[loopA];
  }
  delete [] EMat;
  delete [] DMat;
  delete [] StarMatrix;
}

// ================
// READ FACE FLUXES
// ================
void ShowFaceFluxPatterns(std::string faceFluxFileName, std::string outFileName){
  // NEW SCAN
  MRIScan* MyMRIScan = new MRIScan(0.0);

  // ISOTROPIC
  MyMRIScan->CreateSampleCase(kConstantFlow,5,5,5,1.0,1.0,1.0,0.0,kdirX);

  // ANISOTROPIC
  //MyMRIScan->CreateSampleCase(kConstantFlow,5,5,5,1.0,2.0,3.0,0.0,kdirX);

  // READ FACE FLUXES FROM FILE
  int totalRows = 0;
  int totalCols = 0;
  std::vector<std::vector<double> > faceFluxMat;
  MRIUtils::ReadMatrixFromFile(faceFluxFileName,totalRows,totalCols,faceFluxMat);

  // COPY THE INTERESTING COLUMN
  double faceFluxVec[totalRows];
  for(int loop0=0;loop0<100;loop0++){
    int selectedCol = loop0;
    for(int loopA=0;loopA<totalRows;loopA++){
      faceFluxVec[loopA] = faceFluxMat[loopA][selectedCol];
    }

    // TRASFORM FACE FLUXES IN VELOCITIES
    MyMRIScan->RecoverCellVelocitiesRT0(false,faceFluxVec);
    // UPDATE VELOCITIES
    MyMRIScan->UpdateVelocities();

    // EXPORT TO VTK
    MyMRIScan->ExportToVTK(outFileName + "_" + MRIUtils::IntToStr(loop0) + ".vtk");
  }
}

// ===========================
// Test Expansion Coefficients
// ===========================
void TEST_ExpansionCoefficients(std::string inFileName){

  // Create New Sequence
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // Add File to Sequence
  MRIScan* MyMRIScan = new MRIScan(0.0);
  MyMRIScan->ReadPltFile(inFileName, true);
  MyMRISequence->AddScan(MyMRIScan);

  // SET OPTIONS AND THRESHOLD
  MRIOptions Options(1.5e-2,2000);
  bool useBCFilter = false;
  bool useConstantPatterns = true;
  int thresholdType = 0;
  std::string thresholdTypeString("Conc");
  double thresholdValue = 0.0;
  if (thresholdTypeString == "Conc"){
    thresholdType = kQtyConcentration;
  }else{
    thresholdType = kQtyVelModule;
  }
  MRIThresholdCriteria criteria(kCriterionLessThen,thresholdType,thresholdValue);

  // APPLY FULL FILTER
  MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

  // GET FITRST SCAN
  MRIScan* firstMRIScan = new MRIScan(0.0);
  firstMRIScan = MyMRISequence->GetScan(0);

  // THRESHOLDING ON EXPANSION COEFFICIENTS
  firstMRIScan->expansion->ApplyVortexThreshold(kSoftThreshold,0.5);
  //firstMRIScan->expansion->WriteToFile("Expansion.dat");

  // REBUILD FROM EXPANSION
  MRISequence* ReconstructedSequence = new MRISequence(MyMRISequence);
  MRIScan* currScan = NULL;
  for(int loopA=0;loopA<ReconstructedSequence->GetTotalScans();loopA++){
    currScan = ReconstructedSequence->GetScan(loopA);
    currScan->RebuildFromExpansion(firstMRIScan->expansion,false);
  }

  // COMPARE THE TWO SCANS
  double currDiffNorm = 0.0;
  MRIScan* currScan1 = NULL;
  MRIScan* currScan2 = NULL;
  for(int loopA=0;loopA<MyMRISequence->GetTotalScans();loopA++){
    currScan1 = MyMRISequence->GetScan(loopA);
    currScan2 = ReconstructedSequence->GetScan(loopA);
    // Get Difference Norm
    currDiffNorm += currScan1->GetDiffNorm(currScan2);
  }

  // WRITE OUTPUT FILES TO VTK
  MyMRISequence->ExportToVTK("FilteredSeq.vtk");
  ReconstructedSequence->ExportToVTK("ReconstructedSeq.vtk");

  // PRINT DIFFERENCE NORM
  WriteSchMessage("Difference Norm " + MRIUtils::FloatToStr(currDiffNorm) + "\n");
}

// ==================
// PRINT THRESHOLDING
// ==================
void TEST02_PrintThresholdingToVTK(std::string inFileName){

  // Create New Sequence
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);
  MRISequence* MyRECSequence = new MRISequence(false/*Cyclic Sequence*/);

  // Add File to Sequence
  MRIScan* MyMRIScan = new MRIScan(0.0);
  MyMRIScan->ReadPltFile(inFileName, true);
  MyMRISequence->AddScan(MyMRIScan);

  // SET OPTIONS AND THRESHOLD
  MRIOptions Options(1.0e-4,2000);
  bool useBCFilter = false;
  bool useConstantPatterns = true;
  int thresholdType = 0;
  std::string thresholdTypeString("Conc");
  double thresholdValue = 0.0;
  if (thresholdTypeString == "Conc"){
    thresholdType = kQtyConcentration;
  }else{
    thresholdType = kQtyVelModule;
  }
  MRIThresholdCriteria criteria(kCriterionLessThen,thresholdType,thresholdValue);

  // APPLY FULL FILTER
  MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

  // APPLY SUCCESSIVE THRESHOLD
  for(int loopA=0;loopA<5;loopA++){
    WriteSchMessage("Applying Threshold " + MRIUtils::IntToStr(loopA) + "\n");
    // CREATE NEW EXPANSION
    MRIExpansion* currExp = new MRIExpansion(MyMRISequence->GetScan(0)->expansion);
    // APPLY THRESHOLD
    currExp->ApplyVortexThreshold(kSoftThreshold,(0.95/4.0)*(loopA));
    // WRITE EXPANSION TO FILE
    currExp->WriteToFile("Expansion_" + MRIUtils::IntToStr(loopA) + "\n");
    // GET A SCAN FROM ORIGINAL SEQUENCE
    MRIScan* myScan = new MRIScan(MyMRISequence->GetScan(0));
    myScan->RebuildFromExpansion(currExp,false);
    // ADD SCAN TO RECONSTRUCTION
    MyRECSequence->AddScan(myScan);
  }

  // WRITE OUTPUT FILES TO VTK
  MyMRISequence->ExportToVTK("FilteredSeq");
  MyRECSequence->ExportToVTK("ReconstructedSeq");
}

// ======================
// EVAL REYNOLDS STRESSES
// ======================
void TEST03_EvalReynoldsStresses(std::string inFileName){

  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);
  MRISequence* MyRECSequence = new MRISequence(false/*Cyclic Sequence*/);

  // ADD FILE TO SEQUENCE
  MRIScan* MyMRIScan = new MRIScan(0.0);
  MyMRIScan->ReadPltFile(inFileName, true);
  MyMRISequence->AddScan(MyMRIScan);

  // SET OPTIONS AND THRESHOLD
  MRIOptions Options(5.0e-1,2000);
  bool useBCFilter = false;
  bool useConstantPatterns = true;
  int thresholdType = 0;
  std::string thresholdTypeString("Conc");
  double thresholdValue = 0.0;
  if (thresholdTypeString == "Conc"){
    thresholdType = kQtyConcentration;
  }else{
    thresholdType = kQtyVelModule;
  }
  MRIThresholdCriteria criteria(kCriterionLessThen,thresholdType,thresholdValue);

  // APPLY FULL FILTER
  MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

  WriteSchMessage("Filter Applied!\n");

  // APPLY SUCCESSIVE THRESHOLD
  for(int loopA=0;loopA<2;loopA++){

    // WRITE MESSAGE
    WriteSchMessage("Applying Threshold " + MRIUtils::IntToStr(loopA) + "\n");

    // CREATE NEW EXPANSION
    MRIExpansion* currExp = new MRIExpansion(MyMRISequence->GetScan(0)->expansion);

    // APPLY THRESHOLD
    currExp->ApplyVortexThreshold(kHardThresold,(0.99/1.0)*(loopA));

    // WRITE EXPANSION TO FILE
    currExp->WriteToFile("Expansion_" + MRIUtils::IntToStr(loopA) + "\n");

    // GET A SCAN FROM ORIGINAL SEQUENCE
    MRIScan* myScan = new MRIScan(MyMRISequence->GetScan(0));
    myScan->RebuildFromExpansion(currExp,false);

    // EVAL REYNOLDS STRESSES
    WriteSchMessage("Evaluating Reynolds Stresses...");
    myScan->EvalReynoldsStressComponent();
    WriteSchMessage("Done.\n");

    // ADD SCAN TO RECONSTRUCTION
    MyRECSequence->AddScan(myScan);
  }

  // WRITE OUTPUT FILES TO VTK
  MyRECSequence->ExportToVTK("ReynoldsStressSequence");
}

// ===========================================
// EVAL REYNOLDS STRESSES FROM EXPANSION FILES
// ===========================================
void EvalPressureFromExpansion(std::string inFileName,std::string outFileName,bool exportTECPLOT){

  // SET PARAMETERS
  bool applyThreshold = true;
  double thresholdRatio = 0.0;
  bool doPressureSmoothing = true;

  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // ADD FILE TO SEQUENCE
  MRIScan* MyMRIScan = new MRIScan(0.0);
  //MyMRIScan->ReadFromExpansionFile(inFileName,applyThreshold,thresholdRatio);
  MyMRIScan->ReadPltFile(inFileName,true);
  MyMRISequence->AddScan(MyMRIScan);
  // APPLY MEDIAN FILTER TO VELOCITIES
  //MyMRISequence->GetScan(0)->ApplyMedianFilter(kQtyVelocityX,1);
  //MyMRISequence->GetScan(0)->ApplyMedianFilter(kQtyVelocityY,1);
  //MyMRISequence->GetScan(0)->ApplyMedianFilter(kQtyVelocityZ,1);

  //MyMRISequence->GetScan(0)->ThresholdQuantity(kQtyVelocityZ,1.0e10);

  // EVAL REYNOLDS STRESSES AND PRESSURE GRADIENTS
  MyMRISequence->ComputePressureGradients();

  // APPLY MEDIAN FILTER TO PRESSURE GRADIENT COMPONENTS
  // JET FILIPPO
  //MyMRISequence->GetScan(0)->ThresholdQuantity(kQtyPressGradientMod,2500.0);
  //MyMRISequence->GetScan(0)->ThresholdQuantity(kQtyPressGradientMod,1200.0);
  //MyMRISequence->GetScan(0)->ThresholdQuantity(kQtyPressGradientMod,200.0);

  //MyMRISequence->GetScan(0)->ThresholdQuantity(kQtyPressGradientMod,50000.0);

  //MyMRISequence->GetScan(0)->ApplyMedianFilter(kQtyPressGradientX,5);
  //MyMRISequence->GetScan(0)->ApplyMedianFilter(kQtyPressGradientY,5);
  //MyMRISequence->GetScan(0)->ApplyMedianFilter(kQtyPressGradientZ,5);

  // EVAL LOCATIONS OF NOISY POINTS
  MyMRISequence->GetScan(0)->EvalNoisyPressureGradientPoints();

  // EVAL RELATIVE PRESSURE
  MyMRISequence->ComputeRelativePressure(doPressureSmoothing);

  if(exportTECPLOT){
    // WRITE OUTPUT FILES TO TECPLOT
    MyMRISequence->ExportToTECPLOT(outFileName);
  }else{
    // WRITE OUTPUT FILES TO VTK
    MyMRISequence->ExportToVTK(outFileName);
  }
}

// =========================================
// EVAL PRESSURES FROM EXPANSION COFFICIENTS
// =========================================
void EvalConcentrationGradient(std::string inFileName,std::string outFileName){

  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  bool doPressureSmoothing = false;

  // ADD FILE TO SEQUENCE
  MRIScan* MyMRIScan = new MRIScan(0.0);
  MyMRIScan->ReadPltFile(inFileName,true);
  MyMRIScan->ScalePositions(0.0058);
  MyMRISequence->AddScan(MyMRIScan);

  // EVAL REYNOLDS STRESSES AND PRESSURE GRADIENTS
  MyMRISequence->GetScan(0)->ComputeQuantityGradient(kQtyConcentration);

  // EVAL RELATIVE PRESSURE
  MyMRISequence->ComputeRelativePressure(doPressureSmoothing);

  // WRITE OUTPUT FILES TO VTK
  MyMRISequence->ExportToVTK(outFileName);

}



// ===================
// PERFORM RANDOM TEST
// ===================
void PerformRandomTest(){
  // Set parameters
  int numberMC = 500;
  // SET OPTIONS AND THRESHOLD
  MRIOptions Options(1.0e-10,2000);
  bool useBCFilter = false;
  bool useConstantPatterns = true;
  MRIThresholdCriteria criteria(kCriterionLessThen,kQtyConcentration,1.0);

  // Generating Random Numbers
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator(boost::mt19937(time(0)),boost::normal_distribution<>());

  // Declare Scan
  MRIScan* MyMRIScan;

  // Monte Carlo Loops
  for(int loopA=0;loopA<numberMC;loopA++){

    // Generate Random Velocity Field - Standard Gaussian Distribution
    MyMRIScan = new MRIScan(0.0);
    MyMRIScan->CreateSampleCase(kZeroVelocity,5,5,5,1.0,1.0,1.0,0.0,kdirX);

    // Assign Random Component
    MyMRIScan->AssignRandomComponent(kdirX,generator);
    MyMRIScan->AssignRandomComponent(kdirY,generator);
    MyMRIScan->AssignRandomComponent(kdirZ,generator);

    // Write Generated Velocities To File
    MyMRIScan->ExportVelocitiesToFile("InputVelocities.dat",true);

    // Filter Velocities - NO GLOBAL ATOMS
    MyMRIScan->PerformPhysicsFiltering(Options, false, false, criteria);
    MyMRIScan->UpdateVelocities();

    // Write Resultant Velocities
    MyMRIScan->ExportVelocitiesToFile("OutputVelocities.dat",true);

    // Deallocate
    delete MyMRIScan;
    MyMRIScan = NULL;
  }
}

// ================================
// CROP DICOM AND EVALUATE PRESSURE
// ================================
void CropAndComputeVol(std::string inFileName,std::string outfileName){
  // ================
  // ANALYZE SEQUENCE
  // ================

  // -------------------
  // Create New Sequence
  // -------------------
  MRISequence* MyMRISequence = new MRISequence(true/*Cyclic Sequence*/);

  // -------------------
  // Read From VOL Files
  // -------------------
  MyMRISequence->ReadFromVolSequence(inFileName);

  // ----------------
  // Create Limit Box
  // ----------------
  double limitBox[6] = {0.0};
  limitBox[0] = 135.0;
  limitBox[1] = 160.0;
  limitBox[2] = 130.0;
  limitBox[3] = 165.0;
  limitBox[4] = MyMRISequence->GetScan(0)->domainSizeMin[2];
  limitBox[5] = MyMRISequence->GetScan(0)->domainSizeMax[2];

  // ---------------
  // Reduce Sequence
  // ---------------
  MyMRISequence->Crop(limitBox);

  // ----------------
  // Scale Velocities
  // ----------------
  MyMRISequence->ScaleVelocities(0.001);
  MyMRISequence->ScalePositions(0.001);

  // SET OPTIONS AND THRESHOLD
  MRIOptions Options(1.0e-4,2000);
  bool useBCFilter = false;
  bool useConstantPatterns = true;
  MRIThresholdCriteria criteria(kCriterionLessThen,kQtyConcentration,500.0);

  // ------------------------
  // PRELIMINARY THRESHOLDING
  // ------------------------
  MyMRISequence->ApplyThresholding(criteria);

  // -----------------
  // APPLY FULL FILTER
  // -----------------
  MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

  // --------------------------------------------
  // Apply Boundary Condition Filter to all scans
  // --------------------------------------------
  useBCFilter = true;
  // APPLY FILTER
  MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

  // ---------------------
  // Apply Final Threshold
  // ---------------------
  MyMRISequence->ApplyThresholding(criteria);

  // -------------------------
  // Compute Pressure Gradient
  // -------------------------
  MyMRISequence->ComputePressureGradients();

  // -------------------------
  // Compute Relative Pressure
  // -------------------------
  MyMRISequence->ComputeRelativePressure(false);

  // ---------------------------
  // Export Sequence To VOL File
  // ---------------------------
  MyMRISequence->ExportToTECPLOT(outfileName.c_str());
}

// =========================
// PERFORM STREAMLINE TEST 1
// =========================
void PerformStreamlineTest1(int intValue,std::string inFileName,std::string outfileName){

  // Set Default Options For Streamlines
  // Re   Xmin  Xmax Ymin  Ymax Zmin  Zmax
  // 80   -18.0 15.0 -12.0 21.0 -78.0 75.0
  // 145  -12.0 21.0 -22.0 11.0 -78.0 75.0
  // 190  -14.0 19.0 -3.0  30.0 -78.0 75.0
  // 240  -14.0 19.0 -3.0  30.0 -78.0 75.0
  // 290  -12.0 21.0 -22.0 11.0 -78.0 75.0
  // 458  -12.0 21.0 -22.0 11.0 -78.0 75.0
  // 1390 -15.0 18.0 -12.0 21.0 -78.0 75.0
  // 2740 -15.0 18.0 -12.0 21.0 -78.0 75.0
  MRIStreamlineOptions slOptions(0.0,0.0,0.0,0.0,0.0,0.0);
  switch (intValue){
    case 0:
      // Re 80
      slOptions.setLimits(-18.0,15.0,-12.0,21.0,-78.0,75.0);
      break;
    case 1:
      // Re 145
      slOptions.setLimits(-12.0,21.0,-22.0,11.0,-78.0,75.0);
      break;
    case 2:
      // Re 190
      slOptions.setLimits(-14.0,19.0,-3.0,30.0,-78.0,75.0);
      break;
    case 3:
      // Re 240
      slOptions.setLimits(-14.0,19.0,-3.0,30.0,-78.0,75.0);
      break;
    case 4:
      // Re 290
      slOptions.setLimits(-12.0,21.0,-22.0,11.0,-78.0,75.0);
      break;
    case 5:
      // Re 458
      slOptions.setLimits(-12.0,21.0,-22.0,11.0,-78.0,75.0);
      break;
    case 6:
      // Re 1390
      slOptions.setLimits(-15.0,18.0,-12.0,21.0,-78.0,75.0);
      break;
    case 7:
      // Re 2740
      slOptions.setLimits(-15.0,18.0,-12.0,21.0,-78.0,75.0);
      break;
  }

  // Get First Scan
  MRIScan* myScan = new MRIScan(0.0);

  // Read From File
  //myScan->ReadPltFile(inFileName,true);

  // Export to VTK
  //myScan->ExportToVTK(outfileName+"_Model.vtk");

  // Compute Streamlines
  std::vector<MRIStreamline*> streamlines;
  //myScan->ComputeStreamlines(outfileName+"_grid.dat",slOptions,streamlines);

  // Read Streamlines from File
  MRIUtils::ReadStreamlinesFromLegacyVTK(inFileName,streamlines);

  // Print Streamlines To File
  //MRIUtils::PrintStreamlinesToVTK(streamlines,outfileName+"_DebugStreamlines.vtk");

  // Eval Streamlines Statistics
  myScan->EvalStreamLineStatistics(outfileName,kdirZ,slOptions,streamlines);
}

// =========================
// PERFORM STREAMLINE TEST 2
// =========================
void PerformStreamlineTest2(std::string inFileName,std::string outfileName){

  // Set Default Options For Streamlines
  // Re   Xmin  Xmax Ymin  Ymax Zmin  Zmax
  // 80   -18.0 15.0 -12.0 21.0 -78.0 75.0
  // 145  -12.0 21.0 -22.0 11.0 -78.0 75.0
  // 190  -14.0 19.0 -3.0  30.0 -78.0 75.0
  // 240  -14.0 19.0 -3.0  30.0 -78.0 75.0
  // 290  -12.0 21.0 -22.0 11.0 -78.0 75.0
  // 458  -12.0 21.0 -22.0 11.0 -78.0 75.0
  // 1390 -15.0 18.0 -12.0 21.0 -78.0 75.0
  // 2740 -15.0 18.0 -12.0 21.0 -78.0 75.0
  MRIStreamlineOptions slOptions(0.0,0.0,0.0,0.0,0.0,0.0);
  slOptions.setLimits(-18.0,15.0,-12.0,21.0,-78.0,75.0);

  // Get First Scan
  MRIScan* myScan = new MRIScan(0.0);

  // Read From File
  myScan->ReadPltFile(inFileName,true);

  // Set Velocities constant in Z
  for(int loopA=0;loopA<myScan->totalCellPoints;loopA++){
    myScan->cellPoints[loopA].velocity[0] = 0.01;
    myScan->cellPoints[loopA].velocity[1] = 0.01;
    myScan->cellPoints[loopA].velocity[2] = 1.0;
  }

  // Export to VTK
  myScan->ExportToVTK(outfileName+"_Model.vtk");

  // Compute Streamlines
  std::vector<MRIStreamline*> streamlines;
  myScan->ComputeStreamlines(outfileName+"_grid.dat",slOptions,streamlines);

  // Read Streamlines from File
  //MRIUtils::ReadStreamlinesFromLegacyVTK(inFileName,streamlines);

  // Print Streamlines To File
  MRIUtils::PrintStreamlinesToVTK(streamlines,outfileName+"_DebugStreamlines.vtk");

  // Eval Streamlines Statistics
  myScan->EvalStreamLineStatistics(outfileName,kdirZ,slOptions,streamlines);
}

// ========================================
// BUILD FROM COEFFICIENTS AND WRITE TO PLT
// ========================================
void BuildFromCoeffs(std::string coeffFileName,std::string plotOut,bool performThreshold,int thresholdType,double threshold){

  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // ADD FILE TO SEQUENCE
  MRIScan* MyMRIScan = new MRIScan(0.0);
  MyMRIScan->ReadFromExpansionFile(coeffFileName,performThreshold,thresholdType,threshold);
  MyMRISequence->AddScan(MyMRIScan);

  // APPLY THRESHOLDING
  //MRIThresholdCriteria criteria(kCriterionLessThen,kQtyConcentration,1000.0);
  //MyMRISequence->ApplyThresholding(criteria);

  // EXPORT TO PLT FILE
  //MyMRISequence->ExportToTECPLOT(plotOut);
  MyMRISequence->ExportToVTK(plotOut);
}

// ============================
// EVAL VARIOUS VORTEX CRITERIA
// ============================
void EvalVortexCriteria(std::string inFileName,std::string outFileName){
  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // ADD FILE TO SEQUENCE
  MRIScan* MyMRIScan = new MRIScan(0.0);
  MyMRIScan->ReadPltFile(inFileName,true);
  MyMRISequence->AddScan(MyMRIScan);

  // EVAL VORTEX CRITERIA
  // MyMRISequence->GetScan(0)->EvalVortexCriteria();
  // MyMRISequence->GetScan(0)->EvalVorticity();
  MyMRISequence->GetScan(0)->EvalEnstrophy();


  // WRITE OUTPUT FILES TO VTK
  MyMRISequence->ExportToVTK(outFileName);
}

// ===========================================
// WRITE SPATIAL DISTRIBUTIONS OF COEFFICIENTS
// ===========================================
void WriteSpatialExpansion(std::string expFileName,std::string outFileName){
  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // CREATE NEW SCAN
  MRIScan* MyMRIScan = new MRIScan(0.0);

  // READ FROM EXPANSION COEFFICIENTS
  MyMRIScan->ReadFromExpansionFile(expFileName,false,kSoftThreshold,0.0);

  // ADD TO SEQUENCE
  MyMRISequence->AddScan(MyMRIScan);

  // SPATIALLY EVALUATE VORTEX COEFFICIENTS
  MyMRISequence->GetScan(0)->EvalSMPVortexCriteria(MyMRISequence->GetScan(0)->expansion);

  // EXPORT TO VTK
  MyMRISequence->ExportToVTK(outFileName);
}

// ============
// MAIN PROGRAM
// ============
int main(int argc, char **argv)
{
  // WRITE PROGRAM HEADER
  WriteHeader();
  // CHECK FIRST ARGUMENT IF AVAILABLE
  std::string firstOption("");
  if(argc>1){
    firstOption = argv[1];
  }
  // CHECK VARIOUS OPTIONS FOR THE PROGRAM
  if ((argc == 1)||(firstOption == "-?")||(firstOption == "-help")||(firstOption == "-h")||(firstOption == "")){
    // Write Program Help
    MRIUtils::WriteProgramHelp();
  }else if (firstOption == "-seq"){

    // Get File Names
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);

    // EVAL PRESSURE FROM VELOCITY SEQUENCE
    EvalSequencePressure(inFileName,outfileName);
    
  }else if (firstOption == "-o"){

    // Get File Names
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);

    // EVAL PRESSURE FROM SIGNATURE FLOW
    EvalPressureFromSignatureFlow(inFileName,outfileName);

  }else if (firstOption == "-scan"){

    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);
    // Get Parameters
    double itTol = atof(argv[4]);
    int maxIt = atoi(argv[5]);
    std::string thresholdTypeString(argv[6]);
    double thresholdValue = atof(argv[7]);

    // PROCESS SINGLE SCAN
    ProcessSingleScan(inFileName,outfileName,itTol,maxIt,thresholdTypeString,thresholdValue);

  }else if (firstOption == "-pltToVTK"){

    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);

    // Convert to VTK
    ConvertTECPLOToVTK(inFileName,outfileName);

  }else if (firstOption == "-stat"){

    // Get File Name
    std::string firstFileName(argv[2]);
    std::string secondFileName(argv[3]);
    std::string statFileName(argv[4]);

    // COMPUTE SCAN STATISTICS
    ComputeScanStatistics(firstFileName,secondFileName,statFileName);

    
  }else if (firstOption == "-Mat"){

    // COMPUTE SCAN MATRICES
    ComputeScanMatrices();

  }else if (firstOption == "-randtest"){

    // PERFORM RANDOM TEST
    PerformRandomTest();

  }else if (firstOption == "-reduceVol"){

    // GET FILE NAMES FROM ARGUMENTS
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);

    // CROP AND REDUCE VOL FILE
    CropAndComputeVol(inFileName,outfileName);
    
  }else if (firstOption == "-streamLinesTest1"){

    // GET FILE NAME
    int intValue = atoi(argv[2]);
    std::string inFileName(argv[3]);
    std::string outfileName(argv[4]);

    // PERFORM STREAMLINE TEST 1
    PerformStreamlineTest1(intValue,inFileName,outfileName);
      
  }else if (firstOption == "-streamLinesTest2"){

    // GET FILE NAME
    std::string inFileName(argv[3]);
    std::string outfileName(argv[4]);

    // PERFORM STREAMLINE TEST 2
    PerformStreamlineTest2(inFileName,outfileName);

  }else if (firstOption == "-testExpansion"){

    // Get Input and output File Names
    std::string inFileName(argv[2]);

    // Test Expansion Coefficients
    // TEST_ExpansionCoefficients(inFileName);
    TEST02_PrintThresholdingToVTK(inFileName);


  }else if (firstOption == "-ReynoldsStress"){

    // Get Input and output File Names
    std::string inFileName(argv[2]);

    // Test Expansion Coefficients
    TEST03_EvalReynoldsStresses(inFileName);

  }else if (firstOption == "-readflux"){

    // GET FILE NAMES
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);

    // READ FACE FLUXES FROM FILE AND EXPORT TO VTK
    ShowFaceFluxPatterns(inFileName,outfileName);

  }else if (firstOption == "-buildFromCoeffs"){

    // GET FILE NAMES
    std::string coeffFileName(argv[2]);
    std::string plotOut(argv[3]);
    double threshold = atof(argv[4]);

    // READ FROM COEFFICIENT FILE AND EXPORT TO PLT
    BuildFromCoeffs(coeffFileName,plotOut,true,kSoftThreshold,threshold);

  }else if (firstOption == "-evalPressure"){

    // GET FILE NAMES
    std::string inFileName(argv[2]);
    std::string outFileName(argv[3]);
    int saveAsPLTInt = atoi(argv[4]);
    bool saveAsPLT = (saveAsPLTInt == 0);

    // EVAL PRESSURES FROM EXPANSION COFFICIENTS
    EvalPressureFromExpansion(inFileName,outFileName,saveAsPLT);

  }else if (firstOption == "-evalConcGrad"){

    // GET FILE NAMES
    std::string inFileName(argv[2]);
    std::string outFileName(argv[3]);

    // EVAL PRESSURES FROM EXPANSION COFFICIENTS
    EvalConcentrationGradient(inFileName,outFileName);

  }else if (firstOption == "-vortexCrit"){

      // GET FILE NAMES
      std::string inFileName(argv[2]);
      std::string outFileName(argv[3]);

      // EVAL PRESSURES FROM EXPANSION COFFICIENTS
      EvalVortexCriteria(inFileName,outFileName);

  }else if (firstOption == "-writeSpatialExpansion"){

      // GET FILE NAMES
      std::string inFileName(argv[2]);
      std::string outFileName(argv[3]);

      // EVAL Spatial Expansion Coefficients Distribution
      WriteSpatialExpansion(inFileName,outFileName);

  }else if (firstOption == "-scaleModel"){

    // GET FILE NAMES
    std::string inFileName(argv[2]);
    std::string outFileName(argv[3]);

    // Create New Sequence
    MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

    // Add File to Sequence
    MRIScan* MyMRIScan = new MRIScan(0.0);
    MyMRIScan->ReadPltFile(inFileName, true);
    //MyMRIScan->ScalePositions(0.0058);
    //MyMRIScan->ScaleVelocities(0.5);
    MyMRISequence->AddScan(MyMRIScan);

    double limitBox[6] = {0.0};
    // JET SPERIMENTALE
    //limitBox[0] = 0.0;
    //limitBox[1] = 0.18;
    //limitBox[2] = 0.0213;
    //limitBox[3] = 0.0668;
    //limitBox[4] = 0.00642;
    //limitBox[5] = 0.0499;
    // LUNG - UPPER PART
    //limitBox[0] = -0.025;
    //limitBox[1] = 0.0872;
    //limitBox[2] = -0.0268;
    //limitBox[3] = 0.0434;
    //limitBox[4] = -0.055;
    //limitBox[5] = 0.0546;
    // LUNG - TRACHEA
    //limitBox[0] = 0.0942;
    //limitBox[1] = 0.185;
    //limitBox[2] = -0.00924;
    //limitBox[3] = 0.0294;
    //limitBox[4] = 0.0188;
    //limitBox[5] = 0.0504;
    // LUNG 2 - TRACHEA
    limitBox[0] = -0.0124;
    limitBox[1] = 0.112;
    limitBox[2] = -0.0177;
    limitBox[3] = 0.0259;
    limitBox[4] = -0.00371;
    limitBox[5] = 0.0539;

    // ====
    // CROP
    // ====
    MyMRISequence->Crop(limitBox);

    // EXPORT FILE
    MyMRISequence->ExportToTECPLOT(outFileName);
    //MyMRISequence->ExportToVTK(outFileName);

  }else{

    // INVALID SWITCH
    std::string currMsgs("Error: Invalid Switch. Terminate.\n");
    WriteSchMessage(currMsgs);
    return (1);
  }

  // COMPLETED!
  std::string currMsgs("\nOperation Completed!\n");
  WriteSchMessage(currMsgs);
  return (0);
};

