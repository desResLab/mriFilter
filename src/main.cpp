#include <iostream>

#include "mriScan.h"
#include "mriStructuredScan.h"
#include "mriSequence.h"
#include "mriUtils.h"
#include "mriStatistics.h"
#include "mriConstants.h"
#include "mriOptions.h"
#include "schMessages.h"
#include "mpi.h"
#include "mriCommunicator.h"

#include <boost/random.hpp>

// ======================================
// Read Scan Sequence from Raw Data Files
// ======================================
void ConvertDICOMToVTK(std::string inFileName,std::string outfileName){

  // Add File to Sequence
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // Read from Raw Binary File
  MyMRIScan->ReadRAWFileSequence(inFileName);

  // Export to VTK
  MyMRIScan->ExportToVTK("testVTK.vtk",thresholdCriteria);

}

// ============================
// CONVERT TECPLOT ASCII TO VTK
// ============================
void ConvertTECPLOToVTK(MRIOptions* opts){

  // Add File to Sequence
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // Read from Raw Binary File
  MyMRIScan->ReadPltFile(opts->inputFileName,true);

  WriteSchMessage(std::string("ECCOLO\n"));
  WriteSchMessage(MRIUtils::FloatToStr(MyMRIScan->cellPoints[300].velocity[1]) + std::string("\n"));

  // Export to VTK
  MyMRIScan->ExportToVTK(opts->outputFileName,thresholdCriteria);
  //MyMRIScan->ExportToTECPLOT(outfileName.c_str(),true);
}

// ======================
// SOLVE POISSON EQUATION
// ======================
void SolvePoissonEquation(MRICommunicator* comm){
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
  //MyMRIScan->SolvePoissonEquation(comm);
}

// =====================================
// PROCESS SEQUENCE TO PRODUCE PRESSURES
// =====================================
void EvalSequencePressure(MRIOptions* opts, MRICommunicator* comm){

  // Create New Sequence
  MRISequence* MyMRISequence = new MRISequence(true/*Cyclic Sequence*/);

  // Read From VOL Files
  MyMRISequence->ReadFromVolSequence(opts->inputFileName);

  // PRELIMINARY THRESHOLDING
  MyMRISequence->ApplyThresholding(opts->thresholdCriteria);

  // APPLY FULL FILTER
  MyMRISequence->ApplySMPFilter(opts,false,comm);

  // Apply Boundary Condition Filter to all scans
  if(opts->applyBCFilter){
    // APPLY FILTER
    MyMRISequence->ApplySMPFilter(opts,true,comm);
  }

  // Apply Final Threshold
  MyMRISequence->ApplyThresholding(opts->thresholdCriteria);

  // Compute Pressure Gradient
  MyMRISequence->ComputePressureGradients(opts->thresholdCriteria);

  // Compute Relative Pressure
  MyMRISequence->ComputeRelativePressure(false);

  // Export Sequence To VOL File
  MyMRISequence->ExportToVOL(opts->outputFileName);
}

// ====================================
// PROCESS SIGNATURE FLOW FOR PRESSURES
// ====================================
void EvalPressureFromSignatureFlow(MRIOptions* opts){

  // Set Parameters
  int totalSlides = 1;
  double currentTime = 0.0;
  double timeIncrement = 0.1;

  // Create New Sequence
  MRISequence* MyMRISequence = new MRISequence(true/*Cyclic Sequence*/);

  // Fill Seguence
  currentTime = 0.0;
  for(int loopA=0;loopA<totalSlides;loopA++){
    MRIStructuredScan* MyMRIScan = new MRIStructuredScan(currentTime);
    MyMRIScan->CreateSampleCase(kStagnationFlow,opts->templateParams);
    //MyMRIScan->CreateSampleCase(kTransientFlow,opts->templateParams);
    MyMRISequence->AddScan(MyMRIScan);
    currentTime += timeIncrement;
  }

  // USE DIVERGENCE SMOOTHING FILTER
  //MyMRIScan->ApplySmoothingFilter();

  // APPLY FILTER
  //MyMRIScan->PerformPhysicsFiltering(Options,useBCFilter,useConstantPatterns,criteria);

  // Update Velocities
  //MyMRIScan->UpdateVelocities();

  // =========================
  // Compute Pressure Gradient
  // =========================
  MyMRISequence->ComputePressureGradients(opts->thresholdCriteria);

  // =========================
  // Compute Relative Pressure
  // =========================
  MyMRISequence->ComputeRelativePressure(false);

  // =========================
  // Select the resulting File
  // =========================
  MyMRISequence->ExportToTECPLOT(opts->outputFileName);
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
    MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
    MyMRIScan->ReadPltFile(firstFileName, true);
    MyMRISequence->AddScan(MyMRIScan);

    // Add Second File to Sequence
    MyMRIScan = new MRIStructuredScan(1.0);
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

  bool isIsotropic = true;

  // Print Matrices for Matlab Analysis
  // Create Scan
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
  // Set Template Parameters
  vector<double> params;
  params.push_back(5);
  params.push_back(5);
  params.push_back(5);
  params.push_back(1.0);
  params.push_back(1.0);
  params.push_back(1.0);
  params.push_back(0.0);
  params.push_back(1.0);
  if(isIsotropic){
    // ISOTROPIC
    MyMRIScan->CreateSampleCase(kConstantFlow,params);
  }else{
    // ANISOTROPIC
    params[3] = 1.0;
    params[4] = 2.0;
    params[5] = 3.0;
    MyMRIScan->CreateSampleCase(kConstantFlow,params);
  }
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
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);

  bool isIsotropic = true;

  vector<double> params;
  params.push_back(5.0);
  params.push_back(5.0);
  params.push_back(5.0);
  params.push_back(1.0);
  params.push_back(1.0);
  params.push_back(1.0);
  params.push_back(1.0);
  if(isIsotropic){
    // ISOTROPIC
    MyMRIScan->CreateSampleCase(kConstantFlow,params);
  }else{
   // ANISOTROPIC
   params[3] = 1.0;
   params[4] = 2.0;
   params[5] = 3.0;
   MyMRIScan->CreateSampleCase(kConstantFlow,params);
  }

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

    // Init a Threshold with No quantity
    MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

    // EXPORT TO VTK
    MyMRIScan->ExportToVTK(outFileName + "_" + MRIUtils::IntToStr(loop0) + ".vtk",thresholdCriteria);
  }
}

// ===========================
// Test Expansion Coefficients
// ===========================
void TEST_ExpansionCoefficients(MRIOptions* opts, MRICommunicator* comm){

  // Create New Sequence
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // Add File to Sequence
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
  MyMRIScan->ReadPltFile(opts->inputFileName, true);
  MyMRISequence->AddScan(MyMRIScan);

  // APPLY FULL FILTER
  MyMRISequence->ApplySMPFilter(opts,false,comm);

  // GET FIRST SCAN
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

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // WRITE OUTPUT FILES TO VTK
  MyMRISequence->ExportToVTK("FilteredSeq.vtk",thresholdCriteria);
  ReconstructedSequence->ExportToVTK("ReconstructedSeq.vtk",thresholdCriteria);

  // PRINT DIFFERENCE NORM
  WriteSchMessage("Difference Norm " + MRIUtils::FloatToStr(currDiffNorm) + "\n");
}

// ==================
// PRINT THRESHOLDING
// ==================
void TEST02_PrintThresholdingToVTK(MRIOptions* opts, MRICommunicator* comm){

  // Create New Sequence
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);
  MRISequence* MyRECSequence = new MRISequence(false/*Cyclic Sequence*/);

  // Add File to Sequence
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
  MyMRIScan->ReadPltFile(opts->inputFileName, true);
  MyMRISequence->AddScan(MyMRIScan);

  // APPLY FULL FILTER
  MyMRISequence->ApplySMPFilter(opts,false,comm);

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
    MRIScan* myScan = new MRIScan(*MyMRISequence->GetScan(0));
    myScan->RebuildFromExpansion(currExp,false);
    // ADD SCAN TO RECONSTRUCTION
    MyRECSequence->AddScan(myScan);
  }

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // WRITE OUTPUT FILES TO VTK
  MyMRISequence->ExportToVTK("FilteredSeq",thresholdCriteria);
  MyRECSequence->ExportToVTK("ReconstructedSeq",thresholdCriteria);
}

// ======================
// EVAL REYNOLDS STRESSES
// ======================
void TEST03_EvalReynoldsStresses(MRIOptions* opts, MRICommunicator* comm){

  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);
  MRISequence* MyRECSequence = new MRISequence(false/*Cyclic Sequence*/);

  // ADD FILE TO SEQUENCE
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
  MyMRIScan->ReadPltFile(opts->inputFileName, true);
  MyMRISequence->AddScan(MyMRIScan);

  // APPLY FULL FILTER
  MyMRISequence->ApplySMPFilter(opts,false,comm);

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
    MRIScan* myScan = new MRIScan(*MyMRISequence->GetScan(0));
    myScan->RebuildFromExpansion(currExp,false);

    // EVAL REYNOLDS STRESSES
    WriteSchMessage("Evaluating Reynolds Stresses...");
    myScan->EvalReynoldsStressComponent(opts->thresholdCriteria);
    WriteSchMessage("Done.\n");

    // ADD SCAN TO RECONSTRUCTION
    MyRECSequence->AddScan(myScan);
  }

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // WRITE OUTPUT FILES TO VTK
  MyRECSequence->ExportToVTK("ReynoldsStressSequence",thresholdCriteria);
}

// ===========================================
// EVAL REYNOLDS STRESSES FROM EXPANSION FILES
// ===========================================
void EvalPressureFromExpansion(MRIOptions* opts){

  // SET PARAMETERS
  bool applyThreshold = true;
  double thresholdRatio = 0.0;
  bool doPressureSmoothing = true;

  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // ADD FILE TO SEQUENCE
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
  //MyMRIScan->ReadFromExpansionFile(inFileName,applyThreshold,thresholdRatio);
  MyMRIScan->ReadPltFile(opts->inputFileName,true);
  MyMRISequence->AddScan(MyMRIScan);
  // APPLY MEDIAN FILTER TO VELOCITIES
  //MyMRISequence->GetScan(0)->ApplyMedianFilter(kQtyVelocityX,1);
  //MyMRISequence->GetScan(0)->ApplyMedianFilter(kQtyVelocityY,1);
  //MyMRISequence->GetScan(0)->ApplyMedianFilter(kQtyVelocityZ,1);

  //MyMRISequence->GetScan(0)->ThresholdQuantity(kQtyVelocityZ,1.0e10);

  // EVAL REYNOLDS STRESSES AND PRESSURE GRADIENTS
  MyMRISequence->ComputePressureGradients(opts->thresholdCriteria);

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

  if(opts->outputFormatType == otFILETECPLOT){
    // WRITE OUTPUT FILES TO TECPLOT
    MyMRISequence->ExportToTECPLOT(opts->outputFileName);
  }else{
    // Init a Threshold with No quantity
    MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

    // WRITE OUTPUT FILES TO VTK
    MyMRISequence->ExportToVTK(opts->outputFileName,thresholdCriteria);
  }
}

// =========================================
// EVAL PRESSURES FROM EXPANSION COFFICIENTS
// =========================================
void EvalConcentrationGradient(MRIOptions* opts){

  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  bool doPressureSmoothing = false;

  // ADD FILE TO SEQUENCE
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
  MyMRIScan->ReadPltFile(opts->inputFileName,true);
  MyMRIScan->ScalePositions(0.0058);
  MyMRISequence->AddScan(MyMRIScan);

  // EVAL REYNOLDS STRESSES AND PRESSURE GRADIENTS
  MyMRISequence->GetScan(0)->ComputeQuantityGradient(kQtyConcentration);

  // EVAL RELATIVE PRESSURE
  MyMRISequence->ComputeRelativePressure(doPressureSmoothing);

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // WRITE OUTPUT FILES TO VTK
  MyMRISequence->ExportToVTK(opts->outputFileName,thresholdCriteria);

}

// ===================
// PERFORM RANDOM TEST
// ===================
void PerformRandomTest(MRIOptions* opts,MRICommunicator* comm){
  // Set parameters
  int numberMC = 500;

  // Generating Random Numbers
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator(boost::mt19937(time(0)),boost::normal_distribution<>());

  // Declare Scan
  MRIStructuredScan* MyMRIScan;

  // Monte Carlo Loops
  for(int loopA=0;loopA<numberMC;loopA++){

    // Generate Random Velocity Field - Standard Gaussian Distribution
    MyMRIScan = new MRIStructuredScan(0.0);

    vector<double> params;
    params.push_back(5.0);
    params.push_back(5.0);
    params.push_back(5.0);
    params.push_back(1.0);
    params.push_back(1.0);
    params.push_back(1.0);
    params.push_back(0.0);
    params.push_back(1.0);
    MyMRIScan->CreateSampleCase(kZeroVelocity,params);

    // Assign Random Component
    MyMRIScan->AssignRandomComponent(kdirX,generator);
    MyMRIScan->AssignRandomComponent(kdirY,generator);
    MyMRIScan->AssignRandomComponent(kdirZ,generator);

    // Write Generated Velocities To File
    MyMRIScan->ExportVelocitiesToFile("InputVelocities.dat",true);

    // Filter Velocities - NO GLOBAL ATOMS
    MyMRIScan->applySMPFilter(opts,false,comm);
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
void CropAndComputeVol(MRIOptions* opts, MRICommunicator* comm){
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
  MyMRISequence->ReadFromVolSequence(opts->inputFileName);

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

  // ------------------------
  // PRELIMINARY THRESHOLDING
  // ------------------------
  MyMRISequence->ApplyThresholding(opts->thresholdCriteria);

  // -----------------
  // APPLY FULL FILTER
  // -----------------
  MyMRISequence->ApplySMPFilter(opts,false,comm);

  // --------------------------------------------
  // Apply Boundary Condition Filter to all scans
  // --------------------------------------------
  if(opts->applyBCFilter){
    // APPLY FILTER
    MyMRISequence->ApplySMPFilter(opts,true,comm);
  }

  // ---------------------
  // Apply Final Threshold
  // ---------------------
  MyMRISequence->ApplyThresholding(opts->thresholdCriteria);

  // -------------------------
  // Compute Pressure Gradient
  // -------------------------
  MyMRISequence->ComputePressureGradients(opts->thresholdCriteria);

  // -------------------------
  // Compute Relative Pressure
  // -------------------------
  MyMRISequence->ComputeRelativePressure(false);

  // ---------------------------
  // Export Sequence To VOL File
  // ---------------------------
  MyMRISequence->ExportToTECPLOT(opts->outputFileName);
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
  MRIStructuredScan* myScan = new MRIStructuredScan(0.0);

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
  MRIStructuredScan* myScan = new MRIStructuredScan(0.0);

  // Read From File
  myScan->ReadPltFile(inFileName,true);

  // Set Velocities constant in Z
  for(int loopA=0;loopA<myScan->totalCellPoints;loopA++){
    myScan->cellPoints[loopA].velocity[0] = 0.01;
    myScan->cellPoints[loopA].velocity[1] = 0.01;
    myScan->cellPoints[loopA].velocity[2] = 1.0;
  }

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // Export to VTK
  myScan->ExportToVTK(outfileName+"_Model.vtk",thresholdCriteria);

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
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
  MyMRIScan->ReadFromExpansionFile(coeffFileName,performThreshold,thresholdType,threshold);
  MyMRISequence->AddScan(MyMRIScan);

  // APPLY THRESHOLDING
  //MRIThresholdCriteria criteria(kCriterionLessThen,kQtyConcentration,1000.0);
  //MyMRISequence->ApplyThresholding(criteria);

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // EXPORT TO PLT FILE
  //MyMRISequence->ExportToTECPLOT(plotOut);
  MyMRISequence->ExportToVTK(plotOut, thresholdCriteria);
}

// ============================
// EVAL VARIOUS VORTEX CRITERIA
// ============================
void EvalVortexCriteria(MRIOptions* opts){
  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // ADD FILE TO SEQUENCE
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
  MyMRIScan->ReadPltFile(opts->inputFileName,true);
  MyMRISequence->AddScan(MyMRIScan);

  // EVAL VORTEX CRITERIA
  // MyMRISequence->GetScan(0)->EvalVortexCriteria();
  // MyMRISequence->GetScan(0)->EvalVorticity();
  MyMRISequence->GetScan(0)->EvalEnstrophy(opts->thresholdCriteria);

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // WRITE OUTPUT FILES TO VTK
  MyMRISequence->ExportToVTK(opts->outputFileName,thresholdCriteria);
}

// ===========================================
// WRITE SPATIAL DISTRIBUTIONS OF COEFFICIENTS
// ===========================================
void WriteSpatialExpansion(MRIOptions* opts){
  // CREATE NEW SEQUENCES
  MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

  // CREATE NEW SCAN
  MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);

  // READ FROM EXPANSION COEFFICIENTS
  MyMRIScan->ReadFromExpansionFile(opts->inputFileName,false,kSoftThreshold,0.0);

  // ADD TO SEQUENCE
  MyMRISequence->AddScan(MyMRIScan);

  // SPATIALLY EVALUATE VORTEX COEFFICIENTS
  MyMRISequence->GetScan(0)->EvalSMPVortexCriteria(MyMRISequence->GetScan(0)->expansion);

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // EXPORT TO VTK
  MyMRISequence->ExportToVTK(opts->outputFileName,thresholdCriteria);
}

// ===============================
// RUN APPLICATION IN NORMAL MODE
// ===============================
void runApplication(MRIOptions* opts, MRICommunicator* comm){

  MRISequence* MyMRISequence;

  // MASTER PROCESSOR DOES THE READING
  if(comm->currProc == 0){

      // INIT SEQUENCE
      MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

      // LOOP ON THE NUMBER OF SCANS
      for(size_t loopA=0;loopA<opts->sequenceFileList.size();loopA++){

        // CREATE NEW SCAN
        MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);

        // CHOOSE INPUT FORMAT
        if(opts->inputFormatType == itTEMPLATE){
          // CREATE TEMPLATE
          MyMRIScan->CreateSampleCase(opts->templateType,opts->templateParams);
        }else if(opts->inputFormatType == itEXPANSION){
          // READ FROM EXPANSION COEFFICIENTS
          bool applyThreshold = true;
          int thresholdType = kHardThresold;
          double thresholdRatio = 0.5;
          MyMRIScan->ReadFromExpansionFile(opts->sequenceFileList[loopA],applyThreshold,thresholdType,thresholdRatio);
        }else if (opts->inputFormatType == itFILEVTK){
          // READ FROM FILE
          MyMRIScan->ReadVTKStructuredPoints(opts->sequenceFileList[loopA], true);
        }else if (opts->inputFormatType == itFILETECPLOT){
          // READ FROM FILE
          MyMRIScan->ReadPltFile(opts->sequenceFileList[loopA], true);
        }
        // ADD TO SEQUENCE
        MyMRISequence->AddScan(MyMRIScan);
      }
  }else{
    // LOOP ON THE NUMBER OF SCANS
      // INIT SEQUENCE
      MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);
//    for(size_t loopA=0;loopA<opts->sequenceFileList.size();loopA++){
      // CREATE EMPTY SCAN
      MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
      // ADD TO SEQUENCE
      MyMRISequence->AddScan(MyMRIScan);
//    }
  }

  // All processes are waiting for the root to read the files
  int mpiError = MPI_Barrier(comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);

  // Distribute Sequence Data using MPI
  if(comm->totProc > 1){
    MyMRISequence->DistributeSequenceData(comm);
  }
  if(comm->currProc == 0){
    printf("Data Structure Communication OK.\n");
  }

  // SAVE INITIAL VELOCITIES ON ROOT PROCESSOR
  if(comm->currProc == 0){
    if (opts->saveInitialVel){
      MyMRISequence->saveVelocity();
    }
  }

  // APPLY NOISE
  if (opts->applyNoise){
    MyMRISequence->applyNoise(opts->noiseIntensity);
  }

  // FILTER DATA IF REQUIRED
  if(comm->currProc == 0){
    if (opts->applyMedianFilter){
      MyMRISequence->ApplyMedianFilter(kQtyVelocityX,opts->filterNumIterations,opts->filterUseMedian);
      MyMRISequence->ApplyMedianFilter(kQtyVelocityY,opts->filterNumIterations,opts->filterUseMedian);
      MyMRISequence->ApplyMedianFilter(kQtyVelocityZ,opts->filterNumIterations,opts->filterUseMedian);
    }
  }

  // CLEAN VELOCITY COMPONENTS ON BOUNDARY
  if(comm->currProc == 0){
    if (opts->cleanBoundaryVelocities){
      MyMRISequence->cleanNormalComponentOnBoundary();
    }
  }

  // INTERPOLATE BOUNDARY VELOCITIES WITH POLYNOMIALS
  if(comm->currProc == 0){
    if (opts->interpolateBoundaryVelocities){
      MyMRISequence->InterpolateBoundaryVelocities();
    }
  }

  //if(comm->currProc == 0){
  //  // Open Output File
  //  FILE* outFile;
  //  outFile = fopen("testGauss.log","w");
  //  // Write Header

  //  for(int loopA=0;loopA<MyMRISequence->GetScan(0)->totalCellPoints;loopA++){
  //    fprintf(outFile,"%e %e %e\n",MyMRISequence->GetScan(0)->cellPoints[loopA].velocity[0],MyMRISequence->GetScan(0)->cellPoints[loopA].velocity[1],MyMRISequence->GetScan(0)->cellPoints[loopA].velocity[2]);
  //  }
  //  // Close Output file
  //  fclose(outFile);
  //}


  // APPLY FULL FILTER
  if (opts->applySMPFilter){
    MyMRISequence->ApplySMPFilter(opts,false,comm);
  }

  // APPLY BOUNDARY CONDITION FILTER
  if (opts->applyBCFilter){
    MyMRISequence->ApplySMPFilter(opts,true,comm);
  }

  // APPLY THRESHOLD
  if(comm->currProc == 0){
    MyMRISequence->ApplyThresholding(opts->thresholdCriteria);
  }

  // EVAL VORTEX CRITERIA
  if(comm->currProc == 0){
    if(opts->evalPopVortexCriteria){
      MyMRISequence->EvalVortexCriteria(opts->thresholdCriteria);
      MyMRISequence->EvalVorticity(opts->thresholdCriteria);
      MyMRISequence->EvalEnstrophy(opts->thresholdCriteria);
    }
  }

  if(comm->currProc == 0){
    if(opts->applySMPFilter && opts->evalSMPVortexCriterion){
      // SPATIALLY EVALUATE VORTEX COEFFICIENTS
      MyMRISequence->EvalSMPVortexCriteria();
      // Threshold Vortex Expansion
      //void ApplyVortexThreshold(int thresholdType, double ratio);
      // Eval 2-Norm of Coefficient Vector
      //double Get2Norm(bool onlyVortex);
    }
  }

  if(opts->evalPressure){
    // Compute Pressure Gradient
    MyMRISequence->ComputePressureGradients(opts->thresholdCriteria);

    // Compute Relative Pressure
    MyMRISequence->ComputeRelativePressure(false);
  }

  // SAVE EXPANSION COEFFICIENTS IF REQUESTED
  if(comm->currProc == 0){
    if(opts->applySMPFilter && opts->saveExpansionCoeffs){
      MyMRISequence->WriteExpansionFile(std::string(opts->outputFileName + "_expCoeff"));
    }
  }

  // SAVE FILE FOR POISSON COMPUTATION
  if(comm->currProc == 0){
    if (opts->exportToPoisson){
      MyMRISequence->ExportForPOISSON(opts->poissonFileName,opts->density,opts->viscosity,opts->thresholdCriteria);
    }
  }

  // EXPORT FILE
  if(comm->currProc == 0){
    if(opts->outputFormatType == itFILEVTK){
      // READ FROM FILE
      MyMRISequence->ExportToVTK(opts->outputFileName,opts->thresholdCriteria);
    }else if (opts->outputFormatType == itFILETECPLOT){
      // READ FROM FILE
      MyMRISequence->ExportToTECPLOT(opts->outputFileName);
    }else{
      throw MRIException("ERROR: Invalid output file format.\n");
    }
  }

  // Free Memory
  delete MyMRISequence;
}

// ============
// MAIN PROGRAM
// ============
int main(int argc, char **argv){

  // Init MRI Communicator
  MRICommunicator* comm = new MRICommunicator();

  // Initialize MPI
  MPI::Init();
  int rank; int nproc;
  comm->mpiComm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm->mpiComm, &comm->currProc);
  MPI_Comm_size(comm->mpiComm, &comm->totProc);

  //  Declare
  int val = 0;
  MRIOptions* options;

  // WRITE PROGRAM HEADER - ONLY MASTER NODE
  if(comm->currProc == 0){
    WriteHeader();

    // Create Options
    options = new MRIOptions();

    // Read Options from Command Line
    int res = options->getCommadLineOptions(argc,argv);
    if(res != 0){
      return -1;
    }

    // Read options from command file if required
    if(options->useCommandFile){
      int res = options->getOptionsFromCommandFile(options->commandFileName);
      if(res != 0){
        return -1;
      }
    }

    // Generate Command File if needed
    if(options->generateCommandFile){
      int res = options->writeCommandFilePrototype(options->commandFileName);
      if(res != 0){
        return -1;
      }
      return 0;
    }
  }else{
    options = new MRIOptions();
  }

  // Wait for all processes
  int mpiError = MPI_Barrier(comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);

  // Distribute Options using MPI
  if(comm->totProc > 1){
    options->DistributeProgramOptions(comm);
  }
  if(comm->currProc == 0){
    printf("Program Options Communication OK.\n");
  }
  //string optOut("optionsOut_" + to_string(comm->currProc) + ".out");
  //options->writeOptionsToFile(optOut);

  // Finalize options
  options->finalize();

  // ============
  // MAIN PROGRAM
  // ============
  try{
    // Write Program Help
    switch(options->runMode){
      case rmHELP:
        if(comm->currProc == 0){
          // Write Program Help
          MRIUtils::WriteProgramHelp();
        }
        break;
      // PREFERRED RUNNING MODE
      case rmNORMAL:
      {
        runApplication(options,comm);
        break;
      }
      case rmEVALSEQUENCEPRESSURE:
      {
        // Eval Pressure From Velocity Sequence
        EvalSequencePressure(options,comm);
        break;
      }
      case rmEVALPRESSUREFROMSIGNATUREFLOW:
      {
        // Eval Pressure for signature Flow Cases
        EvalPressureFromSignatureFlow(options);
        break;
      }
      case rmPLTTOVTK:
      {
        // Convert to VTK
        ConvertTECPLOToVTK(options);
        break;
      }
      case rmEVALSCANSTATISTICS:
      {
        // COMPUTE SCAN STATISTICS
        string firstFileName("firstFile.txt");
        string secondFileName("secondFile.txt");
        string statFileName("stats");
        ComputeScanStatistics(firstFileName,secondFileName,statFileName);
        break;
      }
      case rmCOMUTESCANMATRICES:
        // COMPUTE SCAN MATRICES
        ComputeScanMatrices();
        break;
      case rmPERFORMRANDOMTEST:
        // PERFORM RANDOM TEST
        PerformRandomTest(options,comm);
        break;
      case rmCROPANDCOMPUTEVOLUME:
        // CROP AND REDUCE VOL FILE
        CropAndComputeVol(options,comm);
        break;
      case rmSTREAMLINETEST1:
      {
        // PERFORM STREAMLINE TEST 1
        int intValue = 0;
        string inFileName("inputFile.txt");
        string outfileName("outputFile.txt");
        PerformStreamlineTest1(intValue,inFileName,outfileName);
        break;
      }
      case rmSTREAMLINETEST2:
      {
        // PERFORM STREAMLINE TEST 2
        string inFileName("inFile.txt");
        string outFileName("outFile.txt");
        PerformStreamlineTest2(inFileName,outFileName);
        break;
      }
      case rmPRINTTHRESHOLDINGTOVTK:
        // Test Expansion Coefficients
        // TEST_ExpansionCoefficients(inFileName);
        TEST02_PrintThresholdingToVTK(options,comm);
        break;
      case rmEVALREYNOLDSSTRESSES:
        // Test Expansion Coefficients
        TEST03_EvalReynoldsStresses(options,comm);
        break;
      case rmSHOWFACEFLUXPATTERS:
      {
        // READ FACE FLUXES FROM FILE AND EXPORT TO VTK
        string faceFluxFileName("faceFluxFile.txt");
        string outFileName("outFileName.txt");
        ShowFaceFluxPatterns(faceFluxFileName,outFileName);
        break;
      }
      case rmBUILDFROMCOEFFICIENTS:
      {
        // Get File Names
        string coeffFileName("coeffFileName.txt");
        string plotOut("outputFile.txt");
        bool performThreshold = false;
        int thresholdType = 0;
        double thresholdValue = 0.0;
        // READ FROM COEFFICIENT FILE AND EXPORT TO PLT
        BuildFromCoeffs(coeffFileName,plotOut,performThreshold,thresholdType,thresholdValue);
        break;
      }
      case rmEVALPRESSURE:
        // EVAL PRESSURES FROM EXPANSION COFFICIENTS
        EvalPressureFromExpansion(options);
        break;
      case rmEVALCONCGRADIENT:
        // Eval Concetrantion Gradient
        EvalConcentrationGradient(options);
        break;
      case rmEVALVORTEXCRITERIA:
        // Eval Vortex Criterion
        EvalVortexCriteria(options);
        break;
      case rmWRITESPATIALEXPANSION:
        // EVAL Spatial Expansion Coefficients Distribution
        WriteSpatialExpansion(options);
        break;
      case rmSOLVEPOISSON:
        // Solve Poisson Equation to Compute the pressures
        SolvePoissonEquation(comm);
        break;
      default:
        // Invalid Switch
        std::string currMsgs("Error: Invalid Switch. Terminate.\n");
        WriteSchMessage(currMsgs);
        return (1);
        break;
    }
  }catch (std::exception& ex){
    if(comm->currProc == 0){
      WriteSchMessage(std::string(ex.what()));
      WriteSchMessage(std::string("\n"));
      WriteSchMessage(std::string("Program Terminated.\n"));
    }
    // Finalize MPI
    delete options;
    MPI::Finalize();
    return -1;
  }
  if(comm->currProc == 0){
    WriteSchMessage(string("\n"));
    WriteSchMessage(string("Program Completed.\n"));
  }
  // Finalize MPI
  //delete comm;
  delete options;
  MPI::Finalize();
  return 0;
}

