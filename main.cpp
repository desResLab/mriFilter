#include <iostream>

#include "mriSequence.h"
#include "mriScan.h"
#include "mriUtils.h"
#include "mriStatistics.h"
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

  
// ============
// MAIN PROGRAM
// ============
int main(int argc, char **argv)
{
  // Write Program Header
  WriteHeader();
  // Testing Functionalities for the Program
  // TestRandomNumberFacilities();
  // Store Local Option
  std::string firstOption(argv[1]);
  if ((argc == 1)||(firstOption == "-?")||(firstOption == "-help")||(firstOption == "-h")||(firstOption == "")){
    // Write Program Help
    MRIUtils::WriteProgramHelp();
  }else if (firstOption == "-s"){
    // ================
    // ANALYZE SEQUENCE
    // ================
    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);    
    
    // -------------------
    // Create New Sequence
    // -------------------
    MRISequence* MyMRISequence = new MRISequence(true/*Cyclic Sequence*/);
    
    // -------------------
    // Read From VOL Files
    // -------------------
    MyMRISequence->ReadFromVolSequence(inFileName);
    
    // --------------------------------------
    // Apply MP Filter to all Scan Separately
    // --------------------------------------
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
    MyMRISequence->ExportToVOL(outfileName);
    
  }else if (firstOption == "-o"){
    // ----------------
    // READ SINGLE SCAN
    // ----------------
    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);
    
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
    
    // -------------------------
    // Compute Pressure Gradient
    // -------------------------
    MyMRISequence->ComputePressureGradients();
    
    // -------------------------
    // Compute Relative Pressure
    // -------------------------
    MyMRISequence->ComputeRelativePressure(false);
    
    // Select the resulting File
    MyMRISequence->ExportToTECPLOT(outfileName);
    
  }else if (firstOption == "-art"){
    // -------------------------------------------
    // Launch Calculation for article with Filippo
    // -------------------------------------------
    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);
    // Get Parameters
    double itTol = atof(argv[4]);
    std::string thresholdTypeString(argv[5]);
    double thresholdValue = atof(argv[6]);
        
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
    MRIOptions Options(itTol,2000);
    bool useBCFilter = false;
    bool useConstantPatterns = true;
    //MRIThresholdCriteria criteria(kCriterionLessThen,kQtyVelModule,1.0e-4);
    int thresholdType = 0;
    if (thresholdTypeString == "Conc"){
      thresholdType = kQtyConcentration;
    }else{
      thresholdType = kQtyVelModule;
    }
    MRIThresholdCriteria criteria(kCriterionLessThen,thresholdType,thresholdValue);
    
    // APPLY FULL FILTER
    MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);
    
    // Apply Boundary Condition Filter to all scans
    useBCFilter = true;
    // APPLY FILTER
    MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

    // Apply Final Threshold
    MyMRISequence->ApplyThresholding(criteria);
    
    // Evaluate Statistics
    //MyMRISequence->EvalStatistics();
    
    // -------------------------
    // Compute Pressure Gradient
    // -------------------------
    //MyMRISequence->ComputePressureGradients();
    
    // -------------------------
    // Compute Relative Pressure
    // -------------------------
    //MyMRISequence->ComputeRelativePressure(false);
    
    // SELECT THE RESULTING FILE
    MyMRISequence->ExportToTECPLOT(outfileName);
    
  }else if (firstOption == "-filt"){
    // -------------------------------------------
    // Launch Calculation for article with Filippo
    // -------------------------------------------
    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);
    // Get Parameters
    double itTol = atof(argv[4]);

    // Create New Sequence
    MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

    // Add File to Sequence
    MRIScan* MyMRIScan = new MRIScan(0.0);
    MyMRIScan->ReadPltFile(inFileName, true);
    MyMRISequence->AddScan(MyMRIScan);

    // Echo Inserted Parameters
    WriteSchMessage(std::string("--------------------------------------------\n"));
    WriteSchMessage(std::string("MP Iteration tolerance Value: ")+MRIUtils::FloatToStr(itTol)+"\n");
    WriteSchMessage(std::string("--------------------------------------------\n"));

    // SET OPTIONS AND THRESHOLD
    MRIOptions Options(itTol,2000);
    bool useBCFilter = false;
    bool useConstantPatterns = true;
    int thresholdType = kQtyConcentration;
    double thresholdValue = 0.0;
    MRIThresholdCriteria criteria(kCriterionLessThen,thresholdType,thresholdValue);

    // APPLY FULL FILTER
    MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

    // SAVE FILE
    MyMRISequence->ExportToTECPLOT(outfileName);

  }else if (firstOption == "-tvol"){
    
    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);    

    // Create New Sequence
    MRISequence* MyMRISequence = new MRISequence(true/*Cyclic Sequence*/);
    
    // Read VOL Sequence
    MyMRISequence->ReadFromVolSequence(inFileName);
    
    // Select the resulting File
    MyMRISequence->ExportToTECPLOT(outfileName);
    
  }else if (firstOption == "-stats"){
    // =============================
    // EVAL STATISTICS BETWEEN SCANS
    // =============================
    // Get File Name
    std::string firstFileName(argv[2]);
    std::string secondFileName(argv[3]); 
    std::string statFileName(argv[4]); 
    
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
    
  }else if (firstOption == "-Mat"){
    // VAR
    int totalERows = 0;
    int totalECols = 0;
    double** EMat = nullptr;
    int totalDRows = 0;
    int totalDCols = 0;
    double** DMat = nullptr;
    int totalStarRows = 0;
    int totalStarCols = 0;
    double** StarMatrix = nullptr;
    
    // Print Matrices for Matlab Analysis Of Variance
    MRIScan* MyMRIScan = new MRIScan(0.0);
    // ISOTROPIC
    //MyMRIScan->CreateSampleCase(kConstantFlow,5,5,5,1.0,1.0,1.0,0.0,kdirX);  
    // ANISOTROPIC
    MyMRIScan->CreateSampleCase(kConstantFlow,5,5,5,1.0,2.0,3.0,0.0,kdirX);  
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
  }else if (firstOption == "-bctest"){
    // Create New Sequence
    MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);
    
    // Add File to Sequence
    MRIScan* MyMRIScan = new MRIScan(0.0);
    MyMRIScan->CreateSampleCase(kConstantFlowWithStep,19,19,19,1.0,1.0,1.0,0.0,kdirX);
    MyMRISequence->AddScan(MyMRIScan);
    
    // Export As Original
    MyMRISequence->ExportToTECPLOT("BCCube_Original.plt");
    
    // Filter 
    // SET OPTIONS AND THRESHOLD
    MRIOptions Options(1.0e-7,2000);
    bool useBCFilter = false;
    bool useConstantPatterns = true;
    MRIThresholdCriteria criteria(kCriterionLessThen,kQtyConcentration,1.0);

    // APPLY FULL FILTER
    MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);
    
    // Export Filtered
    MyMRISequence->ExportToTECPLOT("BCCube_Filtered.plt");
    
    // APPLY FILTER
    useBCFilter = true;
    MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);
    MyMRISequence->ApplyThresholding(criteria);
    
    // Export Filtered BC
    MyMRISequence->ExportToTECPLOT("BCCube_FilteredBC.plt");
    
  }else if (firstOption == "-randtest"){
    // Set parameters
    int numberMC = 500;
    // SET OPTIONS AND THRESHOLD
    MRIOptions Options(1.0e-10,2000);
    bool useBCFilter = false;
    bool useConstantPatterns = true;
    MRIThresholdCriteria criteria(kCriterionLessThen,kQtyConcentration,1.0);
    
    // Generating Random Numbers
    boost::variate_generator<boost::mt19937, boost::normal_distribution<>> generator(boost::mt19937(time(0)),boost::normal_distribution<>());

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
      MyMRIScan = nullptr;
    }
  }else if (firstOption == "-reduceVol"){
    // ================
    // ANALYZE SEQUENCE
    // ================
    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);    
    
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
    
  }else if (firstOption == "-streamLines"){
    // Get File Name
    int intValue = atoi(argv[2]);
    std::string inFileName(argv[3]);
    std::string outfileName(argv[4]);
    
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
      
  }else if (firstOption == "-streamLinesTest"){
    // Get File Name
    std::string inFileName(argv[3]);
    std::string outfileName(argv[4]);

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

  }else if (firstOption == "-vtk"){
    // ------------------------
    // Read and Write File Only
    // ------------------------
    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);

    // Add File to Sequence
    MRIScan* MyMRIScan = new MRIScan(0.0);
    MyMRIScan->ReadPltFile(inFileName, true);

    // SELECT THE RESULTING FILE
    //MyMRIScan->ExportToTECPLOT(outfileName,true);
    MyMRIScan->ExportToVTK(outfileName);
  }else if (firstOption == "-testEnergy"){

    // Get File Name
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);

    // Create New Sequence
    MRISequence* MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);

    // Add File to Sequence
    MRIScan* MyMRIScan = new MRIScan(0.0);
    MyMRIScan->ReadPltFile(inFileName, true);
    MyMRISequence->AddScan(MyMRIScan);

    // SET OPTIONS AND THRESHOLD
    MRIOptions Options(1.0e-4,2000);
    bool useBCFilter = false;
    bool useConstantPatterns = true;
    //MRIThresholdCriteria criteria(kCriterionLessThen,kQtyVelModule,1.0e-4);
    int thresholdType = 0;
    thresholdType = kQtyConcentration;
    MRIThresholdCriteria criteria(kCriterionLessThen,thresholdType,0.0);

    // APPLY FULL FILTER
    MyMRISequence->ApplyMPFilter(Options,useBCFilter,useConstantPatterns,criteria);

    // SELECT THE RESULTING FILE
    MyMRISequence->ExportToTECPLOT(outfileName);

  }else if (firstOption == "-dico"){
    // Get File Names
    std::string inFileName(argv[2]);
    std::string outfileName(argv[3]);

    // Perform Conversion
    ConvertDICOMToVTK(inFileName,outfileName);

  }else{
    // Invalid switch
    std::string currMsgs("Error: Invalid Switch. Terminate.\n");
    WriteSchMessage(currMsgs);
    return (1);
  }
  // Message: Completed!
  std::string currMsgs("\nOperation Completed!\n");
  WriteSchMessage(currMsgs);
  return (0);
};

