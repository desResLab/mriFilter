# include <iostream>
# include <boost/random.hpp>

# include "mriScan.h"
# include "mriSequence.h"
# include "mriUtils.h"
# include "mriConstants.h"
# include "mriOptions.h"
# include "mriCommunicator.h"

# include "mpi.h"

using namespace std;

// =============================
// EVAL STATISTICS BETWEEN SCANS
// =============================
void ComputeScanStatistics(std::string firstFileName,std::string secondFileName,std::string statFileName){

    // Init File Names
    string statFileNameFirst = statFileName+"_First.dat";
    string statFileNameSecond = statFileName+"_Second.dat";
    string statFileNameDiff = statFileName+"_Diff.dat";

    // Cyclic Sequence
    bool isCyclicSequence = false;

    // Read the File in a Sequence
    // Create New Sequence
    MRISequence* seq = new MRISequence(isCyclicSequence);

    // Add First File to Sequence
    seq->readPLTFile(firstFileName, true);

    // Create New Sequence
    seq->readPLTFile(secondFileName, true);

    // Create Statistics
    bool useBox = false;
    int numberOfBins = 301;
    MRIDoubleVec limitBox(6);
    
    // Set Limits
    limitBox[0] = seq->topology->domainSizeMin[0];
    limitBox[1] = seq->topology->domainSizeMax[0];
    limitBox[2] = seq->topology->domainSizeMin[1];
    limitBox[3] = seq->topology->domainSizeMax[1];
    limitBox[4] = seq->topology->domainSizeMin[2];
    limitBox[5] = seq->topology->domainSizeMax[2];
    
    // Apply Factors to limitBox
    double xFactor = 1.0;
    double yFactor = 1.0;
    double zFactor = 1.0;
    MRIUtils::applyLimitBoxFactors(xFactor,yFactor,zFactor,limitBox);
    
    // Allocate Bin Arrays
    MRIDoubleVec binCenters(numberOfBins);
    MRIDoubleVec binValues(numberOfBins);
    // Eval Single PDFs
    // FIRST
    seq->getScan(0)->evalScanPDF(kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
    MRIUtils::printBinArrayToFile(statFileNameFirst,numberOfBins,binCenters,binValues);

    // SECOND
    seq->getScan(0)->evalScanPDF(kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
    MRIUtils::printBinArrayToFile(statFileNameSecond,numberOfBins,binCenters,binValues);

    // DIFFERENCE
    seq->evalScanDifferencePDF(1,0,kQtyVelModule,numberOfBins,useBox,limitBox,binCenters,binValues);
    MRIUtils::printBinArrayToFile(statFileNameDiff,numberOfBins,binCenters,binValues);

    // Free Sequence
    delete seq;
}

// =============================
// EXPLICITLY EVAL SCAN MATRICES
// =============================
void ComputeScanMatrices(){
  // VAR
  int totalERows = 0;
  int totalECols = 0;
  MRIDoubleMat EMat;
  int totalDRows = 0;
  int totalDCols = 0;
  MRIDoubleMat DMat;
  int totalStarRows = 0;
  int totalStarCols = 0;
  MRIDoubleMat StarMatrix;

  // SET PARAMETERS
  bool isIsotropic = true;

  // Cyclic Sequence
  bool isCyclicSequence = false;
  // Create Scan
  MRISequence* seq = new MRISequence(isCyclicSequence);
  
  // Set Template Parameters
  MRIDoubleVec params(8);
  params[0] = 5;
  params[1] = 5;
  params[2] = 5;
  params[3] = 1.0;
  params[4] = 1.0;
  params[5] = 1.0;
  params[6] = 0.0;
  params[7] = 1.0;
  
  // SET ISOTROPIC OR ANISOTROPIC CASE
  if(isIsotropic){
    // ISOTROPIC
    seq->createSampleCase(kConstantFlow,params);
  }else{
    // ANISOTROPIC
    params[3] = 1.0;
    params[4] = 2.0;
    params[5] = 3.0;
    seq->createSampleCase(kConstantFlow,params);
  }

  // RETRIEVE OPERATORS IN MATRIX FORM
  // ENCODING
  seq->getScan(0)->assembleEncodingMatrix(totalERows,totalECols,EMat);
  MRIUtils::printMatrixToFile("EncodingMat.dat",EMat);
  // DECODING
  seq->getScan(0)->assembleDecodingMatrix(totalDRows,totalDCols,DMat);
  MRIUtils::printMatrixToFile("DecodingMat.dat",DMat);
  // VORTEX FRAME MATRIX
  seq->getScan(0)->assembleStarMatrix(totalStarRows,totalStarCols,StarMatrix);
  MRIUtils::printMatrixToFile("StarMat.dat",StarMatrix);
}

// ================
// READ FACE FLUXES
// ================
void ShowFaceFluxPatterns(std::string faceFluxFileName, std::string outFileName){
  
  // Cyclic Sequence
  bool isCyclicSequence = false;
  MRISequence* seq = new MRISequence(isCyclicSequence);

  bool isIsotropic = true;

  MRIDoubleVec params(7);
  params[0] = 5.0;
  params[1] = 5.0;
  params[2] = 5.0;
  params[3] = 1.0;
  params[4] = 1.0;
  params[5] = 1.0;
  params[6] = 1.0;
  if(isIsotropic){
    // ISOTROPIC
    seq->createSampleCase(kConstantFlow,params);
  }else{
   // ANISOTROPIC
   params[3] = 1.0;
   params[4] = 2.0;
   params[5] = 3.0;
   seq->createSampleCase(kConstantFlow,params);
  }

  // READ FACE FLUXES FROM FILE
  int totalRows = 0;
  int totalCols = 0;
  MRIDoubleMat faceFluxMat;
  MRIUtils::readMatrixFromFile(faceFluxFileName,totalRows,totalCols,faceFluxMat);

  // COPY THE INTERESTING COLUMN
  double faceFluxVec[totalRows];
  for(int loop0=0;loop0<100;loop0++){
    int selectedCol = loop0;
    for(int loopA=0;loopA<totalRows;loopA++){
      faceFluxVec[loopA] = faceFluxMat[loopA][selectedCol];
    }

    // TRASFORM FACE FLUXES IN VELOCITIES
    seq->getScan(0)->recoverCellVelocitiesRT0(false,faceFluxVec);
    
    // UPDATE VELOCITIES
    seq->getScan(0)->updateVelocities();

    // Init a Threshold with No quantity
    MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

    // EXPORT TO VTK
    seq->exportToVTK(outFileName + "_" + to_string(loop0) + ".vtk",thresholdCriteria);
  }
}

// ==================
// PRINT THRESHOLDING
// ==================
void TEST02_PrintThresholdingToVTK(MRIOptions* opts, MRICommunicator* comm){

  // Create New Sequence
  bool isCyclicSequence = false;
  MRISequence* mriSeq = new MRISequence(isCyclicSequence);
  MRISequence* recSeq = new MRISequence(isCyclicSequence);

  // ADD FILE TO SEQUENCE
  mriSeq->readPLTFile(opts->inputFileName, true);

  // APPLY FULL FILTER
  mriSeq->applySMPFilter(opts,false,comm);

  // APPLY SUCCESSIVE THRESHOLD
  for(int loopA=0;loopA<5;loopA++){
    writeSchMessage("Applying Threshold " + MRIUtils::intToStr(loopA) + "\n");
    
    // CREATE NEW EXPANSION
    MRIExpansion* currExp = new MRIExpansion(mriSeq->getScan(0)->expansion);
    
    // APPLY THRESHOLD
    currExp->applyVortexThreshold(kSoftThreshold,(0.95/4.0)*(loopA));
    
    // WRITE EXPANSION TO FILE
    currExp->writeToFile("Expansion_" + MRIUtils::intToStr(loopA) + "\n");
    
    // GET A SCAN FROM ORIGINAL SEQUENCE
    MRIScan* myScan = new MRIScan(*mriSeq->getScan(0));
    myScan->rebuildFromExpansion(currExp,false);
    
    // ADD SCAN TO RECONSTRUCTION
    recSeq->addScan(myScan);
  }

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // WRITE OUTPUT FILES TO VTK
  mriSeq->exportToVTK("FilteredSeq",thresholdCriteria);
  recSeq->exportToVTK("ReconstructedSeq",thresholdCriteria);
}

// ======================
// EVAL REYNOLDS STRESSES
// ======================
void TEST03_EvalReynoldsStresses(MRIOptions* opts, MRICommunicator* comm){

  // CREATE NEW SEQUENCES
  bool isCyclicSequence = false;
  MRISequence* mriSeq = new MRISequence(isCyclicSequence);
  MRISequence* recSeq = new MRISequence(isCyclicSequence);

  // ADD FILE TO SEQUENCE
  mriSeq->readPLTFile(opts->inputFileName, true);

  // APPLY FULL FILTER
  mriSeq->applySMPFilter(opts,false,comm);

  writeSchMessage("Filter Applied!\n");

  // APPLY SUCCESSIVE THRESHOLD
  for(int loopA=0;loopA<2;loopA++){

    // WRITE MESSAGE
    writeSchMessage("Applying Threshold " + MRIUtils::intToStr(loopA) + "\n");

    // CREATE NEW EXPANSION
    MRIExpansion* currExp = new MRIExpansion(mriSeq->getScan(0)->expansion);

    // APPLY THRESHOLD
    currExp->applyVortexThreshold(kHardThresold,(0.99/1.0)*(loopA));

    // WRITE EXPANSION TO FILE
    currExp->writeToFile("Expansion_" + MRIUtils::intToStr(loopA) + "\n");

    // GET A SCAN FROM ORIGINAL SEQUENCE
    MRIScan* myScan = new MRIScan(*mriSeq->getScan(0));
    myScan->rebuildFromExpansion(currExp,false);

    // EVAL REYNOLDS STRESSES
    writeSchMessage("Evaluating Reynolds Stresses...");
    myScan->evalReynoldsStress(opts->thresholdCriteria);
    writeSchMessage("Done.\n");

    // ADD SCAN TO RECONSTRUCTION
    recSeq->addScan(myScan);
  }

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // WRITE OUTPUT FILES TO VTK
  recSeq->exportToVTK("ReynoldsStressSequence",thresholdCriteria);
}

// =========================================
// EVAL PRESSURES FROM EXPANSION COFFICIENTS
// =========================================
void evalConcentrationGradient(MRIOptions* opts){

  // CREATE NEW SEQUENCES
  bool isCyclicSequence = false;
  MRISequence* mriSeq = new MRISequence(isCyclicSequence);

  bool doPressureSmoothing = false;

  // ADD FILE TO SEQUENCE
  mriSeq->readPLTFile(opts->inputFileName,true);
  mriSeq->scalePositions(0.0058);

  // EVAL REYNOLDS STRESSES AND PRESSURE GRADIENTS
  mriSeq->getScan(0)->computeQuantityGradient(kQtyConcentration);

  // EVAL RELATIVE PRESSURE
  mriSeq->computeRelativePressure(doPressureSmoothing);

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // WRITE OUTPUT FILES TO VTK
  mriSeq->exportToVTK(opts->outputFileName,thresholdCriteria);

}

// ===================
// PERFORM RANDOM TEST
// ===================
void performRandomTest(MRIOptions* opts,MRICommunicator* comm){
  
  // Set parameters
  int numberMC = 500;

  // Generating Random Numbers
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> > generator(boost::mt19937(time(0)),boost::normal_distribution<>());

  // Declare Scan
  MRIScan* scan;

  // Monte Carlo Loops
  for(int loopA=0;loopA<numberMC;loopA++){

    // Generate Random Velocity Field - Standard Gaussian Distribution
    scan = new MRIScan(0.0);

    MRIDoubleVec params(8);
    params[0] = 5.0;
    params[1] = 5.0;
    params[2] = 5.0;
    params[3] = 1.0;
    params[4] = 1.0;
    params[5] = 1.0;
    params[6] = 0.0;
    params[7] = 1.0;
    scan->createSampleCase(kZeroVelocity,params);

    // Assign Random Component
    scan->assignRandomComponent(kdirX,generator);
    scan->assignRandomComponent(kdirY,generator);
    scan->assignRandomComponent(kdirZ,generator);

    // Write Generated Velocities To File
    scan->exportVelocitiesToFile("InputVelocities.dat",true);

    // Filter Velocities - NO GLOBAL ATOMS
    scan->applySMPFilter(opts,false,comm);
    scan->updateVelocities();

    // Write Resultant Velocities
    scan->exportVelocitiesToFile("OutputVelocities.dat",true);

    // Deallocate
    delete scan;
    scan = NULL;
  }
}

// ========================================
// BUILD FROM COEFFICIENTS AND WRITE TO PLT
// ========================================
void buildFromCoeffs(std::string coeffFileName,std::string plotOut,bool performThreshold,int thresholdType,double threshold){

  // CREATE NEW SEQUENCES
  MRISequence* seq = new MRISequence(false/*Cyclic Sequence*/);

  // ADD FILE TO SEQUENCE
  seq->readFromExpansionFile(coeffFileName,performThreshold,thresholdType,threshold);

  // INIT A THRESHOLD WITH NO QUANTITY
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // EXPORT TO PLT FILE
  seq->exportToVTK(plotOut, thresholdCriteria);
}

// ============================
// EVAL VARIOUS VORTEX CRITERIA
// ============================
void evalVortexCriteria(MRIOptions* opts){

  int vortexCrit = 0;
  
  // CREATE NEW SEQUENCES
  MRISequence* seq = new MRISequence(false/*Cyclic Sequence*/);

  // ADD FILE TO SEQUENCE
  seq->readPLTFile(opts->inputFileName,true);

  // EVAL VORTEX CRITERIA
  if(vortexCrit == 0){
    seq->getScan(0)->evalVortexCriteria(opts->thresholdCriteria);
  }else if(vortexCrit == 1){
    seq->getScan(0)->evalVorticity(opts->thresholdCriteria);
  }else if(vortexCrit == 2){
    seq->getScan(0)->evalEnstrophy(opts->thresholdCriteria);
  }else{
    throw MRIException("ERROR: Invalid Vortex Criterion.\n");
  }

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // WRITE OUTPUT FILES TO VTK
  seq->exportToVTK(opts->outputFileName,thresholdCriteria);
}

// ===========================================
// WRITE SPATIAL DISTRIBUTIONS OF COEFFICIENTS
// ===========================================
void writeSpatialExpansion(MRIOptions* opts){

  // CREATE NEW SEQUENCES
  MRISequence* seq = new MRISequence(false/*Cyclic Sequence*/);

  // READ FROM EXPANSION COEFFICIENTS
  seq->readFromExpansionFile(opts->inputFileName,false,kSoftThreshold,0.0);

  // SPATIALLY EVALUATE VORTEX COEFFICIENTS
  seq->getScan(0)->evalSMPVortexCriteria(seq->getScan(0)->expansion);

  // Init a Threshold with No quantity
  MRIThresholdCriteria* thresholdCriteria = new MRIThresholdCriteria(kNoQuantity,kCriterionLessThen,0.0);

  // EXPORT TO VTK
  seq->exportToVTK(opts->outputFileName,thresholdCriteria);
}

// ===============================
// RUN APPLICATION IN NORMAL MODE
// ===============================
void runApplication(MRIOptions* opts, MRICommunicator* comm){

  MRISequence* MyMRISequence;

  // MASTER PROCESSOR DOES THE READING
  if(comm->currProc == 0){

      // INIT SEQUENCE
      MyMRISequence = new MRISequence(true/*Cyclic Sequence*/);

      // LOOP ON THE NUMBER OF SCANS
      for(size_t loopA=0;loopA<opts->sequenceFileList.size();loopA++){

        // CREATE NEW SCAN
        MRIStructuredScan* MyMRIScan = new MRIStructuredScan(opts->sequenceFileTimes[loopA]);

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
        MyMRISequence->addScan(MyMRIScan);
      }
  }else{
    // LOOP ON THE NUMBER OF SCANS
      // INIT SEQUENCE
      MyMRISequence = new MRISequence(false/*Cyclic Sequence*/);
//    for(size_t loopA=0;loopA<opts->sequenceFileList.size();loopA++){
      // CREATE EMPTY SCAN
      MRIStructuredScan* MyMRIScan = new MRIStructuredScan(0.0);
      // ADD TO SEQUENCE
      MyMRISequence->addScan(MyMRIScan);
//    }
  }

  // Compute the topology of the sequence
  // The topology of the sequence must be the same
  MyMRISequence->createTopology();

  // All processes are waiting for the root to read the files
  int mpiError = MPI_Barrier(comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);

  // Scale Model if required
  if(comm->currProc == 0){
    if (opts->scaleVelocities){
      MyMRISequence->scaleVelocities(opts->scaleVelocityFactor);
    }
    if (opts->scalePositions){
      MyMRISequence->scalePositions(opts->scalePositionFactor);
    }
  }

  // Distribute Sequence Data using MPI
  if(comm->totProc > 1){
    MyMRISequence->distributeSequenceData(comm);
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
      MyMRISequence->ApplyMedianFilter(kQtyVelocityX,opts->filterNumIterations,opts->medianFilterOrder,opts->medianFilterType,opts->thresholdCriteria);
      MyMRISequence->ApplyMedianFilter(kQtyVelocityY,opts->filterNumIterations,opts->medianFilterOrder,opts->medianFilterType,opts->thresholdCriteria);
      MyMRISequence->ApplyMedianFilter(kQtyVelocityZ,opts->filterNumIterations,opts->medianFilterOrder,opts->medianFilterType,opts->thresholdCriteria);
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
      MyMRISequence->ExportForPoisson(opts->poissonFileName,opts->density,opts->viscosity,opts->thresholdCriteria,
                                      opts->PPE_IncludeAccelerationTerm,
                                      opts->PPE_IncludeAdvectionTerm,
                                      opts->PPE_IncludeDiffusionTerm,
                                      opts->PPE_IncludeReynoldsTerm,
                                      opts->readMuTFromFile,
                                      opts->muTFile,
                                      opts->smagorinskyCoeff);
    }
  }

  // SAVE FILE FOR DISTANCE COMPUTATION
  if(comm->currProc == 0){
    if (opts->exportToDistance){
      MyMRISequence->ExportForDistancing(opts->distanceFileName,opts->thresholdCriteria);
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

