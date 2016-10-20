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

// =================================================
// READ FILES IN VARIOUS FORMATS AND DISTRIBUTE GRID
// =================================================

void readAndDistribute(MRICommunicator* comm, MRIOptions* opts, MRISequence* seq){

  // INIT SEQUENCE
  seq = new MRISequence(true/*Cyclic Sequence*/);

  // LOOP ON THE NUMBER OF SCANS
  for(size_t loopA=0;loopA<opts->sequenceFileList.size();loopA++){

    // CHOOSE INPUT FORMAT
    if(opts->inputFormatType == itTEMPLATE){
      // CREATE TEMPLATE
      seq->createSampleCase(opts->templateType,opts->templateParams);
    }else if(opts->inputFormatType == itEXPANSION){
      // READ FROM EXPANSION COEFFICIENTS
      bool applyThreshold = true;
      int thresholdType = kHardThresold;
      double thresholdRatio = 0.5;
      seq->readFromExpansionFile(opts->sequenceFileList[loopA],applyThreshold,thresholdType,thresholdRatio);
    }else if (opts->inputFormatType == itFILEVTK){
      // READ FROM FILE          
      seq->readVTKStructuredPoints(opts->sequenceFileList[loopA], true);
    }else if (opts->inputFormatType == itFILETECPLOT){
      // READ FROM FILE
      seq->readPLTFile(opts->sequenceFileList[loopA], true);
    }
  }
  
  // Compute the topology for all sequences in all processors
  seq->createTopology();
}

// ============
// WRITE OUTPUT
// ============
void writeOutput(MRICommunicator* comm, MRIOptions* opts, MRISequence* seq){
  // EXPORT FILE FROM ALL PROECESSORS IN ORDER
  if(comm->currProc == 0){
    if(opts->outputFormatType == itFILEVTK){
      // READ FROM FILE
      seq->exportToVTK(opts->outputFileName,opts->thresholdCriteria);
    }else if (opts->outputFormatType == itFILETECPLOT){
      // READ FROM FILE
      seq->exportToTECPLOT(opts->outputFileName);
    }else{
      throw MRIException("ERROR: Invalid output file format.\n");
    }
  }
}

// ===============================
// RUN APPLICATION IN NORMAL MODE
// ===============================
void runApplication(MRIOptions* opts, MRICommunicator* comm){

  // CREATE NEW SEQUENCE
  MRISequence* seq;

  // READ AND DISTRIBUTED MEASUREMENT GRID
  readAndDistribute(comm,opts,seq);

  // SYNC PROCESSES
  int mpiError = MPI_Barrier(comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);

  // PERFOM OPERATIONS ACCORDING TO LIST
  for(int loopA=0;loopA<opts->operationList.size();loopA++){
    opts->operationList[loopA]->processSequence(seq);
  }

  // WRITE TO OUTPUT
  writeOutput(comm,opts,seq);

  // FREE MEMORY
  delete seq;
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
    writeHeader();

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
          MRIUtils::writeProgramHelp();
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
        writeSchMessage(currMsgs);
        return (1);
        break;
    }
  }catch (std::exception& ex){
    if(comm->currProc == 0){
      writeSchMessage(std::string(ex.what()));
      writeSchMessage(std::string("\n"));
      writeSchMessage(std::string("Program Terminated.\n"));
    }
    // Finalize MPI
    delete options;
    MPI::Finalize();
    return -1;
  }
  if(comm->currProc == 0){
    writeSchMessage(string("\n"));
    writeSchMessage(string("Program Completed.\n"));
  }
  // Finalize MPI
  //delete comm;
  delete options;
  MPI::Finalize();
  return 0;
}

