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

// =================================================
// READ FILES IN VARIOUS FORMATS AND DISTRIBUTE GRID
// =================================================
void readAndDistribute(MRICommunicator* comm, MRIOptions* opts, MRISequence* seq){

  // CHOOSE INPUT FORMAT
  if(opts->inputFormatType == itTEMPLATE){
    
    // CREATE TEMPLATE
    seq->createSampleCase(opts->templateType,opts->templateParams);

  }else if(opts->inputFormatType == itEXPANSION){
    
    // READ FROM EXPANSION COEFFICIENTS
    bool applyThreshold = true;
    int thresholdType = kHardThresold;
    double thresholdRatio = 0.5;
    seq->readFromExpansionFiles(opts->sequenceFileList,opts->sequenceFileTimes,    
                                applyThreshold, 
                                thresholdType,
                                thresholdRatio);

  }else if (opts->inputFormatType == itFILEVTK){

    // READ FROM FILE          
    seq->readFromASCIISequence(kInputVTK,opts->sequenceFileList,opts->sequenceFileTimes); 

  }else if (opts->inputFormatType == itFILEPLT){
    
    // READ FROM FILE
    seq->readFromASCIISequence(kInputPLT,opts->sequenceFileList,opts->sequenceFileTimes); 

  }
  
  // Compute the topology for all sequences
  seq->createTopology();
}

// ============
// WRITE OUTPUT
// ============
void writeOutput(MRICommunicator* comm, MRIOptions* opts, MRISequence* seq){  

  // EXPORT FILE FROM ALL PROECESSORS IN ORDER
  if(comm->currProc == 0){
    if(opts->outputFormatType == otFILEVTK){
      // READ FROM FILE      
      seq->exportToVTK(opts->outputFileName,opts->thresholdCriteria);
    }else if (opts->outputFormatType == otFILEPLT){
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

  // INIT SEQUENCE
  seq = new MRISequence(true/*Cyclic Sequence*/);

  // READ AND DISTRIBUTED MEASUREMENT GRID
  readAndDistribute(comm,opts,seq);

  // SYNC PROCESSES
  int mpiError = MPI_Barrier(comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);

  // PERFOM OPERATIONS ACCORDING TO LIST
  for(int loopA=0;loopA<opts->operationList.size();loopA++){
    opts->operationList[loopA]->processSequence(comm,opts->thresholdCriteria,seq);
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

