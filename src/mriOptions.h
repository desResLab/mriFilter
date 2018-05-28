#ifndef MRIPROGRAMOPTIONS_H
#define MRIPROGRAMOPTIONS_H

# include <string>
# include <vector>
# include <stdlib.h>
# include <getopt.h>
# include <boost/algorithm/string.hpp>


# include "mriThresholdCriteria.h"
# include "mriOperation.h"
# include "mriCommunicator.h"
# include "mriConstants.h"
# include "mriUtils.h"
# include "mriException.h"

using namespace std;

  // RUN MODES
  const int rmNORMAL                        = 0;
  const int rmEVALSEQUENCEPRESSURE          = 1;
  const int rmEVALPRESSUREFROMSIGNATUREFLOW = 2;
  const int rmPLTTOVTK                      = 3;
  const int rmEVALSCANSTATISTICS            = 4;
  const int rmCOMUTESCANMATRICES            = 5;
  const int rmPERFORMRANDOMTEST             = 6;
  const int rmCROPANDCOMPUTEVOLUME          = 7;
  const int rmSTREAMLINETEST1               = 8;
  const int rmSTREAMLINETEST2               = 9;
  const int rmPRINTTHRESHOLDINGTOVTK        = 10;
  const int rmEVALREYNOLDSSTRESSES          = 11;
  const int rmSHOWFACEFLUXPATTERS           = 12;
  const int rmBUILDFROMCOEFFICIENTS         = 13;
  const int rmEVALPRESSURE                  = 14;
  const int rmEVALCONCGRADIENT              = 15;
  const int rmEVALVORTEXCRITERIA            = 16;
  const int rmWRITESPATIALEXPANSION         = 17;
  const int rmSOLVEPOISSON                  = 18;
  const int rmHELP                          = 19;

  // INPUT TYPES
  const int itFILEVTK                       = 0;
  const int itFILEPLT                       = 1;
  const int itTEMPLATE                      = 2;
  const int itEXPANSION                     = 3;

  // OUTPUT TYPES
  const int otFILEVTK                       = 0;
  const int otFILEPLT                       = 1;

class MRIOperation;

// CLASS MRIPROGRAMOPTIONS
class MRIOptions{
public:
  // OPTIONS
  int runMode;
  string inputFileName;
  string outputFileName;
  string statFileName;
  // Input Templates
  int templateType;
  MRIDoubleVec templateParams;
  // Command File Options
  bool generateCommandFile;
  bool useCommandFile;
  string commandFileName;
  // Sequence Processing
  bool haveSequence;
  string sequenceFileName;
  MRIStringVec sequenceFileList;
  MRIDoubleVec sequenceFileTimes;
  // Export File Format
  int inputFormatType;
  int outputFormatType;
  // Material Properties
  double density;
  double viscosity;
  // threshold Criteria
  int thresholdQty;
  int thresholdType;
  double thresholdValue;  
  MRIThresholdCriteria* thresholdCriteria;
  // PPE Solver Options
  bool PPE_IncludeAccelerationTerm;
  bool PPE_IncludeAdvectionTerm;
  bool PPE_IncludeDiffusionTerm;
  bool PPE_IncludeReynoldsTerm;  
  bool readMuTFromFile;
  string muTFile;
  double smagorinskyCoeff;
  // List of operations
  vector<MRIOperation*> operationList;

  // CONSTRUCTOR
  MRIOptions();
  
  // Distructor
  ~MRIOptions();

  // Get Command Options from input Arguments
  int getCommadLineOptions(int argc, char **argv);
  
  // Get Options from Command File
  int getOptionsFromCommandFile(string commandFile);
  
  // Write Options to file
  int writeOptionsToFile(string outFile);
  
  // Read Sequence File List
  void readSequenceFileList(string fileName,MRIStringVec& sequenceFileList,MRIDoubleVec& sequenceFileTimes);

  // Process Options
  void finalize();

  // MESSAGE PASSING
  void DistributeProgramOptions(MRICommunicator* comm);

};

#endif // MRIPROGRAMOPTIONS_H
