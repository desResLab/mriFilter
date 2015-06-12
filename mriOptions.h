#ifndef MRIPROGRAMOPTIONS_H
#define MRIPROGRAMOPTIONS_H

# include <string>
# include <vector>

# include "mriThresholdCriteria.h"

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
  const int itFILETECPLOT                   = 1;
  const int itTEMPLATE                      = 2;
  const int itEXPANSION                     = 3;

  // OUTPUT TYPES
  const int otFILEVTK                       = 0;
  const int otFILETECPLOT                   = 1;

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
  vector<double> templateParams;
  // Command File Options
  bool generateCommandFile;
  bool useCommandFile;
  string commandFileName;
  double itTol;
  int maxIt;
  int thresholdQty;
  int thresholdType;
  double thresholdValue; 
  MRIThresholdCriteria* thresholdCriteria;
  // Noise to apply
  bool applyNoise;
  double noiseIntensity;
  // Save Initial Velocities
  bool saveInitialVel;
  // Save Expansion Coefficients
  bool saveExpansionCoeffs;
  // Export File Format
  int inputFormatType;
  int outputFormatType;
  // Apply SMP Filter
  bool applySMPFilter;
  // Apply BC Filter
  bool applyBCFilter;
  bool useConstantPatterns;
  // Sequence Processing
  string sequenceFileName;
  vector<string> sequenceFileList;
  // Post processing
  bool evalPopVortexCriteria;
  bool evalSMPVortexCriterion;
  bool evalPressure;
  // Export to Poisson Solver
  bool exportToPoisson;

  // CONSTRUCTOR
  MRIOptions();

  // Get Command Options from input Arguments
  int getCommadLineOptions(int argc, char **argv);
  // Get Options from Command File
  int getOptionsFromCommandFile(string commandFile);
  // Write Command File Prototype
  int writeCommandFilePrototype(string commandFile);

  // Process Options
  void finalize();
};

#endif // MRIPROGRAMOPTIONS_H
