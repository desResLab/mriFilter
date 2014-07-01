#ifndef MRIPROGRAMOPTIONS_H
#define MRIPROGRAMOPTIONS_H

#include <string>

class mriProgramOptions{
public:
  int runMode;
  std::string inputFileName;
  std::string outputFileName;
  std::string statFileName;
  double itTol;
  int maxIt;
  std::string thresholdType;
  double thresholdValue;
  // Export File Format
  int exportFormat;
  // Constructor
  mriProgramOptions();
  // Get Command Options from input Arguments
  int getCommadLineOptions(int argc, char **argv);
};

#endif // MRIPROGRAMOPTIONS_H
