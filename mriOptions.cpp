#include "mriOptions.h"
#include "mriThresholdCriteria.h"
#include "mriConstants.h"

#include <stdlib.h>
#include <getopt.h>

MRIOptions::MRIOptions(){
  // Set Default Values of the parameters
  runMode = rmHELP;
  // File Names
  inputFileName = "";
  outputFileName = "";
  statFileName = "";
  commandFileName = "";
  // Parameters
  itTol = 1.0e-3;
  maxIt = 2000;
  thresholdType = kNoQuantity;
  thresholdValue = 0.0;
  // Save Initial Velocities
  saveInitialVel = false;
  // Apply Filter
  applySMPFilter = false;
  applyBCFilter = false;
  useConstantPatterns = true;
  // Export Format Type
  inputFormatType = itFILEVTK;
  outputFormatType = otFILEVTK;
  // Default: Process Single Scan
  sequenceFileName = "";
  // Post processing
  evalPopVortexCriteria = false;
  evalSMPVortexCriterion = false;
  evalPressure = false;
}

// =========================
// GET COMMANDLINE ARGUMENTS
// =========================
int MRIOptions::getCommadLineOptions(int argc, char **argv){
  int c;
  printf("--- COMMAND OPTIONS ECHO\n");
  printf("\n");
  while (1){
    static struct option long_options[] =
    { {"input",     required_argument, 0, 'i'},
      {"output",    required_argument, 0, 'o'},
      {"command",   required_argument, 0, 'c'},
      {"filter",        no_argument, 0, 1},
      {"bcfilter",        no_argument, 0, 2},
      {"iterationTol",  required_argument, 0, 3},
      {"maxIterations", required_argument, 0, 4},
      {"thresholdQty",  required_argument, 0, 5},
      {"thresholdType", required_argument, 0, 6},
      {"thresholdValue",required_argument, 0, 7},
      {"iformat",    required_argument, 0, 8},
      {"oformat",    required_argument, 0, 9},
      {"sequence",    required_argument, 0, 10},
      {"noConstPatterns",    no_argument, 0, 11},
      {"evalPOPVortex",    no_argument, 0, 12},
      {"evalSMPVortex",    no_argument, 0, 13},
      {"saveInitialVel",    no_argument, 0, 14},

      {"normal",        no_argument, 0, 22},
      {"writexp",   no_argument, 0, 23},
      {"test",      no_argument, 0, 24},
      {"testMPI",   no_argument, 0, 25},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;    

    c = getopt_long (argc, argv, "c:i:o:",
                     long_options, &option_index);

    /* Detect the end of the options.*/
    if (c == -1){
      break;
    }    
    // Check which option was selected
    switch (c){
      case 0:
        /* If this option sets a flag, do nothing else now. */
        if (long_options[option_index].flag != 0){
          break;
        }
        printf ("option %s", long_options[option_index].name);
        if (optarg){
          printf (" with arg %s", optarg);
        }
        printf ("\n");
        break;
      case 'i':
        // Set input file name
        inputFileName = optarg;
        printf("Input File: %s\n",inputFileName.c_str());
        break;
      case 'o':
        // Set input file name
        outputFileName = optarg;
        printf("Output File: %s\n",outputFileName.c_str());
        break;
      case 'c':
        // Set input file name
        commandFileName = optarg;
        printf("COMMAND FILE MODE\n");
        printf("Command File: %s\n",commandFileName.c_str());
        break;
      case 1:
        applySMPFilter = true;
        printf("SMP Filter Enabled\n");
        break;
      case 2:
        applyBCFilter = true;
        printf("Boundary Filter Enabled\n");
        break;
      case 3:
        itTol = atof(optarg);
        printf("Iteration Tolerance set to %e\n",itTol);
        break;
      case 4:
        maxIt = atoi(optarg);
        printf("Maximum number of iterations set to %d\n",maxIt);
        break;
      case 5:
        thresholdQty = atoi(optarg);
        break;
      case 6:
        thresholdType = atoi(optarg);
        break;
      case 7:
        thresholdValue = atof(optarg);
        break;
      case 8:
        inputFormatType = atoi(optarg);
        printf("Input File Format: %s\n",optarg);
        break;
      case 9:
        outputFormatType = atoi(optarg);
        printf("Output File Format: %s\n",optarg);
        break;
      case 10:
        sequenceFileName = optarg;
        // Read Sequence file names
        printf("Output File Format: %s\n",optarg);
        break;
      case 11:
        useConstantPatterns = false;
        // Read Sequence file names
        printf("Constant Patterns Disabled in BC Filter\n");
        break;
      case 12:
        evalPopVortexCriteria = true;
        printf("Evaluating Vortex Criteria\n");
        break;
      case 13:
        evalSMPVortexCriterion = true;
        printf("Evaluating SMP Vortex Criteria\n");
        break;
      case 14:
        saveInitialVel = true;
        printf("Saving Initial Velocities To Output\n");
        break;
      case 22:
        runMode = rmNORMAL;
        // Read Sequence file names
        printf("Running in NORMAL mode\n");
        break;
      case '?':
        /* getopt_long already printed an error message. */
        break;
      case 16:
        runMode = rmSOLVEPOISSON;
        break;
      default:
        abort ();
      }
    }
    /* Print any remaining command line arguments (not options). */
    /*if (optind < argc){
      printf ("non-option ARGV-elements: ");
      while (optind < argc){
        printf ("%s ", argv[optind++]);
      }
      putchar ('\n');
    }*/
  return 0;
}

// ======================================
// FINALIZE OPTIONS: PROCESS INPUT PARAMS
// ======================================
void MRIOptions::finalize(){
  // CREATE THRESHOLD OBJECT
  thresholdCriteria = new MRIThresholdCriteria(thresholdQty,thresholdType,thresholdValue);
  // FILL FILE LIST SEQUENCE WITH SINGLE FILE
  if(sequenceFileName == ""){
    sequenceFileList.push_back(inputFileName);
  }else{
    // READ FROM SEQUENCE LIST FILE
  }

}
