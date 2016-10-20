#include "mriOptions.h"
#include "mriThresholdCriteria.h"
#include "mriConstants.h"
#include "mriUtils.h"
#include "mriException.h"

#include <boost/algorithm/string.hpp>

#include <stdlib.h>
#include <getopt.h>

using namespace boost::algorithm;


MRIOptions::MRIOptions(){
  // Set Default Values of the parameters
  runMode = rmHELP;
  // File Names
  inputFileName = "";
  outputFileName = "";
  statFileName = "";
  generateCommandFile = false;
  useCommandFile = false;
  commandFileName = "";
  // Parameters
  itTol = 1.0e-3;
  maxIt = 2000;
  thresholdType = kNoQuantity;
  thresholdValue = 0.0;
  // Save Initial Velocities
  saveInitialVel = false;
  saveExpansionCoeffs = false;
  // Apply Filter
  applySMPFilter = false;
  applyBCFilter = false;
  useConstantPatterns = false;
  // Export Format Type
  inputFormatType = itFILEVTK;
  outputFormatType = otFILEVTK;
  // Default: Process Single Scan
  haveSequence = false;
  sequenceFileName = "";
  // Post processing
  evalPopVortexCriteria = false;
  evalSMPVortexCriterion = false;
  evalPressure = false;
  // Pressure Gradient Components To be considered for pressure computation
  PPE_IncludeAccelerationTerm = true;
  PPE_IncludeAdvectionTerm = true;
  PPE_IncludeDiffusionTerm = true;
  PPE_IncludeReynoldsTerm = false;
  // Turbulent Viscosity
  readMuTFromFile = false;
  muTFile = "";
  smagorinskyCoeff = 0.15;
  // Export to Poisson
  exportToPoisson = false;
  poissonFileName = "solution.vtk";
  // Export to wall distance solver
  exportToDistance = false;
  distanceFileName = "solution.vtk";
  // Add Noise
  applyNoise = false;
  noiseIntensity = 0.0;
  // Init with Water Density and Viscosity
  density = 1000.0;
  viscosity = 1.0e-3;
  // Init with No Median Filter
  applyMedianFilter = false;
  filterNumIterations = 1;
  medianFilterType = kMedianFilter;
  medianFilterOrder = 1;
  // Clean Velocity Components
  cleanBoundaryVelocities = false;
  interpolateBoundaryVelocities = false;
  // Model Scaling
  scaleVelocities = false;
  scaleVelocityFactor = 1.0;
  scaleVelocities = false;
  scalePositionFactor = 1.0;
}

MRIOptions::~MRIOptions(){

}

// =========================
// GET COMMANDLINE ARGUMENTS
// =========================
int MRIOptions::getCommadLineOptions(int argc, char **argv){
  int c;
  string threshString;
  printf("--- COMMAND OPTIONS ECHO\n");
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
      {"saveExpansionCoeffs",    no_argument, 0, 15},
      {"poisson",    no_argument, 0, 16},
      {"normal",        no_argument, 0, 22},
      {"writexp",   no_argument, 0, 23},
      {"test",      no_argument, 0, 24},
      {"testMPI",   no_argument, 0, 25},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;    

    c = getopt_long (argc, argv, "c:i:o:g:",
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
        useCommandFile = true;
        commandFileName = optarg;
        printf("COMMAND FILE MODE\n");
        printf("Command File: %s\n",commandFileName.c_str());
        break;
      case 'g':
        // Set input file name
        generateCommandFile = true;
        commandFileName = optarg;
        printf("COMMAND FILE GENERATION MODE\n");
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
        threshString = MRIUtils::getThresholdQtyString(thresholdQty);
        printf("Threshold quantity set: %s\n",threshString.c_str());
        break;
      case 6:
        thresholdType = atoi(optarg);
        threshString = MRIUtils::getThresholdTypeString(thresholdType);
        printf("Threshold type set: %s\n",threshString.c_str());
        break;
      case 7:
        thresholdValue = atof(optarg);
        printf("Threshold value: %f\n",thresholdValue);
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
      case 15:
        saveExpansionCoeffs = true;
        printf("Saving Expansion Coefficients\n");
        break;
      case 16:
        exportToPoisson = true;
        printf("Exporting FE Files for Poisson Solver\n");
        break;
      case 22:
        runMode = rmNORMAL;
        // Read Sequence file names
        printf("Running in NORMAL mode\n");
        break;
      case '?':
        /* getopt_long already printed an error message. */
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
  printf("---\n");
  printf("\n");
  return 0;
}

// ============================
// READ FILE SEQUENCE AND TIMES
// ============================
void MRIOptions::readSequenceFileList(string fileName,MRIStringVec& sequenceFileList,MRIDoubleVec& sequenceFileTimes){
  vector<string> tokenizedString;

  // Clear Outputs
  sequenceFileList.clear();
  sequenceFileTimes.clear();

  // Read Data From File
  std::string buffer;
  std::ifstream infile;
  infile.open(fileName.c_str());
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    // Add Sequence File and time
    try{
      sequenceFileList.push_back(tokenizedString[0]);
      sequenceFileTimes.push_back(atof(tokenizedString[1].c_str()));
    }catch(...){
      infile.close();
      throw MRIException("ERROR: Cannot read sequence file.\n");
    }
  }
  // Close File
  infile.close();
}

// ======================================
// FINALIZE OPTIONS: PROCESS INPUT PARAMS
// ======================================
void MRIOptions::finalize(){
  // CREATE THRESHOLD OBJECT
  thresholdCriteria = new MRIThresholdCriteria(thresholdQty,thresholdType,thresholdValue);
  // FILL FILE LIST SEQUENCE WITH SINGLE FILE
  if(!haveSequence){
    sequenceFileList.push_back(inputFileName);
    sequenceFileTimes.push_back(0.0);
  }else{
    // READ FROM SEQUENCE LIST FILE
    readSequenceFileList(sequenceFileName,sequenceFileList,sequenceFileTimes);
  }
}

// =============================
// GET OPTIONS FROM COMMAND FILE
// =============================
int MRIOptions::getOptionsFromCommandFile(string commandFile){
  // Write Message
  printf("Reading Command file: %s\n",commandFile.c_str());

  // Declare input File
  std::ifstream infile;
  infile.open(commandFile.c_str());

  // Declare
  vector<string> tokenizedString;

  // Read Data From File
  std::string buffer;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(":,"), boost::token_compress_on);
    // Trim All Strings
    for(size_t loopA=0;loopA<tokenizedString.size();loopA++){
      trim(tokenizedString[loopA]);
    }
    // CHECK THE ELEMENT TYPE
    if(boost::to_upper_copy(tokenizedString[0]) == std::string("RUNMODE")){
      // READ RUN MODE
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("NORMAL")){
        runMode =rmNORMAL;
      }else{
        throw MRIException("ERROR: Invalid Run Mode.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("INPUTFILE")){
      try{
        inputFileName = tokenizedString[1];
      }catch(...){
        throw MRIException("ERROR: Invalid Input File.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("OUTPUTFILE")){
      try{
        outputFileName = tokenizedString[1];
      }catch(...){
        throw MRIException("ERROR: Invalid Output File.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("STATFILE")){
      try{
        statFileName = tokenizedString[1];
      }catch(...){
        throw MRIException("ERROR: Invalid Statistics File.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("DENSITY")){
        try{
          density = atof(tokenizedString[1].c_str());
        }catch(...){
          throw MRIException("ERROR: Invalid Density Value.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("APPLYSMOOTHINGFILTER")){
        
        try{
          applyMedianFilter = true;
          filterNumIterations = atoi(tokenizedString[1].c_str());
          if(boost::to_upper_copy(tokenizedString[2]) == std::string("MEDIAN")){
            medianFilterType = kMedianFilter;
          }else if(boost::to_upper_copy(tokenizedString[2]) == std::string("MEAN")){
            medianFilterType = kMeanFilter;
          }else if(boost::to_upper_copy(tokenizedString[2]) == std::string("GAUSSIAN")){
            medianFilterType = kGaussianFilter;
          }else{
            throw MRIException("ERROR: Invalid Filter Type.\n");
          }
          medianFilterOrder = atoi(tokenizedString[3].c_str());
        }catch(...){
          throw MRIException("ERROR: Invalid Median Filter Entry.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("VISCOSITY")){
        try{
          viscosity = atof(tokenizedString[1].c_str());
        }catch(...){
          throw MRIException("ERROR: Invalid Density Value.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SMPITERATIONTOLERANCE")){
      try{
        itTol = atof(tokenizedString[1].c_str());
      }catch(...){
        throw MRIException("ERROR: Invalid SMP Tolerance.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SMPMAXITERATIONS")){
      try{
        maxIt = atoi(tokenizedString[1].c_str());
      }catch(...){
        throw MRIException("ERROR: Invalid Max number of SMP Iterations.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("THRESHOLDQTY")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("POSX")){
        thresholdQty = kQtyPositionX;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("POSY")){
        thresholdQty = kQtyPositionY;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("POSZ")){
        thresholdQty = kQtyPositionZ;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("CONCENTRATION")){
        thresholdQty = kQtyConcentration;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("VELX")){
        thresholdQty = kQtyVelocityX;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("VELY")){
        thresholdQty = kQtyVelocityY;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("VELZ")){
        thresholdQty = kQtyVelocityZ;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("VELMOD")){
        thresholdQty = kQtyVelModule;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("NONE")){
          thresholdQty = kNoQuantity;
      }else{
        throw MRIException("ERROR: Invalid Threshold Quantity.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("THRESHOLDTYPE")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("LT")){
        thresholdType = kCriterionLessThen;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("GT")){
        thresholdType = kCriterionGreaterThen;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("ABSLT")){
        thresholdType = kCriterionABSLessThen;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("ABSGT")){
        thresholdType = kCriterionABSGreaterThen;
      }else{
        throw MRIException("ERROR: Invalid Threshold Type.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("THRESHOLDVALUE")){
      try{
        thresholdValue = atof(tokenizedString[1].c_str());
      }catch(...){
        throw MRIException("ERROR: Invalid Threshold Value.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SAVEINITIALVELOCITIES")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        saveInitialVel = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        saveInitialVel = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for saveInitialVel.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SAVEEXPANSIONCOEFFS")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        saveExpansionCoeffs = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        saveExpansionCoeffs = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for saveExpansionCoeffs.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("INPUTTYPE")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("VTK")){
        inputFormatType = itFILEVTK;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("TECPLOT")){
        inputFormatType = itFILETECPLOT;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("TEMPLATE")){
        inputFormatType = itTEMPLATE;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("EXPANSION")){
        inputFormatType = itEXPANSION;
      }else{
        throw MRIException("ERROR: Invalid input file type.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("OUTPUTTYPE")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("VTK")){
        outputFormatType = otFILEVTK;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("TECPLOT")){
        outputFormatType = otFILETECPLOT;
      }else{
        throw MRIException("ERROR: Invalid output file type.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("USESMPFILTER")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        applySMPFilter = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        applySMPFilter = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for applySMPFilter.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("USEBCFILTER")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        applyBCFilter = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        applyBCFilter = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for applyBCFilter.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("CLEANBOUNDARYVELOCITIES")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        cleanBoundaryVelocities = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        cleanBoundaryVelocities = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for CLEANBOUNDARYVELOCITIES.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("INTERPOLATEBOUNDARYVELOCITY")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        interpolateBoundaryVelocities = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        interpolateBoundaryVelocities = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for CLEANBOUNDARYVELOCITIES.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("USECONSTANTPATTERNS")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        useConstantPatterns = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        useConstantPatterns = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for useConstantPatterns.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SEQUENCEFILENAME")){
      try{
        haveSequence = true;
        sequenceFileName = tokenizedString[1];
      }catch(...){
        throw MRIException("ERROR: Invalid Sequence File Name.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("EVALVORTEXCRITERIA")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        evalPopVortexCriteria = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        evalPopVortexCriteria = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for evalPopVortexCriteria.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("EVALSMPVORTEXCRITERIA")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        evalSMPVortexCriterion = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        evalSMPVortexCriterion = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for evalSMPVortexCriterion.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("EVALPRESSURE")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        evalPressure = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        evalPressure = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for evalPressure.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("EXPORTTOPOISSON")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        exportToPoisson = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        exportToPoisson = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for exportToPoisson.\n");
      }

    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("EXPORTTODISTANCE")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        exportToDistance = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        exportToDistance = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for exportToDistance.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("DISTANCEFILE")){
      distanceFileName = tokenizedString[1];
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("POISSONFILE")){
      poissonFileName = tokenizedString[1];
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("ADDNOISE")){
      applyNoise = true;
      try{
        noiseIntensity = atof(tokenizedString[1].c_str());
      }catch(...){
        throw MRIException("ERROR: Invalid Noise Intensity.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("TEMPLATETYPE")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("ZEROVELOCITY")){
        templateType = kZeroVelocity;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("CONSTANT")){
        templateType = kConstantFlow;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("POISEILLE")){
        templateType = kPoiseilleFlow;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("STAGNATION")){
        templateType = kStagnationFlow;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("CYLINDRICALVOLTEX")){
        templateType = kCylindricalVortex;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("SPHERICALVORTEX")){
        templateType = kSphericalVortex;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("TOROIDALVORTEX")){
          templateType = kToroidalVortex;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRANSIENT")){
          templateType = kTransientFlow;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("CONSTANTWITHSTEP")){
          templateType = kConstantFlowWithStep;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("TAYLORVORTEX")){
          templateType = kTaylorVortex;
      }else{
        throw MRIException("ERROR: Invalid logical value for exportToPoisson.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("TEMPLATEPARAMS")){
      try{
        for(int loopA=1;loopA<tokenizedString.size();loopA++){
          templateParams.push_back(atof(tokenizedString[loopA].c_str()));
        }
      }catch(...){
        throw MRIException("ERROR: Invalid template parameters.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("USETURBVISCOSITYFROMFILE")){
      if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
        readMuTFromFile = true;
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
        readMuTFromFile = false;
      }else{
        throw MRIException("ERROR: Invalid logical value for USETURBVISCOSITYFROMFILE.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("TURBVISCOSITYFILE")){
      muTFile = tokenizedString[1];
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SMAGORINSKYCONSTANT")){
      try{
        smagorinskyCoeff = atof(tokenizedString[1].c_str());
      }catch(...){
        throw MRIException("ERROR: Invalid Smagorinsky Constant.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SCALEVELOCITY")){
        try{
          scaleVelocities = true;
          scaleVelocityFactor = atof(tokenizedString[1].c_str());
        }catch(...){
          throw MRIException("ERROR: Invalid velocity scaling parameters.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("SCALEPOSITION")){
        try{
          scalePositions = true;
          scalePositionFactor = atof(tokenizedString[1].c_str());
        }catch(...){
          throw MRIException("ERROR: Invalid position scaling parameters.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("PRESSUREGRADIENTCOMPONENTS")){
        try{
          // ACCELERATION TERM
          if(tokenizedString[1].compare("Y") == 0){
            PPE_IncludeAccelerationTerm = true;
          }else if(tokenizedString[1].compare("N") == 0){
            PPE_IncludeAccelerationTerm = false;
          }else{
            throw MRIException("ERROR: Invalid token for acceleration term.\n");
          }
          // ADVECTION TERM
          if(tokenizedString[2].compare("Y") == 0){
            PPE_IncludeAdvectionTerm = true;
          }else if(tokenizedString[2].compare("N") == 0){
            PPE_IncludeAdvectionTerm = false;
          }else{
            throw MRIException("ERROR: Invalid token for advection term.\n");
          }
          // DIFFUSION TERM
          if(tokenizedString[3].compare("Y") == 0){
            PPE_IncludeDiffusionTerm = true;
          }else if(tokenizedString[3].compare("N") == 0){
            PPE_IncludeDiffusionTerm = false;
          }else{
            throw MRIException("ERROR: Invalid token for diffusion term.\n");
          }
          // REYNOLDS TERM
          if(tokenizedString[4].compare("Y") == 0){
            PPE_IncludeReynoldsTerm = true;
          }else if(tokenizedString[4].compare("N") == 0){
            PPE_IncludeReynoldsTerm = false;
          }else{
            throw MRIException("ERROR: Invalid token for Reynolds term.\n");
          }
        }catch(...){
          throw MRIException("ERROR: Invalid pressure gradient component inclusion command.\n");
        }
    }else if((tokenizedString[0].empty())||(tokenizedString[0].at(0) == '#')){
      // Comment: Do Nothing
    }else{
      string errorMsg("ERROR: Invalid Token in input File: " + tokenizedString[0] + "\n");
      throw MRIException(errorMsg.c_str());
    }
  }
  // Close File
  infile.close();
  // Return
  return 0;
}

// UTILITY TO PRINT BOOLEAN VARIABLES
string getTrueFalseString(bool value){
  if(value){
    return string("TRUE");
  }else{
    return string("FALSE");
  }
}

// ============================
// WRITE COMMAND FILE PROTOTYPE
// ============================
int MRIOptions::writeCommandFilePrototype(string commandFile){
  // Open Output File
  FILE* f;
  f = fopen(commandFile.c_str(),"w");
  // Run Mode
  fprintf(f,"RUNMODE:");
  switch(runMode){
    case rmNORMAL:
      fprintf(f,"NORMAL\n");
      break;
    case rmEVALSEQUENCEPRESSURE:
      fprintf(f,"EVALSEQUENCEPRESSURE\n");
      break;
    case rmEVALPRESSUREFROMSIGNATUREFLOW:
      fprintf(f,"EVALPRESSUREFROMSIGNATUREFLOW\n");
      break;
    case rmPLTTOVTK:
      fprintf(f,"PLTTOVTK\n");
      break;
    case rmEVALSCANSTATISTICS:
      fprintf(f,"EVALSCANSTATISTICS\n");
      break;
    case rmCOMUTESCANMATRICES:
      fprintf(f,"COMUTESCANMATRICES\n");
      break;
    case rmPERFORMRANDOMTEST:
      fprintf(f,"PERFORMRANDOMTEST\n");
      break;
    case rmCROPANDCOMPUTEVOLUME:
      fprintf(f,"CROPANDCOMPUTEVOLUME\n");
      break;
    case rmSTREAMLINETEST1:
      fprintf(f,"STREAMLINETEST1\n");
      break;
    case rmSTREAMLINETEST2:
      fprintf(f,"STREAMLINETEST2\n");
      break;
    case rmPRINTTHRESHOLDINGTOVTK:
      fprintf(f,"PRINTTHRESHOLDINGTOVTK\n");
      break;
    case rmEVALREYNOLDSSTRESSES:
      fprintf(f,"EVALREYNOLDSSTRESSES\n");
      break;
    case rmSHOWFACEFLUXPATTERS:
      fprintf(f,"SHOWFACEFLUXPATTERS\n");
      break;
    case rmBUILDFROMCOEFFICIENTS:
      fprintf(f,"BUILDFROMCOEFFICIENTS\n");
      break;
    case rmEVALPRESSURE:
      fprintf(f,"EVALPRESSURE\n");
      break;
    case rmEVALCONCGRADIENT:
      fprintf(f,"EVALCONCGRADIENT\n");
      break;
    case rmEVALVORTEXCRITERIA:
      fprintf(f,"EVALVORTEXCRITERIA\n");
      break;
    case rmWRITESPATIALEXPANSION:
      fprintf(f,"WRITESPATIALEXPANSION\n");
      break;
    case rmSOLVEPOISSON:
      fprintf(f,"SOLVEPOISSON\n");
      break;
    case rmHELP:
      fprintf(f,"HELP\n");
      break;
  }
  // PRINT FILE NAMES
  fprintf(f,"INPUTFILENAME:%s\n",inputFileName.c_str());
  fprintf(f,"OUTPUTFILENAME:%s\n",outputFileName.c_str());
  fprintf(f,"STATFILENAME:%s\n",statFileName.c_str());

  // ADD NOISE
  fprintf(f,"ADDNOISE:%e\n",noiseIntensity);

  // SMP PARAMETERS
  fprintf(f,"SMPITERATIONTOLERANCE:%e\n",itTol);
  fprintf(f,"SMPMAXITERATIONS:%d\n",maxIt);
  // THRESHOLD QUANTITY
  fprintf(f,"THRESHOLDQTY:");
  switch (thresholdQty){
    case kQtyPositionX:
      fprintf(f,"POSX\n");
      break;
    case kQtyPositionY:
      fprintf(f,"POSY\n");
      break;
    case kQtyPositionZ:
      fprintf(f,"POSZ\n");
      break;
    case kQtyConcentration:
      fprintf(f,"CONCENTRATION\n");
      break;
    case kQtyVelocityX:
      fprintf(f,"VELX\n");
      break;
    case kQtyVelocityY:
      fprintf(f,"VELY\n");
      break;
    case kQtyVelocityZ:
      fprintf(f,"VELZ\n");
      break;
    case kQtyVelModule:
      fprintf(f,"VELMOD\n");
      break;
  }
  // THRESHOLD TYPE
  fprintf(f,"THRESHOLDTYPE:");
  switch(thresholdType){
    case kCriterionLessThen:
      fprintf(f,"LT\n");
      break;
    case kCriterionGreaterThen:
      fprintf(f,"GT\n");
      break;
    case kCriterionABSLessThen:
      fprintf(f,"ABSLT\n");
      break;
    case kCriterionABSGreaterThen:
      fprintf(f,"ABSGT\n");
      break;
  }
  // THRESHOLD VALUE
  fprintf(f,"THRESHOLDVALUE:%f\n",thresholdValue);
  // SAVE INITIAL VELOCITY
  fprintf(f,"SAVEINITIALVELOCITIES:%s\n",getTrueFalseString(saveInitialVel).c_str());
  // SAVE EXPANSION COEFFICIENTS
  fprintf(f,"SAVEEXPANSIONVOEFFS:%s\n",getTrueFalseString(saveExpansionCoeffs).c_str());
  // SMP FILTER
  fprintf(f,"USESMPFILTER:%s\n",getTrueFalseString(applySMPFilter).c_str());
  fprintf(f,"USEBCFILTER:%s\n",getTrueFalseString(applyBCFilter).c_str());
  fprintf(f,"USECONSTANTPATTERNS:%s\n",getTrueFalseString(useConstantPatterns).c_str());
  // INPUT AND OUTPUT FORMAT TYPE
  fprintf(f,"INPUTTYPE:");
  switch(inputFormatType){
    case itFILEVTK:
      fprintf(f,"VTK\n");
      break;
    case itFILETECPLOT:
      fprintf(f,"TECPLOT\n");
      break;
    case itTEMPLATE:
      fprintf(f,"TEMPLATE\n");
      break;
    case itEXPANSION:
      fprintf(f,"EXPANSION\n");
      break;
  }
  fprintf(f,"OUTPUTTYPE:");
  switch(outputFormatType){
    case otFILEVTK:
      fprintf(f,"VTK\n");
      break;
    case otFILETECPLOT:
      fprintf(f,"TECPLOT\n");
      break;
  }
  // TEMPLATE TO EXPORT
  fprintf(f,"TEMPLATETYPE:");
  switch(templateType){
    case kPoiseilleFlow:
      fprintf(f,"POISEILLE\n");
      break;
    case kCylindricalVortex:
      fprintf(f,"CYLVORTEX\n");
      break;
    case kSphericalVortex:
      fprintf(f,"SPHVORTEX\n");
      break;
    case kToroidalVortex:
      fprintf(f,"TORVORTEX\n");
      break;
    case kTransientFlow:
      fprintf(f,"TRANSIENT\n");
      break;
    case kTaylorVortex:
      fprintf(f,"TAYLORVORTEX\n");
      break;
  }

  // SEQUENCE FILE
  fprintf(f,"SEQUENCEFILENAME:%s\n",sequenceFileName.c_str());
  // POST PROCESSING
  fprintf(f,"EVALVORTEXCRITERIA:%s\n",getTrueFalseString(evalPopVortexCriteria).c_str());
  fprintf(f,"EVALSMPVORTEXCRITERIA:%s\n",getTrueFalseString(evalSMPVortexCriterion).c_str());
  fprintf(f,"EVALPRESSURE:%s\n",getTrueFalseString(evalPressure).c_str());
  // Export to Poisson
  fprintf(f,"EXPORTTOPOISSON:%s\n",getTrueFalseString(exportToPoisson).c_str());

  // Close Output file
  fclose(f);
}

// Distribute Program Options
void MRIOptions::DistributeProgramOptions(MRICommunicator* comm){
  // FORM INTEGER OPTIONS
  int size = 0;
  int mpiError = 0;
  int* intParams = NULL;
  double* doubleParams = NULL;

  // INTEGER OPTIONS
  size = 7;
  intParams = new int[size];
  intParams[0] = runMode;
  intParams[1] = templateType;
  intParams[2] = maxIt;
  intParams[3] = thresholdQty;
  intParams[4] = thresholdType;
  intParams[5] = inputFormatType;
  intParams[6] = outputFormatType;
  mpiError = MPI_Bcast(intParams,size,MPI_INT,0,comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);
  if(comm->currProc > 0){
    runMode = intParams[0];
    templateType = intParams[1];
    maxIt = intParams[2];
    thresholdQty = intParams[3];
    thresholdType = intParams[4];
    inputFormatType = intParams[5];
    outputFormatType = intParams[6];
  }
  delete [] intParams;

  // DOUBLE OPTIONS
  size = 3;
  doubleParams = new double[size];
  doubleParams[0] = itTol;
  doubleParams[1] = thresholdValue;
  doubleParams[2] = noiseIntensity;
  mpiError = MPI_Bcast(doubleParams,size,MPI_DOUBLE,0,comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);
  if(comm->currProc > 0){
    itTol = doubleParams[0];
    thresholdValue = doubleParams[1];
    noiseIntensity = doubleParams[2];
  }
  delete [] doubleParams;

  // BOOLEAN OPTIONS
  size = 12;
  bool* boolParams = new bool[size];
  boolParams[0] = generateCommandFile;
  boolParams[1] =  useCommandFile;
  boolParams[2] =  applyNoise;
  boolParams[3] =  saveInitialVel;
  boolParams[4] =  saveExpansionCoeffs;
  boolParams[5] =  applySMPFilter;
  boolParams[6] =  applyBCFilter;
  boolParams[7] =  useConstantPatterns;
  boolParams[8] =  evalPopVortexCriteria;
  boolParams[9] =  evalSMPVortexCriterion;
  boolParams[10] =  evalPressure;
  boolParams[11] =  exportToPoisson;
  mpiError = MPI_Bcast(boolParams,size,MPI::BOOL,0,comm->mpiComm);
  MRIUtils::checkMpiError(mpiError);
  if(comm->currProc > 0){
    generateCommandFile = boolParams[0];
    useCommandFile = boolParams[1];
    applyNoise = boolParams[2];
    saveInitialVel = boolParams[3];
    saveExpansionCoeffs = boolParams[4];
    applySMPFilter = boolParams[5];
    applyBCFilter = boolParams[6];
    useConstantPatterns = boolParams[7];
    evalPopVortexCriteria = boolParams[8];
    evalSMPVortexCriterion = boolParams[9];
    evalPressure = boolParams[10];
    exportToPoisson = boolParams[11];
  }
  delete [] boolParams;

  // PASS STRINGS
  comm->passString(inputFileName);
  comm->passString(outputFileName);
  comm->passString(statFileName);
  comm->passString(commandFileName);
  comm->passString(sequenceFileName);

  // PASS VECTOR: CAREFULL DO NOT PASS ZERO DIM OBJECTS
  if(templateParams.size() > 0){
    comm->passStdDoubleVector(templateParams);
  }

  // PASS STRING LIST
  //vector<string> sequenceFileList;

}

// Write Options to file
int MRIOptions::writeOptionsToFile(string outFile){
  // Open Output File
  FILE* f;
  f = fopen(outFile.c_str(),"w");
  // Write Options
  fprintf(f,"Run Mode: %d\n",runMode);
  // File Names
  fprintf(f,"Input File: %s\n",inputFileName.c_str());
  fprintf(f,"Output File: %s\n",outputFileName.c_str());
  fprintf(f,"Statistics File: %s\n",statFileName.c_str());
  fprintf(f,"Generate Command File: %s\n",generateCommandFile ? "True" : "False");
  fprintf(f,"Use Command File: %s\n",useCommandFile ? "True" : "False");
  fprintf(f,"Command File Name: %s\n",commandFileName.c_str());
  // Parameters
  fprintf(f,"Iteration Tolerance: %e\n",itTol);
  fprintf(f,"Max Iterations: %d\n",maxIt);
  fprintf(f,"Threshold Type: %d\n",thresholdType);
  fprintf(f,"Threshold Value: %e\n",thresholdValue);
  // Save Initial Velocities
  fprintf(f,"Set Initial Velocities: %s\n",saveInitialVel ? "True" : "False");
  fprintf(f,"Save Expansion Coefficients: %s\n",saveExpansionCoeffs ? "True" : "False");
  // Apply Filter
  fprintf(f,"Apply SMP Filter: %s\n",applySMPFilter ? "True" : "False");
  fprintf(f,"Apply Boundary Filter: %s\n",applyBCFilter ? "True" : "False");
  fprintf(f,"Use Constant Pattern: %s\n",useConstantPatterns ? "True" : "False");
  // Export Format Type
  fprintf(f,"Input Format Type: %d\n",inputFormatType);
  fprintf(f,"Output Format Type: %d\n",outputFormatType);
  // Default: Process Single Scan
  fprintf(f,"Sequence File Name: %s\n",sequenceFileName.c_str());
  // Post processing
  fprintf(f,"Eval Popular Vortex Criteria: %s\n",evalPopVortexCriteria ? "True" : "False");
  fprintf(f,"Eval SMP Vortex Criteria: %s\n",evalSMPVortexCriterion ? "True" : "False");
  fprintf(f,"Eval Pressure: %s\n",evalPressure ? "True" : "False");
  // Export to Poisson
  fprintf(f,"Export To FEM Poisson Solver: %s\n",exportToPoisson ? "True" : "False");
  // Add Noise
  fprintf(f,"Apply Noise: %s\n",applyNoise ? "True" : "False");
  fprintf(f,"Noise Intensity: %e\n",noiseIntensity);
  // Close Output file
  fclose(f);

}





