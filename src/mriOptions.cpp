#include "mriOptions.h"

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
  // Export Format Type
  inputFormatType = itFILEVTK;
  outputFormatType = otFILEVTK;
  // Default: Process Single Scan
  haveSequence = false;
  sequenceFileName = "";
  // Init with Water Density and Viscosity
  density = 1000.0;
  viscosity = 1.0e-3;
  // Init Threshold Value
  thresholdQty = kNoQuantity;
  thresholdType = kCriterionLessThen;
  thresholdValue = 0.0;
  // PPE Options
  PPE_IncludeAccelerationTerm = true;
  PPE_IncludeAdvectionTerm = true;
  PPE_IncludeDiffusionTerm = true;
  PPE_IncludeReynoldsTerm = false;
  readMuTFromFile = false;
  muTFile = "";
  smagorinskyCoeff = 0.15;
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
    { {"command",   required_argument, 0, 'c'},
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
      case 'c':
        // Set input file name
        useCommandFile = true;
        commandFileName = optarg;
        printf("COMMAND FILE MODE\n");
        printf("Command File: %s\n",commandFileName.c_str());
        break;
      case '?':
        /* getopt_long already printed an error message. */
        break;
      default:
        abort ();
      }
    }
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
  // FILL FILE LIST SEQUENCE WITH SINGLE FILE
  if(!haveSequence){
    sequenceFileList.push_back(inputFileName);
    sequenceFileTimes.push_back(0.0);
  }else{
    // READ FROM SEQUENCE LIST FILE
    readSequenceFileList(sequenceFileName,sequenceFileList,sequenceFileTimes);
  }
  // Create Threshold Object
  thresholdCriteria = new MRIThresholdCriteria(thresholdQty, thresholdType, thresholdValue);
}

// =============================
// GET OPTIONS FROM COMMAND FILE
// =============================
int MRIOptions::getOptionsFromCommandFile(string commandFile){
  
  // Write Message
  printf("Reading Command file: %s\n",commandFile.c_str());

  // Local Variables
  bool applyMedianFilter;
  int filterNumIterations;
  int filterType;
  int filterOrder;

  double itTol;
  int maxIt;

  bool saveInitialVel;
  bool saveExpansionCoeffs;

  bool applySMPFilter;
  bool applyBCFilter;
  bool cleanBoundaryVelocities;
  bool interpolateBoundaryVelocities;
  bool useConstantPatterns;

  bool evalPopVortexCriteria;
  bool evalSMPVortexCriterion;
  bool evalPressure;
  bool exportToPoisson;
  bool exportToDistance;
  string distanceFileName;
  string poissonFileName;
  bool applyNoise;
  double noiseIntensity;
  double smagorinskyCoeff;
  bool scaleVelocities;
  double scaleVelocityFactor;
  bool scalePositions;
  double scalePositionFactor;

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
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("VISCOSITY")){
        try{
          viscosity = atof(tokenizedString[1].c_str());
        }catch(...){
          throw MRIException("ERROR: Invalid Density Value.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("APPLYSMOOTHINGFILTER")){
        
        try{
          applyMedianFilter = true;
          filterNumIterations = atoi(tokenizedString[1].c_str());
          if(boost::to_upper_copy(tokenizedString[2]) == std::string("MEDIAN")){
            filterType = kMedianFilter;
          }else if(boost::to_upper_copy(tokenizedString[2]) == std::string("MEAN")){
            filterType = kMeanFilter;
          }else if(boost::to_upper_copy(tokenizedString[2]) == std::string("GAUSSIAN")){
            filterType = kGaussianFilter;
          }else{
            throw MRIException("ERROR: Invalid Filter Type.\n");
          }
          filterOrder = atoi(tokenizedString[3].c_str());

          // Create New Operation
          MRIOperation* op = new MRIOpApplySmoothing(filterNumIterations,filterType,filterOrder);
          // Add to the operation list
          operationList.push_back(op);

        }catch(...){
          throw MRIException("ERROR: Invalid Median Filter Entry.\n");
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
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("PLT")){
        inputFormatType = itFILEPLT;
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
      }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("PLT")){
        outputFormatType = otFILEPLT;
      }else{
        throw MRIException("ERROR: Invalid output file type.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("USESMPFILTER")){

      try{
        // applyBCFilter
        if(boost::to_upper_copy(tokenizedString[1]) == std::string("TRUE")){
          applyBCFilter = true;
        }else if(boost::to_upper_copy(tokenizedString[1]) == std::string("FALSE")){
          applyBCFilter = false;
        }else{
          throw MRIException("ERROR: Invalid logical value for applySMPFilter.\n");
        }
        // useConstantPatterns
        if(boost::to_upper_copy(tokenizedString[2]) == std::string("TRUE")){
          useConstantPatterns = true;
        }else if(boost::to_upper_copy(tokenizedString[2]) == std::string("FALSE")){
          useConstantPatterns = false;
        }else{
          throw MRIException("ERROR: Invalid logical value for applySMPFilter.\n");
        }
        // itTol
        itTol = atof(tokenizedString[3].c_str());
        // maxIt
        maxIt = atoi(tokenizedString[4].c_str());
      }catch(...){
        throw MRIException("ERROR: Invalid Export To Poisson Command Line.\n");
      }    
      // Create Operation
      MRIOperation* op = new MRIOpApplySolenoidalFilter(applyBCFilter,useConstantPatterns,itTol,maxIt);
      // Add to the operation list
      operationList.push_back(op);
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
      try{
        poissonFileName = tokenizedString[1];
      }catch(...){
        throw MRIException("ERROR: Invalid Export To Poisson Command Line.\n");
      }    
      // Create Operation For Poisson Export
      MRIOperation* op = new MRIOpExportForPoissonSolver(poissonFileName,
                                                         density,
                                                         viscosity,
                                                         PPE_IncludeAccelerationTerm,
                                                         PPE_IncludeAdvectionTerm,
                                                         PPE_IncludeDiffusionTerm,
                                                         PPE_IncludeReynoldsTerm,
                                                         readMuTFromFile,
                                                         muTFile,
                                                         smagorinskyCoeff);
      // Add to the operation list
      operationList.push_back(op);
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
    }else if(boost::to_upper_copy(tokenizedString[0]) == std::string("ADDNOISE")){
      applyNoise = true;
      try{
        noiseIntensity = atof(tokenizedString[1].c_str());        
      }catch(...){
        throw MRIException("ERROR: Invalid Noise Intensity.\n");
      }
      // Create Operation
      MRIOperation* op = new MRIOpApplyNoise(noiseIntensity);
      // Add to the operation list
      operationList.push_back(op);
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


// Distribute Program Options
void MRIOptions::DistributeProgramOptions(MRICommunicator* comm){
  // FORM INTEGER OPTIONS
/*  int size = 0;
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
*/
}

