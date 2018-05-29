#include "mriOptions.h"

using namespace boost::algorithm;
using namespace std;

mriOptions::mriOptions(){
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

mriOptions::~mriOptions(){

}

// =========================
// GET COMMANDLINE ARGUMENTS
// =========================
int mriOptions::getCommadLineOptions(int argc, char **argv){
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
void mriOptions::readSequenceFileList(string fileName,mriStringVec& sequenceFileList,mriDoubleVec& sequenceFileTimes){
  vector<string> tokenizedString;

  // Clear Outputs
  sequenceFileList.clear();
  sequenceFileTimes.clear();

  // Read Data From File
  string buffer;
  ifstream infile;
  infile.open(fileName.c_str());
  while (getline(infile,buffer)){
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
      throw mriException("ERROR: Cannot read sequence file.\n");
    }
  }
  // Close File
  infile.close();
}

// ======================================
// FINALIZE OPTIONS: PROCESS INPUT PARAMS
// ======================================
void mriOptions::finalize(){
  // FILL FILE LIST SEQUENCE WITH SINGLE FILE
  if(!haveSequence){
    sequenceFileList.push_back(inputFileName);
    sequenceFileTimes.push_back(0.0);
  }else{
    // READ FROM SEQUENCE LIST FILE
    readSequenceFileList(sequenceFileName,sequenceFileList,sequenceFileTimes);
  }
  // Create Threshold Object
  thresholdCriteria = new mriThresholdCriteria(thresholdQty, thresholdType, thresholdValue);
}

// =============================
// GET OPTIONS FROM COMMAND FILE
// =============================
int mriOptions::getOptionsFromCommandFile(string commandFile){
  
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
  double noiseSeed;
  double smagorinskyCoeff;
  bool scaleVelocities;
  double scaleVelocityFactor;
  bool scalePositions;
  double scalePositionFactor;

  // Declare input File
  ifstream infile;
  infile.open(commandFile.c_str());

  // Declare
  mriStringVec tokenizedString;

  // Read Data From File
  string buffer;
  while (getline(infile,buffer)){
    // printf("%s\n",buffer.c_str());
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(":,"), boost::token_compress_on);
    // Trim All Strings
    for(size_t loopA=0;loopA<tokenizedString.size();loopA++){
      trim(tokenizedString[loopA]);
    }
    // CHECK THE ELEMENT TYPE
    if(boost::to_upper_copy(tokenizedString.at(0)) == string("RUNMODE")){
      // READ RUN MODE
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("NORMAL")){
        runMode =rmNORMAL;
      }else{
        throw mriException("ERROR: Invalid Run Mode.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("INPUTFILE")){
      try{
        inputFileName = tokenizedString[1];
      }catch(...){
        throw mriException("ERROR: Invalid Input File.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("OUTPUTFILE")){
      try{
        outputFileName = tokenizedString[1];
      }catch(...){
        throw mriException("ERROR: Invalid Output File.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("STATFILE")){
      try{
        statFileName = tokenizedString[1];
      }catch(...){
        throw mriException("ERROR: Invalid Statistics File.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("DENSITY")){
        try{
          density = atof(tokenizedString[1].c_str());
        }catch(...){
          throw mriException("ERROR: Invalid Density Value.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("VISCOSITY")){
        try{
          viscosity = atof(tokenizedString[1].c_str());
        }catch(...){
          throw mriException("ERROR: Invalid Density Value.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("APPLYSMOOTHINGFILTER")){
        
        try{
          applyMedianFilter = true;
          filterNumIterations = atoi(tokenizedString.at(1).c_str());
          if(boost::to_upper_copy(tokenizedString.at(2)) == string("MEDIAN")){
            filterType = kMedianFilter;
          }else if(boost::to_upper_copy(tokenizedString.at(2)) == string("MEAN")){
            filterType = kMeanFilter;
          }else if(boost::to_upper_copy(tokenizedString.at(2)) == string("GAUSSIAN")){
            filterType = kGaussianFilter;
          }else{
            throw mriException("ERROR: Invalid Filter Type.\n");
          }
          filterOrder = atoi(tokenizedString.at(3).c_str());

          // Create New Operation
          mriOperation* op = new mriOpApplySmoothing(filterNumIterations,filterType,filterOrder);
          // Add to the operation list
          operationList.push_back(op);

        }catch(...){
          throw mriException("ERROR: Invalid Median Filter Entry.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("SMPITERATIONTOLERANCE")){
      try{
        itTol = atof(tokenizedString.at(1).c_str());
      }catch(...){
        throw mriException("ERROR: Invalid SMP Tolerance.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("SMPMAXITERATIONS")){
      try{
        maxIt = atoi(tokenizedString.at(1).c_str());
      }catch(...){
        throw mriException("ERROR: Invalid Max number of SMP Iterations.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("THRESHOLDQTY")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("POSX")){
        thresholdQty = kQtyPositionX;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("POSY")){
        thresholdQty = kQtyPositionY;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("POSZ")){
        thresholdQty = kQtyPositionZ;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("CONCENTRATION")){
        thresholdQty = kQtyConcentration;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("VELX")){
        thresholdQty = kQtyVelocityX;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("VELY")){
        thresholdQty = kQtyVelocityY;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("VELZ")){
        thresholdQty = kQtyVelocityZ;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("VELMOD")){
        thresholdQty = kQtyVelModule;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("NONE")){
          thresholdQty = kNoQuantity;
      }else{
        throw mriException("ERROR: Invalid Threshold Quantity.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("THRESHOLDTYPE")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("LT")){
        thresholdType = kCriterionLessThen;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("GT")){
        thresholdType = kCriterionGreaterThen;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("ABSLT")){
        thresholdType = kCriterionABSLessThen;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("ABSGT")){
        thresholdType = kCriterionABSGreaterThen;
      }else{
        throw mriException("ERROR: Invalid Threshold Type.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("THRESHOLDVALUE")){
      try{
        thresholdValue = atof(tokenizedString.at(1).c_str());
      }catch(...){
        throw mriException("ERROR: Invalid Threshold Value.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("SAVEINITIALVELOCITIES")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        saveInitialVel = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        saveInitialVel = false;
      }else{
        throw mriException("ERROR: Invalid logical value for saveInitialVel.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("SAVEEXPANSIONCOEFFS")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        saveExpansionCoeffs = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        saveExpansionCoeffs = false;
      }else{
        throw mriException("ERROR: Invalid logical value for saveExpansionCoeffs.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("INPUTTYPE")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("VTK")){
        inputFormatType = itFILEVTK;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("PLT")){
        inputFormatType = itFILEPLT;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("TEMPLATE")){
        inputFormatType = itTEMPLATE;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("EXPANSION")){
        inputFormatType = itEXPANSION;
      }else{
        throw mriException("ERROR: Invalid input file type.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("OUTPUTTYPE")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("VTK")){
        outputFormatType = otFILEVTK;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("PLT")){
        outputFormatType = otFILEPLT;
      }else{
        throw mriException("ERROR: Invalid output file type.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("USESMPFILTER")){

      try{
        // applyBCFilter
        if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
          applyBCFilter = true;
        }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
          applyBCFilter = false;
        }else{
          throw mriException("ERROR: Invalid logical value (applyBCFilter) for USESMPFILTER.\n");
        }
        // useConstantPatterns
        if(boost::to_upper_copy(tokenizedString.at(2)) == string("TRUE")){
          useConstantPatterns = true;
        }else if(boost::to_upper_copy(tokenizedString.at(2)) == string("FALSE")){
          useConstantPatterns = false;
        }else{
          throw mriException("ERROR: Invalid logical value (useConstantPatterns) for USESMPFILTER.\n");
        }
        // itTol
        itTol = atof(tokenizedString.at(3).c_str());
        // maxIt
        maxIt = atoi(tokenizedString.at(4).c_str());
      }catch(...){
        throw mriException("ERROR: Invalid definition of USESMPFILTER.\n");
      }    
      // Create Operation
      mriOperation* op = new mriOpApplySolenoidalFilter(applyBCFilter,useConstantPatterns,itTol,maxIt);
      // Add to the operation list
      operationList.push_back(op);
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("CLEANBOUNDARYVELOCITIES")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        cleanBoundaryVelocities = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        cleanBoundaryVelocities = false;
      }else{
        throw mriException("ERROR: Invalid logical value for CLEANBOUNDARYVELOCITIES.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("INTERPOLATEBOUNDARYVELOCITY")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        interpolateBoundaryVelocities = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        interpolateBoundaryVelocities = false;
      }else{
        throw mriException("ERROR: Invalid logical value for CLEANBOUNDARYVELOCITIES.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("USECONSTANTPATTERNS")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        useConstantPatterns = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        useConstantPatterns = false;
      }else{
        throw mriException("ERROR: Invalid logical value for useConstantPatterns.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("SEQUENCEFILENAME")){
      try{
        haveSequence = true;
        sequenceFileName = tokenizedString.at(1);
      }catch(...){
        throw mriException("ERROR: Invalid Sequence File Name.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("EVALVORTEXCRITERIA")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        evalPopVortexCriteria = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        evalPopVortexCriteria = false;
      }else{
        throw mriException("ERROR: Invalid logical value for evalPopVortexCriteria.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("EVALSMPVORTEXCRITERIA")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        evalSMPVortexCriterion = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        evalSMPVortexCriterion = false;
      }else{
        throw mriException("ERROR: Invalid logical value for evalSMPVortexCriterion.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("EVALPRESSURE")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        evalPressure = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        evalPressure = false;
      }else{
        throw mriException("ERROR: Invalid logical value for evalPressure.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("EXPORTTOPOISSON")){
      try{
        poissonFileName = tokenizedString.at(1);
      }catch(...){
        throw mriException("ERROR: Invalid Export To Poisson Command Line.\n");
      }    
      // Create Operation For Poisson Export
      mriOperation* op = new mriOpExportForPoissonSolver(poissonFileName,
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
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("EXPORTTODISTANCE")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        exportToDistance = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        exportToDistance = false;
      }else{
        throw mriException("ERROR: Invalid logical value for exportToDistance.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("DISTANCEFILE")){
      distanceFileName = tokenizedString.at(1);
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("ADDNOISE")){
      applyNoise = true;
      try{
        noiseIntensity = atof(tokenizedString.at(1).c_str());        
      }catch(...){
        throw mriException("ERROR: Invalid Noise Intensity.\n");
      }
      if(tokenizedString.size() > 2){
	      try{
	        noiseSeed = atof(tokenizedString.at(2).c_str());
        }catch(...){
	        throw mriException("ERROR: Invalid Noise Seed.\n");
        }
      }
      // Create Operation
      mriOperation* op = new mriOpApplyNoise(noiseIntensity,noiseSeed);
      // Add to the operation list
      operationList.push_back(op);
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("TEMPLATETYPE")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("ZEROVELOCITY")){
        templateType = kZeroVelocity;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("CONSTANT")){
        templateType = kConstantFlow;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("POISEILLE")){
        templateType = kPoiseilleFlow;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("STAGNATION")){
        templateType = kStagnationFlow;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("CYLINDRICALVOLTEX")){
        templateType = kCylindricalVortex;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("SPHERICALVORTEX")){
        templateType = kSphericalVortex;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("TOROIDALVORTEX")){
          templateType = kToroidalVortex;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRANSIENT")){
          templateType = kTransientFlow;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("CONSTANTWITHSTEP")){
          templateType = kConstantFlowWithStep;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("TAYLORVORTEX")){
          templateType = kTaylorVortex;
      }else{
        throw mriException("ERROR: Invalid logical value for exportToPoisson.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("TEMPLATEPARAMS")){
      try{
        for(int loopA=1;loopA<tokenizedString.size();loopA++){
          templateParams.push_back(atof(tokenizedString[loopA].c_str()));
        }
      }catch(...){
        throw mriException("ERROR: Invalid template parameters.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("USETURBVISCOSITYFROMFILE")){
      if(boost::to_upper_copy(tokenizedString.at(1)) == string("TRUE")){
        readMuTFromFile = true;
      }else if(boost::to_upper_copy(tokenizedString.at(1)) == string("FALSE")){
        readMuTFromFile = false;
      }else{
        throw mriException("ERROR: Invalid logical value for USETURBVISCOSITYFROMFILE.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("TURBVISCOSITYFILE")){
      muTFile = tokenizedString[1];
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("SMAGORINSKYCONSTANT")){
      try{
        smagorinskyCoeff = atof(tokenizedString.at(1).c_str());
      }catch(...){
        throw mriException("ERROR: Invalid Smagorinsky Constant.\n");
      }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("SCALEVELOCITY")){
        try{
          scaleVelocities = true;
          scaleVelocityFactor = atof(tokenizedString.at(1).c_str());
        }catch(...){
          throw mriException("ERROR: Invalid velocity scaling parameters.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("SCALEPOSITION")){
        try{
          scalePositions = true;
          scalePositionFactor = atof(tokenizedString.at(1).c_str());
        }catch(...){
          throw mriException("ERROR: Invalid position scaling parameters.\n");
        }
    }else if(boost::to_upper_copy(tokenizedString.at(0)) == string("PRESSUREGRADIENTCOMPONENTS")){
        try{
          // ACCELERATION TERM
          if(tokenizedString.at(1).compare("Y") == 0){
            PPE_IncludeAccelerationTerm = true;
          }else if(tokenizedString.at(1).compare("N") == 0){
            PPE_IncludeAccelerationTerm = false;
          }else{
            throw mriException("ERROR: Invalid token for acceleration term.\n");
          }
          // ADVECTION TERM
          if(tokenizedString.at(2).compare("Y") == 0){
            PPE_IncludeAdvectionTerm = true;
          }else if(tokenizedString.at(2).compare("N") == 0){
            PPE_IncludeAdvectionTerm = false;
          }else{
            throw mriException("ERROR: Invalid token for advection term.\n");
          }
          // DIFFUSION TERM
          if(tokenizedString.at(3).compare("Y") == 0){
            PPE_IncludeDiffusionTerm = true;
          }else if(tokenizedString.at(3).compare("N") == 0){
            PPE_IncludeDiffusionTerm = false;
          }else{
            throw mriException("ERROR: Invalid token for diffusion term.\n");
          }
          // REYNOLDS TERM
          if(tokenizedString.at(4).compare("Y") == 0){
            PPE_IncludeReynoldsTerm = true;
          }else if(tokenizedString.at(4).compare("N") == 0){
            PPE_IncludeReynoldsTerm = false;
          }else{
            throw mriException("ERROR: Invalid token for Reynolds term.\n");
          }
        }catch(...){
          throw mriException("ERROR: Invalid pressure gradient component inclusion command.\n");
        }
    }else if((tokenizedString.at(0).empty())||(tokenizedString[0].at(0) == '#')){
      // Comment: Do Nothing
    }else{
      string errorMsg("ERROR: Invalid Token in input File: " + tokenizedString[0] + "\n");
      throw mriException(errorMsg.c_str());
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
void mriOptions::DistributeProgramOptions(mriCommunicator* comm){
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
  mriUtils::checkMpiError(mpiError);
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
  mriUtils::checkMpiError(mpiError);
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
  mriUtils::checkMpiError(mpiError);
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

