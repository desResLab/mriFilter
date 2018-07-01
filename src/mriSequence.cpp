# include "mriSequence.h"

// FORM BIN LIMITS
void mriSequence::formDifferenceBinLimits(int otherScan, int refScan, 
                                          int pdfQuantity, double& currInterval,
                                          const mriDoubleVec& limitBox, 
                                          int numberOfBins, 
                                          mriDoubleVec& binMin, 
                                          mriDoubleVec& binMax, 
                                          mriDoubleVec& binCenter){
  // Get The Scans out of the sequence
  mriScan* scanOther = getScan(otherScan);
  mriScan* scanRef = getScan(refScan);
  // Initialize Limits
  double minRange = std::numeric_limits<double>::max();
  double maxRange = -std::numeric_limits<double>::max();
  double* cellCoord = NULL;
  double otherQuantity = 0.0;
  double refQuantity = 0.0;
  double currValue = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Get quantity
    otherQuantity = scanOther->cells[loopA].getQuantity(pdfQuantity);
    refQuantity = scanRef->cells[loopA].getQuantity(pdfQuantity);
    // Get Value
    currValue = (otherQuantity - refQuantity);
    // Check If Within the Bin 
    if (mriUtils::isPointInsideBox(topology->cellLocations[loopA][0],
                                   topology->cellLocations[loopA][1],
                                   topology->cellLocations[loopA][2],limitBox)){
      // Assign Values
      if(currValue>maxRange) maxRange = currValue;
      if(currValue<minRange) minRange = currValue;
    }
  }
  // If minRange and maxRange are the same than add something
  if (fabs(maxRange - minRange)<kMathZero){
    minRange = minRange - 1.0;
    maxRange = maxRange + 1.0;
  }
  // Fill the bin arrays
  currInterval = ((maxRange - minRange)/(double)numberOfBins);
  double currPtr = minRange;
  for(int loopA=0;loopA<numberOfBins;loopA++){
    binMin[loopA] = currPtr;
    binMax[loopA] = currPtr + currInterval;
    binCenter[loopA] = 0.5*(binMin[loopA] + binMax[loopA]);
    // Update
    currPtr += currInterval;
  }
}

// Constructor
mriSequence::mriSequence(bool cyclic){
  // Set If Cyclic
  isCyclic = cyclic;
  // Initialize Topology
  topology = new mriTopology();
}

// Copy Constructor
mriSequence::mriSequence(mriSequence* copySequence){
  // Copy Cyclic Property
  isCyclic = copySequence->isCyclic;
  // Copy cells totals
  for(int loopA=0;loopA<3;loopA++){
    topology->domainSizeMin[loopA] = copySequence->topology->domainSizeMin[loopA];
    topology->domainSizeMax[loopA] = copySequence->topology->domainSizeMax[loopA];
  }
  // Total number of Cells
  topology->totalCells = copySequence->topology->totalCells;
  // Allocate CellLocations
  topology->cellLocations.resize(topology->totalCells);
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    topology->cellLocations[loopA].resize(3);
  }
  // Initialize
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Copy position
    topology->cellLocations[loopA][0] = copySequence->topology->cellLocations[loopA][0];
    topology->cellLocations[loopA][1] = copySequence->topology->cellLocations[loopA][1];
    topology->cellLocations[loopA][2] = copySequence->topology->cellLocations[loopA][2];
  }  
  // Copy Cell Totals
  topology->cellLengths.resize(3);
  for(int loopA=0;loopA<3;loopA++){
    topology->cellTotals[loopA] = copySequence->topology->cellTotals[loopA];    
    topology->cellLengths[loopA].resize(copySequence->topology->cellLengths[loopA].size());
    // FILL THE LENGTHS
    for(size_t loopB=0;loopB<copySequence->topology->cellLengths[loopA].size();loopB++){
      topology->cellLengths[loopA][loopB] = copySequence->topology->cellLengths[loopA][loopB];
    }
  }

  // Fill with Zero Scans
  for(int loopB=0;loopB<sequence.size();loopB++){
    mriScan* newScan = new mriScan(*copySequence->getScan(loopB));
    sequence.push_back(newScan);
  }
}

// Destructor
mriSequence::~mriSequence(){}

// Print the File List Log
void mriSequence::printSequenceFiles(std::string outFIleName){
  // Open Output File
	FILE* outFile;
	outFile = fopen(outFIleName.c_str(),"w");
  // Write List
  fprintf(outFile,"List of Files in Sequence\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    fprintf(outFile,"File %d: %s\n",loopA+1,fileNames[loopA].c_str());
  }
  // Close Output file
	fclose(outFile);				
}

// Add a Scan to the Sequence
void mriSequence::addScan(mriScan* scan){
  // Add the Scan
  sequence.push_back(scan);
  // Assign the Topology Pointer
  getScan(sequence.size()-1)->topology = this->topology;
}

// Get a Scan Pointer From the Sequence
mriScan* mriSequence::getScan(int scanNumber){
  return sequence[scanNumber];
}

// EXPORT SEQUENCE TO TECPLOT FILE
void mriSequence::exportToTECPLOT(std::string outfileName){
  writeSchMessage(std::string("\n"));
  writeSchMessage(std::string("EXPORTING -------------------------------------\n"));
  for(int loopA=0;loopA<sequence.size();loopA++){
      sequence[loopA]->exportToTECPLOT(outfileName,(loopA == 0));
  }
}

// EXCTRACT SINGKLE POINT CURVE IN TIME
void mriSequence::extractSinglePointTimeCurve(int cellNumber, int exportQty, std::string fileName){
  // Open Output File
	FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");
  // Var
  double currTime = 0.0;
  double currValue = 0.0;
  double currRad = 0.0;
  // Eval Centre Point
  double centrePoint[3] = {0.0};
  double localCoord[3] = {0.0};
  centrePoint[0] = 0.5 * (topology->domainSizeMax[0] + topology->domainSizeMin[0]);
  centrePoint[1] = 0.5 * (topology->domainSizeMax[1] + topology->domainSizeMin[1]);
  centrePoint[2] = 0.5 * (topology->domainSizeMax[2] + topology->domainSizeMin[2]);
  // Write List
  fprintf(outFile,"List of Files in Sequence\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    // Get Time
    currTime = sequence[loopA]->scanTime;
    // Get Local Coordinates
    localCoord[0] = topology->cellLocations[cellNumber][0] - centrePoint[0];
    localCoord[1] = topology->cellLocations[cellNumber][1] - centrePoint[1];
    localCoord[2] = topology->cellLocations[cellNumber][2] - centrePoint[2];
    // Check Quantity
    switch(exportQty){
      case kQtyConcentration:
        currValue = sequence[loopA]->cells[cellNumber].concentration;
        break;
      case kQtyVelocityX:
        currValue = sequence[loopA]->cells[cellNumber].velocity[0];
        break;      
      case kQtyVelocityY:
        currValue = sequence[loopA]->cells[cellNumber].velocity[1];
        break;      
      case kQtyVelocityZ:
        currValue = sequence[loopA]->cells[cellNumber].velocity[2];
        break;      
    }
    currRad = sqrt(localCoord[1]*localCoord[1]+
                   localCoord[2]*localCoord[2]);
    fprintf(outFile,"%e %e %e\n",currTime,currValue,currRad);
  }
  // Close Output file
	fclose(outFile);				  
}

// Export to Poisson Solver
void mriSequence::exportForPoisson(string inputFileName,double density,double viscosity,mriThresholdCriteria* threshold,
                                   bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                                   bool readMuTFromFile, string muTFile, double smagorinskyCoeff){
  string name;
  mriDoubleMat timeDeriv;
  // mriDoubleMat reynoldsDeriv;
  for(int loopA=0;loopA<sequence.size();loopA++){
    name = inputFileName + "_" + mriUtils::floatToStr(loopA);
    if(PPE_IncludeAccelerationTerm){
      printf("Computing Time Derivatives for Scan %d...",loopA);
      evalScanTimeDerivs(loopA,timeDeriv);
      printf("Done.\n");
    }
    sequence[loopA]->exportForPoisson(name,density,viscosity,threshold,timeDeriv,
                                      PPE_IncludeAccelerationTerm,PPE_IncludeAdvectionTerm,PPE_IncludeDiffusionTerm,PPE_IncludeReynoldsTerm,
                                      readMuTFromFile,muTFile,smagorinskyCoeff);
  }
}

// EXPORT TO WALL DISTANCE SOLVER
void mriSequence::exportForDistancing(string inputFileName, mriThresholdCriteria* threshold){
  string name;
  for(int loopA=0;loopA<sequence.size();loopA++){
    name = inputFileName + "_" + mriUtils::floatToStr(loopA);
    sequence[loopA]->exportForDistancing(name,threshold);
  }
}

// GET SEQUENCE OUTPUT FILE NAME
string getSequenceOutputFileName(string outfileName, int loopA){
  // Get Extension
  string ext = outfileName.substr(outfileName.find_last_of('.'),outfileName.size()-1);
  // Get File Name
  string res = outfileName.substr(0,outfileName.find_last_of('.')) + "_" + to_string(loopA) + ext;
  // Return
  return res;
}

// EXPORT TO SEQUENCE OF VTK FILES
void mriSequence::exportToVTK(string outfileName,mriThresholdCriteria* thresholdCriteria){
  // Export All Data
  for(int loopA=0;loopA<sequence.size();loopA++){    
    string outName = getSequenceOutputFileName(outfileName,loopA);
    sequence[loopA]->exportToVTK(outName,thresholdCriteria);
  }
}

// PHYSICS FILTERING FOR ALL SCANS
void mriSequence::applySMPFilter(mriCommunicator* comm, bool isBC, 
                                 mriThresholdCriteria* thresholdCriteria,
                                 double itTol,
                                 int maxIt,
                                 bool useConstantPatterns){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    // Perform Filter
    sequence[loopA]->applySMPFilter(comm,isBC,thresholdCriteria,itTol,maxIt,useConstantPatterns);
    // Update Velocities
    sequence[loopA]->updateVelocities();
  }
}

// SAVE INITIAL VELOCITIES
void mriSequence::saveVelocity(){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    // Perform Filter
    sequence[loopA]->saveVelocity();
  }
}


// APPLY THRESHOLDING TO ALL SCANS
void mriSequence::applyThresholding(mriThresholdCriteria* thresholdCriteria){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    sequence[loopA]->applyThresholding(thresholdCriteria);
  }
}

// EVAL VORTEX CRITERIA
void mriSequence::evalVortexCriteria(mriThresholdCriteria* thresholdCriteria){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    sequence[loopA]->evalVortexCriteria(thresholdCriteria);
  }
}

// EVAL VORTICITY
void mriSequence::evalVorticity(mriThresholdCriteria* thresholdCriteria){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    sequence[loopA]->evalVorticity(thresholdCriteria);
  }
}

// EVAL ENSTROPHY
void mriSequence::evalEnstrophy(mriThresholdCriteria* thresholdCriteria){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    sequence[loopA]->evalEnstrophy(thresholdCriteria);
  }
}

// EVAL SMP VORTEX CRITERION
void mriSequence::evalSMPVortexCriteria(){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    sequence[loopA]->evalSMPVortexCriteria(sequence[loopA]->expansion);
  }
}

// EVAL EXPANSION FILE
void mriSequence::writeExpansionFile(string fileName){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<sequence.size();loopA++){
    sequence[loopA]->writeExpansionFile(fileName);
  }
}

// Scale velocities for all Scans
void mriSequence::scaleVelocities(double factor){
  writeSchMessage(std::string("Scaling Velocities..."));
  for(int loopA=0;loopA<sequence.size();loopA++){
    sequence[loopA]->scaleVelocities(factor);
  }
  writeSchMessage(std::string("Done.\n"));
}

// Scale Positions
void mriSequence::scalePositions(const mriDoubleVec& origin, double factor){
  writeSchMessage(std::string("Scaling Positions..."));
  topology->scalePositions(origin,factor);
  writeSchMessage(std::string("Done.\n"));  
}

// Add noise to measurements
void mriSequence::applyNoise(double noiseIntensity, double seed){
  writeSchMessage(std::string("Applying Noise..."));
  for(int loopA=0;loopA<sequence.size();loopA++){
    sequence[loopA]->applyGaussianNoise(noiseIntensity, seed);
  }
  writeSchMessage(std::string("Done.\n"));
}

// =============================
// DISTRIBUTION OF SEQUENCE DATA
// =============================
void mriSequence::distributeSequenceData(mriCommunicator* comm){
  // Create New Sequence
  for(int loopA=0;loopA<this->sequence.size();loopA++){
    sequence[loopA]->distributeScanData(comm);
  }
}

// =============================
// PERFORM BASIC DATA FILTERING
// =============================
void mriSequence::applyMedianFilter(int qtyID,int maxIt,int order,int filterType,mriThresholdCriteria* threshold){
  // Create New Sequence
  for(int loopA=0;loopA<this->sequence.size();loopA++){
    sequence[loopA]->applyMedianFilter(qtyID,maxIt,order,filterType,threshold);
  }
}

// ===========================
// CLEAN COMPONENT ON BOUNDARY
// ===========================
void mriSequence::cleanNormalComponentOnBoundary(mriThresholdCriteria* threshold){
  for(int loopA=0;loopA<this->sequence.size();loopA++){
    sequence[loopA]->cleanNormalComponentOnBoundary(threshold);
  }
}

// ===============================
// INTERPOLATE BOUNDARY VELOCITIES
// ===============================
void mriSequence::interpolateBoundaryVelocities(mriThresholdCriteria* threshold){
  printf("Interpolating Boundary Velocities...\n");
  for(int loopA=0;loopA<this->sequence.size();loopA++){
    sequence[loopA]->interpolateBoundaryVelocities(threshold);
  }
}

// ===============================
// EVALUATE REYNOLDS STRESS TENSOR
// ===============================
void mriSequence::evalReynoldsStresses(mriThresholdCriteria* threshold){
  for(int loopA=0;loopA<sequence.size();loopA++){
    printf("Evaluating Reynolds Stress Tensor for Scan %d...",loopA);
    sequence[loopA]->evalReynoldsStress(threshold);
    printf("Done.");
  }
}

// =============================
// CREATE SEQUENCE MESH TOPOLOGY
// =============================
void mriSequence::createTopology(){
  // Take Time
  float cellConn_BeginTime,cellConn_TotalTime;
  float faceConn_BeginTime,faceConn_TotalTime;
  float faceArea_BeginTime,faceArea_TotalTime;
  float edgeConn_BeginTime,edgeConn_TotalTime;
  float auxNodes_BeginTime,auxNodes_TotalTime;

  // Build Cell Connections
  writeSchMessage(std::string("Build Cell Connection...\n"));
  cellConn_BeginTime = clock();
  topology->buildCellConnections();  
  cellConn_TotalTime = float( clock () - cellConn_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",cellConn_TotalTime);

  writeSchMessage(std::string("Build Aux Node Coords...\n"));
  auxNodes_BeginTime = clock();
  topology->buildAuxNodesCoords();
  auxNodes_TotalTime = float( clock () - auxNodes_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",auxNodes_TotalTime);

  // Build Face Connections
  writeSchMessage(std::string("Build Face Connection...\n"));
  faceConn_BeginTime = clock();
  topology->buildFaceConnections();
  topology->buildFaceCells();
  faceConn_TotalTime = float( clock () - faceConn_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",faceConn_TotalTime);

  // Build Face Area and Face Normal Vector
  writeSchMessage(std::string("Build Areas and Normals...\n"));
  faceArea_BeginTime = clock();
  topology->buildFaceAreasAndNormals();
  faceArea_TotalTime = float( clock () - faceArea_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",faceArea_TotalTime);

  // Build Edge Connections
  writeSchMessage(std::string("Build Edge Connections...\n"));
  edgeConn_BeginTime = clock();
  topology->buildEdgeConnections();
  edgeConn_TotalTime = float( clock () - edgeConn_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",edgeConn_TotalTime);

  // Completed
  writeSchMessage(std::string("Topology Creation Completed.\n"));
}

// =============
// CROP SEQUENCE
// =============
void mriSequence::crop(const mriDoubleVec& limitBox){
  mriBoolVec indexesToCrop;

  // Crop Topology
  topology->crop(limitBox,indexesToCrop);

  // Crop Scans
  for(int loopA=0;loopA<sequence.size();loopA++){
    printf("Cropping Scan %d...",loopA);
    sequence[loopA]->crop(limitBox,indexesToCrop);
    printf("Done.\n");
  }
}

// Eval The Difference PDF of Scans
void mriSequence::evalScanDifferencePDF(int otherScan, int refScan, const int pdfQuantity, int numberOfBins, bool useBox, mriDoubleVec& limitBox, mriDoubleVec& binCenters, mriDoubleVec& binArray){
  // Get The Scans out of the sequence
  mriScan* scanOther = getScan(otherScan);
  mriScan* scanRef = getScan(refScan);
  // Allocate Quantities
  mriDoubleVec binMin(numberOfBins);
  mriDoubleVec binMax(numberOfBins);
  // Form Bin 
  double currInterval = 0.0;
  formDifferenceBinLimits(otherScan,refScan,pdfQuantity,currInterval,limitBox,numberOfBins,binMin,binMax,binCenters);
  // Intialize bins
  for(int loopA=0;loopA<numberOfBins;loopA++){
    binArray[loopA] = 0.0;
  }
  // Loop through the points
  double sourceValue = 0.0;
  double refValue = 0.0;
  double* cellCoord = NULL;
  double otherQuantity = 0.0;
  double refQuantity = 0.0;
  double currValue = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // Get quantity
    otherQuantity = scanOther->cells[loopA].getQuantity(pdfQuantity);
    refQuantity = scanRef->cells[loopA].getQuantity(pdfQuantity);
    // Get Value
    currValue = (otherQuantity - refQuantity);
    // Assign Value to Bin
    if (mriUtils::isPointInsideBox(topology->cellLocations[loopA][0],
                                   topology->cellLocations[loopA][1],
                                   topology->cellLocations[loopA][2],
                                   limitBox)){
      // COMPLETE
      mriUtils::assignToBin(currValue,numberOfBins,binMin,binMax,binArray);
    }
  }
  // Normalize
  mriUtils::normalizeBinArray(binArray,currInterval);
}

// READ SCAN FROM EXPANSION FILE
void mriSequence::readFromExpansionFiles(const mriStringVec& fileNames, const mriDoubleVec& Times, bool applyThreshold, int thresholdType,double thresholdRatio){

  // Allocate Variables
  mriIntVec tot(3);
  mriDoubleVec lengthX;
  mriDoubleVec lengthY;
  mriDoubleVec lengthZ;
  mriDoubleVec minlimits(3);
  mriDoubleVec maxlimits(3);
  mriExpansion* exp = NULL;

  mriScan* scan = NULL;

  // Loop Through the files
  for(int loopA=0;loopA<fileNames.size();loopA++){

    // Read current filename
    readExpansionFile(fileNames[loopA],tot,lengthX,lengthY,lengthZ,minlimits,maxlimits,exp);

    // Build a New Topology from this file
    mriTopology* topo = new mriTopology(tot,lengthX,lengthY,lengthZ,minlimits,maxlimits);

    if(loopA == 0){
      // Assign as the full sequence topology
      topology = topo;
      // Add Scan
      scan = new mriScan(Times[loopA]);
      scan->rebuildFromExpansion(exp, true);
      addScan(scan);
    }else{
      // Check Compatibility
      if(topology->isCompatibleTopology(topo)){

        // Add Scan From File
        scan = new mriScan(Times[loopA]);
        scan->rebuildFromExpansion(exp, true);
        addScan(scan);

      }else{
        // Skip Scan due to incompatible topology
        printf("WARNING: Skipping Scan, topology is not compatible.\n");
      }
    }

    // Delete expansion
    delete exp;
    exp = NULL;
  }
}

// Make Difference of Scans
// Put the Result Of the Operation in firstScan
void mriSequence::makeScanDifference(int firstScanID, int secondScanID){
  // Get Scans
  if((firstScanID<0)||(firstScanID>sequence.size())){
    throw mriException("Cannot Make Difference. Invalid ID for First Scan.");
  }
  if((secondScanID<0)||(secondScanID>sequence.size())){
    throw mriException("Cannot Make Difference. Invalid ID for Second Scan.");
  }
  mriScan* firstScan = getScan(firstScanID);
  mriScan* secondScan = getScan(secondScanID);
  // set the Tolerance Value
  double DistTol = 1.0e-4;
  // If they Belong to the same sequence they must be compatible
  writeSchMessage("Computing Scan Difference...");  
  // SubTract Velocity Data
  double diffX,diffY,diffZ;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    for(int loopB=0;loopB<3;loopB++){
      firstScan->cells[loopA].velocity[loopB] = firstScan->cells[loopA].velocity[loopB]-secondScan->cells[loopA].velocity[loopB];
    }
  }
}

// Make Average
// Put the Result Of the Operation in firstScan
void mriSequence::makeScanAverage(int numberOfMeasures, int firstScanID, int secondScanID){
  // Get Scans
  if((firstScanID<0)||(firstScanID>sequence.size())){
    throw mriException("Cannot Make Difference. Invalid ID for First Scan.");
  }
  if((secondScanID<0)||(secondScanID>sequence.size())){
    throw mriException("Cannot Make Difference. Invalid ID for Second Scan.");
  }
  mriScan* firstScan = sequence[firstScanID];
  mriScan* secondScan = sequence[secondScanID];
  // Get Distance Tolerance
  double DistTol = 1.0e-4;
  // If they Belong to the same sequence they must be compatible
  writeSchMessage("Computing Scan Average...");
  // SubTract Velocity Data
  double diffX,diffY,diffZ;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    for(int loopB=0;loopB<3;loopB++){
      firstScan->cells[loopA].velocity[loopB] =
      firstScan->cells[loopA].velocity[loopB]*((double)(numberOfMeasures-1)/(double)numberOfMeasures)+
      secondScan->cells[loopA].velocity[loopB]*(1.0/(double)numberOfMeasures);
    }
  }
}

// ============================================
// WRITE SEQUENCE TOPOLOGY STATISTICS TO STDOUT
// ============================================
string mriSequence::writeStatistics(){
  string myresult = "\n";
  myresult += "--------------------------------\n";
  myresult += "SEQUENCE STATISTICS\n";
  myresult += "--------------------------------\n";
  myresult += "Total Number Of Cells: "+mriUtils::intToStr(topology->totalCells)+"\n";
  myresult += "--------------------------------\n";
  myresult += "Total Number Of Coordinate Cells\n";
  myresult += "X Direction: "+mriUtils::intToStr(topology->cellTotals[0])+"\n";
  myresult += "Y Direction: "+mriUtils::intToStr(topology->cellTotals[1])+"\n";
  myresult += "Z Direction: "+mriUtils::intToStr(topology->cellTotals[2])+"\n";
  myresult += "Cells Lengths\n";
  myresult += "X Direction - MIN: "+mriUtils::floatToStr(*min_element(topology->cellLengths[0].begin(),topology->cellLengths[0].end()))+"\n";
  myresult += "X Direction - MAX: "+mriUtils::floatToStr(*max_element(topology->cellLengths[0].begin(),topology->cellLengths[0].end()))+"\n";
  myresult += "Y Direction - MIN: "+mriUtils::floatToStr(*min_element(topology->cellLengths[1].begin(),topology->cellLengths[1].end()))+"\n";
  myresult += "Y Direction - MAX: "+mriUtils::floatToStr(*max_element(topology->cellLengths[1].begin(),topology->cellLengths[1].end()))+"\n";
  myresult += "Z Direction - MIN: "+mriUtils::floatToStr(*min_element(topology->cellLengths[2].begin(),topology->cellLengths[2].end()))+"\n";
  myresult += "Z Direction - MAX: "+mriUtils::floatToStr(*max_element(topology->cellLengths[2].begin(),topology->cellLengths[2].end()))+"\n";
  myresult += "--------------------------------\n";
  myresult += "Domain Size\n";
  myresult += "Minimum X: "+mriUtils::floatToStr(topology->domainSizeMin[0])+"\n";
  myresult += "Maximum X: "+mriUtils::floatToStr(topology->domainSizeMax[0])+"\n";
  myresult += "Minimum Y: "+mriUtils::floatToStr(topology->domainSizeMin[1])+"\n";
  myresult += "Maximum Y: "+mriUtils::floatToStr(topology->domainSizeMax[1])+"\n";
  myresult += "Minimum Z: "+mriUtils::floatToStr(topology->domainSizeMin[2])+"\n";
  myresult += "Maximum Z: "+mriUtils::floatToStr(topology->domainSizeMax[2])+"\n";
  myresult += "--------------------------------\n";
  myresult += "Number of Scans: "+mriUtils::intToStr(sequence.size())+"\n";
  myresult += "--------------------------------\n";
  myresult += "\n";
  for(size_t loopA=0;loopA<sequence.size();loopA++){
    myresult += "--------------------------------\n";
    myresult += "Scan number: "+to_string(loopA)+"\n";
    myresult += "Scan time: "+to_string(sequence[loopA]->scanTime)+"\n";
    myresult += "Maximum Velocity Module: "+to_string(sequence[loopA]->maxVelModule)+"\n";
    myresult += "--------------------------------\n";
    myresult += "\n";
  }
  // Return String
  return myresult;
}

// ====================
// CREATE TEMPLATE CASE
// ====================
void mriSequence::createSampleCase(int sampleType, const mriDoubleVec& params){

  // Loop over the file names
  mriTopology* topo;
  mriScan* scan;

  // Create New Topology from Sample
  topo = new mriTopology();
  topo->createFromTemplate(sampleType,params);  

  // If topology does not exists then assign 
  if(sequence.size() == 0){
    // Assign Current Topology
    topology = topo;
    // Create and Assign Scan
    scan = new mriScan(0.0);
    scan->topology = topo;
    scan->createFromTemplate(sampleType,params);
    // Add to sequence
    addScan(scan);
  }else{
    // Check Compatibility 
    if(topology->isCompatibleTopology(topo)){
      // Create and Assign Scan
      scan = new mriScan(0.0);
      scan->topology = topo;
      scan->createFromTemplate(sampleType,params);
      // Add to sequence
      addScan(scan);
    }
  }

  string CurrentStats = writeStatistics();
  writeSchMessage(CurrentStats);

}

// ======================
// INITIALIZE VTK OPTIONS
// ======================
void resetVTKOptions(vtkStructuredPointsOptionRecord& vtkOptions){
  vtkOptions.dataBlockStart.clear();
  vtkOptions.dataBlockType.clear();
  vtkOptions.dataBlockRead.clear();
}

// ======================
// INITIALIZE PLT OPTIONS
// ======================
void resetPLTOptions(pltOptionRecord& pltOptions){
  pltOptions.i = 0;
  pltOptions.j = 0;
  pltOptions.k = 0;
  pltOptions.N = 0;
  pltOptions.E = 0;
  pltOptions.type = pltUNIFORM;
}

// ==========================
// READ VTK STRUCTURED POINTS
// ==========================
void mriSequence::readFromASCIISequence(int asciiInputType,const mriStringVec& asciiFileNames, const mriDoubleVec& times){

  // Loop over the file names
  mriTopology* topo;
  mriScan* scan;
  vtkStructuredPointsOptionRecord vtkOptions;
  pltOptionRecord pltOptions;
  for(int loopA=0;loopA<asciiFileNames.size();loopA++){

    // Read Topology
    topo = new mriTopology();
    if(asciiInputType == kInputVTK){
      resetVTKOptions(vtkOptions);
      topo->readFromVTK_ASCII(asciiFileNames[loopA],vtkOptions);  
    }else if(asciiInputType == kInputPLT){
      resetPLTOptions(pltOptions);
      topo->readFromPLT_ASCII(asciiFileNames[loopA],pltOptions);  
    }
    
    if(loopA == 0){
      
      // Assign as the full sequence topology
      topology = topo;
      
      // Add Scan
      scan = new mriScan(times[loopA]);
      if(asciiInputType == kInputVTK){
        scan->readFromVTK_ASCII(asciiFileNames[loopA],vtkOptions);
      }else if(asciiInputType == kInputPLT){
        scan->readFromPLT_ASCII(asciiFileNames[loopA],pltOptions);
      }    
      
      // Add scan to sequence
      addScan(scan);
    }else{
      // Check Compatibility
      if(topology->isCompatibleTopology(topo)){

        // Add Scan From File
        scan = new mriScan(times[loopA]);
        if(asciiInputType == kInputVTK){
          scan->readFromVTK_ASCII(asciiFileNames[loopA],vtkOptions);
        }else if(asciiInputType == kInputPLT){
          scan->readFromPLT_ASCII(asciiFileNames[loopA],pltOptions);
        }    
        addScan(scan);

      }else{
        
        // Skip Scan due to incompatible topology
        printf("WARNING: Skipping Scan, topology is not compatible.\n");

      }
      // Delete current topology
      delete topo;
      topo = NULL;
    }
  }

  // WRITE STATISTICS
  string CurrentStats = writeStatistics();
  writeSchMessage(CurrentStats);

}
