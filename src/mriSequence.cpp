# include "mriSequence.h"

// FORM BIN LIMITS
void MRISequence::formDifferenceBinLimits(int otherScan, int refScan, 
                                          int pdfQuantity, double& currInterval,
                                          const MRIDoubleVec& limitBox, 
                                          int numberOfBins, 
                                          MRIDoubleVec& binMin, 
                                          MRIDoubleVec& binMax, 
                                          MRIDoubleVec& binCenter){
  // Get The Scans out of the sequence
  MRIScan* scanOther = getScan(otherScan);
  MRIScan* scanRef = getScan(refScan);
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
    if (MRIUtils::isPointInsideBox(topology->cellLocations[loopA][0],
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
MRISequence::MRISequence(bool cyclic){
  // Initizalize
  totalScans = 0;
  // Set If Cyclic
  isCyclic = cyclic;
  // Initialize Topology
  topology = new MRITopology();
}

// Copy Constructor
MRISequence::MRISequence(MRISequence* copySequence){
  // Copy Scans
  totalScans = copySequence->totalScans;
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
  for(int loopB=0;loopB<totalScans;loopB++){
    MRIScan* newScan = new MRIScan(*copySequence->getScan(loopB));
    sequence.push_back(newScan);
  }
}

// Destructor
MRISequence::~MRISequence(){}

// Print the File List Log
void MRISequence::printSequenceFiles(std::string outFIleName){
  // Open Output File
	FILE* outFile;
	outFile = fopen(outFIleName.c_str(),"w");
  // Write List
  fprintf(outFile,"List of Files in Sequence\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    fprintf(outFile,"File %d: %s\n",loopA+1,fileNames[loopA].c_str());
  }
  // Close Output file
	fclose(outFile);				
}

// Add a Scan to the Sequence
void MRISequence::addScan(MRIScan* scan){
  sequence.push_back(scan);
  totalScans++;
}

// Get a Scan Pointer From the Sequence
MRIScan* MRISequence::getScan(int scanNumber){
  return sequence[scanNumber];
}

// EXPORT SEQUENCE TO TECPLOT FILE
void MRISequence::exportToTECPLOT(std::string outfileName){
  writeSchMessage(std::string("\n"));
  writeSchMessage(std::string("EXPORTING -------------------------------------\n"));
  for(int loopA=0;loopA<totalScans;loopA++){
      sequence[loopA]->exportToTECPLOT(outfileName,(loopA == 0));
  }
}

// EXCTRACT SINGKLE POINT CURVE IN TIME
void MRISequence::extractSinglePointTimeCurve(int cellNumber, int exportQty, std::string fileName){
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
  for(int loopA=0;loopA<totalScans;loopA++){
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

// Compute Relative Pressure
void MRISequence::computeRelativePressure(bool doPressureSmoothing){
  writeSchMessage(std::string("\n"));
  writeSchMessage(std::string("REL PRESSURE COMPUTATION --------------------------\n"));
  // Loop Through Scans
  for(int loopA=0;loopA<totalScans;loopA++){
    // Write Message
    writeSchMessage(std::string("Computing Relative Pressure - Step "+MRIUtils::intToStr(loopA+1)+"/"+MRIUtils::intToStr(totalScans)+"..."));
    // Get The Scan Back
    MRIScan* resultScan = getScan(loopA);
    int startingCell = resultScan->evalCentralCell();
    resultScan->evalRelativePressure(startingCell,0.0);
    if (doPressureSmoothing){
      resultScan->performPressureIterations();
    }
    // Done
    writeSchMessage(std::string("Relative Pressure Computed.\n"));
  }
}

// Export to Poisson Solver
void MRISequence::exportForPoisson(string inputFileName,double density,double viscosity,MRIThresholdCriteria* threshold,
                                   bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                                   bool readMuTFromFile, string muTFile, double smagorinskyCoeff){
  string name;
  MRIDoubleMat timeDeriv;
  // MRIDoubleMat reynoldsDeriv;
  for(int loopA=0;loopA<totalScans;loopA++){
    name = inputFileName + "_" + MRIUtils::floatToStr(loopA);
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
void MRISequence::exportForDistancing(string inputFileName, MRIThresholdCriteria* threshold){
  string name;
  for(int loopA=0;loopA<totalScans;loopA++){
    name = inputFileName + "_" + MRIUtils::floatToStr(loopA);
    sequence[loopA]->exportForDistancing(name,threshold);
  }
}

// EXPORT TO SEQUENCE OF VTK FILES
void MRISequence::exportToVTK(std::string outfileName,MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->exportToVTK(outfileName,thresholdCriteria);
  }
}

// PHYSICS FILTERING FOR ALL SCANS
void MRISequence::applySMPFilter(MRICommunicator* comm, bool isBC, 
                                 MRIThresholdCriteria* thresholdCriteria,
                                 double itTol,
                                 int maxIt,
                                 bool useConstantPatterns){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    // Perform Filter
    sequence[loopA]->applySMPFilter(comm,isBC,thresholdCriteria,itTol,maxIt,useConstantPatterns);
    // Update Velocities
    sequence[loopA]->updateVelocities();
  }
}

// SAVE INITIAL VELOCITIES
void MRISequence::saveVelocity(){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    // Perform Filter
    sequence[loopA]->saveVelocity();
  }
}


// APPLY THRESHOLDING TO ALL SCANS
void MRISequence::applyThresholding(MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->applyThresholding(thresholdCriteria);
  }
}

// EVAL VORTEX CRITERIA
void MRISequence::evalVortexCriteria(MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->evalVortexCriteria(thresholdCriteria);
  }
}

// EVAL VORTICITY
void MRISequence::evalVorticity(MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->evalVorticity(thresholdCriteria);
  }
}

// EVAL ENSTROPHY
void MRISequence::evalEnstrophy(MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->evalEnstrophy(thresholdCriteria);
  }
}

// EVAL SMP VORTEX CRITERION
void MRISequence::evalSMPVortexCriteria(){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->evalSMPVortexCriteria(sequence[loopA]->expansion);
  }
}

// EVAL EXPANSION FILE
void MRISequence::writeExpansionFile(string fileName){
  // Export All Data
  writeSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->writeExpansionFile(fileName);
  }
}

// Scale velocities for all Scans
void MRISequence::scaleVelocities(double factor){
  writeSchMessage(std::string("Scaling Velocities..."));
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->scaleVelocities(factor);
  }
  writeSchMessage(std::string("Done.\n"));
}

// Scale Positions
void MRISequence::scalePositions(const MRIDoubleVec& origin, double factor){
  writeSchMessage(std::string("Scaling Positions..."));
  topology->scalePositions(origin,factor);
  writeSchMessage(std::string("Done.\n"));  
}

// Add noise to measurements
void MRISequence::applyNoise(double noiseIntensity){
  writeSchMessage(std::string("Applying Noise..."));
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->applyGaussianNoise(noiseIntensity);
  }
  writeSchMessage(std::string("Done.\n"));
}

// =============================
// DISTRIBUTION OF SEQUENCE DATA
// =============================
void MRISequence::distributeSequenceData(MRICommunicator* comm){
  // Create New Sequence
  for(int loopA=0;loopA<this->totalScans;loopA++){
    sequence[loopA]->distributeScanData(comm);
  }
}

// =============================
// PERFORM BASIC DATA FILTERING
// =============================
void MRISequence::applyMedianFilter(int qtyID,int maxIt,int order,int filterType,MRIThresholdCriteria* threshold){
  // Create New Sequence
  for(int loopA=0;loopA<this->totalScans;loopA++){
    sequence[loopA]->applyMedianFilter(qtyID,maxIt,order,filterType,threshold);
  }
}

// ===========================
// CLEAN COMPONENT ON BOUNDARY
// ===========================
void MRISequence::cleanNormalComponentOnBoundary(MRIThresholdCriteria* threshold){
  for(int loopA=0;loopA<this->totalScans;loopA++){
    sequence[loopA]->cleanNormalComponentOnBoundary(threshold);
  }
}

// ===============================
// INTERPOLATE BOUNDARY VELOCITIES
// ===============================
void MRISequence::interpolateBoundaryVelocities(MRIThresholdCriteria* threshold){
  printf("Interpolating Boundary Velocities...\n");
  for(int loopA=0;loopA<this->totalScans;loopA++){
    sequence[loopA]->interpolateBoundaryVelocities(threshold);
  }
}

// ===============================
// EVALUATE REYNOLDS STRESS TENSOR
// ===============================
void MRISequence::evalReynoldsStresses(MRIThresholdCriteria* threshold){
  for(int loopA=0;loopA<totalScans;loopA++){
    printf("Evaluating Reynolds Stress Tensor for Scan %d...",loopA);
    sequence[loopA]->evalReynoldsStress(threshold);
    printf("Done.");
  }
}

// =============================
// CREATE SEQUENCE MESH TOPOLOGY
// =============================
void MRISequence::createTopology(){
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
void MRISequence::crop(const MRIDoubleVec& limitBox){
  MRIBoolVec indexesToCrop;

  // Crop Topology
  topology->crop(limitBox,indexesToCrop);

  // Crop Scans
  for(int loopA=0;loopA<totalScans;loopA++){
    printf("Cropping Scan %d...",loopA);
    sequence[loopA]->crop(limitBox,indexesToCrop);
    printf("Done.\n");
  }
}

// Eval The Difference PDF of Scans
void MRISequence::evalScanDifferencePDF(int otherScan, int refScan, const int pdfQuantity, int numberOfBins, bool useBox, MRIDoubleVec& limitBox, MRIDoubleVec& binCenters, MRIDoubleVec& binArray){
  // Get The Scans out of the sequence
  MRIScan* scanOther = getScan(otherScan);
  MRIScan* scanRef = getScan(refScan);
  // Allocate Quantities
  MRIDoubleVec binMin(numberOfBins);
  MRIDoubleVec binMax(numberOfBins);
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
    if (MRIUtils::isPointInsideBox(topology->cellLocations[loopA][0],
                                   topology->cellLocations[loopA][1],
                                   topology->cellLocations[loopA][2],
                                   limitBox)){
      // COMPLETE
      MRIUtils::assignToBin(currValue,numberOfBins,binMin,binMax,binArray);
    }
  }
  // Normalize
  MRIUtils::normalizeBinArray(binArray,currInterval);
}

// READ SCAN FROM EXPANSION FILE
void MRISequence::readFromExpansionFiles(const MRIStringVec& fileNames, const MRIDoubleVec& Times, bool applyThreshold, int thresholdType,double thresholdRatio){

  // Allocate Variables
  MRIIntVec tot(3);
  MRIDoubleVec lengthX;
  MRIDoubleVec lengthY;
  MRIDoubleVec lengthZ;
  MRIDoubleVec minlimits(3);
  MRIDoubleVec maxlimits(3);
  MRIExpansion* exp = NULL;

  MRIScan* scan = NULL;

  // Loop Through the files
  for(int loopA=0;loopA<fileNames.size();loopA++){

    // Read current filename
    readExpansionFile(fileNames[loopA],tot,lengthX,lengthY,lengthZ,minlimits,maxlimits,exp);

    // Build a New Topology from this file
    MRITopology* topo = new MRITopology(tot,lengthX,lengthY,lengthZ,minlimits,maxlimits);

    if(loopA == 0){
      // Assign as the full sequence topology
      topology = topo;
      // Add Scan
      scan = new MRIScan(Times[loopA]);
      scan->rebuildFromExpansion(exp, true);
      addScan(scan);
    }else{
      // Check Compatibility
      if(topology->isCompatibleTopology(topo)){

        // Add Scan From File
        scan = new MRIScan(Times[loopA]);
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
void MRISequence::makeScanDifference(int firstScanID, int secondScanID){
  // Get Scans
  if((firstScanID<0)||(firstScanID>sequence.size())){
    throw MRIException("Cannot Make Difference. Invalid ID for First Scan.");
  }
  if((secondScanID<0)||(secondScanID>sequence.size())){
    throw MRIException("Cannot Make Difference. Invalid ID for Second Scan.");
  }
  MRIScan* firstScan = sequence[firstScanID];
  MRIScan* secondScan = sequence[secondScanID];
  // set the Tolerance Value
  double DistTol = 1.0e-4;
  // Check Compatibility
  bool scansAreCompatible = firstScan->isCompatibleWith(secondScan);
  if(!scansAreCompatible){
    throw MRIException("Scans are not compatible.");
  }
  // Write Progress Message
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
void MRISequence::makeScanAverage(int numberOfMeasures, int firstScanID, int secondScanID){
  // Get Scans
  if((firstScanID<0)||(firstScanID>sequence.size())){
    throw MRIException("Cannot Make Difference. Invalid ID for First Scan.");
  }
  if((secondScanID<0)||(secondScanID>sequence.size())){
    throw MRIException("Cannot Make Difference. Invalid ID for Second Scan.");
  }
  MRIScan* firstScan = sequence[firstScanID];
  MRIScan* secondScan = sequence[secondScanID];
  // Get Distance Tolerance
  double DistTol = 1.0e-4;
  // Check Compatibility
  bool scansAreCompatible = firstScan->isCompatibleWith(secondScan);
  if(!scansAreCompatible){
    throw MRIException("Scans are not compatible.");
  }
  // Write Progress Message
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
string MRISequence::writeStatistics(){
  string myresult = "\n";
  myresult += "--------------------------------\n";
  myresult += "SEQUENCE STATISTICS\n";
  myresult += "--------------------------------\n";
  myresult += "Total Number Of Cells: "+MRIUtils::intToStr(topology->totalCells)+"\n";
  myresult += "--------------------------------\n";
  myresult += "Total Number Of Coordinate Cells\n";
  myresult += "X Direction: "+MRIUtils::intToStr(topology->cellTotals[0])+"\n";
  myresult += "Y Direction: "+MRIUtils::intToStr(topology->cellTotals[1])+"\n";
  myresult += "Z Direction: "+MRIUtils::intToStr(topology->cellTotals[2])+"\n";
  myresult += "Cells Lengths\n";
  myresult += "X Direction - MIN: "+MRIUtils::floatToStr(*min_element(topology->cellLengths[0].begin(),topology->cellLengths[0].end()))+"\n";
  myresult += "X Direction - MAX: "+MRIUtils::floatToStr(*max_element(topology->cellLengths[0].begin(),topology->cellLengths[0].end()))+"\n";
  myresult += "Y Direction - MIN: "+MRIUtils::floatToStr(*min_element(topology->cellLengths[1].begin(),topology->cellLengths[1].end()))+"\n";
  myresult += "Y Direction - MAX: "+MRIUtils::floatToStr(*max_element(topology->cellLengths[1].begin(),topology->cellLengths[1].end()))+"\n";
  myresult += "Z Direction - MIN: "+MRIUtils::floatToStr(*min_element(topology->cellLengths[2].begin(),topology->cellLengths[2].end()))+"\n";
  myresult += "Z Direction - MAX: "+MRIUtils::floatToStr(*max_element(topology->cellLengths[2].begin(),topology->cellLengths[2].end()))+"\n";
  myresult += "--------------------------------\n";
  myresult += "Domain Size\n";
  myresult += "Minimum X: "+MRIUtils::floatToStr(topology->domainSizeMin[0])+"\n";
  myresult += "Maximum X: "+MRIUtils::floatToStr(topology->domainSizeMax[0])+"\n";
  myresult += "Minimum Y: "+MRIUtils::floatToStr(topology->domainSizeMin[1])+"\n";
  myresult += "Maximum Y: "+MRIUtils::floatToStr(topology->domainSizeMax[1])+"\n";
  myresult += "Minimum Z: "+MRIUtils::floatToStr(topology->domainSizeMin[2])+"\n";
  myresult += "Maximum Z: "+MRIUtils::floatToStr(topology->domainSizeMax[2])+"\n";
  myresult += "--------------------------------\n";
  myresult += "Number of Scans: "+MRIUtils::intToStr(sequence.size())+"\n";
  myresult += "--------------------------------\n";
  myresult += "\n";
  // Return String
  return myresult;
}

// ==========================
// READ VTK STRUCTURED POINTS
// ==========================
void MRISequence::readFromASCIISequence(int asciiInputType,const MRIStringVec& asciiFileNames, const MRIDoubleVec& times){

  // Loop over the file names
  MRITopology* topo;
  MRIScan* scan;
  vtkStructuredPointsOptionRecord vtkOptions;
  pltOptionRecord pltOptions;
  for(int loopA=0;loopA<asciiFileNames.size();loopA++){

    // Read Topology
    topo = new MRITopology();
    if(asciiInputType == kInputVTK){
      topo->readFromVTK_ASCII(asciiFileNames[loopA],vtkOptions);  
    }else if(asciiInputType == kInputPLT){
      topo->readFromPLT_ASCII(asciiFileNames[loopA],pltOptions);  
    }
    
    if(loopA == 0){
      // Assign as the full sequence topology
      topology = topo;
      // Add Scan
      scan = new MRIScan(times[loopA]);
      if(asciiInputType == kInputVTK){
        scan->readFromVTK_ASCII(asciiFileNames[loopA],vtkOptions);
      }else if(asciiInputType == kInputPLT){
        scan->readFromPLT_ASCII(asciiFileNames[loopA],pltOptions);
      }    
      addScan(scan);
    }else{
      // Check Compatibility
      if(topology->isCompatibleTopology(topo)){

        // Add Scan From File
        scan = new MRIScan(times[loopA]);
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
    }
  }

  // WRITE STATISTICS
  string CurrentStats = writeStatistics();
  writeSchMessage(CurrentStats);

}
