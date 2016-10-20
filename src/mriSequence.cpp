# include "mriSequence.h"

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

// ==========================
// READ SEQUENCE OF PLT FILES
// ==========================
void MRISequence::readPLTFile(const MRIStringVec& pltFileNames, const MRIDoubleVec& scanTimes , bool DoReorderCells){
  MRITopology* topo;
  MRIScan* scan;
  for(int loopA=0;loopA<pltFileNames.size();loopA++){
    topo = new MRITopology();
    scan = new MRIScan(scanTimes[loopA]);
    readPLTData(pltFileNames[loopA], topo, scan);
    if(loopA == 0){
      // First File: store the topology
      topology = topo;
      // First File: store the scan
      sequence.push_back(scan);
    }else{
      // Other Scan: compare topology
      if(topology->isCompatibleTopology(topo)){
        sequence.push_back(scan);
      }
    }
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
  formDifferenceBinLimits(scanOther,scanRef,pdfQuantity,currInterval,limitBox,numberOfBins,binMin,binMax,binCenters);
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
  for(int loopA=0;loopA<scanRef->totalCellPoints;loopA++){
    // Get Cell Coords
    cellCoord = scanRef->cellPoints[loopA].position;
    // Get quantity
    otherQuantity = scanOther->cellPoints[loopA].getQuantity(pdfQuantity);
    refQuantity = scanRef->cellPoints[loopA].getQuantity(pdfQuantity);
    // Get Value
    currValue = (otherQuantity - refQuantity);
    // Assign Value to Bin
    if (MRIUtils::IsPointInsideBox(cellCoord[0],cellCoord[1],cellCoord[2],limitBox)){
      // COMPLETE
      AssignToBin(currValue,numberOfBins,binMin,binMax,binArray);
    }
  }
  // Normalize
  NormalizeBinArray(numberOfBins,binArray,currInterval);
}

// READ SCAN FROM EXPANSION FILE
void MRISequence::readFromExpansionFile(string fileName,bool applyThreshold, int thresholdType,double thresholdRatio){

  // Allocate Variables
  int tot[3];
  std::vector<double> lengthX;
  std::vector<double> lengthY;
  std::vector<double> lengthZ;
  double minlimits[3];
  double maxlimits[3];
  MRIExpansion* exp = NULL;

  // Read Quantities From File
  readExpansionFile(fileName,tot,lengthX,lengthY,lengthZ,minlimits,maxlimits,exp);

  // SET UP SCAN QUANTITIES
  // CELL TOTALS
  topology->cellTotals[0] = tot[0];
  topology->cellTotals[1] = tot[1];
  topology->cellTotals[2] = tot[2];

  // CELL LENGTHS
  topology->cellLengths.resize(kNumberOfDimensions);
  // X
  for(size_t loopA=0;loopA<lengthX.size();loopA++){
    topology->cellLengths[0].push_back(lengthX[loopA]);
  }
  // Y
  for(size_t loopA=0;loopA<lengthY.size();loopA++){
    topology->cellLengths[1].push_back(lengthY[loopA]);
  }
  // Z
  for(size_t loopA=0;loopA<lengthZ.size();loopA++){
    topology->cellLengths[2].push_back(lengthZ[loopA]);
  }

  // DIMENSIONS
  // MIN
  topology->domainSizeMin[0] = minlimits[0];
  topology->domainSizeMin[1] = minlimits[1];
  topology->domainSizeMin[2] = minlimits[2];
  // MAX
  topology->domainSizeMax[0] = maxlimits[0];
  topology->domainSizeMax[1] = maxlimits[1];
  topology->domainSizeMax[2] = maxlimits[2];
  // MRI EXPANSION
  expansion = new MRIExpansion(exp);

  // APPLY THRESHOLD TO EXPANSION
  if(applyThreshold){
    expansion->ApplyVortexThreshold(thresholdType,thresholdRatio);
  }

  // INITIALIZE VALUES
  hasPressureGradient = false;
  hasRelativePressure = false;
  scanTime = 0.0;
  maxVelModule = 0.0;

  // INITIALIZE SCAN
  topology->totalCells = topology->cellTotals[0]*topology->cellTotals[1]*topology->cellTotals[2];
  // CREATE NEW CELL
  MRICell newCell;
  // INITIALIZE QUANTITIES TO ZERO
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    // ADD IT TO CELL POINTS
    cells.push_back(newCell);
  }

  // INITIALIZE POSITIONS
  int intCoords[3] = {0};
  double Pos[3] = {0.0};
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    mapIndexToCoords(loopA,intCoords);
    mapCoordsToPosition(intCoords,true,Pos);
    cells[loopA].setQuantity(kQtyPositionX,topology->domainSizeMin[0] + Pos[0]);
    cells[loopA].setQuantity(kQtyPositionY,topology->domainSizeMin[1] + Pos[1]);
    cells[loopA].setQuantity(kQtyPositionZ,topology->domainSizeMin[2] + Pos[2]);
  }

  // REORDER MODEL
  reorderScan();

  // REBUILD SCAN
  //RebuildFromExpansion(expansion,true);
  // No Constant Flux
  rebuildFromExpansion(expansion,false);
}

// Make Difference of Scans
// Put the Result Of the Operation in firstScan
void MRISequence::makeScanDifference(int firstScanID, int secondScanID){
  // Get Scans
  if((firstScanID<0)||(firstScanID>sequence.size())){
    throw MRISequenceException("Cannot Make Difference. Invalid ID for First Scan.");
  }
  if((secondScanID<0)||(secondScanID>sequence.size())){
    throw MRISequenceException("Cannot Make Difference. Invalid ID for Second Scan.");
  }
  MRIScan* firstScan = sequence[firstScanID];
  MRIScan* secondScan = sequence[secondScanID];
  // set the Tolerance Value
  double DistTol = 1.0e-4;
  // Check Compatibility
  bool scansAreCompatible = firstScan->isCompatibleWith(secondScan);
  if(!scansAreCompatible){
    throw MRISequenceException("Scans are not compatible.");
  }
  // Write Progress Message
  writeSchMessage("Computing Scan Difference...");  
  // SubTract Velocity Data
  double diffX,diffY,diffZ;
  for(int loopA=0;loopA<firstScan->totalCellPoints;loopA++){
    // Check They have the same Coordinates
    diffX = fabs((firstScan->cellPoints[loopA].position[0]-secondScan->cellPoints[loopA].position[0])/(firstScan->cellPoints[loopA].position[0]));
    diffY = fabs((firstScan->cellPoints[loopA].position[1]-secondScan->cellPoints[loopA].position[1])/(firstScan->cellPoints[loopA].position[1]));
    diffZ = fabs((firstScan->cellPoints[loopA].position[2]-secondScan->cellPoints[loopA].position[2])/(firstScan->cellPoints[loopA].position[2]));
    // Throw exception if the coordinate is different
    if((diffX>DistTol)||(diffY>DistTol)||(diffZ>DistTol)) throw new MRISequenceException("Error In MakeGlobalDataDifference: Different Coords");
    for(int loopB=0;loopB<3;loopB++){
      firstScan->cellPoints[loopA].velocity[loopB] = firstScan->cellPoints[loopA].velocity[loopB]-secondScan->cellPoints[loopA].velocity[loopB];
    }
  }
}

// Make Average
// Put the Result Of the Operation in firstScan
void MRISequence::makeScanAverage(int numberOfMeasures, int firstScanID, int secondScanID){
  // Get Scans
  if((firstScanID<0)||(firstScanID>sequence.size())){
    throw MRISequenceException("Cannot Make Difference. Invalid ID for First Scan.");
  }
  if((secondScanID<0)||(secondScanID>sequence.size())){
    throw MRISequenceException("Cannot Make Difference. Invalid ID for Second Scan.");
  }
  MRIScan* firstScan = sequence[firstScanID];
  MRIScan* secondScan = sequence[secondScanID];
  // Get Distance Tolerance
  double DistTol = 1.0e-4;
  // Check Compatibility
  bool scansAreCompatible = firstScan->isCompatibleWith(secondScan);
  if(!scansAreCompatible){
    throw MRISequenceException("Scans are not compatible.");
  }
  // Write Progress Message
  writeSchMessage("Computing Scan Average...");
  // SubTract Velocity Data
  double diffX,diffY,diffZ;
  for(int loopA=0;loopA<firstScan->totalCellPoints;loopA++){
    // Check They have the same Coordinates
    diffX = fabs((firstScan->cellPoints[loopA].position[0]-secondScan->cellPoints[loopA].position[0])/(firstScan->cellPoints[loopA].position[0]));
    diffY = fabs((firstScan->cellPoints[loopA].position[1]-secondScan->cellPoints[loopA].position[1])/(firstScan->cellPoints[loopA].position[1]));
    diffZ = fabs((firstScan->cellPoints[loopA].position[2]-secondScan->cellPoints[loopA].position[2])/(firstScan->cellPoints[loopA].position[2]));
    if((diffX>DistTol)||(diffY>DistTol)||(diffZ>DistTol))throw new MRISequenceException("Error In MakeGlobalDataAverage: Different Coords");
    for(int loopB=0;loopB<3;loopB++){
      firstScan->cellPoints[loopA].velocity[loopB] =
      firstScan->cellPoints[loopA].velocity[loopB]*((double)(numberOfMeasures-1)/(double)numberOfMeasures)+
      secondScan->cellPoints[loopA].velocity[loopB]*(1.0/(double)numberOfMeasures);
    }
  }
}

// ==========================
// READ VTK STRUCTURED POINTS
// ==========================
void MRISequence::readVTKStructuredPoints(std::string vtkFileName, bool DoReorderCells){

  // Init totalCellPoints
  topology->totalCells = 0;

  // Assign File
  writeSchMessage(std::string("\n"));
  writeSchMessage(std::string("--- READING STRUCTURED POINT FILE\n"));
  writeSchMessage(std::string("\n"));
  std::ifstream vtkFile;
  writeSchMessage(std::string("Open File: ") + vtkFileName + std::string("\n"));
  vtkFile.open(vtkFileName.c_str());

  // Create and initialize vtkOption Vector
  vtkStructuredPointsOptionRecord vtkOptions;
  initVTKStructuredPointsOptions(vtkOptions);

  // Read Through and look for options
  std::vector<std::string> tokenizedString;
  writeSchMessage(std::string("Computing input file size..."));
  int totalLinesInFile = 0;
  std::string Buffer;
  while (std::getline(vtkFile,Buffer)){
    boost::split(tokenizedString, Buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    // Check if you find options
    assignVTKOptions(totalLinesInFile,tokenizedString, vtkOptions);
    // Increase line number
    totalLinesInFile++;
  }
  vtkOptions.dataBlockStart.push_back(totalLinesInFile);
  vtkOptions.dataBlockType.push_back(0);
  vtkOptions.dataBlockRead.push_back(false);

  // Done: Computing Input File Size
  writeSchMessage(std::string("Done.\n"));

  // Check if all properties were defined
  bool fileOK = true;
  for(int loopA=0;loopA<vtkOptions.numDefined;loopA++){
    fileOK = fileOK && vtkOptions.isDefined[loopA];
    if(!fileOK){
      printf("ERROR: DEFINITION %d\n",loopA);
    }
  }
  if(!fileOK){
    writeSchMessage(std::string("ERROR: Invalid VTK File format.\n"));
    writeSchMessage(std::string("\n"));
    exit(1);
  }

  // Print VTK OPTIONS
  printVTKOptions(vtkOptions);

  // Creating Grid Geometry from Options
  writeSchMessage(std::string("Creating Grid Geometry ..."));
  createGridFromVTKStructuredPoints(vtkOptions);
  writeSchMessage(std::string("Done.\n"));

  // Reset File
  vtkFile.clear();
  vtkFile.seekg(0, std::ios::beg);

  // Read Concentration and Velocities
  MRIDoubleVec vtkScalar;
  MRIDoubleMat vtkVector;
  int linesToRead = 0;
  totalLinesInFile = 0;
  for(int loopA=0;loopA<vtkOptions.dataBlockStart.size();loopA++){
    if(vtkOptions.dataBlockRead[loopA]){
      // Go to next block start
      while(totalLinesInFile<vtkOptions.dataBlockStart[loopA]){
        std::getline(vtkFile,Buffer);
        totalLinesInFile++;
      }
      // Read Scalars Or Vectors
      if(vtkOptions.dataBlockType[loopA] == 0){
        // Read Scalars
        writeSchMessage(std::string("Reading Scalars...\n"));
        vtkScalar.clear();
        linesToRead = vtkOptions.dataBlockStart[loopA + 1] - vtkOptions.dataBlockStart[loopA] - 2;
        readVTKScalar(vtkFile,totalLinesInFile,linesToRead,vtkScalar);
        if(vtkScalar.size() != topology->totalCells){
          printf("Scalar Size %d, Total Cells %d\n",(int)vtkScalar.size(), topology->totalCells);
          throw MRIException("ERROR: Total number of scalars differ from number of cells.\n");
        }
      }else{
        // Read Vectors
        writeSchMessage(std::string("Reading Vectors...\n"));
        vtkVector.clear();
        linesToRead = vtkOptions.dataBlockStart[loopA + 1] - vtkOptions.dataBlockStart[loopA];
        readVTKVector(vtkFile,totalLinesInFile,linesToRead,vtkVector);
        if(vtkVector.size() != topology->totalCells){
          printf("Vector Size %d, Total Cells %d\n",(int)vtkVector.size(), topology->totalCells);
          throw MRIException("ERROR: Total number of vectors differs from number of cells.\n");
        }
      }
    }
  }

  // Transfer Scalars and Vectors to Cells
  // Scalars
  for(int loopA=0;loopA<topology->totalCells;loopA++){
    cells[loopA].concentration = vtkScalar[loopA];
  }
  // Vectors
  maxVelModule = 0.0;
  double currModulus = 0.0;
  for(int loopA=0;loopA<topology->totalCells;loopA++){
      cells[loopA].velocity[0] = vtkVector[loopA][0];
      cells[loopA].velocity[1] = vtkVector[loopA][1];
      cells[loopA].velocity[2] = vtkVector[loopA][2];
      currModulus = sqrt((cells[loopA].velocity[0] * cells[loopA].velocity[0]) +
                         (cells[loopA].velocity[1] * cells[loopA].velocity[1]) +
                         (cells[loopA].velocity[2] * cells[loopA].velocity[2]));
      if(currModulus > maxVelModule){
        maxVelModule = currModulus;
      }
  }

  // Finished Reading File
  writeSchMessage(std::string("File reading completed.\n"));

  // Close File
  vtkFile.close();

  // REORDER CELLS
  //if (DoReorderCells){
  //  ReorderScan();
  //}

  // WRITE STATISTICS
  std::string CurrentStats = WriteStatistics();
  writeSchMessage(CurrentStats);

}


