#include <fstream>
#include <string>
#include <cmath>
#include <mriUtils.h>
#include <mriConstants.h>
#include "mriSequence.h"
#include "schMessages.h"

// Constructor
MRISequence::MRISequence(bool cyclic){
  // Initizalize
  totalScans = 0;
  // Set If Cyclic
  isCyclic = cyclic;
}

// Copy Constructor
MRISequence::MRISequence(MRISequence* copySequence){
  // Copy Scans
  totalScans = copySequence->totalScans;
  // Copy Cyclic Property
  isCyclic = copySequence->isCyclic;
  // Fill with Zero Scans
  for(int loopB=0;loopB<totalScans;loopB++){
    MRIScan* newScan = new MRIScan(*copySequence->GetScan(loopB));
    sequence.push_back(newScan);
  }
}

// Destructor
MRISequence::~MRISequence(){}

// Print the File List Log
void MRISequence::PrintSequenceFiles(std::string outFIleName){
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
void MRISequence::AddScan(MRIScan* scan){
  sequence.push_back(scan);
  totalScans++;
}

// Get a Scan Pointer From the Sequence
MRIScan* MRISequence::GetScan(int scanNumber){
  return sequence[scanNumber];
}

// EXPORT SEQUENCE TO TECPLOT FILE
void MRISequence::ExportToTECPLOT(std::string outfileName){
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("EXPORTING -------------------------------------\n"));
  for(int loopA=0;loopA<totalScans;loopA++){
      sequence[loopA]->ExportToTECPLOT(outfileName,(loopA == 0));
  }
}

// EXCTRACT SINGKLE POINT CURVE IN TIME
void MRISequence::ExtractSinglePointTimeCurve(int cellNumber, int exportQty, std::string fileName){
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
  centrePoint[0] = 0.5 * (sequence[0]->domainSizeMax[0]+sequence[0]->domainSizeMin[0]);
  centrePoint[1] = 0.5 * (sequence[0]->domainSizeMax[1]+sequence[0]->domainSizeMin[1]);
  centrePoint[2] = 0.5 * (sequence[0]->domainSizeMax[2]+sequence[0]->domainSizeMin[2]);  
  // Write List
  fprintf(outFile,"List of Files in Sequence\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    // Get Time
    currTime = sequence[loopA]->scanTime;
    // Get Local Coordinates
    localCoord[0] = sequence[loopA]->cellPoints[cellNumber].position[0] - centrePoint[0];
    localCoord[1] = sequence[loopA]->cellPoints[cellNumber].position[1] - centrePoint[1];
    localCoord[2] = sequence[loopA]->cellPoints[cellNumber].position[2] - centrePoint[2];
    // Check Quantity
    switch(exportQty){
      case kQtyConcentration:
        currValue = sequence[loopA]->cellPoints[cellNumber].concentration;
        break;
      case kQtyVelocityX:
        currValue = sequence[loopA]->cellPoints[cellNumber].velocity[0];
        break;      
      case kQtyVelocityY:
        currValue = sequence[loopA]->cellPoints[cellNumber].velocity[1];
        break;      
      case kQtyVelocityZ:
        currValue = sequence[loopA]->cellPoints[cellNumber].velocity[2];
        break;      
      case kQtyPressGradientX:
        currValue = sequence[loopA]->cellPoints[cellNumber].pressGrad[0];
        break;      
      case kQtyPressGradientY:
        currValue = sequence[loopA]->cellPoints[cellNumber].pressGrad[1];
        break;      
      case kQtyPressGradientZ:
        currValue = sequence[loopA]->cellPoints[cellNumber].pressGrad[2];
        break;      
      case kQtyRelPressure:
        currValue = sequence[loopA]->cellPoints[cellNumber].relPressure;
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
void MRISequence::ComputeRelativePressure(bool doPressureSmoothing){
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("REL PRESSURE COMPUTATION --------------------------\n"));
  // Loop Through Scans
  for(int loopA=0;loopA<totalScans;loopA++){
    // Write Message
    WriteSchMessage(std::string("Computing Relative Pressure - Step "+MRIUtils::IntToStr(loopA+1)+"/"+MRIUtils::IntToStr(totalScans)+"..."));
    // Get The Scan Back
    MRIScan* resultScan = GetScan(loopA);
    int startingCell = resultScan->EvalCentralCell();
    resultScan->EvalRelativePressure(startingCell,0.0);
    if (doPressureSmoothing){
      resultScan->PerformPressureIterations();
    }
    // Done
    WriteSchMessage(std::string("Relative Pressure Computed.\n"));
  }
}

// Read From VOL Sequence File
void MRISequence::ReadFromVolSequence(std::string seqfileName){
  // Var
  MRIScan* myScan;
  std::string buffer;
  
  // Open File
  std::ifstream seqFile;
  seqFile.open(seqfileName.c_str());

  // Loop Through Lines
  while (std::getline(seqFile,buffer))
  {
    if (buffer != ""){
      // Tokenize Line
      std::vector<std::string> ResultArray = MRIUtils::ExctractSubStringFromBufferMS(buffer);    
      // Create Scan Objects
      myScan = new MRIScan(atof(ResultArray[3].c_str()));
      // Files for a single scan in the same row
      myScan->ReadScanFromVOLFiles(ResultArray[0],ResultArray[1],ResultArray[2],ResultArray[3]);    
      // Add to Sequence
      AddScan(myScan);
    }
  }
  // Close File
  seqFile.close();
}
    
// Export to VOL File Set
void MRISequence::ExportToVOL(std::string outfileName){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ExportToVOL(outfileName+"_Step"+MRIUtils::IntToStr(loopA));
  }
}

// Export to Poisson Solver
void MRISequence::ExportForPoisson(string inputFileName,double density,double viscosity,MRIThresholdCriteria* threshold,
                                   bool PPE_IncludeAccelerationTerm,bool PPE_IncludeAdvectionTerm,bool PPE_IncludeDiffusionTerm,bool PPE_IncludeReynoldsTerm,
                                   bool readMuTFromFile, string muTFile, double smagorinskyCoeff){
  string name;
  MRIDoubleMat timeDeriv;
  // MRIDoubleMat reynoldsDeriv;
  for(int loopA=0;loopA<totalScans;loopA++){
    name = inputFileName + "_" + MRIUtils::FloatToStr(loopA);
    if(PPE_IncludeAccelerationTerm){
      printf("Computing Time Derivatives for Scan %d...",loopA);
      EvalScanTimeDerivs(loopA,timeDeriv);
      printf("Done.\n");
    }
    sequence[loopA]->ExportForPoisson(name,density,viscosity,threshold,timeDeriv,
                                      PPE_IncludeAccelerationTerm,PPE_IncludeAdvectionTerm,PPE_IncludeDiffusionTerm,PPE_IncludeReynoldsTerm,
                                      readMuTFromFile,muTFile,smagorinskyCoeff);
  }
}

// EXPORT TO WALL DISTANCE SOLVER
void MRISequence::ExportForDistancing(string inputFileName, MRIThresholdCriteria* threshold){
  string name;
  for(int loopA=0;loopA<totalScans;loopA++){
    name = inputFileName + "_" + MRIUtils::FloatToStr(loopA);
    sequence[loopA]->ExportForDistancing(name,threshold);
  }
}

// EXPORT TO SEQUENCE OF VTK FILES
void MRISequence::ExportToVTK(std::string outfileName,MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ExportToVTK(outfileName,thresholdCriteria);
  }
}

// PHYSICS FILTERING FOR ALL SCANS
void MRISequence::ApplySMPFilter(MRIOptions* options, bool isBC, MRICommunicator* comm){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    // Perform Filter
    sequence[loopA]->applySMPFilter(options,isBC,comm);
    // Update Velocities
    sequence[loopA]->UpdateVelocities();
  }
}

// SAVE INITIAL VELOCITIES
void MRISequence::saveVelocity(){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    // Perform Filter
    sequence[loopA]->saveVelocity();
  }
}


// APPLY THRESHOLDING TO ALL SCANS
void MRISequence::ApplyThresholding(MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ApplyThresholding(thresholdCriteria);
  }
}

// EVAL VORTEX CRITERIA
void MRISequence::EvalVortexCriteria(MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->EvalVortexCriteria(thresholdCriteria);
  }
}

// EVAL VORTICITY
void MRISequence::EvalVorticity(MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->EvalVorticity(thresholdCriteria);
  }
}

// EVAL ENSTROPHY
void MRISequence::EvalEnstrophy(MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->EvalEnstrophy(thresholdCriteria);
  }
}

// EVAL SMP VORTEX CRITERION
void MRISequence::EvalSMPVortexCriteria(){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->EvalSMPVortexCriteria(sequence[loopA]->expansion);
  }
}

// EVAL EXPANSION FILE
void MRISequence::WriteExpansionFile(string fileName){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->WriteExpansionFile(fileName);
  }
}

// Crop All Scans in Sequence
void MRISequence::Crop(double* limitBox){
  WriteSchMessage(std::string("Cropping Sequence..."));
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->Crop(limitBox);
  }
  WriteSchMessage(std::string("Done.\n"));
}

// Scale velocities for all Scans
void MRISequence::ScaleVelocities(double factor){
  WriteSchMessage(std::string("Scaling Velocities..."));
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ScaleVelocities(factor);
  }
  WriteSchMessage(std::string("Done.\n"));
}

// Scale Positions
void MRISequence::ScalePositions(double factor){
  WriteSchMessage(std::string("Scaling Positions..."));
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ScalePositions(factor);
  }
  WriteSchMessage(std::string("Done.\n"));
}

// Add noise to measurements
void MRISequence::applyNoise(double noiseIntensity){
  WriteSchMessage(std::string("Applying Noise..."));
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ApplyGaussianNoise(noiseIntensity);
  }
  WriteSchMessage(std::string("Done.\n"));
}

// =============================
// DISTRIBUTION OF SEQUENCE DATA
// =============================
void MRISequence::DistributeSequenceData(MRICommunicator* comm){
  // Create New Sequence
  for(int loopA=0;loopA<this->totalScans;loopA++){
    sequence[loopA]->DistributeScanData(comm);
  }
}

// =============================
// PERFORM BASIC DATA FILTERING
// =============================
void MRISequence::ApplyMedianFilter(int qtyID,int maxIt,int order,int filterType,MRIThresholdCriteria* threshold){
  // Create New Sequence
  for(int loopA=0;loopA<this->totalScans;loopA++){
    sequence[loopA]->ApplyMedianFilter(qtyID,maxIt,order,filterType,threshold);
  }
}

// ===========================
// CLEAN COMPONENT ON BOUNDARY
// ===========================
void MRISequence::cleanNormalComponentOnBoundary(){
  for(int loopA=0;loopA<this->totalScans;loopA++){
    sequence[loopA]->cleanNormalComponentOnBoundary();
  }
}

// ===============================
// INTERPOLATE BOUNDARY VELOCITIES
// ===============================
void MRISequence::InterpolateBoundaryVelocities(){
  printf("Interpolating Boundary Velocities...\n");
  for(int loopA=0;loopA<this->totalScans;loopA++){
    sequence[loopA]->InterpolateBoundaryVelocities();
  }
}

// ===============================
// EVALUATE REYNOLDS STRESS TENSOR
// ===============================
void MRISequence::EvalReynoldsStresses(MRIThresholdCriteria* threshold){
  for(int loopA=0;loopA<totalScans;loopA++){
    printf("Evaluating Reynolds Stress Tensor for Scan %d...",loopA);
    sequence[loopA]->EvalReynoldsStressComponent(threshold);
    printf("Done.");
  }
}

// ==================
// GET CELL FACE AREA
// ==================
void evalCellAreas(int cellNumber,double* Areas){
  int intCoords[3];
  MapIndexToCoords(cellNumber,intCoords);
  // Get the Three Edge Lengths
  double EdgeX = cellLengths[0][intCoords[0]];
  double EdgeY = cellLengths[1][intCoords[1]];
  double EdgeZ = cellLengths[2][intCoords[2]];
  // Write Results
  Areas[0] = EdgeY * EdgeZ;
  Areas[1] = EdgeX * EdgeZ;
  Areas[2] = EdgeX * EdgeY;
}

// ========================================
// BUILD GRID CONNECTIVITY FOR THE SEQUENCE
// ========================================
void buildCellConnections(){
  // Allocate connections
  cellConnections.resize(totalCellPoints);
  // Loop through the cells
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    // Get Coordinates
    int intCoords[3] = {0};
    int totZSilceNodes = 0;
    int zOffset = 0;
    int yOffset = 0;
    int xOffset = 0;
    int node1,node2,node3,node4;
    int node5,node6,node7,node8;
    // Find Integer Coords
    MapIndexToCoords(loopA,intCoords);
    // Get The Nodes
    // Front Nodes in Z
    // Total Nodes in a Z layer
    totZSilceNodes = (cellTotals[0]+1)*(cellTotals[1]+1);
    zOffset = intCoords[2] * totZSilceNodes;
    yOffset = intCoords[1] * (cellTotals[0]+1);
    xOffset = intCoords[0];
    // Add Node 1
    node1 = (zOffset + yOffset + xOffset);
    cellConnections[loopA].push_back(node1);
    // Add Node 2
    node2 = (zOffset + yOffset + xOffset + 1);
    cellConnections[loopA].push_back(node2);
    // Add Node 3
    node3 = (zOffset + yOffset + xOffset + cellTotals[0] + 1);
    cellConnections[loopA].push_back(node3);
    // Add Node 4
    node4 = (zOffset + yOffset + xOffset + cellTotals[0] + 2);
    cellConnections[loopA].push_back(node4);
    // Change zOffset
    zOffset += (cellTotals[0]+1)*(cellTotals[1]+1);
    // Add Node 5
    node5 = (zOffset + yOffset + xOffset);
    cellConnections[loopA].push_back(node5);
    // Add Node 6
    node6 = (zOffset + yOffset + xOffset + 1);
    cellConnections[loopA].push_back(node6);
    // Add Node 7
    node7 = (zOffset + yOffset + xOffset + cellTotals[0] + 1);
    cellConnections[loopA].push_back(node7);
    // Add Node 8
    node8 = (zOffset + yOffset + xOffset + cellTotals[0] + 2);
    cellConnections[loopA].push_back(node8);
  }
}

// ================================
// GET CELL EXTERNAL NORMAL AT FACE
// ================================
void getExternalFaceNormal(int cellID, int localFaceID, double* extNormal){
  // Get Face Nodes
  std::vector<int> faceIds;
  getFaceConnections(localFaceID,cellConnections[cellID],faceIds);

  int node1Coords[3] = {0};
  int node2Coords[3] = {0};
  int node3Coords[3] = {0};
  // Get the integer coordinates for the first three nodes
  MapIndexToAuxNodeCoords(faceIds[0],node1Coords);
  MapIndexToAuxNodeCoords(faceIds[1],node2Coords);
  MapIndexToAuxNodeCoords(faceIds[2],node3Coords);
  // Get the positions for the first three nodes
  double node1Pos[3] = {0.0};
  double node2Pos[3] = {0.0};
  double node3Pos[3] = {0.0};
  double centreCellPos[3] = {0.0};
  MapAuxCoordsToPosition(node1Coords,node1Pos);
  MapAuxCoordsToPosition(node2Coords,node2Pos);
  MapAuxCoordsToPosition(node3Coords,node3Pos);
  centreCellPos[0] = cellPoints[cellID].position[0];
  centreCellPos[1] = cellPoints[cellID].position[1];
  centreCellPos[2] = cellPoints[cellID].position[2];
  // Get the difference
  double diff1[3] = {0.0};
  double diff2[3] = {0.0};
  double normVec[3] = {0.0};
  for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
    diff1[loopB] = node2Pos[loopB] - node1Pos[loopB];
    diff2[loopB] = node3Pos[loopB] - node2Pos[loopB];
    normVec[loopB] = node1Pos[loopB] - centreCellPos[loopB];
  }
  // Get the normal
  MRIUtils::Do3DExternalProduct(diff1,diff2,extNormal);
  MRIUtils::Normalize3DVector(extNormal);
  // Check Sign
  double sign = 0.0;
  for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
    sign += normVec[loopB] * extNormal[loopB];
  }
  if(sign < 0.0){
    extNormal[0] *= -1.0;
    extNormal[1] *= -1.0;
    extNormal[2] *= -1.0;
  }
}

// =====================
// ADD FACE TO FACE LIST
// =====================
int addToFaceConnections(std::vector<std::vector<mriFace* > > &AuxFirstNodeFaceList, std::vector<int> faceIds){
  mriFace* newFace;
  // Get first node in connectivity
  int firstConnectivityNode = MRIUtils::getMinInt(faceIds);
  // Try to find with the first node list
  bool found = false;
  size_t count = 0;
  while((!found)&&(count<AuxFirstNodeFaceList[firstConnectivityNode].size())){
    found = MRIUtils::isSameIntVector(faceIds,AuxFirstNodeFaceList[firstConnectivityNode][count]->connections);
    // Update
    if(!found){
      count++;
    }
  }
  if(!found){
    // Add to Face List
    faceConnections.push_back(faceIds);
    // Add to AuxFirstNodeFaceList
    newFace = new mriFace;
    newFace->number = faceConnections.size()-1;
    for(size_t loopA=0;loopA<faceIds.size();loopA++){
      newFace->connections.push_back(faceIds[loopA]);
    }
    AuxFirstNodeFaceList[firstConnectivityNode].push_back(newFace);
    // Return
    return (faceConnections.size()-1);
  }else{
    return AuxFirstNodeFaceList[firstConnectivityNode][count]->number;
  }
}

// =====================
// ADD EDGE TO FACE LIST
// =====================
int addToEdgeConnections(vector<vector<mriEdge*> > &AuxFirstNodeEdgeList, int* edgeIds){
  mriEdge* newEdge;
  vector<int> tmp;
  tmp.resize(2);
  tmp[0] = edgeIds[0];
  tmp[1] = edgeIds[1];
  int firstConnectivityNode = 0;

  // Get first node in connectivity
  if(edgeIds[0] < edgeIds[1]){
    firstConnectivityNode = edgeIds[0];
  }else{
    firstConnectivityNode = edgeIds[1];
  }

  // Find it in First Node List
  bool found = false;
  size_t count = 0;
  while((!found)&&(count<AuxFirstNodeEdgeList[firstConnectivityNode].size())){
    found = MRIUtils::isSameIntVector(tmp,AuxFirstNodeEdgeList[firstConnectivityNode][count]->connections);
    // Update
    if(!found){
      count++;
    }
  }
  if(!found){
    // Add to Edge List
    edgeConnections.push_back(tmp);
    // Add to AuxFirstNodeEdgeList
    newEdge = new mriEdge;
    newEdge->number = edgeConnections.size()-1;
    for(size_t loopA=0;loopA<2;loopA++){
      newEdge->connections.push_back(edgeIds[loopA]);
    }
    AuxFirstNodeEdgeList[firstConnectivityNode].push_back(newEdge);
    // Return
    return (edgeConnections.size()-1);
  }else{
    return AuxFirstNodeEdgeList[firstConnectivityNode][count]->number;
  }
}

// =======================
// BUILD FACE CONNECTIVITY
// =======================
void buildFaceConnections(){
  std::vector<int> faceIds;
  std::vector<std::vector<mriFace* > > AuxFirstNodeFaceList;
  int currFace = 0;
  cellFaces.resize(totalCellPoints);
  AuxFirstNodeFaceList.resize(getTotalAuxNodes());
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    for(int loopB=0;loopB<k3DNeighbors;loopB++){
      // Get Face Connections
      getFaceConnections(loopB,cellConnections[loopA],faceIds);
      // Add to Face Connections
      currFace = addToFaceConnections(AuxFirstNodeFaceList,faceIds);
      // Add to Cell Faces
      cellFaces[loopA].push_back(currFace);
    }
  }
}

// =======================
// BUILD EDGE CONNECTIVITY
// =======================
void buildFaceCells(){
  faceCells.resize(faceConnections.size());
  int currFace = 0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    for(size_t loopB=0;loopB<cellFaces[loopA].size();loopB++){
      currFace = cellFaces[loopA][loopB];
      faceCells[currFace].push_back(loopA);
    }
  }
}

// =======================
// BUILD EDGE CONNECTIVITY
// =======================
void buildEdgeConnections(){

  int edgeIds[2];
  std::vector<std::vector<mriEdge*> > AuxFirstNodeEdgeList;
  int currEdge = 0;
  double coeff = 0.0;
  faceEdges.resize(faceConnections.size());
  AuxFirstNodeEdgeList.resize(getTotalAuxNodes());
  // Loop on the total number of faces
  for(size_t loopA=0;loopA<faceConnections.size();loopA++){
    for(int loopB=0;loopB<4;loopB++){
      // Get Face Connections
      getEdgeConnections(loopB,faceConnections[loopA],edgeIds);
      // Add to Face Connections
      currEdge = addToEdgeConnections(AuxFirstNodeEdgeList,edgeIds);
      // Add to Face Edges
      faceEdges[loopA].push_back(currEdge);
    }
  }

  // Build edgeFaces
  edgeFaces.resize(edgeConnections.size());
  for(size_t loopA=0;loopA<faceConnections.size();loopA++){
    for(size_t loopB=0;loopB<faceEdges[loopA].size();loopB++){
      currEdge = faceEdges[loopA][loopB];
      // Get Coefficient
      coeff = getEdgeFaceVortexCoeff(currEdge,loopA);
      if(coeff>0.0){
        edgeFaces[currEdge].push_back(loopA+1);
      }else{
        edgeFaces[currEdge].push_back(-(loopA+1));
      }
    }
  }
}

// =============================
// CREATE SEQUENCE MESH TOPOLOGY
// =============================
void MRISequence::CreateTopology(){
  // Take Time
  float cellConn_BeginTime,cellConn_TotalTime;
  float faceConn_BeginTime,faceConn_TotalTime;
  float faceArea_BeginTime,faceArea_TotalTime;
  float edgeConn_BeginTime,edgeConn_TotalTime;
  float auxNodes_BeginTime,auxNodes_TotalTime;

  // Build Cell Connections
  WriteSchMessage(std::string("Build Cell Connection...\n"));
  cellConn_BeginTime = clock();
  buildCellConnections();  
  cellConn_TotalTime = float( clock () - cellConn_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",cellConn_TotalTime);

  WriteSchMessage(std::string("Build Aux Node Coords...\n"));
  auxNodes_BeginTime = clock();
  buildAuxNodesCoords();
  auxNodes_TotalTime = float( clock () - auxNodes_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",auxNodes_TotalTime);

  // Build Face Connections
  WriteSchMessage(std::string("Build Face Connection...\n"));
  faceConn_BeginTime = clock();
  buildFaceConnections();
  buildFaceCells();
  faceConn_TotalTime = float( clock () - faceConn_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",faceConn_TotalTime);

  // Build Face Area and Face Normal Vector
  WriteSchMessage(std::string("Build Areas and Normals...\n"));
  faceArea_BeginTime = clock();
  buildFaceAreasAndNormals();
  faceArea_TotalTime = float( clock () - faceArea_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",faceArea_TotalTime);

  // Build Edge Connections
  WriteSchMessage(std::string("Build Edge Connections...\n"));
  edgeConn_BeginTime = clock();
  buildEdgeConnections();
  edgeConn_TotalTime = float( clock () - edgeConn_BeginTime ) /  CLOCKS_PER_SEC;
  printf("Executed in %f [s]\n",edgeConn_TotalTime);

  // Completed
  WriteSchMessage(std::string("Topology Creation Completed.\n"));
}

