# include "mriTopology.h"

using namespace std;

// Get Local Adjacent Plane
int getLocalAdjacentFace(int localNodeNumber, int totalX, int AdjType){
  switch(AdjType){
    case kfacePlusX:  
      return (((int)localNodeNumber)/((int)totalX)) * (2 * totalX + 1) + (localNodeNumber % totalX) + totalX + 1;
      break;
    case kfaceMinusX: 
      return ((int)localNodeNumber/(int)totalX) * (2*totalX+1) + (localNodeNumber % totalX) + totalX;
    case kfacePlusY : 
      return (((int)localNodeNumber/(int)totalX)+1) * (2*totalX+1) + (localNodeNumber % totalX);
    case kfaceMinusY: 
      return ((int)localNodeNumber/(int)totalX) * (2*totalX+1) + (localNodeNumber % totalX);
    case kfacePlusZ:  
      return localNodeNumber;
    case kfaceMinusZ: 
      return localNodeNumber;
  }
  return -1;
}

// ============================
// READ CELL DATA FROM PLT FILE
// ============================
void readCellsFromPLTFile(string pltFileName, 
                          MRIDoubleVec& domainSizeMin, MRIDoubleVec& domainSizeMax, 
                          double& maxVelModule, 
                          MRIDoubleVec& XCoords,
                          MRIDoubleVec& YCoords,
                          MRIDoubleVec& ZCoords,
                          MRIDoubleMat& gridData){

  ifstream pltFile;
  pltFile.open(pltFileName.c_str());

  // Init Domain Limits
  domainSizeMin[0] =  std::numeric_limits<double>::max();
  domainSizeMin[1] =  std::numeric_limits<double>::max();
  domainSizeMin[2] =  std::numeric_limits<double>::max();
  domainSizeMax[0] = -std::numeric_limits<double>::max();
  domainSizeMax[1] = -std::numeric_limits<double>::max();
  domainSizeMax[2] = -std::numeric_limits<double>::max();  

  // Read Lines
  string Buffer;
  int lineCount = 0;
  MRIStringVec resultArray;
  bool Continue;
  int valueCounter = 0;
  int neededValues = 7;
  MRIDoubleVec LocalVal(7);
  MRIDoubleVec TempVal(7);
  double LocalXCoord = 0.0;
  double LocalYCoord = 0.0;
  double LocalZCoord = 0.0;
  double LocalConc = 0.0;
  double LocalXVel = 0.0;
  double LocalYVel = 0.0;
  double LocalZVel = 0.0;
  double CurrentModule = 0.0;
  MRIDoubleVec tmp;
  while (getline(pltFile,Buffer)){
    // Read Line
    lineCount++;

    // Tokenize Line
    resultArray = MRIUtils::extractSubStringFromBufferMS(Buffer);
    
    // Store Local Structure
    try{
      // Set Continue
      Continue = true;
      // Check Ratio between ResultArray.size, valueCounter, neededValues
      if((int)resultArray.size()+valueCounter < neededValues){
        // Read the whole Result Array
        for(size_t loopA=0;loopA<resultArray.size();loopA++){
          LocalVal[loopA] = atof(resultArray[loopA].c_str());
        }
        Continue = false;
      }else{
        // Read part of the result array
        for(int loopA=0;loopA<neededValues-valueCounter;loopA++){
          LocalVal[loopA] = atof(resultArray[loopA].c_str());
        }
        // Put the rest in temporary array
        for(size_t loopA=0;loopA<resultArray.size()-(neededValues-valueCounter);loopA++){
          TempVal[loopA] = atof(resultArray[loopA].c_str());
        }
        Continue = true;
        // Coords
        LocalXCoord = LocalVal[0];
        LocalYCoord = LocalVal[1];
        LocalZCoord = LocalVal[2];
        // Concentration
        LocalConc = LocalVal[3];
        // Velocity
        LocalXVel = LocalVal[4];
        LocalYVel = LocalVal[5];
        LocalZVel = LocalVal[6];
      }
      Continue = true;
      // Coords
      LocalXCoord = LocalVal[0];
      LocalYCoord = LocalVal[1];
      LocalZCoord = LocalVal[2];
      // Concentration
      LocalConc = LocalVal[3];
      // Velocity
      LocalXVel = LocalVal[4];
      LocalYVel = LocalVal[5];
      LocalZVel = LocalVal[6];
      // Update valueCounter
      valueCounter = ((resultArray.size() + valueCounter) % neededValues);
      // Check Module
      CurrentModule = sqrt((LocalXVel * LocalXVel) + (LocalYVel * LocalYVel) + (LocalZVel * LocalZVel));
      if (CurrentModule>1000.0){
        throw 20;
      }
    }catch (...){
      //Set Continue
      Continue = false;
      std::string outString = "WARNING[*] Error Reading Line: "+MRIUtils::intToStr(lineCount)+"; Line Skipped.\n";
      printf("%s",outString.c_str());
    }
    if(Continue){

      // Update Limits
      // Min
      if (LocalXCoord<domainSizeMin[0]) domainSizeMin[0] = LocalXCoord;
      if (LocalYCoord<domainSizeMin[1]) domainSizeMin[1] = LocalYCoord;
      if (LocalZCoord<domainSizeMin[2]) domainSizeMin[2] = LocalZCoord;
      // Max
      if (LocalXCoord>domainSizeMax[0]) domainSizeMax[0] = LocalXCoord;
      if (LocalYCoord>domainSizeMax[1]) domainSizeMax[1] = LocalYCoord;
      if (LocalZCoord>domainSizeMax[2]) domainSizeMax[2] = LocalZCoord;

      // Update Max Speeds
      if (CurrentModule>maxVelModule) {
        maxVelModule = CurrentModule;
      }

      // Store Node Coords To Find Grid Size
      MRIUtils::insertInList(LocalXCoord,XCoords);
      MRIUtils::insertInList(LocalYCoord,YCoords);
      MRIUtils::insertInList(LocalZCoord,ZCoords);    

      // Store Velocity/Concentrations
      tmp.clear();
      tmp.push_back(LocalXCoord);
      tmp.push_back(LocalYCoord);
      tmp.push_back(LocalZCoord);
      tmp.push_back(LocalConc);
      tmp.push_back(LocalXVel);
      tmp.push_back(LocalYVel);
      tmp.push_back(LocalZVel);
    
      // Add to Vector
      gridData.push_back(tmp);

      // Set Continue
      Continue = true;
    }
  }
}

// ==========================
// GET LOCAL FACE CONNECTIONS
// ==========================
void getFaceConnections(int faceID, vector<int> cellConnections, vector<int> &faceIds){
  faceIds.clear();
  switch(faceID){
    case 0:
      faceIds.push_back(cellConnections[0]);
      faceIds.push_back(cellConnections[1]);
      faceIds.push_back(cellConnections[3]);
      faceIds.push_back(cellConnections[2]);
      break;
    case 1:
      faceIds.push_back(cellConnections[4]);
      faceIds.push_back(cellConnections[6]);
      faceIds.push_back(cellConnections[7]);
      faceIds.push_back(cellConnections[5]);
      break;
    case 2:
      faceIds.push_back(cellConnections[0]);
      faceIds.push_back(cellConnections[2]);
      faceIds.push_back(cellConnections[6]);
      faceIds.push_back(cellConnections[4]);
      break;
    case 3:
      faceIds.push_back(cellConnections[1]);
      faceIds.push_back(cellConnections[5]);
      faceIds.push_back(cellConnections[7]);
      faceIds.push_back(cellConnections[3]);
      break;
    case 4:
      faceIds.push_back(cellConnections[0]);
      faceIds.push_back(cellConnections[4]);
      faceIds.push_back(cellConnections[5]);
      faceIds.push_back(cellConnections[1]);
      break;
    case 5:
      faceIds.push_back(cellConnections[2]);
      faceIds.push_back(cellConnections[3]);
      faceIds.push_back(cellConnections[7]);
      faceIds.push_back(cellConnections[6]);
      break;
  }
}

// ====================================
// GET LOCAL EDGE CONNECTIONS FROM FACE
// ====================================
void getEdgeConnections(int EdgeID, const MRIIntVec& faceConnections, MRIIntVec& edgeIds){
  switch(EdgeID){
    case 0:
      edgeIds[0] = faceConnections[0];
      edgeIds[1] = faceConnections[1];
      break;
    case 1:
      edgeIds[0] = faceConnections[1];
      edgeIds[1] = faceConnections[2];
      break;
    case 2:
      edgeIds[0] = faceConnections[2];
      edgeIds[1] = faceConnections[3];
      break;
    case 3:
      edgeIds[0] = faceConnections[3];
      edgeIds[1] = faceConnections[0];
      break;
  }
}

MRITopology::MRITopology(){
  domainSizeMin.resize(3);
  domainSizeMax.resize(3);
  totalCells = 0;
  cellTotals.resize(3);
  cellLengths.resize(3);
}

// CONSTRUCTOR
MRITopology::MRITopology(const MRIIntVec& totals,
                         const MRIDoubleVec& lengthX,
                         const MRIDoubleVec& lengthY,
                         const MRIDoubleVec& lengthZ,
                         const MRIDoubleVec& minlimits,
                         const MRIDoubleVec& maxlimits){

  // SET UP SCAN QUANTITIES
  // CELL TOTALS
  cellTotals.resize(3);
  cellTotals[0] = totals[0];
  cellTotals[1] = totals[1];
  cellTotals[2] = totals[2];

  // CELL LENGTHS
  cellLengths.resize(kNumberOfDimensions);
  // X
  for(size_t loopA=0;loopA<lengthX.size();loopA++){
    cellLengths[0].push_back(lengthX[loopA]);
  }
  // Y
  for(size_t loopA=0;loopA<lengthY.size();loopA++){
    cellLengths[1].push_back(lengthY[loopA]);
  }
  // Z
  for(size_t loopA=0;loopA<lengthZ.size();loopA++){
    cellLengths[2].push_back(lengthZ[loopA]);
  }

  // DIMENSIONS
  // MIN
  domainSizeMin[0] = minlimits[0];
  domainSizeMin[1] = minlimits[1];
  domainSizeMin[2] = minlimits[2];
  // MAX
  domainSizeMax[0] = maxlimits[0];
  domainSizeMax[1] = maxlimits[1];
  domainSizeMax[2] = maxlimits[2];

  // INITIALIZE SCAN
  totalCells = cellTotals[0]*cellTotals[1]*cellTotals[2];

  // INITIALIZE POSITIONS
  MRIIntVec intCoords(3,0);
  MRIDoubleVec Pos(3,0.0);
  for(int loopA=0;loopA<totalCells;loopA++){
    mapIndexToCoords(loopA,intCoords);
    mapCoordsToPosition(intCoords,true,Pos);
    cellLocations[loopA][0] = domainSizeMin[0] + Pos[0];
    cellLocations[loopA][1] = domainSizeMin[1] + Pos[1];
    cellLocations[loopA][2] = domainSizeMin[2] + Pos[2];
  }
}

// DISTRUCTOR
MRITopology::~MRITopology(){
}

// Map To Cells Coords
void MRITopology::mapIndexToCoords(int index, MRIIntVec& intCoords){
  int CurrentIndex = index;
  intCoords[2] = (int)(CurrentIndex/(cellTotals[0] * cellTotals[1]));
  CurrentIndex = (CurrentIndex-intCoords[2] * cellTotals[0] * cellTotals[1]);
  intCoords[1] = (int)(CurrentIndex / cellTotals[0]);
  CurrentIndex = CurrentIndex-intCoords[1] * cellTotals[0];
  intCoords[0] = CurrentIndex;
}

// Map From Cells Coords
int MRITopology::mapCoordsToIndex(int i, int j, int k){
  // C++ INDEXES ZERO BASED
  return k*(cellTotals[0]*cellTotals[1])+j*(cellTotals[0])+i;
}

// ===================
// GET TOTAL AUX NODES
// ===================
int MRITopology::getTotalAuxNodes(){
  return (cellTotals[0] + 1)*(cellTotals[1] + 1)*(cellTotals[2] + 1);
}

// ===============
// GET TOTAL FACES
// ===============
int MRITopology::getTotalFaces(){
  return cellTotals[0] * cellTotals[1] * (cellTotals[2] + 1) +
         cellTotals[1] * cellTotals[2] * (cellTotals[0] + 1) +
         cellTotals[2] * cellTotals[0] * (cellTotals[1] + 1);
}

// ======================
// GET VORTEX COEFFICIENT
// ======================
double MRITopology::getEdgeFaceVortexCoeff(int edgeID, int faceID){
  
  MRIDoubleVec edgeDirVector(3);
  MRIDoubleVec edgeFaceVector(3);
  MRIDoubleVec resVec(3);
  double res = 0.0;
  
  // Get Vectors
  getEdgeDirection(edgeID,edgeDirVector);
  getEdgeToFaceDirection(edgeID,faceID,edgeFaceVector);

  // Eval Vector Product
  resVec[0] = edgeDirVector[1] * edgeFaceVector[2] - edgeFaceVector[1] * edgeDirVector[2];
  resVec[1] = edgeDirVector[2] * edgeFaceVector[0] - edgeFaceVector[2] * edgeDirVector[0];
  resVec[2] = edgeDirVector[0] * edgeFaceVector[1] - edgeFaceVector[0] * edgeDirVector[1];
  double modulus = (resVec[0] * resVec[0] + resVec[1] * resVec[1] + resVec[2] * resVec[2]);
  resVec[0] = resVec[0]/modulus;
  resVec[1] = resVec[1]/modulus;
  resVec[2] = resVec[2]/modulus;
  // Get Sign
  res = resVec[0] * faceNormal[faceID][0] + 
        resVec[1] * faceNormal[faceID][1] + 
        resVec[2] * faceNormal[faceID][2];
  
  return round(res);
}

// ==================
// GET EDGE DIRECTION
// ==================
void MRITopology::getEdgeDirection(int edgeID, MRIDoubleVec& edgeDirVector){
  int node1 = 0;
  int node2 = 0;
  MRIDoubleVec node1Pos(3,0.0);
  MRIDoubleVec node2Pos(3,0.0);

  // Get The Two Nodes
  node1 = edgeConnections[edgeID][0];
  node2 = edgeConnections[edgeID][1];

  // Eval Auxiliary Node Coordinates
  node1Pos[0] = auxNodesCoords[node1][0];
  node1Pos[1] = auxNodesCoords[node1][1];
  node1Pos[2] = auxNodesCoords[node1][2];
  node2Pos[0] = auxNodesCoords[node2][0];
  node2Pos[1] = auxNodesCoords[node2][1];
  node2Pos[2] = auxNodesCoords[node2][2];

  // Get the versor
  edgeDirVector[0] = node1Pos[0] - node2Pos[0];
  edgeDirVector[1] = node1Pos[1] - node2Pos[1];
  edgeDirVector[2] = node1Pos[2] - node2Pos[2];
  double modulus = (edgeDirVector[0] * edgeDirVector[0] + 
                    edgeDirVector[1] * edgeDirVector[1] + 
                    edgeDirVector[2] * edgeDirVector[2]);

  edgeDirVector[0] = fabs(edgeDirVector[0]/modulus);
  edgeDirVector[1] = fabs(edgeDirVector[1]/modulus);
  edgeDirVector[2] = fabs(edgeDirVector[2]/modulus);
}

// ==========================
// GET EDGE TO FACE DIRECTION
// ==========================
void MRITopology::getEdgeToFaceDirection(int edgeID, int faceID, MRIDoubleVec& edgeFaceVector){
  // Declare
  MRIDoubleVec ec(3,0.0);
  MRIDoubleVec fc(3,0.0);

  // Get Edge Center
  getEdgeCenter(edgeID,ec);

  // Get Face Center
  getFaceCenter(faceID,fc);

  // Get the versor
  edgeFaceVector[0] = fc[0] - ec[0];
  edgeFaceVector[1] = fc[1] - ec[1];
  edgeFaceVector[2] = fc[2] - ec[2];
  double modulus = (edgeFaceVector[0] * edgeFaceVector[0] + edgeFaceVector[1] * edgeFaceVector[1] + edgeFaceVector[2] * edgeFaceVector[2]);
  edgeFaceVector[0] = edgeFaceVector[0]/modulus;
  edgeFaceVector[1] = edgeFaceVector[1]/modulus;
  edgeFaceVector[2] = edgeFaceVector[2]/modulus;
}

// ===============
// GET EDGE CENTER
// ===============
void MRITopology::getEdgeCenter(int edgeID, MRIDoubleVec& ec){
  int node1 = 0;
  int node2 = 0;
  double node1Pos[3] = {0.0};
  double node2Pos[3] = {0.0};

  // Get The Two Nodes
  node1 = edgeConnections[edgeID][0];
  node2 = edgeConnections[edgeID][1];

  // Eval Auxiliary Node Coordinates
  node1Pos[0] = auxNodesCoords[node1][0];
  node1Pos[1] = auxNodesCoords[node1][1];
  node1Pos[2] = auxNodesCoords[node1][2];
  node2Pos[0] = auxNodesCoords[node2][0];
  node2Pos[1] = auxNodesCoords[node2][1];
  node2Pos[2] = auxNodesCoords[node2][2];

  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    ec[loopA] = 0.5*(node1Pos[loopA] + node2Pos[loopA]);
  }
}

// ===============
// GET FACE CENTER
// ===============
void MRITopology::getFaceCenter(int faceID, MRIDoubleVec& fc){
  int currNode = 0;
  MRIDoubleVec pos(3,0.0);

  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    fc[loopA] = 0.0;
  }

  for(size_t loopA=0;loopA<faceConnections[faceID].size();loopA++){
    currNode = faceConnections[faceID][loopA];
    pos[0] = auxNodesCoords[currNode][0];
    pos[1] = auxNodesCoords[currNode][1];
    pos[2] = auxNodesCoords[currNode][2];
    for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
      fc[loopA] += pos[loopA];
    }
  }
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    fc[loopA] /= (double)faceConnections[faceID].size();
  }
}

// ===================================
// CREATE MATRIX WITH AUX NODES COORDS
// ===================================
void MRITopology::buildAuxNodesCoords(){
  MRIDoubleVec nodePos(3,0.0);
  MRIDoubleVec nodePosVec(3);
  int totAuxNodes = getTotalAuxNodes();
  for(int loopA=0;loopA<totAuxNodes;loopA++){
    getAuxNodeCoordinates(loopA,nodePos);
    nodePosVec[0] = nodePos[0];
    nodePosVec[1] = nodePos[1];
    nodePosVec[2] = nodePos[2];
    auxNodesCoords.push_back(nodePosVec);
  }
}

// ====================================
// GET COORDINATEDS FOR AUXILIARY NODES
// ====================================
void MRITopology::getAuxNodeCoordinates(int nodeNum, MRIDoubleVec& pos){
  MRIIntVec intAuxCoords(3,0.0);
  // Map To Integer Coordinates
  mapIndexToAuxNodeCoords(nodeNum,intAuxCoords);
  // Map To Spatial Position
  mapAuxCoordsToPosition(intAuxCoords,pos);
}

// =======================
// MAP TO AUX CELLS COORDS
// =======================
void MRITopology::mapIndexToAuxNodeCoords(int index, MRIIntVec& intCoords){
  int CurrentIndex = index;
  int totalX = cellTotals[0] + 1;
  int totalY = cellTotals[1] + 1;
  intCoords[2] = (int)(CurrentIndex/(totalX*totalY));
  CurrentIndex = (CurrentIndex-intCoords[2]*totalX*totalY);
  intCoords[1] = (int)(CurrentIndex/totalX);
  CurrentIndex = CurrentIndex-intCoords[1]*totalX;
  intCoords[0] = CurrentIndex;
}

// ========================================
// MAP AUXILIARY INTEGER COORDS TO POSITION
// ========================================
void MRITopology::mapAuxCoordsToPosition(const MRIIntVec& auxCoords, MRIDoubleVec& pos){
  // Loop on the three dimensions
  for(int loopA=0;loopA<3;loopA++){
    pos[loopA] = 0.0;
    for(int loopB=0;loopB<auxCoords[loopA];loopB++){
      pos[loopA] += cellLengths[loopA][loopB];
    }
  }
  for(int loopA=0;loopA<3;loopA++){
    pos[loopA] = pos[loopA] + domainSizeMin[loopA] - 0.5*cellLengths[loopA][0];
  }
}

// ========================================
// BUILD GRID CONNECTIVITY FOR THE SEQUENCE
// ========================================
void MRITopology::buildCellConnections(){

  MRIIntVec intCoords(3,0);
  int totZSilceNodes = 0;
  int zOffset = 0;
  int yOffset = 0;
  int xOffset = 0;
  int node1,node2,node3,node4;
  int node5,node6,node7,node8;

  // Allocate connections
  cellConnections.resize(totalCells);
  // Loop through the cells
  for(int loopA=0;loopA<totalCells;loopA++){
    // Find Integer Coords
    mapIndexToCoords(loopA,intCoords);
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
    zOffset += (cellTotals[0] + 1)*(cellTotals[1] + 1);
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

// =======================
// BUILD FACE CONNECTIVITY
// =======================
void MRITopology::buildFaceConnections(){
  MRIIntVec faceIds;
  vector<vector<mriFace* > > AuxFirstNodeFaceList;
  int currFace = 0;
  cellFaces.resize(totalCells);
  AuxFirstNodeFaceList.resize(getTotalAuxNodes());
  for(int loopA=0;loopA<totalCells;loopA++){
    for(int loopB=0;loopB<k3DNeighbors;loopB++){
      // Get Face Connections
      getFaceConnections(loopB,cellConnections[loopA],faceIds);
      // Add to Face Connections
      currFace = addToFaceConnections(faceIds,AuxFirstNodeFaceList);
      // Add to Cell Faces
      cellFaces[loopA].push_back(currFace);
    }
  }
}

// =====================
// ADD FACE TO FACE LIST
// =====================
int MRITopology::addToFaceConnections(const MRIIntVec& faceIds, vector<vector<mriFace* > >& AuxFirstNodeFaceList){
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
int MRITopology::addToEdgeConnections(const MRIIntVec& edgeIds,vector<vector<mriEdge*> >& AuxFirstNodeEdgeList){
  mriEdge* newEdge;
  MRIIntVec tmp;
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
// BUILD EDGE CONNECTIVITY
// =======================
void MRITopology::buildFaceCells(){
  faceCells.resize(faceConnections.size());
  int currFace = 0;
  for(int loopA=0;loopA<totalCells;loopA++){
    for(size_t loopB=0;loopB<cellFaces[loopA].size();loopB++){
      currFace = cellFaces[loopA][loopB];
      faceCells[currFace].push_back(loopA);
    }
  }
}

// =======================
// BUILD EDGE CONNECTIVITY
// =======================
void MRITopology::buildEdgeConnections(){

  MRIIntVec edgeIds(2,0);
  vector<vector<mriEdge*> > AuxFirstNodeEdgeList;
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
      currEdge = addToEdgeConnections(edgeIds,AuxFirstNodeEdgeList);
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

// ======================
// BUILD FACE AREA VECTOR
// ======================
void MRITopology::buildFaceAreasAndNormals(){
  double prod = 0.0;
  faceArea.resize(faceConnections.size());
  faceNormal.resize(faceConnections.size());
  // Declare
  MRIIntVec node1Coords(3,0);
  MRIIntVec node2Coords(3,0);
  MRIIntVec node3Coords(3,0);
  MRIDoubleVec node1Pos(3,0.0);
  MRIDoubleVec node2Pos(3,0.0);
  MRIDoubleVec node3Pos(3,0.0);
  MRIDoubleVec diff1(3,0.0);
  MRIDoubleVec diff2(3,0.0);
  MRIDoubleVec diff3(3,0.0);
  MRIDoubleVec currNormal(3,0.0);
  double d1 = 0.0;
  double d2 = 0.0;
  double innerProd = 0.0;
  for(size_t loopA=0;loopA<faceConnections.size();loopA++){
    // Get the integer coordinates for the first three nodes
    mapIndexToAuxNodeCoords(faceConnections[loopA][0],node1Coords);
    mapIndexToAuxNodeCoords(faceConnections[loopA][1],node2Coords);
    mapIndexToAuxNodeCoords(faceConnections[loopA][2],node3Coords);
    // Get the positions for the first three nodes
    mapAuxCoordsToPosition(node1Coords,node1Pos);
    mapAuxCoordsToPosition(node2Coords,node2Pos);
    mapAuxCoordsToPosition(node3Coords,node3Pos);   
    // Get the difference
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      diff1[loopB] = node2Pos[loopB] - node1Pos[loopB];
      diff2[loopB] = node3Pos[loopB] - node2Pos[loopB];
      diff3[loopB] = node1Pos[loopB] - cellLocations[faceCells[loopA][0]][loopB];
    }

    innerProd = 0.0;
    for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
      innerProd += diff1[loopB] * diff2[loopB];
    }
    if(fabs(innerProd) > kMathZero){
      throw MRIException("ERROR: Face sides are not Orthogonal in buildFaceAreasAndNormals\n");
    }
    d1 = MRIUtils::do3DEucNorm(diff1);
    d2 = MRIUtils::do3DEucNorm(diff2);
    // Evaluate Face Area
    faceArea[loopA] = d1 * d2;
    // Get the normal
    MRIUtils::do3DExternalProduct(diff1,diff2,currNormal);
    MRIUtils::normalize3DVector(currNormal);
    //printf("NODE 1 POS: %f %f %f\n",node1Pos[0],node1Pos[1],node1Pos[2]);
    //printf("CELL CENTRE POS: %f %f %f\n",cellPoints[faceCells[loopA][0]].position[0],cellPoints[faceCells[loopA][0]].position[1],cellPoints[faceCells[loopA][0]].position[2]);
    //getchar();
    if(faceCells[loopA].size() == 1){
      prod = 0.0;
      for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
        prod += currNormal[loopB] * diff3[loopB];
      }
      if(prod > kMathZero){
        //printf("NODE 1 POS: %f %f %f\n",node1Pos[0],node1Pos[1],node1Pos[2]);
        //printf("CELL CENTRE POS: %f %f %f\n",cellPoints[faceCells[loopA][0]].position[0],cellPoints[faceCells[loopA][0]].position[1],cellPoints[faceCells[loopA][0]].position[2]);
        //printf("FLIPPED! Prod: %e, %f %f %f, %f %f %f\n",prod,currNormal[0],currNormal[1],currNormal[2],diff3[0],diff3[1],diff3[2]);
        //getchar();
        currNormal[0] *= - 1.0;
        currNormal[1] *= - 1.0;
        currNormal[2] *= - 1.0;
      }
    }
    faceNormal[loopA].push_back(currNormal[0]);
    faceNormal[loopA].push_back(currNormal[1]);
    faceNormal[loopA].push_back(currNormal[2]);
  }
}

// ===============
// SCALE POSITIONS
// ===============
void MRITopology::scalePositions(const MRIDoubleVec& origin, double factor){
  // SCALE CELL LENGTHS
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    for(size_t loopB=0;loopB<cellLengths[loopA].size();loopB++){
      cellLengths[loopA][loopB] *= factor;
    }
  }
  // SCALE POSITIONS
  for(int loopA=0;loopA<totalCells;loopA++){
    cellLocations[loopA][0] = origin[0] + (cellLocations[loopA][0] - origin[0]) * factor;
    cellLocations[loopA][1] = origin[1] + (cellLocations[loopA][1] - origin[1]) * factor;
    cellLocations[loopA][2] = origin[2] + (cellLocations[loopA][2] - origin[2]) * factor;
  }
  // SCALE FACE AREA
  for(int loopA=0;loopA<faceConnections.size();loopA++){
    faceArea[loopA] = faceArea[loopA] * factor * factor;
  }
  // SCALE DOMAIN DIMENSIONS
  // Max
  domainSizeMax[0] = origin[0] + (domainSizeMax[0] - origin[0]) * factor;
  domainSizeMax[1] = origin[1] + (domainSizeMax[1] - origin[1]) * factor;
  domainSizeMax[2] = origin[2] + (domainSizeMax[2] - origin[2]) * factor;
  // Min
  domainSizeMin[0] = origin[0] + (domainSizeMin[0] - origin[0]) * factor;
  domainSizeMin[1] = origin[1] + (domainSizeMin[1] - origin[1]) * factor;
  domainSizeMin[2] = origin[2] + (domainSizeMin[2] - origin[2]) * factor;
}

// =========
// CROP SCAN
// =========
void MRITopology::crop(const MRIDoubleVec& limitBox, MRIBoolVec& indexes){
  
  // Clear Index Vector
  indexes.clear();

  MRIDoubleVec tmp(3);
  MRIDoubleMat newCellLocations;

  // Count The Number Of Cells Remaining
  int remainingCells = 0;
  for(int loopA=0;loopA<totalCells;loopA++){
    if (MRIUtils::isPointInsideBox(cellLocations[loopA][0],
                                   cellLocations[loopA][1],
                                   cellLocations[loopA][2],
                                   limitBox)){

      tmp[0] = cellLocations[loopA][0];
      tmp[1] = cellLocations[loopA][1];
      tmp[2] = cellLocations[loopA][2];
      newCellLocations.push_back(tmp);

      indexes.push_back(true);
      remainingCells++;
    }else{
      indexes.push_back(false);
    }
  }
  // Allocate New Cellpoints
  cellLocations.resize(remainingCells);
  for(int loopA=0;loopA<remainingCells;loopA++){
    cellLocations[loopA].resize(3);
  }

  // Fill Temporary Cells
  int tempCount = 0;
  for(int loopA=0;loopA<totalCells;loopA++){
    if (indexes[loopA]){
      // Velocity
      cellLocations[tempCount][0] = newCellLocations[loopA][0];
      cellLocations[tempCount][1] = newCellLocations[loopA][1];
      cellLocations[tempCount][2] = newCellLocations[loopA][2];
      // Update Counter
      tempCount++;
    }
  }
}

// ============================================
// Create Grid from VTK Structured Scan Options
// ============================================
void MRITopology::createGridFromVTKStructuredPoints(const vtkStructuredPointsOptionRecord& opts){
  // Assign cell totals
  cellTotals.resize(3);
  cellTotals[0] = opts.dimensions[0];
  cellTotals[1] = opts.dimensions[1];
  cellTotals[2] = opts.dimensions[2];
  // Assign total number of cells
  totalCells = cellTotals[0] * cellTotals[1] * cellTotals[2];
  cellLocations.resize(totalCells);
  for(int loopA=0;loopA<totalCells;loopA++){
    cellLocations[loopA].resize(3);
  }
  // Assign Cell spacing
  cellLengths.resize(3);
  cellLengths[0].resize(cellTotals[0]);
  cellLengths[1].resize(cellTotals[1]);
  cellLengths[2].resize(cellTotals[2]);
  for(int loopA=0;loopA<cellTotals[0];loopA++){
    cellLengths[0][loopA] = opts.spacing[0];
  }
  for(int loopA=0;loopA<cellTotals[1];loopA++){
    cellLengths[1][loopA] = opts.spacing[1];
  }
  for(int loopA=0;loopA<cellTotals[2];loopA++){
    cellLengths[2][loopA] = opts.spacing[2];
  }
  // Set domain size
  // Min
  domainSizeMin.resize(3);
  domainSizeMin[0] = opts.origin[0];
  domainSizeMin[1] = opts.origin[1];
  domainSizeMin[2] = opts.origin[2];
  // Max
  domainSizeMax.resize(3);
  domainSizeMax[0] = opts.origin[0] + (opts.dimensions[0]-1) * opts.spacing[0];
  domainSizeMax[1] = opts.origin[1] + (opts.dimensions[1]-1) * opts.spacing[1];
  domainSizeMax[2] = opts.origin[2] + (opts.dimensions[2]-1) * opts.spacing[2];
  // Allocate the cells
  //cellPoints.reserve(totalCellPoints);
  int count = 0;
  // Fill the position vectors
  double locCoordX = opts.origin[0];
  double locCoordY = opts.origin[1];
  double locCoordZ = opts.origin[2];
  for(int loopA=0;loopA<cellTotals[2];loopA++){
    locCoordY = opts.origin[1];
    for(int loopB=0;loopB<cellTotals[1];loopB++){
      locCoordX = opts.origin[0];
      for(int loopC=0;loopC<cellTotals[0];loopC++){
        // Set Cell positions
        cellLocations[count][0] = locCoordX;
        cellLocations[count][1] = locCoordY;
        cellLocations[count][2] = locCoordZ;
        count++;
        locCoordX += opts.spacing[0];
      }
      locCoordY += opts.spacing[1];
    }
    locCoordZ += opts.spacing[2];
  }
}

// ================================
// GET CELL EXTERNAL NORMAL AT FACE
// ================================
void MRITopology::getExternalFaceNormal(int cellID, int localFaceID, MRIDoubleVec& extNormal){
  // Get Face Nodes
  MRIIntVec faceIds;
  getFaceConnections(localFaceID,cellConnections[cellID],faceIds);

  MRIIntVec node1Coords(3,0);
  MRIIntVec node2Coords(3,0);
  MRIIntVec node3Coords(3,0);
  // Get the integer coordinates for the first three nodes
  mapIndexToAuxNodeCoords(faceIds[0],node1Coords);
  mapIndexToAuxNodeCoords(faceIds[1],node2Coords);
  mapIndexToAuxNodeCoords(faceIds[2],node3Coords);
  
  // Get the positions for the first three nodes
  MRIDoubleVec node1Pos(3,0.0);
  MRIDoubleVec node2Pos(3,0.0);
  MRIDoubleVec node3Pos(3,0.0);
  MRIDoubleVec centreCellPos(3,0.0);
  mapAuxCoordsToPosition(node1Coords,node1Pos);
  mapAuxCoordsToPosition(node2Coords,node2Pos);
  mapAuxCoordsToPosition(node3Coords,node3Pos);
  centreCellPos[0] = cellLocations[cellID][0];
  centreCellPos[1] = cellLocations[cellID][1];
  centreCellPos[2] = cellLocations[cellID][2];
  // Get the difference
  MRIDoubleVec diff1(3,0.0);
  MRIDoubleVec diff2(3,0.0);
  MRIDoubleVec normVec(3,0.0);
  for(int loopB=0;loopB<kNumberOfDimensions;loopB++){
    diff1[loopB] = node2Pos[loopB] - node1Pos[loopB];
    diff2[loopB] = node3Pos[loopB] - node2Pos[loopB];
    normVec[loopB] = node1Pos[loopB] - centreCellPos[loopB];
  }
  // Get the normal
  MRIUtils::do3DExternalProduct(diff1,diff2,extNormal);
  MRIUtils::normalize3DVector(extNormal);
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

// =========================
// CHECK COMPATIBLE TOPOLOGY
// =========================
bool MRITopology::isCompatibleTopology(MRITopology* topo){
  bool result = true;

  result = result && (totalCells == topo->totalCells);

  result = result && (cellTotals[0] == topo->cellTotals[0]);
  result = result && (cellTotals[1] == topo->cellTotals[1]);
  result = result && (cellTotals[2] == topo->cellTotals[2]);

  for(int loopA=0;loopA<cellTotals[0];loopA++){
    result = result && (fabs(cellLengths[0][loopA]-topo->cellLengths[0][loopA]) < kMathZero);
  }

  for(int loopA=0;loopA<cellTotals[1];loopA++){
    result = result && (fabs(cellLengths[1][loopA]-topo->cellLengths[1][loopA]) < kMathZero);
  }

  for(int loopA=0;loopA<cellTotals[2];loopA++){
    result = result && (fabs(cellLengths[2][loopA]-topo->cellLengths[2][loopA]) < kMathZero);
  }

  result = result && (fabs(domainSizeMin[0]-topo->domainSizeMin[0]) < kMathZero);
  result = result && (fabs(domainSizeMin[1]-topo->domainSizeMin[1]) < kMathZero);
  result = result && (fabs(domainSizeMin[2]-topo->domainSizeMin[2]) < kMathZero);

  result = result && (fabs(domainSizeMax[0]-topo->domainSizeMax[0]) < kMathZero);
  result = result && (fabs(domainSizeMax[1]-topo->domainSizeMax[1]) < kMathZero);
  result = result && (fabs(domainSizeMax[2]-topo->domainSizeMax[2]) < kMathZero);

  return result;
}

// ==============================
// MAP INTEGER COORDS TO POSITION
// ==============================
void MRITopology::mapCoordsToPosition(const MRIIntVec& coords, bool addMeshMinima, MRIDoubleVec& pos){
  // Loop on the three dimensions
  for(int loopA=0;loopA<3;loopA++){
    pos[loopA] = 0.0;
    for(int loopB=1;loopB<(coords[loopA]+1);loopB++){
      pos[loopA] += 0.5*(cellLengths[loopA][loopB-1] + cellLengths[loopA][loopB]);
    }
  }
  if(addMeshMinima){
    for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
      pos[loopA] += domainSizeMin[loopA];
    }
  }
}

// ============================
// READ TOPOLOGY FROM ASCII VTK
// ============================
void MRITopology::readFromVTK_ASCII(string vtkFileName, vtkStructuredPointsOptionRecord& vtkOptions){
  // Write Progress
  writeSchMessage(string("Reading Topology From: ") + vtkFileName + string("\n"));

  // Assign File
  ifstream vtkFile;
  vtkFile.open(vtkFileName.c_str());

  // Create and initialize vtkOption Vector
  initVTKStructuredPointsOptions(vtkOptions);

  // Read Through and look for options
  MRIStringVec tokenizedString;
  int totalLinesInFile = 0;
  string Buffer;
  while(getline(vtkFile,Buffer)){
    boost::split(tokenizedString, Buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    // Check if you find options
    assignVTKOptions(totalLinesInFile,tokenizedString, vtkOptions);
    // Increase line number
    totalLinesInFile++;
  }

  vtkOptions.dataBlockStart.push_back(totalLinesInFile);
  vtkOptions.dataBlockType.push_back(0);
  vtkOptions.dataBlockRead.push_back(false);

  // CHECK IF ALL PROPERTIES WERE DEFINED
  bool fileOK = true;
  for(int loopA=0;loopA<vtkOptions.numDefined;loopA++){
    fileOK = fileOK && vtkOptions.isDefined[loopA];
    if(!fileOK){
      printf("ERROR: DEFINITION %d\n",loopA);
    }
  }
  if(!fileOK){
    writeSchMessage(string("ERROR: Invalid VTK File format.\n"));
    writeSchMessage(string("\n"));
    exit(1);
  }

  // Creating Grid Geometry from Options
  createGridFromVTKStructuredPoints(vtkOptions);
}

// ================================
// READ TOPOLOGY FROM ASCII TECPLOT
// ================================
void MRITopology::readFromPLT_ASCII(string pltFileName, pltOptionRecord& pltOptions){
  
  // Write Progress
  writeSchMessage(std::string("Reading topology from : ") + pltFileName + std::string("\n"));

  // Init Line Count
  int lineCount = 0;
  totalCells = 0;

  // READ PLT FILE HEADER
  ifstream pltFile;
  pltFile.open(pltFileName.c_str());
  MRIStringVec tokenizedString;
  bool foundheader = false;
  bool areAllFloats = false;
  int headerCount = 0;
  int totalLinesInFile = 0;
  std::string Buffer;
  while (getline(pltFile,Buffer)){
    if(!foundheader){
      boost::trim(Buffer);
      boost::split(tokenizedString, Buffer, boost::is_any_of("= ,"), boost::token_compress_on);
      areAllFloats = true;
      assignPLTOptions(tokenizedString, pltOptions);
      for(size_t loopA=0;loopA<tokenizedString.size();loopA++){
        areAllFloats = (areAllFloats && (MRIUtils::isFloat(tokenizedString[loopA])));
      }
      foundheader = areAllFloats;
      headerCount++;
    }
    // Increase cell and line number
    totalLinesInFile++;
  }

  cellTotals.resize(3);
  cellTotals[0] = pltOptions.i;
  cellTotals[1] = pltOptions.j;
  cellTotals[2] = pltOptions.k;

  // Read Cells From PLT File
  double maxVelModule;
  MRIDoubleVec XCoords;
  MRIDoubleVec YCoords;
  MRIDoubleVec ZCoords;
  MRIDoubleMat gridData;  
  readCellsFromPLTFile(pltFileName, 
                       domainSizeMin,domainSizeMax,
                       maxVelModule, 
                       XCoords,YCoords,ZCoords,
                       gridData);

  
  // Resize CellLenghts: UNIFORM CASE
  cellLengths.resize(3);
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    cellLengths[loopA].resize(cellTotals[loopA]);
    if(cellTotals[loopA] == 1){
      cellLengths[loopA][0] = 1.0;
    }else{
      for(int loopB=0;loopB<cellTotals[loopA];loopB++){
        cellLengths[loopA][loopB] = fabs(domainSizeMax[loopA]-domainSizeMin[loopA])/(cellTotals[loopA]-1);
      }
    }
  }

  // Close File
  pltFile.close();
}

// =========================
// Get Global Adjacent Plane
// =========================
int MRITopology::getAdjacentFace(int globalNodeNumber /*Already Ordered Globally x-y-z*/, int AdjType){
  
  // Get The Z Coord
  double currentZCoord = cellLocations[globalNodeNumber][2] - domainSizeMin[2];
  
  // Find The Node Number in The Current Plane
  int ZCompleteLevels = MRIUtils::findHowMany(currentZCoord,cellLengths[2]);
  int localNodeNumber = globalNodeNumber - ZCompleteLevels * cellTotals[0] * cellTotals[1];

  // Find The Adjacent face in the Current Plane
  int localFaceNumber = getLocalAdjacentFace(localNodeNumber, cellTotals[0], AdjType);

  // Map to Global Face Numbering: Differentiate if in Z
  if(AdjType == kfacePlusZ){
    return localFaceNumber+(ZCompleteLevels+1)*cellTotals[0]*cellTotals[1]+
                           (ZCompleteLevels+1)*(cellTotals[0]*(cellTotals[1]+1)+
                           (cellTotals[0]+1)*cellTotals[1]);
    
  }else if(AdjType == kfaceMinusZ){
    return localFaceNumber+(ZCompleteLevels)*cellTotals[0]*cellTotals[1]+
                           (ZCompleteLevels)*(cellTotals[0]*(cellTotals[1]+1)+
                           (cellTotals[0]+1)*cellTotals[1]);
    
  }else{
    return localFaceNumber+(ZCompleteLevels+1)*cellTotals[0]*cellTotals[1]+
                           (ZCompleteLevels)*(cellTotals[0]*(cellTotals[1]+1)+
                           (cellTotals[0]+1)*cellTotals[1]);
  }
}

// ================================
// GET VORTEXES ASSOCIATED TO CELLS
// ================================
void MRITopology::getNeighborVortexes(int cellNumber,int dim,MRIIntVec& idx){
  // Loop through the edges
  MRIIntVec ElEdgeList;
  MRIDoubleVec currEdgeDirVector(3);
  int currFace = 0;
  int currEdge = 0;
  for(int loopA=0;loopA<cellFaces[cellNumber].size();loopA++){
    currFace = cellFaces[cellNumber][loopA];
    for(int loopB=0;loopB<faceEdges[currFace].size();loopB++){
      currEdge = faceEdges[currFace][loopB];
      MRIUtils::insertInList(currEdge,ElEdgeList);
    }
  }
  // Find the Edges Aligned with the Selected Dimension
  idx.clear();
  int currDir = 0;
  for(int loopA=0;loopA<ElEdgeList.size();loopA++){
    getEdgeDirection(ElEdgeList[loopA],currEdgeDirVector);
    if((fabs(currEdgeDirVector[1])<kMathZero)&&(fabs(currEdgeDirVector[2])<kMathZero)){
      currDir = 0;
    }else if((fabs(currEdgeDirVector[0])<kMathZero)&&(fabs(currEdgeDirVector[2])<kMathZero)){
      currDir = 1;
    }else if((fabs(currEdgeDirVector[0])<kMathZero)&&(fabs(currEdgeDirVector[1])<kMathZero)){
      currDir = 2;
    }else{
      throw MRIException("ERROR: Invalid Edge Direction in getNeighborVortexes.\n");
    }
    if(currDir == dim){
      idx.push_back(ElEdgeList[loopA]);
    }
  }
}

// ===================
// CREATE SAMPLE FLOWS
// ===================
void MRITopology::createFromTemplate(MRISamples sampleType,const MRIDoubleVec& params){

  // Store Parameter Values
  int sizeX = int(params[0]);
  int sizeY = int(params[1]);
  int sizeZ = int(params[2]);
  double distX = params[3];
  double distY = params[4];
  double distZ = params[5];
  double currTime = params[6];

  // Template Orientation
  int dir = 0;
  int direction = int(params[7]);
  if(direction == 0){
    dir = kdirX;
  }else if(direction == 1){
    dir = kdirY;
  }else if(direction == 2){
    dir = kdirZ;
  }else{
    throw MRIException("ERROR: Invalid template direction in CreateSampleCase.\n");
  }

  MRIIntVec currentCoords(kNumberOfDimensions);
  // Set Cells Totals
  cellTotals.resize(3);
  cellTotals[0] = sizeX;
  cellTotals[1] = sizeY;
  cellTotals[2] = sizeZ;

  cellLengths.resize(3);
  cellLengths[0].resize(sizeX);
  cellLengths[1].resize(sizeY);
  cellLengths[2].resize(sizeZ);  

  // Set Cell Lengths
  for(int loopA=0;loopA<kNumberOfDimensions;loopA++){
    for(int loopB=0;loopB<cellTotals[loopA];loopB++){
      switch(loopA){
        case 0:
          cellLengths[loopA][loopB] = distX;
          break;
        case 1:
          cellLengths[loopA][loopB] = distY;
          break;
        case 2:
          cellLengths[loopA][loopB] = distZ;
          break;
      }
    }
  }

  // Set Global Dimensions
  // Min
  domainSizeMin.resize(3);
  domainSizeMin[0] = 0.0;
  domainSizeMin[1] = 0.0;
  domainSizeMin[2] = 0.0;
  // Max
  domainSizeMax.resize(3);
  domainSizeMax[0] = (sizeX-1) * distX;
  domainSizeMax[1] = (sizeY-1) * distY;
  domainSizeMax[2] = (sizeZ-1) * distZ;
  
  // Set Total Cells
  totalCells = sizeX * sizeY * sizeZ;

  // Allocate CellLocations
  cellLocations.resize(totalCells);
  for(int loopA=0;loopA<totalCells;loopA++){
    cellLocations[loopA].resize(totalCells);
  }

  // Assign Coordinates
  for(int loopA=0;loopA<totalCells;loopA++){
    mapIndexToCoords(loopA,currentCoords);
    cellLocations[loopA][0] = currentCoords[0] * distX;
    cellLocations[loopA][1] = currentCoords[1] * distY;
    cellLocations[loopA][2] = currentCoords[2] * distZ;
  }
}
