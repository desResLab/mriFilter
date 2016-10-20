# include "mriTopology.h"

// ==========================
// GET LOCAL FACE CONNECTIONS
// ==========================
void getFaceConnections(int faceID, std::vector<int> cellConnections, std::vector<int> &faceIds){
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
  cellTotals[0] = opts.dimensions[0];
  cellTotals[1] = opts.dimensions[1];
  cellTotals[2] = opts.dimensions[2];
  // Assign total number of cells
  totalCells = cellTotals[0] * cellTotals[1] * cellTotals[2];
  cellLocations.resize(totalCells);
  for(int loopA=0;loopA<totalCells;loopA++){
    cellLocations[loopA].resize(totalCells);
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
  domainSizeMin[0] = opts.origin[0];
  domainSizeMin[1] = opts.origin[1];
  domainSizeMin[2] = opts.origin[2];
  // Max
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

// CHECK COMPATIBLE TOPOLOGY
bool MRITopology::isCompatibleTopology(MRITopology* topo){
  bool result = true;

  result = result && (totalCells != topo->totalCells);

  result = result && (cellTotals[0] != topo->cellTotals[0]);
  result = result && (cellTotals[1] != topo->cellTotals[1]);
  result = result && (cellTotals[2] != topo->cellTotals[2]);

  for(int loopA=0;loopA<cellTotals[0];loopA++){
    result = result && (fabs(cellLengths[0][loopA]-topo->cellLengths[0][loopA]) > kMathZero);
  }

  for(int loopA=0;loopA<cellTotals[1];loopA++){
    result = result && (fabs(cellLengths[1][loopA]-topo->cellLengths[1][loopA]) > kMathZero);
  }

  for(int loopA=0;loopA<cellTotals[2];loopA++){
    result = result && (fabs(cellLengths[2][loopA]-topo->cellLengths[2][loopA]) > kMathZero);
  }

  result = result && (fabs(domainSizeMin[0]-topo->domainSizeMin[0]) > kMathZero);
  result = result && (fabs(domainSizeMin[1]-topo->domainSizeMin[1]) > kMathZero);
  result = result && (fabs(domainSizeMin[2]-topo->domainSizeMin[2]) > kMathZero);

  result = result && (fabs(domainSizeMax[0]-topo->domainSizeMax[0]) > kMathZero);
  result = result && (fabs(domainSizeMax[1]-topo->domainSizeMax[1]) > kMathZero);
  result = result && (fabs(domainSizeMax[2]-topo->domainSizeMax[2]) > kMathZero);

  return result;
}