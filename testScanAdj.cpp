#ifndef TESTRANDOM_CPP
#define TESTRANDOM_CPP

#include "mriScan.h"

void MRIScan::TestScanAdjacency(std::string fileName){
  FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");
  // Write mesh Dimensions
  fprintf(outFile,"Mesh Dims, X: %d, Y: %d, Z: %d\n",cellTotals[0],cellTotals[1],cellTotals[2]);
  // TEST Face ADJACENCY
  int currFace = 0;
  int faceLocation = 0;
  for(int loopA=0;loopA<totalCellPoints;loopA++){
    fprintf(outFile,"Face Adjacency for Cell: %d\n",loopA);
    for(int loopB=0;loopB<k3DNeighbors;loopB++){
      switch(loopB){
        case 0:
          faceLocation = kfacePlusX;
          break;
        case 1:
          faceLocation = kfaceMinusX;
          break;
        case 2:
          faceLocation = kfacePlusY;
          break;
        case 3:
          faceLocation = kfaceMinusY;
          break;
        case 4:
          faceLocation = kfacePlusZ;
          break;
        case 5:
          faceLocation = kfaceMinusZ;
          break;
      }        
      currFace = GetAdjacentFace(loopA,faceLocation);  
      fprintf(outFile,"%d\n",currFace);
    }
  }
	// Close Output file
	fclose(outFile);
}

#endif //TESTRANDOM_CPP
