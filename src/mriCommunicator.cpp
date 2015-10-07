#include <stddef.h>

#include "mriCommunicator.h"
#include "mriUtils.h"
#include "mriTypes.h"
#include "mriSequence.h"
#include "mriStructuredScan.h"

MRICommunicator::MRICommunicator(){
}

MRICommunicator::~MRICommunicator(){

}

// ===========================
// PASS STD MATRIX OF INTEGERS
// ===========================
void MRICommunicator::passStdIntMatrix(MRIIntMat& matrix){
  int source = 0;
  int tag = 0;
  int size = 0;
  int totSize = 0;
  int unrolledSize = 0;
  int* sizeVec = NULL;
  int mpiError = 0;
  int* unrolledVec = NULL;
  int count = 0;
  MPI_Status status;
  if(currProc == 0){
    // UNROLL SIZES
    int matSize = matrix.size();
    int* sizeVec = new int[matSize];

    int mpiError = 0;
    for(int loopA=0;loopA<matrix.size();loopA++){
      sizeVec[loopA] = matrix[loopA].size();
      totSize +=  matrix[loopA].size();
    }
    // SEND SIZES
    for(int loopDest=1;loopDest<totProc;loopDest++){
      mpiError = MPI_Send(sizeVec,matSize, MPI_INT, loopDest, tag, mpiComm);
      MRIUtils::checkMpiError(mpiError);
    }    
  }else{
    // Probe Sizes Lengths
    mpiError = MPI_Probe(source,tag, mpiComm, &status);
    MRIUtils::checkMpiError(mpiError);
    mpiError = MPI_Get_count(&status,MPI_INT,&size);
    MRIUtils::checkMpiError(mpiError);
    // Receive Lengths Vector
    if(size > 0){
      // Receive Vector with Sizes
      sizeVec = new int[size];
      mpiError = MPI_Recv(sizeVec,size,MPI_INT,source,tag,mpiComm, &status);
      MRIUtils::checkMpiError(mpiError);
    }
  }

  if(currProc == 0){
    // SEND MATRIX VALUES
    int* intMatrixToSend = new int[totSize];
    int count = 0;
    for(int loopA=0;loopA<matrix.size();loopA++){
      for(int loopB=0;loopB<matrix[loopA].size();loopB++){
        intMatrixToSend[count] = matrix[loopA][loopB];
        // Update Count
        count++;
      }
    }
    // SEND VALUES
    for(int loopDest=1;loopDest<totProc;loopDest++){
      mpiError = MPI_Send(intMatrixToSend,totSize, MPI_INT, loopDest, tag, mpiComm);
      MRIUtils::checkMpiError(mpiError);
    }
    delete [] intMatrixToSend;
  }else{
    // Probe the Unrolled Vector Size
    mpiError = MPI_Probe(source,tag,mpiComm,&status);
    MRIUtils::checkMpiError(mpiError);
    mpiError = MPI_Get_count(&status,MPI_CHAR,&unrolledSize);
    MRIUtils::checkMpiError(mpiError);
    if(unrolledSize > 0){
      // Receive Vector with Sizes
      unrolledVec = new int[unrolledSize];
      mpiError = MPI_Recv(unrolledVec,unrolledSize,MPI_INT,source,tag,mpiComm, &status);
      MRIUtils::checkMpiError(mpiError);
    }
    // Allocate Matrix
    matrix.clear();
    matrix.resize(size);
    for(int loopA=0;loopA<size;loopA++){
      matrix[loopA].resize(sizeVec[loopA]);
    }
    // Fill Matrix
    count = 0;
    //printf("RECV MAT\n");
    for(int loopA=0;loopA<size;loopA++){
      for(int loopB=0;loopB<sizeVec[loopA];loopB++){
        matrix[loopA][loopB] = unrolledVec[count];
        //printf("%d ",matrix[loopA][loopB]);
        // Update Counter
        count++;
      }
      //printf("\n");
    }
    // Delete Temporary Vectors
    delete [] sizeVec;
    delete [] unrolledVec;
  }
}

// =======================
// BCAST ARRAY OF INTEGERS
// =======================
void MRICommunicator::passStdDoubleMatrix(MRIDoubleMat& matrix){
  int mpiError = 0;
  int tag = 0;
  int source = 0;
  int size = 0;
  int totSize = 0;
  int unrolledSize = 0;
  int count = 0;
  int* sizeVec = NULL;
  MPI_Status status;
  if(currProc == 0){
    // UNROLL SIZES
    int matSize = matrix.size();
    int* sizeVec = new int[matSize];
    for(int loopA=0;loopA<matSize;loopA++){
      sizeVec[loopA] = matrix[loopA].size();
      totSize +=  matrix[loopA].size();
    }
    // SEND SIZES
    for(int loopDest=1;loopDest<totProc;loopDest++){
      mpiError = MPI_Send(sizeVec,matSize, MPI_INT, loopDest, tag, mpiComm);
      MRIUtils::checkMpiError(mpiError);
    }
    if(size > 0){
      delete [] sizeVec;
    }
  }else{
    // Probe Sizes Lengths
    mpiError = MPI_Probe(source,tag,mpiComm,&status);
    MRIUtils::checkMpiError(mpiError);
    mpiError = MPI_Get_count(&status,MPI_INT,&size);
    MRIUtils::checkMpiError(mpiError);
    // Receive Lengths Vector
    if(size > 0){
      // Receive Vector with Sizes
      sizeVec = new int[size];
      mpiError = MPI_Recv(sizeVec,size,MPI_INT,source,tag,mpiComm, &status);
      MRIUtils::checkMpiError(mpiError);
    }
  }

  // SEND MATRIX VALUES
  if(currProc == 0){
    double* doubleMatrixToSend = new double[totSize];
    int count = 0;
    for(int loopA=0;loopA<matrix.size();loopA++){
      for(int loopB=0;loopB<matrix[loopA].size();loopB++){
        doubleMatrixToSend[count] = matrix[loopA][loopB];
        // Update Count
        count++;
      }
    }
    // SEND VALUES
    for(int loopDest=1;loopDest<totProc;loopDest++){
      mpiError = MPI_Send(doubleMatrixToSend,totSize, MPI_DOUBLE, loopDest, tag, mpiComm);
      MRIUtils::checkMpiError(mpiError);
    }
    delete [] doubleMatrixToSend;
  }else{
    double* unrolledVec;
    // Probe the Unrolled Vector Size
    mpiError = MPI_Probe(source,tag,mpiComm,&status);
    MRIUtils::checkMpiError(mpiError);
    mpiError = MPI_Get_count(&status,MPI_INT,&unrolledSize);
    MRIUtils::checkMpiError(mpiError);
    if(unrolledSize > 0){
      // Receive Vector with Sizes
      unrolledVec = new double[unrolledSize];
      mpiError = MPI_Recv(unrolledVec,unrolledSize,MPI_DOUBLE,source,tag,mpiComm, &status);
      MRIUtils::checkMpiError(mpiError);
    }
    // Allocate Matrix
    matrix.clear();
    matrix.resize(size);
    for(int loopA=0;loopA<size;loopA++){
      matrix[loopA].resize(sizeVec[loopA]);
    }
    // Fill Matrix
    count = 0;
    //printf("RECV MAT\n");
    for(int loopA=0;loopA<size;loopA++){
      for(int loopB=0;loopB<sizeVec[loopA];loopB++){
        matrix[loopA][loopB] = unrolledVec[count];
        //printf("%f ",matrix[loopA][loopB]);
        // Update Counter
        count++;
      }
      //printf("\n");
    }
    // Delete Temporary Vectors
    if(unrolledSize > 0){
      delete [] unrolledVec;
    }
  }
}

// ======================
// PASS ARRAY OF INTEGERS
// ======================
void MRICommunicator::passStdIntVector(MRIIntVec& vector){
  int mpiError = 0;
  int tag = 0;
  int source = 0;
  if(currProc == 0){
    int size = vector.size();
    int* buf = new int[size];
    for(int loopA=0;loopA<size;loopA++){
      buf[loopA] = vector[loopA];
    }
    for(int loopDest=1;loopDest<totProc;loopDest++){
      mpiError = MPI_Send(buf,size,MPI_INT,loopDest,tag,mpiComm);
    }
  }else{
    int vecSize = 0;
    int* buf = NULL;
    MPI_Status status;
    // Probe Vector Length
    mpiError = MPI_Probe(source,tag,mpiComm,&status);
    MRIUtils::checkMpiError(mpiError);
    mpiError = MPI_Get_count(&status,MPI_INT,&vecSize);
    MRIUtils::checkMpiError(mpiError);
    // Receive String
    if(vecSize > 0){
      buf = new int[vecSize];
      mpiError = MPI_Recv(buf,vecSize,MPI_INT,source,tag,mpiComm,&status);
      MRIUtils::checkMpiError(mpiError);
      vector.clear();
      for(int loopA=0;loopA<vecSize;loopA++){
        vector.push_back(buf[loopA]);
      }
      //for(int loopA=0;loopA<vector.size();loopA++){
      //  printf("Rec Vector %d %d\n",loopA,vector[loopA]);
      //}
      delete [] buf;
    }
  }
}

// =====================
// PASS ARRAY OF DOUBLES
// =====================
void MRICommunicator::passStdDoubleVector(MRIDoubleVec& vector){
  int source = 0;
  int tag = 0;
  int mpiError = 0;
  if(currProc == 0){
    int vecSize = vector.size();
    double* vecToSend = new double[vecSize];
    for(int loopA=0;loopA<vecSize;loopA++){
      vecToSend[loopA] = vector[loopA];
    }
    for(int loopDest=1;loopDest<totProc;loopDest++){
      mpiError = MPI_Send(vecToSend,vecSize,MPI_DOUBLE,loopDest,tag,mpiComm);
      MRIUtils::checkMpiError(mpiError);
    }
    delete [] vecToSend;
  }else{
    int mpiError = 0;
    int vecSize = 0;
    double* buf = NULL;
    MPI_Status status;
    // Probe Vector Length
    mpiError = MPI_Probe(source,tag,mpiComm,&status);
    MRIUtils::checkMpiError(mpiError);
    mpiError = MPI_Get_count(&status,MPI_DOUBLE,&vecSize);
    MRIUtils::checkMpiError(mpiError);
    // Receive String
    if(vecSize > 0){
      buf = new double[vecSize];
      mpiError = MPI_Recv(buf,vecSize,MPI_DOUBLE,source,tag,mpiComm,&status);
      MRIUtils::checkMpiError(mpiError);
      vector.clear();
      for(int loopA=0;loopA<vecSize;loopA++){
        vector.push_back(buf[loopA]);
        //printf("%f\n",buf[loopA]);
      }
      //printf("\n");
      delete [] buf;
    }
  }
}

// ===========
// PASS STRING
// ===========
void MRICommunicator::passString(string& msg){
  int source = 0;
  int tag = 0;
  if(currProc == 0){
    int strSize = msg.length() + 1;
    for(int loopDest=1;loopDest<totProc;loopDest++){
      int mpiError = MPI_Send(msg.c_str(),strSize,MPI_CHAR,loopDest,tag,mpiComm);
      MRIUtils::checkMpiError(mpiError);
    }
  }else{
    int mpiError = 0;
    int strSize = 0;
    char* buf = NULL;
    MPI_Status status;
    // Probe Vector Length
    mpiError = MPI_Probe(source,tag,mpiComm,&status);
    MRIUtils::checkMpiError(mpiError);
    mpiError = MPI_Get_count(&status,MPI_CHAR,&strSize);
    MRIUtils::checkMpiError(mpiError);
    // Receive String
    if(strSize > 0){
      buf = new char[strSize];
      mpiError = MPI_Recv(buf,strSize,MPI_CHAR,source,tag,mpiComm,&status);
      MRIUtils::checkMpiError(mpiError);
      msg = string(buf);
      delete [] buf;
    }
  }
}

// ==============
// PASS CELL DATA
// ==============
void MRICommunicator::passCellData(int& totalCellPoints,vector<MRICell>& cellPoints){
  MRIDoubleMat storeMat;
  MRIDoubleVec temp;  
  if(currProc == 0){
    temp.resize(7);
    // Create Matrices
    for(int loopA=0;loopA<cellPoints.size();loopA++){
      //printf("PASSING CELLS: %d\n",loopA);
      temp[0] = cellPoints[loopA].concentration;
      temp[1] = cellPoints[loopA].velocity[0];
      temp[2] = cellPoints[loopA].velocity[1];
      temp[3] = cellPoints[loopA].velocity[2];
      temp[4] = cellPoints[loopA].position[0];
      temp[5] = cellPoints[loopA].position[1];
      temp[6] = cellPoints[loopA].position[2];
      //printf("BEFORE STORE_PUSH\n");
      storeMat.push_back(temp);
      //printf("AFTER STORE_PUSH\n");
    }
  }

  // Pass Storage Matrix to other Processes
  passStdDoubleMatrix(storeMat);

  // Rebuild Cell Data in other processes
  if(currProc != 0){
    totalCellPoints = storeMat.size();
    // Resize
    cellPoints.resize(totalCellPoints);
    // Create Matrices
    for(int loopA=0;loopA<totalCellPoints;loopA++){
      cellPoints[loopA].concentration = storeMat[loopA][0];
      cellPoints[loopA].velocity[0] = storeMat[loopA][1];
      cellPoints[loopA].velocity[1] = storeMat[loopA][2];
      cellPoints[loopA].velocity[2] = storeMat[loopA][3];
      cellPoints[loopA].position[0] = storeMat[loopA][4];
      cellPoints[loopA].position[1] = storeMat[loopA][5];
      cellPoints[loopA].position[2] = storeMat[loopA][6];
    }
  }
}




