#include <stddef.h>

#include "mriCommunicator.h"
#include "mriUtils.h"
#include "mriTypes.h"
#include "mriSequence.h"
#include "mriStructuredScan.h"

MRICommunicator::MRICommunicator(){
}

// ===========================
// SEND STD MATRIX OF INTEGERS
// ===========================
void MRICommunicator::sendStdIntMatrix(MRIIntMat matrix,int dest,int tag1, int tag2){
  // UNROLL SIZES
  int* sizeVec = new int(matrix.size());
  int totSize = 0;
  int mpiError = 0;
  for(int loopA=0;loopA<matrix.size();loopA++){
    sizeVec[loopA] = matrix[loopA].size();
    totSize +=  matrix[loopA].size();
  }
  // SEND SIZES
  mpiError = MPI_Send(&sizeVec[0],matrix.size(), MPI_INT, dest, tag1, mpiComm);
  MRIUtils::checkMpiError(mpiError);
  printf("Sending Int Matrix Size Array of size: %d, tag %d\n",matrix.size(),tag1);

  delete [] sizeVec;
  // SEND MATRIX VALUES
  int* intMatrixToSend = new int(totSize);
  int count = 0;
  for(int loopA=0;loopA<matrix.size();loopA++){
    for(int loopB=0;loopB<matrix[loopA].size();loopB++){
      intMatrixToSend[count] = matrix[loopA][loopB];
      // Update Count
      count++;
    }
  }
  // SEND VALUES
  mpiError = MPI_Send(&intMatrixToSend[0],totSize, MPI_INT, dest, tag2, mpiComm);
  MRIUtils::checkMpiError(mpiError);
  delete [] intMatrixToSend;
  printf("Sending Int Matrix of size: %d, tag %d\n",totSize,tag2);

}

// ===========================
// BCAST STD MATRIX OF INTEGERS
// ===========================
void MRICommunicator::bcasStdIntMatrix(MRIIntMat matrix){
  // UNROLL SIZES
  int matSize = matrix.size();
  int* sizeVec = new int[matSize];
  int totSize = 0;
  int mpiError = 0;
  for(int loopA=0;loopA<matrix.size();loopA++){
    sizeVec[loopA] = matrix[loopA].size();
    totSize +=  matrix[loopA].size();
  }
  // SEND SIZES
  mpiError = MPI_Bcast(sizeVec,matSize, MPI_INT,0,mpiComm);
  MRIUtils::checkMpiError(mpiError);
  printf("Sending Int Matrix Size Array of size: %d.\n",matrix.size());

  delete [] sizeVec;
  // SEND MATRIX VALUES
  int* intMatrixToSend = new int(totSize);
  int count = 0;
  for(int loopA=0;loopA<matrix.size();loopA++){
    for(int loopB=0;loopB<matrix[loopA].size();loopB++){
      intMatrixToSend[count] = matrix[loopA][loopB];
      // Update Count
      count++;
    }
  }
  // SEND VALUES
  mpiError = MPI_Bcast(intMatrixToSend,totSize, MPI_INT,0,mpiComm);
  MRIUtils::checkMpiError(mpiError);
  delete [] intMatrixToSend;
  printf("Sending Int Matrix of size: %d\n",totSize);
}


// ==============================
// RECEIVE STD MATRIX OF INTEGERS
// ==============================
void MRICommunicator::recvStdIntMatrix(MRIIntMat& matrix,int source,int tag1, int tag2){
  MPI_Status status;
  int size = 0;
  int unrolledSize = 0;
  int* sizeVec = NULL;
  int* unrolledVec = NULL;
  int mpiError = 0;
  int count = 0;
  // Probe Sizes Lengths
  mpiError = MPI_Probe(source,tag1, mpiComm, &status);
  MRIUtils::checkMpiError(mpiError);
  mpiError = MPI_Get_count(&status,MPI_INT,&size);
  MRIUtils::checkMpiError(mpiError);
  // Receive Lengths Vector
  if(size > 0){
    // Receive Vector with Sizes
    sizeVec = new int[size];
    mpiError = MPI_Recv(sizeVec,size,MPI_INT,source,tag1,mpiComm, &status);
    MRIUtils::checkMpiError(mpiError);
    printf("Receiving Int Matrix sizeArray of size: %d, tag %d\n",size,tag1);

    // Probe the Unrolled Vector Size
    mpiError = MPI_Probe(source,tag2,mpiComm,&status);
    MRIUtils::checkMpiError(mpiError);
    mpiError = MPI_Get_count(&status,MPI_CHAR,&unrolledSize);
    MRIUtils::checkMpiError(mpiError);
    if(unrolledSize > 0){
      // Receive Vector with Sizes
      unrolledVec = new int[unrolledSize];
      mpiError = MPI_Recv(unrolledVec,unrolledSize,MPI_INT,source,tag2,mpiComm, &status);
      MRIUtils::checkMpiError(mpiError);
      printf("Receiving Int Matrix of unrolled size: %d, tag %d\n",unrolledSize,tag2);
    }
    // Allocate Matrix
    matrix.clear();
    matrix.resize(size);
    for(int loopA=0;loopA<size;loopA++){
      matrix[loopA].resize(sizeVec[loopA]);
    }
    // Fill Matrix
    count = 0;
    for(int loopA=0;loopA<size;loopA++){
      for(int loopB=0;loopB<sizeVec[loopA];loopB++){
        matrix[loopA][loopB] = unrolledVec[count];
        // Update Counter
        count++;
      }
    }
    // Delete Temporary Vectors
    delete [] sizeVec;
    delete [] unrolledVec;
  }
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
    delete [] sizeVec;
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
      mpiError = MPI_Send(&intMatrixToSend[0],totSize, MPI_INT, loopDest, tag, mpiComm);
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
    printf("RECV MAT\n");
    for(int loopA=0;loopA<size;loopA++){
      for(int loopB=0;loopB<sizeVec[loopA];loopB++){
        matrix[loopA][loopB] = unrolledVec[count];
        printf("%d ",matrix[loopA][loopB]);
        // Update Counter
        count++;
      }
      printf("\n");
    }
    // Delete Temporary Vectors
    delete [] unrolledVec;
  }
}

// ==========================
// SEND STD MATRIX OF DOUBLES
// ==========================
void MRICommunicator::sendStdDoubleMatrix(MRIDoubleMat matrix,int dest,int tag1, int tag2){
  int mpiError = 0;
  // UNROLL SIZES
  int* sizeVec = new int(matrix.size());
  int totSize = 0;
  for(int loopA=0;loopA<matrix.size();loopA++){
    sizeVec[loopA] = matrix[loopA].size();
    totSize +=  matrix[loopA].size();
  }

  // SEND SIZES
  mpiError = MPI_Send(sizeVec,matrix.size(), MPI_INT, dest, tag1, mpiComm);
  MRIUtils::checkMpiError(mpiError);
  printf("Send Double Matrix of sizeArray %d, tag %d\n",matrix.size(),tag1);
  delete [] sizeVec;

  // SEND MATRIX VALUES
  double* doubleMatrixToSend = new double(totSize);
  int count = 0;
  for(int loopA=0;loopA<matrix.size();loopA++){
    for(int loopB=0;loopB<matrix[loopA].size();loopB++){
      doubleMatrixToSend[count] = matrix[loopA][loopB];
      // Update Count
      count++;
    }
  }
  // SEND VALUES
  mpiError = MPI_Send(doubleMatrixToSend,totSize, MPI_DOUBLE, dest, tag2, mpiComm);
  MRIUtils::checkMpiError(mpiError);
  delete [] doubleMatrixToSend;
  printf("Send Double Matrix of size %d, tag %d\n",totSize,tag2);
}

// ===========================
// BCAST STD MATRIX OF DOUBLES
// ===========================
void MRICommunicator::bcasStdDoubleMatrix(MRIDoubleMat matrix){
  int mpiError = 0;
  // UNROLL SIZES
  int matSize = matrix.size();
  int* sizeVec = new int[matSize];
  int totSize = 0;
  for(int loopA=0;loopA<matrix.size();loopA++){
    sizeVec[loopA] = matrix[loopA].size();
    totSize +=  matrix[loopA].size();
  }

  // SEND SIZES
  mpiError = MPI_Bcast(sizeVec,matSize, MPI_INT,0,mpiComm);
  MRIUtils::checkMpiError(mpiError);
  printf("Send Double Matrix of sizeArray %d\n",matrix.size());
  delete [] sizeVec;

  // SEND MATRIX VALUES
  double* doubleMatrixToSend = new double(totSize);
  int count = 0;
  for(int loopA=0;loopA<matrix.size();loopA++){
    for(int loopB=0;loopB<matrix[loopA].size();loopB++){
      doubleMatrixToSend[count] = matrix[loopA][loopB];
      // Update Count
      count++;
    }
  }
  // SEND VALUES
  mpiError = MPI_Bcast(doubleMatrixToSend,totSize, MPI_INT,0,mpiComm);
  MRIUtils::checkMpiError(mpiError);
  delete [] doubleMatrixToSend;
  printf("Send Double Matrix of size %d.\n",totSize);
}

// =============================
// RECEIVE STD MATRIX OF DOUBLES
// =============================
void MRICommunicator::recvStdDoubleMatrix(MRIDoubleMat& matrix,int source,int tag1, int tag2){
  int mpiError = 0;
  MPI_Status status;
  int size = 0;
  int unrolledSize = 0;
  int* sizeVec = NULL;
  double* unrolledVec = NULL;
  int count = 0;
  // Probe Sizes Lengths
  mpiError = MPI_Probe(source,tag1,mpiComm,&status);
  MRIUtils::checkMpiError(mpiError);
  mpiError = MPI_Get_count(&status,MPI_INT,&size);
  MRIUtils::checkMpiError(mpiError);
  printf("Received Double Matrix sizeInfo of size %d, tag %d\n",size,tag1);
  // Receive Lengths Vector
  if(size > 0){
    // Receive Vector with Sizes
    sizeVec = new int(size);
    mpiError = MPI_Recv(sizeVec,size,MPI_INT,source,tag1,mpiComm, &status);
    MRIUtils::checkMpiError(mpiError);
    printf("Received Double Matrix sizeArray of size %d, tag %d\n",size,tag1);

    // Probe the Unrolled Vector Size
    mpiError = MPI_Probe(source,tag2,mpiComm,&status);
    MRIUtils::checkMpiError(mpiError);
    mpiError = MPI_Get_count(&status,MPI_INT,&unrolledSize);
    MRIUtils::checkMpiError(mpiError);
    printf("RECV 1: Unrolled Size: %d\n",unrolledSize);
    if(unrolledSize > 0){
      // Receive Vector with Sizes
      unrolledVec = new double[unrolledSize];
      mpiError = MPI_Recv(unrolledVec,unrolledSize,MPI_DOUBLE,source,tag2,mpiComm, &status);
      MRIUtils::checkMpiError(mpiError);
      printf("Received Double Matrix of size %d, tag %d\n",unrolledSize,tag2);
    }
    // Allocate Matrix
    matrix.clear();
    matrix.resize(size);
    for(int loopA=0;loopA<size;loopA++){
      matrix[loopA].resize(sizeVec[loopA]);
    }
    // Fill Matrix
    count = 0;
    for(int loopA=0;loopA<size;loopA++){
      for(int loopB=0;loopB<sizeVec[loopA];loopB++){
        matrix[loopA][loopB] = unrolledVec[count];
        // Update Counter
        count++;
      }
    }
    // Delete Temporary Vectors
    if(size > 0){
      delete [] sizeVec;
    }
    if(unrolledSize > 0){
      delete [] unrolledVec;
    }
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
    printf("RECV 1: Unrolled Size: %d\n",unrolledSize);
    if(unrolledSize > 0){
      // Receive Vector with Sizes
      unrolledVec = new double[unrolledSize];
      mpiError = MPI_Recv(unrolledVec,unrolledSize,MPI_DOUBLE,source,tag,mpiComm, &status);
      MRIUtils::checkMpiError(mpiError);
      printf("Received Double Matrix of size %d, tag %d\n",unrolledSize,tag);
    }
    // Allocate Matrix
    matrix.clear();
    matrix.resize(size);
    for(int loopA=0;loopA<size;loopA++){
      matrix[loopA].resize(sizeVec[loopA]);
    }
    // Fill Matrix
    count = 0;
    printf("RECV MAT\n");
    for(int loopA=0;loopA<size;loopA++){
      for(int loopB=0;loopB<sizeVec[loopA];loopB++){
        matrix[loopA][loopB] = unrolledVec[count];
        printf("%f ",matrix[loopA][loopB]);
        // Update Counter
        count++;
      }
      printf("\n");
    }
    // Delete Temporary Vectors
    if(unrolledSize > 0){
      delete [] unrolledVec;
    }
  }
}

// =======================
// BCAST ARRAY OF INTEGERS
// =======================
void MRICommunicator::sendStdIntVector(MRIIntVec vector,int dest,int tag){
  int size = vector.size();
  int* buf = new int(size);
  for(int loopA=0;loopA<size;loopA++){
    buf[loopA] = vector[loopA];
  }
  printf("Sending Int Vector of Size: %d\n",size);
  int mpiError = MPI_Send(buf,size,MPI_INT,dest,tag,mpiComm);
  MRIUtils::checkMpiError(mpiError);
  delete [] buf;
}

// ======================
// SEND ARRAY OF INTEGERS
// ======================
void MRICommunicator::bcasStdIntVector(MRIIntVec vector){
  int size = vector.size();
  int mpiError = 0;
  int* buf = new int[size];
  for(int loopA=0;loopA<size;loopA++){
    buf[loopA] = vector[loopA];
  }
  printf("Sending Int Vector of Size: %d\n",size);
  mpiError = MPI_Bcast(buf,size,MPI_INT,0,mpiComm);
  MRIUtils::checkMpiError(mpiError);
  delete [] buf;
}

// ======================
// PASS ARRAY OF INTEGERS
// ======================
void MRICommunicator::passStdIntVector(MRIIntVec vector){
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
    printf("Received Integer Vector Size: %d\n",vecSize);
    // Receive String
    if(vecSize > 0){
      buf = new int[vecSize];
      mpiError = MPI_Recv(buf,vecSize,MPI_INT,source,tag,mpiComm,&status);
      MRIUtils::checkMpiError(mpiError);
      vector.clear();
      for(int loopA=0;loopA<vecSize;loopA++){
        vector.push_back(buf[loopA]);
      }
      for(int loopA=0;loopA<vector.size();loopA++){
        printf("Rec Vector %d %d\n",loopA,vector[loopA]);
      }
      delete [] buf;
    }
  }
}

// =========================
// RECEIVE ARRAY OF INTEGERS
// =========================
void MRICommunicator::recvStdIntVector(MRIIntVec& vector,int source,int tag){
  int mpiError = 0;
  int vecSize = 0;
  int* buf = NULL;
  MPI_Status status;
  // Probe Vector Length
  mpiError = MPI_Probe(source,tag,mpiComm,&status);
  MRIUtils::checkMpiError(mpiError);
  mpiError = MPI_Get_count(&status,MPI_INT,&vecSize);
  MRIUtils::checkMpiError(mpiError);
  printf("Received Integer Vector Size: %d\n",vecSize);
  // Receive String
  if(vecSize > 0){
    buf = new int[vecSize];
    mpiError = MPI_Recv(buf,vecSize,MPI_INT,source,tag,mpiComm,&status);
    MRIUtils::checkMpiError(mpiError);
    vector.clear();
    for(int loopA=0;loopA<vecSize;loopA++){
      vector.push_back(buf[loopA]);
    }
    delete [] buf;
  }
}

// =====================
// SEND ARRAY OF DOUBLES
// =====================
void MRICommunicator::sendStdDoubleVector(MRIDoubleVec vector,int dest,int tag){
  int mpiError = MPI_Send(&vector[0],vector.size(),MPI_DOUBLE,dest,tag,mpiComm);
  MRIUtils::checkMpiError(mpiError);
  // Message
  printf("Double Vector Sent From %d, size %d with tag %d\n",currProc,vector.size(),tag);
}

// ======================
// BCAST ARRAY OF DOUBLES
// ======================
void MRICommunicator::bcasStdDoubleVector(MRIDoubleVec vector){
  int size = vector.size();
  int* vecToSend = new int[size];
  int mpiError = MPI_Bcast(vecToSend,size,MPI_INT,0,mpiComm);
  MRIUtils::checkMpiError(mpiError);
  // Message
  printf("Double Vector Sent From %d, size %d\n",currProc,vector.size());
}

// ========================
// RECEIVE ARRAY OF DOUBLES
// ========================
void MRICommunicator::recvStdDoubleVector(MRIDoubleVec& vector,int source,int tag){
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
    }
    delete [] buf;
  }
  // Message
  printf("Double Vector Received From %d, size %d with tag %d\n",currProc,vecSize,tag);

}

// =====================
// PASS ARRAY OF DOUBLES
// =====================
void MRICommunicator::passStdDoubleVector(MRIDoubleVec& vector){
  int source = 0;
  int tag = 0;
  if(currProc == 0){
    int vecSize = vector.size();
    double* vecToSend = new double[vecSize];
    for(int loopA=0;loopA<vecSize;loopA++){
      vecToSend[loopA] = vector[loopA];
    }
    for(int loopDest=1;loopDest<totProc;loopDest++){
      int mpiError = MPI_Send(vecToSend,vecSize,MPI_DOUBLE,loopDest,tag,mpiComm);
      MRIUtils::checkMpiError(mpiError);
    }
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
      printf("RECV DOUBLE VEC\n");
      for(int loopA=0;loopA<vecSize;loopA++){
        vector.push_back(buf[loopA]);
        printf("%f\n",buf[loopA]);
      }
      printf("\n");
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
    int strSize = msg.length();
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



