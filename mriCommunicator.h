#ifndef MRICOMMUNICATOR_H
#define MRICOMMUNICATOR_H

# include "string"

# include "mpi.h"
# include "mriTypes.h"

using namespace std;

class MRICommunicator{
public:
  // Main Communicator
  MPI_Comm mpiComm;
  // Data Members
  int currProc;
  int totProc;
  // Constructor
  MRICommunicator();

  // SEND AND RECEIVE STD MATRICES AND VECTORS
  // Int Mat
  void sendStdIntMatrix(MRIIntMat matrix,int dest,int tag1, int tag2);
  void bcasStdIntMatrix(MRIIntMat matrix);
  void recvStdIntMatrix(MRIIntMat& matrix,int source,int tag1, int tag2);
  void passStdIntMatrix(MRIIntMat& matrix);
  // Double Mat
  void sendStdDoubleMatrix(MRIDoubleMat matrix,int dest,int tag1, int tag2);
  void bcasStdDoubleMatrix(MRIDoubleMat matrix);
  void recvStdDoubleMatrix(MRIDoubleMat& matrix,int source,int tag1, int tag2);
  void passStdDoubleMatrix(MRIDoubleMat& matrix);
  // Int Vector
  void sendStdIntVector(MRIIntVec vector,int dest,int tag);
  void bcasStdIntVector(MRIIntVec vector);
  void recvStdIntVector(MRIIntVec& vector,int source,int tag);
  void passStdIntVector(MRIIntVec vector);
  // Double Vector
  void sendStdDoubleVector(MRIDoubleVec vector,int dest,int tag);
  void bcasStdDoubleVector(MRIDoubleVec vector);
  void recvStdDoubleVector(MRIDoubleVec& vector,int source,int tag);
  void passStdDoubleVector(MRIDoubleVec& vector);
  // String
  void passString(string& msg);

};

#endif // MRICOMMUNICATOR_H
