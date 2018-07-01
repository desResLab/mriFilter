#ifndef MRICOMMUNICATOR_H
#define MRICOMMUNICATOR_H

# include <stddef.h>
# include <string>

# include "mriCell.h"
# include "mriUtils.h"
# include "mriTypes.h"
# include "mriScan.h"
# include "mriSequence.h"

# include "mpi.h"

using namespace std;

class mriCommunicator{
public:
  // Main Communicator
  MPI_Comm mpiComm;
  // Data Members
  int currProc;
  int totProc;
  // Constructor
  mriCommunicator();
  virtual ~mriCommunicator();

  // PASS CELL DATA
  void passCellData(int& totalCellPoints,vector<mriCell>& cellPoints);

  // SEND AND RECEIVE STD MATRICES AND VECTORS
  // Int Mat
  void passStdIntMatrix(mriIntMat& matrix);
  // Double Mat
  void passStdDoubleMatrix(mriDoubleMat& matrix);
  // Int Vector
  void passStdIntVector(mriIntVec& vector);
  // Double Vector
  void passStdDoubleVector(mriDoubleVec& vector);
  // String
  void passString(string& msg);

};

#endif // MRICOMMUNICATOR_H
