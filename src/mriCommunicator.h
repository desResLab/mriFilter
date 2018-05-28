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

class MRICommunicator{
public:
  // Main Communicator
  MPI_Comm mpiComm;
  // Data Members
  int currProc;
  int totProc;
  // Constructor
  MRICommunicator();
  ~MRICommunicator();

  // PASS CELL DATA
  void passCellData(int& totalCellPoints,vector<MRICell>& cellPoints);

  // SEND AND RECEIVE STD MATRICES AND VECTORS
  // Int Mat
  void passStdIntMatrix(MRIIntMat& matrix);
  // Double Mat
  void passStdDoubleMatrix(MRIDoubleMat& matrix);
  // Int Vector
  void passStdIntVector(MRIIntVec& vector);
  // Double Vector
  void passStdDoubleVector(MRIDoubleVec& vector);
  // String
  void passString(string& msg);

};

#endif // MRICOMMUNICATOR_H
