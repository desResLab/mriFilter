#ifndef MRICOMMUNICATOR_H
#define MRICOMMUNICATOR_H

# include "mpi.h"

class MRICommunicator{
public:
  // Main Communicator
  MPI_Comm mpiComm;
  // Data Members
  int currProc;
  int totProc;
  // Constructor
  MRICommunicator();
};

#endif // MRICOMMUNICATOR_H
