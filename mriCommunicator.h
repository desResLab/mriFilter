#ifndef MRICOMMUNICATOR_H
#define MRICOMMUNICATOR_H

class MRICommunicator{
public:
  // Main Communicator
  int mpiComm;
  // Data Members
  int currProc;
  int totProc;
  // Constructor
  MRICommunicator();
};

#endif // MRICOMMUNICATOR_H
