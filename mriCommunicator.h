#ifndef MRICOMMUNICATOR_H
#define MRICOMMUNICATOR_H

class mriCommunicator{
public:
  // Main Communicator
  int mpiComm;
  // Data Members
  int currProc;
  int totProc;
  // Constructor
  mriCommunicator();
};

#endif // MRICOMMUNICATOR_H
