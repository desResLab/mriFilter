#include <iostream>
#include "mriStructuredScan.h"
#include "mpi.h"
#include "parmetis.h"
#include "schMessages.h"
#include "mriCommunicator.h"

using namespace std;

// ===========================================
// SOLVE POISSON EQUATION WITH FINITE ELEMENTS
// ===========================================
int MRIStructuredScan::SolvePoissonEquation(mriCommunicator* comm){

  WriteSchMessage(std::string("Sono Dentro Poisson\n"));

  // Declare Vectors to store the partition
  int* elmdist = nullptr;
  int *eptr = nullptr;
  int *eind = nullptr;
  int *elmwgt = nullptr;
  // Unweighted Graph
  int wgtflag = 0;
  // C-Style Numbering Flag
  // Array Indexes starting from 0
  int numflag = 0;
  // Number of Constraints (1 constraint)
  int ncon = 1;
  // Common Nodes ??
  int ncommonnodes = 0;
  // Number of Partitions (inependent on the number of processors)
  int nparts = 0;
  // Total weight in every single partition
  double *tpwgts = nullptr;
  tpwgts = new double(ncon*nparts);
  // Maximum degree of imbalance (size ncon)
  double *ubvec = nullptr;
  ubvec = new double(ncon);
  // Option vector
  int options[3] = {0};
  // Parmetis Results
  int *edgecut;
  int *part;

  // Set the number of partition to the
  // number of processors
  nparts = comm->totProc;

  // Set equal weights to all parts
  for(int loopA=0;loopA<nparts;loopA++){
    tpwgts[loopA] = 1.0/(double) nparts;
  }

  // Set the maximum imbalance
  for(int loopA=0;loopA<ncon;loopA++){
    ubvec[loopA] = 1.05;
  }

  if(comm->currProc == 0){
    // Build the node coordinate vector


    // Build the vectors needed for partioning
    buildMetisConnectivities(elmdist,eptr,eind);
  }

  // Partition Problem
  int resFlag = ParMETIS_V3_PartMeshKway(
    elmdist,eptr,eind,elmwgt,
    &wgtflag,&numflag,&ncon,&ncommonnodes,&nparts,
    tpwgts,ubvec,options,
    edgecut,part,
    &comm->mpiComm);

  // Check Execution
  if(resFlag == METIS_OK){
    WriteSchMessage(std::string("Parmetis excution successful.\n"));
  }else{
    WriteSchMessage(std::string("Parmetis excution terminated.\n"));
  }

  // Assemble Partition Matrix and Load Vector

  // Solve with CG

  // Return OK
  return 0;
}


int test(){

  int result;
// Needed by parmetis

  idx_t *vtxdist=NULL;
  idx_t *xadj=NULL;
  idx_t *adjncy=NULL;
  idx_t *vwgt=NULL, *adjwgt=NULL;
  idx_t wgtflag=0;
  idx_t numflag=0;
  idx_t ncon=1;
  idx_t nparts=3;
  real_t *tpwgts=NULL, ubvec;
  idx_t options[4], edgecut;
  idx_t part[5];

// For AdaptiveRepart
  real_t itr;
  idx_t *vsize=NULL;

// Start Comm
  // MPI VARIABLES
  int MPI_PROC_ID;
  int MPI_PROC_TOTAL_NUM;
  MPI_Comm comm;
  double TOTAL_TIME_ELAPSED;
  int ierr;
  MPI::Init ();
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &MPI_PROC_TOTAL_NUM );
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &MPI_PROC_ID );
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  if ( MPI_PROC_ID == 0 ){
    cout << " Parmetis example from LiberLocus."<< '\n';
  }
  cout << "I am Proc " << MPI_PROC_ID  << '\n';
// Common for every processor
    vtxdist = new idx_t[4];
// For AdaptiveRepart
    itr = 1000.0;

    vtxdist[0] = 0;
    vtxdist[1] = 5;
    vtxdist[2] = 10;
    vtxdist[3] = 15;

    ubvec = 1.05;

    options[0] = 0;
    options[1] = 0;
    options[2] = 0;
    options[3] = 0;

    part[0] = MPI_PROC_ID;
    part[1] = MPI_PROC_ID;
    part[2] = MPI_PROC_ID;
    part[3] = MPI_PROC_ID;
    part[4] = MPI_PROC_ID;

    tpwgts = new real_t[3];
//    tpwgts[0] = static_cast<float>(ncon) / static_cast<float>(nparts);
    tpwgts[0] = 1.0/3.0;
    tpwgts[1] = 1.0/3.0;
    tpwgts[2] = 1.0/3.0;
// Dependent on each processor
    if ( MPI_PROC_ID == 0 ){
    xadj = new idx_t[6];
    adjncy = new idx_t[13];

    xadj[0] = 0;
    xadj[1] = 2;
    xadj[2] = 5;
    xadj[3] = 8;
    xadj[4] = 11;
    xadj[5] = 13;

    adjncy[0] = 1;
    adjncy[1] = 5;
    adjncy[2] = 0;
    adjncy[3] = 2;
    adjncy[4] = 6;
    adjncy[5] = 1;
    adjncy[6] = 3;
    adjncy[7] = 7;
    adjncy[8] = 2;
    adjncy[9] = 4;
    adjncy[10] = 8;
    adjncy[11] = 3;
    adjncy[12] = 9;

    }
    else if ( MPI_PROC_ID == 1 ){
    xadj = new idx_t[6];
    adjncy = new idx_t[18];

    xadj[0] = 0;
    xadj[1] = 3;
    xadj[2] = 7;
    xadj[3] = 11;
    xadj[4] = 15;
    xadj[5] = 18;

    adjncy[0] = 0;
    adjncy[1] = 6;
    adjncy[2] = 10;
    adjncy[3] = 1;
    adjncy[4] = 5;
    adjncy[5] = 7;
    adjncy[6] = 11;
    adjncy[7] = 2;
    adjncy[8] = 6;
    adjncy[9] = 8;
    adjncy[10] = 12;
    adjncy[11] = 3;
    adjncy[12] = 7;
    adjncy[13] = 9;
    adjncy[14] = 13;
    adjncy[15] = 4;
    adjncy[16] = 8;
    adjncy[17] = 14;


    }
    else if ( MPI_PROC_ID == 2 ){
    xadj = new idx_t[6];
    adjncy = new idx_t[13];

    xadj[0] = 0;
    xadj[1] = 2;
    xadj[2] = 5;
    xadj[3] = 8;
    xadj[4] = 11;
    xadj[5] = 13;

    adjncy[0] = 5;
    adjncy[1] = 11;
    adjncy[2] = 6;
    adjncy[3] = 10;
    adjncy[4] = 12;
    adjncy[5] = 7;
    adjncy[6] = 11;
    adjncy[7] = 13;
    adjncy[8] = 8;
    adjncy[9] = 12;
    adjncy[10] = 14;
    adjncy[11] = 9;
    adjncy[12] = 13;

    }
    if ( MPI_PROC_ID == 0 ){
      cout << "parmetis initialized." << '\n';
    }

//  result = ParMETIS_V3_PartKway( vtxdist, xadj, adjncy, vwgt, adjwgt,
//                                 &wgtflag, &numflag, &ncon,
//                                 &nparts, tpwgts, &ubvec, options,
//                                 &edgecut, part, &comm );
  result = ParMETIS_V3_AdaptiveRepart( vtxdist, xadj, adjncy, vwgt, vsize,
                                 adjwgt, &wgtflag, &numflag, &ncon,
                                 &nparts, tpwgts, &ubvec, &itr, options,
                                 &edgecut, part, &comm );

    if ( MPI_PROC_ID == 0 ){
      cout << "parmetis finalized." << '\n';
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if ( MPI_PROC_ID == 0 ){
      cout << MPI_PROC_ID << " edgecut " << edgecut << '\n';
      for(int i=0; i<5; i++){
        cout << MPI_PROC_ID << " " << part[i] << '\n';
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if ( MPI_PROC_ID == 1 ){
      cout << MPI_PROC_ID << " edgecut " << edgecut << '\n';
      for(int i=0; i<5; i++){
        cout << MPI_PROC_ID << " " << part[i] << '\n';
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if ( MPI_PROC_ID == 2 ){
      cout << MPI_PROC_ID << " edgecut " << edgecut << '\n';
      for(int i=0; i<5; i++){
        cout << MPI_PROC_ID << " " << part[i] << '\n';
      }
    }
    delete vtxdist;
    delete xadj;
    delete adjncy;
    delete tpwgts;

// Finish Comm

  ierr = MPI_Finalize();

  return 0;

}
