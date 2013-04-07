#ifndef MRIOPTIONS_H
#define MRIOPTIONS_H

class MRIOptions
{
  public:
	  // Data Members
  	double tolerance;
    int maxIterations;		
	  // Constructor and Distructors
	  MRIOptions(double tol, int maxIt);
	  ~MRIOptions(){};
};
	
#endif // MRIOPTIONS_H
