#include <string>
#include <stdio.h>
#include "mriScan.h"
#include "mriSamplingOptions.h"
#include "schMessages.h"
#include "mriUtils.h"

//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_int_distribution.hpp>

// Print Bins to File
void PrintBinsToFile(std::string fileName, int totalSamples, double binInterval, int totalBins, int* bins){
	// Open Output File
	FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");

  // Print Bins and Corresponding Interval
  for(int loopA=0;loopA<totalBins;loopA++){
    fprintf(outFile,"%e %e\n",(loopA-1)*binInterval+0.5*binInterval,((double)bins[loopA]/(double)totalSamples));
  }
	// Close Output file
	fclose(outFile);				  
}

// Sample Velocities
void MRIScan::SampleVelocities(MRISamplingOptions SamplingOptions){
  // Initialize Generator
  //boost::random::mt19937 gen;
  
  // Set Parameters For PDF File
  int* bins = new int[SamplingOptions.numberOfBins];
  for(int loopA=0;loopA<SamplingOptions.numberOfBins;loopA++) bins[loopA] = 0;

  // Set End Index
  int endIndex = 0;
  if(SamplingOptions.totalSamples>totalCellPoints){
    // Write Message
    WriteSchMessage("Number of Samples exceeding the total number of cells. Reducing Samples.\n");
    endIndex = totalCellPoints;  
  }else{
    endIndex = SamplingOptions.totalSamples;  
  }
  double binInterval = ((double)maxVelModule/(double)SamplingOptions.numberOfBins);
  int count = 0;
  int currentCell = 0;
  double currentXCoord,currentYCoord,currentZCoord;
  double currentModule;
  while(count<endIndex){
    // Draw a Random Number Between 1 and GlobalData.TotalCellPoints
    WriteSchMessage("ERROR: ADD LIBRARY TO GENERATE RANDOM NUMBERS!!!\n");
    //currentCell = MRIUtils::GenerateUniformIntegers(0, totalCellPoints-1);

    // Find cell centroid
    currentXCoord = cellPoints[currentCell].position[0];
    currentYCoord = cellPoints[currentCell].position[1];
    currentZCoord = cellPoints[currentCell].position[2];

    // Eval Cell Velocity Module
    currentModule = sqrt((cellPoints[currentCell].velocity[0]*cellPoints[currentCell].velocity[0])+
                         (cellPoints[currentCell].velocity[1]*cellPoints[currentCell].velocity[1])+
                         (cellPoints[currentCell].velocity[2]*cellPoints[currentCell].velocity[2]));

    // Increment Counter
    if((SamplingOptions.useBox)&&(MRIUtils::IsPointInsideBox(currentXCoord,currentYCoord,currentZCoord,SamplingOptions.limitBox))){
      // Increment The Number of Samples
      count++;
      // Assign to the Right Bin
      bins[(int)(currentModule/binInterval)+1]++;
    }
  }
}
