# include <string>
# include <stdio.h>
# include <random>

# include "mriScan.h"
# include "mriSamplingOptions.h"
# include "schMessages.h"
# include "mriUtils.h"

using namespace std;

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
void MRIScan::sampleVelocities(MRISamplingOptions SamplingOptions, MRIIntVec& bins){
  random_device rand_dev;
  mt19937 generator(rand_dev());
  uniform_int_distribution<int>  distr(0, topology->totalCells - 1);

  // Set Parameters For PDF File
  bins.resize(SamplingOptions.numberOfBins);
  for(int loopA=0;loopA<bins.size();loopA++){
    bins[loopA] = 0.0;
  }
  MRIBoolVec visited(topology->totalCells - 1, false);

  // Set End Index
  int endIndex = 0;
  if(SamplingOptions.totalSamples>topology->totalCells){
    // Write Message
    WriteSchMessage("WARNING: Number of Samples exceeding the total number of cells. Reducing Samples.\n");
    endIndex = topology->totalCells;  
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
    currentCell = distr(generator);
    
    if(!visited[currentCell]){

      // Find cell centroid
      currentXCoord = topology->cellLocations[currentCell][0];
      currentYCoord = topology->cellLocations[currentCell][1];
      currentZCoord = topology->cellLocations[currentCell][2];

      // Eval Cell Velocity Module
      currentModule = sqrt((cells[currentCell].velocity[0]*cells[currentCell].velocity[0])+
                           (cells[currentCell].velocity[1]*cells[currentCell].velocity[1])+
                           (cells[currentCell].velocity[2]*cells[currentCell].velocity[2]));

      // Increment Counter
      if((SamplingOptions.useBox)&&(MRIUtils::isPointInsideBox(currentXCoord,currentYCoord,currentZCoord,SamplingOptions.limitBox))){
        // Increment The Number of Samples
        count++;
        // Assign to the Right Bin
        bins[(int)(currentModule/binInterval) + 1]++;
      }

    }
    visited[currentCell] = true;
  }
}
