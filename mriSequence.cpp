#include <fstream>
#include <cmath>
#include <mriUtils.h>
#include <mriConstants.h>
#include "mriSequence.h"
#include "schMessages.h"

// Constructor
MRISequence::MRISequence(bool cyclic){
  // Initizalize
  totalScans = 0;
  // Set If Cyclic
  isCyclic = cyclic;
}

// Copy Constructor
MRISequence::MRISequence(MRISequence* copySequence){
  // Copy Scans
  totalScans = copySequence->totalScans;
  // Copy Cyclic Property
  isCyclic = copySequence->isCyclic;
  // Fill with Zero Scans
  for(int loopB=0;loopB<totalScans;loopB++){
    MRIScan* newScan = new MRIScan(*copySequence->GetScan(loopB));
    sequence.push_back(newScan);
  }
}

// Destructor
MRISequence::~MRISequence(){}

// Print the File List Log
void MRISequence::PrintSequenceFiles(std::string outFIleName){
  // Open Output File
	FILE* outFile;
	outFile = fopen(outFIleName.c_str(),"w");
  // Write List
  fprintf(outFile,"List of Files in Sequence\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    fprintf(outFile,"File %d: %s\n",loopA+1,fileNames[loopA].c_str());
  }
  // Close Output file
	fclose(outFile);				
}

// Add a Scan to the Sequence
void MRISequence::AddScan(MRIScan* scan){
  sequence.push_back(scan);
  totalScans++;
}

// Get a Scan Pointer From the Sequence
MRIScan* MRISequence::GetScan(int scanNumber){
  return sequence[scanNumber];
}

// EXPORT SEQUENCE TO TECPLOT FILE
void MRISequence::ExportToTECPLOT(std::string outfileName){
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("EXPORTING -------------------------------------\n"));
  for(int loopA=0;loopA<totalScans;loopA++){
      sequence[loopA]->ExportToTECPLOT(outfileName,(loopA == 0));
  }
}

// EXCTRACT SINGKLE POINT CURVE IN TIME
void MRISequence::ExtractSinglePointTimeCurve(int cellNumber, int exportQty, std::string fileName){
  // Open Output File
	FILE* outFile;
	outFile = fopen(fileName.c_str(),"w");
  // Var
  double currTime = 0.0;
  double currValue = 0.0;
  double currRad = 0.0;
  // Eval Centre Point
  double centrePoint[3] = {0.0};
  double localCoord[3] = {0.0};
  centrePoint[0] = 0.5 * (sequence[0]->domainSizeMax[0]+sequence[0]->domainSizeMin[0]);
  centrePoint[1] = 0.5 * (sequence[0]->domainSizeMax[1]+sequence[0]->domainSizeMin[1]);
  centrePoint[2] = 0.5 * (sequence[0]->domainSizeMax[2]+sequence[0]->domainSizeMin[2]);  
  // Write List
  fprintf(outFile,"List of Files in Sequence\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    // Get Time
    currTime = sequence[loopA]->scanTime;
    // Get Local Coordinates
    localCoord[0] = sequence[loopA]->cellPoints[cellNumber].position[0] - centrePoint[0];
    localCoord[1] = sequence[loopA]->cellPoints[cellNumber].position[1] - centrePoint[1];
    localCoord[2] = sequence[loopA]->cellPoints[cellNumber].position[2] - centrePoint[2];
    // Check Quantity
    switch(exportQty){
      case kQtyConcentration:
        currValue = sequence[loopA]->cellPoints[cellNumber].concentration;
        break;
      case kQtyVelocityX:
        currValue = sequence[loopA]->cellPoints[cellNumber].velocity[0];
        break;      
      case kQtyVelocityY:
        currValue = sequence[loopA]->cellPoints[cellNumber].velocity[1];
        break;      
      case kQtyVelocityZ:
        currValue = sequence[loopA]->cellPoints[cellNumber].velocity[2];
        break;      
      case kQtyPressGradientX:
        currValue = sequence[loopA]->cellPoints[cellNumber].pressGrad[0];
        break;      
      case kQtyPressGradientY:
        currValue = sequence[loopA]->cellPoints[cellNumber].pressGrad[1];
        break;      
      case kQtyPressGradientZ:
        currValue = sequence[loopA]->cellPoints[cellNumber].pressGrad[2];
        break;      
      case kQtyRelPressure:
        currValue = sequence[loopA]->cellPoints[cellNumber].relPressure;
        break;      
    }
    currRad = sqrt(localCoord[1]*localCoord[1]+
                   localCoord[2]*localCoord[2]);
    fprintf(outFile,"%e %e %e\n",currTime,currValue,currRad);
  }
  // Close Output file
	fclose(outFile);				  
}

// Compute Relative Pressure
void MRISequence::ComputeRelativePressure(bool doPressureSmoothing){
  WriteSchMessage(std::string("\n"));
  WriteSchMessage(std::string("REL PRESSURE COMPUTATION --------------------------\n"));
  // Loop Through Scans
  for(int loopA=0;loopA<totalScans;loopA++){
    // Write Message
    WriteSchMessage(std::string("Computing Relative Pressure - Step "+MRIUtils::IntToStr(loopA+1)+"/"+MRIUtils::IntToStr(totalScans)+"..."));
    // Get The Scan Back
    MRIScan* resultScan = GetScan(loopA);
    int startingCell = resultScan->EvalCentralCell();
    resultScan->EvalRelativePressure(startingCell,0.0);
    if (doPressureSmoothing){
      resultScan->PerformPressureIterations();
    }
    // Done
    WriteSchMessage(std::string("Relative Pressure Computed.\n"));
  }
}

// Read From VOL Sequence File
void MRISequence::ReadFromVolSequence(std::string seqfileName){
  // Var
  MRIScan* myScan;
  std::string buffer;
  
  // Open File
  std::ifstream seqFile;
  seqFile.open(seqfileName.c_str());

  // Loop Through Lines
  while (std::getline(seqFile,buffer))
  {
    if (buffer != ""){
      // Tokenize Line
      std::vector<std::string> ResultArray = MRIUtils::ExctractSubStringFromBufferMS(buffer);    
      // Create Scan Objects
      myScan = new MRIScan(atof(ResultArray[3].c_str()));
      // Files for a single scan in the same row
      myScan->ReadScanFromVOLFiles(ResultArray[0],ResultArray[1],ResultArray[2],ResultArray[3]);    
      // Add to Sequence
      AddScan(myScan);
    }
  }
  // Close File
  seqFile.close();
}
    
// Export to VOL File Set
void MRISequence::ExportToVOL(std::string outfileName){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ExportToVOL(outfileName+"_Step"+MRIUtils::IntToStr(loopA));
  }
}

// EXPORT TO SEQUENCE OF VTK FILES
void MRISequence::ExportToVTK(std::string outfileName){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ExportToVTK(outfileName+"_Step"+MRIUtils::IntToStr(loopA));
  }
}

// PHYSICS FILTERING FOR ALL SCANS
void MRISequence::ApplySMPFilter(MRIOptions* options, MRICommunicator* comm){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    // Perform Filter
    sequence[loopA]->applySMPFilter(options);
    // Update Velocities
    sequence[loopA]->UpdateVelocities();
  }
}

// APPLY THRESHOLDING TO ALL SCANS
void MRISequence::ApplyThresholding(MRIThresholdCriteria* thresholdCriteria){
  // Export All Data
  WriteSchMessage("\n");
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ApplyThresholding(thresholdCriteria);
  }
}

// Crop All Scans in Sequence
void MRISequence::Crop(double* limitBox){
  WriteSchMessage(std::string("Cropping Sequence..."));
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->Crop(limitBox);
  }
  WriteSchMessage(std::string("Done.\n"));
}

// Scale velocities for all Scans
void MRISequence::ScaleVelocities(double factor){
  WriteSchMessage(std::string("Scaling Velocities..."));
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ScaleVelocities(factor);
  }
  WriteSchMessage(std::string("Done.\n"));
}

// Scale Positions
void MRISequence::ScalePositions(double factor){
  WriteSchMessage(std::string("Scaling Positions..."));
  for(int loopA=0;loopA<totalScans;loopA++){
    sequence[loopA]->ScalePositions(factor);
  }
  WriteSchMessage(std::string("Done.\n"));
}









