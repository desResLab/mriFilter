#ifndef MRIIO_H
#define MRIIO_H

# include "mriTypes.h"
# include "mriScan.h"
# include "mriTopology.h"
# include "mriExpansion.h"

class mriScan;
class mriTopology;

void readPLTData(string PltFileName, mriTopology* topo, mriScan* scan);

void readExpansionFile(string fileName,
                       mriIntVec& tot,
                       mriDoubleVec& lengthX,
                       mriDoubleVec& lengthY,
                       mriDoubleVec& lengthZ,
                       mriDoubleVec& minlimits,
                       mriDoubleVec& maxlimits,
                       mriExpansion* exp);

void initVTKStructuredPointsOptions(vtkStructuredPointsOptionRecord &opts);

void assignVTKOptions(int lineNum, const mriStringVec& tokens, vtkStructuredPointsOptionRecord &vtkOptions);

void assignPLTOptions(const mriStringVec& tokens, pltOptionRecord& pltOptions);


#endif // MRIIO_H
