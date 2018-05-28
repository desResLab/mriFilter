#ifndef MRIIO_H
#define MRIIO_H

# include "mriTypes.h"
# include "mriScan.h"
# include "mriTopology.h"
# include "mriExpansion.h"

class MRIScan;
class MRITopology;

void readPLTData(string PltFileName, MRITopology* topo, MRIScan* scan);

void readExpansionFile(string fileName,
                       MRIIntVec& tot,
                       MRIDoubleVec& lengthX,
                       MRIDoubleVec& lengthY,
                       MRIDoubleVec& lengthZ,
                       MRIDoubleVec& minlimits,
                       MRIDoubleVec& maxlimits,
                       MRIExpansion* exp);

void initVTKStructuredPointsOptions(vtkStructuredPointsOptionRecord &opts);

void assignVTKOptions(int lineNum, const MRIStringVec& tokens, vtkStructuredPointsOptionRecord &vtkOptions);

void assignPLTOptions(const MRIStringVec& tokens, pltOptionRecord& pltOptions);


#endif // MRIIO_H
