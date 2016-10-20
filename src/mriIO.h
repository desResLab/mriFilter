#ifndef MRIIO_H
#define MRIIO_H

# include "mriTypes.h"
# include "mriScan.h"
# include "mriTopology.h"

class MRIScan;

void readPLTData(string PltFileName, MRITopology* topo, MRIScan* scan);

#endif // MRIIO_H
