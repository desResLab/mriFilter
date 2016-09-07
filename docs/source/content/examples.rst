Examples
========

Template Poiseille Flow
^^^^^^^^^^^^^^^^^^^^^^^

This example simply creates a 20x20x25 grid with spacing of 0.002,0.002 and 0.004 along the x,y and z axis respectively. A parabolic velocity profile is assigned with velocities along the z axis (direction parameter set to 2). Gaussian noise is also added to the velocity field with zero mean and standard deviation equal to 5% of the maximum velocity module in the grid. 

Command file: ::

  # Running in Normal Mode
  RUNMODE: NORMAL

  # Assign density and dynamic viscosity
  DENSITY: 1060.0
  VISCOSITY: 4.0e-3

  # Add 5% Gaussian Noise
  ADDNOISE: 5.0

  # Input TEMPLATE
  INPUTTYPE: TEMPLATE
  TEMPLATETYPE: POISEILLE
  TEMPLATEPARAMS: 20.0,20.0,25.0,0.002,0.002,0.004,0.0,2.0

  # Output File and Format
  OUTPUTTYPE: VTK
  OUTPUTFILE: Poiseille.vtk

  # Threshold - Select what belongs to the wall
  THRESHOLDQTY: CONCENTRATION
  THRESHOLDTYPE: LT
  THRESHOLDVALUE: 0.5

Template Stagnation Flow
^^^^^^^^^^^^^^^^^^^^^^^^

Command file: ::

  # Running in Normal Mode
  RUNMODE: NORMAL

  # Input TEMPLATE
  INPUTTYPE: TEMPLATE
  TEMPLATETYPE: STAGNATION
  TEMPLATEPARAMS: 20.0,20.0,20.0,0.01,0.01,0.01,0.0,2.0

  # Output File and Format
  OUTPUTTYPE: VTK
  OUTPUTFILE: Stagnation_NoNoise_NoFilter.vtk

  # Add 10% Gaussian Noise
  ADDNOISE: 10.0

  # Use SMP Filters
  USESMPFILTER: TRUE
  USEBCFILTER: FALSE
  USECONSTANTPATTERNS: TRUE

  # Set Filter Tolrances
  SMPITERATIONTOLERANCE: 1.0e-03
  SMPMAXITERATIONS: 100

  # Threshold - Select what to remove
  THRESHOLDQTY: CONCENTRATION
  THRESHOLDTYPE: LT
  THRESHOLDVALUE: 0.5

