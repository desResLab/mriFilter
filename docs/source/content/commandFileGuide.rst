Command File Guide for smpFilter
================================

Running Mode
^^^^^^^^^^^^

Only the normal running mode is supported at the moment ::

  RUNMODE: NORMAL

Input/Output
^^^^^^^^^^^^

Input file format
"""""""""""""""""

This token allows the user to specify various input file formats. The available formats are:

1. **VTK**. VTK File format. A structured_point dataset need to be specified in Ascii format.
2. **TECPLOT**. ASCII Tecplot file format. 
3. **TEMPLATE**. 
4. **EXPANSION**. 

Example input: ::

  INPUTTYPE: VTK

An example VTK file is ::

  # vtk DataFile Version 2.0
  Example dataset
  ASCII
  DATASET STRUCTURED_POINTS
  DIMENSIONS 10 10 10
  ORIGIN 0.0 0.0 0.0
  SPACING 1.000000 1.000000 1.000000
  POINT_DATA 1000
  VECTORS velocity float
  vx1 vy1 vz1
  ...
  ...
  ...

Input file name
"""""""""""""""

Example input: ::
  
  INPUTFILE: inputFile.vtk

Using pre-defined templates
"""""""""""""""""""""""""""

Using the TEMPLATETYPE token it is possible to create a pre-defined flow configuration 

1. **ZEROVELOCITY**. 
2. **CONSTANT**. 
3. **POISEILLE**. 
4. **STAGNATION**. 
5. **CYLINDRICALVOLTEX**. 
6. **SPHERICALVORTEX**. 
7. **TOROIDALVORTEX**. 
8. **TRANSIENT**. 
9. **CONSTANTWITHSTEP**. 
10. **TAYLORVORTEX**. 

Example input: ::

Using pre-defined templates
"""""""""""""""""""""""""""

Example input: ::

TEMPLATEPARAMS

Output file type
""""""""""""""""

The output file format can be specified using the following formats:

1. **VTK**. 
2. **TECPLOT**. ASCII Tecplot file format. 

Example input: ::

  OUTPUTTYPE: VTK

Output file name
""""""""""""""""

The name of the output file is specified through the following token:

  OUTPUTFILE: outputFile.vtk

Gaussian/Median Filter
^^^^^^^^^^^^^^^^^^^^^^

It is possible to filter the velocity data using either a Gaussian or Median filter approach as follows:

1. **MEDIAN**. 
2. **GAUSSIAN**. 

Example input: ::

  APPLYMEDIANFILTER: MEDIAN

SMP Filter
^^^^^^^^^^

SMP Filter activation
"""""""""""""""""""""

Use the USESMPFILTER token to activate/deactivate the SMP filter. 

Two possible inputs can be specified:

1. **TRUE**. The filter is active.
2. **FALSE**. The filter is inactive.

Example input: :: 

  USESMPFILTER: TRUE


Use Constant flow waveforms
"""""""""""""""""""""""""""

The USECONSTANTPATTERNS will include three constant waveform at each SMP iteration. This helps to speed up the convergence especially for flows characterized by a strong average component.

Example input: ::

  USECONSTANTPATTERNS: TRUE

Iteration tolerance and number of iterations
""""""""""""""""""""""""""""""""""""""""""""

The SMP convergence tolerance can be specified using the SMPITERATIONTOLERANCE token. This is the tolerance for the relative change in the 2-norm of the residual between successive iterations.

Example input: ::

  SMPITERATIONTOLERANCE: 1.0e-4

Example input: ::

  SMPMAXITERATIONS: 

Adding Noise
^^^^^^^^^^^^

In some situations you may want to add Gaussianly distributed, component independent noise, to an input velocity field. To do so, the ADDNOISE token allows to enter the intensity of the noise as a percent of the maximum velocity module.

Example input: ::

  ADDNOISE: 10.0

This will use 10\% of the maximum velocity module in the field as the standard deviation of the Gaussian noise intensity.

Physical Constants
^^^^^^^^^^^^^^^^^^

Two tokens, DENSITY and VISCOSITY can be used to specify physical constants. These constant are mainly use for pressure estimation through the Pressure Poisson Equation method or other methods.  

Example input: ::

  DENSITY: 1060.0
  VISCOSITY: 4.0e-3

Wall thresholds
^^^^^^^^^^^^^^^

When applying the boundary condition filter, the subset of the computational domain occupied by solid walls is specified through the THRESHOLDQTY, THRESHOLDTYPE and THRESHOLDVALUE tokens.

Threshold Quantity
""""""""""""""""""

The following options specify the quantity used to define the threshold:

1. **POSX**. X coordinate. 
2. **POSY**. Y coordinate. 
3. **POSZ**. Z coordinate. 
4. **CONCENTRATION**. The concentration 
5. **VELX**. Velocity in the X direction. 
6. **VELY**. Velocity in the Y direction. 
7. **VELZ**. Velocity in the Z direction. 
8. **VELMOD**. Velocity module. 
9. **NONE**. No thresholding.

Threshold Type
""""""""""""""

Various threshold criteria can be used:

1. **LT**. Less than.
2. **GT**. Greater than.
3. **ABSLT**. Less than in absolute value.
4. **ABSGT**. Greater than in absolute value.

Threshold Value
"""""""""""""""
The THRESHOLDVALUE token is used to specify the numerical value of the threshold. 

Example input: ::

  THRESHOLDQTY: CONCENTRATION
  THRESHOLDTYPE: GT
  THRESHOLDVALUE: 0.5

This means that all the cells with concentration greater than 0.5 will be considered as walls. 

Export to PPE Poisson Solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The EXPORTTOPOISSON and POISSONFILE tell the application to export a finite element input file for successive solution with a PPE solver. 

Example input: ::

  EXPORTTOPOISSON: TRUE
  POISSONFILE: poissonInputFile.txt

Vortex Criteria
^^^^^^^^^^^^^^^

Traditional Vortex Criteria
"""""""""""""""""""""""""""

The EVALVORTEXCRITERIA is responsible to add three criteria to the Q, L2 and Delta result file (this works only for the VTK output file format).

Example input: ::

  EVALVORTEXCRITERIA: TRUE

An additional vortex criteria based on the vertex frame representation can also be plotted using the EVALSMPVORTEXCRITERIA token. 

Example input: ::

  EVALSMPVORTEXCRITERIA: TRUE  

Other Options
^^^^^^^^^^^^^

Save Initial Velocities
"""""""""""""""""""""""

In some cases, the user may want to save the original velocity field before any manipulation is performed. This is accomplished through the SAVEINITIALVELOCITIES token.

Example input: ::

  SAVEINITIALVELOCITIES: TRUE


Save Expansion Coefficients
"""""""""""""""""""""""""""

The coefficient representation computed using the SMP filter can be saved and a flow field can be restored by an expansion file. See also the INPUTTYPE token above for instructions on how to read an expansion file. 

Example input: ::

  SAVEEXPANSIONCOEFFS: TRUE
