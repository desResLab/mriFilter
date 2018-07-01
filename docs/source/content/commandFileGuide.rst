Command File Guide for smpFilter
================================

Running Mode
^^^^^^^^^^^^

The normal running mode is the typical execution mode (other modes are experimental and will be integrated in future releases) ::

  RUNMODE: NORMAL

Input/Output
^^^^^^^^^^^^

Input file format
"""""""""""""""""

This token allows the user to specify various input file formats. The available formats are:

1. **VTK**, VTK File format. A structured_point dataset need to be specified in ASCII format.
2. **TECPLOT**, ASCII Tecplot file format. 
3. **TEMPLATE**, Creates a template model using the parameters specified in **TEMPLATEPARAMS**. 
4. **EXPANSION**, Reads the model from a file containing the vortex frame expansion coefficients. 

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
  SCALARS concentration double 1
  LOOKUP_TABLE default
  ...
  ...
  VECTORS velocity float
  vx1 vy1 vz1
  ...
  ...
  ...

Note that the dataset type is **STRUCTURED_POINTS**, the file format is ASCII and the there must be a *SCALAR* quantity called **concentration** (i.e., the magnitude of the acquired field or a mask) and a *VECTOR* quantity called **velocity**.

An example of TECPLOT file is: ::

  TITLE     = "Volume 50"
  VARIABLES = "X"
  "Y"
  "Z"
  "Magnitude<4000"
  "Vx"
  "Vy"
  "Vz"
  ZONE T="SubZone"
  STRANDID=0, SOLUTIONTIME=0
  I=2, J=61, K=115, ZONETYPE=Ordered
  DATAPACKING=POINT
  DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )
  -9.350372553E-001 -1.990501761E+000 5.128377533E+001 1.000000000E+000 0.000000000E+000 0.000000000E+000 0.000000000E+000
  -5.867734528E+001 -1.990501761E+000 5.128377533E+001 1.000000000E+000 0.000000000E+000 0.000000000E+000 0.000000000E+000
  -9.350372553E-001 -1.090502143E+000 5.128377533E+001 1.000000000E+000 0.000000000E+000 0.000000000E+000 0.000000000E+000
  -5.867734528E+001 -1.090502143E+000 5.128377533E+001 1.000000000E+000 0.000000000E+000 0.000000000E+000 0.000000000E+000
  -9.350372553E-001 -1.905025542E-001 5.128377533E+001 1.000000000E+000 0.000000000E+000 0.000000000E+000 0.000000000E+000
  -5.867734528E+001 -1.905025542E-001 5.128377533E+001 1.000000000E+000 0.000000000E+000 0.000000000E+000 0.000000000E+000
  -9.350372553E-001 7.094970942E-001 5.128377533E+001 1.000000000E+000 0.000000000E+000 0.000000000E+000 0.000000000E+000
  ...
  ...

Note that the datapacking type is **POINT** and the file format is ASCII. Note also that there are **7** columns where the first three contain the **x,y and z** cell-center coordinates, a column with the **concentration** and the last three columns containing the **x,y and z** components of the velocity. 

Input file name
"""""""""""""""

Example input: ::
  
  INPUTFILE: inputFile.vtk

Using pre-defined templates
"""""""""""""""""""""""""""

Using the TEMPLATETYPE token it is possible to create a pre-defined flow configuration. The following hardcoded templates are available: 

1. **ZEROVELOCITY**, a velocity field with uniformly zero velocities. 
2. **CONSTANT**, a constant flow along one of the coordinate axis, x,y or z. 
3. **POISEUILLE**, a parabolic velocity profile in a cilinder.
4. **STAGNATION**, stagnation point flow. 
5. **CYLINDRICALVOLTEX**. 
6. **SPHERICALVORTEX**. 
7. **TOROIDALVORTEX**. 
8. **TRANSIENT**. 
9. **CONSTANTWITHSTEP**. 
10. **ROTATINGVORTEX** Flow in a rotating cylinder.

Specification of the selected template is completed by specifing a few parameters using **TEMPLATEPARAMS**. The first six parameters specify the geometry of the cartesian acquisition grid with uniform spacing. These are followed by a time and a direction parameter.

Example input: ::

  TEMPLATEPARAMS: sizeX,sizeY,sizeZ,distX,distY,distZ,time,direction,auxParams1,auxParams2,auxParams3... 

These parameters are:

1. **sizeX**, the number of measurement cells in the x direction.
2. **sizeY**, the number of measurement cells in the y direction.
3. **sizeZ**, the number of measurement cells in the z direction.
4. **distX**, cell spacing in the x direction.
5. **distY**, cell spacing in the y direction.
6. **distZ**, cell spacing in the z direction.
7. **time**, time parameter (used only for the **TRANSIENT** template)
8. **direction**, orientation parameter.
9. **auxParamters**, template specific parameters.

The auxiliary parameters for the various templates are:

1. **ZEROVELOCITY**, no additional parameters. 
2. **CONSTANT**, no additional parameters.
3. **POISEUILLE**, two addtional parameters:

   - Radius of the cylindrical domain.
   - Peak velocity.

4. **STAGNATION**, no additional parameters.
5. **CYLINDRICALVOLTEX**, no additional parameters. 
6. **SPHERICALVORTEX**, no additional parameters. 
7. **TOROIDALVORTEX**, no additional parameters. 
8. **TRANSIENT**, no additional parameters. 
9. **CONSTANTWITHSTEP**, no additional parameters. 
10. **ROTATINGVORTEX**, no additional parameters.

Output file type
""""""""""""""""

The output file format can be specified using the following formats:

1. **VTK**, VTK Legacy file format, either STRUCTURED_POINTS or STRUCTURED_GRID (form non uniformly spaced grids). 
2. **TECPLOT**, Tecplot ASCII file format. 

Example input: ::

  OUTPUTTYPE: VTK

Output file name
""""""""""""""""

The name of the output file is specified through the following token:

  OUTPUTFILE: outputFile.vtk

Mean/Median/Gaussian Filter
^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to filter the velocity data using either a Mean, Median or Gaussian filter approach as follows:

1. **MEAN**, applies a mean filter to the neighbors of every cell (direct neighbors for filter of order one). 
2. **MEDIAN**, applies a median filter. 
3. **GAUSSIAN**, applied a Gaussian filter. 

Example input: ::

  APPLYSMOOTHINGFILTER: 2,MEDIAN,1

The **first integer parameter** defines the number of filter iterations. If it is greater the one the filter is applied multiple times.
The **third integer parameter** is the neighboring order that defines how many cells are included when computing the filtered velocities. 

Scaling Geometry and Velocities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes when the units used for the geometry of the acquisition grid and the velocities are not consistent it may be usefull to scale indipendently these two quantities. This is achived by the **SCALEPOSITION** and **SCALEVELOCITY** tokens.

Example input: ::

  # Scale from mm to m 
  SCALEPOSITION: 1.0e-3

  # Scale from mm/s to m/s 
  SCALEVELOCITY: 1.0e-3

Note that when scaling the position, the measurement grid is **translated in space** and the minimum coordinates are set equal to the **origin** of the axis system (0.0,0.0,0.0).

Solenoidal Filter
^^^^^^^^^^^^^^^^^

Use the USESMPFILTER token to apply a solenoidal filter. Four possible options can be specified:

1. Use **boundary filter**. This activates a boundary condition filter after the main filter. This provides improvement in the situation where solid walls are present in the domain and one wants to both have a solenoidal velocity field plus satisfy the boundary conditions at the interface between fluid and walls. 
2. Use **constant flow waveforms**. This includes three constant waveform at each SMP iteration. This helps to speed up the convergence especially for flows characterized by a strong average component.
3. The SMP **convergence tolerance**. This is the tolerance for the relative change in the 2-norm of the residual between successive iterations.
4. The **number of iterations**. Maximum number of iterations for the iterative solver.

Example input: :: 

  USESMPFILTER: TRUE,TRUE,1.0e-3,2000

Adding Noise
^^^^^^^^^^^^

In some situations you may want to add Gaussianly distributed, component independent noise, to an input velocity field. To do so, the ADDNOISE token allows to enter the intensity of the noise as a percent of the maximum velocity module together with the seed for the random number generator

Example input: ::

  ADDNOISE: 10.0,1234

This will use 10\% of the maximum velocity module in the field as the standard deviation of the Gaussian noise intensity and a random seed equal to **1234**.

Physical Constants
^^^^^^^^^^^^^^^^^^

Two tokens, **DENSITY** and **VISCOSITY** can be used to specify physical constants. These constant are mainly use for pressure estimation through the Pressure Poisson Equation (PPE) or alternative approaches.  

Example input: ::

  # Apply density and viscosity of blood in SI units
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

The **EXPORTTOPOISSON** token tells the application to export a finite element input file for successive solution with a PPE solver. 

Example input: ::

  EXPORTTOPOISSON: poissonInputFile.txt

Components of the pressure gradient to evaluate
"""""""""""""""""""""""""""""""""""""""""""""""

The **PRESSUREGRADIENTCOMPONENTS** token tell the program which components of the pressure gradient need to be included in the computation of the relative pressures by solving th PPE. The following components can be activated/deactivated:

1. Acceleration (or transient) term ("Y" to activate and "N" to deactivate).
2. Advection (or convection) term ("Y" to activate and "N" to deactivate).
3. Diffusion term ("Y" to activate and "N" to deactivate).
4. Reynolds Stress term ("Y" to activate and "N" to deactivate). The eddy viscosity approximation is activated by default. Therefore by specifying "N" the turbulent viscosity is neglected. 

By default, if the **PRESSUREGRADIENTCOMPONENTS** is not specified all but the Reynolds Stress terms are included. 

The example below shows how to include only the acceleration and advection components of the pressure gradient and to exclude the diffusion and Reynolds Stress term: ::

  PRESSUREGRADIENTCOMPONENTS: Y,Y,N,N

Computing the Turbulent Viscosity
"""""""""""""""""""""""""""""""""

The token **USETURBVISCOSITYFROMFILE** allows to read the value of the turbulent viscosity from file. The file contains only one column where the value of the turbulent viscosity is specified at each row. The order should reflect the order of the grid points in the velocity measurements (i.e., the first row should contain the turbulent viscosity for the grid point one, the second row for point two, and so on) ::

  muT_Point_1
  muT_Point_2
  muT_Point_3
  muT_Point_4
  ...

The token **TURBVISCOSITYFILE** specifies the location of the file above. 

In case the turbulent viscosity is not read from file, the token **SMAGORINSKYCONSTANT** allows the user to specify the Smagorinsky constant in the associated subgrid scale model. If this constant is not specified a default value of 0.15 is used. 

Export to Laplace wall distancing solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The **EXPORTTODISTANCE** token tell the application to export a finite element input file for wall distance evaluation.

Example input: ::

  EXPORTTODISTANCE: distanceInputFile.txt

Vortex Criteria
^^^^^^^^^^^^^^^

Traditional Vortex Criteria
"""""""""""""""""""""""""""

The EVALVORTEXCRITERIA is responsible to add three criteria to the Q, L2 and Delta result file.

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

The coefficient representation computed using the SMP filter can be saved and a flow field can be restored by an expansion file. See also the **INPUTTYPE** token above for instructions on how to read an expansion file. 

Example input: ::

  SAVEEXPANSIONCOEFFS: TRUE
