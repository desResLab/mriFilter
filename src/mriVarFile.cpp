/*  TApplicationModes = (amInitial,amFileOpen);
  VolDataTypes =  (vdAnatomy,vdVelocityX,vdVelocityY,vdVelocityZ);
  GlobalAdjTypes = (gaPlusX,gaMinusX,gaPlusY,gaMinusY,gaPlusZ,gaMinusZ);
  CellPointItem = Packed Record
    Position: Array3Single;
    Concentration: Single;
    Velocity: Array3Single;
    FilteredVel: Array3Single;
    PressGrad: Array3Single;
    RelPressure: Single;
  End;
  CellPointArray = Array Of CellPointItem;
  {Record with Global Data}
  GlobalDataRecord = Packed Record
    {Cells Totals}
    CellTotals: Array3Integer4;
    CellLength: Array3Single;
    {Domain Dimension}
    DomainSizeMin: Array3Single;
    DomainSizeMax: Array3Single;
    {Envelope Velocities}
    MaxVelModule: Double;
    {Velocities And Concentrations for all Measure Points}
    TotalCellPoints: Integer4;
    CellPoints: CellPointArray;
  End;
  {Option for MRI Reconstruction Algorithm}
  MRIOptionRecord = Packed Record
    Tolerance: Double;
    MaxIterations: Integer4;
  End;
  {Record For VOL Files}
  VolDataRecord = Packed Record
    GridX: Integer4;
    GridY: Integer4;
    GridZ: Integer4;
    SpaceX: Single;
    SpaceY: Single;
    SpaceSlice: Single;
    SpaceThick: Single;
    Voxels: Integer2Array;
  End;
  {Threshold Application Criteria}
  ThresholdTypes = (tcLessThen,tcGreaterThen,tcABSLessThen,tcABSGreaterThen);
  ThresholdCritRecord = Packed Record
    ThresholdQty: Integer2;
    ThresholdValue: Double;
    ThresholdType: ThresholdTypes;
  End;


Const
  VersionString = '0.0.4';
  NumberOfDimensions = 3;

Var
  {Global Data Array}
  GlobalData: GlobalDataRecord;
  {Check If File Is Loaded}
  GridLoaded: Boolean;
  {Check If Filter Is Applied}
  AreVelocityFiltered: Boolean;
  {PLT File Header}
  PltFileHeader: StringArray;
  {Stop Iterations}
  StopIterations: Boolean;
  {Pressure Gradient Calculated}
  IsPressureGradientAvailable: Boolean;
  {Relative Pressure Available}
  IsRelativePressureAvailable: Boolean;
  {Is Relative Pressure Smoothed}
  IsRelativePressureSmoothed: Boolean;

  {Reset Global Data Structure}
  Procedure ResetGlobalData(Var GlobalData: GlobalDataRecord);

  {Check If The Criterion Is Meet}
  Function MeetCriterion(ThresholdCriteria: ThresholdCritRecord;CurrentValue: Double): Boolean;

implementation

{Reset Global Data Record}
Procedure ResetGlobalData(Var GlobalData: GlobalDataRecord);
Begin
  With GlobalData Do
  Begin
    {Cells Totals}
    CellTotals[1]:=0;
    CellTotals[2]:=0;
    CellTotals[3]:=0;
    {Cells Lengths}
    CellLength[1]:=0.0;
    CellLength[2]:=0.0;
    CellLength[3]:=0.0;
    {Domain Dimension}
    {Min}
    DomainSizeMin[1]:=MaxDouble;
    DomainSizeMin[2]:=MaxDouble;
    DomainSizeMin[3]:=MaxDouble;
    {Max}
    DomainSizeMax[1]:=-MaxDouble;
    DomainSizeMax[2]:=-MaxDouble;
    DomainSizeMax[3]:=-MaxDouble;
    {Envelope Velocities}
    MaxVelModule:=0.0;
    {Velocities And Concentrations for all Measure Points}
    TotalCellPoints:=0;
    SetLength(CellPoints,0);
  End;
End;

{Check if Criterion is Meet}
Function MeetCriterion(ThresholdCriteria: ThresholdCritRecord;
                       CurrentValue: Double): Boolean;
Begin
  Case ThresholdCriteria.ThresholdType Of
    tcLessThen: Begin
                  Result:=(CurrentValue<ThresholdCriteria.ThresholdValue);
                End;
    tcGreaterThen: Begin
                     Result:=(CurrentValue>ThresholdCriteria.ThresholdValue);
                   End;
    tcABSLessThen: Begin
                     Result:=(ABS(CurrentValue)<ThresholdCriteria.ThresholdValue);
                   End;
    tcABSGreaterThen: Begin
                        Result:=(ABS(CurrentValue)>ThresholdCriteria.ThresholdValue);
                      End;
  End;
End;*/
