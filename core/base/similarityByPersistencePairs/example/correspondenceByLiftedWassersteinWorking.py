# state file generated using paraview version 5.8.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.8.1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [773, 530]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [42.0, 42.0, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [34.052854108481576, 54.03090381692208, 229.03886399896956]
renderView1.CameraFocalPoint = [42.0, 42.0, 0.0]
renderView1.CameraViewUp = [0.1929299752044011, 0.9801895823046299, -0.044792939277352437]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 59.39696961966999

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [627, 530]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.CenterOfRotation = [31.5, 31.5, 31.5]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [-133.19493765754288, -91.38311706750326, 261.78483859066694]
renderView2.CameraFocalPoint = [31.499999999999957, 31.49999999999996, 31.500000000000004]
renderView2.CameraViewUp = [-0.8093661788602131, 0.496316346497576, -0.31400075273240297]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 54.559600438419636
renderView2.BackEnd = 'OSPRay raycaster'

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.447740)
layout1.AssignView(1, renderView2)
layout1.AssignView(2, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'TTK CinemaReader'
tTKCinemaReader1 = TTKCinemaReader(DatabasePath='F:/Dev/Thesis/tracking/PN_Translate.cdb')

# create a new 'TTK CinemaQuery'
tTKCinemaQuery1 = TTKCinemaQuery(InputTable=tTKCinemaReader1)
tTKCinemaQuery1.SQLStatement = """SELECT * FROM InputTable0
ORDER BY Time
LIMIT 5"""

# create a new 'TTK CinemaProductReader'
tTKCinemaProductReader1 = TTKCinemaProductReader(Input=tTKCinemaQuery1)

# create a new 'TTK ForEach'
tTKForEach1 = TTKForEach(Input=tTKCinemaProductReader1)
tTKForEach1.IterationMode = 'Block'
tTKForEach1.InputArray = ['POINTS', 'Field']
tTKForEach1.OutputType = 'vtkImageData'
tTKForEach1.ImageExtent = [0, 63, 0, 63, 0, 63]

# create a new 'TTK PersistenceDiagram'
tTKPersistenceDiagram1 = TTKPersistenceDiagram(Input=tTKForEach1)
tTKPersistenceDiagram1.ScalarField = ['POINTS', 'Field']
tTKPersistenceDiagram1.InputOffsetField = ['POINTS', '1']
tTKPersistenceDiagram1.EmbedinDomain = 1

# create a new 'Threshold'
threshold3 = Threshold(Input=tTKPersistenceDiagram1)
threshold3.Scalars = ['CELLS', 'Persistence']
threshold3.ThresholdRange = [0.1907828253507614, 0.9084896445274353]

# create a new 'TTK CorrespondenceByPersistencePairs'
tTKCorrespondenceByPersistencePairs1 = TTKCorrespondenceByPersistencePairs(Diagrams=threshold3)
tTKCorrespondenceByPersistencePairs1.Extremumweight = 0.91
tTKCorrespondenceByPersistencePairs1.Labels = ['POINTS', 'CriticalType']

# create a new 'TTK BlockAggregator'
tTKBlockAggregator4 = TTKBlockAggregator(Input=tTKCorrespondenceByPersistencePairs1)

# create a new 'TTK Extract'
tTKExtract1 = TTKExtract(Input=tTKPersistenceDiagram1)
tTKExtract1.ExtractionMode = 'Geometry'
tTKExtract1.Expression = '0.01'
tTKExtract1.ImageExtent = [2147483647, -2147483647, 2147483647, -2147483647, 2147483647, -2147483647]
tTKExtract1.ValidationMode = '> '
tTKExtract1.InputArray = ['CELLS', 'Persistence']

# create a new 'TTK TopologicalSimplification'
tTKTopologicalSimplification1 = TTKTopologicalSimplification(Domain=tTKForEach1,
    Constraints=tTKExtract1)
tTKTopologicalSimplification1.ScalarField = ['POINTS', 'Field']
tTKTopologicalSimplification1.InputOffsetField = ['POINTS', '2']
tTKTopologicalSimplification1.VertexIdentifierField = ['POINTS', '1']

# create a new 'Calculator'
calculator1 = Calculator(Input=tTKTopologicalSimplification1)
calculator1.ResultArrayName = 'Mask'
calculator1.Function = 'Field>0.7'
calculator1.ResultArrayType = 'Signed Char'

# create a new 'TTK ConnectedComponents'
tTKConnectedComponents1 = TTKConnectedComponents(Segmentation=calculator1)
tTKConnectedComponents1.FeatureMask = ['POINTS', 'Mask']

# find source
tTKConnectedComponents1_1 = FindSource('TTKConnectedComponents1')

# create a new 'TTK BlockAggregator'
tTKBlockAggregator2 = TTKBlockAggregator(Input=OutputPort(tTKConnectedComponents1_1,1))

# create a new 'TTK MorphologicalOperators'
tTKMorphologicalOperators1 = TTKMorphologicalOperators(Input=tTKConnectedComponents1)
tTKMorphologicalOperators1.Labels = ['POINTS', 'ComponentId']
tTKMorphologicalOperators1.Iterations = 0
tTKMorphologicalOperators1.UseGrayscaleOperators = 1

# create a new 'TTK CorrespondenceByOverlap'
tTKCorrespondenceByOverlap1 = TTKCorrespondenceByOverlap(Segmentations=tTKMorphologicalOperators1)
tTKCorrespondenceByOverlap1.Labels = ['POINTS', 'ComponentId']

# create a new 'TTK BlockAggregator'
tTKBlockAggregator1 = TTKBlockAggregator(Input=tTKCorrespondenceByOverlap1)

# create a new 'TTK BlockAggregator'
tTKBlockAggregator3 = TTKBlockAggregator(Input=[tTKBlockAggregator1, tTKBlockAggregator2, tTKBlockAggregator4, ])
tTKBlockAggregator3.ForceReset = 1
tTKBlockAggregator3.FlattenInput = 0

# create a new 'TTK EndFor'
tTKEndFor1 = TTKEndFor(Data=tTKBlockAggregator3,
    For=tTKForEach1)

# create a new 'TTK Extract'
tTKExtract2 = TTKExtract(Input=tTKEndFor1)
tTKExtract2.Expression = '0'
tTKExtract2.OutputType = 'vtkMultiBlockDataSet'
tTKExtract2.ImageExtent = [0, 84, 0, 84, 0, 0]
tTKExtract2.InputArray = ['POINTS', 'Label']

# create a new 'TTK Extract'
tTKExtract3 = TTKExtract(Input=tTKEndFor1)
tTKExtract3.Expression = '1'
tTKExtract3.OutputType = 'vtkMultiBlockDataSet'
tTKExtract3.ImageExtent = [0, 84, 0, 84, 0, 0]
tTKExtract3.InputArray = ['POINTS', 'Label']

# create a new 'TTK Extract'
tTKExtract4 = TTKExtract(Input=tTKEndFor1)
tTKExtract4.Expression = '2'
tTKExtract4.OutputType = 'vtkMultiBlockDataSet'

# create a new 'Threshold'
threshold1 = Threshold(Input=tTKExtract3)
threshold1.Scalars = ['POINTS', 'Size']
threshold1.ThresholdRange = [3.0, 877.0]

# create a new 'TTK TrackingGraph'
tTKTrackingGraph1 = TTKTrackingGraph(Correspondences=tTKExtract2,
    Features=threshold1)
tTKTrackingGraph1.Matrix = ['POINTS', 'Overlap']
tTKTrackingGraph1.Labels = ['POINTS', 'ComponentId']

# create a new 'Threshold'
threshold2 = Threshold(Input=tTKTrackingGraph1)
threshold2.Scalars = ['POINTS', 'Size']
threshold2.ThresholdRange = [0.0, 877.0]

# create a new 'TTK IcospheresFromPoints'
tTKIcospheresFromPoints1 = TTKIcospheresFromPoints(Input=threshold2)
tTKIcospheresFromPoints1.Radius = 0.2

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from tTKCorrespondenceByPersistencePairs1
tTKCorrespondenceByPersistencePairs1Display = Show(tTKCorrespondenceByPersistencePairs1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'LiftedWassersteinDistance'
liftedWassersteinDistanceLUT = GetColorTransferFunction('LiftedWassersteinDistance')
liftedWassersteinDistanceLUT.RGBPoints = [nan, 0.231373, 0.298039, 0.752941, nan, 0.865003, 0.865003, 0.865003, nan, 0.705882, 0.0156863, 0.14902]
liftedWassersteinDistanceLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
tTKCorrespondenceByPersistencePairs1Display.Representation = 'Points'
tTKCorrespondenceByPersistencePairs1Display.ColorArrayName = ['POINTS', 'LiftedWassersteinDistance']
tTKCorrespondenceByPersistencePairs1Display.LookupTable = liftedWassersteinDistanceLUT
tTKCorrespondenceByPersistencePairs1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKCorrespondenceByPersistencePairs1Display.SelectOrientationVectors = 'None'
tTKCorrespondenceByPersistencePairs1Display.ScaleFactor = -2.0000000000000002e+298
tTKCorrespondenceByPersistencePairs1Display.SelectScaleArray = 'None'
tTKCorrespondenceByPersistencePairs1Display.GlyphType = 'Arrow'
tTKCorrespondenceByPersistencePairs1Display.GlyphTableIndexArray = 'None'
tTKCorrespondenceByPersistencePairs1Display.GaussianRadius = -1e+297
tTKCorrespondenceByPersistencePairs1Display.SetScaleArray = [None, '']
tTKCorrespondenceByPersistencePairs1Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKCorrespondenceByPersistencePairs1Display.OpacityArray = [None, '']
tTKCorrespondenceByPersistencePairs1Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKCorrespondenceByPersistencePairs1Display.DataAxesGrid = 'GridAxesRepresentation'
tTKCorrespondenceByPersistencePairs1Display.PolarAxes = 'PolarAxesRepresentation'

# setup the color legend parameters for each legend in this view

# get color legend/bar for liftedWassersteinDistanceLUT in view renderView1
liftedWassersteinDistanceLUTColorBar = GetScalarBar(liftedWassersteinDistanceLUT, renderView1)
liftedWassersteinDistanceLUTColorBar.Title = 'LiftedWassersteinDistance'
liftedWassersteinDistanceLUTColorBar.ComponentTitle = ''

# set color bar visibility
liftedWassersteinDistanceLUTColorBar.Visibility = 1

# show color legend
tTKCorrespondenceByPersistencePairs1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from tTKTopologicalSimplification1
tTKTopologicalSimplification1Display = Show(tTKTopologicalSimplification1, renderView2, 'UniformGridRepresentation')

# get color transfer function/color map for 'Field'
fieldLUT = GetColorTransferFunction('Field')
fieldLUT.RGBPoints = [0.03379811346530914, 0.231373, 0.298039, 0.752941, 0.4880429282784462, 0.865003, 0.865003, 0.865003, 0.9422877430915833, 0.705882, 0.0156863, 0.14902]
fieldLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Field'
fieldPWF = GetOpacityTransferFunction('Field')
fieldPWF.Points = [0.03379811346530914, 0.0, 0.5, 0.0, 0.6627524495124817, 0.0, 0.5, 0.0, 0.9422877430915833, 0.680497944355011, 0.5, 0.0]
fieldPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
tTKTopologicalSimplification1Display.Representation = 'Volume'
tTKTopologicalSimplification1Display.ColorArrayName = ['POINTS', 'Field']
tTKTopologicalSimplification1Display.LookupTable = fieldLUT
tTKTopologicalSimplification1Display.OSPRayScaleArray = 'Field'
tTKTopologicalSimplification1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKTopologicalSimplification1Display.SelectOrientationVectors = 'None'
tTKTopologicalSimplification1Display.ScaleFactor = 6.300000000000001
tTKTopologicalSimplification1Display.SelectScaleArray = 'Field'
tTKTopologicalSimplification1Display.GlyphType = 'Arrow'
tTKTopologicalSimplification1Display.GlyphTableIndexArray = 'Field'
tTKTopologicalSimplification1Display.GaussianRadius = 0.315
tTKTopologicalSimplification1Display.SetScaleArray = ['POINTS', 'Field']
tTKTopologicalSimplification1Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKTopologicalSimplification1Display.OpacityArray = ['POINTS', 'Field']
tTKTopologicalSimplification1Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKTopologicalSimplification1Display.DataAxesGrid = 'GridAxesRepresentation'
tTKTopologicalSimplification1Display.PolarAxes = 'PolarAxesRepresentation'
tTKTopologicalSimplification1Display.ScalarOpacityUnitDistance = 1.7320508075688774
tTKTopologicalSimplification1Display.ScalarOpacityFunction = fieldPWF
tTKTopologicalSimplification1Display.SliceFunction = 'Plane'
tTKTopologicalSimplification1Display.Slice = 31

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKTopologicalSimplification1Display.ScaleTransferFunction.Points = [0.06124124675989151, 0.0, 0.5, 0.0, 0.9422877430915833, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKTopologicalSimplification1Display.OpacityTransferFunction.Points = [0.06124124675989151, 0.0, 0.5, 0.0, 0.9422877430915833, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
tTKTopologicalSimplification1Display.SliceFunction.Origin = [31.5, 31.5, 31.5]

# show data from calculator1
calculator1Display = Show(calculator1, renderView2, 'UniformGridRepresentation')

# trace defaults for the display properties.
calculator1Display.Representation = 'Outline'
calculator1Display.ColorArrayName = ['POINTS', '']
calculator1Display.LookupTable = fieldLUT
calculator1Display.OSPRayScaleArray = 'Label'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 6.300000000000001
calculator1Display.SelectScaleArray = 'Label'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'Label'
calculator1Display.GaussianRadius = 0.315
calculator1Display.SetScaleArray = ['POINTS', 'Label']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'Label']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.ScalarOpacityUnitDistance = 1.7320508075688774
calculator1Display.ScalarOpacityFunction = fieldPWF
calculator1Display.IsosurfaceValues = [-1.0]
calculator1Display.SliceFunction = 'Plane'
calculator1Display.Slice = 31

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
calculator1Display.SliceFunction.Origin = [31.5, 31.5, 31.5]

# show data from tTKIcospheresFromPoints1
tTKIcospheresFromPoints1Display = Show(tTKIcospheresFromPoints1, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'Time'
timeLUT = GetColorTransferFunction('Time')
timeLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 2.00048828125, 0.865003, 0.865003, 0.865003, 4.0009765625, 0.705882, 0.0156863, 0.14902]
timeLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
tTKIcospheresFromPoints1Display.Representation = 'Surface'
tTKIcospheresFromPoints1Display.ColorArrayName = ['POINTS', 'Time']
tTKIcospheresFromPoints1Display.LookupTable = timeLUT
tTKIcospheresFromPoints1Display.OSPRayScaleArray = 'Label'
tTKIcospheresFromPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints1Display.SelectOrientationVectors = 'None'
tTKIcospheresFromPoints1Display.ScaleFactor = 6.4
tTKIcospheresFromPoints1Display.SelectScaleArray = 'Label'
tTKIcospheresFromPoints1Display.GlyphType = 'Arrow'
tTKIcospheresFromPoints1Display.GlyphTableIndexArray = 'Label'
tTKIcospheresFromPoints1Display.GaussianRadius = 0.32
tTKIcospheresFromPoints1Display.SetScaleArray = ['POINTS', 'Label']
tTKIcospheresFromPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints1Display.OpacityArray = ['POINTS', 'Label']
tTKIcospheresFromPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'
tTKIcospheresFromPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tTKIcospheresFromPoints1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tTKIcospheresFromPoints1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 84.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tTKIcospheresFromPoints1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 84.0, 1.0, 0.5, 0.0]

# show data from threshold2
threshold2Display = Show(threshold2, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold2Display.Representation = 'Surface'
threshold2Display.ColorArrayName = ['POINTS', '']
threshold2Display.OSPRayScaleArray = 'Label'
threshold2Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold2Display.SelectOrientationVectors = 'None'
threshold2Display.ScaleFactor = 6.300000000000001
threshold2Display.SelectScaleArray = 'Label'
threshold2Display.GlyphType = 'Arrow'
threshold2Display.GlyphTableIndexArray = 'Label'
threshold2Display.GaussianRadius = 0.315
threshold2Display.SetScaleArray = ['POINTS', 'Label']
threshold2Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold2Display.OpacityArray = ['POINTS', 'Label']
threshold2Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold2Display.DataAxesGrid = 'GridAxesRepresentation'
threshold2Display.PolarAxes = 'PolarAxesRepresentation'
threshold2Display.ScalarOpacityUnitDistance = 12.494150575757386

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
threshold2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 216.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
threshold2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 216.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'LiftedWassersteinDistance'
liftedWassersteinDistancePWF = GetOpacityTransferFunction('LiftedWassersteinDistance')
liftedWassersteinDistancePWF.Points = [nan, 0.0, 0.5, 0.0, nan, 1.0, 0.5, 0.0]
liftedWassersteinDistancePWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'Time'
timePWF = GetOpacityTransferFunction('Time')
timePWF.Points = [0.0, 0.0, 0.5, 0.0, 4.0009765625, 1.0, 0.5, 0.0]
timePWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(tTKEndFor1)
# ----------------------------------------------------------------