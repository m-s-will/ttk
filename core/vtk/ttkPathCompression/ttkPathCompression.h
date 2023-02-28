/// \ingroup vtk
/// \class ttkPathCompression
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK VTK-filter that wraps the pathCompression processing package.
///
/// TTK module for the computation of Morse-Smale complexes.
/// Morse-Smale complexes are useful topological abstractions of scalar
/// fields for data segmentation, feature extraction, etc.
///
/// \b Related \b publication \n
/// "Parallel Computation of 3D Morse-Smale Complexes" \n
/// Nithin Shivashankar, Vijay Natarajan \n
/// Proc. of EuroVis 2012. \n
/// Computer Graphics Forum, 2012.
///
/// \param Input Input scalar field, defined as a point data scalar field
/// attached to a geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output0 Output critical points (vtkPolyData)
/// \param Output1 Output 1-separatrices (vtkPolyData)
/// \param Output2 Output 2-separatrices (vtkPolyData)
/// \param Output3 Output data segmentation (vtkDataSet)
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// The optional offset array can be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 1 (FIXED: the second array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the offset array)
/// \note: To use this optional array, `ForceInputOffsetScalarField` needs to be
/// enabled with the setter `setForceInputOffsetScalarField()'.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::PathCompression
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearning/">1-Manifold
///   Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearningCircles/">1-Manifold
///   Learning Circles example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/2manifoldLearning/">
///   2-Manifold Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/imageProcessing/">Image
///   Processing example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/karhunenLoveDigits64Dimensions/">Karhunen-Love
///   Digits 64-Dimensions example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseMolecule/">Morse
///   molecule example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morsePersistence/">Morse
///   Persistence example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 0 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering1/">Persistence
///   clustering 1 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering2/">Persistence
///   clustering 2 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering3/">Persistence
///   clustering 3 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering4/">Persistence
///   clustering 4 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_at/">Persistent
///   Generators AT example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistentGenerators_darkSky/">Persistent
///   Generators DarkSky example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tectonicPuzzle/">Tectonic
///   Puzzle example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tribute/">Tribute
///   example</a> \n
///

#pragma once

// VTK Module
#include <ttkPathCompressionModule.h>

// ttk code includes
#include <PathCompression.h>
#include <ttkAlgorithm.h>
#include <ttkMacros.h>

class vtkPolyData;

class TTKPATHCOMPRESSION_EXPORT ttkPathCompression
  : public ttkAlgorithm,
    protected ttk::PathCompression {

public:
  static ttkPathCompression *New();

  vtkTypeMacro(ttkPathCompression, ttkAlgorithm);

  vtkSetMacro(ComputeAscendingSegmentation, bool);
  vtkGetMacro(ComputeAscendingSegmentation, bool);

  vtkSetMacro(ComputeDescendingSegmentation, bool);
  vtkGetMacro(ComputeDescendingSegmentation, bool);

  vtkSetMacro(ComputeFinalSegmentation, bool);
  vtkGetMacro(ComputeFinalSegmentation, bool);

  vtkSetMacro(ForceInputOffsetScalarField, bool);
  vtkGetMacro(ForceInputOffsetScalarField, bool);

protected:
  template <typename scalarType, typename triangulationType>
  int dispatch(vtkDataArray *const inputScalars,
               const SimplexId *const inputOffsets,
               const triangulationType &triangulation);

  ttkPathCompression();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  bool ForceInputOffsetScalarField{};
  OutputManifold segmentations_{};
};
