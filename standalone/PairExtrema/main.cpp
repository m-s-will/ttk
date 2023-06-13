/// \author Julien Tierny <julien.tierny@sorbonne-universite.fr>.
/// \date February 2017.
///
/// \brief Command line program for critical point computation.

// TTK Includes
#include <CommandLineParser.h>
#include <ttkFTMTree.h>
#include <ttkPairExtrema.h>
#include <ttkPathCompression.h>
#include <ttkScalarFieldCriticalPoints2.h>

// VTK Includes
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLDataObjectWriter.h>
#include <vtkXMLGenericDataObjectReader.h>

int main(int argc, char **argv) {

  // ---------------------------------------------------------------------------
  // Program variables
  // ---------------------------------------------------------------------------
  std::string inputFilePath = "";
  std::string inputArrayName = "";
  std::string outputPathPrefix{"output"};
  bool listArrays{false};
  bool compare{false};
  int threadNumber{112};
  int repetitions{1};

  // ---------------------------------------------------------------------------
  // Set program variables based on command line arguments
  // ---------------------------------------------------------------------------
  {
    ttk::CommandLineParser parser;

    // -------------------------------------------------------------------------
    // Standard options and arguments
    // -------------------------------------------------------------------------
    parser.setArgument(
      "i", &inputFilePath, "Input data-set (*.vti, *vtu, *vtp)", false);
    parser.setArgument("a", &inputArrayName, "Input array name", false);
    parser.setArgument("r", &repetitions,
                       "Number of times you want to run the algorithms", true);
    parser.setOption("l", &listArrays, "List available arrays");
    parser.setOption("c", &compare, "Also build the TTK FTM Tree to compare");
    parser.setArgument(
      "n", &threadNumber, "The number of OMP Threads to run the filters", true);

    parser.parse(argc, argv);
  }

  // ---------------------------------------------------------------------------
  // Command line output messages.
  // ---------------------------------------------------------------------------
  ttk::Debug msg;
  msg.setDebugMsgPrefix("PairExtrema");

  // ---------------------------------------------------------------------------
  // Initialize ttkPairExtrema module (adjust parameters)
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Read input vtkDataObjects (optionally: print available arrays)
  // ---------------------------------------------------------------------------
  vtkDataArray *defaultArray = nullptr;
  // init a reader that can parse any vtkDataObject stored in xml format
  auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
  reader->SetFileName(inputFilePath.data());
  reader->Update();

  // check if input vtkDataObject was successfully read
  auto inputDataObject = reader->GetOutput();
  if(!inputDataObject) {
    msg.printErr("Unable to read input file `" + inputFilePath + "' :(");
    return 1;
  }

  auto inputAsVtkDataSet = vtkDataSet::SafeDownCast(inputDataObject);

  // if requested print list of arrays, otherwise proceed with execution
  if(listArrays) {
    msg.printMsg(inputFilePath + ":");
    if(inputAsVtkDataSet) {
      // Point Data
      msg.printMsg("  PointData:");
      auto pointData = inputAsVtkDataSet->GetPointData();
      for(int j = 0; j < pointData->GetNumberOfArrays(); j++)
        msg.printMsg("    - " + std::string(pointData->GetArrayName(j)));

      // Cell Data
      msg.printMsg("  CellData:");
      auto cellData = inputAsVtkDataSet->GetCellData();
      for(int j = 0; j < cellData->GetNumberOfArrays(); j++)
        msg.printMsg("    - " + std::string(cellData->GetArrayName(j)));
      return 0;
    } else {
      msg.printErr("Unable to list arrays on file `" + inputFilePath + "'");
      return 1;
    }
  }

  // ---------------------------------------------------------------------------
  // Specify which arrays of the input vtkDataObjects will be processed
  // ---------------------------------------------------------------------------
  if(!defaultArray) {
    defaultArray = inputAsVtkDataSet->GetPointData()->GetArray(0);
  }
  if(inputArrayName == "") {
    if(defaultArray)
      inputArrayName = defaultArray->GetName();
  }

  msg.setDebugMsgPrefix("PathCompression");
  for(int i = 0; i < repetitions; i++) {
    auto pathCompression = vtkSmartPointer<ttkPathCompression>::New();
    auto criticalPoints = vtkSmartPointer<ttkScalarFieldCriticalPoints2>::New();
    auto pairExtrema = vtkSmartPointer<ttkPairExtrema>::New();

    pathCompression->SetInputDataObject(0, reader->GetOutput());
    pathCompression->SetComputeFinalSegmentation(0);
    pathCompression->SetInputArrayToProcess(0, 0, 0, 0, inputArrayName.data());
    pathCompression->SetUseAllCores(false);
    pathCompression->SetThreadNumber(threadNumber);
    pathCompression->Update();

    msg.setDebugMsgPrefix("CriticalPoints");
    criticalPoints->SetInputDataObject(0, pathCompression->GetOutput());
    criticalPoints->SetInputArrayToProcess(0, 0, 0, 0, inputArrayName.data());
    criticalPoints->SetUseAllCores(false);
    criticalPoints->SetThreadNumber(threadNumber);
    criticalPoints->Update();

    pairExtrema->SetInputDataObject(0, pathCompression->GetOutput());
    pairExtrema->SetInputDataObject(1, criticalPoints->GetOutput());

    pairExtrema->SetInputArrayToProcess(0, 0, 0, 0, inputArrayName.data());

    // ---------------------------------------------------------------------------
    // Execute ttkPairExtrema filter
    // ---------------------------------------------------------------------------
    pairExtrema->SetUseAllCores(false);
    pairExtrema->SetThreadNumber(threadNumber);
    pairExtrema->Modified();
    pairExtrema->Update();
  }

  if(compare) {
    msg.setDebugMsgPrefix("FTMTree");
    auto ftmTree = vtkSmartPointer<ttkFTMTree>::New();
    ftmTree->SetInputDataObject(0, reader->GetOutput());
    ftmTree->SetInputArrayToProcess(0, 0, 0, 0, inputArrayName.data());
    ftmTree->SetTreeType(1); // Split tree
    ftmTree->SetUseAllCores(false);
    ftmTree->SetThreadNumber(threadNumber);
    for(int i = 0; i < repetitions; i++) {
      ftmTree->Modified();
      ftmTree->Update();
    }
  }
  // ---------------------------------------------------------------------------
  // If output prefix is specified then write all output objects to disk
  // ---------------------------------------------------------------------------
  /*if(!outputPathPrefix.empty()) {
    for(int i = 0; i < pairExtrema->GetNumberOfOutputPorts();
        i++) {
      auto output = pairExtrema->GetOutputDataObject(i);
      auto writer = vtkSmartPointer<vtkXMLWriter>::Take(
        vtkXMLDataObjectWriter::NewWriter(output->GetDataObjectType()));

      std::string outputFileName = outputPathPrefix + "_port_"
                                   + std::to_string(i) + "."
                                   + writer->GetDefaultFileExtension();
      msg.printMsg("Writing output file `" + outputFileName + "'...");
      writer->SetInputDataObject(output);
      writer->SetFileName(outputFileName.data());
      writer->Update();
    }
  }*/

  return 0;
}
