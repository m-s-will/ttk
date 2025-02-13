#include <ttkTemporalMergeTreeMap.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkImageData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkTemporalMergeTreeMap);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkTemporalMergeTreeMap::ttkTemporalMergeTreeMap() {
  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(1);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkTemporalMergeTreeMap::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  if(port == 2) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

/**
 * TODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkTemporalMergeTreeMap::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

void dfs_linearization(
  int curr_node,
  std::vector<double> &lin,
  std::vector<int> &seg,
  std::vector<int> &bar,
  std::vector<std::vector<int>> &memiChildren,
  std::vector<std::vector<double>> &memiSegmentScalars,
  std::vector<int> &memiSizes,
  std::vector<int> &memiSegs,
  std::vector<int> &branchNodeIDs,
  std::vector<double> &memiScalars,
  std::vector<double> &memiOrdering){
      if(memiChildren[curr_node].size()==0){
        auto curr_size = memiSizes[curr_node];
        auto segment = memiSegmentScalars[memiSegs[curr_node]];
        for(int s=0; s<segment.size(); s+=2){
          lin.push_back(segment[s]);
          seg.push_back(memiSegs[curr_node]);
          bar.push_back(branchNodeIDs[curr_node]);
        }
        for(int s=segment.size()-1-(segment.size()%2); s>=0; s-=2){
        // for(int s=1; s<segment.size(); s+=2){
          lin.push_back(segment[s]);
          seg.push_back(memiSegs[curr_node]);
          bar.push_back(branchNodeIDs[curr_node]);
        }
      }
      else{
        std::vector<int> curr_children = memiChildren[curr_node];
        std::sort(curr_children.begin(),curr_children.end(),[memiOrdering](int x1,int x2)->bool{return memiOrdering[x1]<memiOrdering[x1];});
        // std::vector<int> curr_children;
        auto curr_size = int(memiSizes[curr_node]);
        auto segment = memiSegmentScalars[memiSegs[curr_node]];
        auto cs1 = int(curr_size/2);
        auto cs2 = curr_size-cs1;
        for(int s=0; s<segment.size()-1; s+=2){
          lin.push_back(segment[s]);
          seg.push_back(memiSegs[curr_node]);
          bar.push_back(branchNodeIDs[curr_node]);
        }
        if(curr_children.size()>2){
          std::cout << "!! " << curr_children.size() << std::endl;
        }
        for(int ci=0; ci<curr_children.size(); ci++){
          auto c = curr_children[ci];
          dfs_linearization(c,lin,seg,bar,memiChildren,memiSegmentScalars,memiSizes,memiSegs,branchNodeIDs,memiScalars,memiOrdering);
          if(ci < curr_children.size()-1){
            lin.push_back(memiScalars[curr_node]);
            seg.push_back(memiSegs[curr_node]);
            bar.push_back(branchNodeIDs[curr_node]);
          }
        }
        for(int s=segment.size()-2-(1-(segment.size()%2)); s>=0; s-=2){
        // for(int s=1; s<segment.size()-1; s+=2){
          lin.push_back(segment[s]);
          seg.push_back(memiSegs[curr_node]);
          bar.push_back(branchNodeIDs[curr_node]);
        }
      }
}

/**
 * TODO 10: Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkTemporalMergeTreeMap::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  auto members = vtkMultiBlockDataSet::GetData(inputVector[0]);
  if(!members)
    return 0;
  auto barycenter = vtkMultiBlockDataSet::GetData(inputVector[1]);
  if(!barycenter)
    return 0;
  auto domains = vtkMultiBlockDataSet::GetData(inputVector[2]);
  if(!domains)
    return 0;
  // std::cout << members->GetNumberOfBlocks() <<std::endl;
  // std::cout << barycenter->GetNumberOfBlocks() <<std::endl;
  // std::cout << domains->GetNumberOfBlocks() <<std::endl;

  auto memberNodes = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(0));
  auto memberArcs = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(1));
  auto baryNodes_ = vtkMultiBlockDataSet::SafeDownCast(barycenter->GetBlock(0));
  auto baryArcs_ = vtkMultiBlockDataSet::SafeDownCast(barycenter->GetBlock(1));
  auto baryNodes = vtkUnstructuredGrid::SafeDownCast(baryNodes_->GetBlock(0));
  auto baryArcs = vtkUnstructuredGrid::SafeDownCast(baryArcs_->GetBlock(0));

  // std::cout << memberNodes->GetNumberOfBlocks() <<std::endl;
  // std::cout << memberArcs->GetNumberOfBlocks() <<std::endl;
  // std::cout << baryNodes_->GetNumberOfBlocks() <<std::endl;
  // std::cout << baryArcs_->GetNumberOfBlocks() <<std::endl;
  // std::cout << baryNodes->GetNumberOfPoints() <<std::endl;
  // std::cout << baryNodes->GetNumberOfCells() <<std::endl;
  // std::cout << baryArcs->GetNumberOfPoints() <<std::endl;
  // std::cout << baryArcs->GetNumberOfCells() <<std::endl;

  // compute ordering of barycenter branchIDs (from tree layout/actual positions)
  auto ordering_branches = std::vector<double>(baryNodes->GetNumberOfPoints(),-1);
  for(int i=0; i<baryNodes->GetNumberOfPoints(); i++){
    auto posX = baryNodes->GetPoint(i)[0];
    // auto posY = baryNodes->GetPoint(i)[1];
    // auto nId = baryNodes->GetPointData()->GetArray("NodeId")->GetComponent(i,0);
    auto bId = vtkIntArray::SafeDownCast(baryNodes->GetPointData()->GetArray("BranchNodeID"))->GetValue(i);
    ordering_branches[bId] = posX;
  }

  std::vector<std::vector<double>> linearizations;
  std::vector<std::vector<int>> segmentations;
  std::vector<std::vector<int>> barycenterRefs;
  int maxlen = 0;

  for(int blockIdx=0; blockIdx < memberNodes->GetNumberOfBlocks(); blockIdx++){
    int totalSize = 0;
    std::vector<int> arcRegions;
    auto memiNodes = vtkUnstructuredGrid::SafeDownCast(memberNodes->GetBlock(blockIdx));
    auto memiArcs = vtkUnstructuredGrid::SafeDownCast(memberArcs->GetBlock(blockIdx));
    auto memiDomain = vtkDataSet::SafeDownCast(domains->GetBlock(blockIdx));
    // # print(" ")

    // prepare scalar values in segments
    std::vector<std::vector<double>> memiSegmentScalars(memiArcs->GetNumberOfCells());
    for(int j=0; j<memiDomain->GetNumberOfPoints(); j++){
      auto scalar = memiDomain->GetPointData()->GetArray("nrrd")->GetComponent(j,0);
      auto segId = vtkIntArray::SafeDownCast(memiDomain->GetPointData()->GetArray("SegmentationId"))->GetValue(j);
      memiSegmentScalars[segId].push_back(scalar);
    }
    for(int i=0; i<memiSegmentScalars.size(); i++){
      auto l = memiSegmentScalars[i];
      std::sort(memiSegmentScalars[i].begin(),memiSegmentScalars[i].end());
    }

    // get node properties of member tree
    int numNodesi = 0;
    std::vector<int> memiNodeIsDummy(memiNodes->GetNumberOfPoints());
    std::vector<double> memiScalars(memiNodes->GetNumberOfPoints());
    for(int i=0; i<memiNodes->GetNumberOfPoints(); i++){
      auto nId = vtkIntArray::SafeDownCast(memiNodes->GetPointData()->GetArray("NodeId"))->GetValue(i);
      auto isDummy = vtkIntArray::SafeDownCast(memiNodes->GetPointData()->GetArray("isDummyNode"))->GetValue(i);
      auto scalar = memiNodes->GetPointData()->GetArray("Scalar")->GetComponent(i,0);
      if(!isDummy){
        memiNodeIsDummy[nId] = 1;
        memiScalars[nId] = scalar;
        numNodesi += 1;
      }
    }

    // create tree structure of member tree (parent pointers)
    std::vector<int> memiParents(memiNodes->GetNumberOfPoints(),-1);
    for(int i=0; i<memiArcs->GetNumberOfCells(); i++){
      auto childId = vtkIntArray::SafeDownCast(memiArcs->GetCellData()->GetArray("downNodeId"))->GetValue(i);
      auto parentId = vtkIntArray::SafeDownCast(memiArcs->GetCellData()->GetArray("upNodeId"))->GetValue(i);
      memiParents[childId] = parentId;
    }

    //  create tree structure of member tree (children lists)
    std::vector<std::vector<int>> memiChildren(memiNodes->GetNumberOfPoints());
    int root = -1;
    for(int i=0; i<memiParents.size(); i++){
      auto parentId = memiParents[i];
      if(parentId >= 0)
        memiChildren[parentId].push_back(i);
      if(parentId==-1 && memiNodeIsDummy[i])
        root = i;
    }
    
    // get barycenter branchIDs of member tree nodes
    std::vector<int> branchNodeIDs(memiNodes->GetNumberOfPoints(),-1);
    for(int j=0; j<memiNodes->GetNumberOfPoints(); j++){
      auto nId = vtkIntArray::SafeDownCast(memiNodes->GetPointData()->GetArray("NodeId"))->GetValue(j);
      auto bId = vtkIntArray::SafeDownCast(memiNodes->GetPointData()->GetArray("BranchBaryNodeID"))->GetValue(j);
      branchNodeIDs[nId] = bId;
    }

    // compute ordering of member tree nodes for layout (derived from barycenter branchIDs of leaves through dfs)
    std::vector<double> memiOrdering(memiNodes->GetNumberOfPoints(),-1);
    std::stack<int> s1;
    std::stack<int> s2;
    std::vector<int> postorder;
    std::vector<int> preorder;
    s1.push(root);
    while(!s1.empty()){
      auto node = s1.top();
      s1.pop();
      s2.push(node);
      preorder.push_back(node);
      for(auto child : memiChildren[node]){
        s1.push(child);
      }
    }
    while(!s2.empty()){
      auto node = s2.top();
      s2.pop();
      postorder.push_back(node);
    }
    for(auto node : postorder){
      if(memiChildren[node].size()==0){
        memiOrdering[node] = ordering_branches[branchNodeIDs[node]];
        // return memiOrdering[node];
      }
      else{
        double oidx = std::numeric_limits<double>::infinity();
        for(auto child : memiChildren[node]){
          auto cidx = memiOrdering[child];
          if(cidx < oidx)
            oidx = cidx;
        }
        memiOrdering[node] = oidx;
        // return oidx;
      }
    }

    // # get area sizes and segmentation ids for member tree nodes (from parent edge)
    std::vector<int> memiSizes(memiNodes->GetNumberOfPoints());
    std::vector<int> memiSegs(memiNodes->GetNumberOfPoints());
    for(int j=0; j<memiArcs->GetNumberOfCells(); j++){
      auto rS = memiArcs->GetCellData()->GetArray("RegionSize")->GetComponent(j,0);
      auto isDummy = vtkIntArray::SafeDownCast(memiArcs->GetCellData()->GetArray("isDummyArc"))->GetValue(j);
      auto downNodeId = vtkIntArray::SafeDownCast(memiArcs->GetCellData()->GetArray("downNodeId"))->GetValue(j);
      int segId = memiArcs->GetCellData()->GetArray("SegmentationId")->GetComponent(j,0);
      memiSizes[downNodeId] = int(rS);
      memiSegs[downNodeId] = int(segId);
    }

    // # compute linearization based on odering
    std::vector<double> linearization;
    std::vector<int> segmentation;
    std::vector<int> barycenterRef;
    // def dfs_linearization(curr_node,lin,seg,bar):
    //   if len(memiChildren[curr_node])==0:
    //     curr_size = memiSizes[curr_node]
    //     segment = memiSegmentScalars[memiSegs[curr_node]]
    //     lin += [segment[i] for i in range(0,len(segment),2)]
    //     seg += [memiSegs[curr_node] for i in range(0,len(segment),2)]
    //     bar += [branchNodeIDs[curr_node] for i in range(0,len(segment),2)]
    //     lin += [segment[i] for i in reversed(range(1,len(segment),2))]
    //     seg += [memiSegs[curr_node] for i in reversed(range(1,len(segment),2))]
    //     bar += [branchNodeIDs[curr_node] for i in reversed(range(1,len(segment),2))]
    //   else:
    //     curr_children = sorted(memiChildren[curr_node],key=lambda x: memiOrdering[x])
    //     curr_size = int(memiSizes[curr_node])
    //     segment = memiSegmentScalars[memiSegs[curr_node]]
    //     cs1 = int(curr_size/2)
    //     cs2 = curr_size-cs1
    //     lin += [segment[i] for i in range(0,len(segment)-1,2)]
    //     seg += [memiSegs[curr_node] for i in range(0,len(segment)-1,2)]
    //     bar += [branchNodeIDs[curr_node] for i in range(0,len(segment)-1,2)]
    //     if len(curr_children)>2:
    //       print("!!",len(curr_children))
    //     for ci in range(len(curr_children)):
    //       c = curr_children[ci]
    //       dfs_linearization(c,lin,seg,bar)
    //       if ci < len(curr_children)-1:
    //         lin.append(memiScalars[curr_node])
    //         seg.append(memiSegs[curr_node])
    //         bar.append(branchNodeIDs[curr_node])
    //     lin += [segment[i] for i in reversed(range(1,len(segment)-1,2))]
    //     seg += [memiSegs[curr_node] for i in reversed(range(1,len(segment)-1,2))]
    //     bar += [branchNodeIDs[curr_node] for i in reversed(range(1,len(segment)-1,2))]
    dfs_linearization(memiChildren[root][0],linearization,segmentation,barycenterRef,memiChildren,memiSegmentScalars,memiSizes,memiSegs,branchNodeIDs,memiScalars,memiOrdering);
    linearizations.push_back(linearization);
    segmentations.push_back(segmentation);
    barycenterRefs.push_back(barycenterRef);
    // print(len(linearization),len(segmentation))
    if (linearization.size() > maxlen){
      maxlen = linearization.size();
    }
    std::cout << linearization.size() << " " << memiDomain->GetNumberOfPoints() << " " << memiNodes->GetNumberOfPoints() << std::endl;
  }
  std::cout << maxlen << std::endl;

  // write linearization to output vti

  auto outputmb = vtkMultiBlockDataSet::GetData(outputVector,0);
  outputmb->SetNumberOfBlocks(1);
  vtkNew<vtkImageData> tmtm;
  tmtm->SetDimensions(linearizations.size()+1,maxlen+1,1);
  tmtm->SetSpacing(1024,1,1);
  tmtm->SetOrigin(0,0,0);

  vtkNew<vtkFloatArray> linArray{};
  vtkNew<vtkIntArray> segArray{};
  vtkNew<vtkIntArray> blockArray{};
  vtkNew<vtkIntArray> barArray{};

  linArray->SetName("Scalar");
  linArray->SetNumberOfComponents(1);
  linArray->SetNumberOfTuples(linearizations.size()*maxlen);

  segArray->SetName("SegmentationId");
  segArray->SetNumberOfComponents(1);
  segArray->SetNumberOfTuples(linearizations.size()*maxlen);

  blockArray->SetName("BlockId");
  blockArray->SetNumberOfComponents(1);
  blockArray->SetNumberOfTuples(linearizations.size()*maxlen);

  barArray->SetName("BarycenterBranchID");
  barArray->SetNumberOfComponents(1);
  barArray->SetNumberOfTuples(linearizations.size()*maxlen);

  int k = 0;
  for(int j=0; j<maxlen; j++){
    for(int i=0; i<linearizations.size(); i++){
      auto v = j < linearizations[i].size() ? linearizations[i][j] : 0;
      auto s = j < linearizations[i].size() ? segmentations[i][j]: 0;
      auto b = j < linearizations[i].size() ? barycenterRefs[i][j]: 0;
      linArray->SetValue(k,v);
      segArray->SetValue(k,s);
      blockArray->SetValue(k,i);
      barArray->SetValue(k,b);
      k += 1;
    }
  }
  tmtm->GetCellData()->AddArray(linArray);
  tmtm->GetCellData()->AddArray(segArray);
  tmtm->GetCellData()->AddArray(blockArray);
  tmtm->GetCellData()->AddArray(barArray);

  outputmb->SetBlock(0,tmtm);

  // Get input array that will be processed
  //
  // Note: VTK provides abstract functionality to handle array selections, but
  //       this essential functionality is unfortunately not well documented.
  //       Before you read further, please keep in mind the the TTK developer
  //       team is not responsible for the existing VTK Api ;-)
  //
  //       In a nutshell, prior to the RequestData execution one has to call
  //
  //           SetInputArrayToProcess (
  //               int idx,
  //               int port,
  //               int connection,
  //               int fieldAssociation,
  //               const char *name
  //            )
  //
  //       The parameter 'idx' is often misunderstood: lets say the filter
  //       requires n arrays, then idx enumerates them from 0 to n-1.
  //
  //       The 'port' is the input port index at which the object is connected
  //       from which we want to get the array.
  //
  //       The 'connection' is the connection index at that port (we have to
  //       specify this because VTK allows multiple connections at the same
  //       input port).
  //
  //       The 'fieldAssociation' integer specifies if the array should be taken
  //       from 0: point data, 1: cell data, or 2: field data.
  //
  //       The final parameter is the 'name' of the array.
  //
  //       Example: SetInputArrayToProcess(3,1,0,1,"EdgeLength") will store that
  //                for the 3rd array the filter needs the cell data array named
  //                "EdgeLength" that it will retrieve from the vtkDataObject
  //                at input port 1 (first connection). During the RequestData
  //                method one can then actually retrieve the 3rd array it
  //                requires for its computation by calling
  //                GetInputArrayToProcess(3, inputVector)
  //
  //       If this filter is run within ParaView, then the UI will automatically
  //       call SetInputArrayToProcess (see TemporalMergeTreeMap.xml file).
  //
  //       During the RequestData execution one can then retrieve an actual
  //       array with the method "GetInputArrayToProcess".
  // vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);
  // if(!inputArray) {
  //   this->printErr("Unable to retrieve input array.");
  //   return 0;
  // }

  // // To make sure that the selected array can be processed by this filter,
  // // one should also check that the array association and format is correct.
  // if(this->GetInputArrayAssociation(0, inputVector) != 0) {
  //   this->printErr("Input array needs to be a point data array.");
  //   return 0;
  // }
  // if(inputArray->GetNumberOfComponents() != 1) {
  //   this->printErr("Input array needs to be a scalar array.");
  //   return 0;
  // }

  // // If all checks pass then log which array is going to be processed.
  // this->printMsg("Starting computation...");
  // this->printMsg("  Scalar Array: " + std::string(inputArray->GetName()));

  // // Create an output array that has the same data type as the input array
  // // Note: vtkSmartPointers are well documented
  // //       (https://vtk.org/Wiki/VTK/Tutorials/SmartPointers)
  // vtkSmartPointer<vtkDataArray> const outputArray
  //   = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  // outputArray->SetName(this->OutputArrayName.data()); // set array name
  // outputArray->SetNumberOfComponents(1); // only one component per tuple
  // outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());

  // // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // // not exist already).
  // ttk::Triangulation *triangulation
  //   = ttkAlgorithm::GetTriangulation(inputDataSet);
  // if(!triangulation)
  //   return 0;

  // // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  // this->preconditionTriangulation(triangulation); // implemented in base class

  // // Templatize over the different input array data types and call the base code
  // int status = 0; // this integer checks if the base code returns an error
  // ttkVtkTemplateMacro(inputArray->GetDataType(), triangulation->getType(),
  //                     (status = this->computeAverages<VTK_TT, TTK_TT>(
  //                        (VTK_TT *)ttkUtils::GetVoidPointer(outputArray),
  //                        (VTK_TT *)ttkUtils::GetVoidPointer(inputArray),
  //                        (TTK_TT *)triangulation->getData())));

  // // On error cancel filter execution
  // if(status != 1)
  //   return 0;

  // // Get output vtkDataSet (which was already instantiated based on the
  // // information provided by FillOutputPortInformation)
  // vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // // make a SHALLOW copy of the input
  // outputDataSet->ShallowCopy(inputDataSet);

  // // add to the output point data the computed output array
  // outputDataSet->GetPointData()->AddArray(outputArray);

  // return success
  return 1;
}
