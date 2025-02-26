#include <ttkTemporalMergeTreeMap2.h>

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
#include <Debug.h>
#include <vtkCellData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <FTMTreeUtils.h>
#include <ttkMergeTreeUtils.h>
#include <ttkMergeTreeClustering.h>
#include <MergeTreeBarycenter.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkTemporalMergeTreeMap2);

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
ttkTemporalMergeTreeMap2::ttkTemporalMergeTreeMap2() {
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
int ttkTemporalMergeTreeMap2::FillInputPortInformation(int port, vtkInformation *info) {
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
int ttkTemporalMergeTreeMap2::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

void ttkTemporalMergeTreeMap2::dfs_linearization(
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
        std::sort(curr_children.begin(),curr_children.end(),[memiOrdering](int x1,int x2)->bool{return memiOrdering[x1]<memiOrdering[x2];});
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
        //if(curr_children.size()>2){
        //  std::cout << "!! " << curr_children.size() << std::endl;
        //}
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

void ttkTemporalMergeTreeMap2::computeBaryBranchOrdering(vtkMultiBlockDataSet* mtmb, vtkMultiBlockDataSet* members,std::vector<double> &ordering_branches){
  vtkNew<ttkMergeTreeClustering> c;
  c->SetInputDataObject(0,mtmb);
  c->SetComputeBarycenter(true);
  c->SetImportantPairs(0);
  c->SetEpsilonTree1(0);
  c->SetEpsilon2Tree1(100);
  c->SetEpsilon3Tree1(100);
  c->SetDeterministic(true);
  c->SetUseFixedInit(true);
  c->SetPlanarLayout(true);
  c->SetDebugLevel(4);
  c->SetBarycenterSizeLimitPercent(barycenterSize);
  c->Update();

  members->DeepCopy(vtkMultiBlockDataSet::SafeDownCast(c->GetOutputDataObject(0)));
  auto barycenter = vtkMultiBlockDataSet::SafeDownCast(c->GetOutputDataObject(1));

  auto baryNodes_ = vtkMultiBlockDataSet::SafeDownCast(barycenter->GetBlock(0));
  auto baryArcs_ = vtkMultiBlockDataSet::SafeDownCast(barycenter->GetBlock(1));
  auto baryNodes = vtkUnstructuredGrid::SafeDownCast(baryNodes_->GetBlock(0));
  auto baryArcs = vtkUnstructuredGrid::SafeDownCast(baryArcs_->GetBlock(0));

  int baryNumNodes = 0;
  std::vector<int> baryNodeIsDummy(baryNodes->GetNumberOfPoints());
  for(int i=0; i<baryNodes->GetNumberOfPoints(); i++){
    auto nId = baryNodes->GetPointData()->GetArray("NodeId")->GetComponent(i,0);
    auto isDummy = baryNodes->GetPointData()->GetArray("isDummyNode")->GetComponent(i,0);
    if(!isDummy){
      baryNodeIsDummy[nId] = 1;
      baryNumNodes += 1;
    }
  }
  //std::cout << "survived the loop" << std::endl;

  std::vector<int> baryParents(baryNodes->GetNumberOfPoints(),-1);
  for(int i=0; i<baryArcs->GetNumberOfCells(); i++){
    auto childId = (baryArcs->GetCellData()->GetArray("downNodeId"))->GetComponent(i,0);
    auto parentId = (baryArcs->GetCellData()->GetArray("upNodeId"))->GetComponent(i,0);
    baryParents[childId] = parentId;
  }

  std::vector<std::vector<int>> baryChildren(baryNodes->GetNumberOfPoints());
  int root = -1;
  for(int i=0; i<baryParents.size(); i++){
    auto parentId = baryParents[i];
    if(parentId >= 0)
      baryChildren[parentId].push_back(i);
    if(parentId==-1 && baryNodeIsDummy[i])
      root = i;
  }

  std::stack<int> s;
  s.push(root);
  int idx = 0;
  auto ordering_baryNodes = std::vector<double>(baryNodes->GetNumberOfPoints(),-1);
  while(!s.empty()){
    auto node = s.top();
    s.pop();
    ordering_baryNodes[node] = idx;
    idx++;
    for(auto child : baryChildren[node]){
      s.push(child);
    }
  }

  // compute ordering of barycenter branchIDs
  ordering_branches = std::vector<double>(baryNodes->GetNumberOfPoints(),-1);
  auto bary_branches = std::vector<double>(baryNodes->GetNumberOfPoints(),-1);
  auto bary_scalars = std::vector<double>(baryNodes->GetNumberOfPoints(),-1);
  for(int i=0; i<baryNodes->GetNumberOfPoints(); i++){
    auto posX = baryNodes->GetPoint(i)[0];
    // auto posY = baryNodes->GetPoint(i)[1];
    int nId = baryNodes->GetPointData()->GetArray("NodeId")->GetComponent(i,0);
    int scalar = baryNodes->GetPointData()->GetArray("Scalar")->GetComponent(i,0);
    // if(!baryChildren[nId].empty()) continue;
    auto bId = (baryNodes->GetPointData()->GetArray("BranchNodeID"))->GetComponent(i,0);
    ordering_branches[bId] = posX; //ordering_baryNodes[nId];
    bary_branches[nId] = bId;
    bary_scalars[nId] = scalar;
  }

  // std::cout << "bary" << "-----\n  ";
  // for(int j=0; j<bary_scalars.size(); j++){
  //   std::cout << j << ":";
  //   std::cout << baryParents[j] << "  ";
  // }
  // std::cout << "\n  ";
  // for(int j=0; j<bary_scalars.size(); j++){
  //   std::cout << bary_branches[j] << "/";
  //   std::cout << ordering_branches[bary_branches[j]] << "  ";
  // }
  // // std::cout << "\n  ";
  // // for(int j=0; j<bary_scalars.size(); j++){
  // //   std::cout << std::setprecision(2) << bary_scalars[j] << "/";
  // //   std::cout << memiOrdering[j] << "  ";
  // // }
  // std::cout << "\n-----" << std::endl;
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
int ttkTemporalMergeTreeMap2::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {


  //---------------------------------------------------------------
  // internal base layer filter approach

  // auto inputNodes = vtkMultiBlockDataSet::GetData(inputVector[0]);
  // if(!inputNodes)
  //   return 0;
  // auto inputArcs = vtkMultiBlockDataSet::GetData(inputVector[1]);
  // if(!inputArcs)
  //   return 0;
  // auto domains = vtkMultiBlockDataSet::GetData(inputVector[2]);
  // if(!domains)
  //   return 0;

  // std::vector<ttk::ftm::MergeTree<double>> trees(inputNodes->GetNumberOfBlocks());
  // for(int i = 0; i < inputNodes->GetNumberOfBlocks(); i++) {
  //   auto treeNodes = vtkUnstructuredGrid::SafeDownCast(inputNodes->GetBlock(i));
  //   auto treeArcs = vtkUnstructuredGrid::SafeDownCast(inputArcs->GetBlock(i));
  //   trees[i] = ttk::ftm::makeTree<double>(treeNodes, treeArcs);
  //   std::cout << trees[i].tree.getNumberOfNodes() << std::endl;
  // }

  // ttk::MergeTreeBarycenter mergeTreeBarycenter;
  // mergeTreeBarycenter.setThreadNumber(this->threadNumber_);
  // mergeTreeBarycenter.setDebugLevel(this->debugLevel_);

  // std::vector<std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, double>>> matchingBary;
  // std::vector<std::vector<std::pair<std::pair<ttk::ftm::idNode, ttk::ftm::idNode>, std::pair<ttk::ftm::idNode, ttk::ftm::idNode>>>> matchingPath;
  // ttk::ftm::MergeTree<double> baryMT;

  // mergeTreeBarycenter.execute<double>(
  //   trees, matchingBary, matchingPath, baryMT);
  // // trees1NodeCorrMesh = mergeTreeBarycenter.getTreesNodeCorr();
  // // finalDistances = mergeTreeBarycenter.getFinalDistances();

  // std::cout << trees.size() << baryMT.tree.getNumberOfNodes() << std::endl;

  //------------------------------------------------------------
  // internal vtk filter approach
  ttk::Timer completeTimer;
  auto inputNodes = vtkMultiBlockDataSet::GetData(inputVector[0]);
  if(!inputNodes)
    return 0;
  auto inputArcs = vtkMultiBlockDataSet::GetData(inputVector[1]);
  if(!inputArcs)
    return 0;
  auto domains = vtkMultiBlockDataSet::GetData(inputVector[2]);
  if(!domains)
    return 0;

  bool isJoinTree = false;
  std::cout << inputNodes->GetNumberOfBlocks() << std::endl;
  if(inputNodes->GetNumberOfBlocks()>0){
    auto inputNodes0 = vtkUnstructuredGrid::SafeDownCast(inputNodes->GetBlock(0));
    isJoinTree = inputNodes0->GetPointData()->GetArray("Scalar")->GetComponent(0,0) < inputNodes0->GetPointData()->GetArray("Scalar")->GetComponent(1,0);
  }
  std::cout << "isJoinTree: " << isJoinTree << std::endl;

  // bool sliding_window = true;
  // int windowSize = 5;
  std::vector<double> ordering_branches;
  vtkNew<vtkMultiBlockDataSet> members;

  if(!this->useSlidingWindow){
    ttk::Timer clusteringTimer;
    this->printMsg("Computing global barycenter", 0,0);
    vtkNew<vtkMultiBlockDataSet> mtmb;
    mtmb->SetNumberOfBlocks(2);
    mtmb->SetBlock(0,inputNodes);
    mtmb->SetBlock(1,inputArcs);

    computeBaryBranchOrdering(mtmb.GetPointer(),members.GetPointer(),ordering_branches);
    this->printMsg("Computed barycenter", 0.5,clusteringTimer.getElapsedTime());
  }

  auto memberNodes = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(0));
  auto memberArcs = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(1));

  // std::cout << members->GetNumberOfBlocks() << " " << barycenter->GetNumberOfBlocks() << std::endl;

  // return 1;

  //---------------------------------------------------------
  // barycenter input approach

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  // auto members = vtkMultiBlockDataSet::GetData(inputVector[0]);
  // if(!members)
  //   return 0;
  // auto barycenter = vtkMultiBlockDataSet::GetData(inputVector[1]);
  // if(!barycenter)
  //   return 0;
  // auto domains = vtkMultiBlockDataSet::GetData(inputVector[2]);
  // if(!domains)
  //   return 0;

  std::vector<std::vector<double>> linearizations;
  std::vector<std::vector<int>> segmentations;
  std::vector<std::vector<int>> barycenterRefs;
  int maxlen = 0;

  for(int blockIdx=0; blockIdx < inputNodes->GetNumberOfBlocks(); blockIdx++){
    int memberIdx = blockIdx;
    if(this->useSlidingWindow){
      ttk::Timer clusteringTimer;
      vtkNew<vtkMultiBlockDataSet> mtmb;
      mtmb->SetNumberOfBlocks(2);
      mtmb->SetBlock(0,vtkNew<vtkMultiBlockDataSet>());
      mtmb->SetBlock(1,vtkNew<vtkMultiBlockDataSet>());
      int sb = std::max(blockIdx-this->windowSize,0);
      memberIdx = this->windowSize;
      if(blockIdx-this->windowSize<0) memberIdx += blockIdx-this->windowSize;
      int eb = std::min(blockIdx+this->windowSize,(int)inputNodes->GetNumberOfBlocks()-1);
      this->printMsg("Computing barycenter from blocks " + std::to_string(sb) + " to " + std::to_string(eb) , 0.5 + 0.5*(blockIdx)/(inputNodes->GetNumberOfBlocks()),0);
      vtkMultiBlockDataSet::SafeDownCast(mtmb->GetBlock(0))->SetNumberOfBlocks(eb-sb+1);
      vtkMultiBlockDataSet::SafeDownCast(mtmb->GetBlock(1))->SetNumberOfBlocks(eb-sb+1);
      int bidx = 0;
      for(int b=sb; b<=eb; b++){
        vtkMultiBlockDataSet::SafeDownCast(mtmb->GetBlock(0))->SetBlock(bidx,inputNodes->GetBlock(b));
        vtkMultiBlockDataSet::SafeDownCast(mtmb->GetBlock(1))->SetBlock(bidx,inputArcs->GetBlock(b));
        bidx++;
      }
      computeBaryBranchOrdering(mtmb.GetPointer(),members.GetPointer(),ordering_branches);
      this->printMsg("Computed barycenter from blocks " + std::to_string(sb) + " to " + std::to_string(eb) , 0.5 + 0.5*(blockIdx+1)/(inputNodes->GetNumberOfBlocks()),0);

      memberNodes = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(0));
      memberArcs = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(1));
    }

    int totalSize = 0;
    std::vector<int> arcRegions;
    auto memiNodes = vtkUnstructuredGrid::SafeDownCast(memberNodes->GetBlock(memberIdx));
    auto memiArcs = vtkUnstructuredGrid::SafeDownCast(memberArcs->GetBlock(memberIdx));
    auto memiDomain = vtkDataSet::SafeDownCast(domains->GetBlock(blockIdx));

    auto scalarArrayDomain = this->GetInputArrayToProcess(0, memiDomain);
    // # print(" ")

    // prepare scalar values in segments
    std::vector<std::vector<double>> memiSegmentScalars(memiArcs->GetNumberOfCells());
    for(int j=0; j<memiDomain->GetNumberOfPoints(); j++){
      auto scalar = scalarArrayDomain->GetComponent(j,0);
      auto segId = (memiDomain->GetPointData()->GetArray("SegmentationId"))->GetComponent(j,0);
      memiSegmentScalars[segId].push_back(scalar);
    }
    for(int i=0; i<memiSegmentScalars.size(); i++){
      auto l = memiSegmentScalars[i];
      // sort in descending order if we're a join tree, ascending if we're a split tree
      if(isJoinTree){
        std::sort(l.begin(),l.end(),std::greater<double>());
      }
      else{
        std::sort(l.begin(),l.end());
      }
    }

    // get node properties of member tree
    int numNodesi = 0;
    int nnmti = vtkUnstructuredGrid::SafeDownCast(inputArcs->GetBlock(blockIdx))->GetNumberOfPoints();
    std::vector<int> memiNodeIsDummy(memiNodes->GetNumberOfPoints());
    std::vector<double> memiScalars(memiNodes->GetNumberOfPoints());
    for(int i=0; i<memiNodes->GetNumberOfPoints(); i++){
      auto nId = (memiNodes->GetPointData()->GetArray("NodeId"))->GetComponent(i,0);
      auto isDummy = (memiNodes->GetPointData()->GetArray("isDummyNode"))->GetComponent(i,0);
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
      auto childId = (memiArcs->GetCellData()->GetArray("downNodeId"))->GetComponent(i,0);
      auto parentId = (memiArcs->GetCellData()->GetArray("upNodeId"))->GetComponent(i,0);
      memiParents[childId] = parentId;
    }

    //  create tree structure of member tree (children lists)
    std::vector<std::vector<int>> memiChildren(memiNodes->GetNumberOfPoints());
    int root = -1;
    int maxDegree = 0;
    for(int i=0; i<memiParents.size(); i++){
      auto parentId = memiParents[i];
      if(parentId >= 0){
        memiChildren[parentId].push_back(i);
        maxDegree = std::max(maxDegree,(int)memiChildren[parentId].size());
        // std::cout << i << "/" << memiChildren[parentId].size() << " ; ";
      }
      if(parentId==-1 && memiNodeIsDummy[i])
        root = i;
    }

    //std::cout << blockIdx << ": " << maxDegree << std::endl;
    //if(maxDegree>2){
    //  std::cout << "  !!" << maxDegree << std::endl;
    //}

    // get barycenter branchIDs of member tree nodes
    std::vector<int> branchNodeIDs(memiNodes->GetNumberOfPoints(),-1);
    for(int j=0; j<memiNodes->GetNumberOfPoints(); j++){
      auto nId = (memiNodes->GetPointData()->GetArray("NodeId"))->GetComponent(j,0);
      auto bId = (memiNodes->GetPointData()->GetArray("BranchBaryNodeID"))->GetComponent(j,0);
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
        // std:: cout << memiOrdering[node] << std::endl;
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
        // std:: cout << memiOrdering[node] << std::endl;
        // return oidx;
      }
    }

    // # get area sizes and segmentation ids for member tree nodes (from parent edge)
    std::vector<int> memiSizes(memiNodes->GetNumberOfPoints());
    std::vector<int> memiSegs(memiNodes->GetNumberOfPoints());
    for(int j=0; j<memiArcs->GetNumberOfCells(); j++){
      auto rS = memiArcs->GetCellData()->GetArray("RegionSize")->GetComponent(j,0);
      auto isDummy = (memiArcs->GetCellData()->GetArray("isDummyArc"))->GetComponent(j,0);
      auto downNodeId = (memiArcs->GetCellData()->GetArray("downNodeId"))->GetComponent(j,0);
      int segId = memiArcs->GetCellData()->GetArray("SegmentationId")->GetComponent(j,0);
      memiSizes[downNodeId] = int(rS);
      memiSegs[downNodeId] = int(segId);
    }

    // std::cout << blockIdx << "-----\n  ";
    // for(int j=0; j<memiScalars.size(); j++){
    //   std::cout << j << ":";
    //   std::cout << memiParents[j] << "  ";
    // }
    // std::cout << "\n  ";
    // for(int j=0; j<memiScalars.size(); j++){
    //   std::cout << branchNodeIDs[j] << "/";
    //   std::cout << ordering_branches[branchNodeIDs[j]] << "  ";
    // }
    // std::cout << "\n  ";
    // for(int j=0; j<memiScalars.size(); j++){
    //   std::cout << std::setprecision(2) << memiScalars[j] << "/";
    //   std::cout << memiOrdering[j] << "  ";
    // }
    // std::cout << "\n-----" << std::endl;

    // # compute linearization based on odering
    std::vector<double> linearization;
    std::vector<int> segmentation;
    std::vector<int> barycenterRef;
    dfs_linearization(memiChildren[root][0],linearization,segmentation,barycenterRef,memiChildren,memiSegmentScalars,memiSizes,memiSegs,branchNodeIDs,memiScalars,memiOrdering);
    linearizations.push_back(linearization);
    segmentations.push_back(segmentation);
    barycenterRefs.push_back(barycenterRef);
    // print(len(linearization),len(segmentation))
    if (linearization.size() > maxlen){
      maxlen = linearization.size();
    }
    // std::cout << linearization.size() << " " << memiDomain->GetNumberOfPoints() << " " << memiNodes->GetNumberOfPoints() << std::endl;
  }
  // std::cout << maxlen << std::endl;

  // write linearization to output vti

  auto outputmb = vtkMultiBlockDataSet::GetData(outputVector,0);
  outputmb->SetNumberOfBlocks(1);
  vtkNew<vtkImageData> tmtm;
  tmtm->SetDimensions(linearizations.size()+1,maxlen+1,1);
  tmtm->SetSpacing(std::ceil(maxlen/linearizations.size())*2,1,1);
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
  this->printMsg("Computed temporal merge tree map", 1, completeTimer.getElapsedTime());
  return 1;
}
