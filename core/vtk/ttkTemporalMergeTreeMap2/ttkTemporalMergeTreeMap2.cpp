#include <ttkTemporalMergeTreeMap2.h>

#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <Debug.h>
#include <vtkCellData.h>

#include <FTMTreeUtils.h>
#include <MergeTreeBarycenter.h>
#include <ttkMacros.h>
#include <ttkMergeTreeClustering.h>
#include <ttkMergeTreeFeatureTracking.h>
#include <ttkMergeTreeUtils.h>
#include <ttkUtils.h>

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
int ttkTemporalMergeTreeMap2::FillInputPortInformation(int port,
                                                       vtkInformation *info) {
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
int ttkTemporalMergeTreeMap2::FillOutputPortInformation(int port,
                                                        vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

void ttkTemporalMergeTreeMap2::dfs_linearization(
  ttk::SimplexId curr_node,
  std::vector<double> &lin,
  std::vector<ttk::SimplexId> &seg,
  std::vector<ttk::SimplexId> &bar,
  std::vector<double> &nodePositions,
  std::vector<std::vector<ttk::SimplexId>> &memiChildren,
  std::vector<std::vector<double>> &memiSegmentScalars,
  std::vector<ttk::SimplexId> &memiSizes,
  std::vector<ttk::SimplexId> &memiSegs,
  std::vector<ttk::SimplexId> &branchNodeIDs,
  std::vector<double> &memiScalars,
  std::vector<double> &memiOrdering,
  std::vector<ttk::SimplexId> prevMatching,
  std::vector<double> prevOrdering) {
  if(memiChildren[curr_node].size() == 0) {
    auto curr_size = memiSizes[curr_node];
    auto segment = memiSegmentScalars[memiSegs[curr_node]];
    for(ttk::SimplexId s = 0; s < segment.size(); s += 2) {
      lin.push_back(segment[s]);
      seg.push_back(memiSegs[curr_node]);
      bar.push_back(branchNodeIDs.empty() ? -1 : branchNodeIDs[curr_node]);
    }
    nodePositions[curr_node] = lin.size();
    for(ttk::SimplexId s = segment.size() - 1 - (segment.size() % 2); s >= 0;
        s -= 2) {
      // for(ttk::SimplexId s=1; s<segment.size(); s+=2){
      lin.push_back(segment[s]);
      seg.push_back(memiSegs[curr_node]);
      bar.push_back(branchNodeIDs.empty() ? -1 : branchNodeIDs[curr_node]);
    }
  } else {
    std::vector<ttk::SimplexId> curr_children = memiChildren[curr_node];
    if(!memiOrdering.empty()){
      std::sort(curr_children.begin(), curr_children.end(),
                [memiOrdering](ttk::SimplexId x1, ttk::SimplexId x2) -> bool {
                  return memiOrdering[x1] < memiOrdering[x2];
                });
      std::cout << "sorting by barycenter layout" << std::endl;
    }
    else if(!prevMatching.empty()){
      std::sort(curr_children.begin(), curr_children.end(),
                [prevMatching,prevOrdering](ttk::SimplexId x1, ttk::SimplexId x2) -> bool {
                  if(prevMatching[x1]<0) return true;
                  if(prevMatching[x2]<0) return false;
                  return prevOrdering[prevMatching[x1]] < prevOrdering[prevMatching[x2]];
                });
      std::cout << "sorting by previous step" << std::endl;
    }
    else{
      std::cout << "no sorting at all" << std::endl;
    }
    // std::vector<ttk::SimplexId> curr_children;
    auto curr_size = ttk::SimplexId(memiSizes[curr_node]);
    auto segment = memiSegmentScalars[memiSegs[curr_node]];
    auto cs1 = ttk::SimplexId(curr_size / 2);
    auto cs2 = curr_size - cs1;
    for(ttk::SimplexId s = 0; s < segment.size() - 1; s += 2) {
      lin.push_back(segment[s]);
      seg.push_back(memiSegs[curr_node]);
      bar.push_back(branchNodeIDs.empty() ? -1 : branchNodeIDs[curr_node]);
    }
    if(curr_children.size() > 2) {
      std::cout << "!! " << curr_children.size() << std::endl;
    }
    for(ttk::SimplexId ci = 0; ci < curr_children.size(); ci++) {
      auto c = curr_children[ci];
      dfs_linearization(c, lin, seg, bar, nodePositions, memiChildren, memiSegmentScalars,
                        memiSizes, memiSegs, branchNodeIDs, memiScalars,
                        memiOrdering,prevMatching,prevOrdering);
      if(ci < curr_children.size() - 1) {
        lin.push_back(memiScalars[curr_node]);
        seg.push_back(memiSegs[curr_node]);
        bar.push_back(branchNodeIDs.empty() ? -1 : branchNodeIDs[curr_node]);
      }
      if(ci == 0) {
        nodePositions[curr_node] = lin.size();
      }
    }
    for(ttk::SimplexId s = segment.size() - 2 - (1 - (segment.size() % 2));
        s >= 0; s -= 2) {
      // for(ttk::SimplexId s=1; s<segment.size()-1; s+=2){
      lin.push_back(segment[s]);
      seg.push_back(memiSegs[curr_node]);
      bar.push_back(branchNodeIDs.empty() ? -1 : branchNodeIDs[curr_node]);
    }
  }
}

void ttkTemporalMergeTreeMap2::computeBaryBranchOrdering(
  vtkMultiBlockDataSet *mtmb,
  vtkMultiBlockDataSet *members,
  std::vector<double> &ordering_branches) {

  vtkNew<ttkMergeTreeClustering> c;
  c->SetInputDataObject(0, mtmb);
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

  members->DeepCopy(
    vtkMultiBlockDataSet::SafeDownCast(c->GetOutputDataObject(0)));
  auto barycenter
    = vtkMultiBlockDataSet::SafeDownCast(c->GetOutputDataObject(1));

  auto baryNodes_ = vtkMultiBlockDataSet::SafeDownCast(barycenter->GetBlock(0));
  auto baryArcs_ = vtkMultiBlockDataSet::SafeDownCast(barycenter->GetBlock(1));
  auto baryNodes = vtkUnstructuredGrid::SafeDownCast(baryNodes_->GetBlock(0));
  auto baryArcs = vtkUnstructuredGrid::SafeDownCast(baryArcs_->GetBlock(0));

  ttk::SimplexId baryNumNodes = 0;
  std::vector<ttk::SimplexId> baryNodeIsDummy(baryNodes->GetNumberOfPoints());
  for(ttk::SimplexId i = 0; i < baryNodes->GetNumberOfPoints(); i++) {
    auto nId
      = vtkIntArray::SafeDownCast(baryNodes->GetPointData()->GetArray("NodeId"))
          ->GetValue(i);
    auto isDummy = vtkIntArray::SafeDownCast(
                     baryNodes->GetPointData()->GetArray("isDummyNode"))
                     ->GetValue(i);
    if(!isDummy) {
      baryNodeIsDummy[nId] = 1;
      baryNumNodes += 1;
    }
  }
  //std::cout << "survived the loop" << std::endl;

  std::vector<ttk::SimplexId> baryParents(baryNodes->GetNumberOfPoints(), -1);
  for(ttk::SimplexId i = 0; i < baryArcs->GetNumberOfCells(); i++) {
    auto childId = vtkIntArray::SafeDownCast(
                     baryArcs->GetCellData()->GetArray("downNodeId"))
                     ->GetValue(i);
    auto parentId
      = vtkIntArray::SafeDownCast(baryArcs->GetCellData()->GetArray("upNodeId"))
          ->GetValue(i);
    baryParents[childId] = parentId;
  }

  std::vector<std::vector<ttk::SimplexId>> baryChildren(
    baryNodes->GetNumberOfPoints());
  ttk::SimplexId root = -1;
  for(ttk::SimplexId i = 0; i < baryParents.size(); i++) {
    auto parentId = baryParents[i];
    if(parentId >= 0)
      baryChildren[parentId].push_back(i);
    if(parentId == -1 && baryNodeIsDummy[i])
      root = i;
  }

  std::stack<ttk::SimplexId> s;
  s.push(root);
  ttk::SimplexId idx = 0;
  auto ordering_baryNodes
    = std::vector<double>(baryNodes->GetNumberOfPoints(), -1);
  while(!s.empty()) {
    auto node = s.top();
    s.pop();
    ordering_baryNodes[node] = idx;
    idx++;
    for(auto child : baryChildren[node]) {
      s.push(child);
    }
  }

  // compute ordering of barycenter branchIDs
  ordering_branches = std::vector<double>(baryNodes->GetNumberOfPoints(), -1);
  auto bary_branches = std::vector<double>(baryNodes->GetNumberOfPoints(), -1);
  auto bary_scalars = std::vector<double>(baryNodes->GetNumberOfPoints(), -1);
  for(ttk::SimplexId i = 0; i < baryNodes->GetNumberOfPoints(); i++) {
    auto posX = baryNodes->GetPoint(i)[0];
    // auto posY = baryNodes->GetPoint(i)[1];
    ttk::SimplexId nId
      = baryNodes->GetPointData()->GetArray("NodeId")->GetComponent(i, 0);
    ttk::SimplexId scalar
      = baryNodes->GetPointData()->GetArray("Scalar")->GetComponent(i, 0);
    // if(!baryChildren[nId].empty()) continue;
    auto bId = vtkIntArray::SafeDownCast(
                 baryNodes->GetPointData()->GetArray("BranchNodeID"))
                 ->GetValue(i);
    ordering_branches[bId] = posX; // ordering_baryNodes[nId];
    bary_branches[nId] = bId;
    bary_scalars[nId] = scalar;
  }

  // std::cout << "bary" << "-----\n  ";
  // for(ttk::SimplexId j=0; j<bary_scalars.size(); j++){
  //   std::cout << j << ":";
  //   std::cout << baryParents[j] << "  ";
  // }
  // std::cout << "\n  ";
  // for(ttk::SimplexId j=0; j<bary_scalars.size(); j++){
  //   std::cout << bary_branches[j] << "/";
  //   std::cout << ordering_branches[bary_branches[j]] << "  ";
  // }
  // // std::cout << "\n  ";
  // // for(ttk::SimplexId j=0; j<bary_scalars.size(); j++){
  // //   std::cout << std::setprecision(2) << bary_scalars[j] << "/";
  // //   std::cout << memiOrdering[j] << "  ";
  // // }
  // std::cout << "\n-----" << std::endl;
}

// void ttkTemporalMergeTreeMap2::computeOrderings(
//   vtkMultiBlockDataSet *mtmb,
//   vtkMultiBlockDataSet *members,
//   std::vector<std::vector<double>> &orderings) {

//   vtkNew<ttkMergeTreeFeatureTracking> ft;
//   ft->SetInputDataObject(0,mtmb);
//   ft->SetImportantPairs(0);
//   ft->SetEpsilonTree1(0);
//   ft->SetEpsilon2Tree1(1);
//   ft->SetEpsilon3Tree1(1);
//   ft->SetPlanarLayout(true);
//   ft->Update();

//   members->DeepCopy(
//     vtkMultiBlockDataSet::SafeDownCast(ft->GetOutputDataObject(0)));
//   auto matchings
//     = vtkMultiBlockDataSet::SafeDownCast(ft->GetOutputDataObject(1));
//   auto memberNodes = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(0));
//   auto memberArcs = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(1));

//   auto firstMemberNodes = vtkUnstructuredGrid::SafeDownCast(memberNodes->GetBlock(0));
//   auto firstMemberArcs = vtkUnstructuredGrid::SafeDownCast(memberArcs->GetBlock(0));

//   ttk::SimplexId firstMemberNumNodes = 0;
//   std::vector<ttk::SimplexId> firstMemberNodeIsDummy(firstMemberNodes->GetNumberOfPoints());
//   for(ttk::SimplexId i = 0; i < firstMemberNodes->GetNumberOfPoints(); i++) {
//     auto nId
//       = vtkIntArray::SafeDownCast(firstMemberNodes->GetPointData()->GetArray("NodeId"))
//           ->GetValue(i);
//     ttk::SimplexId isDummy = firstMemberNodes->GetPointData()->GetArray("isDummyNode")->GetComponent(i,0);
//     if(!isDummy) {
//       firstMemberNodeIsDummy[nId] = 1;
//       firstMemberNumNodes += 1;
//     }
//   }

//   std::vector<ttk::SimplexId> firstMemberParents(firstMemberNodes->GetNumberOfPoints(), -1);
//   for(ttk::SimplexId i = 0; i < firstMemberArcs->GetNumberOfCells(); i++) {
//     auto childId = vtkIntArray::SafeDownCast(
//                      firstMemberArcs->GetCellData()->GetArray("downNodeId"))
//                      ->GetValue(i);
//     auto parentId
//       = vtkIntArray::SafeDownCast(firstMemberArcs->GetCellData()->GetArray("upNodeId"))
//           ->GetValue(i);
//     firstMemberParents[childId] = parentId;
//   }

//   std::vector<std::vector<ttk::SimplexId>> firstMemberChildren(
//     firstMemberNodes->GetNumberOfPoints());
//   ttk::SimplexId root = -1;
//   for(ttk::SimplexId i = 0; i < firstMemberParents.size(); i++) {
//     auto parentId = firstMemberParents[i];
//     if(parentId >= 0)
//       firstMemberChildren[parentId].push_back(i);
//     if(parentId == -1 && firstMemberNodeIsDummy[i])
//       root = i;
//   }

//   std::stack<ttk::SimplexId> s;
//   s.push(root);
//   ttk::SimplexId idx = 0;
//   auto ordering_firstMemberNodes
//     = std::vector<double>(firstMemberNodes->GetNumberOfPoints(), -1);
//   while(!s.empty()) {
//     auto node = s.top();
//     s.pop();
//     ordering_firstMemberNodes[node] = idx;
//     idx++;
//     for(auto child : firstMemberChildren[node]) {
//       s.push(child);
//     }
//   }

//   auto firstMember_ordering = std::vector<double>(firstMemberNodes->GetNumberOfPoints(), -1);
//   auto firstMember_branches = std::vector<double>(firstMemberNodes->GetNumberOfPoints(), -1);
//   auto firstMember_scalars = std::vector<double>(firstMemberNodes->GetNumberOfPoints(), -1);
//   for(ttk::SimplexId i = 0; i < firstMemberNodes->GetNumberOfPoints(); i++) {
//     auto posX = firstMemberNodes->GetPoint(i)[0];
//     // auto posY = firstMemberNodes->GetPoint(i)[1];
//     ttk::SimplexId nId
//       = firstMemberNodes->GetPointData()->GetArray("NodeId")->GetComponent(i, 0);
//     ttk::SimplexId scalar
//       = firstMemberNodes->GetPointData()->GetArray("Scalar")->GetComponent(i, 0);
//     // if(!firstMemberChildren[nId].empty()) continue;
//     auto bId = vtkIntArray::SafeDownCast(
//                  firstMemberNodes->GetPointData()->GetArray("BranchNodeID"))
//                  ->GetValue(i);
//     firstMember_ordering[bId] = posX; // ordering_firstMemberNodes[nId];
//     firstMember_branches[nId] = bId;
//     firstMember_scalars[nId] = scalar;
//   }

//   auto last_branches = firstMember_branches;
//   auto last_ordering = firstMember_ordering;
//   for(ttk::SimplexId blockIdx = 0; blockIdx < members->GetNumberOfBlocks(); blockIdx++) {

//     ttk::SimplexId totalSize = 0;
//     std::vector<ttk::SimplexId> arcRegions;
//     auto memiNodes
//       = vtkUnstructuredGrid::SafeDownCast(memberNodes->GetBlock(blockIdx));
//     auto memiArcs
//       = vtkUnstructuredGrid::SafeDownCast(memberArcs->GetBlock(blockIdx));

//     auto matchingbranchesi = std::vector<ttk::SimplexId>(memiNodes->GetNumberOfPoints(),-1);
//     if(blockIdx==0){
//       for(ttk::SimplexId i = 0; i < memiNodes->GetNumberOfPoints(); i++) {
//         matchingbranchesi[i] = last_branches[i];
//       }
//     }
//     else{
//       auto matching_vtk = blockIdx == 0 ? nullptr : vtkUnstructuredGrid::SafeDownCast(matchings->GetBlock(blockIdx-1));
//       for(ttk::SimplexId cellIdx = 0; cellIdx < matching_vtk->GetNumberOfCells(); cellIdx++) {
//         ttk::SimplexId n1 = matching_vtk->GetCellData()->GetArray("mergeTree1NodeId")->GetComponent(cellIdx,0);
//         ttk::SimplexId n2 = matching_vtk->GetCellData()->GetArray("mergeTree2NodeId")->GetComponent(cellIdx,0);
//         matchingbranchesi[n2] = n1;
//       }
//     }
//     auto last_branch_positions = std::set<double>();
//     for(auto n1 : matchingbranchesi){
//       if(n1>=0) last_branch_positions.insert(last_ordering[n1]);
//     }
//     auto last_branch_positions_ordered = std::vector<double>();
//     std::copy(last_branch_positions.begin(), last_branch_positions.end(), std::back_inserter(last_branch_positions_ordered));
//     double min_diff = std::numeric_limits<double>::infinity();
//     for(size_t b=1; b<last_branch_positions_ordered.size(); b++){
//       min_diff = std::min(std::abs(last_branch_positions_ordered[b-1]-last_branch_positions_ordered[b]),min_diff);
//     }

//     // get node properties of member tree
//     ttk::SimplexId numNodesi = 0;
//     ttk::SimplexId nnmti
//       = vtkUnstructuredGrid::SafeDownCast(memberArcs->GetBlock(blockIdx))
//           ->GetNumberOfPoints();
//     std::vector<ttk::SimplexId> memiNodeIsDummy(memiNodes->GetNumberOfPoints());
//     std::vector<double> memiScalars(memiNodes->GetNumberOfPoints());
//     for(ttk::SimplexId i = 0; i < memiNodes->GetNumberOfPoints(); i++) {
//       auto nId = vtkIntArray::SafeDownCast(
//                    memiNodes->GetPointData()->GetArray("NodeId"))
//                    ->GetValue(i);
//       auto isDummy = vtkIntArray::SafeDownCast(
//                        memiNodes->GetPointData()->GetArray("isDummyNode"))
//                        ->GetValue(i);
//       auto scalar
//         = memiNodes->GetPointData()->GetArray("Scalar")->GetComponent(i, 0);
//       if(!isDummy) {
//         memiNodeIsDummy[nId] = 1;
//         memiScalars[nId] = scalar;
//         numNodesi += 1;
//       }
//     }

//     // create tree structure of member tree (parent pointers)
//     std::vector<ttk::SimplexId> memiParents(memiNodes->GetNumberOfPoints(), -1);
//     for(ttk::SimplexId i = 0; i < memiArcs->GetNumberOfCells(); i++) {
//       auto childId = vtkIntArray::SafeDownCast(
//                        memiArcs->GetCellData()->GetArray("downNodeId"))
//                        ->GetValue(i);
//       auto parentId = vtkIntArray::SafeDownCast(
//                         memiArcs->GetCellData()->GetArray("upNodeId"))
//                         ->GetValue(i);
//       memiParents[childId] = parentId;
//     }

//     //  create tree structure of member tree (children lists)
//     std::vector<std::vector<ttk::SimplexId>> memiChildren(
//       memiNodes->GetNumberOfPoints());
//     ttk::SimplexId root = -1;
//     ttk::SimplexId maxDegree = 0;
//     for(ttk::SimplexId i = 0; i < memiParents.size(); i++) {
//       auto parentId = memiParents[i];
//       if(parentId >= 0) {
//         memiChildren[parentId].push_back(i);
//         maxDegree
//           = std::max(maxDegree, (ttk::SimplexId)memiChildren[parentId].size());
//         // std::cout << i << "/" << memiChildren[parentId].size() << " ; ";
//       }
//       if(parentId == -1 && memiNodeIsDummy[i])
//         root = i;
//     }

//     std::cout << blockIdx << ": " << maxDegree << std::endl;
//     if(maxDegree > 2) {
//       std::cout << "  !!" << maxDegree << std::endl;
//     }

//     // get barycenter branchIDs of member tree nodes
//     std::vector<ttk::SimplexId> branchNodeIDs(
//       memiNodes->GetNumberOfPoints(), -1);
//     for(ttk::SimplexId j = 0; j < memiNodes->GetNumberOfPoints(); j++) {
//       auto nId = vtkIntArray::SafeDownCast(
//                    memiNodes->GetPointData()->GetArray("NodeId"))
//                    ->GetValue(j);
//       auto bId = vtkIntArray::SafeDownCast(
//                    memiNodes->GetPointData()->GetArray("BranchBaryNodeID"))
//                    ->GetValue(j);
//       branchNodeIDs[nId] = bId;
//     }

//     // compute ordering of member tree nodes for layout (derived from barycenter
//     // branchIDs of leaves through dfs)
//     std::vector<double> memiOrdering(memiNodes->GetNumberOfPoints(), -1);
//     std::stack<ttk::SimplexId> s1;
//     std::stack<ttk::SimplexId> s2;
//     std::vector<ttk::SimplexId> postorder;
//     std::vector<ttk::SimplexId> preorder;
//     s1.push(root);
//     while(!s1.empty()) {
//       auto node = s1.top();
//       s1.pop();
//       s2.push(node);
//       preorder.push_back(node);
//       for(auto child : memiChildren[node]) {
//         s1.push(child);
//       }
//     }
//     while(!s2.empty()) {
//       auto node = s2.top();
//       s2.pop();
//       postorder.push_back(node);
//     }
//     for(auto node : postorder) {
//       if(memiChildren[node].size() == 0) {
//         if(matchingbranchesi[node]>=0){
//           memiOrdering[node] = last_ordering[last_branches[matchingbranchesi[node]]];
//         }
//         else{
//           memiOrdering[node] =
//         }
//         // std:: cout << memiOrdering[node] << std::endl;
//         // return memiOrdering[node];
//       } else {
//         double oidx = std::numeric_limits<double>::infinity();
//         for(auto child : memiChildren[node]) {
//           auto cidx = memiOrdering[child];
//           if(cidx < oidx)
//             oidx = cidx;
//         }
//         memiOrdering[node] = oidx;
//         // std:: cout << memiOrdering[node] << std::endl;
//         // return oidx;
//       }
//     }

//     // # get area sizes and segmentation ids for member tree nodes (from parent
//     // edge)
//     std::vector<ttk::SimplexId> memiSizes(memiNodes->GetNumberOfPoints());
//     std::vector<ttk::SimplexId> memiSegs(memiNodes->GetNumberOfPoints());
//     for(ttk::SimplexId j = 0; j < memiArcs->GetNumberOfCells(); j++) {
//       auto rS
//         = memiArcs->GetCellData()->GetArray("RegionSize")->GetComponent(j, 0);
//       auto isDummy = vtkIntArray::SafeDownCast(
//                        memiArcs->GetCellData()->GetArray("isDummyArc"))
//                        ->GetValue(j);
//       auto downNodeId = vtkIntArray::SafeDownCast(
//                           memiArcs->GetCellData()->GetArray("downNodeId"))
//                           ->GetValue(j);
//       ttk::SimplexId segId = memiArcs->GetCellData()
//                                ->GetArray("SegmentationId")
//                                ->GetComponent(j, 0);
//       memiSizes[downNodeId] = ttk::SimplexId(rS);
//       memiSegs[downNodeId] = ttk::SimplexId(segId);
//     }

//     // std::cout << blockIdx << "-----\n  ";
//     // for(ttk::SimplexId j=0; j<memiScalars.size(); j++){
//     //   std::cout << j << ":";
//     //   std::cout << memiParents[j] << "  ";
//     // }
//     // std::cout << "\n  ";
//     // for(ttk::SimplexId j=0; j<memiScalars.size(); j++){
//     //   std::cout << branchNodeIDs[j] << "/";
//     //   std::cout << ordering_branches[branchNodeIDs[j]] << "  ";
//     // }
//     // std::cout << "\n  ";
//     // for(ttk::SimplexId j=0; j<memiScalars.size(); j++){
//     //   std::cout << std::setprecision(2) << memiScalars[j] << "/";
//     //   std::cout << memiOrdering[j] << "  ";
//     // }
//     // std::cout << "\n-----" << std::endl;

//     // # compute linearization based on odering
//     std::vector<double> linearization;
//     std::vector<ttk::SimplexId> segmentation;
//     std::vector<ttk::SimplexId> barycenterRef;
//     dfs_linearization(memiChildren[root][0], linearization, segmentation,
//                       barycenterRef, memiChildren, memiSegmentScalars,
//                       memiSizes, memiSegs, branchNodeIDs, memiScalars,
//                       memiOrdering);
//     linearizations.push_back(linearization);
//     segmentations.push_back(segmentation);
//     barycenterRefs.push_back(barycenterRef);
//     // print(len(linearization),len(segmentation))
//     if(linearization.size() > maxlen) {
//       maxlen = linearization.size();
//     }
//     // std::cout << linearization.size() << " " <<
//     // memiDomain->GetNumberOfPoints() << " " << memiNodes->GetNumberOfPoints()
//     // << std::endl;
//   }
// }

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
  // for(ttk::SimplexId i = 0; i < inputNodes->GetNumberOfBlocks(); i++) {
  //   auto treeNodes = vtkUnstructuredGrid::SafeDownCast(inputNodes->GetBlock(i));
  //   auto treeArcs = vtkUnstructuredGrid::SafeDownCast(inputArcs->GetBlock(i));
  //   trees[i] = ttk::ftm::makeTree<double>(treeNodes, treeArcs);
  //   std::cout << trees[i].tree.getNumberOfNodes() << std::endl;

  // std::vector<ttk::ftm::MergeTree<double>>
  // trees(inputNodes->GetNumberOfBlocks()); for(ttk::SimplexId i = 0; i <
  // inputNodes->GetNumberOfBlocks(); i++) {
  //   auto treeNodes =
  //   vtkUnstructuredGrid::SafeDownCast(inputNodes->GetBlock(i)); auto treeArcs
  //   = vtkUnstructuredGrid::SafeDownCast(inputArcs->GetBlock(i)); trees[i] =
  //   ttk::ftm::makeTree<double>(treeNodes, treeArcs); std::cout <<
  //   trees[i].tree.getNumberOfNodes() << std::endl;
  // }

  // ttk::MergeTreeBarycenter mergeTreeBarycenter;
  // mergeTreeBarycenter.setThreadNumber(this->threadNumber_);
  // mergeTreeBarycenter.setDebugLevel(this->debugLevel_);

  // std::vector<std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode,
  // double>>> matchingBary;
  // std::vector<std::vector<std::pair<std::pair<ttk::ftm::idNode,
  // ttk::ftm::idNode>, std::pair<ttk::ftm::idNode, ttk::ftm::idNode>>>>
  // matchingPath; ttk::ftm::MergeTree<double> baryMT;

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
  // ttk::SimplexId windowSize = 5;
  std::vector<double> ordering_branches;
  vtkNew<vtkMultiBlockDataSet> members;
  std::vector<std::vector<ttk::SimplexId>> matchings;

  if(this->layoutMode==2 && !this->useSlidingWindow) {
    ttk::Timer clusteringTimer;
    this->printMsg("Computing global barycenter", 0,0);
    vtkNew<vtkMultiBlockDataSet> mtmb;
    mtmb->SetNumberOfBlocks(2);
    mtmb->SetBlock(0, inputNodes);
    mtmb->SetBlock(1, inputArcs);
    vtkNew<ttkMergeTreeFeatureTracking> ft;
    ft->SetInputDataObject(0,mtmb.GetPointer());
    ft->SetImportantPairs(0);
    ft->SetEpsilonTree1(0);
    ft->SetEpsilon2Tree1(100);
    ft->SetEpsilon3Tree1(100);
    ft->SetPlanarLayout(true);
    ft->Update();
    this->printMsg("Computed barycenter", 0.5,clusteringTimer.getElapsedTime());
    members->DeepCopy(
      vtkMultiBlockDataSet::SafeDownCast(ft->GetOutputDataObject(0)));

    auto parents = std::vector<std::vector<ttk::SimplexId>>(inputNodes->GetNumberOfBlocks());
    auto children = std::vector<std::vector<std::vector<ttk::SimplexId>>>(inputNodes->GetNumberOfBlocks());
    auto isDummy = std::vector<std::vector<ttk::SimplexId>>(inputNodes->GetNumberOfBlocks());
    for(ttk::SimplexId blockIdx = 0; blockIdx < inputNodes->GetNumberOfBlocks();
      blockIdx++) {

      auto memberNodes = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(0));
      auto memberArcs = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(1));
      auto memiNodes
        = vtkUnstructuredGrid::SafeDownCast(memberNodes->GetBlock(blockIdx));
      auto memiArcs
        = vtkUnstructuredGrid::SafeDownCast(memberArcs->GetBlock(blockIdx));

      // get node properties of member tree
      std::vector<ttk::SimplexId> memiNodeIsDummy(memiNodes->GetNumberOfPoints());
      std::vector<double> memiScalars(memiNodes->GetNumberOfPoints());
      for(ttk::SimplexId i = 0; i < memiNodes->GetNumberOfPoints(); i++) {
        auto nId = vtkIntArray::SafeDownCast(
                    memiNodes->GetPointData()->GetArray("NodeId"))
                    ->GetValue(i);
        auto isDummy = vtkIntArray::SafeDownCast(
                        memiNodes->GetPointData()->GetArray("isDummyNode"))
                        ->GetValue(i);
        auto scalar
          = memiNodes->GetPointData()->GetArray("Scalar")->GetComponent(i, 0);
        if(!isDummy) {
          memiNodeIsDummy[nId] = 1;
          memiScalars[nId] = scalar;
        }
      }

      // create tree structure of member tree (parent pointers)
      std::vector<ttk::SimplexId> memiParents(memiNodes->GetNumberOfPoints(), -1);
      for(ttk::SimplexId i = 0; i < memiArcs->GetNumberOfCells(); i++) {
        auto childId = vtkIntArray::SafeDownCast(
                        memiArcs->GetCellData()->GetArray("downNodeId"))
                        ->GetValue(i);
        auto parentId = vtkIntArray::SafeDownCast(
                          memiArcs->GetCellData()->GetArray("upNodeId"))
                          ->GetValue(i);
        memiParents[childId] = parentId;
      }

      //  create tree structure of member tree (children lists)
      std::vector<std::vector<ttk::SimplexId>> memiChildren(
        memiNodes->GetNumberOfPoints());
      ttk::SimplexId root = -1;
      ttk::SimplexId maxDegree = 0;
      for(ttk::SimplexId i = 0; i < memiParents.size(); i++) {
        auto parentId = memiParents[i];
        if(parentId >= 0) {
          memiChildren[parentId].push_back(i);
          maxDegree
            = std::max(maxDegree, (ttk::SimplexId)memiChildren[parentId].size());
          // std::cout << i << "/" << memiChildren[parentId].size() << " ; ";
        }
        if(parentId == -1 && memiNodeIsDummy[i])
          root = i;
      }
      parents[blockIdx] = memiParents;
      children[blockIdx] = memiChildren;
      isDummy[blockIdx] = memiNodeIsDummy;
    }

    auto matchings_mb = vtkMultiBlockDataSet::SafeDownCast(ft->GetOutputDataObject(1));
    matchings = std::vector<std::vector<ttk::SimplexId>>(matchings_mb->GetNumberOfBlocks());
    for(ttk::SimplexId i=0; i<matchings.size(); i++){
      auto matchingi_vtk = vtkUnstructuredGrid::SafeDownCast(matchings_mb->GetBlock(i));
      auto memberNodes = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(0));
      auto currmemberNodes = vtkUnstructuredGrid::SafeDownCast(memberNodes->GetBlock(i+1));
      auto prevmemberNodes = vtkUnstructuredGrid::SafeDownCast(memberNodes->GetBlock(i));
      matchings[i] = std::vector<ttk::SimplexId>(currmemberNodes->GetNumberOfPoints(),-1);
      for(ttk::SimplexId cellIdx = 0; cellIdx < matchingi_vtk->GetNumberOfCells(); cellIdx++) {
        ttk::SimplexId id1 = matchingi_vtk->GetCellData()->GetArray("tree1NodeId")->GetComponent(cellIdx,0);
        ttk::SimplexId id2 = matchingi_vtk->GetCellData()->GetArray("tree2NodeId")->GetComponent(cellIdx,0);
        ttk::SimplexId n1 = prevmemberNodes->GetPointData()->GetArray("NodeId")->GetComponent(id1,0);
        ttk::SimplexId n2 = currmemberNodes->GetPointData()->GetArray("NodeId")->GetComponent(id2,0);
        matchings[i][n2] = n1;
      }
    }
  }

  auto memberNodes = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(0));
  auto memberArcs = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(1));

  std::vector<std::vector<double>> linearizations;
  std::vector<std::vector<ttk::SimplexId>> segmentations;
  std::vector<std::vector<ttk::SimplexId>> barycenterRefs;
  ttk::SimplexId maxlen = 0;

  auto prevMatching = std::vector<ttk::SimplexId>();
  auto prevOrdering = std::vector<double>();
  for(ttk::SimplexId blockIdx = 0; blockIdx < inputNodes->GetNumberOfBlocks();
      blockIdx++) {

    ttk::SimplexId memberIdx = blockIdx;
    if(this->layoutMode==2 && this->useSlidingWindow) {
      ttk::Timer clusteringTimer;
      vtkNew<vtkMultiBlockDataSet> mtmb;
      mtmb->SetNumberOfBlocks(2);
      mtmb->SetBlock(0, vtkNew<vtkMultiBlockDataSet>());
      mtmb->SetBlock(1, vtkNew<vtkMultiBlockDataSet>());
      ttk::SimplexId sb = std::max((int)(blockIdx - this->windowSize), 0);
      memberIdx = this->windowSize;
      if(blockIdx - this->windowSize < 0)
        memberIdx += blockIdx - this->windowSize;
      ttk::SimplexId eb
        = std::min(blockIdx + this->windowSize,
                   (ttk::SimplexId)inputNodes->GetNumberOfBlocks() - 1);
      this->printMsg("Computing barycenter from blocks " + std::to_string(sb) + " to " + std::to_string(eb) , 0.5 + 0.5*(blockIdx)/(inputNodes->GetNumberOfBlocks()),0);
      vtkMultiBlockDataSet::SafeDownCast(mtmb->GetBlock(0))
        ->SetNumberOfBlocks(eb - sb + 1);
      vtkMultiBlockDataSet::SafeDownCast(mtmb->GetBlock(1))
        ->SetNumberOfBlocks(eb - sb + 1);
      ttk::SimplexId bidx = 0;
      for(ttk::SimplexId b = sb; b <= eb; b++) {
        vtkMultiBlockDataSet::SafeDownCast(mtmb->GetBlock(0))
          ->SetBlock(bidx, inputNodes->GetBlock(b));
        vtkMultiBlockDataSet::SafeDownCast(mtmb->GetBlock(1))
          ->SetBlock(bidx, inputArcs->GetBlock(b));
        bidx++;
      }
      computeBaryBranchOrdering(
        mtmb.GetPointer(), members.GetPointer(), ordering_branches);
      this->printMsg("Computed barycenter from blocks " + std::to_string(sb) + " to " + std::to_string(eb) , 0.5 + 0.5*(blockIdx+1)/(inputNodes->GetNumberOfBlocks()),0);

      memberNodes = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(0));
      memberArcs = vtkMultiBlockDataSet::SafeDownCast(members->GetBlock(1));
    }

    ttk::SimplexId totalSize = 0;
    std::vector<ttk::SimplexId> arcRegions;
    auto memiNodes
      = vtkUnstructuredGrid::SafeDownCast(memberNodes->GetBlock(memberIdx));
    auto memiArcs
      = vtkUnstructuredGrid::SafeDownCast(memberArcs->GetBlock(memberIdx));
    auto memiDomain = vtkDataSet::SafeDownCast(domains->GetBlock(blockIdx));

    auto scalarArrayDomain = this->GetInputArrayToProcess(0, memiDomain);
    // # print(" ")

    auto nodePositions = std::vector<double>(memiNodes->GetNumberOfPoints(),-1);
    if(this->layoutMode==1 && blockIdx>0){
      prevMatching = matchings[blockIdx-1];
    }

    // prepare scalar values in segments
    std::vector<std::vector<double>> memiSegmentScalars(
      memiArcs->GetNumberOfCells());
    for(ttk::SimplexId j = 0; j < memiDomain->GetNumberOfPoints(); j++) {
      auto scalar = scalarArrayDomain->GetComponent(j, 0);
      auto segId = vtkIntArray::SafeDownCast(
                     memiDomain->GetPointData()->GetArray("SegmentationId"))
                     ->GetValue(j);
      memiSegmentScalars[segId].push_back(scalar);
    }
    for(ttk::SimplexId i = 0; i < memiSegmentScalars.size(); i++) {
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
    ttk::SimplexId numNodesi = 0;
    ttk::SimplexId nnmti
      = vtkUnstructuredGrid::SafeDownCast(inputArcs->GetBlock(blockIdx))
          ->GetNumberOfPoints();
    std::vector<ttk::SimplexId> memiNodeIsDummy(memiNodes->GetNumberOfPoints());
    std::vector<double> memiScalars(memiNodes->GetNumberOfPoints());
    for(ttk::SimplexId i = 0; i < memiNodes->GetNumberOfPoints(); i++) {
      auto nId = vtkIntArray::SafeDownCast(
                   memiNodes->GetPointData()->GetArray("NodeId"))
                   ->GetValue(i);
      auto isDummy = vtkIntArray::SafeDownCast(
                       memiNodes->GetPointData()->GetArray("isDummyNode"))
                       ->GetValue(i);
      auto scalar
        = memiNodes->GetPointData()->GetArray("Scalar")->GetComponent(i, 0);
      if(!isDummy) {
        memiNodeIsDummy[nId] = 1;
        memiScalars[nId] = scalar;
        numNodesi += 1;
      }
    }

    // create tree structure of member tree (parent pointers)
    std::vector<ttk::SimplexId> memiParents(memiNodes->GetNumberOfPoints(), -1);
    for(ttk::SimplexId i = 0; i < memiArcs->GetNumberOfCells(); i++) {
      auto childId = vtkIntArray::SafeDownCast(
                       memiArcs->GetCellData()->GetArray("downNodeId"))
                       ->GetValue(i);
      auto parentId = vtkIntArray::SafeDownCast(
                        memiArcs->GetCellData()->GetArray("upNodeId"))
                        ->GetValue(i);
      memiParents[childId] = parentId;
    }

    //  create tree structure of member tree (children lists)
    std::vector<std::vector<ttk::SimplexId>> memiChildren(
      memiNodes->GetNumberOfPoints());
    ttk::SimplexId root = -1;
    ttk::SimplexId maxDegree = 0;
    for(ttk::SimplexId i = 0; i < memiParents.size(); i++) {
      auto parentId = memiParents[i];
      if(parentId >= 0) {
        memiChildren[parentId].push_back(i);
        maxDegree
          = std::max(maxDegree, (ttk::SimplexId)memiChildren[parentId].size());
        // std::cout << i << "/" << memiChildren[parentId].size() << " ; ";
      }

      if(parentId == -1 && memiNodeIsDummy[i])
        root = i;
    }


    // get barycenter branchIDs of member tree nodes
    std::vector<ttk::SimplexId> branchNodeIDs;
    if(this->layoutMode == 2){
      branchNodeIDs = std::vector<ttk::SimplexId>(memiNodes->GetNumberOfPoints(), -1);
      for(ttk::SimplexId j=0; j<memiNodes->GetNumberOfPoints(); j++){
        auto nId = (memiNodes->GetPointData()->GetArray("NodeId"))->GetComponent(j,0);
        auto bId = (memiNodes->GetPointData()->GetArray("BranchBaryNodeID"))->GetComponent(j,0);
        branchNodeIDs[nId] = bId;
      }
    }

    // compute ordering of member tree nodes for layout (derived from barycenter
    // branchIDs of leaves through dfs)
    std::vector<double> memiOrdering;
    if(this->layoutMode == 2){
      memiOrdering = std::vector<double>(memiNodes->GetNumberOfPoints(), -1);
    }
    if(this->layoutMode == 1 && blockIdx>0){
      memiOrdering = std::vector<double>(memiNodes->GetNumberOfPoints(), std::numeric_limits<double>::infinity());
    }
    std::stack<ttk::SimplexId> s1;
    std::stack<ttk::SimplexId> s2;
    std::vector<ttk::SimplexId> postorder;
    std::vector<ttk::SimplexId> preorder;
    s1.push(root);
    while(!s1.empty()) {
      auto node = s1.top();
      s1.pop();
      s2.push(node);
      preorder.push_back(node);
      for(auto child : memiChildren[node]) {
        s1.push(child);
      }
    }
    while(!s2.empty()) {
      auto node = s2.top();
      s2.pop();
      postorder.push_back(node);
    }
    for(auto node : postorder) {
      if(memiChildren[node].size() == 0) {
        if(this->layoutMode==2){
          memiOrdering[node] = ordering_branches[branchNodeIDs[node]];
        }
        if(this->layoutMode==1 && blockIdx>0 && prevMatching[node]>=0){
          memiOrdering[node] = prevOrdering[prevMatching[node]];
        }
        // std:: cout << memiOrdering[node] << std::endl;
        // return memiOrdering[node];
      } else {
        if(this->layoutMode==2 || blockIdx>0){
          double oidx = std::numeric_limits<double>::infinity();
          for(auto child : memiChildren[node]) {
            auto cidx = memiOrdering[child];
            if(cidx < oidx)
              oidx = cidx;
          }
          memiOrdering[node] = oidx;
          // std:: cout << memiOrdering[node] << std::endl;
          // return oidx;
        }
      }
    }

    // # get area sizes and segmentation ids for member tree nodes (from parent
    // edge)
    std::vector<ttk::SimplexId> memiSizes(memiNodes->GetNumberOfPoints());
    std::vector<ttk::SimplexId> memiSegs(memiNodes->GetNumberOfPoints());
    for(ttk::SimplexId j = 0; j < memiArcs->GetNumberOfCells(); j++) {
      auto rS
        = memiArcs->GetCellData()->GetArray("RegionSize")->GetComponent(j, 0);
      auto isDummy = vtkIntArray::SafeDownCast(
                       memiArcs->GetCellData()->GetArray("isDummyArc"))
                       ->GetValue(j);
      auto downNodeId = vtkIntArray::SafeDownCast(
                          memiArcs->GetCellData()->GetArray("downNodeId"))
                          ->GetValue(j);
      ttk::SimplexId segId = memiArcs->GetCellData()
                               ->GetArray("SegmentationId")
                               ->GetComponent(j, 0);
      memiSizes[downNodeId] = ttk::SimplexId(rS);
      memiSegs[downNodeId] = ttk::SimplexId(segId);
    }

    // std::cout << blockIdx << "-----\n  ";
    // for(ttk::SimplexId j=0; j<memiScalars.size(); j++){
    //   std::cout << j << ":";
    //   std::cout << memiParents[j] << "  ";
    // }
    // std::cout << "\n  ";
    // for(ttk::SimplexId j=0; j<memiScalars.size(); j++){
    //   std::cout << branchNodeIDs[j] << "/";
    //   std::cout << ordering_branches[branchNodeIDs[j]] << "  ";
    // }
    // std::cout << "\n  ";
    // for(ttk::SimplexId j=0; j<memiScalars.size(); j++){
    //   std::cout << std::setprecision(2) << memiScalars[j] << "/";
    //   std::cout << memiOrdering[j] << "  ";
    // }
    // std::cout << "\n-----" << std::endl;

    // # compute linearization based on odering
    std::vector<double> linearization;
    std::vector<ttk::SimplexId> segmentation;
    std::vector<ttk::SimplexId> barycenterRef;
    dfs_linearization(memiChildren[root][0], linearization, segmentation,
                      barycenterRef, nodePositions, memiChildren, memiSegmentScalars,
                      memiSizes, memiSegs, branchNodeIDs, memiScalars,
                      memiOrdering,prevMatching,prevOrdering);
    linearizations.push_back(linearization);
    segmentations.push_back(segmentation);
    barycenterRefs.push_back(barycenterRef);
    // print(len(linearization),len(segmentation))
    if(linearization.size() > maxlen) {
      maxlen = linearization.size();
    }
    prevOrdering = nodePositions;
    // std::cout << linearization.size() << " " <<
    // memiDomain->GetNumberOfPoints() << " " << memiNodes->GetNumberOfPoints()
    // << std::endl;
  }
  // std::cout << maxlen << std::endl;

  // write linearization to output vti

  auto outputmb = vtkMultiBlockDataSet::GetData(outputVector, 0);
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
  linArray->SetNumberOfTuples(linearizations.size() * maxlen);

  segArray->SetName("SegmentationId");
  segArray->SetNumberOfComponents(1);
  segArray->SetNumberOfTuples(linearizations.size() * maxlen);

  blockArray->SetName("BlockId");
  blockArray->SetNumberOfComponents(1);
  blockArray->SetNumberOfTuples(linearizations.size() * maxlen);

  barArray->SetName("BarycenterBranchID");
  barArray->SetNumberOfComponents(1);
  barArray->SetNumberOfTuples(linearizations.size() * maxlen);

  ttk::SimplexId k = 0;
  for(ttk::SimplexId j = 0; j < maxlen; j++) {
    for(ttk::SimplexId i = 0; i < linearizations.size(); i++) {
      auto v = j < linearizations[i].size() ? linearizations[i][j] : 0;
      auto s = j < linearizations[i].size() ? segmentations[i][j] : 0;
      auto b = j < linearizations[i].size() ? barycenterRefs[i][j] : 0;
      linArray->SetValue(k, v);
      segArray->SetValue(k, s);
      blockArray->SetValue(k, i);
      barArray->SetValue(k, b);
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
