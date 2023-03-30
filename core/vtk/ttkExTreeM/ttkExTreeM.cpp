#include <ttkExTreeM.h>

#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkSignedCharArray.h>
#include <vtkUnsignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <PathCompression.h>
#include <ScalarFieldCriticalPoints2.h>

const std::array<ttk::SimplexId, 64*14*3> offsetsLUT{
0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,-1,0,1,0,0,1,-1,0,0,-1,1,0,0,1,0,-1,1,1,0,1,1,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,0,-1,0,-1,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,1,0,0,1,0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,1,0,0,-1,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,-1,0,0,-1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,0,0,-1,0,1,0,0,1,0,-1,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,-1,0,1,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,-1,0,1,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,0,0,-1,0,1,0,0,1,1,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,-1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,0,0,-1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,-1,0,0,-1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,-1,1,0,0,1,0,0,-1,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,-1,1,0,0,1,0,1,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,-1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,-1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,-1,0,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
std::array<ttk::SimplexId, 64*14> offsetsLUT2;

const std::array<ttk::SimplexId,64> nNeighborsLUT{
14,10,10,0,10,6,8,0,10,8,6,0,0,0,0,0,10,6,8,0,8,4,7,0,6,4,4,0,0,0,0,0,10,8,6,0,6,4,4,0,8,7,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

struct MyImplicitTriangulation {

  ttk::SimplexId dim[3];
  ttk::SimplexId dimM1[3];
  float dimM1F[3];
  ttk::SimplexId dimXY;

  void setDimension(int* dim_){
    this->dim[0] = dim_[0];
    this->dim[1] = dim_[1];
    this->dim[2] = dim_[2];
    this->dimM1[0] = dim_[0]-1;
    this->dimM1[1] = dim_[1]-1;
    this->dimM1[2] = dim_[2]-1;

    this->dimXY = this->dim[0]*this->dim[1];
  }

  void preconditionVertexNeighbors() {
    for(int i=0,j=0; i<64*14; i++,j+=3){
      const auto& dx =  offsetsLUT[j+0];
      const auto& dy =  offsetsLUT[j+1];
      const auto& dz =  offsetsLUT[j+2];

      offsetsLUT2[i] = dx + dy*this->dim[0] + dz*this->dimXY;
    }
  }

  ttk::SimplexId getNumberOfVertices() const {
    return this->dim[0]*this->dim[1]*this->dim[2];
  }

  ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId idx) const {
    const ttk::SimplexId z = idx / this->dimXY;
    const ttk::SimplexId idx2 = idx - (z * this->dimXY);
    const ttk::SimplexId y = idx2 / this->dim[0];
    const ttk::SimplexId x = idx2 % this->dim[0];

    int key =
      (x==0?1:x==this->dimM1[0]?2:0)+
      (y==0?4:y==this->dimM1[1]?8:0)+
      (z==0?16:z==this->dimM1[2]?32:0)
    ;

    return nNeighborsLUT[key];
  }

  void getVertexNeighbor(const ttk::SimplexId idx, const ttk::SimplexId n, ttk::SimplexId& nIdx) const {
    const ttk::SimplexId z = idx / this->dimXY;
    const ttk::SimplexId idx2 = idx - (z * this->dimXY);
    const ttk::SimplexId y = idx2 / this->dim[0];
    const ttk::SimplexId x = idx2 % this->dim[0];

    int key =
      (x==0?1:x==this->dimM1[0]?2:0)+
      (y==0?4:y==this->dimM1[1]?8:0)+
      (z==0?16:z==this->dimM1[2]?32:0)
    ;

    nIdx = idx + offsetsLUT2[key*14+n];
  }
};


vtkStandardNewMacro(ttkExTreeM);

ttkExTreeM::ttkExTreeM() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(3);
}

ttkExTreeM::~ttkExTreeM() {
}

int ttkExTreeM::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkExTreeM::FillOutputPortInformation(int port, vtkInformation *info) {
  switch(port) {
    case 0:
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      return 1;
    case 2:
      info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
      return 1;
    default:
      return 0;
  }
}

template <class triangulationType>
int ttkExTreeM::getMergeTree(vtkUnstructuredGrid *outputSkeletonArcs,
                             std::vector<ExTreeM::Branch> &mergeTree,
                             const triangulationType *triangulation) {
  vtkNew<vtkUnstructuredGrid> skeletonArcs{};
  ttk::SimplexId pointIds[2];
  ttk::SimplexId pointOrders[2];
  vtkNew<vtkPoints> points{};
  vtkNew<vtkLongLongArray> data{};
  data->SetNumberOfComponents(1);
  data->SetName("Order");
  vtkNew<vtkLongLongArray> gIdArray{};
  gIdArray->SetNumberOfComponents(1);
  gIdArray->SetName("GlobalPointIds");
  float point[3];
  std::map<ttk::SimplexId, ttk::SimplexId> addedPoints;
  ttk::SimplexId currentId = 0;
  for(auto const &b : mergeTree) {
    auto &vertices = b.vertices;
    for(size_t p = 0; p < vertices.size() - 1; p++) {
      pointIds[0] = vertices[p].second;
      pointIds[1] = vertices[p + 1].second;
      pointOrders[0] = vertices[p].first;
      pointOrders[1] = vertices[p + 1].first;
      // add each point only once to the vtkPoints
      // addedPoints.insert(x).second inserts x and is true if x was not in
      // addedPoints beforehand
      if(addedPoints.insert({pointIds[0], currentId}).second) {
        // this->printMsg("point " + std::to_string(pointIds[0]));
        triangulation->getVertexPoint(
          pointIds[0], point[0], point[1], point[2]);
        points->InsertNextPoint(point);
        data->InsertNextTuple1(pointOrders[0]);
        gIdArray->InsertNextTuple1(pointIds[0]);
        currentId++;
      }
      if(addedPoints.insert({pointIds[1], currentId}).second) {
        // this->printMsg("point " + std::to_string(pointIds[1]));
        triangulation->getVertexPoint(
          pointIds[1], point[0], point[1], point[2]);
        points->InsertNextPoint(point);
        data->InsertNextTuple1(pointOrders[1]);
        gIdArray->InsertNextTuple1(pointIds[1]);
        currentId++;
      }
      // this->printMsg("Join Tree Arc: " + std::to_string(pointIds[0]) + " "
      //                + std::to_string(pointIds[1]));
      pointIds[0] = addedPoints.at(pointIds[0]);
      pointIds[1] = addedPoints.at(pointIds[1]);
      skeletonArcs->InsertNextCell(VTK_LINE, 2, pointIds);
    }
  }
  skeletonArcs->SetPoints(points);
  outputSkeletonArcs->ShallowCopy(skeletonArcs);
  outputSkeletonArcs->GetPointData()->AddArray(data);
  outputSkeletonArcs->GetPointData()->AddArray(gIdArray);

  return 1;
}

int ttkExTreeM::RequestData(vtkInformation *,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {
  // Get the input
  auto input = vtkImageData::GetData(inputVector[0]);
  if(!input) {
    this->printErr("Unable to retrieve input data object.");
    return 0;
  }
  const size_t nVertices = input->GetNumberOfPoints();

  // Get triangulation of the input object
  auto triangulation2 = ttkAlgorithm::GetTriangulation(input);
  this->printMsg("#Edges in original domain: " + std::to_string(triangulation2->getNumberOfEdges()));
  if(!triangulation2)
    return 0;
  // // Precondition triangulation
  // this->preconditionTriangulation(triangulation);
  MyImplicitTriangulation triangulation;
  int dim[3];
  input->GetDimensions(dim);
  triangulation.setDimension(dim);
  triangulation.preconditionVertexNeighbors();

  // Get input array
  auto scalarArray = this->GetInputArrayToProcess(0, inputVector);
  if(!scalarArray) {
    this->printErr("Unable to retrieve scalar array.");
    return 0;
  }

  // Order Array
  auto orderArray = this->GetOrderArray(input, 0);
  auto orderArrayData = ttkUtils::GetPointer<ttk::SimplexId>(orderArray);

  vtkNew<ttkSimplexIdTypeArray> ascendingManifold{};
  ascendingManifold->SetNumberOfComponents(1);
  ascendingManifold->SetNumberOfTuples(nVertices);
  ascendingManifold->SetName(ttk::MorseSmaleAscendingName);

  vtkNew<ttkSimplexIdTypeArray> descendingManifold{};
  descendingManifold->SetNumberOfComponents(1);
  descendingManifold->SetNumberOfTuples(nVertices);
  descendingManifold->SetName(ttk::MorseSmaleDescendingName);

  vtkNew<ttkSimplexIdTypeArray> joinSegmentationId{};
  joinSegmentationId->SetNumberOfComponents(1);
  joinSegmentationId->SetNumberOfTuples(nVertices);
  joinSegmentationId->SetName("JoinSegmentationId");

  vtkNew<ttkSimplexIdTypeArray> splitSegmentationId{};
  splitSegmentationId->SetNumberOfComponents(1);
  splitSegmentationId->SetNumberOfTuples(nVertices);
  splitSegmentationId->SetName("SplitSegmentationId");

  vtkNew<vtkUnsignedCharArray> isSplitLeaf{};
  isSplitLeaf->SetNumberOfComponents(1);
  isSplitLeaf->SetNumberOfTuples(nVertices);
  isSplitLeaf->SetName("isSplitLeaf");

  vtkNew<vtkUnsignedCharArray> isJoinLeaf{};
  isJoinLeaf->SetNumberOfComponents(1);
  isJoinLeaf->SetNumberOfTuples(nVertices);
  isJoinLeaf->SetName("isJoinLeaf");

  // compute path compression
  {
    ttk::PathCompression subModule;
    subModule.setThreadNumber(this->threadNumber_);
    subModule.setDebugLevel(this->debugLevel_);
    subModule.setComputeSegmentation(true,true,false);

    ttk::PathCompression::OutputManifold om{
      ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
      ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
      nullptr
    };

    int status = 0;
    // ttkTypeMacroT(
    //   triangulation->getType(),
    //   (status = subModule.execute<T0>(
    //     om,
    //     ttkUtils::GetPointer<const ttk::SimplexId>(orderArray),
    //     *static_cast<T0*>(triangulation->getData())
    //   ))
    // );
    status = subModule.execute<MyImplicitTriangulation>(
      om,
      ttkUtils::GetPointer<const ttk::SimplexId>(orderArray),
      triangulation
    );
    if(status!=0)
      return 0;
  }

  // compute critical points
  // type -> thread -> cp
  std::array<std::vector<ttk::SimplexId>,4> criticalPoints;
  {
    std::array<std::vector<std::vector<ttk::SimplexId>>,4> criticalPoints_;
    ttk::ScalarFieldCriticalPoints2 subModule;
    subModule.setThreadNumber(this->threadNumber_);
    subModule.setDebugLevel(this->debugLevel_);

    int status = 0;
    status = subModule.computeCritialPoints<MyImplicitTriangulation>(
      criticalPoints_,
      ttkUtils::GetPointer<ttk::SimplexId>(orderArray),
      ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
      ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
      &triangulation
    );
    if(!status)
      return 0;

    status = subModule.mergeCriticalPointVectors(criticalPoints,criticalPoints_);
    if(!status)
      return 0;
  }

  // compute joinTree
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> persistencePairsJoin{};
  std::vector<ExTreeM::Branch> mergeTreeJoin{};
  {
    int status = 0;
#pragma omp parallel for num_threads(this->threadNumber_)
    for(size_t i = 0; i < nVertices; i++) {
      orderArrayData[i] = nVertices - orderArrayData[i] - 1;
    }

    status = this->computePairs<MyImplicitTriangulation>(
      persistencePairsJoin, mergeTreeJoin,
      ttkUtils::GetPointer<ttk::SimplexId>(joinSegmentationId),
      ttkUtils::GetPointer<unsigned char>(isJoinLeaf),
      criticalPoints[3].data(), criticalPoints[1].data(),
      criticalPoints[0].data(),
      ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold),
      ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold), orderArrayData,
      &triangulation, criticalPoints[3].size(), criticalPoints[1].size(),
      criticalPoints[0].size());

    if(status != 1)
      return 0;
  }

  // compute splitTree
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>>
    persistencePairsSplit{};
  std::vector<ExTreeM::Branch> mergeTreeSplit{};
  {
    int status = 0;
#pragma omp parallel for num_threads(this->threadNumber_)
    for(size_t i = 0; i < nVertices; i++) {
      orderArrayData[i] = nVertices - orderArrayData[i] - 1;
    }

    status = this->computePairs<MyImplicitTriangulation>(
      persistencePairsSplit, mergeTreeSplit,
      ttkUtils::GetPointer<ttk::SimplexId>(splitSegmentationId),
      ttkUtils::GetPointer<unsigned char>(isSplitLeaf),
      criticalPoints[0].data(), // minima
      criticalPoints[2].data(), // 2-saddles
      criticalPoints[3].data(), // maxima
      ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
      ttkUtils::GetPointer<ttk::SimplexId>(ascendingManifold), orderArrayData,
      &triangulation, criticalPoints[0].size(), criticalPoints[2].size(),
      criticalPoints[3].size());

    if(status != 1)
      return 0;
  }

  // Finalize Output
  {
    ttk::Timer timer;
    this->printMsg(
      "Generating Output Data Objects", 0, 0, ttk::debug::LineMode::REPLACE);

    {
      auto outputMergeTreeJoin = vtkUnstructuredGrid::GetData(outputVector, 0);
      ttkTypeMacroT(triangulation2->getType(),
                    getMergeTree<T0>(outputMergeTreeJoin, mergeTreeJoin,
                                     (T0 *)triangulation2->getData()));
    }

    {
      auto outputMergeTreeSplit = vtkUnstructuredGrid::GetData(outputVector, 1);
      ttkTypeMacroT(triangulation2->getType(),
                    getMergeTree<T0>(outputMergeTreeSplit, mergeTreeSplit,
                                     (T0 *)triangulation2->getData()));
    }
    // Create segmentation output
    {
      auto segmentation = vtkDataSet::GetData(outputVector, 2);
      segmentation->ShallowCopy(input);

      auto segmentationPD = segmentation->GetPointData();
      segmentationPD->AddArray(ascendingManifold);
      segmentationPD->AddArray(descendingManifold);
      segmentationPD->AddArray(joinSegmentationId);
      segmentationPD->AddArray(splitSegmentationId);
      segmentationPD->AddArray(isSplitLeaf);
      segmentationPD->AddArray(isJoinLeaf);
    }

    this->printMsg("Generating Output Data Objects", 1, timer.getElapsedTime());
    this->printMsg(ttk::debug::Separator::L1);
  }

  return 1;
}
