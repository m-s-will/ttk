#include <ttkPairExtrema.h>

#include <vtkInformation.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkPairExtrema);

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
ttkPairExtrema::ttkPairExtrema() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(3);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkPairExtrema::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
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
int ttkPairExtrema::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  } else if(port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

struct MyImplicitTriangulation {

  const std::array<ttk::SimplexId, 64*14*3> offsetsLUT{
  0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,-1,0,1,0,0,1,-1,0,0,-1,1,0,0,1,0,-1,1,1,0,1,1,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,0,-1,0,-1,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,1,0,0,1,0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,1,0,0,-1,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,-1,0,0,-1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,0,0,-1,0,1,0,0,1,0,-1,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,-1,0,1,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,-1,0,1,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,1,0,0,1,0,-1,0,1,0,0,1,-1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,0,0,-1,0,1,0,0,1,1,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,-1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,0,0,-1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,-1,0,0,-1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,-1,1,0,0,1,0,0,-1,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,-1,1,0,0,1,0,1,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,-1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,-1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,-1,-1,0,0,-1,1,0,-1,0,-1,0,1,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,-1,0,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  };


  const std::array<ttk::SimplexId,64> nNeighborsLUT{
  14,10,10,0,10,6,8,0,10,8,6,0,0,0,0,0,10,6,8,0,8,4,7,0,6,4,4,0,0,0,0,0,10,8,6,0,6,4,4,0,8,7,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  };

  std::array<ttk::SimplexId, 64*14> offsetsLUT2;

  ttk::SimplexId dim[3];
  ttk::SimplexId dimM1[3];
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

  ttk::SimplexId getNumberOfVertices() const {
    return this->dim[0]*this->dim[1]*this->dim[2];
  }

  // bool onBoundary(ttk::SimplexId* xyz){
  //   return x==0 || x+1==this->dim[0] || y==0 || y+1==this->dim[1] || z==0 || z+1==this->dim[2];
  // }
  inline ttk::SimplexId toIndex(const ttk::SimplexId& x, const ttk::SimplexId& y, const ttk::SimplexId& z) const {
    return (z * this->dimXY) + (y * this->dim[0]) + x;
  }
  // ttk::SimplexId toIndexSave(const ttk::SimplexId& x, const ttk::SimplexId& y, const ttk::SimplexId& z){
  //   if(x<0 || x>=this->dim[0] || y<0 || y>=this->dim[1] || z<0 || z>=this->dim[2])
  //     return -1;
  //   return (z * this->dim[0] * this->dim[1]) + (y * this->dim[0]) + x;
  // }

  inline void toXYZ(ttk::SimplexId* xyz, const ttk::SimplexId idx) const {
    xyz[2] = idx / this->dimXY;
    const auto idx2 = idx - (xyz[2] * this->dimXY);
    xyz[1] = idx2 / this->dim[0];
    xyz[0] = idx2 % this->dim[0];
  }

  inline ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId idx) const {
    ttk::SimplexId xyz[3];
    ttk::SimplexId& x = xyz[0];
    ttk::SimplexId& y = xyz[1];
    ttk::SimplexId& z = xyz[2];
    xyz[2] = idx / this->dimXY;
    const auto idx2 = idx - (xyz[2] * this->dimXY);
    xyz[1] = idx2 / this->dim[0];
    xyz[0] = idx2 % this->dim[0];

    // build case key
    int key =
      (x==0?1:x==this->dimM1[0]?2:0)+
      (y==0?4:y==this->dimM1[1]?8:0)+
      (z==0?16:z==this->dimM1[2]?32:0)
    ;

    return nNeighborsLUT[key];
  }

  inline void getVertexNeighbor(const ttk::SimplexId idx, const ttk::SimplexId n, ttk::SimplexId& nIdx) const {
    ttk::SimplexId xyz[3];
    ttk::SimplexId& x = xyz[0];
    ttk::SimplexId& y = xyz[1];
    ttk::SimplexId& z = xyz[2];
    xyz[2] = idx / this->dimXY;
    const auto idx2 = idx - (xyz[2] * this->dimXY);
    xyz[1] = idx2 / this->dim[0];
    xyz[0] = idx2 % this->dim[0];

    int key =
      (x==0?1:x==this->dimM1[0]?2:0)+
      (y==0?4:y==this->dimM1[1]?8:0)+
      (z==0?16:z==this->dimM1[2]?32:0)
    ;

    if(key==0){
      // const ttk::SimplexId n3 = n*3;
      nIdx = idx + offsetsLUT2[n];
      // nIdx = this->toIndex(x+offsetsLUT[n3],y+offsetsLUT[n3+1],z+offsetsLUT[n3+2]);
    } else {
      nIdx = idx + offsetsLUT2[key*14+n];
      // const ttk::SimplexId n3 = key*14*3+n*3;
      // nIdx = this->toIndex(x+offsetsLUT[n3],y+offsetsLUT[n3+1],z+offsetsLUT[n3+2]);
    }
  }

  void preconditionVertexNeighbors(){
    for(int i=0,j=0; i<64*14; i++,j+=3){
      const auto& dx =  offsetsLUT[j+0];
      const auto& dy =  offsetsLUT[j+1];
      const auto& dz =  offsetsLUT[j+2];

      offsetsLUT2[i] = dx + dy*this->dim[0] + dz*this->dimXY;
    }
  }
};

template <class triangulationType>
int ttkPairExtrema::getSkeletonArcs(
  vtkUnstructuredGrid *outputSkeletonArcs,
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> &persistencePairs,
  const ttk::SimplexId *order,
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
  for(auto const &p : persistencePairs) {
    pointIds[0] = p.first;
    pointIds[1] = p.second;
    pointOrders[0] = order[p.first];
    pointOrders[1] = order[p.second];
    // add each point only once to the vtkPoints
    // addedPoints.insert(x).second inserts x and is true if x was not in
    // addedPoints beforehand
    if(addedPoints.insert({pointIds[0], currentId}).second) {
      // this->printMsg("point " + std::to_string(pointIds[0]));
      triangulation->getVertexPoint(pointIds[0], point[0], point[1], point[2]);
      points->InsertNextPoint(point);
      data->InsertNextTuple1(pointOrders[0]);
      gIdArray->InsertNextTuple1(pointIds[0]);
      currentId++;
    }
    if(addedPoints.insert({pointIds[1], currentId}).second) {
      // this->printMsg("point " + std::to_string(pointIds[1]));
      triangulation->getVertexPoint(pointIds[1], point[0], point[1], point[2]);
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
  skeletonArcs->SetPoints(points);
  outputSkeletonArcs->ShallowCopy(skeletonArcs);
  outputSkeletonArcs->GetPointData()->AddArray(data);
  outputSkeletonArcs->GetPointData()->AddArray(gIdArray);

  return 1;
}

template <class triangulationType>
int ttkPairExtrema::getMergeTree(vtkUnstructuredGrid *outputSkeletonArcs,
                                 std::vector<PairExtrema::Branch> &mergeTree,
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
int ttkPairExtrema::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation

  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  vtkPointSet *inputCriticalPoints = vtkPointSet::GetData(inputVector[1]);
  if(!inputDataSet || !inputCriticalPoints)
    return 0;

  auto order = ttkAlgorithm::GetOrderArray(inputDataSet, 0);
  auto manifoldPointData = inputDataSet->GetPointData();
  auto descendingManifold
    = manifoldPointData->GetArray(ttk::MorseSmaleDescendingName);
  auto tempArray = manifoldPointData->GetArray(ttk::MorseSmaleAscendingName);
  auto criticalFieldData = inputCriticalPoints->GetFieldData();
  auto minimaIds = criticalFieldData->GetArray("cp0id");
  auto saddle2Ids = criticalFieldData->GetArray("cp2id");
  auto maximaIds = criticalFieldData->GetArray("cp3id");
  ttk::SimplexId nSaddle2 = saddle2Ids->GetNumberOfTuples();
  ttk::SimplexId nMaxima = maximaIds->GetNumberOfTuples();
  ttk::SimplexId nMinima = minimaIds->GetNumberOfTuples();

  if(!descendingManifold | !order | !tempArray | !minimaIds | !saddle2Ids
     | !maximaIds) {
    this->printErr("Unable to retrieve input arrays.");
    return 0;
  }

  // To make sure that the selected array can be processed by this filter,
  // one should also check that the array association and format is correct.

  if((descendingManifold->GetNumberOfComponents() != 1)
     | (order->GetNumberOfComponents() != 1)) {
    this->printErr("Input arrays needs to be a scalar arrays.");
    return 0;
  }

  vtkSmartPointer<vtkDataArray> segmentation
    = vtkSmartPointer<vtkDataArray>::Take(order->NewInstance());
  segmentation->SetName("SegmentationId"); // set array name
  segmentation->SetNumberOfComponents(1); // only one component per tuple
  segmentation->SetNumberOfTuples(order->GetNumberOfTuples());

  // If all checks pass then log which array is going to be processed.
  this->printMsg("Starting computation...");

  // ttk::SimplexId nCriticalPoints = criticalType->GetNumberOfTuples();
  std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> persistencePairs{};
  std::vector<PairExtrema::Branch> mergeTree{};

  // Get ttk::triangulation of the input vtkDataSet (will create one if one does
  // not exist already).
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  this->preconditionTriangulation(triangulation); // implemented in base class

  // Templatize over the different input array data types and call the base code
  int status = 0; // this integer checks if the base code returns an error
  // construct the tree
  /*ttkTypeMacroT(triangulation->getType(),
                (status = this->computePairs<T0>(
                   persistencePairs, mergeTree,
                   ttkUtils::GetPointer<ttk::SimplexId>(segmentation),
                   ttkUtils::GetPointer<ttk::SimplexId>(minimaIds),
                   ttkUtils::GetPointer<ttk::SimplexId>(saddle2Ids),
                   ttkUtils::GetPointer<ttk::SimplexId>(maximaIds),
                   ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
                   ttkUtils::GetPointer<ttk::SimplexId>(tempArray),
                   ttkUtils::GetPointer<ttk::SimplexId>(order),
                   (T0 *)triangulation->getData(), nMinima, nSaddle2, nMaxima)));*/
  MyImplicitTriangulation trian;
  int dim[3];
  ((vtkImageData*)inputDataSet)->GetDimensions(dim);
  trian.setDimension(dim);
  trian.preconditionVertexNeighbors();

  status = this->computePairs<MyImplicitTriangulation>(
      persistencePairs, mergeTree,
      ttkUtils::GetPointer<ttk::SimplexId>(segmentation),
      ttkUtils::GetPointer<ttk::SimplexId>(minimaIds),
      ttkUtils::GetPointer<ttk::SimplexId>(saddle2Ids),
      ttkUtils::GetPointer<ttk::SimplexId>(maximaIds),
      ttkUtils::GetPointer<ttk::SimplexId>(descendingManifold),
      ttkUtils::GetPointer<ttk::SimplexId>(tempArray),
      ttkUtils::GetPointer<ttk::SimplexId>(order),
      &trian, nMinima, nSaddle2, nMaxima);

  // On error cancel filter execution
  if(status != 1)
    return 0;

  //  Construct output, see ttkFTMTree.cpp
  auto outputSkeletonArcs = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputMergeTree = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputSegmentation = vtkDataSet::GetData(outputVector, 2);

  ttkTypeMacroT(
    triangulation->getType(),
    status = getSkeletonArcs<T0>(outputSkeletonArcs, persistencePairs,
                                 ttkUtils::GetPointer<ttk::SimplexId>(order),
                                 (T0 *)triangulation->getData()));

  ttkTypeMacroT(triangulation->getType(),
                status = getMergeTree<T0>(
                  outputMergeTree, mergeTree, (T0 *)triangulation->getData()));

  outputSegmentation->ShallowCopy(inputDataSet);
  outputSegmentation->GetPointData()->AddArray(segmentation);

  // return success
  return 1;
}
