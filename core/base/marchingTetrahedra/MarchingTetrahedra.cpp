#include <MarchingTetrahedra.h>

ttk::MarchingTetrahedra::MarchingTetrahedra() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("MarchingTetrahedra");
  #ifdef TTK_ENABLE_MPI
      hasMPISupport_ = true;
  #endif
}
