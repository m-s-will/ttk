#include <SimplifiedMT.h>

ttk::SimplifiedMT::SimplifiedMT() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("SimplifiedMT");
  #ifdef TTK_ENABLE_MPI
      hasMPISupport_ = true;
  #endif
}
