#include <ConnectedComponentsPC.h>

ttk::ConnectedComponentsPC::ConnectedComponentsPC() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("ConnectedComponentsPC");
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}
