#include <PathCompressionLocal.h>

ttk::PathCompressionLocal::PathCompressionLocal() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PathCompressionLocal");
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}
