#include <PathCompression.h>

ttk::PathCompression::PathCompression() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("PathCompression");
#ifdef TTK_ENABLE_MPI
  hasMPISupport_ = true;
#endif
}
