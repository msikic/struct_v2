
#include <exception>
#include <sstream>
#include <string>
using namespace std;

#include "mlexception.h"

static char * newMsg(string str) {
  size_t len = str.size();
  char * msg = new char[len + 1];
  str.copy(msg, len);
  msg[len] = '\0';
  return msg;
}

MLException::MLException(const char *msg) {
  this->msg = newMsg(string(msg));
  retval = 1;
}

MLException::MLException(const char *msg, int retval) {
  this->msg = newMsg(string(msg));
  this->retval = retval;
}

MLException::MLException(ostringstream& oss) {
  this->msg = newMsg(oss.str());
  retval = 1;
}

MLException::~MLException() throw() {
  delete msg;
}

const char* MLException::what() {
  return msg;
}
