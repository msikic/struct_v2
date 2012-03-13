
class MLException: public exception {
  char * msg;
 public:
  int retval;
  MLException(const char *msg);
  MLException(const char *msg, int retval);
  MLException(ostringstream& oss);
  ~MLException() throw();
  const char* what();
};
