
class Line {
  int n_lines_max;
  Line** next;
  double* tsd;
 public:
  int bgn;
  int end;
  Vec pos;
  Vec dir;

  Line();
  Line(int n_lines_max, int bgn, int end);
  ~Line();

  //Line(const Line& line);
  //Line& operator=(const Line& other);

  double getTSD(int n_lines);
  Line*  getNext(int n_lines);
  void setTSD(int n_lines, double rmsd);
  void setNext(int n_lines, Line* lines);

  Line* append(int n_lines, Line* lines) throw(MLException);

  void print();
};

//void free(int n_lines, Line* lines);
void copyToArray(Line* array, int n_lines, Line* lines);
