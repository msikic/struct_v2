
class LineFit {
  int min_line_len;
  int n_lines_max;
  //
  int n_points;
  Vec *points;
  //
  Line** list_cache;
  Line** seg_cache;
  //
  void checkCachedList(int n_lines, int bgn) throw(MLException);
  Line*  getCachedList(int n_lines, int bgn) throw(MLException);
  Line*  setCachedList(int n_lines, int bgn, Line* lines) throw(MLException);
  void checkCachedSeg(int bgn, int end) throw(MLException);
  Line*  getCachedSeg(int bgn, int end) throw(MLException);
  Line*  setCachedSeg(int bgn, int end, Line* line) throw(MLException);
  //
  Line* calcSeg(int bgn, int end) throw(MLException);
  Line* findSeg(int bgn, int end) throw(MLException);
  //
  Line* findList(int bgn1, int n_lines) throw (MLException);
  //
 public:
  LineFit(int min_line_len, int n_points, Vec points[]) throw(MLException);
  ~LineFit();

  Line* find(int* n_lines, int n_points, Vec points[]) throw(MLException);
};
