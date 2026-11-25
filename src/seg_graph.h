#ifndef SEG_GRAPH_H
#define SEG_GRAPH_H

namespace seg_spatial {
class GraphSearch {
  struct Impl;
  std::unique_ptr<Impl> impl;
public:
  explicit GraphSearch(Rcpp::List network);
  ~GraphSearch();
  std::vector<Neighbor> find(int focal, double threshold, bool knn,
                            const std::vector<const double*>& populations);
};
} // namespace seg_spatial
#endif
