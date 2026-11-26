#ifndef SEG_GRID_H
#define SEG_GRID_H

#include <cstdint>
#include <unordered_map>

namespace seg_spatial {

// Translating these offsets replaces a new spatial query at each grid location.
class GridSearch {
  struct Offset { int col, row; double d; };
  struct Ring { std::size_t begin, end; double d; };
  static constexpr std::size_t cache_limit = 512ULL * 1024 * 1024;
  Rcpp::IntegerVector col, row;
  Rcpp::NumericVector origin;
  int ncol, nrow;
  double width, height;
  std::vector<int> heads, next;
  std::unordered_map<std::int64_t, int> sparse_heads;
  std::vector<Offset> offsets;
  std::vector<Ring> rings;
  std::vector<int> candidate_ring;
  std::vector<double> fixed_weights;
  std::vector<std::vector<double>> band_weights;
  Rcpp::NumericVector bands;
  double power;
  int weighting;
  bool normalized;
  std::size_t weight_values = 0;

  std::int64_t cell(int c, int r) const {
    return static_cast<std::int64_t>(r) * ncol + c;
  }

  int head(int c, int r) const {
    if (c < 0 || c >= ncol || r < 0 || r >= nrow) return -1;
    const auto id = cell(c, r);
    if (!heads.empty()) return heads[id];
    const auto found = sparse_heads.find(id);
    return found == sparse_heads.end() ? -1 : found->second;
  }

  void append_ring(std::vector<Neighbor>& out, int focal, std::size_t ring) {
    const std::size_t begin = out.size();
    for (std::size_t k = rings[ring].begin; k < rings[ring].end; ++k) {
      const auto& offset = offsets[k];
      const std::int64_t c = static_cast<std::int64_t>(col[focal]) + offset.col;
      const std::int64_t r = static_cast<std::int64_t>(row[focal]) + offset.row;
      if (c < 0 || c >= ncol || r < 0 || r >= nrow) continue;
      for (int j = head(static_cast<int>(c), static_cast<int>(r)); j >= 0; j = next[j])
        out.push_back({j, rings[ring].d});
    }
    // Match the existing backends' row ordering within equidistant shells.
    std::sort(out.begin() + begin, out.end(),
      [](const Neighbor& a, const Neighbor& b) { return a.id < b.id; });
    candidate_ring.resize(out.size(), static_cast<int>(ring));
  }

public:
  GridSearch(const Rcpp::List& layout, const arma::vec& X, const arma::vec& Y,
      Rcpp::NumericVector bands_, double power_, int weighting_, bool normalized_,
      bool knn, bool fixed)
    : col(layout["col"]), row(layout["row"]), origin(layout["origin"]), next(X.n_elem, -1),
      bands(bands_), power(power_), weighting(weighting_), normalized(normalized_) {
    const Rcpp::IntegerVector dimensions = layout["dimensions"];
    const Rcpp::NumericVector spacing = layout["size"];
    ncol = dimensions[0]; nrow = dimensions[1];
    width = spacing[0]; height = spacing[1];
    const std::int64_t cells = static_cast<std::int64_t>(ncol) * nrow;
    if (cells <= std::max<std::int64_t>(1024, 16LL * X.n_elem)) heads.assign(cells, -1);
    else sparse_heads.reserve(X.n_elem);
    for (int i = static_cast<int>(X.n_elem) - 1; i >= 0; --i) {
      const auto id = cell(col[i], row[i]);
      const int prior = head(col[i], row[i]);
      if (prior >= 0 && (X[prior] != X[i] || Y[prior] != Y[i]))
        Rcpp::stop("Distinct coordinates map to the same grid cell; use 'boost' or 'brute'");
      next[i] = prior;
      if (!heads.empty()) heads[id] = i;
      else sparse_heads[id] = i;
    }
    const double largest = *std::max_element(bands.begin(), bands.end());
    const double radius = std::nextafter(largest, std::numeric_limits<double>::infinity());
    const int nc = knn ? ncol - 1 : static_cast<int>(std::min<double>(ncol - 1, std::floor(radius / width)));
    const int nr = knn ? nrow - 1 : static_cast<int>(std::min<double>(nrow - 1, std::floor(radius / height)));
    const long double bound = (2.0L * nc + 1) * (2.0L * nr + 1);
    if (bound > cache_limit / sizeof(Offset))
      Rcpp::stop("Grid offsets exceed the 512 MiB cache budget; use smaller radius bands or search = 'boost'");
    offsets.reserve(static_cast<std::size_t>(bound));
    for (int r = -nr; r <= nr; ++r) {
      if ((r + nr) % 256 == 0) Rcpp::checkUserInterrupt();
      for (int c = -nc; c <= nc; ++c) {
        const double dx = c * width, dy = r * height;
        const double d = std::sqrt(dx * dx + dy * dy);
        if (knn || d <= largest) offsets.push_back({c, r, d});
      }
    }
    std::sort(offsets.begin(), offsets.end(), [](const Offset& a, const Offset& b) {
      return a.d != b.d ? a.d < b.d : a.row != b.row ? a.row < b.row : a.col < b.col;
    });
    for (std::size_t begin = 0; begin < offsets.size();) {
      std::size_t end = begin + 1;
      while (end < offsets.size() && offsets[end].d == offsets[begin].d) ++end;
      rings.push_back({begin, end, offsets[begin].d});
      begin = end;
    }
    const std::size_t geometry_bytes = offsets.capacity() * sizeof(Offset) + rings.capacity() * sizeof(Ring);
    std::size_t required = fixed ? rings.size() : 0;
    if (!knn && !fixed) for (double band : bands)
      required += std::upper_bound(rings.begin(), rings.end(), band,
        [](double d, const Ring& ring) { return d < ring.d; }) - rings.begin();
    if (geometry_bytes > cache_limit || required > (cache_limit - geometry_bytes) / sizeof(double))
      Rcpp::stop("Grid offsets and weights exceed the 512 MiB cache budget; use fewer/smaller radius bands or search = 'boost'");
    weight_values = required;
    if (fixed) {
      fixed_weights.reserve(rings.size());
      for (const auto& ring : rings)
        fixed_weights.push_back(seg_weight(ring.d, 1.0, power, weighting));
    }
    else if (!knn) for (double band : bands) {
      const double bw = normalized && band > 0 ? band : 1.0;
      band_weights.emplace_back();
      const auto end = std::upper_bound(rings.begin(), rings.end(), band,
        [](double d, const Ring& ring) { return d < ring.d; });
      band_weights.back().reserve(end - rings.begin());
      for (const auto& ring : rings) {
        if (ring.d > band) break;
        band_weights.back().push_back(weighting == Inverse && ring.d == 0 ? 1.0 :
          seg_weight(ring.d, bw, power, weighting));
      }
    }
  }

  std::vector<Neighbor> find(int focal, double target, bool knn,
      const std::vector<const double*>& populations) {
    candidate_ring.clear();
    if (knn && std::all_of(populations.begin(), populations.end(),
        [&](const double* p) { return p[focal] >= target; })) {
      candidate_ring.push_back(0);
      return {{focal, 0}};
    }
    std::vector<Neighbor> out;
    std::vector<PopulationSum> sums(populations.size());
    for (std::size_t ring = 0; ring < rings.size(); ++ring) {
      if (ring % 4096 == 0) Rcpp::checkUserInterrupt();
      const std::size_t begin = out.size();
      append_ring(out, focal, ring);
      if (knn) {
        for (std::size_t j = begin; j < out.size(); ++j)
          for (std::size_t g = 0; g < sums.size(); ++g) sums[g].add(populations[g][out[j].id]);
        if (seg_populations_reached(sums, target)) break;
      }
    }
    if (knn) {
      auto focal_row = std::find_if(out.begin(), out.end(),
        [focal](const Neighbor& nb) { return nb.id == focal; });
      if (focal_row != out.end()) std::rotate(out.begin(), focal_row, focal_row + 1);
    }
    return out;
  }

  double fixed_weight(std::size_t candidate) const {
    return fixed_weights[candidate_ring[candidate]];
  }

  void radius_weights(NeighborhoodWeights& out, const std::vector<Neighbor>& candidates,
      int band, double nearest_positive) const {
    const auto end = std::upper_bound(candidates.begin(), candidates.end(), bands[band],
      [](double d, const Neighbor& nb) { return d < nb.d; });
    const std::size_t count = end - candidates.begin();
    out.ids.resize(count); out.weights.resize(count);
    const double bw = normalized && bands[band] > 0 ? bands[band] : 1.0;
    const double tolerance = nearest_positive <= bands[band] ? nearest_positive / 2 : bw;
    const double zero_weight = weighting == Inverse ? seg_weight(tolerance, bw, power, weighting) : 1.0;
    for (std::size_t j = 0; j < count; ++j) {
      out.ids[j] = candidates[j].id;
      out.weights[j] = weighting == Inverse && candidates[j].d == 0 ? zero_weight :
        band_weights[band][candidate_ring[j]];
    }
  }

  Rcpp::List diagnostics() const {
    return Rcpp::List::create(Rcpp::Named("cell_size") = Rcpp::NumericVector::create(width, height),
      Rcpp::Named("origin") = origin,
      Rcpp::Named("dimensions") = Rcpp::IntegerVector::create(ncol, nrow),
      Rcpp::Named("offsets") = static_cast<double>(offsets.size()),
      Rcpp::Named("distance_rings") = static_cast<double>(rings.size()),
      Rcpp::Named("cached_weights") = static_cast<double>(weight_values),
      Rcpp::Named("lookup") = heads.empty() ? "sparse" : "dense");
  }
};
}
#endif
