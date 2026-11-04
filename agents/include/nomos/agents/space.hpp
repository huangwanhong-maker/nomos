// Spaces for agents: lattices with their neighbourhoods, graphs, and the plane with a grid for radius queries.
//
// - Lattice: cells (x, y) of a w × h grid numbered y w + x, bounded or a torus; the Moore neighbourhood (the eight
//   cells around) or von Neumann's (the four sharing a side), of radius 1.
// - Graph: undirected adjacency lists from a list of edges (morphe's graphs are exact objects; an agent's space
//   needs only who is next to whom).
// - Plane: points in a box, binned in square cells of a given side, so that the points within a radius of a query
//   are found by looking at the cells the radius covers (a uniform grid; Bentley and Friedman 1979, "Data
//   structures for range searching", ACM Comput. Surv. 11, 397–409).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nomos::agents {

enum class Neighbourhood { moore, von_neumann };

struct Lattice {
  std::size_t width = 0, height = 0;
  bool torus = true;
  Neighbourhood kind = Neighbourhood::moore;

  std::size_t cells() const { return width * height; }
  std::size_t at(std::size_t x, std::size_t y) const { return y * width + x; }

  std::vector<std::size_t> neighbours(std::size_t cell) const {
    const long x = static_cast<long>(cell % width), y = static_cast<long>(cell / width);
    const long W = static_cast<long>(width), H = static_cast<long>(height);
    std::vector<std::size_t> out;
    for (long dy = -1; dy <= 1; dy++)
      for (long dx = -1; dx <= 1; dx++) {
        if (dx == 0 && dy == 0) continue;
        if (kind == Neighbourhood::von_neumann && dx != 0 && dy != 0) continue;
        long nx = x + dx, ny = y + dy;
        if (torus) nx = (nx + W) % W, ny = (ny + H) % H;
        else if (nx < 0 || ny < 0 || nx >= W || ny >= H) continue;
        out.push_back(static_cast<std::size_t>(ny * W + nx));
      }
    return out;
  }
};

struct Graph {
  std::vector<std::vector<std::size_t>> adjacency;

  Graph() = default;
  Graph(std::size_t n, const std::vector<std::pair<std::size_t, std::size_t>>& edges) : adjacency(n) {
    for (const auto& [a, b] : edges) {
      if (a >= n || b >= n) throw std::invalid_argument("nomos::agents::Graph: an edge names a vertex beyond the graph");
      if (a == b) continue;
      adjacency[a].push_back(b), adjacency[b].push_back(a);
    }
    for (auto& v : adjacency) std::sort(v.begin(), v.end()), v.erase(std::unique(v.begin(), v.end()), v.end());
  }
  // The complete graph on n vertices.
  static Graph complete(std::size_t n) {
    Graph g;
    g.adjacency.resize(n);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++)
        if (i != j) g.adjacency[i].push_back(j);
    return g;
  }
  std::size_t size() const { return adjacency.size(); }
  const std::vector<std::size_t>& neighbours(std::size_t v) const { return adjacency[v]; }
};

struct Plane {
  double x0 = 0, y0 = 0, x1 = 1, y1 = 1, cell = 1;
  std::size_t nx = 1, ny = 1;
  std::vector<std::vector<std::size_t>> bins;
  std::vector<std::pair<double, double>> points;

  Plane() = default;
  Plane(double xmin, double ymin, double xmax, double ymax, double side) : x0(xmin), y0(ymin), x1(xmax), y1(ymax), cell(side) {
    if (!(side > 0) || !(xmax > xmin) || !(ymax > ymin)) throw std::invalid_argument("nomos::agents::Plane: a box and a positive cell side");
    nx = static_cast<std::size_t>(std::ceil((xmax - xmin) / side)), ny = static_cast<std::size_t>(std::ceil((ymax - ymin) / side));
    bins.assign(nx * ny, {});
  }
  std::size_t bin(double x, double y) const {
    const auto cx = static_cast<std::size_t>(std::clamp((x - x0) / cell, 0.0, static_cast<double>(nx - 1)));
    const auto cy = static_cast<std::size_t>(std::clamp((y - y0) / cell, 0.0, static_cast<double>(ny - 1)));
    return cy * nx + cx;
  }
  std::size_t insert(double x, double y) {
    points.push_back({x, y});
    bins[bin(x, y)].push_back(points.size() - 1);
    return points.size() - 1;
  }
  // The points within r of (x, y), in increasing index order.
  std::vector<std::size_t> within(double x, double y, double r) const {
    std::vector<std::size_t> out;
    const long cx0 = static_cast<long>(std::floor((x - r - x0) / cell)), cx1 = static_cast<long>(std::floor((x + r - x0) / cell));
    const long cy0 = static_cast<long>(std::floor((y - r - y0) / cell)), cy1 = static_cast<long>(std::floor((y + r - y0) / cell));
    for (long cy = std::max(0L, cy0); cy <= std::min(static_cast<long>(ny) - 1, cy1); cy++)
      for (long cx = std::max(0L, cx0); cx <= std::min(static_cast<long>(nx) - 1, cx1); cx++)
        for (const std::size_t k : bins[static_cast<std::size_t>(cy) * nx + static_cast<std::size_t>(cx)]) {
          const double dx = points[k].first - x, dy = points[k].second - y;
          if (dx * dx + dy * dy <= r * r) out.push_back(k);
        }
    std::sort(out.begin(), out.end());
    return out;
  }
};

}  // namespace nomos::agents
