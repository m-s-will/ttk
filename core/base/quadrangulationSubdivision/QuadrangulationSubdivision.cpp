#include <Dijkstra.h>
#include <Geometry.h>
#include <QuadrangulationSubdivision.h>
#include <array>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stack>

#define MODULE_S "[QuadrangulationSubdivision] "

ttk::SimplexId
  ttk::QuadrangulationSubdivision::findEdgeMiddle(const size_t a,
                                                  const size_t b) const {

  std::vector<float> sum(vertexDistance_[a].size());

  // euclidian barycenter of a and b
  Point edgeEuclBary = (outputPoints_[a] + outputPoints_[b]) * 0.5F;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sum.size(); ++i) {
    float m = vertexDistance_[a][i];
    float n = vertexDistance_[b][i];
    // stay on the shortest path between a and b
    sum[i] = m + n;
    if(m != std::numeric_limits<float>::infinity()
       && n != std::numeric_limits<float>::infinity()) {
      // try to get the middle of the shortest path
      sum[i] += std::abs(m - n);
    }

    // get the euclidian distance to AB
    Point curr{};
    triangulation_->getVertexPoint(i, curr.x, curr.y, curr.z);
    // try to minimize the euclidian distance to AB too
    sum[i] += Geometry::distance(&curr.x, &edgeEuclBary.x);
  }

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

ttk::SimplexId ttk::QuadrangulationSubdivision::findQuadBary(
  const std::vector<size_t> &quadVertices) const {

  std::vector<float> sum(vertexDistance_[*quadVertices.begin()].size(),
                         std::numeric_limits<float>::infinity());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < sum.size(); ++i) {

    // skip following computation if too far from any parent quad vertex
    bool skip = std::any_of(
      quadVertices.begin(), quadVertices.end(), [&](const size_t &a) {
        return vertexDistance_[a][i] == std::numeric_limits<float>::infinity();
      });

    if(skip) {
      continue;
    }

    float m = vertexDistance_[quadVertices[0]][i];
    float n = vertexDistance_[quadVertices[1]][i];
    float o = vertexDistance_[quadVertices[2]][i];
    float p = vertexDistance_[quadVertices[3]][i];

    // try to be "near" the four parent vertices
    sum[i] = m + n + o + p;

    // try to be on the diagonals intersection
    sum[i] += std::abs(m - o);
    sum[i] += std::abs(n - p);
  }

  return std::min_element(sum.begin(), sum.end()) - sum.begin();
}

int ttk::QuadrangulationSubdivision::subdivise() {

  using edgeType = std::pair<long long, long long>;
  using vertexType = std::pair<long long, Point>;
  using std::make_pair;
  std::map<edgeType, vertexType> processedEdges;

  // deep copy of coarse input quads
  auto prevQuads(outputQuads_);

  // clear input quads buffer before re-writing it
  outputQuads_.clear();

  Timer t;

  // avoid reallocation in loop, causing invalid pointers
  outputPoints_.reserve(outputPoints_.size() * 5);

  vertexDistance_.resize(outputPoints_.size());

  // get all other vertices sharing a quad
  quadNeighbors_.clear();
  quadNeighbors_.resize(outputPoints_.size());
  getQuadNeighbors(prevQuads, quadNeighbors_, true);

  // compute shortest distance from every vertex to all other that share a quad
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_.size(); ++i) {

    // skip if already computed on a coarser subdivision
    if(vertexDistance_[i].empty()) {

      // do not propagate on the whole mesh
      std::vector<SimplexId> bounds;
      for(auto &p : quadNeighbors_[i]) {
        bounds.emplace_back(nearestVertexIdentifier_[p]);
      }

      Dijkstra::shortestPath(nearestVertexIdentifier_[i], *triangulation_,
                             vertexDistance_[i], bounds);
    }
  }

  for(auto &q : prevQuads) {
    assert(q.n == 4); // magic number...

    auto i = static_cast<size_t>(q.i);
    auto j = static_cast<size_t>(q.j);
    auto k = static_cast<size_t>(q.k);
    auto l = static_cast<size_t>(q.l);

    // middles of edges
    auto ijid = findEdgeMiddle(i, j);
    auto jkid = findEdgeMiddle(j, k);
    auto klid = findEdgeMiddle(k, l);
    auto liid = findEdgeMiddle(l, i);

    Point midij{};
    triangulation_->getVertexPoint(ijid, midij.x, midij.y, midij.z);
    Point midjk{};
    triangulation_->getVertexPoint(jkid, midjk.x, midjk.y, midjk.z);
    Point midkl{};
    triangulation_->getVertexPoint(klid, midkl.x, midkl.y, midkl.z);
    Point midli{};
    triangulation_->getVertexPoint(liid, midli.x, midli.y, midli.z);

    std::vector<size_t> quadVertices{i, j, k, l};
    // barycenter TTK identifier
    auto baryid = findQuadBary(quadVertices);
    // barycenter 3D coordinates
    Point bary{};
    triangulation_->getVertexPoint(baryid, bary.x, bary.y, bary.z);

    // order edges to avoid duplicates (ij vs. ji)
    auto ij = make_pair(std::min(q.i, q.j), std::max(q.i, q.j));
    auto jk = make_pair(std::min(q.j, q.k), std::max(q.j, q.k));
    auto kl = make_pair(std::min(q.k, q.l), std::max(q.k, q.l));
    auto li = make_pair(std::min(q.l, q.i), std::max(q.l, q.i));

    auto process_edge_middle = [&](const std::pair<long long, long long> &pair,
                                   const Point &pt, const SimplexId id) {
      /* check if edge already processed by a neighbor quad */
      if(processedEdges.find(pair) == processedEdges.end()) {
        processedEdges.insert(
          make_pair(pair, make_pair(outputPoints_.size(), pt)));
        /* add new point 3d coordinates to vector of output points */
        outputPoints_.emplace_back(pt);
        /* new point is an edge middle */
        outputVertType_.emplace_back(1);
        /* store also TTK identifier of triangular mesh vertex */
        nearestVertexIdentifier_.emplace_back(id);
      }
    };

    process_edge_middle(ij, midij, ijid);
    process_edge_middle(jk, midjk, jkid);
    process_edge_middle(kl, midkl, klid);
    process_edge_middle(li, midli, liid);

    // barycenter index in outputPoints_
    auto baryIdx = static_cast<long long>(outputPoints_.size());
    outputPoints_.emplace_back(bary);
    outputVertType_.emplace_back(2);
    nearestVertexIdentifier_.emplace_back(baryid);

    // add the four new quads
    outputQuads_.emplace_back(Quad{
      4, q.i, processedEdges[ij].first, baryIdx, processedEdges[li].first});
    outputQuads_.emplace_back(Quad{
      4, q.j, processedEdges[jk].first, baryIdx, processedEdges[ij].first});
    outputQuads_.emplace_back(Quad{
      4, q.k, processedEdges[kl].first, baryIdx, processedEdges[jk].first});
    outputQuads_.emplace_back(Quad{
      4, q.l, processedEdges[li].first, baryIdx, processedEdges[kl].first});
  }

  // output subdivision level
  auto currSubd = outputSubdivision_.back() + 1;
  auto subdBeg = outputSubdivision_.size();
  outputSubdivision_.resize(outputPoints_.size());
  std::fill(
    outputSubdivision_.begin() + subdBeg, outputSubdivision_.end(), currSubd);

  {
    std::stringstream msg;
    msg << MODULE_S "Subdivised " << prevQuads.size() << " quads into "
        << outputQuads_.size() << " new quads (" << outputPoints_.size()
        << " points) in " << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  return 0;
}

ttk::QuadrangulationSubdivision::Point
  ttk::QuadrangulationSubdivision::getQuadNormal(
    const size_t a, const std::vector<Point> &inputPoints) const {
  Point quadNormal{};

  // current vertex 3d coordinates
  Point pa = inputPoints[a];

  // find all quads that have a as vertex
  std::vector<Quad> quads{};
  for(auto &q : outputQuads_) {
    auto _a = static_cast<long long>(a);
    if(q.i == _a || q.j == _a || q.k == _a || q.l == _a) {
      quads.emplace_back(q);
    }
  }

  // store for each quad around a its two direct neighbors
  std::set<std::array<size_t, 2>> couples{};

  // find couple of neighbors of a sharing a quad
  for(auto &q : quads) {
    auto _a = static_cast<long long>(a);
    std::array<size_t, 2> tmp{};
    if(q.i == _a) {
      tmp[0] = static_cast<size_t>(q.l);
      tmp[1] = static_cast<size_t>(q.j);
    } else if(q.j == _a) {
      tmp[0] = static_cast<size_t>(q.i);
      tmp[1] = static_cast<size_t>(q.k);
    } else if(q.k == _a) {
      tmp[0] = static_cast<size_t>(q.j);
      tmp[1] = static_cast<size_t>(q.l);
    } else if(q.l == _a) {
      tmp[0] = static_cast<size_t>(q.k);
      tmp[1] = static_cast<size_t>(q.i);
    }
    couples.insert(tmp);
  }

  // triangle normals around current quadrangulation vertex
  std::vector<Point> normals{};

  for(auto &t : couples) {
    Point pb = inputPoints[t[0]];
    Point pc = inputPoints[t[1]];

    // triangle normal: cross product of two edges
    Point crossP{};
    // ab, ac vectors
    Point ab = pb - pa;
    Point ac = pc - pa;
    // compute ab ^ ac
    Geometry::crossProduct(&ab.x, &ac.x, &crossP.x);
    // magnitude
    auto mag = Geometry::magnitude(&crossP.x);
    // ensure normal not null
    if(mag > powf(10, -FLT_DIG)) {
      // unitary normal vector
      normals.emplace_back(crossP / mag);
    }
  }

  // ensure normals have same direction
  for(size_t i = 1; i < normals.size(); ++i) {
    auto dotProd = Geometry::dotProduct(&normals[0].x, &normals[i].x);
    if(dotProd < 0.0F) {
      normals[i] = normals[i] * -1.0F;
    }
  }

  if(!normals.empty()) {
    // compute mean of normals
    quadNormal = std::accumulate(
      normals.begin(), normals.end(), Point{}, std::plus<Point>());
    quadNormal = quadNormal / static_cast<float>(normals.size());
  } else {
    // set error value directly in output variable...
    quadNormal.x = NAN;
  }

  return quadNormal;
}

ttk::QuadrangulationSubdivision::Point
  ttk::QuadrangulationSubdivision::findProjection(
    const size_t a,
    const std::vector<Point> &inputPoints,
    const bool lastIter) {

  // current vertex 3d coordinates
  Point pa = inputPoints[a];

  Point res{};

  // fallback to euclidian projection code if no normals
  bool doReverseProj = reverseProjection_;

  // quad normal in a
  Point normalsMean{};

  if(doReverseProj) {
    // compute mean of normals
    normalsMean = getQuadNormal(a, inputPoints);

    doReverseProj = !std::isnan(normalsMean.x);
  }

  // found a projection in one triangle
  bool success = false;
  // list of triangle IDs to test to find a potential projection
  std::stack<SimplexId> trianglesToTest;
  // list of triangle IDs already tested
  // (takes more memory to reduce computation time)
  std::vector<bool> trianglesTested(
    triangulation_->getNumberOfTriangles(), false);
  // vertex in triangle with highest barycentric coordinate
  SimplexId nearestVertex = nearestVertexIdentifier_[a];

  // number of triangles around nearest vertex
  SimplexId triangleNumber
    = triangulation_->getVertexTriangleNumber(nearestVertex);
  // init pipeline by checking in every triangle around selected vertex
  for(SimplexId j = 0; j < triangleNumber; j++) {
    SimplexId ntid;
    triangulation_->getVertexTriangle(nearestVertex, j, ntid);
    trianglesToTest.push(ntid);
  }

  while(!trianglesToTest.empty()) {
    SimplexId i = trianglesToTest.top();
    trianglesToTest.pop();

    // skip if already tested
    if(trianglesTested[i]) {
      continue;
    }

    // get triangle vertices
    std::array<SimplexId, 3> tverts{};
    triangulation_->getTriangleVertex(i, 0, tverts[0]);
    triangulation_->getTriangleVertex(i, 1, tverts[1]);
    triangulation_->getTriangleVertex(i, 2, tverts[2]);

    // get coordinates of triangle vertices
    Point pm{}, pn{}, po{};
    triangulation_->getVertexPoint(tverts[0], pm.x, pm.y, pm.z);
    triangulation_->getVertexPoint(tverts[1], pn.x, pn.y, pn.z);
    triangulation_->getVertexPoint(tverts[2], po.x, po.y, po.z);

    // triangle normal: cross product of two edges
    Point crossP{};
    // mn, mo vectors
    Point mn = pn - pm;
    Point mo = po - pm;
    // compute mn ^ mo
    Geometry::crossProduct(&mn.x, &mo.x, &crossP.x);
    // unitary normal vector
    Point normTri = crossP / Geometry::magnitude(&crossP.x);

    // compute intersection of triangle plane and line (a, normalsMean)
    if(doReverseProj) {

      auto denom = Geometry::dotProduct(&normalsMean.x, &normTri.x);

      // check if triangle plane is parallel to quad normal
      if(std::abs(denom) < powf(10, -FLT_DIG)) {
        // skip this iteration after filling pipeline
        trianglesTested[i] = true;
        // fill pipeline with neighboring triangles
        for(auto &vert : tverts) {
          auto ntr = triangulation_->getVertexTriangleNumber(vert);
          for(SimplexId j = 0; j < ntr; ++j) {
            SimplexId tid;
            triangulation_->getVertexTriangle(vert, j, tid);
            if(tid != i) {
              trianglesToTest.push(tid);
            }
          }
        }
        continue;
      }

      // use formula from Wikipedia: line-plane intersection
      auto tmp = pm - pa;
      auto alpha = Geometry::dotProduct(&tmp.x, &normTri.x) / denom;

      // intersection
      res = pa + normalsMean * alpha;

    }
    // compute euclidian projection of a in triangle plane
    else {

      auto tmp = pa - pm;
      // projection
      res = pa - normTri * Geometry::dotProduct(&normTri.x, &tmp.x);
    }

    // compute barycentric coords of projection
    std::vector<float> baryCoords;
    Geometry::computeBarycentricCoordinates(
      &pm.x, &pn.x, &po.x, &res.x, baryCoords);

    // check if projection in triangle
    bool inTriangle = true;
    for(auto &coord : baryCoords) {
      if(coord < powf(10, -FLT_DIG)) {
        inTriangle = false;
      }
      if(coord > 1 + powf(10, -FLT_DIG)) {
        inTriangle = false;
      }
    }

    // mark triangle as tested
    trianglesTested[i] = true;

    if(inTriangle) {
      success = true;
      // should we check if we have the nearest triangle?
      break;
    }

    // extrema values in baryCoords
    auto extrema = std::minmax_element(baryCoords.begin(), baryCoords.end());

    // find the nearest triangle vertices (with the highest/positive
    // values in baryCoords) from proj
    std::vector<SimplexId> vertices(2);
    vertices[0] = tverts[extrema.second - baryCoords.begin()];
    for(size_t j = 0; j < baryCoords.size(); j++) {
      if(j != static_cast<size_t>(extrema.first - baryCoords.begin())
         && j != static_cast<size_t>(extrema.second - baryCoords.begin())) {
        vertices[1] = tverts[j];
        break;
      }
    }

    // store vertex with highest barycentric coordinate
    nearestVertex = vertices[0];

    // triangles around vertices[0] and vertices[1]
    std::array<std::set<SimplexId>, 2> vertsTriangles{};

    // get triangles around vertices
    for(size_t j = 0; j < vertices.size(); ++j) {
      SimplexId tnum = triangulation_->getVertexTriangleNumber(vertices[j]);
      for(SimplexId k = 0; k < tnum; k++) {
        SimplexId tid;
        triangulation_->getVertexTriangle(vertices[j], k, tid);
        if(tid == i) {
          continue;
        }
        vertsTriangles[j].insert(tid);
      }
    }

    // triangles to test next
    std::vector<SimplexId> common_triangles;

    // look for triangles sharing the vertices with max values in baryCoords
    std::set_intersection(vertsTriangles[0].begin(), vertsTriangles[0].end(),
                          vertsTriangles[1].begin(), vertsTriangles[1].end(),
                          std::back_inserter(common_triangles));

    for(auto &ntid : common_triangles) {
      if(!trianglesTested[ntid]) {
        trianglesToTest.push(ntid);
      }
    }
  }

  size_t trChecked
    = std::count(trianglesTested.begin(), trianglesTested.end(), true);
  const size_t maxTrChecked = 100;

  if(success && trChecked > maxTrChecked) {
    success = false;
  }

  if(success) {
    // replace in nearestVertexIdentifier_ by nearest vertex in projected
    // triangle
    nearestVertexIdentifier_[a] = nearestVertex;
  }

  if(!success) {
    if(!reverseProjection_) {
      reverseProjection_ = true;
      res = findProjection(a, inputPoints, lastIter);
      reverseProjection_ = false;
    } else {
      // replace proj by the nearest vertex?
      std::vector<float> dists(vertexNumber_);
      for(SimplexId i = 0; i < vertexNumber_; ++i) {
        Point pv{};
        triangulation_->getVertexPoint(i, pv.x, pv.y, pv.z);
        dists[i] = Geometry::distance(&pa.x, &pv.x);
      }
      auto min = std::min_element(dists.begin(), dists.end()) - dists.begin();
      triangulation_->getVertexPoint(min, res.x, res.y, res.z);
      nearestVertexIdentifier_[a] = min;
    }
  }

  // fill in debug info
  if(lastIter) {
    trianglesChecked_[a] = trChecked;
    projSucceeded_[a] = success ? (doReverseProj ? 1 : 2) : 3;
  }

  return res;
}

int ttk::QuadrangulationSubdivision::project(const std::set<size_t> &filtered,
                                             const bool lastIter) {
  Timer t;

  if(lastIter) {
    trianglesChecked_.clear();
    projSucceeded_.clear();
    trianglesChecked_.resize(outputPoints_.size());
    projSucceeded_.resize(outputPoints_.size());
  }

  // outputPoints_ deep copy to avoid OpenMP data races in
  // findTriangleQuadNormal
  auto inputPoints(outputPoints_);

  // main loop
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_.size(); i++) {

    // skip computation if i in filtered
    if(filtered.find(i) != filtered.end()) {
      continue;
    }

    // replace curr in outputPoints_ by its projection
    outputPoints_[i] = findProjection(i, inputPoints, lastIter);
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Projected " << outputPoints_.size() - filtered.size()
        << " points in " << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  return 0;
}

int ttk::QuadrangulationSubdivision::getQuadNeighbors(
  const std::vector<Quad> &quads,
  std::vector<std::set<size_t>> &neighbors,
  const bool secondNeighbors) const {
  Timer t;

  for(auto &q : quads) {
    auto i = static_cast<size_t>(q.i);
    auto j = static_cast<size_t>(q.j);
    auto k = static_cast<size_t>(q.k);
    auto l = static_cast<size_t>(q.l);
    if(secondNeighbors) {
      neighbors[i].insert(j);
      neighbors[i].insert(k);
      neighbors[i].insert(l);
      neighbors[j].insert(i);
      neighbors[j].insert(k);
      neighbors[j].insert(l);
      neighbors[k].insert(i);
      neighbors[k].insert(j);
      neighbors[k].insert(l);
      neighbors[l].insert(i);
      neighbors[l].insert(j);
      neighbors[l].insert(k);
    } else {
      neighbors[i].insert(j);
      neighbors[i].insert(l);
      neighbors[k].insert(j);
      neighbors[k].insert(l);
      neighbors[j].insert(k);
      neighbors[j].insert(i);
      neighbors[l].insert(k);
      neighbors[l].insert(i);
    }
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Computed neighbors mapping of " << outputPoints_.size()
        << " points in " << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

int ttk::QuadrangulationSubdivision::relax(const std::set<size_t> &filtered) {
  Timer t;

  // outputPoints_ deep copy
  auto tmp(outputPoints_);

  // loop over output points, do not touch input MSC critical points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_.size(); i++) {

    // skip computation if i in filtered
    if(filtered.find(i) != filtered.end()) {
      continue;
    }

    // barycenter of curr neighbors
    Point relax{};
    for(auto &neigh : quadNeighbors_[i]) {
      relax = relax + tmp[neigh];
    }
    relax = relax * (1.0F / static_cast<float>(quadNeighbors_[i].size()));

    outputPoints_[i] = relax;
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Relaxed " << outputPoints_.size() - inputVertexNumber_
        << " points in " << t.getElapsedTime() << " s." << std::endl;
    dMsg(std::cout, msg.str(), detailedInfoMsg);
  }

  return 0;
}

int ttk::QuadrangulationSubdivision::findExtraordinaryVertices(
  std::set<size_t> &output) {

  // clear output
  output.clear();

  // hold input quads in a vector
  std::vector<Quad> inputQuads(inputQuadNumber_);

  std::vector<std::set<size_t>> neighbors(inputVertexNumber_);

  // use outputQuads_ here because it contains the input quadrangles before
  // subdivision wrt function call position in execute()
  getQuadNeighbors(outputQuads_, neighbors);

  const size_t NORMAL_VALENCE = 4;

  for(size_t i = 0; i < neighbors.size(); ++i) {
    if(neighbors[i].size() != NORMAL_VALENCE) {
      output.insert(i);
    }
  }

  return 0;
}

bool ttk::QuadrangulationSubdivision::checkBadProjectionTube() const {
  // Compute the minimum euclidian distance between every point. If
  // several non-neighboring points are nearer than direct neighbors,
  // there is some overlab.
  Timer t;

  std::vector<std::set<size_t>> neighbors(outputPoints_.size());
  getQuadNeighbors(outputQuads_, neighbors);
  std::vector<std::set<size_t>> secondNeighbors(outputPoints_.size());
  getQuadNeighbors(outputQuads_, secondNeighbors, true);

  std::map<std::pair<size_t, size_t>, float> quadVertsDists{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  // compute euclidian distance between every quadrangle point
  for(size_t i = 0; i < outputPoints_.size(); ++i) {
    for(size_t j = i + 1; j < outputPoints_.size(); ++j) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp critical
#endif // TTK_ENABLE_OPENMP
      quadVertsDists[std::make_pair(i, j)]
        = Geometry::distance(&outputPoints_[i].x, &outputPoints_[j].x);
    }
  }

  size_t count{};

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < outputPoints_.size(); ++i) {
    // compute minimal distance to neighbors
    float minDistNeigh{std::numeric_limits<float>::infinity()};
    for(const auto n : neighbors[i]) {
      std::pair<size_t, size_t> in{};
      if(i < n) {
        in = std::make_pair(i, n);
      } else {
        in = std::make_pair(n, i);
      }
      if(quadVertsDists[in] < minDistNeigh) {
        minDistNeigh = quadVertsDists[in];
      }
    }

    // compute minimal distance to points that aren't neighbors
    float minDistPoint{std::numeric_limits<float>::infinity()};
    for(size_t j = 0; j < outputPoints_.size(); ++j) {
      // skip itself and neighbors
      if(j == i || secondNeighbors[i].find(j) != secondNeighbors[i].end()) {
        continue;
      }
      std::pair<size_t, size_t> ij{};
      if(i < j) {
        ij = std::make_pair(i, j);
      } else {
        ij = std::make_pair(j, i);
      }
      if(quadVertsDists[ij] < minDistPoint) {
        minDistPoint = quadVertsDists[ij];
      }
    }

    if(minDistPoint < minDistNeigh) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic update
#endif // TTK_ENABLE_OPENMP
      count++;
    }
  }
  const size_t maxNearPoints{10};

  std::stringstream msg;
  msg << MODULE_S "Points nearer than direct neighbors: " << count << " (in "
      << t.getElapsedTime() << "s)" << std::endl;
  dMsg(std::cout, msg.str(), detailedInfoMsg);

  return !(count > maxNearPoints);
}

void ttk::QuadrangulationSubdivision::clearData() {
  outputQuads_.clear();
  outputPoints_.clear();
  outputValences_.clear();
  outputVertType_.clear();
  outputSubdivision_.clear();
  quadNeighbors_.clear();
  vertexDistance_.clear();
  trianglesChecked_.clear();
  projSucceeded_.clear();
}

// main routine
int ttk::QuadrangulationSubdivision::execute() {

  Timer t;

  // clear output variables
  clearData();

  // store input points (MSC critical points)
  for(size_t i = 0; i < inputVertexNumber_; i++) {
    outputPoints_.emplace_back(inputVertices_[i]);
  }

  // copy input quads into vector
  for(size_t i = 0; i < inputQuadNumber_; i++) {
    outputQuads_.emplace_back(inputQuads_[i]);
  }

  // fill outputInfos_ with input data (critical points)
  outputVertType_.resize(outputPoints_.size());
  std::fill(outputVertType_.begin(), outputVertType_.end(), 0);

  // fill outputSubdivision with input data
  outputSubdivision_.resize(outputPoints_.size());
  std::fill(outputSubdivision_.begin(), outputSubdivision_.end(), 0);

  // vertices to filter from relaxation, projection
  std::set<size_t> filtered{};
  if(!lockAllInputVertices) {
    if(lockInputExtrema) {
      // get extraordinary vertices
      findExtraordinaryVertices(filtered);
    }
  } else {
    // fill vector with all input points indices from 0 to inputVertexNumber_
    for(size_t i = 0; i < inputVertexNumber_; ++i) {
      filtered.insert(i);
    }
  }

  // main loop
  for(size_t i = 0; i < subdivisionLevel_; i++) {
    // subdivise each quadrangle by creating five new points, at the
    // center of each edge (4) and at the barycenter of the four
    // vertices (1).
    subdivise();
  }

  // retrieve mapping between every vertex and its neighbors
  quadNeighbors_.clear();
  quadNeighbors_.resize(outputPoints_.size());
  getQuadNeighbors(outputQuads_, quadNeighbors_);

  // "relax" the new points, i.e. replace it by the barycenter of its
  // four neighbors
  for(size_t i = 0; i < relaxationIterations_; i++) {
    relax(filtered);

    // project all points on the nearest triangle (except MSC critical
    // points)
    project(filtered, (i == relaxationIterations_ - 1));
  }

  // compute valence of every quadrangle vertex
  outputValences_.resize(outputPoints_.size());
  std::transform(
    quadNeighbors_.begin(), quadNeighbors_.end(), outputValences_.begin(),
    [&](const std::set<size_t> &neighbors) { return neighbors.size(); });

  if(!checkBadProjectionTube()) {
    clearData();
    return 1;
  }

  {
    std::stringstream msg;
    msg << MODULE_S "Produced " << outputQuads_.size() << " quadrangles with "
        << outputPoints_.size() << " points in " << t.getElapsedTime()
        << " s. (" << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), infoMsg);
  }

  return 0;
}
