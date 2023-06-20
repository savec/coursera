#include <cstddef>
#include <iterator>
#include <limits>
#include <ostream>
#include <vector>
#include <iostream>
#include <random>
#include <iomanip>
#include <map>

using namespace std;
// class Graph - implements graph ADT, using adjacency matrix representation (for sake of simplicity)
// public methods:
//    size_t get_n_vertices() const
//    size_t get_n_edges() const
//    double get_edge(size_t x, size_t y) const
//    void set_edge(size_t x, size_t y, double value)
//    void delete_edge(size_t x, size_t y)
//    bool adjacent(size_t x, size_t y) const
//    vector<size_t> get_neighbors(size_t x) const

class Graph {
public:
  Graph(size_t n_v, double density, double min_distance, double max_distance);
  // default constructor implementation
  Graph(size_t n_v=10)
  : m_n_v(n_v)
  , m_n_e(0) 
  , m_matrix(n_v, vector<double>(n_v, -1)) {}
  // copy constructor implementation
  Graph(const Graph& other)
  : m_n_v(other.m_n_v)
  , m_n_e(other.m_n_e) 
  , m_matrix(other.m_n_v, vector<double>(other.m_n_v, -1)) {
    for (size_t i = 0; i < m_n_v; ++i) {
      for (size_t j = 0; j < m_n_v; ++j) {
        m_matrix[i][j] = other.m_matrix[i][j];
      }
    }
  }
  // move constructor implementation (just in case :))
  Graph(Graph&& other) noexcept
  : m_n_v(other.m_n_v)
  , m_n_e(other.m_n_e) 
  , m_matrix(std::move(other.m_matrix)) {
    other.m_n_v = 0;
    other.m_n_e = 0;
  }
  // we don't need to specify a custom destructor for the class since the default one 
  // will invoke vectors' destructor automatically
  // ~Graph() {}

  // get number of vertices
  size_t get_n_vertices() const {
    return m_n_v;
  }
  // get number of edges
  size_t get_n_edges() const {
    return m_n_e;
  }
  // get edge value from vertex x to vertex y, negative value means there is no edge
  double get_edge(size_t x, size_t y) const {
    return m_matrix[x][y];
  }
  // set edge value from vertex x to vertex y, negative value means there is no edge
  void set_edge(size_t x, size_t y, double value) {
    m_matrix[x][y] = m_matrix[y][x] = value;
  }
  // an alias for set_edge(x, y, -1)
  void delete_edge(size_t x, size_t y) {
    set_edge(x, y, -1);
  }
  // check if there is an edge between vertices x and y 
  bool adjacent(size_t x, size_t y) const {
    return get_edge(x, y) >= 0;
  }
  // get list of neighbors of the vertex x
  vector<size_t> get_neighbors(size_t x) const {
    vector<size_t> neighbors;
    for (size_t i = 0; i < m_n_v; ++i) {
      if (m_matrix[x][i] >= 0) {
        neighbors.push_back(i);
      }
    }
    return neighbors;
  }

private:
  // output to a stream operator overload 
  friend ostream& operator<<(ostream& os, Graph const& g) {
    for (int i = 0; i < g.m_n_v; ++i) {
      for (int j = 0; j < g.m_n_v; ++j) {
        const auto v = g.m_matrix[i][j]; 
        if (v == -1.0) {
          os << "[---] ";
        } else {
          os << "[" << fixed << setprecision(1) << g.m_matrix[i][j] << "] ";
        }
      }
      os << endl;
    }
    return os;
  }

  vector<vector<double>> m_matrix;  
  size_t m_n_v;
  size_t m_n_e;
};
// the constructor creates a weighted undirected graph with n_v vertices and randomly builds 
// edges so to correspond the desired density; the weight of each edge is choosen randomly from 
// the range (min_distance, max_distance) 
Graph::Graph(size_t n_v, double density, double min_distance=1.0, double max_distance=10.0)
: m_n_v(n_v)
, m_n_e(0) 
, m_matrix(n_v, vector<double>(n_v, -1)){
  // Create a random number generator engine
  random_device rd;
  mt19937 engine(rd());

  // Define a distribution for double values between 0.0 and 1.0
  uniform_real_distribution<double> edge_distribution(0.0, 1.0);    
  uniform_real_distribution<double> dist_distribution(min_distance, max_distance);    
  for (int i = 0; i < m_n_v; ++i) {
    for (int j = 0; j < m_n_v; ++j) {
      if (i == j) {
        continue;
      }     
      if (edge_distribution(engine) < density) {
        // the graph is undirected
        m_matrix[i][j] = m_matrix[j][i] = dist_distribution(engine);
        m_n_e ++;
      }
    }
  }
}
// class PriorityQueue - a helper class, implements a priority queue
// in contrast to std::priority_queue provides methods for updating random elements
// public methods:
//    bool update(size_t vertex, double distance)
//    double distance(size_t vertex)
//    bool contains(size_t vertex) const
//    size_t size() const
//    pair<size_t, double> top() const
//    pair<size_t, double> pop()
class PriorityQueue {
public:
  PriorityQueue(size_t size) {
    for (size_t i = 0; i < size; ++i) {
      m_storage.insert({i, numeric_limits<double>::infinity()});
    }
  }
  // updates the given vertex' distance only if new distance is less than the old
  // returns true if distance was updated, false otherwise
  bool update(size_t vertex, double distance) {
    if (contains(vertex) && m_storage[vertex] > distance) {
      m_storage[vertex] = distance;
      return true;
    }
    return false;
  }
  // returns the distance of the given vertex
  double distance(size_t vertex) {
    if (contains(vertex)) {
      return m_storage[vertex];
    }
    return numeric_limits<double>::infinity();
  }
  // returns true if the vertex present in the queue, false otherwise
  bool contains(size_t vertex) const {
    auto found = m_storage.find(vertex);
    return found != m_storage.end();
  }
  // returns size of the queue
  size_t size() const {
    return m_storage.size();
  }
  // returns (using peek semantic) an element with the least distance
  pair<size_t, double> top() const {
    pair<size_t, double> result = {0, numeric_limits<double>::infinity()}; 
    for (auto it = m_storage.begin(); it != m_storage.end(); ++it) {
      if (it->second < result.second) {
        result = *it;
      }
    }
    return result;
  }
  // returns an element with the least distance and removes it from the queue
  pair<size_t, double> pop() {
    const auto t = top(); 
    m_storage.erase(t.first);
    return t;
  }

private:
  map<size_t, double> m_storage;
};

// Implements Dijkstra's shortest path algorithm
// public methods: 
//    ShortestPath(const Graph& graph, size_t from, size_t to)
//    double find()
class ShortestPath {
public:
  // Constructor: expects a graph to find the shortest path for, start (from) and destination (to)
  // vertices
  ShortestPath(const Graph& graph, size_t from, size_t to)
  : m_graph(graph) 
  , m_to_vertex(to) 
  , m_from_vertex(from) 
  , m_unexplored(graph.get_n_vertices())
  , m_trace(graph.get_n_vertices(), to) { 
    // the base case: distance from the start vertex to its own is 0
    m_unexplored.update(from, 0);
  }  
  // Implements Dijkstra's shortest path algorithm
  double find();

private:
  // implements the 1st stage of the algorithm: updates estimated distances for each neighbor 
  // of the given vertex; beside of it updates m_trace array for future route tracing
  void update_estimates(size_t current_vertex, double current_distance) {
    for (const auto& neighbor: m_graph.get_neighbors(current_vertex)) {
      if (m_unexplored.update(neighbor, current_distance + m_graph.get_edge(current_vertex, neighbor))) {
        m_trace[neighbor] = current_vertex;
      }
    }
  }
  // implements the 2nd stage of the algorithm: takes a vertex with the least distance and removes it 
  // from the m_unexplored priority queue
  pair<size_t, double> choose_next_vertex() {
    return m_unexplored.pop();
  }
  // output to a stream operator overload 
  friend ostream& operator<<(ostream& os, ShortestPath& p) {
    // if either m_to_vertex or m_from_vertex lays outside of the graph the route is none
    if (p.m_trace.size() <= p.m_to_vertex || 
        p.m_trace.size() <= p.m_from_vertex) {
      os << "none";
      return os;
    }
    // if start and destination are the same vertex the rout is just a vertex name
    if (p.m_from_vertex == p.m_to_vertex) {
      os << p.m_from_vertex;
      return os;
    }
    // if destination was never visited the route is none
    if (p.m_trace[p.m_to_vertex] == p.m_to_vertex) {
      os << "none";
      return os;
    }
    // traverse the trace array back from the end to the beginning, recovering the shortest route
    vector<size_t> route;
    for (size_t  v = p.m_to_vertex; v != p.m_from_vertex; v = p.m_trace[v]) {
      route.push_back(v);
    }
    route.push_back(p.m_from_vertex);
    // output the recovered route
    auto v = route.rbegin();
    os << *v++;
    for (; v != route.rend(); ++v) {
      os << "->" << *v; 
    }
    return os;
  }

  Graph m_graph;
  size_t m_to_vertex;
  size_t m_from_vertex;
  PriorityQueue m_unexplored;
  vector<size_t> m_trace;
};

double ShortestPath::find() {
  // if either m_to_vertex or m_from_vertex lays outside of the graph the distance is infinity
  if(m_to_vertex >= m_graph.get_n_vertices() ||
     m_from_vertex >= m_graph.get_n_vertices()) {
    return numeric_limits<double>::infinity();
  }
  // begin with the start vertex and zero distance
  auto current_vertex = m_from_vertex;
  auto current_distance = double(0);
  while(true) {
    // 1st step: update distances for each neighbor
    update_estimates(current_vertex, current_distance);
    // 2nd step: choose the next vertex - one with the least distance
    const auto next = choose_next_vertex();
    // reset current vertex
    current_vertex = next.first;
    // reset current distance
    current_distance = next.second;
    // if we just removed the destination vertex from the unexplored queue then we're done
    if (current_vertex == m_to_vertex) {
      return current_distance;
      // if there're no unexplored vertices left then we're in trouble
    } else if (current_distance == numeric_limits<double>::infinity()) {
      break;
    }
  }
  return numeric_limits<double>::infinity();
}

int main() {
  Graph g1(20, 0.1);
  cout << g1 << endl;
  ShortestPath path(g1, 0, 19);
  cout << path.find() << ": " << path << endl;
  return 0;
}