#include <cstddef>
#include <ostream>
#include <vector>
#include <iostream>
#include <random>
#include <iomanip>

using namespace std;

class Graph {
public:
  Graph(size_t n_v, double density, double min_distance, double max_distance);

  Graph(size_t n_v=10)
  : m_n_v(n_v)
  , m_n_e(0) 
  , m_matrix(n_v, vector<double>(n_v, -1)) {}

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
  
  size_t n_v() const {
    return m_n_v;
  }

  size_t n_e() const {
    return m_n_e;
  }

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

private:
  vector<vector<double>> m_matrix;  
  size_t m_n_v;
  size_t m_n_e;
};

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
        m_matrix[i][j] = m_matrix[j][i] = dist_distribution(engine);
        m_n_e ++;
      }
    }
  }
}

int main() {
  Graph g1(10, 0.2);
  cout << g1 << endl;
  Graph g2 = g1;
  cout << g2 << endl;
  return 0;
}