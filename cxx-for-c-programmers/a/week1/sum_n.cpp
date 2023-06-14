// Convert this program to C++
// change to C++ io
// change to one line comments
// change defines of constants to const
// change array to vector<>
// inline any short function

#include <iostream>
#include <vector>

using namespace std;

const int N = 40;

template <class summable>
inline void sum(summable& p, const vector<summable> d) {
  // don't try to set p to 'zero' here since it's exact type is unknown
  // caller shall take care of it
  for(summable i = 0; i < d.size(); ++i) {
    p += d[i];
  }
}

int main() {
  vector<int> data;
  for(int i = 0; i < N; ++i) {
    data.push_back(i);
  }
  int accum = 0;

  sum(accum, data);
  cout << "sum is " << accum << endl;
  return 0;
}