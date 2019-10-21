#ifndef LWGA_LIBTSP_READER_H_
#define LWGA_LIBTSP_READER_H_

#include <fstream>
#include <string>
#include <vector>

namespace libtsp_reader {

template <typename T>
using Matrix = std::vector<std::vector<T>>;

template <typename T>
using Vector = std::vector<T>;

enum class ProblemType {
  TSP,
  ATSP,
  SOP,
  HCP,
  CVRP,
  TOUR,
  POP,
};

enum class FileSection {
  HEADER,
  DATA,
}

class LIBTSPReader {
 public:
  LIBTSPReader(std::string file_name);
  ~LIBTSPReader();
  bool Initialise();
  Matrix<double_t> GetCostMatrix();
  Vector<double_t> GetProbabilityVector();
  Vector<double_t> GetRewardsVector();
  double_t GetCapacity();

 private:
  void ParseFile();
  std::string file_name_;
  std::string problem_name_;
  ProblemType problem_type_;
  FileSection file_section_;
  std::ifstream file_stream_;
  Matrix<double_t> cost_mat_;
  Vector<double_t> probabilities_;
  Vector<double_t> rewards_;
  double_t capacity_;
};
}  // namespace libtsp_reader
#endif  // LWGA_LIBTSP_READER_H_
