#include "libtsp_reader.h"

#include <sstream>

namespace libtsp_reader {

LIBTSPReader::LIBTSPReader(std::string file_name) { file_name_ = file_name; }

LIBTSPReader::~LIBTSPReader() { file_stream_.close(); }

bool LIBTSPReader::Initialise() {
  file_stream_.open(file_name_);
  if (!file_stream_.is_open()) {
    return false;
  }
  file_section_ = FileSection::HEADER;
  return true;
}
Matrix<double_t> LIBTSPReader::GetCostMatrix() { return cost_mat_; }

Vector<double_t> LIBTSPReader::GetProbabilityVector() { return probabilties_; }

Vector<double_t> LIBTSPReader::GetRewardsVector() { return rewards_; }
double_t LIBTSPReader::GetCapacity() { return capacity_; }
void LIBTSPReader::ParseFile() {
  std::string line;
  while (std::getline(file_stream_, line)) {
    if (file_section == FileSection::HEADER) {
      HandleHeaderEntry(line);
    } else {
      HandleDataEntry(line);
    }
  }
}
void LibTSPReader::HandleHeaderEntry(const std::string &line) {
  std::istringstream iss(line);
  std::string key;
  iss >> key;

  // Remove the trailing colon if in the same word
  if (key.back() == ':') {
    key.pop_back();
  }
}
void LibTSPReader::HandleDataEntry(const std::string &line) {}
}  // namespace libtsp_reader
