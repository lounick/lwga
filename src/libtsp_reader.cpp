#include "libtsp_reader.h"

#include <sstream>

namespace libtsp_reader {

void RemoveSpaces(std::string &str) {
  str.erase(std::remove_if(str.begin(), str.end(),
                           [](unsigned char x) { return std::isspace(x); }),
            str.end());
}

LIBTSPReader::LIBTSPReader(std::string file_name) {
  file_name_ = file_name;
  eof_ = false;
}

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
  while (std::getline(file_stream_, line) && !eof_) {
    if (file_section == FileSection::HEADER) {
      HandleHeaderEntry(line);
    } else {
      HandleDataEntry(line);
    }
  }
}
bool LIBTSPReader::HandleHeaderEntry(const std::string &line) {
  std::string::size_type pos;
  pos = line.find(kDelimeter);

  // If you don't find the delimeter try to handle it as data
  if (pos == std::string::npos) {
    file_section_ = FileSection::DATA;
    HandleDataEntry(line);
  } else {
    std::string key = line.substr(0, pos);
    std::string value = line.substr(pos + 1, std::string::npos);
    RemoveSpaces(key);
    RemoveSpaces(value);
    if (key.compare(kNameIdentifier) == 0) {
      problem_name_ = value;
    } else if (key.compare(kTypeIdentifier) == 0) {
      return HandleProblemType(value);
    } else if (key.compare(kCommentIdentifier) == 0) {
      return true;  // Ignore the comment
    } else if (key.compare(kDimensionIdentifier) == 0) {
      dimension_ = std::stoi(value);
    } else if (key.compare(kCapacityIdentifier) == 0) {
      capacity_ = std::stoi(value);
    } else if (key.compare(kEdgeWeightTypeIdentifier) == 0) {
      return HandleEdgeWeightType(value);
    } else if (key.compare(kEdgeWeightFormatIdentifier) == 0) {
      return HandleEdgeWeightFormat(value);
    } else if (key.compare(kEdgeDataFormatIdentifier) == 0) {
      return HandleEdgeDataFormat(value);
    } else if (key.compare(kNodeCoordTypeIdentifier) == 0) {
      return HandleNodeCoordType(value);
    } else if (key.compare(kDisplayDataTypeIdentifier) == 0) {
      return HandleDispleyDataType(value);
    } else if (key.compare(kEOFIdentifier) == 0) {
      eof_ = true;
    } else {
      std::cout << "Processing header entry. " << key << " is not a valid key!"
                << std::endl;
      return false;
    }
  }
}

bool LIBTSPReader::HandleDataEntry(const std::string &line) {}

bool LIBTSPReader::HandleProblemType(const std::string &value) {
  if (value.compare(kTSPIdentifier) == 0) {
    problem_type_ = ProblemType::TSP;
  } else if (value.compare(kATSPIdentifier) == 0) {
    problem_type_ = ProblemType::ATSP;
  } else if (value.compare(kSOPIdentifier) == 0) {
    problem_type_ = ProblemType::SOP;
  } else if (value.compare(kHCPIdentifier) == 0) {
    problem_type_ = ProblemType::HCP;
  } else if (value.compare(kCVRPIdentifier) == 0) {
    problem_type_ = ProblemType::CVRP;
  } else if (value.compare(kTOURIdentifier) == 0) {
    problem_type_ = ProblemType::TOUR;
  } else if (value.compare(kPOPIdentifier) == 0) {
    problem_type_ = ProblemType::POP;
  } else {
    std::cout << "Processing problem type. " << value
              << " is not a valid value!" << std::endl;
    return false;
  }
  return true;
}

bool LIBTSPReader::HandleEdgeWeightType(const std::string &value) {
  if (value.compare(kExplicitIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::EXPLICIT;
  } else if (value.compare(kEuc2DIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::EUC_2D;
  } else if (value.compare(kEuc3DIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::EUC_3D;
  } else if (value.compare(kMax2DIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::MAX_2D;
  } else if (value.compare(kMax3DIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::MAX_3D;
  } else if (value.compare(kMan2DIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::MAN_2D;
  } else if (value.compare(kMan3DIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::MAN_3D;
  } else if (value.compare(kCeil2DIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::CEIL_2D;
  } else if (value.compare(kGeopIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::GEO;
  } else if (value.compare(kAttIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::ATT;
  } else if (value.compare(kXRay1Identifier) == 0) {
    edge_weight_type_ = EdgeWeightType::XRAY1;
  } else if (value.compare(kXRay2Identifier) == 0) {
    edge_weight_type_ = EdgeWeightType::XRAY2;
  } else if (value.compare(kSpecialIdentifier) == 0) {
    edge_weight_type_ = EdgeWeightType::SPECIAL;
  } else {
    std::cout << "Processing edge weight type. " << value
              << " is not a valid value!" << std::endl;
    return false;
  }
  return true;
}

bool LIBTSPReader::HandleEdgeWeightFormat(const std::string &value) {
  if (value.compare(kFunctionIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::FUNCTION;
  } else if (value.compare(kFullMatrixIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::FULL_MATRIX;
  } else if (value.compare(kUpperRowIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::UPPER_ROW;
  } else if (value.compare(kLowerRowIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::LOWER_ROW;
  } else if (value.compare(kUpperDiagRowIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::UPPER_DIAG_ROW;
  } else if (value.compare(kLowerDiagRowIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::LOWER_DIAG_ROW;
  } else if (value.compare(kUpperColIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::UPPER_COL;
  } else if (value.compare(kLowerColIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::LOWER_COL;
  } else if (value.compare(kUpperDiagColIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::UPPER_DIAG_COL;
  } else if (value.compare(kLowerDiagColIdentifier) == 0) {
    edge_weight_format_ = WeightFormat::LOWER_DIAG_COL;
  } else {
    std::cout << "Processing edge weight format. " << value
              << " is not a valid value!" << std::endl;
    return false;
  }
  return true;
}

bool LIBTSPReader::HandleEdgeDataFormat(const std::string &value) {
  if (value.compare(kEdgeListIdentifier) == 0) {
    edge_data_format_ = EdgeDataFormat::EDGE_LIST;
  } else if (value.compare(kAdjListIdentifier) == 0) {
    edge_data_format_ = EdgeDataFormat::ADJ_LIST;
  } else {
    std::cout << "Processing edge data format. " << value
              << " is not a valid value!" << std::endl;
    return false;
  }
  return true;
}

bool LIBTSPReader::HandleNodeCoordType(const std::string &value) {
  if (value.compare(kTwoDIdentifier) == 0) {
    node_coord_type_ = NodeCoordType::TWOD_COORDS;
  } else if (value.compare(kThreeDIdentifier) == 0) {
    node_coord_type_ = NodeCoordType::THREED_COORDS;
  } else if (value.compare(kNoCoordsIdentifier) == 0) {
    node_coord_type_ = NodeCoordType::NO_COORDS;
  } else {
    std::cout << "Processing node coord type. " << value
              << " is not a valid value!" << std::endl;
    return false;
  }
  return true;
}

bool LIBTSPReader::HandleDispleyDataType(const std::string &value) {
  if (value.compare(kCoordDisplayIdentifier) == 0) {
    display_data_type_ = DisplayDataType::COORD_DISPLAY;
  } else if (value.compare(kTwoDDisplayIdentifier) == 0) {
    display_data_type_ = DisplayDataType::TWOD_DISPLAY;
  } else if (value.compare(kNoDisplayIdentifier) == 0) {
    display_data_type_ = DisplayDataType::NO_DISPLAY;
  } else {
    std::cout << "Processing display data type. " << value
              << " is not a valid value!" << std::endl;
    return false;
  }
  return true;
}
}  // namespace libtsp_reader
