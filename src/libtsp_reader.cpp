#include "libtsp_reader.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

namespace libtsp_reader {

void RemoveSpaces(std::string &str) {
  str.erase(std::remove_if(str.begin(), str.end(),
                           [](unsigned char x) { return std::isspace(x); }),
            str.end());
}

void ResizeIntMat(Matrix<size_t> &mat, size_t size) {
  mat.resize(size);
  for (Vector<size_t> &v : mat) {
    v.resize(size);
  }
}

size_t CalculateEuclideanDistance(const Node2D &a, const Node2D &b) {
  double_t xd = a.x - b.x;
  double_t yd = a.y - b.y;
  double_t dist = std::sqrt(std::pow(xd, 2) + std::pow(yd, 2)) + 0.5;
  return std::round(dist);
}

size_t CalculateEuclideanDistance(const Node3D &a, const Node3D &b) {
  double_t xd = a.x - b.x;
  double_t yd = a.y - b.y;
  double_t zd = a.z - b.z;
  double_t dist =
      std::sqrt(std::pow(xd, 2) + std::pow(yd, 2) + std::pow(zd, 2)) + 0.5;
  return std::round(dist);
}

size_t CalculateManhattanDistance(const Node2D &a, const Node2D &b) {
  double_t xd = std::abs(a.x - b.x);
  double_t yd = std::abs(a.y - b.y);
  double_t dist = xd + yd + 0.5;
  return std::round(dist);
}

size_t CalculateManhattanDistance(const Node3D &a, const Node3D &b) {
  double_t xd = std::abs(a.x - b.x);
  double_t yd = std::abs(a.y - b.y);
  double_t zd = std::abs(a.z - b.z);
  double_t dist = xd + yd + zd + 0.5;
  return std::round(dist);
}

size_t CalculateMaximumDistance(const Node2D &a, const Node2D &b) {
  double_t xd = std::abs(a.x - b.x);
  double_t yd = std::abs(a.y - b.y);
  return std::max(std::round(xd + 0.5), std::round(yd + 0.5));
}

size_t CalculateMaximumDistance(const Node3D &a, const Node3D &b) {
  double_t xd = std::abs(a.x - b.x);
  double_t yd = std::abs(a.y - b.y);
  double_t zd = std::abs(a.z - b.z);
  return std::max(std::max(std::round(xd + 0.5), std::round(yd + 0.5)),
                  std::round(zd + 0.5));
}

size_t CalculateGeographicalDistance(const Node2D &a, const Node2D &b) {
  double_t PI = 3.141592;

  size_t deg = std::round(a.x);
  double_t min = a.x - deg;
  double_t lat_a = PI * (deg + ((5.0 * min) / 3.0)) / 180.0;

  deg = std::round(a.y);
  min = a.y - deg;
  double_t lon_a = PI * (deg + ((5.0 * min) / 3.0)) / 180.0;

  deg = std::round(b.x);
  min = b.x - deg;
  double_t lat_b = PI * (deg + ((5.0 * min) / 3.0)) / 180.0;

  deg = std::round(b.y);
  min = b.y - deg;
  double_t lon_b = PI * (deg + ((5.0 * min) / 3.0)) / 180.0;

  double_t RRR = 6378.388;

  double_t q1 = std::cos(lon_a - lon_b);
  double_t q2 = std::cos(lat_a - lat_b);
  double_t q3 = std::cos(lat_a + lat_b);

  return std::round(RRR * std::acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) +
                    1.0);
}

size_t CalculatePseudoEuclideanDistance(const Node2D &a, const Node2D &b) {
  double_t xd = a.x - b.x;
  double_t yd = a.y - b.y;
  double_t r = std::sqrt((std::pow(xd, 2) + std::pow(yd, 2)) / 10.0);
  size_t t = std::round(r);
  if (t < r) {
    t += 1;
  }
  return t;
}

size_t CalculateiCeil2DDistance(const Node2D &a, const Node2D &b) {
  double_t xd = a.x - b.x;
  double_t yd = a.y - b.y;
  double_t dist = std::sqrt(std::pow(xd, 2) + std::pow(yd, 2)) + 0.5;
  return std::round(std::ceil(dist));
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
  edge_weight_type_ = EdgeWeightType::UNDEFINED;
  return true;
}
Matrix<size_t> LIBTSPReader::GetCostMatrix() { return cost_mat_; }

Vector<double_t> LIBTSPReader::GetProbabilityVector() { return probabilities_; }

Vector<double_t> LIBTSPReader::GetRewardsVector() { return rewards_; }

double_t LIBTSPReader::GetCapacity() { return capacity_; }

double_t LIBTSPReader::GetMaxCost() { return max_cost_; }
double_t LIBTSPReader::GetMaxReward() { return max_reward_; }
size_t LIBTSPReader::GetStartVertexID() { return start_vertex_id_; }
size_t LIBTSPReader::GetEndVertexID() { return end_vertex_id_; }

bool LIBTSPReader::ParseFile() {
  std::string line;
  while (std::getline(file_stream_, line) && !eof_) {
    if (file_section_ == FileSection::HEADER) {
      if (!HandleHeaderEntry(line)) {
        std::cout << "An error occured. Parsing stopping." << std::endl;
        return false;
      }
    } else {
      if (!HandleDataEntry(line)) {
        std::cout << "An error occured. Parsing stopping." << std::endl;
        return false;
      }
    }
  }
  return GenerateCostMatrix();
}

bool LIBTSPReader::HandleHeaderEntry(const std::string &line) {
  std::string::size_type pos;
  pos = line.find(kDelimeter);

  // If you don't find the delimeter try to handle it as data
  if (pos == std::string::npos) {
    file_section_ = FileSection::DATA;
    return HandleDataEntry(line);
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
    } else if (key.compare(kTMaxIdendtifier) == 0) {
      max_cost_ = std::stod(value);
    } else if (key.compare(kTPrizeIdendtifier) == 0) {
      max_reward_ = std::stod(value);
    } else if (key.compare(kOriginIdendtifier) == 0) {
      start_vertex_id_ = std::stoi(value);
    } else if (key.compare(kDestinationIdendtifier) == 0) {
      end_vertex_id_ = std::stoi(value);
    } else {
      std::cerr << "Processing header entry. " << key << " is not a valid key!"
                << std::endl;
      return false;
    }
  }
  return true;
}

bool LIBTSPReader::HandleDataEntry(const std::string &line) {
  if (line.compare(kNodeCoordSectionIdentifier) == 0) {
    data_section_ = DataSection::NODE_COORD_SECTION;
  } else if (line.compare(kDepotSectionIdentifier) == 0) {
    data_section_ = DataSection::DEPOT_SECTION;
  } else if (line.compare(kDemandSectionIdentifier) == 0) {
    data_section_ = DataSection::DEMAND_SECTION;
  } else if (line.compare(kEdgeDataSectionIdentifier) == 0) {
    data_section_ = DataSection::EDGE_DATA_SECTION;
  } else if (line.compare(kFixedEdgesSectionIdentifier) == 0) {
    data_section_ = DataSection::FIXED_EDGES_SECTION;
  } else if (line.compare(kDisplayDataSectionIdentifier) == 0) {
    data_section_ = DataSection::DISPLAY_DATA_SECTION;
  } else if (line.compare(kTourSectionIdentifier) == 0) {
    data_section_ = DataSection::TOUR_SECTION;
  } else if (line.compare(kEdgeWeightSectionIdentifier) == 0) {
    data_section_ = DataSection::EDGE_WEIGHT_SECTION;
  } else if (line.compare(kNodePrizeProbabilityIdentifier) == 0) {
    data_section_ = DataSection::NODE_PRIZE_PROBABILITY;
  } else if (line.compare(kEOFIdentifier) == 0) {
    eof_ = true;
  } else {
    return HandleSectionEntry(line);
  }
  return true;
}

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
    edge_weight_format_ = EdgeWeightFormat::FUNCTION;
  } else if (value.compare(kFullMatrixIdentifier) == 0) {
    edge_weight_format_ = EdgeWeightFormat::FULL_MATRIX;
  } else if (value.compare(kUpperRowIdentifier) == 0) {
    edge_weight_format_ = EdgeWeightFormat::UPPER_ROW;
  } else if (value.compare(kLowerRowIdentifier) == 0) {
    edge_weight_format_ = EdgeWeightFormat::LOWER_ROW;
  } else if (value.compare(kUpperDiagRowIdentifier) == 0) {
    edge_weight_format_ = EdgeWeightFormat::UPPER_DIAG_ROW;
  } else if (value.compare(kLowerDiagRowIdentifier) == 0) {
    edge_weight_format_ = EdgeWeightFormat::LOWER_DIAG_ROW;
  } else if (value.compare(kUpperColIdentifier) == 0) {
    edge_weight_format_ = EdgeWeightFormat::UPPER_COL;
  } else if (value.compare(kLowerColIdentifier) == 0) {
    edge_weight_format_ = EdgeWeightFormat::LOWER_COL;
  } else if (value.compare(kUpperDiagColIdentifier) == 0) {
    edge_weight_format_ = EdgeWeightFormat::UPPER_DIAG_COL;
  } else if (value.compare(kLowerDiagColIdentifier) == 0) {
    edge_weight_format_ = EdgeWeightFormat::LOWER_DIAG_COL;
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

bool LIBTSPReader::HandleSectionEntry(const std::string &line) {
  if (data_section_ == DataSection::NODE_COORD_SECTION) {
    return HandleNodeCoord(line);
  } else if (data_section_ == DataSection::DEPOT_SECTION) {
    return HandleDepot(line);
  } else if (data_section_ == DataSection::DEMAND_SECTION) {
    return HandleDemand(line);
  } else if (data_section_ == DataSection::EDGE_DATA_SECTION) {
    return HandleEdgeData(line);
  } else if (data_section_ == DataSection::FIXED_EDGES_SECTION) {
    return HandleFixedEdges(line);
  } else if (data_section_ == DataSection::DISPLAY_DATA_SECTION) {
    return HandleDisplayData(line);
  } else if (data_section_ == DataSection::TOUR_SECTION) {
    return HandleTour(line);
  } else if (data_section_ == DataSection::EDGE_WEIGHT_SECTION) {
    return HandleEdgeWeight(line);
  } else if (data_section_ == DataSection::NODE_PRIZE_PROBABILITY) {
    return HandleNodePrizeProbability(line);
  } else {
    std::cout << "Can't process unknown data section entry" << std::endl;
    return false;
  }
  return true;
}

bool LIBTSPReader::HandleNodeCoord(const std::string &line) {
  std::stringstream iss(line);
  if (node_coord_type_ == NodeCoordType::TWOD_COORDS) {
    nodes_2d_.resize(dimension_);
    size_t id;
    Node2D node;
    iss >> id >> node.x >> node.y;
    nodes_2d_[id - 1] = node;
  } else if (node_coord_type_ == NodeCoordType::THREED_COORDS) {
    nodes_3d_.resize(dimension_);
    size_t id;
    Node3D node;
    iss >> id >> node.x >> node.y >> node.z;
    nodes_3d_[id - 1] = node;
  } else if (node_coord_type_ == NodeCoordType::NO_COORDS) {
    std::cout << "Node coordinate type set to NO_COORDS" << std::endl;
    return true;
  } else {
    std::cout << "No coordinate types set. Don't know what to do." << std::endl;
    return false;
  }
  return true;
}

bool LIBTSPReader::HandleDepot(const std::string &line) {
  std::stringstream iss(line);
  size_t depot;
  iss >> depot;
  if (depot > 0) {
    depots_.push_back(depot - 1);
  }
  return true;
}

bool LIBTSPReader::HandleDemand(const std::string &line) {
  std::stringstream iss(line);
  demands_.resize(dimension_);
  size_t id, demand;
  iss >> id >> demand;
  demands_[id - 1] = demand;
  return true;
}

bool LIBTSPReader::HandleEdgeData(const std::string &line) {
  std::stringstream iss(line);
  if (edge_data_format_ == EdgeDataFormat::EDGE_LIST) {
    Edge e;
    iss >> e.start >> e.end;
    if (e.start > 0) {
      --e.start;
      --e.end;
      edge_list_.push_back(e);
    }
  } else if (edge_data_format_ == EdgeDataFormat::ADJ_LIST) {
    adj_list_.resize(dimension_);
    int id;
    iss >> id;
    --id;
    int adj_id;
    iss >> adj_id;
    while (adj_id > 0) {
      --adj_id;
      adj_list_[id].push_back(adj_id);
      iss >> adj_id;
    }
  } else {
    std::cout << "Edge data format not defined." << std::endl;
    return false;
  }
  return true;
}

bool LIBTSPReader::HandleFixedEdges(const std::string &line) {
  std::stringstream iss(line);
  Edge e;
  iss >> e.start >> e.end;
  if (e.start > 0) {
    --e.start;
    --e.end;
    fixed_edges_.push_back(e);
  }
  return true;
}

bool LIBTSPReader::HandleDisplayData(const std::string &line) {
  std::stringstream iss(line);
  if (display_data_type_ == DisplayDataType::TWOD_DISPLAY) {
    display_nodes_.resize(dimension_);
    size_t id;
    Node2D node;
    iss >> id >> node.x >> node.y;
    display_nodes_[id - 1] = node;
  } else {
    std::cout << "Only TWOD_DISPLAY is supported" << std::endl;
  }
  return true;
}

bool LIBTSPReader::HandleTour(const std::string &line) {
  std::stringstream iss(line);
  Vector<int> tour;
  int node;
  iss >> node;
  if (node > 0) {
    while (node > 0) {
      tour.push_back(--node);
      iss >> node;
    }
    tours_.push_back(tour);
  }
  return true;
}

bool LIBTSPReader::HandleEdgeWeight(const std::string &line) {
  std::stringstream iss(line);
  size_t cost;
  while (iss >> cost) {
    edge_weights_.push_back(cost);
  }
}

bool LIBTSPReader::HandleNodePrizeProbability(const std::string &line) {
  std::stringstream iss(line);
  rewards_.resize(dimension_);
  probabilities_.resize(dimension_);
  size_t id;
  double_t reward;
  double_t probability;
  iss >> id >> reward >> probability;
  rewards_[id] = reward;
  probabilities_[id] = probability;
}

bool LIBTSPReader::GenerateCostMatrix() {
  // After reading all the file we need to generate the cost matrix.
  // See if you have nodes or edge weights
  // If Edge weights fill the matrix based on the format
  // if nodes use the correct algorithm to fill the matrix
  if (edge_weight_type_ == EdgeWeightType::EXPLICIT) {
    return GenerateCostMatrixFromEdges();
  } else if (edge_weight_type_ == EdgeWeightType::UNDEFINED) {
    std::cerr << "Undefined edge weight type." << std::endl;
    return false;
  }
  return GenerateCostMatrixFromNodes();
}

// TODO: Break to smaller functions
bool LIBTSPReader::GenerateCostMatrixFromEdges() {
  if (edge_weight_format_ == EdgeWeightFormat::FULL_MATRIX) {
    if (edge_weights_.size() < dimension_ * dimension_) {
      std::cout << "Edge weights are not enough for a full matrix."
                << std::endl;
      return false;
    }
    ResizeIntMat(cost_mat_, dimension_);
    size_t idx = 0;
    for (size_t row = 0; row < dimension_; ++row) {
      for (size_t col = 0; col < dimension_; ++col) {
        cost_mat_[row][col] = edge_weights_[idx++];
      }
    }
  } else if (edge_weight_format_ == EdgeWeightFormat::UPPER_ROW) {
    size_t num_elements = dimension_ * (dimension_ - 1) / 2;
    if (edge_weights_.size() < num_elements) {
      std::cout << "Edge weights are not enough for upper row matrix ."
                << std::endl;
      return false;
    }
    ResizeIntMat(cost_mat_, dimension_);
    size_t idx = 0;
    for (size_t row = 0; row < dimension_; ++row) {
      cost_mat_[row][row] = 0.0;
      for (size_t col = row + 1; col < dimension_; ++col) {
        cost_mat_[row][col] = edge_weights_[idx++];
        cost_mat_[col][row] = cost_mat_[row][col];
      }
    }
  } else if (edge_weight_format_ == EdgeWeightFormat::LOWER_ROW) {
    size_t num_elements = dimension_ * (dimension_ - 1) / 2;
    if (edge_weights_.size() < num_elements) {
      std::cout << "Edge weights are not enough for lower row matrix ."
                << std::endl;
      return false;
    }
    ResizeIntMat(cost_mat_, dimension_);
    size_t idx = 0;
    for (size_t row = 0; row < dimension_; ++row) {
      cost_mat_[row][row] = 0.0;
      for (int col = 0; col < row; ++col) {
        cost_mat_[row][col] = edge_weights_[idx++];
        cost_mat_[col][row] = cost_mat_[row][col];
      }
    }
  } else if (edge_weight_format_ == EdgeWeightFormat::UPPER_DIAG_ROW) {
    size_t num_elements = (dimension_ * (dimension_ - 1) / 2) + dimension_;
    if (edge_weights_.size() < num_elements) {
      std::cout << "Edge weights are not enough for upper diag row matrix ."
                << std::endl;
      return false;
    }
    ResizeIntMat(cost_mat_, dimension_);
    size_t idx = 0;
    for (size_t row = 0; row < dimension_; ++row) {
      for (size_t col = row; col < dimension_; ++col) {
        cost_mat_[row][col] = edge_weights_[idx++];
        cost_mat_[col][row] = cost_mat_[row][col];
      }
    }
  } else if (edge_weight_format_ == EdgeWeightFormat::LOWER_DIAG_ROW) {
    size_t num_elements = (dimension_ * (dimension_ - 1) / 2) + dimension_;
    if (edge_weights_.size() < num_elements) {
      std::cout << "Edge weights are not enough for lower diag row matrix ."
                << std::endl;
      return false;
    }
    ResizeIntMat(cost_mat_, dimension_);
    size_t idx = 0;
    for (size_t row = 0; row < dimension_; ++row) {
      for (size_t col = 0; col < row + 1; ++col) {
        cost_mat_[row][col] = edge_weights_[idx++];
        cost_mat_[col][row] = cost_mat_[row][col];
      }
    }
  } else if (edge_weight_format_ == EdgeWeightFormat::UPPER_COL) {
    std::cout << "UPPER_COL weight format not yet supported" << std::endl;
  } else if (edge_weight_format_ == EdgeWeightFormat::LOWER_COL) {
    std::cout << "LOWER_COL weight format not yet supported" << std::endl;
  } else if (edge_weight_format_ == EdgeWeightFormat::UPPER_DIAG_COL) {
    std::cout << "UPPER_DIAG_COL weight format not yet supported" << std::endl;
  } else if (edge_weight_format_ == EdgeWeightFormat::LOWER_DIAG_COL) {
    std::cout << "LOWER_DIAG_COL weight format not yet supported" << std::endl;
  } else {
    std::cout << "This edge weight format is not supported." << std::endl;
  }
}

bool LIBTSPReader::GenerateCostMatrixFromNodes() {
  if (edge_weight_type_ == EdgeWeightType::EUC_2D) {
    return HandleEuc2DCost();
  } else if (edge_weight_type_ == EdgeWeightType::EUC_3D) {
    return HandleEuc3DCost();
  } else if (edge_weight_type_ == EdgeWeightType::MAX_2D) {
    return HandleMax2DCost();
  } else if (edge_weight_type_ == EdgeWeightType::MAX_3D) {
    return HandleMax3DCost();
  } else if (edge_weight_type_ == EdgeWeightType::MAN_2D) {
    return HandleMan2DCost();
  } else if (edge_weight_type_ == EdgeWeightType::MAN_3D) {
    return HandleMan3DCost();
  } else if (edge_weight_type_ == EdgeWeightType::CEIL_2D) {
    return HandleCeil2DCost();
  } else if (edge_weight_type_ == EdgeWeightType::GEO) {
    return HandleGeoCost();
  } else if (edge_weight_type_ == EdgeWeightType::ATT) {
    return HandleAttCost();
  } else if (edge_weight_type_ == EdgeWeightType::XRAY1) {
    std::cerr << "XRAY1 weight type is not yet supported." << std::endl;
    return false;
  } else if (edge_weight_type_ == EdgeWeightType::XRAY2) {
    std::cerr << "XRAY2 weight type is not yet supported." << std::endl;
    return false;
  } else {
    std::cerr << "This edge weight type is not supported for node generation."
              << std::endl;
    return false;
  }
  return true;
}

bool LIBTSPReader::HandleEuc2DCost() {
  if (nodes_2d_.empty()) {
    std::cout << "EUC_2D edge weight requires 2D nodes." << std::endl;
    return false;
  }

  if (nodes_2d_.size() != dimension_) {
    std::cout << "Nodes are not equal to the defined dimension." << std::endl;
    return false;
  }
  ResizeIntMat(cost_mat_, dimension_);
  for (size_t row = 0; row < dimension_; ++row) {
    for (size_t col = row; col < dimension_; ++col) {
      cost_mat_[row][col] =
          CalculateEuclideanDistance(nodes_2d_[row], nodes_2d_[col]);
      cost_mat_[col][row] = cost_mat_[row][col];
    }
  }
  return true;
}

bool LIBTSPReader::HandleEuc3DCost() {
  if (nodes_3d_.empty()) {
    std::cout << "EUC_3D edge weight requires 3D nodes." << std::endl;
    return false;
  }

  if (nodes_3d_.size() != dimension_) {
    std::cout << "Nodes are not equal to the defined dimension." << std::endl;
    return false;
  }
  ResizeIntMat(cost_mat_, dimension_);
  for (size_t row = 0; row < dimension_; ++row) {
    for (size_t col = row; col < dimension_; ++col) {
      cost_mat_[row][col] =
          CalculateEuclideanDistance(nodes_3d_[row], nodes_3d_[col]);
      cost_mat_[col][row] = cost_mat_[row][col];
    }
  }
  return true;
}

bool LIBTSPReader::HandleMax2DCost() {
  if (nodes_2d_.empty()) {
    std::cout << "MAX_2D edge weight requires 2D nodes." << std::endl;
    return false;
  }

  if (nodes_2d_.size() != dimension_) {
    std::cout << "Nodes are not equal to the defined dimension." << std::endl;
    return false;
  }
  ResizeIntMat(cost_mat_, dimension_);
  for (size_t row = 0; row < dimension_; ++row) {
    for (size_t col = row; col < dimension_; ++col) {
      cost_mat_[row][col] =
          CalculateMaximumDistance(nodes_2d_[row], nodes_2d_[col]);
      cost_mat_[col][row] = cost_mat_[row][col];
    }
  }
  return true;
}

bool LIBTSPReader::HandleMax3DCost() {
  if (nodes_3d_.empty()) {
    std::cout << "MAX_3D edge weight requires 3D nodes." << std::endl;
    return false;
  }

  if (nodes_3d_.size() != dimension_) {
    std::cout << "Nodes are not equal to the defined dimension." << std::endl;
    return false;
  }
  ResizeIntMat(cost_mat_, dimension_);
  for (size_t row = 0; row < dimension_; ++row) {
    for (size_t col = row; col < dimension_; ++col) {
      cost_mat_[row][col] =
          CalculateEuclideanDistance(nodes_3d_[row], nodes_3d_[col]);
      cost_mat_[col][row] = cost_mat_[row][col];
    }
  }
  return true;
}

bool LIBTSPReader::HandleMan2DCost() {
  if (nodes_2d_.empty()) {
    std::cout << "MAN_2D edge weight requires 2D nodes." << std::endl;
    return false;
  }

  if (nodes_2d_.size() != dimension_) {
    std::cout << "Nodes are not equal to the defined dimension." << std::endl;
    return false;
  }
  ResizeIntMat(cost_mat_, dimension_);
  for (size_t row = 0; row < dimension_; ++row) {
    for (size_t col = row; col < dimension_; ++col) {
      cost_mat_[row][col] =
          CalculateManhattanDistance(nodes_2d_[row], nodes_2d_[col]);
      cost_mat_[col][row] = cost_mat_[row][col];
    }
  }
  return true;
}

bool LIBTSPReader::HandleMan3DCost() {
  if (nodes_3d_.empty()) {
    std::cout << "MAN_3D edge weight requires 3D nodes." << std::endl;
    return false;
  }

  if (nodes_3d_.size() != dimension_) {
    std::cout << "Nodes are not equal to the defined dimension." << std::endl;
    return false;
  }
  ResizeIntMat(cost_mat_, dimension_);
  for (size_t row = 0; row < dimension_; ++row) {
    for (size_t col = row; col < dimension_; ++col) {
      cost_mat_[row][col] =
          CalculateEuclideanDistance(nodes_3d_[row], nodes_3d_[col]);
      cost_mat_[col][row] = cost_mat_[row][col];
    }
  }
  return true;
}

bool LIBTSPReader::HandleCeil2DCost() {
  if (nodes_2d_.empty()) {
    std::cout << "CEIL_2D edge weight requires 2D nodes." << std::endl;
    return false;
  }

  if (nodes_2d_.size() != dimension_) {
    std::cout << "Nodes are not equal to the defined dimension." << std::endl;
    return false;
  }
  ResizeIntMat(cost_mat_, dimension_);
  for (size_t row = 0; row < dimension_; ++row) {
    for (size_t col = row; col < dimension_; ++col) {
      cost_mat_[row][col] =
          CalculateiCeil2DDistance(nodes_2d_[row], nodes_2d_[col]);
      cost_mat_[col][row] = cost_mat_[row][col];
    }
  }
  return true;
}

bool LIBTSPReader::HandleGeoCost() {
  if (nodes_2d_.empty()) {
    std::cout << "GEO edge weight requires 2D nodes." << std::endl;
    return false;
  }

  if (nodes_2d_.size() != dimension_) {
    std::cout << "Nodes are not equal to the defined dimension." << std::endl;
    return false;
  }
  ResizeIntMat(cost_mat_, dimension_);
  for (size_t row = 0; row < dimension_; ++row) {
    for (size_t col = row; col < dimension_; ++col) {
      cost_mat_[row][col] =
          CalculateGeographicalDistance(nodes_2d_[row], nodes_2d_[col]);
      cost_mat_[col][row] = cost_mat_[row][col];
    }
  }
  return true;
}

bool LIBTSPReader::HandleAttCost() {
  if (nodes_2d_.empty()) {
    std::cout << "ATT edge weight requires 2D nodes." << std::endl;
    return false;
  }

  if (nodes_2d_.size() != dimension_) {
    std::cout << "Nodes are not equal to the defined dimension." << std::endl;
    return false;
  }
  ResizeIntMat(cost_mat_, dimension_);
  for (size_t row = 0; row < dimension_; ++row) {
    for (size_t col = row; col < dimension_; ++col) {
      cost_mat_[row][col] =
          CalculatePseudoEuclideanDistance(nodes_2d_[row], nodes_2d_[col]);
      cost_mat_[col][row] = cost_mat_[row][col];
    }
  }
  return true;
}

}  // namespace libtsp_reader
