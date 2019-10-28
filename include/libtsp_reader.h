#ifndef LWGA_LIBTSP_READER_H_
#define LWGA_LIBTSP_READER_H_

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace libtsp_reader {

template <typename T>
using Matrix = std::vector<std::vector<T>>;

template <typename T>
using Vector = std::vector<T>;

struct Node2D {
  double_t x;
  double_t y;
};

struct Node3D {
  double_t x;
  double_t y;
  double_t z;
};

struct Edge {
  int start;
  int end;
};

enum class ProblemType {
  TSP,
  ATSP,
  SOP,
  HCP,
  CVRP,
  TOUR,
  POP,
};

enum class EdgeWeightType {
  EXPLICIT,
  EUC_2D,
  EUC_3D,
  MAX_2D,
  MAX_3D,
  MAN_2D,
  MAN_3D,
  CEIL_2D,
  GEO,
  ATT,
  XRAY1,
  XRAY2,
  SPECIAL,
};

enum class EdgeWeightFormat {
  FUNCTION,
  FULL_MATRIX,
  UPPER_ROW,
  LOWER_ROW,
  UPPER_DIAG_ROW,
  LOWER_DIAG_ROW,
  UPPER_COL,
  LOWER_COL,
  UPPER_DIAG_COL,
  LOWER_DIAG_COL,
};

enum class EdgeDataFormat {
  EDGE_LIST,
  ADJ_LIST,
};

enum class NodeCoordType {
  TWOD_COORDS,
  THREED_COORDS,
  NO_COORDS,
};

enum class DisplayDataType {
  COORD_DISPLAY,
  TWOD_DISPLAY,
  NO_DISPLAY,
};

enum class FileSection {
  HEADER,
  DATA,
};

enum class DataSection {
  NODE_COORD_SECTION,
  DEPOT_SECTION,
  DEMAND_SECTION,
  EDGE_DATA_SECTION,
  FIXED_EDGES_SECTION,
  DISPLAY_DATA_SECTION,
  TOUR_SECTION,
  EDGE_WEIGHT_SECTION,
  NODE_PRIZE_PROBABILITY,
};

class LIBTSPReader {
 public:
  LIBTSPReader(std::string file_name);
  ~LIBTSPReader();
  bool Initialise();
  Matrix<size_t> GetCostMatrix();
  Vector<double_t> GetProbabilityVector();
  Vector<double_t> GetRewardsVector();
  double_t GetMaxCost();
  double_t GetMaxReward();
  size_t GetStartVertexID();
  size_t GetEndVertexID();
  double_t GetCapacity();
  bool ParseFile();

 private:
  bool HandleHeaderEntry(const std::string &line);
  bool HandleDataEntry(const std::string &line);
  bool HandleProblemType(const std::string &value);
  bool HandleEdgeWeightType(const std::string &value);
  bool HandleEdgeWeightFormat(const std::string &value);
  bool HandleEdgeDataFormat(const std::string &value);
  bool HandleNodeCoordType(const std::string &value);
  bool HandleDispleyDataType(const std::string &value);
  bool HandleSectionEntry(const std::string &line);
  bool HandleNodeCoord(const std::string &line);
  bool HandleDepot(const std::string &line);
  bool HandleDemand(const std::string &line);
  bool HandleEdgeData(const std::string &line);
  bool HandleFixedEdges(const std::string &line);
  bool HandleDisplayData(const std::string &line);
  bool HandleTour(const std::string &line);
  bool HandleEdgeWeight(const std::string &line);
  bool HandleNodePrizeProbability(const std::string &line);
  bool GenerateCostMatrix();
  bool GenerateCostMatrixFromEdges();
  bool GenerateCostMatrixFromNodes();
  bool HandleEuc2DCost();
  bool HandleEuc3DCost();
  bool HandleMax2DCost();
  bool HandleMax3DCost();
  bool HandleMan2DCost();
  bool HandleMan3DCost();
  bool HandleCeil2DCost();
  bool HandleGeoCost();
  bool HandleAttCost();

  std::string file_name_;
  std::string problem_name_;
  ProblemType problem_type_;
  FileSection file_section_;
  DataSection data_section_;
  EdgeWeightType edge_weight_type_;
  EdgeWeightFormat edge_weight_format_;
  EdgeDataFormat edge_data_format_;
  NodeCoordType node_coord_type_;
  DisplayDataType display_data_type_;
  std::ifstream file_stream_;
  Matrix<size_t> cost_mat_;
  Vector<double_t> probabilities_;
  Vector<double_t> rewards_;
  double_t max_cost_;
  double_t max_reward_;
  size_t start_vertex_id_;
  size_t end_vertex_id_;
  int capacity_;
  int dimension_;
  bool eof_;

  Vector<Node2D> nodes_2d_;
  Vector<Node3D> nodes_3d_;
  Vector<int> depots_;
  Vector<size_t> demands_;
  Vector<Edge> edge_list_;
  Matrix<int> adj_list_;
  Vector<Edge> fixed_edges_;
  Vector<Node2D> display_nodes_;
  Matrix<int> tours_;
  Vector<size_t> edge_weights_;

  const std::string kNameIdentifier = "NAME";
  const std::string kTypeIdentifier = "TYPE";
  const std::string kCommentIdentifier = "COMMENT";
  const std::string kDimensionIdentifier = "DIMENSION";
  const std::string kCapacityIdentifier = "CAPACITY";
  const std::string kEdgeWeightTypeIdentifier = "EDGE_WEIGHT_TYPE";
  const std::string kEdgeWeightFormatIdentifier = "EDGE_WEIGHT_FORMAT";
  const std::string kEdgeDataFormatIdentifier = "EDGE_DATA_FORMAT";
  const std::string kNodeCoordTypeIdentifier = "NODE_COORD_TYPE";
  const std::string kDisplayDataTypeIdentifier = "DISPLAY_DATA_TYPE";
  const std::string kEOFIdentifier = "EOF";
  const std::string kTMaxIdendtifier = "TMAX";
  const std::string kTPrizeIdendtifier = "TPRIZE";
  const std::string kOriginIdendtifier = "ORIGIN";
  const std::string kDestinationIdendtifier = "DESTINATION";

  const std::string kTSPIdentifier = "TSP";
  const std::string kATSPIdentifier = "ATSP";
  const std::string kSOPIdentifier = "SOP";
  const std::string kHCPIdentifier = "HCP";
  const std::string kCVRPIdentifier = "CVRP";
  const std::string kTOURIdentifier = "TOUR";
  const std::string kPOPIdentifier = "POP";

  const std::string kExplicitIdentifier = "EXPLICIT";
  const std::string kEuc2DIdentifier = "EUC_2D";
  const std::string kEuc3DIdentifier = "EUC_3D";
  const std::string kMax2DIdentifier = "MAX_2D";
  const std::string kMax3DIdentifier = "MAX_3D";
  const std::string kMan2DIdentifier = "MAN_2D";
  const std::string kMan3DIdentifier = "MAN_3D";
  const std::string kCeil2DIdentifier = "CEIL_2D";
  const std::string kGeopIdentifier = "GEO";
  const std::string kAttIdentifier = "ATT";
  const std::string kXRay1Identifier = "XRAY1";
  const std::string kXRay2Identifier = "XRAY2";
  const std::string kSpecialIdentifier = "SPECIAL";

  const std::string kFunctionIdentifier = "FUNCTION";
  const std::string kFullMatrixIdentifier = "FULL_MATRIX";
  const std::string kUpperRowIdentifier = "UPPER_ROW";
  const std::string kLowerRowIdentifier = "LOWER_ROW";
  const std::string kUpperDiagRowIdentifier = "UPPER_DIAG_ROW";
  const std::string kLowerDiagRowIdentifier = "LOWER_DIAG_ROW";
  const std::string kUpperColIdentifier = "UPPER_COL";
  const std::string kLowerColIdentifier = "LOWER_COL";
  const std::string kUpperDiagColIdentifier = "UPPER_DIAG_COL";
  const std::string kLowerDiagColIdentifier = "LOWER_DIAG_COL";

  const std::string kEdgeListIdentifier = "EDGE_LIST";
  const std::string kAdjListIdentifier = "ADJ_LIST";

  const std::string kTwoDIdentifier = "TWOD_COORDS";
  const std::string kThreeDIdentifier = "THREED_COORDS";
  const std::string kNoCoordsIdentifier = "NO_COORDS";

  const std::string kCoordDisplayIdentifier = "COORD_DISPLAY";
  const std::string kTwoDDisplayIdentifier = "TWOD_DISPLAY";
  const std::string kNoDisplayIdentifier = "NO_DISPLAY";

  const std::string kDelimeter = ":";

  const std::string kNodeCoordSectionIdentifier = "NODE_COORD_SECTION";
  const std::string kDepotSectionIdentifier = "DEPOT_SECTION";
  const std::string kDemandSectionIdentifier = "DEMAND_SECTION";
  const std::string kEdgeDataSectionIdentifier = "EDGE_DATA_SECTION";
  const std::string kFixedEdgesSectionIdentifier = "FIXED_EDGES_SECTION";
  const std::string kDisplayDataSectionIdentifier = "DISPLAY_DATA_SECTION";
  const std::string kTourSectionIdentifier = "TOUR_SECTION";
  const std::string kEdgeWeightSectionIdentifier = "EDGE_WEIGHT_SECTION";
  const std::string kNodePrizeProbabilityIdentifier = "NODE_PRIZE_PROBABILITY";
};
}  // namespace libtsp_reader
#endif  // LWGA_LIBTSP_READER_H_
