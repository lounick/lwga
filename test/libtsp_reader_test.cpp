#include <iostream>
#include <fstream>

#include "gtest/gtest.h"
#include "libtsp_reader.h"

bool WriteToFile(const std::string &filename, const std::string &content) {
  std::ofstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Could not open output file" << std::endl;
    return false;
  }

  file << content;
  file.close();
  return true;
}

TEST(LibTSPTest, nonExistingFileTest) {
  std::string filename = "/tmp/bad_input.txt";
  libtsp_reader::LIBTSPReader reader(filename);
  EXPECT_EQ(reader.Initialise(), false);
}

TEST(LibTSPTest, existingFileTest) {
  std::string filename = "/tmp/input.txt";

  libtsp_reader::LIBTSPReader reader(filename);

  EXPECT_EQ(WriteToFile(filename, ""), true);
  EXPECT_EQ(reader.Initialise(), true);
}

TEST(LibTSPTest, parseEmptyFile) {
  std::string filename = "/tmp/input.txt";

  libtsp_reader::LIBTSPReader reader(filename);

  EXPECT_EQ(WriteToFile(filename, ""), true);
  EXPECT_EQ(reader.Initialise(), true);
  EXPECT_EQ(reader.ParseFile(), false);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
