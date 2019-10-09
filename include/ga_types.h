//
// Created by nick on 31/03/16.
//

#ifndef LWGA_GA_TYPES_H
#define LWGA_GA_TYPES_H

#include <vector>
#include <limits>
#include <cmath>
#include <cstdint>
#include <memory>
#include <iostream>

template<typename T> using Matrix = std::vector<std::vector<T> >;
template<typename T> using Vector = std::vector<T>;
using VertexId = uint_fast32_t;
using Path = Vector<VertexId>;
using Point2D = std::pair<double, double>;

class SamplingArea {
 public:
  SamplingArea(Point2D origin, uint_fast32_t length, uint_fast32_t width)
      : origin_(origin), length_(length), width_(width) {
    min_x_ = origin_.first - (length_ - 1) / 2.0;
    max_x_ = origin_.first + (length_ - 1) / 2.0;
    min_y_ = origin_.second - (width_ - 1) / 2.0;
    max_y_ = origin_.second + (width_ - 1) / 2.0;
//    start_ = {min_x_ - 1, origin_.second};
//    finish_ = {max_x_ + 1, origin_.second};
  }

//  SamplingArea(Point2D origin, Vector<Point2D> &nodes) :
//      SamplingArea(origin, nodes.front(), nodes.back(), nodes) {}
//
//  SamplingArea(Point2D origin, Point2D start, Point2D finish,
//               const Vector<Point2D> &nodes) :
//      SamplingArea(origin, start, finish, uint_fast32_t(sqrt(nodes.size())),
//                   uint_fast32_t(sqrt(nodes.size())), nodes) {}
//
//  SamplingArea(Point2D origin, Point2D start, Point2D finish,
//               uint_fast32_t length, uint_fast32_t width,
//               const Vector<Point2D> &nodes) : origin_(origin), start_(start),
//                                               finish_(finish), length_(length),
//                                               width_(width) {
//    nodes_ = std::make_shared<Vector<Point2D>>(Vector<Point2D>(nodes));
//  }

  std::shared_ptr<Vector<Point2D>> nodes() const { return nodes_; }
  void nodes(const Vector<Point2D> &nodes) {
    nodes_ = std::make_shared<Vector<Point2D>>(Vector<Point2D>(nodes));
  }
  Point2D origin() const { return origin_; }
  void origin(Point2D origin) { origin_ = origin; }

  Point2D start() const { return start_; }
  void start(Point2D start) { start_ = start; }

  Point2D finish() const { return finish_; }
  void finish(Point2D finish) { finish_ = finish; }

  uint_fast32_t length() const { return length_; }
  void length(uint_fast32_t length) { length_ = length; }

  uint_fast32_t width() const { return width_; }
  void width(uint_fast32_t width) { width_ = width; }

  double_t minX() const { return min_x_; }
  double_t minY() const { return min_y_; }
  double_t maxX() const { return max_x_; }
  double_t maxY() const { return max_y_; }

  bool overlapping(SamplingArea a) const {
    bool ret = false;
    if (max_x_ > a.minX() && a.maxX() > min_x_) {
      if (max_y_ > a.minY() && a.maxY() > min_y_) {
        ret = true;
      }
    }
    return ret;
  }

  Point2D pointFromIdx(uint_fast32_t idx) const {
    return {min_x_ + (idx / width_), min_y_ + (idx % width_)};
  }

  uint_fast32_t idxFromPoint(Point2D p) const {
    return (p.first - min_x_) * width_ + (p.second - min_y_);
  }

  void print(bool nodes = false, std::ostream &out = std::cout) {
    out << "Origin: {" << origin_.first << "," << origin_.second << "}\n";
    out << "Length: " << length_ << "\n";
    out << "Width: " << width_ << "\n";
    out << "Min: {" << min_x_ << "," << min_y_ << "}\n";
    out << "Max: {" << max_x_ << "," << max_y_ << "}\n";
    out << "Start: {" << start_.first << "," << start_.second << "}\n";
    out << "Finish: {" << finish_.first << "," << finish_.second << "}\n";
    out << "Nodes: ";
    if (nodes) {
      out << "[";
      std::string comma = "";
      for (Vector<Point2D>::iterator it = nodes_->begin(); it != nodes_->end();
           ++it) {
        out << comma << "[" << it->first << "," << it->second << "]";
        comma = ", ";
      }
      out << "]" << std::endl;
    }
  }
 private:
  Point2D origin_;
  uint_fast32_t length_;
  uint_fast32_t width_;
  double_t min_x_, max_x_, min_y_, max_y_;
  Point2D start_;
  Point2D finish_;
  std::shared_ptr<Vector<Point2D>> nodes_;

};

#endif //LWGA_GA_TYPES_H
