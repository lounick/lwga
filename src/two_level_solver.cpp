//
// Created by nick on 17/03/18.
//

#include <two_level_solver.h>

double_t getAreaStart(Point2D point, SamplingArea &area) {
  Point2D intersection;
  Point2D area_start, area_end;
  Point2D prev[2];
  Point2D next[2];
  Point2D line_start = point;
  Point2D line_end = area.origin();
  double_t dx = line_end.first - line_start.first;
  double_t dy = line_end.second - line_start.second;

  if (dx != 0.0) {
    if (std::abs(dy / dx) < 1.0) {
      // Front or back based on dx
      if (dx > 0) {
        // Back
        Point2D area_segment_start = {area.minX(), area.minY()};
        Point2D area_segment_end = {area.minX(), area.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      } else {
        // Front
        Point2D area_segment_start = {area.maxX(), area.minY()};
        Point2D area_segment_end = {area.maxX(), area.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      }
      double_t diff = intersection.second - area.minY();
      double_t inter_1 = area.minY() + std::floor(diff);
      double_t inter_2 = area.minY() + std::ceil(diff);
      next[0] = {intersection.first, inter_1};
      next[1] = {intersection.first, inter_2};
      if (find_distance(line_start, next[0])
          < find_distance(line_start, next[1])) {
        area_start = next[0];
      } else {
        area_start = next[1];
      }
    } else if (std::abs(dy / dx) > 1.0) {
      // Top or bottom based on dy
      if (dy > 0) {
        // Bottom
        Point2D area_segment_start = {area.minX(), area.minY()};
        Point2D area_segment_end = {area.maxX(), area.minY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      } else {
        // Top
        Point2D area_segment_start = {area.minX(), area.maxY()};
        Point2D area_segment_end = {area.maxX(), area.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      }
      double_t diff = intersection.first - area.minX();
      double_t inter_1 = area.minX() + std::floor(diff);
      double_t inter_2 = area.minX() + std::ceil(diff);
      next[0] = {inter_1, intersection.second};
      next[1] = {inter_2, intersection.second};
      if (find_distance(line_start, next[0])
          < find_distance(line_start, next[1])) {
        area_start = next[0];
      } else {
        area_start = next[1];
      }
    } else {
      // Corner based on dx and dy
      if (dx > 0) {
        if (dy > 0) {
          // Bottom left corner
          area_start = {area.minX(), area.minY()};
        } else {
          // Top left corner
          area_start = {area.minX(), area.maxY()};
        }
      } else {
        if (dy > 0) {
          // Bottom right corner
          area_start = {area.maxX(), area.minY()};
        } else {
          // Top right corner
          area_start = {area.maxX(), area.maxY()};
        }
      }
    }
  } else {
    // Top or bottom based on dy
    if (dy > 0) {
      // Bottom
      Point2D area_segment_start = {area.minX(), area.minY()};
      Point2D area_segment_end = {area.maxX(), area.minY()};
      LineIntersect(line_start.first, line_start.second,
                    line_end.first, line_end.second,
                    area_segment_start.first, area_segment_start.second,
                    area_segment_end.first, area_segment_end.second,
                    intersection.first, intersection.second);
    } else {
      // Top
      Point2D area_segment_start = {area.minX(), area.maxY()};
      Point2D area_segment_end = {area.maxX(), area.maxY()};
      LineIntersect(line_start.first, line_start.second,
                    line_end.first, line_end.second,
                    area_segment_start.first, area_segment_start.second,
                    area_segment_end.first, area_segment_end.second,
                    intersection.first, intersection.second);
    }
    double_t diff = intersection.first - area.minX();
    double_t inter_1 = area.minX() + std::floor(diff);
    double_t inter_2 = area.minX() + std::ceil(diff);
    next[0] = {inter_1, intersection.second};
    next[1] = {inter_2, intersection.second};
    if (find_distance(line_start, next[0])
        < find_distance(line_start, next[1])) {
      area_start = next[0];
    } else {
      area_start = next[1];
    }
  }
  if (area_start.first < area.minX() || area_start.first > area.maxX()) {
    std::cout << "area_start x out of bounds" << std::endl;
    // TODO: This shouldn't happen. Fix it.
    if (area_start.first < area.minX()) {
      area_start.first = area.minX();
    } else {
      area_start.first = area.maxX();
    }
    std::cout << "(" << point.first << ", " << point.second << ")->("
              << area.origin().first << ", " << area.origin().second << ")"
              << std::endl;
  }
  if (area_start.second < area.minY() || area_start.second > area.maxY()) {
    std::cout << "area_start y out of bounds" << std::endl;
    // TODO: This shouldn't happen. Fix it.
    if (area_start.second < area.minY()) {
      area_start.second = area.minY();
    } else {
      area_start.second = area.maxY();
    }
    std::cout << "(" << point.first << ", " << point.second << ")->("
              << area.origin().first << ", " << area.origin().second << ")"
              << std::endl;
  }
  area.start(area_start);
  return find_distance(line_start, area_start);
}

double_t getAreaStartEnd(SamplingArea &from, SamplingArea &to) {
  Point2D intersection;
  Point2D area_start, area_end;
  Point2D prev[2];
  Point2D next[2];
  double_t min_dist;
  Point2D line_start = from.origin();
  Point2D line_end = to.origin();
  double_t dx = line_end.first - line_start.first;
  double_t dy = line_end.second - line_start.second;

  if (dx != 0.0) {
    // For point exiting "from" area
    if (std::abs(dy / dx) < 1.0) {
      // Front or back based on dx
      if (dx > 0) {
        // Front
        Point2D area_segment_start = {from.maxX(), from.minY()};
        Point2D area_segment_end = {from.maxX(), from.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      } else {
        // Back
        Point2D area_segment_start = {from.minX(), from.minY()};
        Point2D area_segment_end = {from.minX(), from.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      }
      double_t diff = intersection.second - from.minY();
      double_t inter_1 = from.minY() + std::floor(diff);
      double_t inter_2 = from.minY() + std::ceil(diff);
      prev[0] = {intersection.first, inter_1};
      prev[1] = {intersection.first, inter_2};
    } else if (std::abs(dy / dx) > 1.0) {
      // Top or bottom based on dy
      if (dy > 0) {
        // Top
        Point2D area_segment_start = {from.minX(), from.maxY()};
        Point2D area_segment_end = {from.maxX(), from.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      } else {
        // Bottom
        Point2D area_segment_start = {from.minX(), from.minY()};
        Point2D area_segment_end = {from.maxX(), from.minY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      }
      double_t diff = intersection.first - from.minX();
      double_t inter_1 = from.minX() + std::floor(diff);
      double_t inter_2 = from.minX() + std::ceil(diff);
      prev[0] = {inter_1, intersection.second};
      prev[1] = {inter_2, intersection.second};
    } else {
      // Corner based on dx and dy
      if (dx > 0) {
        if (dy > 0) {
          // Top right
          area_end = {from.maxX(), from.maxY()};
        } else {
          // Bottom right
          area_end = {from.maxX(), from.minY()};
        }
      } else {
        if (dy > 0) {
          // Top left
          area_end = {from.minX(), from.maxY()};
        } else {
          // Bottom left
          area_end = {from.minX(), from.minY()};
        }
      }
    }
    // For point entering "to" area
    if (std::abs(dy / dx) < 1.0) {
      // Front or back based on dx
      if (dx > 0) {
        // Back
        Point2D area_segment_start = {to.minX(), to.minY()};
        Point2D area_segment_end = {to.minX(), to.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      } else {
        // Front
        Point2D area_segment_start = {to.maxX(), to.minY()};
        Point2D area_segment_end = {to.maxX(), to.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      }
      double_t diff = intersection.second - to.minY();
      double_t inter_1 = to.minY() + std::floor(diff);
      double_t inter_2 = to.minY() + std::ceil(diff);
      next[0] = {intersection.first, inter_1};
      next[1] = {intersection.first, inter_2};
      // Calculate min distance
      min_dist = std::numeric_limits<double_t>::max();
      for (uint_fast32_t i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          double_t dist = find_distance(prev[i], next[j]);
          if (dist < min_dist) {
            min_dist = dist;
            area_end = prev[i];
            area_start = next[i];
          }
        }
      }
    } else if (std::abs(dy / dx) > 1.0) {
      // Top or bottom based on dy
      if (dy > 0) {
        // Bottom
        Point2D area_segment_start = {to.minX(), to.minY()};
        Point2D area_segment_end = {to.maxX(), to.minY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      } else {
        // Top
        Point2D area_segment_start = {to.minX(), to.maxY()};
        Point2D area_segment_end = {to.maxX(), to.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      }
      double_t diff = intersection.first - to.minX();
      double_t inter_1 = to.minX() + std::floor(diff);
      double_t inter_2 = to.minX() + std::ceil(diff);
      next[0] = {inter_1, intersection.second};
      next[1] = {inter_2, intersection.second};
      // Calculate min distance
      min_dist = std::numeric_limits<double_t>::max();
      for (uint_fast32_t i = 0; i < 2; ++i) {
        for (uint_fast32_t j = 0; j < 2; ++j) {
          double_t dist = find_distance(prev[i], next[j]);
          if (dist < min_dist) {
            min_dist = dist;
            area_end = prev[i];
            area_start = next[i];
          }
        }
      }
    } else {
      // Corner based on dx and dy
      if (dx > 0) {
        if (dy > 0) {
          // Bottom left corner
          area_start = {to.minX(), to.minY()};
        } else {
          // Top left corner
          area_start = {to.minX(), to.maxY()};
        }
      } else {
        if (dy > 0) {
          // Bottom right corner
          area_start = {to.maxX(), to.minY()};
        } else {
          // Top right corner
          area_start = {to.maxX(), to.maxY()};
        }
      }
    }
  } else {
    // Top or bottom based on dy
    // For point exiting "from" area
    if (dy > 0) {
      // Top
      Point2D area_segment_start = {from.minX(), from.maxY()};
      Point2D area_segment_end = {from.maxX(), from.maxY()};
      LineIntersect(line_start.first, line_start.second,
                    line_end.first, line_end.second,
                    area_segment_start.first, area_segment_start.second,
                    area_segment_end.first, area_segment_end.second,
                    intersection.first, intersection.second);
    } else {
      // Bottom
      Point2D area_segment_start = {from.minX(), from.minY()};
      Point2D area_segment_end = {from.maxX(), from.minY()};
      LineIntersect(line_start.first, line_start.second,
                    line_end.first, line_end.second,
                    area_segment_start.first, area_segment_start.second,
                    area_segment_end.first, area_segment_end.second,
                    intersection.first, intersection.second);
    }
    double_t diff = intersection.first - from.minX();
    double_t inter_1 = from.minX() + std::floor(diff);
    double_t inter_2 = from.minX() + std::ceil(diff);
    prev[0] = {inter_1, intersection.second};
    prev[1] = {inter_2, intersection.second};
    // For point entering "to" area
    if (dy > 0) {
      // Bottom
      Point2D area_segment_start = {to.minX(), to.minY()};
      Point2D area_segment_end = {to.maxX(), to.minY()};
      LineIntersect(line_start.first, line_start.second,
                    line_end.first, line_end.second,
                    area_segment_start.first, area_segment_start.second,
                    area_segment_end.first, area_segment_end.second,
                    intersection.first, intersection.second);
    } else {
      // Top
      Point2D area_segment_start = {to.minX(), to.maxY()};
      Point2D area_segment_end = {to.maxX(), to.maxY()};
      LineIntersect(line_start.first, line_start.second,
                    line_end.first, line_end.second,
                    area_segment_start.first, area_segment_start.second,
                    area_segment_end.first, area_segment_end.second,
                    intersection.first, intersection.second);
    }
    diff = intersection.first - to.minX();
    inter_1 = to.minX() + std::floor(diff);
    inter_2 = to.minX() + std::ceil(diff);
    next[0] = {inter_1, intersection.second};
    next[1] = {inter_2, intersection.second};
    // Calculate min distance
    min_dist = std::numeric_limits<double_t>::max();
    for (uint_fast32_t i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        double_t dist = find_distance(prev[i], next[j]);
        if (dist < min_dist) {
          min_dist = dist;
          area_end = prev[i];
          area_start = next[i];
        }
      }
    }
  }
  if (area_start.first < to.minX() || area_start.first > to.maxX()) {
    std::cout << "area_start x out of bounds" << std::endl;
    // TODO: This shouldn't happen. Fix it.
    if (area_start.first < to.minX()) {
      area_start.first = to.minX();
    } else {
      area_start.first = to.maxX();
    }
    std::cout << "(" <<from.origin().first << ", " << from.origin().second
              << ")->(" << to.origin().first << ", " << to.origin().second
              << ")" << std::endl;
  }
  if (area_start.second < to.minY() || area_start.second > to.maxY()) {
    std::cout << "area_start y out of bounds" << std::endl;
    // TODO: This shouldn't happen. Fix it.
    if (area_start.second < to.minY()) {
      area_start.second = to.minY();
    } else {
      area_start.second = to.maxY();
    }
    std::cout << "(" <<from.origin().first << ", " << from.origin().second
              << ")->(" << to.origin().first << ", " << to.origin().second
              << ")" << std::endl;
  }

  if (area_end.first < from.minX() || area_end.first > from.maxX()) {
    std::cout << "area_end x out of bounds" << std::endl;
    // TODO: This shouldn't happen. Fix it.
    if (area_end.first < from.minX()) {
      area_end.first = from.minX();
    } else {
      area_end.first = from.maxX();
    }
    std::cout << "(" <<from.origin().first << ", " << from.origin().second
              << ")->(" << to.origin().first << ", " << to.origin().second
              << ")" << std::endl;
  }
  if (area_end.second < from.minY() || area_end.second > from.maxY()) {
    std::cout << "area_end y out of bounds" << std::endl;
    // TODO: This shouldn't happen. Fix it.
    if (area_end.second < from.minY()) {
      area_end.second = from.minY();
    } else {
      area_end.second = from.maxY();
    }
    std::cout << "(" <<from.origin().first << ", " << from.origin().second
              << ")->(" << to.origin().first << ", " << to.origin().second
              << ")" << std::endl;
  }
  from.finish(area_end);
  to.start(area_start);
  return min_dist;
}

double_t getAreaEnd(Point2D point, SamplingArea &area) {
  Point2D intersection;
  Point2D area_start, area_end;
  Point2D prev[2];
  Point2D next[2];
  Point2D line_start = area.origin();
  Point2D line_end = point;
  double_t dx = line_end.first - line_start.first;
  double_t dy = line_end.second - line_start.second;

  if (dx != 0.0) {
    if (std::abs(dy / dx) < 1.0) {
      // Front or back based on dx
      if (dx > 0) {
        // Front
        Point2D area_segment_start = {area.maxX(), area.minY()};
        Point2D area_segment_end = {area.maxX(), area.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      } else {
        // Back
        Point2D area_segment_start = {area.minX(), area.minY()};
        Point2D area_segment_end = {area.minX(), area.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      }
      double_t diff = intersection.second - area.minY();
      double_t inter_1 = area.minY() + std::floor(diff);
      double_t inter_2 = area.minY() + std::ceil(diff);
      prev[0] = {intersection.first, inter_1};
      prev[1] = {intersection.first, inter_2};
      if (find_distance(prev[0], line_end)
          < find_distance(prev[1], line_end)) {
        area_end = prev[0];
      } else {
        area_end = prev[1];
      }
    } else if (std::abs(dy / dx) > 1.0) {
      // Top or bottom based on dy
      if (dy > 0) {
        // Top
        Point2D area_segment_start = {area.minX(), area.maxY()};
        Point2D area_segment_end = {area.maxX(), area.maxY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      } else {
        // Bottom
        Point2D area_segment_start = {area.minX(), area.minY()};
        Point2D area_segment_end = {area.maxX(), area.minY()};
        LineIntersect(line_start.first, line_start.second,
                      line_end.first, line_end.second,
                      area_segment_start.first, area_segment_start.second,
                      area_segment_end.first, area_segment_end.second,
                      intersection.first, intersection.second);
      }
      double_t diff = intersection.first - area.minX();
      double_t inter_1 = area.minX() + std::floor(diff);
      double_t inter_2 = area.minX() + std::ceil(diff);
      prev[0] = {inter_1, intersection.second};
      prev[1] = {inter_2, intersection.second};
      if (find_distance(prev[0], line_end)
          < find_distance(prev[1], line_end)) {
        area_end = prev[0];
      } else {
        area_end = prev[1];
      }
    } else {
      // Corner based on dx and dy
      if (dx > 0) {
        if (dy > 0) {
          // Top right
          area_end = {area.maxX(), area.maxY()};
        } else {
          // Bottom right
          area_end = {area.maxX(), area.minY()};
        }
      } else {
        if (dy > 0) {
          // Top left
          area_end = {area.minX(), area.maxY()};
        } else {
          // Bottom left
          area_end = {area.minX(), area.minY()};
        }
      }
    }
  } else {
    // Top or bottom based on dy
    if (dy > 0) {
      // Top
      Point2D area_segment_start = {area.minX(), area.maxY()};
      Point2D area_segment_end = {area.maxX(), area.maxY()};
      LineIntersect(line_start.first, line_start.second,
                    line_end.first, line_end.second,
                    area_segment_start.first, area_segment_start.second,
                    area_segment_end.first, area_segment_end.second,
                    intersection.first, intersection.second);
    } else {
      // Bottom
      Point2D area_segment_start = {area.minX(), area.minY()};
      Point2D area_segment_end = {area.maxX(), area.minY()};
      LineIntersect(line_start.first, line_start.second,
                    line_end.first, line_end.second,
                    area_segment_start.first, area_segment_start.second,
                    area_segment_end.first, area_segment_end.second,
                    intersection.first, intersection.second);
    }
    double_t diff = intersection.first - area.minX();
    double_t inter_1 = area.minX() + std::floor(diff);
    double_t inter_2 = area.minX() + std::ceil(diff);
    prev[0] = {inter_1, intersection.second};
    prev[1] = {inter_2, intersection.second};
    if (find_distance(prev[0], line_end)
        < find_distance(prev[1], line_end)) {
      area_end = prev[0];
    } else {
      area_end = prev[1];
    }
  }

  if (area_end.first < area.minX() || area_end.first > area.maxX()) {
    std::cout << "area_end x out of bounds" << std::endl;
    // TODO: This shouldn't happen. Fix it.
    if (area_end.first < area.minX()) {
      area_end.first = area.minX();
    } else {
      area_end.first = area.maxX();
    }
    std::cout << "(" <<area.origin().first << ", " << area.origin().second
              << ")->(" << point.first << ", " << point.second << ")"
              << std::endl;
  }
  if (area_end.second < area.minY() || area_end.second > area.maxY()) {
    std::cout << "area_end y out of bounds" << std::endl;
    // TODO: This shouldn't happen. Fix it.
    if (area_end.second < area.minY()) {
      area_end.second = area.minY();
    } else {
      area_end.second = area.maxY();
    }
    std::cout << "(" <<area.origin().first << ", " << area.origin().second
              << ")->(" << point.first << ", " << point.second << ")"
              << std::endl;
  }
  area.finish(area_end);
  return find_distance(area_end, line_end);
}

two_level_solution solveTwoLevel(Point2D start,
                                 Point2D finish,
                                 Vector<SamplingArea> &problem,
                                 uint_fast32_t num_vehicles,
                                 const Vector<double_t> &budget,
                                 bool add_area_cost,
                                 bool fair) {
  /**
   * 1. Solve high level.
   * 2. For each vehicle.
   * * 2.1. Compute available cost per area.
   * * 2.2. For each area solve 10 GA instances and get the best.
   * 3. Return solution.
   */
  two_level_solution ret;
  ret.total_time = 0.0;
  ret.total_utility = 0.0;

  std::random_device rd;
  std::mt19937 gen(rd());

  // Setup high level problem
  Vector<Point2D> depots;
  depots.reserve(num_vehicles + 1);
  depots.push_back(finish);
  for (size_t i = 0; i < num_vehicles; ++i) {
    depots.push_back(start);
  }

  Vector<Point2D> nodes;
  nodes.reserve(problem.size());
  for (const SamplingArea &a:problem) {
    nodes.push_back(a.origin());
  }

  Vector<Point2D> hl_problem = depots;
  hl_problem.insert(std::end(hl_problem), std::begin(nodes), std::end(nodes));

  Matrix<double_t> cost_mat;
  for (int i = 0; i < hl_problem.size(); i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < hl_problem.size(); ++i) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < hl_problem.size(); ++j) {
      double dist = find_distance(hl_problem[i], hl_problem[j]);

      cost_mat[i].push_back(dist);
      if (j >= depots.size() && add_area_cost) {
        cost_mat[i].back() += problem[j - depots.size()].nodes()->size();
      }
      cost_mat[j].push_back(dist);
      if (i >= depots.size() && add_area_cost) {
        cost_mat[j].back() += problem[i - depots.size()].nodes()->size();
      }
    }
  }
  auto timer_start = std::chrono::high_resolution_clock::now();
  Matrix<uint_fast32_t> tours;
  double_t success = runMDMTSP(
      num_vehicles, nodes.size(), depots.size(), cost_mat, fair, tours);

  if (success > 0.0) {
//    for (const Vector<uint_fast32_t> &tour: tours) {
//      print_vector(tour, std::cout);
//    }
    size_t tour_idx = 0;
    for (const Vector<uint_fast32_t> &tour: tours) {
      double_t travel_cost = 0.0;
      Vector<uint_fast32_t> path;
      path.push_back(0);

      Vector<SamplingArea> areas;
      Vector<size_t> indices;

      size_t total_nodes = 0;
      Vector<uint_fast32_t>::const_iterator idx_iter;
      for (idx_iter = tour.begin() + 1; idx_iter != tour.end(); ++idx_iter) {
        size_t prob_idx = *idx_iter - depots.size();
        if (idx_iter == tour.begin() + 1) {
          //First area. Start calculated based on global start.
          total_nodes +=
              problem[prob_idx].width()
                  * problem[prob_idx].length();
          travel_cost +=
              getAreaStart(start, problem[prob_idx]);
          indices.push_back(prob_idx);
        } else if (idx_iter == tour.end() - 1) {
          travel_cost +=
              getAreaEnd(finish, problem[*(idx_iter - 1) - depots.size()]);
          indices.push_back(prob_idx);
        } else {
          travel_cost +=
              getAreaStartEnd(problem[*(idx_iter - 1) - depots.size()],
                              problem[prob_idx]);
          total_nodes +=
              problem[prob_idx].width()
                  * problem[prob_idx].length();
          indices.push_back(prob_idx);
        }
      }

      double_t available_budget = budget[tour_idx] - travel_cost;
      std::cout << "Budget: " << budget[tour_idx] << " travel cost: " << travel_cost << std::endl;
      std::cout << "Available budget: " << available_budget << std::endl;
      double_t used_budget = 0.0;
      size_t area_idx = 0;
      size_t ga_iterations = 1;
      const uint_fast32_t pop_size = 250;
      const uint_fast32_t num_gen = 25;
      const uint_fast8_t tournament_size = 3;
      const double_t cx_rate = 0.5;
      const double_t mut_rate = 0.8;
      const double_t elite_rate = 0.17;
      for (idx_iter = tour.begin() + 1; idx_iter != tour.end() - 1;
           ++idx_iter) {
//      for (const SamplingArea &a:areas) {
        size_t prob_idx = *idx_iter - depots.size();
        size_t area_nodes_size =
            problem[prob_idx].width()
                * problem[prob_idx].length();
        double_t area_budget =
            available_budget * double_t(area_nodes_size)
                / double_t(total_nodes);
        Matrix<double_t> area_cost_mat;
        area_cost_mat.clear();
        std::shared_ptr<Vector<Point2D>>
            area_nodes = problem[prob_idx].nodes();
        if (area_budget > area_nodes->size()) {
          area_budget = area_nodes->size()+1; //Todo: check if helps
        }
        used_budget += area_budget;
//        std::cout << "Area budget: " << area_budget << std::endl;
        for (int i = 0; i < area_nodes->size(); i++) {
          std::vector<double> tmp_vec;
          area_cost_mat.push_back(tmp_vec);
        }

        for (size_t i = 0; i < area_nodes->size(); ++i) {
          area_cost_mat[i].push_back(0.0);
          for (size_t j = i + 1; j < area_nodes->size(); ++j) {
            double dist = find_distance(area_nodes->at(i), area_nodes->at(j));
            area_cost_mat[i].push_back(dist);
            area_cost_mat[j].push_back(dist);
          }
        }

        Vector<double_t> rewards(area_nodes->size(), 1);
//        std::cout << "{" << problem[prob_idx].start().first << ", " << problem[prob_idx].start().second << "} " << problem[prob_idx].idxFromPoint(problem[prob_idx].start()) << std::endl;
//        std::cout << "{" << problem[prob_idx].finish().first << ", " << problem[prob_idx].finish().second << "} " << problem[prob_idx].idxFromPoint(problem[prob_idx].finish()) << std::endl;
//        std::cout << "Start end cost: " << area_cost_mat[problem[prob_idx].idxFromPoint(problem[prob_idx].start())][problem[prob_idx].idxFromPoint(problem[prob_idx].finish())] << std::endl;
//        std::cout << "Area size: "
//                  << problem[prob_idx].width() * problem[prob_idx].length()
//                  << std::endl;
        cop_ga::Chromosome best;
        best.fitness = 0.0;
        best.cost = 1.0;
        for (size_t iter = 0; iter < ga_iterations; ++iter) {
          cop_ga::Chromosome c = cop_ga::ga_cop(
              area_cost_mat,
              rewards,
              area_budget,
              problem[prob_idx].idxFromPoint(problem[prob_idx].start()),
              problem[prob_idx].idxFromPoint(problem[prob_idx].finish()),
              gen, pop_size, num_gen, tournament_size, cx_rate, mut_rate,
              elite_rate);
          if (c.fitness > best.fitness) {
            best = c;
          }
        }
//        std::cout << best.cost << " " << best.fitness << " "
//                  << std::cbrt(best.fitness * best.cost) << std::endl;
//        if (std::cbrt(best.fitness * best.cost)
//            > problem[prob_idx].width() * problem[prob_idx].length()) {
//          print_vector(best.path);
//        }

//        path.resize(path.size()+best.path.size());
        size_t prev_area_size = 0;
        for (size_t i = 0; i < indices[area_idx]; ++i) {
          prev_area_size += problem[i].nodes()->size();
        }
        for (const uint_fast32_t &idx:best.path) {
          path.push_back(prev_area_size + idx + 1);
        }
        ret.total_utility += cbrt(best.fitness * best.cost);
        ++area_idx;
      }
      size_t total_problem_size = 0;
      for (const SamplingArea &a:problem) {
        total_problem_size += a.nodes()->size();
      }
      path.push_back(total_problem_size + 1);
      ret.paths.push_back(path);
      ++tour_idx;
      std::cout << "Used budget: " << used_budget << std::endl;
    }
    auto timer_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = timer_end - timer_start;
    ret.total_time = diff.count();
  } else {
    std::cout << "Problem was infeasible" << std::endl;
  }
  return ret;
}
