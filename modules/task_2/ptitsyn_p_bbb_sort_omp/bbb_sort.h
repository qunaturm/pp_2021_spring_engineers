// Copyright 2021 Ptitsyn Pavel
#ifndef MODULES_TASK_2_PTITSYN_P_BBB_SORT_OMP_BBB_SORT_H_
#define MODULES_TASK_2_PTITSYN_P_BBB_SORT_OMP_BBB_SORT_H_

#include <bitset>
#include <random>
#include <vector>

class RandomDouble {
 public:
  static double Next();

 private:
  RandomDouble();
  std::mt19937_64 gen_;
  std::uniform_real_distribution<double> dist_;
};

std::vector<double> random_double_array(size_t size);

std::vector<double> transit(std::vector<double> inputNumbers,
                            std::vector<double> loc, const int size, int add);
std::vector<double> loc_sort(std::vector<double> inputNumbers, const int size);
std::vector<double> BatcherSort(const std::vector<double>& vec);
std::vector<double> uns(const std::vector<double>& vec, int l, int r);
std::vector<double> sh(const std::vector<double>& vec, int l, int r);
std::vector<double> BatcherSortOMP(const std::vector<double>& vec);
std::vector<double> unsOMP(const std::vector<double>& vec, int l, int r);

#endif  // MODULES_TASK_2_PTITSYN_P_BBB_SORT_OMP_BBB_SORT_H_
