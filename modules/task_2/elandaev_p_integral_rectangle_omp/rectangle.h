// Copyright 2021 Elandaev Pavel
#ifndef MODULES_TASK_2_ELANDAEV_P_INTEGRAL_RECTANGLE_OMP_RECTANGLE_H_
#define MODULES_TASK_2_ELANDAEV_P_INTEGRAL_RECTANGLE_OMP_RECTANGLE_H_

#include <utility>
#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

void iterplus(std::vector<int> *B, int it, const std::vector<std::vector<int>> &p);
double RecInt(std::vector<double> start,
              std::vector<double> end,
              std::function<double(std::vector<double>)> f,
              double step);
double RecIntOmp(std::vector<double> start,
              std::vector<double> end,
              std::function<double(std::vector<double>)> f,
              double step);

double RecIntT(std::vector<double> start,
    std::vector<double> end,
    std::function<double(std::vector<double>)> f,
    double step);
double RecIntOmpT(std::vector<double> start,
    std::vector<double> end,
    std::function<double(std::vector<double>)> f,
    double step);
#endif  // MODULES_TASK_2_ELANDAEV_P_INTEGRAL_RECTANGLE_OMP_RECTANGLE_H_
