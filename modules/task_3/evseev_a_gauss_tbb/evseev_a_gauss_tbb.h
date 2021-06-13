// Copyright 2021 Evseev Alexander
#ifndef MODULES_TASK_3_EVSEEV_A_GAUSS_TBB_EVSEEV_A_GAUSS_TBB_H_
#define MODULES_TASK_3_EVSEEV_A_GAUSS_TBB_EVSEEV_A_GAUSS_TBB_H_


#include <tbb/parallel_for.h>
#include <tbb/tick_count.h>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

using Matrix = std::vector<std::vector<double>>;
Matrix RandMatrix(int height, int width);
int clamp(int value, int max, int min);
Matrix GaussKernel(int R, double sigma);
Matrix SeqGauss(const Matrix& mat, int height, int width, int R, double sigma);
Matrix TbbGauss(const Matrix& mat, int height, int width, int R, double sigma);
void printMatrix(Matrix mat, int height, int width);

#endif  // MODULES_TASK_3_EVSEEV_A_GAUSS_TBB_EVSEEV_A_GAUSS_TBB_H_
