// Copyright 2021 Khismatulina Karina
#include "../../../modules/task_2/khismatulina_k_gradient/omp.h"

#include <omp.h>

#include <cassert>
#include <iostream>
#include <random>
#include <vector>

double MyAbs(double x) {
  if (x < 0) return -x;
  return x;
}

std::vector<double> getRandomVectorOMP(int size) {
  if (size < 0) throw("size ouhgt to be > 0");
  std::random_device rd;
  std::mt19937 mersenne(rd());
  std::vector<double> vec(size, 1);
  for (int i = 0; i < size; ++i) {
    vec[i] = mersenne() * mersenne() % 7;
  }
  return vec;
}

std::vector<double> getRandomMatrixOMP(int size) {
  assert(size > 0);
  std::random_device rd;
  std::mt19937 mersenne(rd());
  std::vector<double> matrix(size * size, 1);
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      matrix[i * size + j] = matrix[j * size + i] = mersenne() * mersenne() % 7;
    }
  }
  return matrix;
}

double multVVOMP(std::vector<double> A, std::vector<double> B) {
  assert(A.size() == B.size());
  double result = 0;
  int size_A = A.size();
#pragma omp parallel for reduction(+ : result)
  for (int i = 0; i < size_A; ++i) {
    result += A[i] * B[i];
  }
  return result;
}

std::vector<double> multMVOMP(std::vector<double> m, std::vector<double> v) {
  std::vector<double> result(m.size() / v.size());
  int size_r = result.size();
#pragma omp parallel for
  for (int i = 0; i < size_r; ++i) {
    for (int j = 0; j < size_r; ++j) {
      result[i] += m[i * v.size() + j] * v[j];
    }
  }
  return result;
}

bool gradientParOMP(const std::vector<double>& matrix,
                    const std::vector<double>& vector, int size) {
  assert(size > 0);
  int iters = 0;
  double eps = 0.1, beta = 0.0, alpha = 0.0, check = 0.0;
  int i, k = 0;
  std::vector<double> res(size, 1);
  std::vector<double> Ah = multMVOMP(matrix, res);
  std::vector<double> discrepancyCurrent(size), discrepancyNext(size);
#pragma omp parallel for private(i)
  for (i = 0; i < size; ++i) {
    discrepancyCurrent[i] = vector[i] - Ah[i];
  }
  std::vector<double> h(discrepancyCurrent);
  do {
    iters++;
    Ah = multMVOMP(matrix, h);
    alpha = multVVOMP(discrepancyCurrent, discrepancyNext) / multVVOMP(h, Ah);
#pragma omp parallel for private(i)
    for (i = 0; i < size; ++i) {
      res[i] += alpha * h[i];
      discrepancyNext[i] = discrepancyCurrent[i] - alpha * Ah[i];
    }
    beta = multVVOMP(discrepancyNext, discrepancyNext) /
           multVVOMP(discrepancyCurrent, discrepancyCurrent);
    check = sqrt(multVVOMP(discrepancyNext, discrepancyNext));
#pragma omp parallel for private(k)
    for (k = 0; k < size; ++k) {
      h[k] = discrepancyNext[k] + beta * h[k];
    }
    std::vector<double> swap = discrepancyNext;
    discrepancyNext = discrepancyCurrent;
    discrepancyCurrent = swap;
  } while ((check > eps) && (iters <= size));

  std::vector<double> A_check(size);
  int correct = 0;
  int F = size / 10;
  if (size / 10 == 0) {
    F = 1;
  }
  double fail = 0.1 * F;
  for (int i = 0; i < size; ++i) {
    A_check[i] = 0;
    for (int j = 0; j < size; ++j) {
      A_check[i] += matrix[i * size + j] * res[j];
    }
    if ((MyAbs(A_check[i]) > MyAbs(vector[i]) &&
         (MyAbs(A_check[i]) - MyAbs(vector[i]) > fail)) ||
        ((MyAbs(A_check[i]) <= MyAbs(vector[i])) &&
         (MyAbs(vector[i]) - MyAbs(A_check[i]) > fail))) {
      correct = 1;
      break;
    }
  }
  return correct;
}
