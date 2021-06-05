// Copyright 2021 Khismatulina Karina
#include <gtest/gtest.h>
#include <omp.h>

#include <algorithm>
#include <iostream>

#include "../../../modules/task_2/khismatulina_k_gradient/omp.h"

TEST(gradient_omp, create_vec) {
  std::vector<double> vec = getRandomVectorOMP(1);
  size_t check = 1;
  ASSERT_EQ(check, vec.size());
}

TEST(gradient_omp, create_matrix) {
  std::vector<double> matrix = getRandomMatrixOMP(2);
  size_t check = 4;
  ASSERT_EQ(check, matrix.size());
}

TEST(gradient_omp, gradient_5_proc_9x9) {
  std::vector<double> s = getRandomMatrixOMP(3);
  std::vector<double> v = getRandomVectorOMP(9);
  bool resPar = gradientParOMP(s, v, 3);
  bool check = 1;
  ASSERT_EQ(check, resPar);
}

TEST(gradient_omp, gradient_9_proc_16x16) {
  std::vector<double> s = getRandomMatrixOMP(4);
  std::vector<double> p = s;
  std::vector<double> v = getRandomVectorOMP(16);
  bool resPar = gradientParOMP(p, v, 4);
  bool check = 1;
  ASSERT_EQ(check, resPar);
}

TEST(gradient_omp, gradient_16_proc_2500x2500) {
  std::vector<double> s = getRandomMatrixOMP(10);
  std::vector<double> v = getRandomVectorOMP(100);
  bool resPar = gradientParOMP(s, v, 10);
  bool check = 1;
  ASSERT_EQ(check, resPar);
}

TEST(gradient_omp, gradient_25) {
  std::vector<double> s = getRandomMatrixOMP(25);
  std::vector<double> v = getRandomVectorOMP(625);
  bool resPar = gradientParOMP(s, v, 25);
  bool check = 1;
  ASSERT_EQ(check, resPar);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
