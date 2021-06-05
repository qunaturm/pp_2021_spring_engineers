// Copyright 2021 Khismatulina Karina
#ifndef MODULES_TASK_2_KHISMATULINA_K_GRADIENT_OMP_H_
#define MODULES_TASK_2_KHISMATULINA_K_GRADIENT_OMP_H_
#include <vector>

std::vector<double> getRandomVectorOMP(int size);
std::vector<double> getRandomMatrixOMP(int size);
double multVVOMP(std::vector<double> A, std::vector<double> B);
std::vector<double> multMVOMP(std::vector<double> m, std::vector<double> v);
bool gradientParOMP(const std::vector<double>& matrix, const std::vector<double>& vector, int size, int proc);

#endif  // MODULES_TASK_2_KHISMATULINA_K_GRADIENT_OMP_H_
