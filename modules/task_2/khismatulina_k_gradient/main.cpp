// Copyright 2021 Khismatulina Karina
#include <gtest/gtest.h>
#include <omp.h>
#include <iostream>
#include <algorithm>
#include "../../../modules/task_2/khismatulina_k_gradient/omp.h"

TEST(gradient_omp, create_vec) {
    std::vector<double> vec = getRandomVectorOMP(1);
    int check = 1;
    ASSERT_EQ(check, vec.size());
}

TEST(gradient_omp, create_matrix) {
    std::vector<double> matrix = getRandomMatrixOMP(2);
    ASSERT_EQ(4, matrix.size());
}

TEST(gradient_omp, gradient_1_proc_4x4) {
    std::vector<double> s = getRandomMatrixOMP(2);
    std::vector<double> v = getRandomVectorOMP(4);
    bool resPar = gradientParOMP(s, v, 2, 1);
    bool check = 1;
    ASSERT_EQ(check, resPar);
}

TEST(gradient_omp, gradient_5_proc_9x9) {
    std::vector<double> s = getRandomMatrixOMP(3);
    std::vector<double> v = getRandomVectorOMP(9);
    bool resPar = gradientParOMP(s, v, 3, 5);
    bool check = 1;
    ASSERT_EQ(check, resPar);
}

TEST(gradient_omp, gradient_9_proc_16x16) {
    std::vector<double> s = getRandomMatrixOMP(4);
     std::vector<double> p = s;
    std::vector<double> v = getRandomVectorOMP(16);
    bool resPar = gradientParOMP(p, v, 4, 9);
    bool check = 1;
    ASSERT_EQ(check, resPar);
}

TEST(gradient_omp, gradient_16_proc_2500x2500) {
    std::vector<double> s = getRandomMatrixOMP(10);
    std::vector<double> v = getRandomVectorOMP(100);
    bool resPar = gradientParOMP(s, v, 10, 16);
    bool check = 1;
    ASSERT_EQ(check, resPar);
}

TEST(gradient_omp, gradient_16_proc_many_items) {
    std::vector<double> s = getRandomMatrixOMP(1000);
    std::vector<double> v = getRandomVectorOMP(1000000);
    bool resPar = gradientParOMP(s, v, 1000, 4);
    bool check = 1;
    ASSERT_EQ(check, resPar);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
