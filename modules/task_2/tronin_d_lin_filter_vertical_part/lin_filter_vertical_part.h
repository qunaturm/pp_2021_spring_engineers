// Copyright 2020 Tronin Dmitry
#ifndef MODULES_TASK_2_TRONIN_D_LIN_FILTER_VERTICAL_PART_LIN_FILTER_VERTICAL_PART_H_
#define MODULES_TASK_2_TRONIN_D_LIN_FILTER_VERTICAL_PART_LIN_FILTER_VERTICAL_PART_H_

#include <cstdint>
#include <vector>

using std::size_t;

std::vector<double> CalculateGaussFilterP(size_t size, double sigma);
std::vector<uint8_t> ApplyFilterParallel(const std::vector<double> &filter,
                                 const std::vector<uint8_t> &image,
                                 size_t width,
                                 size_t height,
                                 size_t filter_size);

template<class T>
inline T ClampP(T value, T low, T high) {
    if (value < low) {
        return low;
    }
    if (value > high) {
        return high;
    }
    return value;
}

#endif  // MODULES_TASK_2_TRONIN_D_LIN_FILTER_VERTICAL_PART_LIN_FILTER_VERTICAL_PART_H_
