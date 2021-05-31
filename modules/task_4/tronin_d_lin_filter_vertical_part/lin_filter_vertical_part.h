// Copyright 2021 Tronin Dmitry
#ifndef MODULES_TASK_4_TRONIN_D_LIN_FILTER_VERTICAL_PART_LIN_FILTER_VERTICAL_PART_H_
#define MODULES_TASK_4_TRONIN_D_LIN_FILTER_VERTICAL_PART_LIN_FILTER_VERTICAL_PART_H_

#include <stddef.h>
#include <cstdint>
#include <vector>

std::vector<double> CalculateGaussFilterSTD(size_t size, double sigma);
std::vector<uint8_t> ApplyFilterSTD(const std::vector<double> &filter,
                                    const std::vector<uint8_t> &image,
                                    size_t width,
                                    size_t height,
                                    size_t filter_size, int number_of_threads);

template<class T>
inline T ClampSTD(T value, T low, T high) {
    if (value < low) {
        return low;
    }
    if (value > high) {
        return high;
    }
    return value;
}

#endif  // MODULES_TASK_4_TRONIN_D_LIN_FILTER_VERTICAL_PART_LIN_FILTER_VERTICAL_PART_H_
