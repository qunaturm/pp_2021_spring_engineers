// Copyright 2020 Tronin Dmitry
#include "../../modules/task_2/tronin_d_lin_filter_vertical_part/lin_filter_vertical_part.h"
#include <omp.h>
#include <cmath>
#include <stdexcept>
#include <iostream>

static const double kPi = 3.141592653589793;

std::vector<double> CalculateGaussFilterP(size_t size, double sigma) {
    if (size % 2 == 0) {
        throw std::invalid_argument("Size can't be even");
    }
    if (size <= 0) {
        throw std::invalid_argument("Size should be positive");
    }
    if (sigma <= 0) {
        throw std::invalid_argument("Sigma should be positive");
    }
    double sum = 0;

    std::vector<double> filter(size * size, 0);
    for (int row = -static_cast<int>(size) / 2; row <= static_cast<int>(size) / 2; ++row) {
        for (int col = -static_cast<int>(size) / 2; col <= static_cast<int>(size) / 2; ++col) {
            filter[(row + static_cast<int>(size) / 2) * size + col + static_cast<int>(size) / 2] =
                1. / (2 * kPi * sigma * sigma)
                    * exp(-static_cast<double>(row * row + col * col) / (2 * sigma * sigma));
            sum += filter[(row + static_cast<int>(size) / 2) * size + col
                + static_cast<int>(size) / 2];
        }
    }

    double a = 0;

    for (int row = -static_cast<int>(size) / 2; row <= static_cast<int>(size) / 2; ++row) {
        for (int col = -static_cast<int>(size) / 2; col <= static_cast<int>(size) / 2; ++col) {
            filter[(row + static_cast<int>(size) / 2) * size + col + static_cast<int>(size) / 2] /=
                sum;
            a += filter[(row + static_cast<int>(size) / 2) * size + col
                + static_cast<int>(size) / 2];
        }
    }

    return filter;
}

std::vector<uint8_t> ApplyFilterParallel(const std::vector<double> &filter,
                                         const std::vector<uint8_t> &image,
                                         size_t width,
                                         size_t height,
                                         size_t filter_size) {
    if (image.size() != width * height) {
        throw std::invalid_argument("Incorrect image size");
    }
    if (filter.size() != filter_size * filter_size) {
        throw std::invalid_argument("Incorrect filter size");
    }
    if (width <= 0 || height <= 0 || filter_size <= 0) {
        throw std::invalid_argument("Arguments should be positive");
    }
    std::vector<uint8_t> result_image(image.size(), 1);

    double pixel_value, test;
    int image_row, image_col;

#pragma omp parallel default(none) shared(image, width, height, filter_size, filter, result_image, std::cout)
    {
#pragma omp for private(pixel_value, test, image_row, image_col)
        for (int col = 0; col < static_cast<int>(width); ++col) {
            for (int row = 0; row < static_cast<int>(height); ++row) {
                pixel_value = 0;
                for (int filter_row = -static_cast<int>(filter_size) / 2;
                     filter_row <= static_cast<int>(filter_size) / 2; ++filter_row) {
                    for (int filter_col = -static_cast<int>(filter_size) / 2;
                         filter_col <= static_cast<int>(filter_size) / 2; ++filter_col) {
                        image_row =
                            static_cast<size_t>(ClampP(static_cast<int>(row) + filter_row,
                                                       0,
                                                       static_cast<int>(height) - 1));
                        image_col =
                            static_cast<size_t>(ClampP(static_cast<int>(col) + filter_col,
                                                       0,
                                                       static_cast<int>(width) - 1));
                        test = image[image_row * width + image_col]
                            * filter[
                                (filter_row + static_cast<int>(filter_size) / 2) * filter_size
                                    + filter_col + static_cast<int>(filter_size) / 2];
                        pixel_value += test;
                    }
                }
                result_image[row * width + col] =
                    static_cast<uint8_t>(ClampP(pixel_value, 0., 255.));
            }
        }
    }

    return result_image;
}
