// Copyright 2021 Elandaev Pavel

#include "../../modules/task_2/elandaev_p_integral_rectangle_omp/rectangle.h"
#include <omp.h>

// double RectangleIntegral(double start, double end, std::function<double(double)> f, double step) {
//    if (end < start)
//        throw "wrong interval";
//    if (step <= 0)
//        throw "wrong splitting: smallest";
//    if (step > end - start)
//        throw "step is biggest";
//    int countstep = static_cast<int>(end - start) / step;
//    double x;
//    double sum = 0.0;
//    for (double i = 0.0; i < countstep; i ++) {
//        x = start + i * step;
//        sum += f(x);
//    }
//    return (step * sum);
// }
void iterplus(std::vector<int> *B, int it, const std::vector<std::vector<double>> &p) {
    int size = p[it].size();
    if (B->operator[](it) == (size - 1)) {
        if ((it-1) < 0)
            return;
        B->operator[](it) = 0;
        iterplus(B, --it, p);
    } else {
        B->operator[](it)++;
    }
}
double RecInt(std::vector<double> start,
              std::vector<double> end,
              std::function<double(std::vector<double>)> f,
              double step) {
    std::vector<int> countstep;
    int size = start.size();
    for (int i = 0; i < size; i++) {
        if (end[i] < start[i])
            throw "wrong segment";
        if (step > end[i] - start[i])
            throw "step is biggest";
        countstep.push_back(static_cast<int>((end[i] - start[i]) / step));
    }
    if (step <= 0)
        throw "wrong splitting: smallest";
    std::vector<std::vector<double> > point;
    point.resize(start.size());
    int countTrial = 1;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < countstep[i]; j++) {
            point[i].push_back(start[i] + j * step);
        }
        countTrial *= point[i].size();
    }
    double sum = 0.0;
    std::vector<int> B(start.size());
    for (int i = 0; i < size; i++) {
        B[i] = 0;
    }
    std::vector<double> Trial(start.size());
    int dim = static_cast<int>(start.size() - 1);
    for (int i = 0; i < countTrial; i++) {
        for (int j = 0; j < size; j++) {
            Trial[j] = point[j][B[j]];
        }
        iterplus(&B, dim, point);
        sum += f(Trial);
    }
    return sum * (std::pow(step, dim+1));
}

double RecIntOmp(std::vector<double> start,
              std::vector<double> end,
              std::function<double(std::vector<double>)> f,
              double step) {
    std::vector<int> countstep;
    int size = start.size();
    for (int i = 0; i < size; i++) {
        if (end[i] < start[i])
            throw "wrong segment";
        if (step > end[i] - start[i])
            throw "step is biggest";
        countstep.push_back(static_cast<int>((end[i] - start[i]) / step));
    }
    if (step <= 0)
        throw "wrong splitting: smallest";
    std::vector<std::vector<double> > point;
    point.resize(start.size());
    int countTrial = 1;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < countstep[i]; j++) {
            point[i].push_back(start[i] + j * step);
        }
        countTrial *= point[i].size();
    }
    double sum = 0.0;
    std::vector<int> B(size);
    for (int i = 0; i < size; i++) {
        B[i] = 0;
    }
    std::vector<double> Trial(size);
    int dim = size - 1;
    int flag = 1;
#pragma omp parallel for schedule(static) reduction(+:sum) firstprivate(Trial, B, flag) shared(point, dim)
    for (int i = 0; i < countTrial; i++) {
        if (flag == 1) {
            int k = i;
            int count = 0;
            while (k > 0) {
                B[dim-count] = k % point[dim - count].size();
                k = k / point[dim - count].size();
                count++;
            }
            flag = 0;
        }
        for (int j = 0; j < size; j++) {
            Trial[j] = point[j][B[j]];
        }
        iterplus(&B, dim, point);
        sum += f(Trial);
    }
    return sum * (std::pow(step, size));
}

double RecIntT(std::vector<double> start,
    std::vector<double> end,
    std::function<double(std::vector<double>)> f,
    double step) {
    std::vector<int> countstep;
    int size = start.size();
    for (int i = 0; i < size; i++) {
        if (end[i] < start[i])
            throw "wrong segment";
        if (step > end[i] - start[i])
            throw "step is biggest";
        countstep.push_back(static_cast<int>((end[i] - start[i]) / step));
    }
    if (step <= 0)
        throw "wrong splitting: smallest";
    std::vector<std::vector<double> > point;
    point.resize(start.size());
    int countTrial = 1;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < countstep[i]; j++) {
            point[i].push_back(start[i] + j * step);
        }
        countTrial *= point[i].size();
    }
    double sum = 0.0;
    std::vector<int> B(start.size());
    std::vector<double> Trial(start.size());
    int dim = static_cast<int>(start.size() - 1);
    for (int i = 0; i < countTrial; i++) {
        for (int j = 0; j < size; j++) {
            int k = i;
            while (k > 0) {
                int count = 0;
                B[dim - count] = k % point[dim - count].size();
                k = k / point[dim - count].size();
                count++;
            }
            Trial[j] = point[j][B[j]];
        }
        sum += f(Trial);
    }
    return sum * (std::pow(step, dim + 1));
}
