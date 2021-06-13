// Copyright 2021 Ptitsyn Pavel
#include "../../../modules/task_2/ptitsyn_p_bbb_sort_omp/bbb_sort.h"

#include <math.h>
#include <omp.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

double RandomDouble::Next() {
  static RandomDouble rand = RandomDouble();
  return rand.dist_(rand.gen_);
}

RandomDouble::RandomDouble() : gen_(std::random_device()()), dist_(0, 0e3) {}

std::vector<double> random_double_array(size_t size) {
  std::vector<double> arr(size);

  for (size_t i = 0; i < size; ++i) {
    arr[i] = RandomDouble::Next();
  }

  return arr;
}

std::vector<double> transit(std::vector<double> inputNumbers,
                            std::vector<double> loc, const int size, int add) {
  unsigned char* pmem = reinterpret_cast<unsigned char*>(inputNumbers.data());
  int counters[256] = {0};
  int sum = 0;

  for (int i = 0; i < size; i++) {
    counters[pmem[8 * i + add]]++;
  }
  for (int i = 0; i < 256; i++) {
    int tmp = counters[i];
    counters[i] = sum;
    sum += tmp;
  }
  for (int i = 0; i < size; i++) {
    int index = 8 * i + add;
    loc[counters[pmem[index]]] = inputNumbers[i];
    counters[pmem[index]]++;
  }
  return loc;
}

std::vector<double> loc_sort(std::vector<double> inputNumbers, const int size) {
  std::vector<double> loc(size, 0);
  std::vector<double> tmp(size, 0);

  for (int i = 0; i < 8; i++) {
    loc = transit(inputNumbers, loc, size, i);
    tmp = inputNumbers;
    inputNumbers = loc;
    loc = tmp;
  }
  return inputNumbers;
}

std::vector<double> BatcherSort(const std::vector<double>& vec) {
  int size = vec.size();

  if (vec.size() == 2) loc_sort(vec, size);
  if (vec.size() == 1) return vec;
  std::vector<double> vecAuxLeft, vecAuxRight, res;
  res = vec;
  vecAuxLeft.reserve(vec.size() / 2);
  vecAuxRight.reserve(vec.size() / 2);
  res = unsOMP(res, 0, size);
  vecAuxLeft.insert(vecAuxLeft.end(), res.begin(),
                    res.begin() + res.size() / 2);
  vecAuxRight.insert(vecAuxRight.end(), res.begin() + res.size() / 2,
                     res.end());
  vecAuxLeft = BatcherSort(vecAuxLeft);
  vecAuxRight = BatcherSort(vecAuxRight);
  res.erase(res.begin(), res.end());
  res.shrink_to_fit();
  res.insert(res.end(), vecAuxLeft.begin(), vecAuxLeft.end());
  res.insert(res.end(), vecAuxRight.begin(), vecAuxRight.end());
  res = sh(res, 0, size);
  return res;
}

std::vector<double> uns(const std::vector<double>& vec, int l, int r) {
  std::vector<double> auxLeft, auxRight, res;
  res = vec;
  if (res.size() % 2 != 0) {
    int lValue = res.at(res.size() - 1), lPos = res.size() - 1;

    res.erase(res.end() - 1);

    for (int j = l; j < r - 1; j += 2) {
      auxLeft.push_back(res.at(j));
      auxRight.push_back(res.at(j + 1));
    }

    if (lPos % 2 == 0)
      auxLeft.push_back(lValue);
    else
      auxRight.push_back(lValue);
  } else {
    for (int j = l; j < r; j += 2) {
      auxLeft.push_back(res.at(j));
      auxRight.push_back(res.at(j + 1));
    }
  }

  res.erase(res.begin(), res.end());
  res.shrink_to_fit();

  res.reserve(auxLeft.size() + auxRight.size());
  res.insert(res.end(), auxLeft.begin(), auxLeft.end());
  res.insert(res.end(), auxRight.begin(), auxRight.end());

  return res;
}

std::vector<double> sh(const std::vector<double>& vec, int l, int r) {
  std::vector<double> aux, res;
  res = vec;
  int m = (l + r) / 2;

  if (res.size() % 2 != 0) {
    int lValue = res.at(res.size() - 1);

    res.erase(res.end() - 1);

    for (int j = l; j < (r - 1) / 2; j++) {
      aux.push_back(res.at(j));
      aux.push_back(res.at(j + m));
    }

    aux.push_back(lValue);
  } else {
    for (int j = l; j < r / 2; j++) {
      aux.push_back(res.at(j));
      aux.push_back(res.at(j + m));
    }
  }

  res.erase(res.begin(), res.end());
  res.shrink_to_fit();

  for (int i = l; i < r; i++) {
    res.push_back(aux[i]);
  }
  res = loc_sort(res, res.size());

  return res;
}

std::vector<double> BatcherSortOMP(const std::vector<double>& vec) {
  std::vector<double> res(vec);
  std::vector<double> auxLeft, auxRight;
  int size = vec.size();
  int r = size / 2;
  int i;

  if (vec.size() == 2) loc_sort(res, size);
  if (vec.size() == 1) return vec;

  res = unsOMP(res, 0, size);

#pragma omp parallel sections num_threads(1)
  {
#pragma omp section
    { auxLeft = loc_sort(res, size); }
#pragma omp section
    { auxRight = loc_sort(res, size); }
  }

#pragma omp parallel shared(res)
  {
#pragma omp for
    for (i = 0; i < r; i++) {
      res[i] = auxLeft[i];
    }
  }

#pragma omp parallel shared(res)
  {
#pragma omp for
    for (i = r; i < size; i++) {
      res[i] = auxRight[i];
    }
  }

  res = sh(res, 0, size);

  return res;
}

std::vector<double> unsOMP(const std::vector<double>& vec, int l, int r) {
  std::vector<double> auxLeft, auxRight, res;
  res = vec;
  if (res.size() % 2 != 0) {
    int lValue = res.at(res.size() - 1), lPos = res.size() - 1;

    res.erase(res.end() - 1);
    for (int j = l; j < r - 1; j += 2) {
      auxLeft.push_back(res.at(j));
      auxRight.push_back(res.at(j + 1));
    }

    if (lPos % 2 == 0)
      auxLeft.push_back(lValue);
    else
      auxRight.push_back(lValue);
  } else {
    for (int j = l; j < r; j += 2) {
      auxLeft.push_back(res.at(j));
      auxRight.push_back(res.at(j + 1));
    }
  }

  res.erase(res.begin(), res.end());
  res.shrink_to_fit();

  res.reserve(auxLeft.size() + auxRight.size());
  res.insert(res.end(), auxLeft.begin(), auxLeft.end());
  res.insert(res.end(), auxRight.begin(), auxRight.end());

  return res;
}
