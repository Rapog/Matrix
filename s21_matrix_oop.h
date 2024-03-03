#ifndef S21_MATRIX_OOP_S21_MATRIX_OOP_H_
#define S21_MATRIX_OOP_S21_MATRIX_OOP_H_

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iterator>

namespace S21 {
class S21Matrix {
 public:
  S21Matrix() noexcept;
  explicit S21Matrix(size_t rows, size_t cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other) noexcept;
  ~S21Matrix();
  //
  bool EqMatrix(const S21Matrix& other) const noexcept;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num) noexcept;
  void MulMatrix(const S21Matrix& other);
  //
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix operator*(const double num);
  S21Matrix& operator*=(const double num);
  friend S21Matrix operator*(const double num, S21Matrix& other);
  double& operator()(const size_t i, const size_t j);
  const double& operator()(const size_t i, const size_t j) const;
  bool operator==(const S21Matrix& other) const noexcept;
  bool operator!=(const S21Matrix& other) const noexcept;

  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  S21Matrix Transform(const S21Matrix& other, double (*foo)(double, double));
  S21Matrix Transform(double (*foo)(double));

  size_t GetRows() const noexcept { return rows_; }
  size_t GetCols() const noexcept { return cols_; }
  void Resize(size_t newrows, size_t newcols);
  double GetNum(size_t i, size_t j) const noexcept { return matrix_[i][j]; }

  template <typename InputIt, typename Func>
  static void ForEach(InputIt first, InputIt last, Func foo) {
    for (InputIt it = first; it != last; ++it) {
      foo(*it);
    }
  }

  template <typename T>
  static typename std::remove_reference<T>::type&& move(T&& t) {
    return static_cast<typename std::remove_reference<T>::type&&>(t);
  }

  double* begin() { return matrix_[0]; }
  double* end() { return matrix_[0] + rows_ * cols_; }
  const double* begin() const { return matrix_[0]; }
  const double* end() const { return matrix_[0] + rows_ * cols_; }
  // void PrintMatrix() noexcept;

  friend std::ostream& operator<<(std::ostream& out, const S21Matrix& other);
  friend std::istream& operator>>(std::istream& in, S21Matrix& other);

 private:
  size_t rows_;
  size_t cols_;
  double** matrix_;

 private:
  int SwapLines(size_t row, size_t col);
  void GaussDiv(size_t row, size_t col);
  void GaussSubLines(size_t row, size_t col);
  // void CutMatxForMinor(S21Matrix* big,size_t row, size_t col);
  // double DetForCalcComp(size_t row, size_t col);
};

};  // namespace S21

#endif  // S21_MATRIX_OOP_S21_MATRIX_OOP_H_