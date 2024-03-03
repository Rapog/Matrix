#include "s21_matrix_oop.h"

#include <vector>

namespace S21 {
S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0), matrix_(nullptr) {}
S21Matrix::S21Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
  if (!rows_ && !cols_) {
    matrix_ = nullptr;
  } else if (!rows_ || !cols_) {
    throw std::logic_error("Creation error: Invalid matrix dimensions");
  } else {
    matrix_ = new double*[rows_];
    matrix_[0] = new double[rows_ * cols_]();
    for (size_t i = 1; i < rows_; ++i) {
      matrix_[i] = matrix_[0] + cols_ * i;
    }
  }
}
S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  if (matrix_)
    std::memcpy(matrix_[0], other.matrix_[0], sizeof(double) * rows_ * cols_);
}
S21Matrix::S21Matrix(S21Matrix&& other) noexcept : S21Matrix() {
  // rows_(other.rows_);
  // cols_(other.cols_);
  // matrix_(other.matrix_);
  // other.rows_ = 0;
  // other.cols_ = 0;
  // other.matrix_ = nullptr:
  std::swap(other.rows_, rows_);
  std::swap(other.cols_, cols_);
  std::swap(other.matrix_, matrix_);
}
S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (&other == this) return *this;
  S21Matrix boba(other);  // boba = other; this
  *this = move(boba);
  return *this;
}
S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
  if (&other == this) return *this;
  if (matrix_) {
    delete[] matrix_[0];
  }
  delete[] matrix_;
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
  // std::swap(other.rows_, rows_);
  // std::swap(other.cols_, cols_);
  // std::swap(other.matrix_, matrix_);
  return *this;
}

S21Matrix::~S21Matrix() {
  if (matrix_) {
    delete[] matrix_[0];
  }
  delete[] matrix_;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const noexcept {
  bool eq = true;
  if (rows_ != other.rows_ || cols_ != other.cols_) return false;
  for (size_t i = 0; i < rows_ && eq; ++i)
    for (size_t j = 0; j < cols_ && eq; ++j)
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7) return eq = false;
  return eq;
}
void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::logic_error("SumMatrix: Different matrix dimensions");
  std::transform(begin(), end(), other.begin(), begin(), std::plus<>());
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::logic_error("SubMatrix: Different matrix dimensions");
  std::transform(begin(), end(), other.begin(), begin(), std::minus<>());
}

void S21Matrix::MulNumber(const double num) noexcept {
  std::for_each(begin(), end(), [num](double& x) { x *= num; });
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_)
    throw std::logic_error("MulMatrix: Different matrix dimensions");
  *this = (*this * other);
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_)
    throw std::logic_error("Op+: Different matrix dimensions");
  // S21Matrix res(rows_, cols_);
  // std::transform(begin(), end(), other.begin(), res.begin(), std::plus<>());
  auto fun = [](double x, double y) { return x + y; };
  return Transform(other, fun);
}
S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_)
    throw std::logic_error("Op+=: Different matrix dimensions");
  SumMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::logic_error("Op-: Different matrix dimensions");
  S21Matrix res(rows_, cols_);
  std::transform(begin(), end(), other.begin(), res.begin(), std::minus<>());
  return res;
}
S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_)
    throw std::logic_error("Op-=: Different matrix dimensions");
  SubMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  if (cols_ != other.rows_)
    throw std::logic_error("Op*: Different matrix dimensions");
  S21Matrix tmp(rows_, cols_);
  for (size_t i = 0; i < rows_; ++i)
    for (size_t j = 0; j < other.cols_; ++j)
      for (size_t k = 0; k < cols_; ++k)
        tmp.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
  return tmp;
}
S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}
S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix res(rows_, cols_);
  std::transform(begin(), end(), res.begin(),
                 [num](double x) { return x * num; });
  return res;
}
S21Matrix& S21Matrix::operator*=(const double num) {
  std::for_each(matrix_[0], matrix_[0] + rows_ * cols_,
                [num](double& x) { x *= num; });
  return *this;
}
double& S21Matrix::operator()(const size_t i, const size_t j) {
  if (i > rows_ || j > cols_)
    throw std::out_of_range("Op(): Index is outside the matrix");
  return matrix_[i][j];
}
const double& S21Matrix::operator()(const size_t i, const size_t j) const {
  if (i > rows_ || j > cols_)
    throw std::out_of_range("Op const(): Index is outside the matrix");
  return matrix_[i][j];
}
bool S21Matrix::operator==(const S21Matrix& other) const noexcept {
  return EqMatrix(other);
}
bool S21Matrix::operator!=(const S21Matrix& other) const noexcept {
  return !EqMatrix(other);
}
S21Matrix S21Matrix::Transpose() {
  S21Matrix res(cols_, rows_);
  for (size_t i = 0; i < cols_; ++i)
    for (size_t j = 0; j < rows_; ++j) res.matrix_[i][j] = matrix_[j][i];
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix res(rows_, cols_);
  S21Matrix tmp(rows_ - 1, cols_ - 1);
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      size_t minor_row = 0;
      for (size_t k = 0; k < rows_; ++k) {
        bool num_was_written = false;
        size_t minor_col = 0;
        for (size_t l = 0; l < cols_; ++l) {
          if (k != i && l != j) {
            tmp.matrix_[minor_row][minor_col] = matrix_[k][l];
            ++minor_col;
            num_was_written = true;
          }
        }
        if (num_was_written) ++minor_row;
      }
      res.matrix_[i][j] = tmp.Determinant() * ((i + j) % 2 == 0 ? 1 : -1);
    }
  }
  return res;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_)
    throw std::logic_error("Determinant: The matrix is not square");
  double det = 1;
  int sign_det = 1;
  S21Matrix res(*this);
  size_t row = 0, col = 0;
  for (size_t i = 0; i < rows_ && det != 0; ++i) {
    if (std::fabs(res.matrix_[row][col]) < 1e-7) {
      sign_det = res.SwapLines(row, col);
    }
    det *= res.matrix_[row][col] * sign_det;
    res.GaussDiv(row, col);
    res.GaussSubLines(row, col);
    ++row;
    ++col;
  }
  return det * sign_det;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  if (fabs(det) < 1e-7)
    throw std::logic_error("InversMatrix: Matrix determinant is 0");
  S21Matrix res(rows_, cols_);
  res = CalcComplements();
  res = res.Transpose();
  res *= (1 / det);
  return res;
}

S21Matrix S21Matrix::Transform(const S21Matrix& other,
                               double (*foo)(double, double)) {
  S21Matrix res(rows_, cols_);
  auto start = matrix_[0];
  auto start_other = other.matrix_[0];
  auto start_res = res.matrix_[0];

  auto end = matrix_[0] + rows_ * cols_;
  while (start != end) {
    *start_res = foo(*start, *start_other);
    ++start_other;
    ++start;
    ++start_res;
  }
  return res;
}

S21Matrix S21Matrix::Transform(double (*foo)(double)) {
  S21Matrix res(rows_, cols_);
  double* start = matrix_[0];
  double* start_res = res.matrix_[0];
  double* edn = matrix_[0] + rows_ * cols_;

  for (; start != edn; ++start, ++start_res) {
    *start_res = foo(*start);
  }
  return res;
}

void S21Matrix::Resize(size_t newrows, size_t newcols) {
  if (newrows <= 0 || newcols <= 0)
    throw std::logic_error("Resize: Dimension can't be equal or less 0");
  if (newrows == rows_ && newcols == cols_) return;
  S21Matrix tmp(newrows, newcols);
  size_t rows_iter = rows_ < newrows ? rows_ : newrows;
  size_t cols_iter = cols_ < newcols ? cols_ : newcols;
  for (size_t i = 0; i < rows_iter; ++i)
    for (size_t j = 0; j < cols_iter; ++j) tmp(i, j) = matrix_[i][j];
  *this = std::move(tmp);
}

int S21Matrix::SwapLines(size_t row, size_t col) {
  int sign_det = 1;
  for (size_t i = 1; i < rows_ - row && std::fabs(matrix_[row][col]) < 1e-7;
       ++i) {
    if (std::fabs(matrix_[row + i][col]) > 1e-7) {
      double temp = matrix_[row][col];
      // std::cout << "temp - "<<temp<<std::endl;
      for (size_t j = col; j < cols_; ++j) {
        matrix_[row][j] = matrix_[row + i][j];
        matrix_[row + i][j] = temp;
        temp = matrix_[row][col + (j - col + 1)];
      }
      sign_det = -1;
    }
  }
  return sign_det;
}

void S21Matrix::GaussDiv(size_t row, size_t col) {
  double divr = matrix_[row][col];
  for (size_t i = col; i < cols_; ++i) {
    matrix_[row][i] /= divr;
  }
}

void S21Matrix::GaussSubLines(size_t row, size_t col) {
  for (size_t i = row + 1; i < rows_; ++i) {
    double mul = matrix_[i][col];
    for (size_t j = col; j < cols_; ++j) {
      matrix_[i][j] -= matrix_[row][j] * mul;
    }
  }
}

S21Matrix operator*(const double num, S21Matrix& other) { return other * num; }

std::ostream& operator<<(std::ostream& out, const S21Matrix& other) {
  for (size_t i = 0; i < other.rows_; ++i) {
    for (size_t j = 0; j < other.cols_; ++j) {
      out << other.matrix_[i][j] << "\t";
    }
    out << std::endl;
  }
  return out;
}

std::istream& operator>>(std::istream& in, S21Matrix& other) {
  for (auto& v : other) {
    in >> v;
  }
  return in;
}
};  // namespace S21
