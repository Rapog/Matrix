#include <gtest/gtest.h>

#include "S21_matrix_oop.h"

int main() {
  S21::S21Matrix mat(2, 2);
  mat(0, 0) = 3;
  mat(0, 1) = 2;
  mat(1, 0) = -6;
  mat(1, 1) = 0;
  S21::S21Matrix mat1(2, 2);
  mat1(0, 0) = -7;
  mat1(0, 1) = 0;
  mat1(1, 0) = -3;
  mat1(1, 1) = 2;
  S21::S21Matrix res(2, 2);
  // std::cin >> res;
  res = mat + mat1;
  // mat.SubMatrix(mat1);
  std::cout << mat;
  std::cout << std::endl;
  std::cout << mat1;
  std::cout << std::endl;
  std::cout << res;
  // // mat.SetNum(2.3, 0, 0);
  // // mat1.SetNum(2.3, 0, 0);
  // // mat.SetNum(2.3, 1, 1);
  // // mat1.SetNum(4.3, 1, 1);
  // int count = 0;
  // mat.Resize(2,2);
  // for(size_t i = 0; i < mat.GetRows(); ++i)
  //   for(size_t j = 0; j < mat.GetCols(); ++j)
  //     {
  //       mat(i,j) = count;
  //       // mat1.SetNum(i+j, i, j);
  //       ++count;
  //     }
  // mat(0,0) = 234;
  // std::cout << std::endl;
  // std::cout << mat;
  // mat.Resize(1,1);
  // std::cout << std::endl;
  // std::cout << mat;
  // mat(1,0) = 0;
  // mat(2,0) = 0;
  // mat(1,1) = 0;
  // mat.SumMatrix(mat1);
  // mat.SubMatrix(mat1);
  // std::cin >> mat;
  // std::stringstream input_date;
  // input_date << "1";
  // input_date >> mat;
  // mat.PrintMatrix();
  // std::cout << std::endl<<std::endl<<std::endl;
  // std::cout << mat;
  // mat1.PrintMatrix();
  // mat.MulNumber(2);
  // res = mat * mat1;
  // res = mat.Transpose();
  // double rdet = mat.Determinant();
  // std::cout << rdet << std::endl;
  // std::cout << res
  // mat *= 2;
  // mat1 = mat;
  // mat.PrintMatrix();
  // std::iterator_traits<double*>::value_type a;
  // if((typeid(a) == typeid(mat1.GetNum(1, 1)))) std:: cout << "FUCK TEMPLATE
  // DEDUCTION"; auto foo = [](double& x){ x = x*x; };
  // S21::S21Matrix::ForEach(mat.begin(), mat.end(), foo);
  // mat.PrintMatrix();
  // mat *= mat1;
  // res = mat * mat1;
  // mat(1,1) = 3243.54;
  // std::cout << mat(1,1) << std::endl;
  // std::cout << mat1.GetNum(1,1) << std::endl;
  // std::cout << res.GetNum(1,1) << std::endl;
  // mat.PrintMatrix();
  // mat1.PrintMatrix();
  // res.PrintMatrix();
  // foo(a, std::move(mat));

  // std::cout << mat.EqMatrix(mat1) << std::endl;
  // std::for_each
  return 0;
}

// 0 1
// 1 2