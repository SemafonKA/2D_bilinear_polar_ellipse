#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <string>
#include <format>

std::vector<double> ReadVecFromFile(size_t size, const std::string& path);

/// <summary>
/// Класс объектов матриц, хранящихся в разреженном строчно-столбцовом виде
/// <para> Точность хранения элементов - double </para>
/// </summary>
class SparseMatrix {
// Переменные матрицы
public:
   /// <summary>
   /// Массив индексов строк/столбцов, вида 0, 0, 0 + k2, ..., 0+k2+...+kn, где ki - число элементов в i cтроке/столбце
   /// <para> Помимо этого первый элемент i строки можно найти как ggl[ig[i]] </para>
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ig: { 0, 0, 1, 2 } </para>
   /// </summary>
   std::vector<uint32_t> ig;

   /// <summary>
   /// Массив индексов столбцов/строк элементов (ставит индекс в соответствие элементу)
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> jg: { 0, 1 } </para>
   /// </summary>
   std::vector<uint16_t> jg;

   /// <summary>
   /// Массив элементов нижнего треугольника матрицы
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ggl: { 3, 2 } </para>
   /// </summary>
   std::vector<double> ggl;

   /// <summary>
   /// Массив элементов верхнего треугольника матрицы
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> ggu: { 2, 1 } </para>
   /// </summary>
   std::vector<double> ggu;

   /// <summary>
   /// Массив элементов диагонали матрицы
   /// <para> Пример массива для матрицы 3х3: </para>
   ///
   /// <para> Матрица: </para>
   /// <para> | 1 2 0 | </para>
   /// <para> | 3 8 1 | </para>
   /// <para> | 0 2 4 | </para>
   /// <para> di: { 1, 8, 4 } </para>
   /// </summary>
   std::vector<double> di;

// Методы матрицы
public:
   size_t Size() const;

   /// <summary>
   /// Умножение матрицы на вектор
   /// </summary>
   std::vector<double> MultToVec(const std::vector<double>& right) const;
   std::vector<double>& MultToVec(const std::vector<double>& right, std::vector<double>& result) const;

   std::vector<double> operator*(const std::vector<double>& right) const;

   /// <summary>
   /// Умножение транспонированной матрицы на вектор
   /// </summary>
   std::vector<double> TranspMultToVec(const std::vector<double>& right) const;
   std::vector<double>& TranspMultToVec(const std::vector<double>& right, std::vector<double>& result) const;

   SparseMatrix& operator= (SparseMatrix&& right) noexcept;

   double val(uint16_t row, uint16_t column) {
      // if element lay on diagonal
      if (row == column) return di[row];

      auto& v = row > column ? ggl : ggu;
      if (row < column) std::swap(row, column);

      // find element by his pos
      auto i = ig[row];
      while (i < ig[row + 1ll] && jg[i] < column) i++;

      // if element exists, return him
      if (i < ig[row + 1ll] && jg[i] == column) return v[i];

      // else return zero
      return 0.0;
   }

   std::string toStringAsDense() noexcept {
      std::string out = "[ ";
      auto size = Size();

      for (auto i = 0; i < size; i++)
      {
         if (i != 0) out += "  ";
         out += "[ ";
         for (auto j = 0; j < size; j++)
         {
            out += std::format("{: .5f}", val(i, j)); //std::to_string(val(i, j));
            if (j + 1 < size) out += ", ";
         }
         out += " ]";
         if (i + 1 < size) out += "\n";
      }
      out += " ]";

      return out;
   }

// Конструкторы матрицы
public:
   SparseMatrix();

   SparseMatrix(SparseMatrix& right);

   // Конструктор перемещения (нужен для метода ReadFromFiles)
   SparseMatrix(SparseMatrix&& right) noexcept;

// Статические методы матрицы
public:
   static SparseMatrix ReadFromFiles(uint16_t matrixSize, const std::string& igP, const std::string& jgP, const std::string& gglP, const std::string& gguP, const std::string& diP);
};