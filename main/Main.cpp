#include <iostream>
#include <fstream>
#include <vector>

#include "ThreeSteppers/Headers/IterSolvers.h"

// Файл, содержащий в себе пути до файлов, функции f, lambda и gamma
#include "Constants.h"

using namespace std;

#pragma region GLOBAL_OBJECTS
// Глобальная разреженная матрица системы
SparseMatrix global_mat;
// Глобальный вектор системы
vector<double> global_b;
// Локальная матрица, будь то M, G, и прочие матрицы одного размера
Matrix local_mat;
// Локальный вектор, в пару локальным матрицам
vector<double> local_b;
// Массив прямоугольников
vector<Rectangle> rectangles;
// Массив узлов
vector<Node> nodes;
// Массив сопоставления узлов и первых краевых
vector<S1_node> s1_nodes;
// Массив сопоставления рёбер и вторых краевых
vector<S23_edge> s2_edges;
// Массив сопоставления рёбер и третьих краевых
vector<S23_edge> s3_edges;
#pragma endregion GLOBAL_OBJECTS

void readDataFromFiles() {
   // Считывание данных для структуры узлов nodes
   auto nodesFile = ifstream(GlobalPaths::nodesPath);
   if (!nodesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::nodesPath);
   int size;
   nodesFile >> size;
   nodes.resize(size);
   for (auto& node : nodes)
   {
      nodesFile >> node.r >> node.phi;
   }
   nodesFile.close();

   // Считывание данных для структуры прямоугольников rectangles
   auto rectanglesFile = ifstream(GlobalPaths::rectanglesPath);
   if (!rectanglesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::rectanglesPath);
   rectanglesFile >> size;
   rectangles.resize(size);
   for (auto& rect : rectangles)
   {
      rectanglesFile >> rect.a >> rect.b >> rect.c >> rect.d >> rect.regionNum;
   }
   rectanglesFile.close();

   // Считывание данных для первых краевых условий s1_nodes
   auto s1_nodesFile = ifstream(GlobalPaths::s1_nodesPath);
   if (!s1_nodesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::s1_nodesPath);
   s1_nodesFile >> size;
   s1_nodes.resize(size);
   for (auto& s1 : s1_nodes)
   {
      rectanglesFile >> s1.node >> s1.funcNum;
   }
   s1_nodesFile.close();

   // Считывание данных для вторых краевых условий s2_edges
   auto s2_edgesFile = ifstream(GlobalPaths::s2_edgesPath);
   if (!s2_edgesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::s2_edgesPath);
   s2_edgesFile >> size;
   s2_edges.resize(size);
   for (auto& s2 : s2_edges)
   {
      s2_edgesFile >> s2.node1 >> s2.node2 >> s2.funcNum;
   }
   s2_edgesFile.close();

   // Считывание данных для третьих краевых условий s3_edges
   auto s3_edgesFile = ifstream(GlobalPaths::s3_edgesPath);
   if (!s3_edgesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::s3_edgesPath);
   s3_edgesFile >> size;
   s3_edges.resize(size);
   for (auto& s3 : s3_edges)
   {
      s3_edgesFile >> s3.node1 >> s3.node2 >> s3.funcNum;
   }
   s3_edgesFile.close();
}

void generatePortrait() {
   global_mat.di.resize(nodes.size());
   global_mat.ig.resize(nodes.size() + 1);

   for (auto& rect : rectangles)
   {
      int elems[4] = { rect.a, rect.b, rect.c, rect.d };
      for (int i = 0; i < 4; i++)
      {
         for (int k = 0; k < 4; k++)
         {
            // Если элемент на диагонали или в верхнем треугольнике, то скипаем
            if (i == k || elems[k] > elems[i])
            {
               continue;
            }

            bool isExist = false;
            // Пробегаем по всей строке для проверки, существует ли такой элемент
            for (int it = global_mat.ig[elems[i]]; it < global_mat.ig[elems[i] + 1]; it++)
            {
               if (global_mat.jg[it] == elems[k]) isExist = true;
            }
            if (!isExist)
            {
               // Ищем, куда вставить элемент портрета
               int it = global_mat.ig[elems[i]];
               while (it < global_mat.ig[elems[i] + 1] && global_mat.jg[it] < elems[k]) it++;
               // Для вставки нужно взять итератор массива от начала, так что...
               global_mat.jg.insert(global_mat.jg.begin() + it, elems[k]);
               // Добавляем всем элементам ig с позиции elems[i]+1 один элемент
               for (int j = elems[i] + 1; j < global_mat.ig.size(); j++)
                  global_mat.ig[j]++;
            }
         }
      }
   }
   global_mat.ggl.resize(global_mat.jg.size());
   global_mat.ggu.resize(global_mat.jg.size());
}

void addLocalG(const Rectangle& rect) {
   double rp = nodes[rect.a].r;
   double hr = nodes[rect.b].r - nodes[rect.a].r;
   double phi_s = nodes[rect.a].phi;
   double h_phi = nodes[rect.c].phi - nodes[rect.a].phi;
   double lambda = (lambda_value(rect.regionNum, nodes[rect.a]) +
      lambda_value(rect.regionNum, nodes[rect.b]) +
      lambda_value(rect.regionNum, nodes[rect.c]) +
      lambda_value(rect.regionNum, nodes[rect.d])) / 4;

   double integrals[9];
   integrals[0] = rp / hr + 0.5;
   integrals[1] = -integrals[0];
   integrals[2] = (((rp * rp) / (hr * hr)) + (2 * rp / hr) + 1) * std::log((rp + hr) / rp) - (rp / hr) - 1.5;
   integrals[3] = -(rp / hr) * std::log((rp + hr) / rp) + 1;
   integrals[4] = -(rp / hr + (rp * rp) / (hr * hr)) * std::log((rp + hr) / rp) + (rp / hr) + 0.5;
   integrals[5] = h_phi / 3;
   integrals[6] = h_phi / 6;
   integrals[7] = 1 / h_phi;
   integrals[8] = -integrals[7];


   local_mat[0][0] = lambda * (integrals[0] * integrals[5] + integrals[2] * integrals[7]);
   local_mat[0][1] = lambda * (integrals[1] * integrals[5] + integrals[4] * integrals[7]);
   local_mat[0][2] = lambda * (integrals[0] * integrals[6] + integrals[2] * integrals[8]);
   local_mat[0][3] = lambda * (integrals[1] * integrals[6] + integrals[4] * integrals[8]);
   local_mat[1][0] = local_mat[0][1];
   local_mat[1][1] = lambda * (integrals[0] * integrals[5] + integrals[3] * integrals[7]);
   local_mat[1][2] = local_mat[0][3];
   local_mat[1][3] = lambda * (integrals[0] * integrals[6] + integrals[3] * integrals[8]);
   local_mat[2][0] = local_mat[0][2];
   local_mat[2][1] = local_mat[1][2];
   local_mat[2][2] = lambda * (integrals[0] * integrals[5] + integrals[2] * integrals[7]);
   local_mat[2][3] = lambda * (integrals[1] * integrals[5] + integrals[4] * integrals[7]);
   local_mat[3][0] = local_mat[0][3];
   local_mat[3][1] = local_mat[1][3];
   local_mat[3][2] = local_mat[2][3];
   local_mat[3][3] = lambda * (integrals[0] * integrals[5] + integrals[3] * integrals[7]);
}

void addLocalM(const Rectangle& rect) {
   double rp = nodes[rect.a].r;
   double hr = nodes[rect.b].r - nodes[rect.a].r;
   double phi_s = nodes[rect.a].phi;
   double h_phi = nodes[rect.c].phi - nodes[rect.a].phi;
   double gamma[4] = {
      gamma_value(rect.regionNum, nodes[rect.a]),
      gamma_value(rect.regionNum, nodes[rect.b]),
      gamma_value(rect.regionNum, nodes[rect.c]),
      gamma_value(rect.regionNum, nodes[rect.d])
   };
   double tmp[4][4];

   tmp[0][0] = hr * h_phi * (gamma[0] * (rp / 4 + hr / 20) / 4
      + gamma[1] * (rp / 12 + hr / 30) / 4
      + gamma[2] * (rp / 4 + hr / 20) / 12
      + gamma[3] * (rp / 12 + hr / 30)  / 12
      );
   tmp[0][1] = hr * h_phi * (gamma[0] * (rp / 12 + hr / 30) / 4
      + gamma[1] * (rp / 12 + hr / 20) / 4
      + gamma[2] * (rp / 12 + hr / 30) / 12
      + gamma[3] * (rp / 12 + hr / 20)  / 12
      );
   tmp[0][2] = hr * h_phi * (gamma[0] * (rp / 4 + hr / 20) / 12
      + gamma[1] * (rp / 12 + hr / 30) / 12
      + gamma[2] * (rp / 4 + hr / 20) / 12
      + gamma[3] * (rp / 12 + hr / 30)  / 12
      );
   tmp[0][3] = hr * h_phi * (gamma[0] * (rp / 12 + hr / 30) / 12
      + gamma[1] * (rp / 12 + hr / 20) / 12
      + gamma[2] * (rp / 12 + hr / 30) / 12
      + gamma[3] * (rp / 12 + hr / 20)  / 12
      );
   tmp[1][0] = tmp[0][1];
   tmp[1][1] = hr * h_phi * (gamma[0] * (rp / 12 + hr / 20) / 4
      + gamma[1] * (rp / 4 + hr / 5) / 4
      + gamma[2] * (rp / 12 + hr / 20) / 12
      + gamma[3] * (rp / 4 + hr / 5) / 12
      );
   tmp[1][2] = hr * h_phi * (gamma[0] * (rp / 12 + hr / 30) / 12
      + gamma[1] * (rp / 12 + hr / 20) / 12
      + gamma[2] * (rp / 12 + hr / 30) / 12
      + gamma[3] * (rp / 12 + hr / 20) / 12
      );
   tmp[1][3] = hr * h_phi * (gamma[0] * (rp / 12 + hr / 20) / 12
      + gamma[1] * (rp / 4 + hr / 5) / 12
      + gamma[2] * (rp / 12 + hr / 20) / 12
      + gamma[3] * (rp / 4 + hr / 5) / 12
      );
   tmp[2][0] = tmp[0][2];
   tmp[2][1] = tmp[1][2];
   tmp[2][2] = hr * h_phi * (gamma[0] * (rp / 4 + hr / 20) / 12
      + gamma[1] * (rp / 12 + hr / 30) / 12
      + gamma[2] * (rp / 4 + hr / 20) / 4
      + gamma[3] * (rp / 12 + hr / 30) / 4
      );
   tmp[2][3] = hr * h_phi * (gamma[0] * (rp / 12 + hr / 30) / 12
      + gamma[1] * (rp / 12 + hr / 20) / 12
      + gamma[2] * (rp / 12 + hr / 30) / 4
      + gamma[3] * (rp / 12 + hr / 20) / 4
      );
   tmp[3][0] = tmp[0][3];
   tmp[3][1] = tmp[1][3];
   tmp[3][2] = tmp[2][3];
   tmp[3][3] = hr * h_phi * (gamma[0] * (rp / 12 + hr / 20) / 12
      + gamma[1] * (rp / 4 + hr / 5) / 12
      + gamma[2] * (rp / 12 + hr / 20) / 4
      + gamma[3] * (rp / 4 + hr / 5) / 4
      );
   for (int i = 0; i < 4; i++)
   {
      for (int k = 0; k < 4; k++)
      {
         local_mat[i][k] += tmp[i][k];
      }
   }
}

void addLocalsToGlobal(const Rectangle& rect) {
   addLocalG(rect);
   addLocalM(rect);
   //addLocalB(rect);
   //addAllToGlobal(rect);
}

void main() {
   readDataFromFiles();
   generatePortrait();

   local_mat.resize(4);
   for (auto& vec : local_mat)
      vec.resize(4);
   local_b.resize(4);

   for (const auto& rect : rectangles)
   {
      addLocalsToGlobal(rect);
   }
   //include_s3();
   //include_s2();
   //include_s1();

   // дохуя важные расчёты исходной функции в узлах
   // дохуя важные расчёты невязки полученной и исходной
   // вывод дохуя важной информации

   return;
}