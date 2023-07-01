#include <fstream>
#include <iostream>
#include <vector>

#include "ThreeSteppers/Headers/IterSolvers.h"
#include "gaussian_quadrature.h"

// Файл, содержащий в себе пути до файлов, функции f, lambda и gamma
#include "Constants.h"

using namespace std;

#pragma region GLOBAL_OBJECTS
// Глобальная разреженная матрица системы
SparseMatrix global_mat;
// Глобальный вектор системы
vector<double> global_b;
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

double func_R(int ind, double rp, double hr, double r) {
   ind = ind % 2;
   if (ind == 0) {
      return (rp + hr - r) / hr;
   }
   return (r - rp) / hr;
}

double func_R_dif(int ind, double rp, double hr, double r) {
   ind = ind % 2; 
   if (ind == 0) {
      return -1.0 / hr;
   }
   return 1.0 / hr;
}

double func_Phi(int ind, double phi_s, double h_phi, double phi) {
   ind = (ind / 2) % 2;
   if (ind == 0) {
      return (phi_s + h_phi - phi) / h_phi;
   }
   return (phi - phi_s) / h_phi;
}

double func_Phi_dif(int ind, double phi_s, double h_phi, double phi) {
   ind = (ind / 2) % 2;
   if (ind == 0) {
      return -1.0 / h_phi;
   }
   return 1.0 / h_phi;
}

void readDataFromFiles() {
   // Считывание данных для структуры узлов nodes
   auto nodesFile = ifstream(GlobalPaths::nodesPath);
   if (!nodesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::nodesPath);
   int size;
   nodesFile >> size;
   nodes.resize(size);
   for (auto& node : nodes) {
      nodesFile >> node.r >> node.phi;
   }
   nodesFile.close();

   // Считывание данных для структуры прямоугольников rectangles
   auto rectanglesFile = ifstream(GlobalPaths::rectanglesPath);
   if (!rectanglesFile.is_open())
      throw runtime_error("Не удалось открыть файл " +
         GlobalPaths::rectanglesPath);
   rectanglesFile >> size;
   rectangles.resize(size);
   for (auto& rect : rectangles) {
      rectanglesFile >> rect.a >> rect.b >> rect.c >> rect.d >> rect.regionNum;
   }
   rectanglesFile.close();

   // Считывание данных для первых краевых условий s1_nodes
   auto s1_nodesFile = ifstream(GlobalPaths::s1_nodesPath);
   if (!s1_nodesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::s1_nodesPath);
   s1_nodesFile >> size;
   s1_nodes.resize(size);
   for (auto& s1 : s1_nodes) {
      s1_nodesFile >> s1.node >> s1.funcNum;
   }
   s1_nodesFile.close();

   // Считывание данных для вторых краевых условий s2_edges
   auto s2_edgesFile = ifstream(GlobalPaths::s2_edgesPath);
   if (!s2_edgesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::s2_edgesPath);
   s2_edgesFile >> size;
   s2_edges.resize(size);
   for (auto& s2 : s2_edges) {
      s2_edgesFile >> s2.node1 >> s2.node2 >> s2.funcNum;
   }
   s2_edgesFile.close();

   // Считывание данных для третьих краевых условий s3_edges
   auto s3_edgesFile = ifstream(GlobalPaths::s3_edgesPath);
   if (!s3_edgesFile.is_open())
      throw runtime_error("Не удалось открыть файл " + GlobalPaths::s3_edgesPath);
   s3_edgesFile >> size;
   s3_edges.resize(size);
   for (auto& s3 : s3_edges) {
      s3_edgesFile >> s3.node1 >> s3.node2 >> s3.funcNum;
   }
   s3_edgesFile.close();
}

void generatePortrait() {
   global_mat.di.resize(nodes.size());
   global_mat.ig.resize(nodes.size() + 1);

   for (auto& rect : rectangles) {
      const int elems[4] = { rect.a, rect.b, rect.c, rect.d };
      for (int i = 0; i < 4; i++) {
         for (int k = 0; k < i; k++) {
            // Если элемент в верхнем прямоугольнике, то скипаем
            if (elems[k] > elems[i]) continue;

            bool isExist = false;
            // Пробегаем по всей строке для проверки, существует ли такой элемент
            for (auto it = global_mat.ig[elems[i]];
               it < global_mat.ig[elems[i] + 1ll]; it++) {
               if (global_mat.jg[it] == elems[k]) {
                  isExist = true;
                  break;
               }
            }
            if (!isExist) {
               // Ищем, куда вставить элемент портрета
               auto it = global_mat.ig[elems[i]];
               while (it < global_mat.ig[elems[i] + 1ll] &&
                  global_mat.jg[it] < elems[k])
                  it++;

               // Для вставки нужно взять итератор массива от начала, так что...
               global_mat.jg.insert(global_mat.jg.begin() + it, elems[k]);

               // Добавляем всем элементам ig с позиции elems[i]+1 один элемент
               for (auto j = elems[i] + 1; j < global_mat.ig.size(); j++)
                  global_mat.ig[j]++;
            }
         }
      }
   }
   global_mat.ggl.resize(global_mat.jg.size());
   global_mat.ggu.resize(global_mat.jg.size());
}

Matrix getLocalG(const Rectangle& rect) {
   Matrix g = {};

   double rp = nodes[rect.a].r;
   double hr = abs(nodes[rect.b].r - nodes[rect.a].r);
   double phi_s = nodes[rect.a].phi;
   double h_phi = abs(nodes[rect.c].phi - nodes[rect.a].phi);

   int i, j;

   // [i] & [j] variables are linked to [solverFunc] function
   auto solverFunc = [&](double r, double phi) {
      double ans = 0.0;
      ans += func_Phi(i, phi_s, h_phi, phi) * func_R_dif(i, rp, hr, r) *
             func_Phi(j, phi_s, h_phi, phi) * func_R_dif(j, rp, hr, r);
      ans += (1.0 / (r * r)) *
         func_R(i, rp, hr, r) * func_Phi_dif(i, phi_s, h_phi, phi) * 
         func_R(j, rp, hr, r) * func_Phi_dif(j, phi_s, h_phi, phi);
      ans *= lambda_value(rect.regionNum, r, phi) * r;
      return ans;
   };

   auto solver = Gaussian_4p::TwoDimentionalSolver::withStep(rp, hr, phi_s,
      h_phi, solverFunc);

   for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
         g[i][j] = solver.compute(); // here i and j are linked to solver by lambda
      }
   }

   // debug output
#ifndef NDEBUG
   cout << "Local_G:" << endl;
   for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
         cout << std::format(" {: .5f}", g[i][j]);
      }
      cout << endl;
   }
   cout << endl;

   // All tmp must be almost equal zero
   cout << "All this values must be zero:\n";
   for (i = 0; i < 4; i++) {
      double tmp = 0;
      for (j = 0; j < 4; j++) {
         tmp += g[i][j];
      }
      cout << " " << (tmp > 5e-14 ? tmp : 0);
   }
   cout << endl << endl;
#endif

   return g;
}

Matrix getLocalM(const Rectangle& rect, bool getWithoutGamma = false) {
   Matrix m = {};
   std::function<double(int, double, double)> maybeGamma;
   if (getWithoutGamma == true) {
      maybeGamma = [](int reg, double r, double phi) {
         return 1.0;
      };
   }
   else {
      maybeGamma = [](int reg, double r, double phi) {
         return gamma_value(reg, r, phi);
      };
   }

   double rp = nodes[rect.a].r;
   double hr = abs(nodes[rect.b].r - nodes[rect.a].r);
   double phi_s = nodes[rect.a].phi;
   double h_phi = abs(nodes[rect.c].phi - nodes[rect.a].phi);

   int i, j;

   // [i] & [j] variables are linked to [solverFunc] function
   auto solverFunc = [&](double r, double phi) {
      double res = 1.0;
      res *= func_R(i, rp, hr, r) * func_Phi(i, phi_s, h_phi, phi);
      res *= func_R(j, rp, hr, r) * func_Phi(j, phi_s, h_phi, phi);
      res *= r * maybeGamma(rect.regionNum, r, phi);
      return res;
   };

   auto solver = Gaussian_4p::TwoDimentionalSolver::withStep(rp, hr, phi_s,
      h_phi, solverFunc);

   for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
         m[i][j] = solver.compute(); // [i] & [j] variables are linked to [solverFunc] function
      }
   }

#ifndef NDEBUG

   cout << "Local_M" << endl;

   for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
         cout << std::format(" {: .5f}", m[i][j]);
      }
      cout << endl;
   }
   cout << endl;

#endif

   return m;
}

std::vector<double> getLocalB(const Rectangle& rect) {
   std::vector<double> b(4);

   double rp = nodes[rect.a].r;
   double hr = abs(nodes[rect.b].r - nodes[rect.a].r);
   double phi_s = nodes[rect.a].phi;
   double h_phi = abs(nodes[rect.c].phi - nodes[rect.a].phi);

   int i;

   // Данный код по хорошему должен был работать лучше, поскольку использует
   // более точный интеграл, но по факту - нифига подобного
   // 
   //auto solverFunc = [&](double r, double phi) {
   //   double res = 1.0;
   //   res *= func_R(i, rp, hr, r) * func_Phi(i, phi_s, h_phi, phi);
   //   res *= r * f_value(rect.regionNum, r, phi);
   //   return res;
   //};

   //auto solver = Gaussian_4p::TwoDimentionalSolver::withStep(rp, hr, phi_s,
   //   h_phi, solverFunc);

   //for (i = 0; i < 4; i++) {
   //   b[i] = solver.compute();
   //}

   auto M = getLocalM(rect, true);
   for (int i = 0; i < 4; i++) {
      b[i] = 0;
      b[i] += M[i][0] * f_value(rect.regionNum, nodes[rect.a]);
      b[i] += M[i][1] * f_value(rect.regionNum, nodes[rect.b]);
      b[i] += M[i][2] * f_value(rect.regionNum, nodes[rect.c]);
      b[i] += M[i][3] * f_value(rect.regionNum, nodes[rect.d]);
   }

#ifndef NDEBUG
   // debug output
   cout << "Local b:\n";
   for (i = 0; i < 4; i++) {
      cout << std::format(" {: .5f}", b[i]);
   }
   cout << endl << endl;

#endif

   return b;
}

void addLocalMatrixToGlobal(const Rectangle& rect, SparseMatrix& globalMat, const Matrix& localMat) {
   const int elems[4] = { rect.a, rect.b, rect.c, rect.d };
   for (int i = 0; i < 4; i++) {
      // добавляем все внедиагональные элементы на строке elems[i]
      for (int k = 0; k < i; k++) {
         // Если элемент в верхнем прямоугольнике, то скипаем
         if (elems[k] > elems[i]) {
            continue;
         }

         auto id = globalMat.ig[elems[i]];
         for (id; id < globalMat.ig[elems[i] + 1ll] && globalMat.jg[id] != elems[k]; id++);

         globalMat.ggl[id] += localMat[i][k];
         globalMat.ggu[id] += localMat[k][i];
      }
      // добавляем диагональные элементы
      globalMat.di[elems[i]] += localMat[i][i];
   }
}

void addLocalbToGlobal(const Rectangle& rect, std::vector<double>& globalVec, const std::vector<double>& localVec) {
   const int elems[4] = { rect.a, rect.b, rect.c, rect.d };
   for (int i = 0; i < 4; i++) {
      globalVec[elems[i]] += localVec[i];
   }
}

void addLocalsToGlobal(const Rectangle& rect) {
#ifndef NDEBUG

   cout << "Current rect:\n";
   cout << rect.toString() << endl << endl;

#endif

   auto localG = getLocalG(rect);
   addLocalMatrixToGlobal(rect, global_mat, localG);

   auto localM = getLocalM(rect);
   addLocalMatrixToGlobal(rect, global_mat, localM);

   auto localVec = getLocalB(rect);
   addLocalbToGlobal(rect, global_b, localVec);

   // debug output
#ifndef NDEBUG

   cout << "global matrix at this step: " << endl << endl;
   cout << global_mat.toStringAsDense() << endl << endl;
   cout << "Global vector at this step: " << endl;
   cout << "[";
   for (auto i = 0; i < global_b.size(); i++) {
      cout << std::format(" {: .5f}", global_b[i]);
   }
   cout << " ]\n\n\n";

#endif
}

void include_s3() {
   for (const auto& edge : s3_edges) {
      double M[2][2] = {};
      double b[2] = {};
      double beta = (s3_beta_value(edge.funcNum, nodes[edge.node1]) +
         s3_beta_value(edge.funcNum, nodes[edge.node2])) /
         2;
      double u[2] = { s3_u_value(edge.funcNum, nodes[edge.node1]),
                     s3_u_value(edge.funcNum, nodes[edge.node2]) };
      double rp = nodes[edge.node1].r;
      double hr = nodes[edge.node2].r - nodes[edge.node1].r;
      double phi_s = nodes[edge.node1].phi;
      double h_phi = nodes[edge.node2].phi - nodes[edge.node1].phi;
      // Если краевое задано вдоль оси r
      if (hr > 1e-7) {
         M[0][0] = beta * ((hr * hr) / 12 + (rp * hr) / 3);
         M[0][1] = beta * ((hr * rp) / 6 + (hr * hr) / 12);
         M[1][0] = M[0][1];
         M[1][1] = beta * (rp * hr / 3 + hr * hr / 4);
         b[0] = beta * (u[0] * (hr * rp / 3 + hr * hr / 12) +
            u[1] * (hr * rp / 6 + hr * hr / 12));
         b[1] = beta * (u[0] * (hr * rp / 6 + hr * hr / 12) +
            u[1] * (hr * rp / 3 + hr * hr / 4));
      }
      // Если краевое задано вдоль оси phi
      else {
         M[0][0] = beta * rp * h_phi / 3;
         M[0][1] = beta * rp * h_phi / 6;
         M[1][0] = M[0][1];
         M[1][1] = beta * rp * h_phi / 3;
         b[0] = beta * rp * (u[0] * h_phi / 3 + u[1] * h_phi / 6);
         b[1] = beta * rp * (u[0] * h_phi / 6 + u[1] * h_phi / 3);
      }

      // добавляем полученный результат в глобальную матрицу
      int elems[2] = { edge.node1, edge.node2 };
      for (int i = 0; i < 2; i++) {
         // добавляем все внедиагональные элементы на строке elems[i]
         for (int k = 0; k < i; k++) {
            auto id = global_mat.ig[elems[i]];
            for (id; id < global_mat.ig[elems[i] + 1ll] &&
               global_mat.jg[id] != elems[k];
               id++)
               ;
            global_mat.ggl[id] += M[i][k];
            global_mat.ggu[id] += M[i][k];
         }
         // добавляем диагональные элементы и вектор b
         global_mat.di[elems[i]] += M[i][i];
         global_b[elems[i]] += b[i];
      }
   }
}

void include_s2() {
   for (const auto& edge : s2_edges) {
      double b[2] = {};
      double theta[2] = { s2_theta_value(edge.funcNum, nodes[edge.node1]),
                         s2_theta_value(edge.funcNum, nodes[edge.node2]) };
      double rp = nodes[edge.node1].r;
      double hr = nodes[edge.node2].r - nodes[edge.node1].r;
      double phi_s = nodes[edge.node1].phi;
      double h_phi = nodes[edge.node2].phi - nodes[edge.node1].phi;
      // Если краевое задано вдоль оси r
      if (hr > 1e-7) {
         b[0] = (theta[0] * (hr * rp / 3 + hr * hr / 12) +
            theta[1] * (hr * rp / 6 + hr * hr / 12));
         b[1] = (theta[0] * (hr * rp / 6 + hr * hr / 12) +
            theta[1] * (hr * rp / 3 + hr * hr / 4));
      }
      // Если краевое задано вдоль оси phi
      else {
         b[0] = rp * (theta[0] * h_phi / 3 + theta[1] * h_phi / 6);
         b[1] = rp * (theta[0] * h_phi / 6 + theta[1] * h_phi / 3);
      }

      // добавляем полученный результат в глобальную матрицу
      int elems[2] = { edge.node1, edge.node2 };
      for (int i = 0; i < 2; i++) {
         // добавляем вектор b
         global_b[elems[i]] += b[i];
      }
   }
}

void include_s1() {
   for (const auto& node : s1_nodes) {
      double u = s1_u_value(node.funcNum, nodes[node.node]);

      // ставим на диагональ значение 1
      global_mat.di[node.node] = 1;
      // ставим в соответствующую ячейку вектора b значение u
      global_b[node.node] = u;
      // зануляем строку в нижнем треугольнике
      for (auto j = global_mat.ig[node.node]; j < global_mat.ig[node.node + 1ll];
         j++) {
         global_mat.ggl[j] = 0;
      }
      // зануляем строку в верхнем треугольнике
      for (int i = node.node + 1; i < global_mat.Size(); i++) {
         for (auto j = global_mat.ig[i]; j < global_mat.ig[i + 1ll]; j++) {
            if (global_mat.jg[j] == node.node) {
               global_mat.ggu[j] = 0;
               break;
            }
         }
      }
   }
}

int main() {
   setlocale(LC_ALL, "ru-RU");
   readDataFromFiles();
   generatePortrait();
   global_b.resize(global_mat.Size());

   for (const auto& rect : rectangles) {
      addLocalsToGlobal(rect);
   }
   include_s3();
   include_s2();
   include_s1();

   // debug output
#ifndef NDEBUG

   cout << "After include edges:\n\n Matrix:\n\n";
   cout << global_mat.toStringAsDense() << endl << endl;
   cout << "Vector b:\n\n";
   cout << "[";
   for (auto i = 0; i < global_b.size(); i++) {
      cout << std::format(" {: .5f}", global_b[i]);
   }
   cout << " ]\n\n\n";

#endif

   vector<double> q;
   q.resize(global_mat.Size());
   //IterSolvers::MSG_Assimetric::Init_Default(q.size());
   IterSolvers::LOS::Init_LuPrecond(q.size(), global_mat);
   IterSolvers::minEps = 1e-20;
   double eps;
   //IterSolvers::MSG_Assimetric::Default(global_mat, global_b, q, eps);
   IterSolvers::LOS::LuPrecond(global_mat, global_b, q, eps);
   IterSolvers::Destruct();

   cout << "Полученное решение: " << endl;
   for (auto elem : q) {
      cout << format("{: .14f}", elem) << endl;
   }

   double nev = 0;
   cout << "Погрешность полученного решения: " << endl;
   for (auto i = 0; i < q.size(); i++) {
      double delta = s1_u_value(0, nodes.at(i)) - q.at(i);
      nev += std::pow(delta, 2);
      cout << format("{: .14f}", delta) << endl;
   }
   nev = std::sqrt(nev);
   cout << "Погрешность решения: " << nev << endl;

   return 0;
}