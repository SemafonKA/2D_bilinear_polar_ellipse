#include <iostream>
#include <fstream>
#include <vector>

#include "ThreeSteppers/Headers/IterSolvers.h"

// Файл, содержащий в себе пути до файлов, функции f, lambda и gamma
#include "Constants.h"

using namespace std;

#pragma region TYPEDEFINES
using Matrix = vector<vector<double>>;

struct Rectangle {
   int a = 0;
   int b = 0;
   int c = 0;
   int d = 0;
   int regionNum = 0;
};

struct Node {
   double r = 0.0;
   double phi = 0.0;
};

struct S1_node {
   int node = 0;
   int funcNum = 0;
};

struct S23_edge {
   int node1 = 0;
   int node2 = 0;
   int funcNum = 0;
};
#pragma endregion TYPEDEFINES

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

void main() {
   readDataFromFiles();
}