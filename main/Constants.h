/*
   Файл, содержащий в себе только вынесенные константы и константные функции для main.cpp
   Ни в коем случае не добавлять его никуда, кроме main.cpp!
*/

#pragma once
#include <string>
#include <stdexcept>
#include <array>

namespace GlobalPaths {
   // Пути файлов:
   const std::string filesPath = "./iofiles/";
   const std::string nodesPath = filesPath + "nodes.txt";
   const std::string rectanglesPath = filesPath + "rectangles.txt";
   const std::string s1_nodesPath = filesPath + "s1_nodes.txt";
   const std::string s2_edgesPath = filesPath + "s2_edges.txt";
   const std::string s3_edgesPath = filesPath + "s3_edges.txt";
}

#pragma region TYPEDEFINES
using Matrix = std::array<std::array<double, 4>, 4>;

/// <summary>
/// Структура прямоугольника, имеет 4 номера вершины: [a], [b], [c], [d], а также
/// номер области, в которой находится сам прямоугольник, [region]
/// </summary>
struct Rectangle {
   int a = 0;
   int b = 0;
   int c = 0;
   int d = 0;
   int regionNum = 0;

   std::string toString() const {
      std::string out = "( ";
      out += "a: " + std::to_string(a);
      out += ", b: " + std::to_string(b);
      out += ", c: " + std::to_string(c);
      out += ", d: " + std::to_string(d);
      out += ", region: " + std::to_string(regionNum);
      out += " )";

      return out;
   }
};

/// <summary>
/// Структура описания узла сетки. Содержит координаты этого узла [r] и [phi]
/// </summary>
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

double f_value(int regionNum, double r, double phi) {
   switch (regionNum)
   {
   case 0: {
      return phi*phi + r*r - 8*r*r - 2.0;
   }

   default:
      throw std::runtime_error("Значения функции f для региона с номером " + std::to_string(regionNum) + " не найдено.");
   }
}
double f_value(int regionNum, Node node) {
   return f_value(regionNum, node.r, node.phi);
}

double lambda_value(int regionNum, double r, double phi) {
   switch (regionNum)
   {
   case 0: {
      return r*r;
   }

   default:
      throw std::runtime_error("Значения функции lambda для региона с номером " + std::to_string(regionNum) + " не найдено.");
   }

}
double lambda_value(int regionNum, Node node) {
   return lambda_value(regionNum, node.r, node.phi);
}

double gamma_value(int regionNum, double r, double phi) {
   switch (regionNum)
   {
   case 0: {
      return 1;
   }

   default:
      throw std::runtime_error("Значения функции gamma для региона с номером " + std::to_string(regionNum) + " не найдено.");
   }
}
double gamma_value(int regionNum, Node node) {
   return gamma_value(regionNum, node.r, node.phi);
}

double s1_u_value(int s1_funcNum, double r, double phi) {
   switch (s1_funcNum) {
   case 0: {
      return phi*phi + r*r;
   }

   default:
      throw std::runtime_error("Значения функции u для s1-краевого с номером " + std::to_string(s1_funcNum) + " не найдено.");
   }
}
double s1_u_value(int s1_funcNum, Node node) {
   return s1_u_value(s1_funcNum, node.r, node.phi);
}

double s3_beta_value(int s3_funcNum, Node node) {
   double ans = 0.0;
   switch (s3_funcNum)
   {
      //case 0: ans = 1; break;
      //case 1: ans = 1.0; break;
   default:
      throw std::runtime_error("Значения функции beta для s3-краевого с номером " + std::to_string(s3_funcNum) + " не найдено.");
      break;
   }
   return ans;
}

double s3_u_value(int s3_funcNum, Node node) {
   double ans = 0.0;
   switch (s3_funcNum)
   {
      //case 0: ans = node.r - 1; break;
      //case 1: ans = node.r; break;
   default:
      throw std::runtime_error("Значения функции U_beta для s3-краевого с номером " + std::to_string(s3_funcNum) + " не найдено.");
      break;
   }
   return ans;
}

double s2_theta_value(int s2_funcNum, Node node) {
   double ans = 0.0;
   switch (s2_funcNum)
   {
      //case 0: ans = 0; break;
      //case 1: ans = 0.0; break;
   default:
      throw std::runtime_error("Значения функции theta для s2-краевого с номером " + std::to_string(s2_funcNum) + " не найдено.");
      break;
   }
   return ans;
}
