/*
   Файл, содержащий в себе только вынесенные константы и константные функции для main.cpp
   Ни в коем случае не добавлять его никуда, кроме main.cpp!
*/

#pragma once
#include <string>
#include <stdexcept>

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
using Matrix = std::vector<std::vector<double>>;

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

double f_value(int regionNum, double r, double phi) {
   double ans = 0.0;
   switch (regionNum)
   {
      case 0: ans = 8.0; break;
      //case 1: ans = 1.0; break;
      default:
         throw std::runtime_error("Значения функции f для региона с номером " + std::to_string(regionNum) + " не найдено.");
         break;
   }
   return ans;
}
double f_value(int regionNum, Node node) {
   return f_value(regionNum, node.r, node.phi);
}

double lambda_value(int regionNum, double r, double phi) {
   double ans = 0.0;
   switch (regionNum)
   {
      case 0: ans = 1.0; break;
      //case 1: ans = 1.0; break;
      default:
         throw std::runtime_error("Значения функции lambda для региона с номером " + std::to_string(regionNum) + " не найдено.");
         break;
   }
   return ans;
}
double lambda_value(int regionNum, Node node) {
   return lambda_value(regionNum, node.r, node.phi);
}

double gamma_value(int regionNum, double r, double phi) {
   double ans = 0.0;
   switch (regionNum)
   {
      case 0: ans = 2.0; break;
      //case 1: ans = 1.0; break;
      default:
         throw std::runtime_error("Значения функции gamma для региона с номером " + std::to_string(regionNum) + " не найдено.");
         break;
   }
   return ans;
}
double gamma_value(int regionNum, Node node) {
   return gamma_value(regionNum, node.r, node.phi);
}

double s3_beta_value(int s3_funcNum, Node node) {
   double ans = 0.0;
   switch (s3_funcNum)
   {
      case 0: ans = 1.0; break;
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
      case 0: ans = 4.0; break;
      //case 1: ans = 1.0; break;
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
      case 0: ans = 0.0; break;
      case 1: ans = 0.0; break;
      default:
         throw std::runtime_error("Значения функции theta для s2-краевого с номером " + std::to_string(s2_funcNum) + " не найдено.");
         break;
   }
   return ans;
}

double s1_u_value(int s1_funcNum, Node node) {
   double ans = 0.0;
   switch (s1_funcNum)
   {
      case 0: ans = 4.0; break;
      //case 1: ans = 1.0; break;
      default:
         throw std::runtime_error("Значения функции u для s1-краевого с номером " + std::to_string(s1_funcNum) + " не найдено.");
         break;
   }
   return ans;
}