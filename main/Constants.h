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

double f_value(int regionNum, double r, double phi) {
   double ans = 0.0;
   switch (regionNum)
   {
      //case 0: ans = 0.0; break;
      //case 1: ans = 1.0; break;
      default:
         throw std::runtime_error("Значения функции f для региона с номером " + std::to_string(regionNum) + " не найдено.");
         break;
   }
   return ans;
}

double lambda_value(int regionNum, double r, double phi) {
   double ans = 0.0;
   switch (regionNum)
   {
      //case 0: ans = 0.0; break;
      //case 1: ans = 1.0; break;
      default:
         throw std::runtime_error("Значения функции lambda для региона с номером " + std::to_string(regionNum) + " не найдено.");
         break;
   }
   return ans;
}

double gamma_value(int regionNum, double r, double phi) {
   double ans = 0.0;
   switch (regionNum)
   {
      //case 0: ans = 0.0; break;
      //case 1: ans = 1.0; break;
      default:
         throw std::runtime_error("Значения функции gamma для региона с номером " + std::to_string(regionNum) + " не найдено.");
         break;
   }
   return ans;
}