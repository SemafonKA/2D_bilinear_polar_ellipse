#include "../Headers/IterSolvers.h"
using namespace std;

namespace Vec {
   inline double Scalar(const vector<double>& l, const vector<double>& r) {
      if (l.size() != r.size()) throw runtime_error("Размеры векторов не совпадают");

      double res = 0.0;
      for (size_t i = 0; i < l.size(); i++)
      {
         res += l[i] * r[i];
      }

      return res;
   }

   // l or r may be similar vectors to ans
   inline void Mult(const vector<double>& l, const vector<double>& r, vector<double>& ans) {
      if (ans.size() != l.size() || ans.size() != r.size()) throw runtime_error("Ошибка: размеры векторов должны совпадать.");

      for (size_t i = 0; i < ans.size(); i++)
      {
         ans[i] = l[i] * r[i];
      }
   }
   inline vector<double> Mult(const vector<double>& l, const vector<double>& r) {
      if (r.size() != l.size()) throw runtime_error("Ошибка: размеры векторов должны совпадать.");
      vector<double> ans(l.size());

      for (size_t i = 0; i < ans.size(); i++)
      {
         ans[i] = l[i] * r[i];
      }
      return ans;
   }
}

namespace IterSolvers {
   double minEps = 1e-8;
   size_t maxIter = 2000;
   bool globalDebugOutput = false;

   std::vector<double>* _tmp1 = nullptr, * _tmp2 = nullptr,
      * _tmp3 = nullptr, * _tmp4 = nullptr, * _tmp5 = nullptr;
   LU* _lu_mat = nullptr;

   inline void VecInit(vector<double>*& vec, size_t size) {
      if (vec == nullptr)
      {
         vec = new vector<double>(size);
      }
      else if (vec->size() != size)
      {
         vec->resize(size);
      }
   }

   namespace LOS {
      size_t resetIter = 10;

      void Init_LuPrecond(size_t diSize, const SparseMatrix& A) {
         VecInit(_tmp1, diSize); // Массив для вектора r метода
         VecInit(_tmp2, diSize); // Массив для вектора z
         VecInit(_tmp3, diSize); // Массив для вектора p
         VecInit(_tmp4, diSize); // Массив для вектора Ar
         VecInit(_tmp5, diSize); // Массив для вектора tmp

         if (_lu_mat == nullptr)
         {
            _lu_mat = new LU(A);
         }
         else
         {
            _lu_mat->MakeLuFor(A);
         }
      }

      size_t LuPrecond(const SparseMatrix& A, const vector<double>& f, vector<double>& x, double& eps, bool debugOutput) {
         uint16_t size = x.size();

         LU& lu = *_lu_mat;

         vector<double>& tmp = *_tmp5;
         vector<double>& r = *_tmp1;               // r0 = L^-1 * (f - A * x)
         A.MultToVec(x, tmp);
         for (uint16_t i = 0; i < size; i++) tmp[i] = f[i] - tmp[i];
         lu.LSlauSolve(tmp, r);

         vector<double>& z = *_tmp2;               // z0 = U^-1 * r
         lu.USlauSolve(r, z);

         vector<double>& p = *_tmp3;               // p0 = L^-1 * A * z0
         A.MultToVec(z, tmp);
         lu.LSlauSolve(tmp, p);

         vector<double>& Ar = *_tmp4;              // Ar = L^-1 * A * U^-1 * r

         double ppScalar;
         double nev = Vec::Scalar(r, r);
         double ffScalar = Vec::Scalar(f, f);
         eps = nev / ffScalar;
         double a;                  // alpha
         double b;                  // beta
         size_t iter;

         for (iter = 1; iter <= maxIter && eps > minEps; iter++)
         {
            ppScalar = Vec::Scalar(p, p);          // (p_k-1, p_k-1)
            a = Vec::Scalar(p, r) / ppScalar;      // (p_k-1, r_k-1) / (p_k-1, p_k-1)

            for (uint16_t i = 0; i < size; i++)
            {
               x[i] += a * z[i];                   // [x_k] = [x_k-1] + a*z_k-1
               r[i] -= a * p[i];                   // [r_k] = [r_k-1] - a*p_k-1
            }

            lu.USlauSolve(r, Ar);
            A.MultToVec(Ar, tmp);
            lu.LSlauSolve(tmp, Ar);                // Ar = L^-1 * A * U^-1 * r
            //Vec::Mult(D, r, tmp);
            //A.MultToVec(tmp, Ar);
            //Vec::Mult(D, Ar, Ar);                  

            b = -Vec::Scalar(p, Ar) / ppScalar;    // b = - (p_k-1, L^-1 * A * U^-1 * r_k) / (p_k-1, p_k-1)
            lu.USlauSolve(r, tmp);                 // tmp = U^-1 * r_k

            for (uint16_t i = 0; i < size; i++)
            {
               z[i] = tmp[i] + b * z[i];            // [z_k] = U^-1 * r_k + b * [z_k-1]
               p[i] = Ar[i] + b * p[i];             // [p_k] = A * r_k + b * [p_k-1]
            }

            if (iter % resetIter == 0)
            {
               A.MultToVec(x, tmp);
               for (uint16_t i = 0; i < size; i++) tmp[i] = f[i] - tmp[i];
               lu.LSlauSolve(tmp, r);

               lu.USlauSolve(r, z);

               A.MultToVec(z, tmp);
               lu.LSlauSolve(tmp, p);
            }
            nev = Vec::Scalar(r, r);
            eps = sqrt(nev / ffScalar);

            // Выводим на то же место, что и раньше (со сдвигом каретки)
            if (debugOutput)
            {
               //cout << format("Итерация: {0:<10} относительная невязка: {1:<15.3e}\n", iter, eps);
               cout << format("\rИтерация: {0:<10} относительная невязка: {1:<15.3e}", iter, eps);
            }
            if (isinf(eps))
            {
               break;
            }
         }

         if (debugOutput)
         {
            cout << endl;
            if (isinf(eps))
            {
               cout << "Выход по переполнению метода" << endl << endl;
            }
            else if (iter > maxIter)
            {
               cout << "Выход по числу итераций" << endl << endl;
            }
            else
            {
               cout << "Выход по относительной невязке" << endl << endl;
            }
         }

         return iter - 1;
      }
   }

   void Destruct() {
      delete _tmp1, _tmp2, _tmp3, _tmp4, _tmp5, _tmp6;
      _tmp1 = _tmp2 = _tmp3 = _tmp4 = _tmp5 = _tmp6 = nullptr;
      delete _lu_mat;
      _lu_mat = nullptr;
   }
};
