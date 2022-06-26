#include <iostream>
#include <matvec/matvec.h>
#include "icgs.h"

#include <gnu_gama/adj/adj_gso.h>

#include <chrono>
#include <memory>
#include <random>

using std::cout;
using std::endl;

const double maxn = 1000;

#define NUMERICAL_NOISE
#ifdef NUMERICAL_NOISE

static std::random_device rd;
static std::mt19937 gen(rd());

const double eps=1e-12;
static std::uniform_real_distribution<> uniform_noise(-eps, eps);

const int noise_dim = (2*maxn+1)*(maxn+1);
double noise[noise_dim];

#endif



const double rho = 1e-9;

GNU_gama::Mat<> Lauchli(int n)
{
  GNU_gama::Mat<> t(2*n+1,n+1);

  for (int i=1; i<=2*n+1; i++)
    for (int j=1; j<=n+1; j++)
      {
        t(i,j) = 0;
      }

  t(1,n+1) = n*rho;
  for (int c=1; c<=n; c++) t(1,c) = 1;
  for (int r=2; r<=n+1; r++) t(r,r-1) = 1;
  for (int r=2; r<=n+1; r++) t(r,n+1) = rho;
  for (int r=1; r<=n; r++)   t(n+r+1,r) = 1;

#ifdef NUMERICAL_NOISE
  int m = 0;
  for (int i=1; i<=2*n+1; i++)
    for (int j=1; j<=n+1; j++)
      {
        t(i,j) += noise[m++];
      }
#endif

  return t;
}

int test_icgs(int n)
{
  auto A = Lauchli(n);
  // cout << A;

  ICGS icgs;
  icgs.reset(A, n+1, n);
  icgs.icgs1();

#if 0
  const double* rb = icgs.residuals_begin();
  const double* re = icgs.residuals_end();
  cout << "residuals : "; while (rb < re) cout << *rb++ << " "; cout << endl;

  const double* ub = icgs.unknowns_begin();
  const double* ue = icgs.unknowns_end();
  cout << "unknowns  : "; while (ub < ue) cout << *ub++ << " "; cout << endl;
#endif

  double max_error {0};
  const double* errb = icgs.unknowns_begin();
  const double* erre = icgs.unknowns_end();
  while (errb < erre)
    {
      double t = std::abs(rho + *errb++);
      if (t > max_error) max_error = t;
    }

  cout << "n = " << n << "  max error : " << max_error << endl;
  return 0;
}

int test_gama(int n)
{
  auto L = Lauchli(n);

  using namespace GNU_gama;
  Mat<> A(n+1, n);
  Vec<> b(n+1);

  for (int r=1; r<=n+1; r++)
    {
      for (int c=1; c<=n; c++) A(r,c) = L(r,c);

      b(r) = L(r, n+1);
    }

  GNU_gama::AdjGSO<double, int, GNU_gama::Exception::matvec> gama;

  gama.reset(A,b);
  gama.min_x();
  gama.solve();

  auto unk = gama.unknowns();
  //cout << trans(unk) << trans(gama.residuals());

  double max_error {0};
  for (int i=1; i<=n; i++)
    {
      double t = std::abs(rho - unk(i));
      if (t > max_error) max_error = t;
    }

  cout << "n = " << n << "  max error : " << max_error << endl;

  return 0;
}

extern "C" {
  void mtgso1_ (double A[],
                int*  adim,
                int* m1, int* n1, int* n2, int* n4, int* mode,
                bool* sc,
                double* eps, double* tol,
                int* nx, int* df,
                // EXTERNAL PROC ... removed
                int* ier, int indx[], int inds[],
                double cp[]);
}

int test_mtgso1(int n)
{


  auto A = Lauchli(n);
  // cout << "********************************  " << A;

  int adim = (2*n+1)*(n+1);
  std::shared_ptr<double[]> A_mtgso1 ( new double[adim] );
  std::shared_ptr<double[]> cp (new double[2*n+1] );  // one base indexing
  std::shared_ptr<int[]> indx (new int[n] );
  std::shared_ptr<int[]> inds (new int[n] );

  int df, ier;

  int m1 = n+1; // number of project equations
  int n1 = n;   // number of unknowns
  int n2 = 1;   // number of rhs vectors
  int n4 = 0;   // number of apriori orthogonal columns
                //    in the left part of the matrix of project equations
  int mode = 2; // solution variant
  bool sc = 1;  // inputn normalization
  double eps = 1e-16; // relative precision of the computer
  double tol = 1e-8;  // tolerance for detection of linearly
                      // dependent columns in the orthonslizes
                      // matrices
  int nx = 0;   // number of elements of the subvector of unknowns
                // whose length is to be minimized

  for (int k=1; k<=n; k++) indx[nx++] = k;
  double* a1 = A_mtgso1.get();
  for (int j=1; j<=n+1; j++)
    for (int i=1; i<=2*n+1; i++)
      {
        // cout << i << " " << j << " " << A(i,j) << endl;
        *a1++ = A(i,j);
      }

  mtgso1_(A_mtgso1.get(), &adim, &m1, &n1, &n2, &n4, &mode, &sc, &eps,
          &tol, &nx,  &df, &ier, indx.get(), inds.get(), cp.get());




  double max_error {0};
  double* x = A_mtgso1.get() + n*(2*n+1) + (n+1);
  for (int i=1; i<=n; i++)
    {
      //cout << "i/x_i " << i << "/" << *x << " ";
      double t = std::abs(rho + *x++);
      if (t > max_error) max_error = t;
    }

  cout << "n = " << n << "  max error : " << max_error << endl;
  return 0;
}

int main(int argc, char *argv[])
{
#ifdef NUMERICAL_NOISE
  {
    int n = 0;
    while (n < noise_dim) noise[n++] = uniform_noise(gen);
  }
#endif

  cout << "\nICGS\n\n";
  auto start = std::chrono::steady_clock::now();
  int n = maxn;
  while (n >= 1)
    {
      test_icgs(n);
      n *= 0.38;
    }
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";

#if 1
  cout << "\nGNU_gama\n\n";
  n = maxn;
  start = std::chrono::steady_clock::now();
  while (n >= 1)
    {
      test_gama(n);
      n *= 0.38;
    }
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";
#endif

  cout << "\nmtgso1\n\n";
  n = maxn;
  start = std::chrono::steady_clock::now();
  while (n >= 1)
    {
      test_mtgso1(n);
      n *= 0.38;
    }
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";


  return 0;
}
