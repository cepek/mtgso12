#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <matvec/matvec.h>
#include "icgs.h"

#include <gnu_gama/adj/adj_gso.h>

#include <chrono>
#include <memory>

using std::cout;
using std::endl;
using std::setw;

std::vector<GNU_gama::Mat<>> VB, VA;
std::vector<GNU_gama::Vec<>> Vb;
std::vector<std::string>     Vfiles;

int round_digits = 7;

double round(double d)
{
#if 1
  std::ostringstream ostr;
  ostr.precision(round_digits);
  ostr << std::scientific << d;
  return std::stod(ostr.str());
#else
  return d;
#endif
}

GNU_gama::Mat<> randMat(std::string filename)
{
  GNU_gama::Mat<> t;

  std::ifstream inp("/tmp/" + filename);
  inp >> t;

  std::cout << t.rows() << " " << t.cols() << "\n\n" << t << "\n";

  return t;
}

int test_icgs(int index)
{
  auto A = VB[index];
  //std::cout << "file " << filename << "  " << A;

  int n = A.cols() - 1;
  int m = A.rows() - n;

  //std::cout << "M = " << m << "   N = " << n << "\n\n";

  ICGS icgs;
  icgs.reset(A, m, n);
  icgs.icgs1();
  icgs.icgs2();

#if 0
  std::cout << "\nfile      : " << Vfiles[index] << "\n";
  const double* rb = icgs.residuals_begin();
  const double* re = icgs.residuals_end();
  std::cout << "residuals : "; while (rb < re) std::cout << *rb++ << " "; std::cout << endl;

  const double* ub = icgs.unknowns_begin();
  const double* ue = icgs.unknowns_end();
  std::cout << "unknowns  : "; while (ub < ue) std::cout << *ub++ << " "; std::cout << endl;
#endif

  double max_error {0};
  const double* errb = icgs.unknowns_begin();
  const double* erre = icgs.unknowns_end();
  while (errb < erre)
    {
      double t = std::abs(1 - *errb++);
      if (t > max_error) max_error = t;
    }

  auto file = Vfiles[index];
  std::cout << file.substr(file.rfind('/')+1);
  std::cout << "  n = " << setw(2) << n << "  m = " << setw(2) << m ;
  if (icgs.defect())
    std::cout << "  defect " << icgs.defect() << endl;
  else
    std::cout << "  max error : " << std::scientific << max_error << endl;

  return 0;
}

#if 1
int test_gama(int index)
{
  //to A = VB[index];

  using namespace GNU_gama;
  Mat<> A = VA[index];
  Vec<> b = Vb[index];

  int n = A.cols();
  int m = A.rows();

  GNU_gama::AdjGSO<double, int, GNU_gama::Exception::matvec> gama;
  gama.reset(A,b);
  gama.min_x();
  gama.solve();

  auto unk = gama.unknowns();
  //cout << trans(unk) << trans(gama.residuals());

  double max_error {0};
  for (int i=1; i<=n; i++)
    {
      double t = std::abs(1 - unk(i));
      if (t > max_error) max_error = t;
    }

  auto file = Vfiles[index];
  std::cout << file.substr(file.rfind('/')+1);
  std::cout << "  n = " << setw(2) << n << "  m = " << setw(2) << m ;
  if (gama.defect())
    std::cout << "  defect " << gama.defect() << endl;
  else
    std::cout << "  max error : " << std::scientific << max_error << endl;


  return 0;
}
#endif

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

#if 1
int test_mtgso1(int index)
{
  auto A = VB[index];

  int n = A.cols() - 1;
  int m = A.rows() - n;

  int adim = (m+n)*(n+1);
  std::shared_ptr<double[]> A_mtgso1 ( new double[adim] );
  std::shared_ptr<double[]> cp (new double[m+n+1] );  // one base indexing
  std::shared_ptr<int[]> indx (new int[n+1] );
  std::shared_ptr<int[]> inds (new int[n+1] );

  int df, ier;

  int m1 = m;   // number of project equations
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
  for (int c=1; c<=n+1; c++)
    for (int r=1; r<=m+n; r++)

      {
        *a1++ = A(r,c);
      }

  mtgso1_(A_mtgso1.get(), &adim, &m1, &n1, &n2, &n4, &mode, &sc, &eps,
          &tol, &nx,  &df, &ier, indx.get(), inds.get(), cp.get());


  double max_error {0};
  double* x = A_mtgso1.get() + (m+n)*n + m;
  for (int i=1; i<=n; i++)
    {
      double t = std::abs(1 - *x++);
      if (t > max_error) max_error = t;
    }

  auto file = Vfiles[index];
  std::cout << file.substr(file.rfind('/')+1);
  cout << "  n = " << setw(2) << n << "  m = " << setw(2) << m ;
  if (df)
    std::cout << "  defect " << df << endl;
  else
    std::cout << "  max error : " << std::scientific << max_error << endl;

  return 0;
}
#endif

int main(int argc, char *argv[])
{
  std::cout.precision(9);

  if (argc != 3)
    {
      std::cerr << "Usage: director files_list.txt\n";
      return 1;
    }

  std::string   dir  (argv[1]);
  if (dir.back() != '/') dir.append("/");
  std::ifstream files(dir + argv[2]);
  std::string filename;
  while (files >> filename) Vfiles.push_back(filename);
  std::sort(Vfiles.begin(), Vfiles.end());

  cout << "rhs elements rounded to " << round_digits << " digits\n\n";

  for (const auto& filename : Vfiles) {
      std::cout << "reading " << dir+filename << "\n" << std::flush;
      std::ifstream input_data(dir + filename);
      GNU_gama::Mat<> A, B;
      GNU_gama::Vec<> b;
      input_data >> B >> A >> b;
      for (int r=1; r<=B.rows()-B.cols()+1; r++)
        {
          B(r,B.cols()) = round(B(r,B.cols()));
          b(r) = round(b(r));
        }
      VB.push_back(B);
      VA.push_back(A);
      Vb.push_back(b);
    }

  if (VB.size() == 0)
    {
      std::cerr << "no input data read *******************\n";
      return 2;
    }

#if 1
  cout << "\nICGS\n\n";
  auto start = std::chrono::steady_clock::now();
  for (int i=0; i<VB.size(); i++)
    {
      test_icgs(i);
    }
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";
#endif

#if 1
  cout << "\nGNU_gama\n\n";

  start = std::chrono::steady_clock::now();
  for (int i=0; i<VB.size(); i++)
    {
      test_gama(i);
    }
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";
#endif

#if 1
  cout << "\nmtgso1\n\n";

  start = std::chrono::steady_clock::now();
  for (int i=0; i<VB.size(); i++)
    {
      test_mtgso1(i);
    }
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";
#endif

//  cout << "round example " << round(1.23456789) << endl;
//  cout << "round example " << round(0.123456789) << endl;
//  cout << "round example " << round(12345.6789) << endl;

  return 0;
}
