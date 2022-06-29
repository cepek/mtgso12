#include <iostream>
#include <iomanip>
#include <random>
#include <set>
#include <vector>
#include <algorithm>
#include <memory>
#include <matvec/matvec.h>
#include <gnu_gama/adj/adj_gso.h>
#include <gnu_gama/adj/icgs.h>

using std::cout;
using std::endl;

/* Compare results of four implementations of Gram-Schmidt Orthogonalization
 * -------------------------------------------------------------------------
 *
 * Two fortran subroutines mtgso1 and mtgso2 implements modified
 * Gram-Schmidt, the second one is optimized for 2D surveying networks
 * and swap colums corresponding to point coordinates x and y to improve
 * numerical stability.
 *
 * C++ templates class GSO from the project GNU Gama
 * https://www.gnu.org/software/gama/. Clone Gama project from git server:
 * git clone git clone https://git.savannah.gnu.org/git/gama.git
 *
 * C++ class ICGS from this directory implementing "Iterated Classical
 * Gram-Schmidt algorithm" (ICGS)
 */


/* Parameters of fortran subroutines
   .................................

      SUBROUTINE MTGSO1(A,ADIM,M1,N1,N2,N4,MODE,SC,EPS,TOL,
     /                  NX,DF,PROC,IER,INDX,INDS,CP)
      IMPLICIT INTEGER (A-Z)
      INTEGER  ADIM
      INTEGER  M1,N1,N2,N4,MODE,NX,DF,IER,INDX(1),INDS(1)
      LOGICAL  SC
      REAL     A(ADIM),EPS,TOL
      DOUBLE PRECISION CP(1)
      EXTERNAL PROC

  Input/output paramaters: A, DF, INDS
  Output parameter: IER
  Input parameters: ADIM, M1, N1, N2, N4, MODE, SC, EPS, TOL NX, INDX
  Other parameters: CP
*/

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

  void mtgso2_ (double A[],
                int*  adim,
                int* m1, int* n1, int* n2, int* n4, int* mode,
                bool* sc,
                double* eps, double* tol,
                int* nx, int* df,
                // EXTERNAL PROC ... removed
                int* ier, int indx[], int inds[],
                double cp[]);
}

GNU_gama::Vec<> unknowns(double* A, int m, int n)
{
  GNU_gama::Vec<> v(n);

  double* a = A + n*(m+n)+m;
  for (int i=1; i<=n; i++)
    v(i) = *a++;

  return v;
}

GNU_gama::Vec<> residuals(double* A, int m, int n)
{
  GNU_gama::Vec<> v(m);

  double* a = A + n*(m+n);
  for (int i=1; i<=m; i++)
    v(i) = *a++;

  return v;
}

GNU_gama::Mat<> blockMat(double* A, int m, int n)
{
  GNU_gama::Mat<> t(m+n, n+1);

  double* a = A;
  for (int j=1; j<=n+1; j++)
    for (int i=1; i<=m+n; i++)
      t(i,j) = *a++;

  return t;
}

GNU_gama::Mat<> BMmatrix(GNU_gama::Mat<> A, int m, int n)
{
  GNU_gama::Mat<> t(m,n);

  for (int i=1; i<=m; i++)
    for (int j=1; j<=n; j++)
      {
        t(i,j) = A(i,j);
      }

  return t;
}


GNU_gama:: Vec<> BMvector(GNU_gama::Mat<> A, int m, int n)
{
  GNU_gama::Vec<> v(m);

  for (int i=1; i<=m; i++)
      {
        v(i) = A(i,n+1);
      }

  return v;
}

int random_main()
{
  cout << "\n\n"
            << "**************************************************\n"
            << "*                                                *\n"
            << "*   Compare GNU Gama and Geodet.pc GSO results   *\n"
            << "*                                                *\n"
            << "**************************************************\n"
            << "\n";

  using namespace GNU_gama;

    static std::random_device rd;
    static std::mt19937 gen(rd());

#if 1
    static std::uniform_int_distribution<> uniformM(60, 70);
    static std::uniform_int_distribution<> uniformN(40, 55);
#else
    static std::uniform_int_distribution<> uniformM(6, 9);
    static std::uniform_int_distribution<> uniformN(4, 6);
#endif
    int M = uniformM(gen);
    int N = uniformN(gen);

    while (N > M)
      {
        M = uniformM(gen);
        N = uniformN(gen);
      }

    static std::uniform_real_distribution<> uniformA(-1, 1);

    Mat<> BM_gso(M+N,N+1);
    std::set<int> BM_indx;

    // GSO block matrix BM_gso
    {
      for (int i=1; i<=M; i++)
        for (int j=1; j<=N+1; j++)
          {
            BM_gso(i,j) = uniformA(gen);
          }

      // random defect
      static std::uniform_int_distribution<> uniformD(0,  6);
      int defect = uniformD(gen);

      std::vector<int> rind;
      for (int i=1; i<=N; i++) rind.push_back(i);
      std::random_shuffle(rind.begin(), rind.end());

      cout << "A size = " << M << " x " << N << "\n";

      if (defect) cout << endl;

      for (int d=1; d<=defect && rind.size() >= 3; d++)
        {
          if (rind.size() < 3) break;
          int i = rind[0];      // column to be rewritten
          int j = rind[1];
          int k = rind[2];
          BM_indx.insert(j);    // j, k regularization unknowns
          BM_indx.insert(k);

          std:: cout << "A[*,"<< i <<"] = A[*,"
                     <<  j <<"] + A[*, "<< k <<"]\n";

          for (int r=1; r<=M; r++)
            {
              BM_gso(r,i) = BM_gso(r,j) + BM_gso(r,k);
            }
          rind.erase(rind.begin(), rind.begin() + 3);
        }

      for (int i=1; i<=N; i++)
        for (int j=1; j<=N+1; j++)
          {
            BM_gso(M+i,j) = i==j ? 1 : 0;
          }
    } // GSO block matrix BM_gso

    bool full_print = false;
    /**********************/

#if 0
    {
      full_print = true;
      M = 5;   N = 3;
      Mat<> Demo1
               {{19, 12, 13, 14},
                {21, 29, 23, 24},
                {31, 32, 39, 34},
                {41, 42, 43, 44},
                {51, 52, 53, 54},
                { 1,  0,  0,  0},
                { 0,  1,  0,  0},
                { 0,  0,  1,  0}};
      BM_gso = Demo1;
    }
#endif

    if (full_print)
      cout << "\nVstupni blokova matice BM_gso\n"
           <<   "=============================\n\n"
           << std::setw(8) << std::setprecision(3) << std::fixed
           << BM_gso;

  Mat<> BM_icgs(BM_gso.rows(), BM_gso.cols());
  for (int i=1; i<=BM_gso.rows(); i++)
    for (int j=1; j<=BM_gso.cols(); j++)
      BM_icgs(i,j) = BM_gso(i,j);

  int adim = (M+N)*(N+1);
  std::shared_ptr<double[]> A_gso1 ( new double[adim] );
  std::shared_ptr<double[]> A_gso2 ( new double[adim] );

  int m1 = M;   // number of project equations
  int n1 = N;   // number of unknowns
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

  int df1, df2;
  int ier1, ier2;

  std::shared_ptr<double[]> cp ( new double[M+N+1] );
  std::shared_ptr<int[]>  inds1( new int[N+1] );
  std::shared_ptr<int[]>  inds2( new int[N+1] );
  std::shared_ptr<int[]>  indx ( new int[N+1] );

  nx = 0;
  for (int k : BM_indx) indx[nx++] = k;
  {
    double* a1 = A_gso1.get();
    double* a2 = A_gso2.get();
    for (int j=1; j<=N+1; j++)
      for (int i=1; i<=M+N; i++)
        {
          *a1++ = BM_gso(i,j);
          *a2++ = BM_gso(i,j);
        }
  }
  mtgso1_(A_gso1.get(), &adim, &m1, &n1, &n2, &n4, &mode, &sc, &eps,
          &tol, &nx,  &df1, &ier1, indx.get(), inds1.get(), cp.get());

  n4 = N%2;
  mtgso2_(A_gso2.get(), &adim, &m1, &n1, &n2, &n4, &mode, &sc, &eps,
          &tol, &nx,  &df2, &ier2, indx.get(), inds2.get(), cp.get());

    auto unk_gso1 = unknowns (A_gso1.get(), M, N);
    auto res_gso1 = residuals(A_gso1.get(), M, N);

    auto unk_gso2 = unknowns (A_gso2.get(), M, N);
    auto res_gso2 = residuals(A_gso2.get(), M, N);

    //cout << "*** residuals 1 " << res_gso1
    //cout << "*** unknowns  1 " << unk_gso1;

    auto res_diff12 = res_gso1 - res_gso2;
    auto unk_diff12 = unk_gso1 - unk_gso2;

    //cout << trans(res_diff12) << trans(unk_diff12);

    if (df1 > 0)
      {
         for (int i=0; i<df1; i++)
          {
            double* a = A_gso1.get();
            int k = inds1[i] - 1;
            a += (M+N)*k;
            for (int j=1; j<=(M+N); j++) *a++ = 0;
          }
       }
    if (df2 > 0)
      {
         for (int i=0; i<df2; i++)
          {
            double* a = A_gso2.get();
            int k = inds2[i] - 1;
            a += (M+N)*k;
            for (int j=1; j<=(M+N); j++) *a++ = 0;
          }
      }
#if 0
    if (df1 > 0)
      {
        cout << "\ndf1 = " << df1 << " : ";
        for (int i=0; i<df1; i++)
          {
            cout << inds1[i] << " ";
          }
        cout << endl;
      }
    if (df2 > 0)
      {
        cout << "df2 = " << df2 << " : ";
        for (int i=0; i<df2; i++)
          {
            cout << inds2[i] << " ";
          }
        cout << endl;
      }
#endif

#if 0
    cout << "\nBlock matrices A_gso1 and A_gso2 after ortogonalization\n"
              <<   "--------------------------------------------------\n\n";
    cout << std::setw(8) << std::setprecision(3) << std::fixed
              << blockMat(A_gso1.get(), M, N) << "\n"
              << std::setw(8) << std::setprecision(3) << std::fixed
              << blockMat(A_gso2.get(), M, N) << "\n\n"
              << endl;
#endif

    //cout << BMmatrix(BM_gso, M, N);
    //cout << BMvector(BM_gso, M, N);
    auto bmm = BMmatrix(BM_gso, M, N);
    auto bmv = BMvector(BM_gso, M, N);

    // cout << "Matice vstup gama " << bmm << bmv;

    GNU_gama::AdjGSO<double, int, GNU_gama::Exception::matvec> gama;
    gama.reset(bmm, bmv);
    gama.min_x(nx, indx.get());
    // cout << "nx = " << nx <<" gama.minx = ";
    // for (int i=0; i<nx; i++) cout << indx[i] << " ";
    // cout << endl;
    gama.solve();

    //if (full_print)
    //  cout << "gama.solve() "
    //       << std::setw(8) << std::setprecision(13) << std::fixed << gama.A_;

    auto unk_diff1g = unk_gso1 + gama.unknowns();
    auto unk_diff2g = unk_gso2 + gama.unknowns();

    auto res_diff1g = res_gso1 + gama.residuals();
    auto res_diff2g = res_gso2 + gama.residuals();



    //---cout << "----------------------------------------------\n";
    //---cout << "BM_icgs " << BM_icgs;
    ICGS icgs; //(BM_icgs);

    std::shared_ptr<double[]> temp ( new double[(M+N)*(N+1)] );
    static int icgs_count {0};
    if (++icgs_count % 2)
      {
        cout << "\nicgs uses external memory\n";
        double* a = temp.get();
        for (int c=1; c<=BM_icgs.cols(); c++)
          for (int r=1; r<=BM_icgs.rows(); r++)
            *a++ = BM_icgs(r,c);

        icgs.reset(temp.get(), M, N);
        icgs.min_x(nx, indx.get());
        icgs.icgs1();
        icgs.icgs2();

        a = temp.get();
        for (int c=1; c<=BM_icgs.cols(); c++)
          for (int r=1; r<=BM_icgs.rows(); r++)
            BM_icgs(r,c) = *a++;
      }
    else
      {
        cout << "\nicgs uses internal memory\n";
        icgs.reset(BM_icgs,BM_icgs.rows()-BM_icgs.cols()+1, BM_icgs.cols()-1);
        icgs.min_x(nx, indx.get());
        icgs.icgs1();
        icgs.icgs2();
      }


    //cout << "icgs defect " << icgs.defect() << "\n";
    //cout << "icgs error  " << icgs.error() << "\n";


    Vec<> icgs_res(M);
    Vec<> icgs_unk(N);
    {
      int i = 1;
      for (const double* r=icgs.residuals_begin(); r<icgs.residuals_end(); r++)
        {
          icgs_res(i++) = *r;
        }

      i = 1;
      for (const double* x=icgs.unknowns_begin(); x<icgs.unknowns_end(); x++)
        {
          icgs_unk(i++) = *x;
        }
    }
    //cout << trans(icgs_res) << trans(icgs_unk);

    auto unk_diffcg = icgs_unk + gama.unknowns();
    auto res_diffcg = icgs_res + gama.residuals();

    cout << "==============================================\n";

#if 0
    int i = 1;
    while (res.first != res.second) icgs_res(i++) = *res.first++;
    auto unk = icgs.unknowns();
    i = 1;
    while (unk.first != unk.second) icgs_unk(i++) = *res.first++;
    cout << trans(icgs_res);// << trans(gama.residuals());
    //cout << trans(icgs_unk);// << trans(gama.unknowns());
#endif


    auto print = [](const char* text, Vec<> v, bool f)
    {
        cout << text << v.norm_L1();
        if (f) cout << "  ... "<<  trans(v);
        else   cout << "\n";
    };

    cout << "\n   RESIDUALS\n";
    print("diff 1-2  L1 = ", res_diff12, full_print);
    print("diff 1-g  L1 = ", res_diff1g, full_print);
    print("diff 2-g  L1 = ", res_diff2g, full_print);
    print("diff c-g  L1 = ", res_diffcg, full_print);

    cout << "\n   UNKNOWNS\n";
    print("diff 1-2  L1 = ", unk_diff12, full_print);
    print("diff 1-g  L1 = ", unk_diff1g, full_print);
    print("diff 2-g  L1 = ", unk_diff2g, full_print);
    print("diff c-g  L1 = ", unk_diffcg, full_print);
    //print("diff 2-c  L1 = ", unk_diff2c, full_print);

    if (gama.defect())
      {
        const auto& x_gso1 = unk_gso1;
        const auto& x_gso2 = unk_gso1;
        const auto& x_gama = gama.unknowns();

        Vec<> x_icgs(x_gama.dim());
        auto ptr = icgs.unknowns_begin();
        auto end = icgs.unknowns_end();
        int k=1;
        while (ptr < end) x_icgs(k++) = *ptr++;

        auto B1 = blockMat(A_gso1.get(), M, N);
        auto B2 = blockMat(A_gso2.get(), M, N);
        double t, norm_gso1 = 0, norm_gso2 = 0, norm_gama = 0, norm_icgs = 0;
        double sum_qxx_gso1 = 0, sum_qxx_gso2 = 0, sum_qxx_gama = 0;

        double sum_qxx_icgs = 0;
        Mat<> Bc;
        icgs.getMat(Bc);

        cout << "\nRegularization: nx = " << nx << "  indx {";
        for (int j = 0; j<nx; j++)
          {
            int i = indx[j];
            cout << " " << i;

            t = x_gso1(i);
            norm_gso1 += t*t;
            t = x_gso2(i);
            norm_gso2 += t*t;
            t = x_gama(i);
            norm_gama += t*t;
            t = x_icgs(i);
            norm_icgs += t*t;

            for (int s=1; s<=N; s++)
              {
                t = B1(M+i,s);
                sum_qxx_gso1 += t*t;
                t = B2(M+i,s);
                sum_qxx_gso2 += t*t;
                t = Bc(M+i,s);
                sum_qxx_icgs += t*t;
              }
            sum_qxx_gama += gama.q_xx(i,i);
          }

        norm_gso1 = std::sqrt(norm_gso1);
        norm_gso2 = std::sqrt(norm_gso2);
        norm_gama = std::sqrt(norm_gama);
        norm_icgs = std::sqrt(norm_icgs);
        cout << " }\n\n";

        auto prec = cout.precision();
        cout.precision(6);
        cout << "norm_gso1 dx   " << norm_gso1
             << "      sum qx   " << sum_qxx_gso1 << endl;
        cout << "norm_gso2 dx   " << norm_gso2
             << "      sum qx   " << sum_qxx_gso2 << endl;
        cout << "norm_gama dx   " << norm_gama
             << "      sum qx   " << sum_qxx_gama << endl;
        cout << "norm_icgs dx   " << norm_icgs
             << "      sum qx   " << sum_qxx_icgs << endl;
        cout << endl;
        cout.precision(prec);
      }

    {
      if (gama.defect() == 0)  cout << endl;

      //cout << BMmatrix(BM_gso, M, N);
      cout << "gama    defect " << gama.defect() << endl;
      cout << "icgs    defect " << icgs.defect() << endl;
      cout << "gso 1 2 defect " << df1  << " " << df2  << "\n";
      cout << "gso 1 2 ier    " << ier1 << " " << ier2 << "\n\n";

        if (df1 || df2)
          {
            auto s1 = inds1.get();
            cout << "inds1 lin.dep. { ";
            for (int i=0; i<df1; i++) cout << *s1++ << " ";
            cout << "}\n";

            auto s2 = inds2.get();
            cout << "inds2          { ";
            for (int i=0; i<df2; i++) cout << *s2++ << " ";
            cout << "}\n";

            cout << "gama           { ";
            for (int i=1; i<=N; i++) if (gama.lindep(i)) cout << i << " ";
            cout << "}\n";

            cout << "icgs           { ";
            for (int n : icgs.lindep_columns) cout << n << " ";
            cout << "}\n";
          }
    }

    double l1 = std::max(1.0, gama.unknowns().norm_L1());
    const double tol1 = 1e-10;

    if (res_diff1g.norm_L1()/l1 > tol1)
      {
        cout << "\ndiff gso1 gama relative error "
             << std::scientific << res_diff1g.norm_L1() << " / " << l1 << endl;
        return 1;
      }
    if (false && unk_diff2g.norm_L1()/l1 > tol1)
      {
        cout << "\ndiff gso2 gama relative error "
             << std::scientific << unk_diff2g.norm_L1() << " / " << l1 << endl;
        return 2;
      }

    int repeat = 10;

    if (res_diff12.norm_L1() > 1e-8) repeat = 0;
    if (unk_diff12.norm_L1() > 1e-8) repeat = 0;

    return repeat;
}

int main()
{
  int repeat = 100;
  do {
    } while (repeat-- && random_main());

  cout << "\nexit status = " << repeat+1 << endl;

  return 0;
}
