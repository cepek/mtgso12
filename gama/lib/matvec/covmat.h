/*
  C++ Matrix/Vector templates (GNU Gama / matvec)
  Copyright (C) 2002, 2005, 2007, 2018  Ales Cepek <cepek@gnu.org>

  This file is part of the GNU Gama C++ Matrix/Vector template library.

  This library is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with GNU Gama.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GNU_gama_CovMat_gMatVec_Symmetric_Band_Matrix_2_h
#define GNU_gama_CovMat_gMatVec_Symmetric_Band_Matrix_2_h

#include <matvec/matvec.h>
#include <matvec/choldec.h>
#include <algorithm>

/*
 * CovMat - Covariance Matrix (symmetric band matrix)
 * ==================================================
 *
 * Bandwidth is defined as max{ |i-j| | a_ij != 0 }
 *
 * Upper triangular part of the matrix is stored by rows, ie
 *
 *        d*(b+1) - b*(b+1)/2   of unzero elements,
 *
 * where `d' is the matrix dimension and `b' the bandwidth
 *
 */


namespace GNU_gama {   /** \brief Covariance Matrix (symmetric band matrix) */


  template <typename Float=double,
            typename Index=int,
            typename Exc=Exception::matvec>
  class CovMat : public MatBase<Float, Index, Exc>,
                 public CholDecLD<Float, Index, Exc> {
  public:

    CovMat() : band_(0), band_1(0), dim_b(0)
    {
    }
    CovMat(Index d, Index b)
      : MatBase<Float, Index, Exc>(d,d,d*(b+1) - b*(b+1)/2), band_(b),
      band_1(b+1), dim_b(d-b)
    {
    }

    void reset()
    {
      this->row_ = this->col_ = band_ = band_1 = dim_b = 0;
      this->resize(0);
    }

    void   reset    (Index d, Index b) override;
    Index  dim      () const { return this->row_; }
    Index  bandWidth() const { return band_; }
    Float  operator ()(Index, Index) const override;
    Float& operator ()(Index, Index) override;
    void   cholDec  () override;
    void   solve    (Vec<Float, Index, Exc>&) const override;

    const Float* operator[](Index row) const
    {
      const Float* a_ = this->begin() + --row*band_1;
      if (row > dim_b) {
        const Index i_  = row - dim_b;
        a_ -= i_*(i_+1)/2;
      }
      return a_;
    }

    Float* operator[](Index row)
    {
      Float* a_ = this->begin() + --row*band_1;
      if (row > dim_b)
        {
          const Index i_  = row - dim_b;
          a_ -= i_*(i_+1)/2;
        }
      return a_;
    }

    Vec<Float, Index, Exc> operator*(const Vec<Float, Index, Exc>&) const;

    std::istream&  read (std::istream&) override;
    std::ostream&  write(std::ostream&) const override;

  private:

    Index   band_, band_1, dim_b;

  };      /* class CovMat */


  template <typename Float, typename Index, typename Exc>
  void
  CovMat<Float, Index, Exc>::reset(Index d, Index b)
  {
    if (dim() != d || band_ != b)
      {
        this->row_ = this->col_ = d;
        band_  = b;
        band_1 = b+1;
        dim_b  = d-b;
        this->resize(d*(b+1) - b*(b+1)/2);
      }
  }


  template <typename Float, typename Index, typename Exc>
  Float
  CovMat<Float, Index, Exc>::operator()(Index r, Index s) const
  {
    if (r > s) {
      Index t = r;
      r = s;
      s = t;
    }

    if (s > r+band_)
      return 0;

    s -= r;
    return *(operator[](r) + s);
  }


  template <typename Float, typename Index, typename Exc>
  Float&
  CovMat<Float, Index, Exc>::operator()(Index r, Index s)
  {
    if (r > s) {
      Index t = r;
      r = s;
      s = t;
    }

    if (s > r+band_)
      throw Exc(Exception::BadIndex,
                "Float& CovMat::operator()(Index r, Index s)");

    s -= r;
    return *(operator[](r) + s);
  }


  template <typename Float, typename Index, typename Exc>
  void
  CovMat<Float, Index, Exc>::cholDec()
  {
    /*
     * Cholesky factorization of positive definite matrix A = L*D*trans(L)
     *
     * L is lower triangular matrix with unity diagonal; D is diagonal matrix.
     * Matrices L and D replace factored band symmetric matrix `in situ'.
     */
    using namespace std;

    Float *B = this->begin();
    Index  N = dim();
    Index  W = bandWidth();

    const  Float  Tol = this->Abs(*B*this->cholTol());
    Float *p;
    Index  row, k, l, n;
    Float  pivot, q;

    if (N == 0)
      throw Exc(Exception::BadRank,
                "CovMat::cholDec(Float  tol) - zero dim matrix");

    for (row=1; row<=N; row++)
      {
        if((pivot = *B) < Tol)
          throw Exc(Exception::NonPositiveDefinite,
                    "CovMat::cholDec(Float  tol) - "
                    "Matrix is not positive definite");

        k = min(W, N-row);             // number of of-diagonal elements
        p = B+k;                       // next row address - 1
        for (n=1; n<=k; n++)
          {
            q = B[n]/pivot;
            for (l=n; l<=k; l++) p[l] -= q*B[l];
            p += min(W, N-row-n);
          }

        B++;                           // *B++ = pivot = sqrt(pivot);
        for (; k; k--) *B++ /= pivot;
      }
  }


  template <typename Float, typename Index, typename Exc>
  void
  CovMat<Float, Index, Exc>::solve(Vec<Float, Index, Exc>& rhs) const
  {
    using namespace std;
    Index i, j, k;
    Float s;
    const Float *m;

    // forward substitution
    for (i=2; i<=dim(); i++)
      {
        s = 0;
        for (j = i>band_ ? i-band_ : 1; j<i; j++) s += operator()(i,j)*rhs(j);
        rhs(i) -= s;
      }

    // inverse of diagonal
    for (i=1; i<=dim(); i++) rhs(i) /= *operator[](i);

    // backward substituiton
    for (i=dim()-1; i>0; i--)
      {
        s = 0;
        m = operator[](i) + 1;
        for (k=i+1; k<=min(i+band_,dim()); k++) s += *m++ * rhs(k);
        rhs(i) -= s;
      }
  }


  template <typename Float, typename Index, typename Exc>
  Vec<Float, Index, Exc>
  CovMat<Float, Index, Exc>::operator*(const Vec<Float, Index, Exc>& v) const
  {
    Vec<Float, Index, Exc> T(dim());
    Index i, j, k;
    Float s;

    for (i=1; i<=dim(); i++)
      {
        s = 0;

        k = i+band_;
        if (k > dim()) k = dim();
        for (j = i>band_ ? i-band_ : 1; j<=k; j++)
          {
            s += operator()(i,j)*v(j);
          }

        T(i) = s;
      }

    return T;
  }


  template <typename Float, typename Index, typename Exc>
  std::istream&
  CovMat<Float, Index, Exc>::read(std::istream& inp)
  {
    int inpd, inpb;
    inp >> inpd >> inpb;
    reset(inpd, inpb);

    Float  *b = this->begin();
    for (Index i=1; i<=dim(); i++)
      for (Index j=i; j<=i+band_; j++)
        if (j <= dim())
          inp >> *b++;

    return inp;
  }


  template <typename Float, typename Index, typename Exc>
  std::ostream&
  CovMat<Float, Index, Exc>::write(std::ostream& out) const
  {
    auto w = out.width();
    out.width(w);
    out << dim() << ' ';
    out.width(w);
    out << band_ << "\n\n";

    const Float  *b = this->begin();
    for (Index i=1; i<=dim(); i++, out << '\n')
      for (Index j=i; j<=i+band_; j++)
        if (j <= dim())
          {
            out.width(w);
            out << *b++ << ' ';
          }

    return out;
  }


}      //  namespace GNU_gama

#endif
