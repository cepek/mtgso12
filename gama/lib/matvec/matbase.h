/*
  C++ Matrix/Vector templates (GNU Gama / matvec)
  Copyright (C) 1999, 2007, 2018  Ales Cepek <cepek@gnu.org>

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

#ifndef GNU_gama_gMatVec_MatBase_h
#define GNU_gama_gMatVec_MatBase_h

#include <iostream>
#include <matvec/matvecbase.h>


namespace GNU_gama {   /** \brief Base matrix class */


  template <typename Float=double,
            typename Index=int,
            typename Exc=Exception::matvec>
  class MatBase : public MatVecBase<Float, Index, Exc>
  {
  protected:

    Index row_;
    Index col_;

    MatBase() : row_(0), col_(0) {}
    MatBase(Index r, Index c, Index nsz)
      : MatVecBase<Float, Index, Exc>(nsz), row_(r), col_(c) {}
    MatBase(Index r, Index c, const MatBase& m)
      : MatVecBase<Float, Index, Exc>(m), row_(r), col_(c) {}
    virtual ~MatBase() = default;

  public:

    using iterator = typename MatVecBase<Float, Index, Exc>::iterator;
    using const_iterator = typename MatVecBase<Float, Index, Exc>::const_iterator;

    Index rows() const { return row_; }
    Index cols() const { return col_; }

    Index min_rc() const { return row_ > col_ ? col_ : row_; }
    Index max_rc() const { return row_ < col_ ? col_ : row_; }

    virtual Float& operator()(Index r, Index c) = 0;
    virtual Float  operator()(Index r, Index c) const = 0;

    void reset() { row_ = col_ = 0; this->resize(0); }
    virtual void reset(Index r, Index c) {
      if (r != row_ || c != col_) {
        row_ = r; col_ = c; this->resize(r*c);
      }
    }

    void set_diagonal(Float d)
    {
      this->set_zero();
      for (Index i=1; i<=min_rc(); i++)
        this->operator()(i,i) = d;
    }
    void set_identity() { set_diagonal(Float(1.0)); }

    virtual void transpose()
    {
      throw Exc(Exception::NotImplemented, "MatBase::transpose()");
    }

    virtual void invert()
    {
      throw Exc(Exception::NotImplemented, "MatBase::invert()");
    }

    virtual std::istream& read(std::istream& inp)
    {
      Index r, c;
      if (inp >> r >> c)
        {
          reset(r, c);
          for (Index i=1; i<=r; i++)
            for (Index j=1; j<=c; j++)
              inp >> operator()(i,j);
        }
      return inp;
    }
    virtual std::ostream& write(std::ostream& out) const
    {

      const auto fw = out.width();
      out.width(fw);
      out << rows() << " ";
      out.width(fw);
      out << cols() << "\n\n";
      for (Index i=1; i<=rows(); i++)
        {
          for (Index j=1; j<=cols(); j++) {
            out.width(fw);
            out << operator()(i,j) << " ";
          }
          out << '\n';
        }
      return out;
    }

  };


  template <typename Float, typename Index, typename Exc>
  std::istream& operator>>(std::istream& inp, MatBase<Float, Index, Exc>& M)
  {
    return M.read(inp);
  }


  template <typename Float, typename Index, typename Exc>
  std::ostream&
  operator<<(std::ostream& out, const MatBase<Float, Index, Exc>& M)
  {
    return M.write(out);
  }


}   // namespace GNU_gama

#endif
