/*************************************************************************** \
 * Copyright (C) by University Paris-Est - MISS team
 * Matrix4x4.hpp created in 09 2008.
 * Mail : biri@univ-mlv.fr
 *
 * This file is part of the OpenKraken-math.
 *
 * The OpenKraken-math is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * The OpenKraken-math is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.      If not, see <http://www.gnu.org/licenses/>.
 *
 \***************************************************************************/

/*
 * Anti-doublon
 */
#ifndef __OPENKN_MATH__MATRIX4X4_HPP__
#define __OPENKN_MATH__MATRIX4X4_HPP__

/*
 * Internal Includes
 */
#include "Matrix.hpp"
#include "Vector4.hpp"

namespace kn {

  template<class T>
    class Matrix4x4 : public Matrix<T> {

      protected:

      public:

        Matrix4x4(void)
          : Matrix<T>(4) {
            setIdentity();
          }

        Matrix4x4(const Matrix4x4<T>& m)
          : Matrix<T>(4) {
            std::copy(m.begin_, m.end_, this->begin_);
          }

        Matrix4x4(const Matrix<T>& m);

        Matrix4x4(const Matrix4x4<T>* m);

        Matrix4x4(const T& d)
          : Matrix<T>(4) {
            std::fill(this->begin_, this->end_, d);
          }

        Matrix4x4(const T* d);

        Matrix4x4(const Vector<T>& v,
            const bool& setasrows = true);

        Matrix4x4(const Vector4<T>* v,
            const bool& setasrows = true);

        ~Matrix4x4(void) {
        }

      public:


        using Matrix<T>::operator*;

        Vector4<T> operator*(const Vector4<T>& v) const;

        inline size_t rows(void) const {
          return 4;
        }

        inline size_t columns(void) const {
          return 4;
        }

        Matrix4x4<T>& transpose(void);

        void setIdentity(void);

        inline bool isSquare(void) const {
          return true;
        }

        Matrix4x4<T>& operator=(const Matrix4x4<T>& m);

        Matrix4x4<T> operator+(const Matrix4x4<T>& m) const;

        Matrix4x4<T> operator-(const Matrix4x4<T>& m) const;

        Matrix4x4<T> operator/(const T& d) const;

        Matrix4x4<T> operator*(const Matrix4x4<T>& m) const;

        Matrix4x4<T> operator*(const T& d) const;

        Matrix4x4<T> operator-(void) const;

        Matrix4x4<T>& operator+=(const Matrix4x4<T>& m);

        Matrix4x4<T>& operator-=(const Matrix4x4<T>& m);

        Matrix4x4<T>& operator/=(const T& d);

        Matrix4x4<T>& operator*=(const Matrix4x4<T>& m);

        Matrix4x4<T>& operator*=(const T& d);

        Vector4<T> getRow(const unsigned int& row) const;

        Vector4<T> getColumn(const unsigned int& column) const;

        Vector4<T> getDiagonal(void) const;

        inline Matrix4x4<T> getTranspose(void) const
        {
          return Matrix4x4<T>(*this).transpose();
        }

        Matrix4x4<T>& power(const unsigned int& p);

        T trace(void) const;

    };


  template<class T>
    Matrix4x4<T>::Matrix4x4(const Matrix4x4<T>* m)
    : Matrix<T>(4)
    {
      if(m==0)
        throw MathException("Pointer is null");
      std::copy(m->begin_, m->end_, this->begin_);
    }



  template<class T>
    Matrix4x4<T>::Matrix4x4(const T* a)
    : Matrix<T>(4)
    {
      if(a==0)
        throw MathException("Pointer is null");
      std::copy(a, a + 16, this->begin_);
    }


  template<class T>
    Matrix4x4<T>::Matrix4x4(const Vector<T>& v,
        const bool& setasrows)
    : Matrix<T>(4)
    {
      if(v.size()!=16)
        throw MathException("Vector size is different from Matrix size");

      if(setasrows)
        for(unsigned int i = 0; i < 16; ++i)
          (this->data)[i] = v[i];
      else
        for(unsigned int j = 0; j < 4; ++j)
          for(unsigned int i = 0; i < 4; ++i)
            (this->accessor)[i][j] = v[j*4+i];
    }


  template<class T>
    Matrix4x4<T>::Matrix4x4(const Vector4<T>* v,
        const bool& setasrows)
    : Matrix<T>(4)
    {

      if(v==0)
        throw MathException("Pointer is null");

      if(setasrows)
        for(unsigned int i = 0; i < 4; ++i)
          for(unsigned int j = 0; j < 4; ++j)
            (this->accessor)[i][j] = v[i][j];
      else
        for(unsigned int j = 0; j < 4; ++j)
          for(unsigned int i = 0; i < 4; ++i)
            (this->accessor)[i][j] = v[j][i];

    }


  template<class T>
    Matrix4x4<T>::Matrix4x4(const Matrix<T>& m)
    : Matrix<T>(4) {
      if ((m.rows() != 4) || (m.columns() != 4))
        throw MathException("Matrix4x4 Constructor","Matrix sizes are different");
      std::copy(m.getMatrixArray(), m.getMatrixArray()+m.rows()*m.columns(), this->begin_);
    }


  template<class T>
    Vector4<T> Matrix4x4<T>::operator*(const Vector4<T>& v) const {

      Vector<T> vtmp(4);
      vtmp[0] = v[0]*(this->accessor)[0][0] + v[1]*(this->accessor)[0][1] + v[2]*(this->accessor)[0][2] + v[3]*(this->accessor)[0][3];
      vtmp[1] = v[0]*(this->accessor)[1][0] + v[1]*(this->accessor)[1][1] + v[2]*(this->accessor)[1][2] + v[3]*(this->accessor)[1][3];
      vtmp[2] = v[0]*(this->accessor)[2][0] + v[1]*(this->accessor)[2][1] + v[2]*(this->accessor)[2][2] + v[3]*(this->accessor)[2][3];
      vtmp[3] = v[0]*(this->accessor)[3][0] + v[1]*(this->accessor)[3][1] + v[2]*(this->accessor)[3][2] + v[3]*(this->accessor)[3][3];
      return vtmp;
    }



  template<class T>
    Matrix4x4<T>& Matrix4x4<T>::transpose(void) {

      for(unsigned int i = 0; i < 4; ++i)
        for(unsigned int j = 0; j < i; ++j)
          std::swap((this->accessor)[i][j],(this->accessor)[j][i]);

      return *this;
    }

  template<class T>
    void Matrix4x4<T>::setIdentity(void){
      this->setZero();
      (this->accessor)[0][0] = (this->accessor)[1][1] = (this->accessor)[2][2] = (this->accessor)[3][3] = T(1.0);
    }




  template<class U>
    Vector4<U> operator* (const Vector4<U>& v,
        const Matrix4x4<U>& m){

      Vector4<U> result;
      result[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2] + m[3][0]*v[3];
      result[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2] + m[3][1]*v[3];
      result[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2] + m[3][2]*v[3];
      result[3] = m[0][3]*v[0] + m[1][3]*v[1] + m[2][3]*v[2] + m[3][3]*v[3];
      return result;
    }


  template<typename T>
    Matrix4x4<T>& Matrix4x4<T>::operator=(const Matrix4x4<T>& m){
      if(&m == this) return *this;

      std::copy(m.begin_,m.end_,this->begin_);

      return *this;
    }


  template<typename T>
    Matrix4x4<T> Matrix4x4<T>::operator+(const Matrix4x4<T>& m) const{
      Matrix4x4<T> result = *this;
      result += m;
      return result;
    }


  template<typename T>
    Matrix4x4<T> Matrix4x4<T>::operator-(const Matrix4x4<T>& m) const{
      Matrix4x4<T> result = *this;
      result -= m;
      return result;
    }


  template<typename T>
    Matrix4x4<T> Matrix4x4<T>::operator/(const T& d) const{
      Matrix4x4<T> result = *this;
      result /= d;
      return result;
    }


  template<typename T>
    Matrix4x4<T> Matrix4x4<T>::operator*(const Matrix4x4<T>& m) const{

      Matrix4x4<T> result;
      T sum;
      for(unsigned int i = 0; i < 4; ++i)
        for(unsigned int j = 0; j < 4; ++j){
          sum = T(0.0);

          for(unsigned int k = 0; k < 4; ++k){
            sum += (*this)[i][k]*m[k][j];
          }
          result[i][j] = sum;
        }

      return result;
    }

  template<typename T>
    Matrix4x4<T> Matrix4x4<T>::operator*(const T& d) const{
      Matrix4x4<T> result = *this;
      result *= d;
      return result;
    }


  template<typename T>
    Matrix4x4<T> Matrix4x4<T>::operator-(void) const{
      Matrix4x4<T> result = *this;
      std::transform(this->begin_, this->end_, result.begin_, std::negate<T>());
      return result;
    }


  template<typename T>
    Matrix4x4<T>& Matrix4x4<T>::operator+=(const Matrix4x4<T>& m){
      std::transform(this->begin_, this->end_, m.begin_, this->begin_, std::plus<T>());
      return *this;
    }


  template<typename T>
    Matrix4x4<T>& Matrix4x4<T>::operator-=(const Matrix4x4<T>& m){
      std::transform(this->begin_, this->end_, m.begin_, this->begin_, std::minus<T>());
      return *this;
    }

  template<typename T>
    Matrix4x4<T>& Matrix4x4<T>::operator/=(const T& d){
      if(d==0)
        throw MathException("Value is null");
      //std::transform(this->begin_, this->end_, this->begin_, std::bind2nd(std::divides<T>(),d));
      auto lambda = std::bind(std::divides<T>(), std::placeholders::_1, d);
      std::transform(this->begin_, this->end_, this->begin_, lambda);
      return *this;
    }


  template<typename T>
    Matrix4x4<T>& Matrix4x4<T>::operator*=(const Matrix4x4<T>& m){
      Matrix4x4<T> tmp = *this;
      T sum;
      for(unsigned int i = 0; i < 4; ++i)
        for(unsigned int j = 0; j < 4; ++j){
          sum = T(0.0);

          for(unsigned int k = 0; k < 4; ++k){
            sum += tmp[i][k]*m[k][j];
          }
          this->accessor[i][j] = sum;
        }
      return *this;
    }

  template<typename T>
    Matrix4x4<T>& Matrix4x4<T>::operator*=(const T& d){
      //std::transform(this->begin_, this->end_, this->begin_, std::bind2nd(std::multiplies<T>(),d));
        auto lambda = std::bind(std::multiplies<T>(), std::placeholders::_1, d);
        std::transform(this->begin_, this->end_, this->begin_, lambda);
      return *this;
    }


  template<typename T>
    Matrix4x4<T>& Matrix4x4<T>::power(const unsigned int& p){
      if(p == 0){
        setIdentity();
        return *this;
      }
      Matrix4x4<T> tmp(*this);

      for(unsigned int i=1; i < p; ++i)
        *this *= tmp;

      return *this;
    }

  template<typename T>
    Vector4<T> Matrix4x4<T>::getDiagonal(void) const{
      Vector4<T> v;
      for(unsigned int i = 0; i < 4; ++i)
        v[i] = this->accessor[i][i];
      return v;
    }

  template<typename T>
    Vector4<T> Matrix4x4<T>::getRow(const unsigned int& row) const{
      if(row >= this->rowsMatrix)
        throw MathException("Out of bounds");
      Vector4<T> v;
      for(unsigned int i = 0; i < this->columnsMatrix; ++i)
        v[i] = this->accessor[row][i];
      return v;
    }

  template<typename T>
    Vector4<T> Matrix4x4<T>::getColumn(const unsigned int& column) const{
      if(column >= this->columnsMatrix)
        throw MathException("Out of bounds");
      Vector4<T> v;
      for(unsigned int i = 0; i < this->columnsMatrix; ++i)
        v[i] = this->accessor[i][column];
      return v;
    }

  template<typename T>
    T Matrix4x4<T>::trace(void) const{
      return (this->accessor[0][0] + this->accessor[1][1] + this->accessor[2][2] + this->accessor[3][3]);
    }


  template<class U> Matrix4x4<U> operator* (const U& d,
      const Matrix4x4<U>& m){
    return m*d;
  }


  /*
   * Type definition
   */

  typedef Matrix4x4<float>  Matrix4x4f;
  typedef Matrix4x4<double> Matrix4x4d;
  typedef Matrix4x4<int>    Matrix4x4i;

}


#endif
