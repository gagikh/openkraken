/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * Matrix.hpp created in 09 2008.
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
#ifndef __OPENKN_MATH__MATRIX_HPP__
#define __OPENKN_MATH__MATRIX_HPP__

/*
 * External Includes
 */
#include <iostream>
#include <algorithm>
#include <functional>
//#include <cmath>

/*
 * Internal Includes
 */
#include "Vector.hpp"
#include "MathException.hpp"

/*
 * Namespace
 */
namespace kn {

  /*
   * Class definition
   */

  template<class T> class Matrix {

    protected:

      size_t rowsMatrix;

      size_t columnsMatrix;

      T* data;

      T* begin_;

      T* end_;

      T** accessor;

    public:

      /*
       * Constructor
       */

      Matrix();

      Matrix(const Matrix<T>& m);

      Matrix(const Matrix<T>* m);

      Matrix(const size_t& n,
          const size_t& m);

      explicit Matrix(const size_t& n);

      Matrix(const size_t& n,
          const size_t& m,
          const T& d);

      Matrix(const size_t& n,
          const size_t& m,
          T* a);


      Matrix(const size_t& n, T* a);

      Matrix(const size_t& n,
          const size_t& m,
          const Vector<T>& v,
          const bool& setasrows = true);

      Matrix(const size_t& n,
          const Vector<T>* v,
          const bool& setasrows = true);

      /*
       * Destructor
       */

      virtual ~Matrix(void){
        if(data!=0){
          delete[] accessor;
          delete[] data;
        }
      }

    protected:

      void allocate(void);

    public:

       inline virtual T  * begin() const {return begin_;}

       inline virtual T  * end() const {return end_;}


      T& at(const unsigned int& i, const unsigned int& j);

      const T& at(const unsigned int& i, const unsigned int& j) const;

      /*
       * Operators redefinition
       */

      bool operator==(const Matrix<T>& m) const;

      inline bool operator!=(const Matrix<T>& m) const
      {
        return !(*this == m);
      }

      Matrix<T>& operator=(const Matrix<T>& m);

      Matrix<T>& operator=(const T & value);

      Matrix<T> operator+(const Matrix<T>& m) const;

      Matrix<T> operator-(const Matrix<T>& m) const;

      Matrix<T> operator/(const T& d) const;

      Matrix<T> operator*(const Matrix<T>& m) const;

      Vector<T> operator*(const Vector<T>& v) const;

      Matrix<T> operator*(const T& d) const;

      void times(const Matrix<T>& m1,const Matrix<T>& m2);

      Matrix<T> operator-(void) const;

      Matrix<T>& operator+=(const Matrix<T>& m);

      Matrix<T>& operator+=(const T& d);

      Matrix<T>& operator-=(const Matrix<T>& m);

      Matrix<T>& operator-=(const T& d);

      Matrix<T>& operator/=(const T& d);

      Matrix<T>& operator*=(const Matrix<T>& m);

      Matrix<T>& operator*=(const T& d);

      T* operator[](const unsigned int& i) const;

      T& operator()(const unsigned int& i,
          const unsigned int& j);

      const T& operator()(const unsigned int& i,
          const unsigned int& j) const;

      /*
       * Accessors
       */

      inline size_t rows(void) const{
        return rowsMatrix;
      }

      inline size_t columns(void) const{
        return columnsMatrix;
      }

      double getNorm(void) const;

      const T* getMatrixArray(void) const{
        return data;
      }


      Matrix<T>& transpose(void);

      inline Matrix<T> getTranspose(void) const
      {
        return Matrix<T>(*this).transpose();
      }

      Matrix<T>& power(const unsigned int& p);

      void setIdentity(void);

      inline void setZero(void)
      {
        std::fill(begin_,end_,T(0.0));
      }

      inline void fill(const T &d){
        std::fill(begin_,end_,d);
      }

      void roundZero(const double& d = 1.0e-14);

      void setSubMatrix(const unsigned int& row,
          const unsigned int& column,
          const Matrix<T>& m);

      Matrix<T> getSubMatrix(const unsigned int& firstrow,const unsigned int& row,
        const unsigned int& firstcolumn,const unsigned int& column);
      void setRow(const unsigned int& row,
          const Vector<T>& v);

      void setRow(const unsigned int& row,
          T* a);

      void setRow(const unsigned int& row,
          const T& d);

      void setColumn(const unsigned int& column,
          const Vector<T>& v);

      void setColumn(const unsigned int& column,
          const T* a);

      void setColumn(const unsigned int& column,
          const T& d);

      Vector<T> getRow(const unsigned int& row) const;

      Vector<T> getColumn(const unsigned int& column) const;

      void swapRows(const unsigned int& row1,
          const unsigned int& row2);

      void swapColumns(const unsigned int& column1,
          const unsigned int& column2);

      void setDiagonal(const Vector<T>& v);

      void setDiagonal(const T* a);

      void setDiagonal(const T& d);

      Vector<T> getDiagonal(void) const;

      inline bool isSquare(void) const{
        return rowsMatrix == columnsMatrix;
      }

      T trace(void) const;

      void cross3x3(const kn::Vector<T> v);

  };

  template<typename T>
    Matrix<T>:: Matrix()
    : rowsMatrix(0), columnsMatrix(0){
      data     = NULL;
      begin_   = NULL;
      end_     = NULL;
      accessor = NULL;
    }

  template<typename T>
    Matrix<T>:: Matrix(const Matrix<T>& m)
    : rowsMatrix(m.rows()), columnsMatrix(m.columns()){
      allocate();
      std::copy(m.begin_, m.end_, begin_);
    }


  template<typename T>
    Matrix<T>::Matrix(const Matrix<T>* m)
    : rowsMatrix(m->rows()), columnsMatrix(m->columns()){
      if(m==0)
        throw MathException("Pointer is null");
      allocate();
      std::copy(m->begin_, m->end_, begin_);
    }


  template<typename T>
    Matrix<T>::Matrix(const size_t& n,
        const size_t& m)
    :rowsMatrix(n), columnsMatrix(m){
      if(m==0 || n==0)
        throw MathException("Matrix size is null");
      allocate();
    }


  template<typename T>
    Matrix<T>::Matrix(const size_t& n)
    :rowsMatrix(n), columnsMatrix(n){
      if(n==0)
        throw MathException("Matrix size is null");
      allocate();
    }


  template<typename T>
    Matrix<T>::Matrix(const size_t& n,
        const size_t& m,
        const T& d)
    :rowsMatrix(n), columnsMatrix(m){
      if(m==0 || n==0)
        throw MathException("Matrix size is null");
      allocate();
      std::fill(begin_, end_, d);
    }


  template<typename T>
    Matrix<T>::Matrix(const size_t& n,
        const size_t& m,
        T* a)
    :rowsMatrix(n), columnsMatrix(m){
      if(m==0 || n==0)
        throw MathException("Matrix size is null");
      if(a==0)
        throw MathException("Pointer is null");
      allocate();
      std::copy(a, a + n*m, this->begin_);
    }



  template<typename T>
    Matrix<T>::Matrix(const size_t& n, T* a)
    :rowsMatrix(n), columnsMatrix(n){
      if(n==0)
        throw MathException("Matrix size is null");
      if(a==0)
        throw MathException("Pointer is null");
      allocate();
      std::copy(a, a + n*n, this->begin_);
    }


  template<typename T>
    Matrix<T>::Matrix(const size_t& n,
        const size_t& m,
        const Vector<T>& v,
        const bool& setasrows)
    :rowsMatrix(n), columnsMatrix(m){
      if(m==0 || n==0)
        throw MathException("Matrix size is null");
      if(v.size()!=n*m)
        throw MathException("Vector size is different from Matrix size");

      allocate();

      if(setasrows)
        for(unsigned int i = 0; i < rowsMatrix*columnsMatrix; ++i)
          data[i] = v[i];
      else
        for(unsigned int j = 0; j < columnsMatrix; ++j)
          for(unsigned int i = 0; i < rowsMatrix; ++i)
            accessor[i][j] = v[j*columnsMatrix+i];
    }


  template<typename T> Matrix<T>::Matrix(const size_t& n,
      const Vector<T>* v,
      const bool& setasrows){

    if(n==0)
      throw MathException("Matrix size is null");
    if(v==0)
      throw MathException("Pointer is null");
    unsigned int checksize = 0;
    for(unsigned int i = 0; i < n; ++i){
      if(v[i].size()==0)
        throw MathException("One vector size is null");
      if(i==0)
        checksize = v[i].size();
      else
        if(v[i].size() != checksize)
          throw MathException("Vectors size are different");
    }

    if(setasrows){
      rowsMatrix = n;
      columnsMatrix = v[0].size();
      allocate();
      for(unsigned int i = 0; i < n; ++i)
        for(unsigned int j = 0; j < (v[i]).size(); ++j)
          accessor[i][j] = v[i][j];
    }else{
      rowsMatrix = v[0].size();
      columnsMatrix = n;
      allocate();
      for(unsigned int j = 0; j < n; ++j)
        for(unsigned int i = 0; i < (v[i]).size(); ++i)
          accessor[i][j] = v[j][i];
    }
  }



  template<typename T>
    void Matrix<T>::allocate(void){
      this->data = this->begin_ = this->end_ = 0;
      this->accessor = 0;
      this->data = new T[rowsMatrix*columnsMatrix];
      this->begin_ = this->data;
      this->end_ = this->begin_ + this->rowsMatrix*this->columnsMatrix;
      this->accessor = new T*[this->rowsMatrix];
      for(unsigned int i = 0; i < this->rowsMatrix; ++i)
        this->accessor[i] = this->begin_ + i*this->columnsMatrix;
    }

  template<typename T>
    T& Matrix<T>::at(const unsigned int& i, const unsigned int& j) {
      if(i>=rowsMatrix || j>=columnsMatrix)
        throw MathException("Out of bounds");

      return accessor[i][j];
    }

  template<typename T>
    const T& Matrix<T>::at(const unsigned int& i, const unsigned int& j) const {
      if(i>=rowsMatrix || j>=columnsMatrix)
        throw MathException("Out of bounds");

      return accessor[i][j];
    }


  template<typename T>
    bool Matrix<T>::operator==(const Matrix<T>& m) const{
      if(m.rowsMatrix != rowsMatrix || m.columnsMatrix != columnsMatrix) return false;
      return std::equal(begin_, end_, m.begin_);
    }


  template<typename T>
    Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m){
      if(&m == this) return *this;

      if(rowsMatrix == 0 && columnsMatrix == 0){
        rowsMatrix = m.rowsMatrix;
        columnsMatrix = m.columnsMatrix;
        allocate();
      }

      if(m.rowsMatrix == rowsMatrix && m.columnsMatrix == columnsMatrix)
        std::copy(m.begin_,m.end_,begin_);
      else
        throw  MathException("Matrices'size are different");

      return *this;
    }

  template<typename T>
    Matrix<T>& Matrix<T>::operator=(const T & value) {
      std::fill(begin_,end_,T(value));

      return *this;
    }


  template<typename T>
    Matrix<T> Matrix<T>::operator+(const Matrix<T>& m) const{
      Matrix<T> result = *this;
      result += m;
      return result;
    }


  template<typename T>
    Matrix<T> Matrix<T>::operator-(const Matrix<T>& m) const{
      Matrix<T> result = *this;
      result -= m;
      return result;
    }


  template<typename T>
    Matrix<T> Matrix<T>::operator/(const T& d) const{
      Matrix<T> result = *this;
      result /= d;
      return result;
    }


  template<typename T>
    Matrix<T> Matrix<T>::operator*(const Matrix<T>& m) const{

      if(m.rowsMatrix != columnsMatrix)
        throw MathException("Matrices'size is incorrect");

      Matrix<T> result(rowsMatrix,m.columnsMatrix);
      result.times(this,m);
      return result;
    }

    template<typename T>
    void Matrix<T>::times(const Matrix<T>& m1,const Matrix<T>& m2){

      if(m1.columnsMatrix != m2.rowsMatrix)
        throw MathException("Matrices'size is incorrect");

      if(m1.rowsMatrix != rowsMatrix)
        throw MathException("Matrices' rows size is incorrect");

      if(m2.columnsMatrix != columnsMatrix)
        throw MathException("Matrices' columns size is incorrect");

      T sum;
      for(unsigned int i = 0; i < rowsMatrix; ++i)
        for(unsigned int j = 0; j < columnsMatrix; ++j){
          sum = T(0.0);
          for(unsigned int k = 0; k < m1.columnsMatrix; ++k){
            sum += m1[i][k]*m2[k][j];
          }
          accessor[i][j] = sum;
        }
    }



  template<typename T>
    Vector<T> Matrix<T>::operator*(const Vector<T>& v) const{
      if(v.size() != columnsMatrix)
        throw MathException("Matrice rows and vector size are different");

      Vector<T> vtmp(rowsMatrix);
      T sum;
      for(unsigned int i = 0; i < rowsMatrix; ++i){
        sum = T(0);
        for(unsigned int j = 0; j < columnsMatrix; ++j){
          sum += accessor[i][j]*v[j];
        }
        vtmp[i] = sum;
      }
      return vtmp;
    }


  template<typename T>
    Matrix<T> Matrix<T>::operator*(const T& d) const{
      Matrix<T> result = *this;
      result *= d;
      return result;
    }


  template<typename T>
    Matrix<T> Matrix<T>::operator-(void) const{
      Matrix<T> result = *this;
      std::transform(begin_, end_, result.begin_, std::negate<T>());
      return result;
    }


  template<typename T>
    Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m){
      if(m.rowsMatrix != rowsMatrix || m.columnsMatrix != columnsMatrix)
        throw MathException("Matrices'size is incorrect");
      std::transform(begin_, end_, m.begin_, begin_, std::plus<T>());
      return *this;
    }


  template<typename T>
    Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m){
      if(m.rowsMatrix != rowsMatrix || m.columnsMatrix != columnsMatrix)
        throw MathException("Matrices'size is incorrect");
      std::transform(begin_, end_, m.begin_, begin_, std::minus<T>());
      return *this;
    }

  template<typename T>
    Matrix<T>& Matrix<T>::operator/=(const T& d){
      if(d==0)
        throw MathException("Value is null");
      //std::transform(begin_, end_, begin_, std::bind2nd(std::divides<T>(),d));
      auto lambda = std::bind(std::divides<T>(), std::placeholders::_1, d);
      std::transform(begin_, end_, begin_, lambda);
      return *this;
    }

  template<typename T>
    Matrix<T>& Matrix<T>::operator+=(const T& d){
      //std::transform(begin_, end_, begin_, std::bind2nd(std::plus<T>(),d));
      auto lambda = std::bind(std::plus<T>(), std::placeholders::_1, d);
      std::transform(begin_, end_, begin_, lambda);
      return *this;
    }

  template<typename T>
    Matrix<T>& Matrix<T>::operator-=(const T& d){
      //std::transform(begin_, end_, begin_, std::bind2nd(std::minus<T>(),d));
      auto lambda = std::bind(std::minus<T>(), std::placeholders::_1, d);
      std::transform(begin_, end_, begin_, lambda);
      return *this;
    }

  template<typename T>
    Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& m){

      if(rowsMatrix != columnsMatrix)
        throw MathException("Matrices'size is incorrect");

      if(m.rowsMatrix != m.columnsMatrix)
        throw MathException("Matrices'size is incorrect");

      if(m.rowsMatrix != columnsMatrix)
        throw MathException("Matrices'size is incorrect");

      Matrix<T> tmp = *this;
      T sum;
      for(unsigned int i = 0; i < rowsMatrix; ++i)
        for(unsigned int j = 0; j < m.columnsMatrix; ++j){
          sum = T(0.0);

          for(unsigned int k = 0; k < columnsMatrix; ++k){
            sum += tmp[i][k]*m[k][j];
          }
          accessor[i][j] = sum;
        }
      return *this;
    }

  template<typename T>
    Matrix<T>& Matrix<T>::operator*=(const T& d){
      //std::transform(begin_, end_, begin_, std::bind2nd(std::multiplies<T>(),d));
      auto lambda = std::bind(std::multiplies<T>(), std::placeholders::_1, d);
      std::transform(begin_, end_, begin_, lambda);
      return *this;
    }


  template<typename T>
    T* Matrix<T>::operator[](const unsigned int& i) const{
      return accessor[i];
    }

  template<typename T>
    T& Matrix<T>::operator()(const unsigned int& i, const unsigned int& j) {
      return accessor[i][j];
    }

  template<typename T>
    const T& Matrix<T>::operator()(const unsigned int& i, const unsigned int& j) const{
      return accessor[i][j];
    }




  template<typename T>
    double Matrix<T>::getNorm(void) const{
      double norm = 0.0;
      for(unsigned int i = 0; i < rowsMatrix*columnsMatrix; ++i)
        norm += data[i]*data[i];
      return sqrt(norm);
    }



  template<typename T>
    Matrix<T>& Matrix<T>::transpose(void){
      Matrix<T> tmp(*this);
      std::swap(rowsMatrix,columnsMatrix);

      if(data!=0){
        delete[] accessor;
        delete[] data;
      }
      allocate();
      for(unsigned int i = 0; i < rowsMatrix; ++i)
        for(unsigned int j = 0; j < columnsMatrix; ++j)
          accessor[i][j] = tmp[j][i];

      return *this;
    }


  template<typename T>
    Matrix<T>& Matrix<T>::power(const unsigned int& p){
      if(rowsMatrix != columnsMatrix)
        throw MathException("Matrix is not square");

      if(p == 0){
        setIdentity();
        return *this;
      }
      Matrix<T> tmp(*this);

      for(unsigned int i=1; i < p; ++i)
        *this *= tmp;

      return *this;
    }

  template<typename T>
    void Matrix<T>::setIdentity(void){
      setZero();
      setDiagonal(T(1.0));
    }


  template<typename T>
    void Matrix<T>::roundZero(const double& d){
      T *iter = begin_;

      while(iter != end_){
        if(std::abs((*iter)) < d) (*iter) = T(0.0);
        ++iter;
      }
    }

  template<typename T>
    void Matrix<T>::setSubMatrix(const unsigned int& row,
        const unsigned int& column,
        const Matrix<T>& m){

      if(row + m.rowsMatrix > rowsMatrix || column + m.columnsMatrix > columnsMatrix)
        throw MathException("Input matrix size is incorrect");

      for(unsigned int i = 0; i < m.rowsMatrix; ++i)
        std::copy(m.accessor[i], m.accessor[i]+m.columnsMatrix,accessor[i+row]+column);

    }

  template<typename T>
    Matrix<T> Matrix<T>::getSubMatrix(const unsigned int& firstrow,const unsigned int& row,
        const unsigned int& firstcolumn,const unsigned int& column){

      if(firstrow + row > rowsMatrix || firstcolumn + column > columnsMatrix)
        throw MathException("Input parameter row/column size is incorrect");

      Matrix<T> tmp(row,column);

      for(unsigned int i = 0; i < row; ++i)
        std::copy(accessor[firstrow+i]+firstcolumn, accessor[firstrow+i]+firstcolumn+column,tmp.accessor[i]);

      return tmp;
    }

  template<typename T>
    void Matrix<T>::setRow(const unsigned int& row,
        const Vector<T>& v){
      if(row >= rowsMatrix)
        throw MathException("Out of bounds");
      if (v.size() > columnsMatrix)
        throw MathException("Input vector size is incorrect");

      for(unsigned int i = 0; i < columnsMatrix; ++i)
        accessor[row][i] = v[i];
    }

  template<typename T>
    void Matrix<T>::setRow(const unsigned int& row,
        T* a){
      if(row >= rowsMatrix)
        throw MathException("Out of bounds");
      std::copy(a, a+columnsMatrix,accessor[row]);
    }

  template<typename T>
    void Matrix<T>::setRow(const unsigned int& row,
        const T& d){
      if(row >= rowsMatrix)
        throw MathException("Out of bounds");
      std::fill(accessor[row], accessor[row]+columnsMatrix,d);
    }


  template<typename T>
    void Matrix<T>::setColumn(const unsigned int& column,
        const Vector<T>& v){
      if(column >= columnsMatrix)
        throw MathException("Out of bounds");
      if (v.size() > rowsMatrix)
        throw MathException("Input vector size is incorrect");

      for(unsigned int i = 0; i < rowsMatrix; ++i)
        accessor[i][column] = v[i];
    }

  template<typename T>
    void Matrix<T>::setColumn(const unsigned int& column,
        const T* a){
      if(column >= columnsMatrix)
        throw MathException("Out of bounds");

      for(unsigned int i = 0; i < rowsMatrix; ++i)
        accessor[i][column] = a[i];
    }


  template<typename T>
    void Matrix<T>::setColumn(const unsigned int& column,
        const T& d){
      if(column >= columnsMatrix)
        throw MathException("Out of bounds");

      for(unsigned int i = 0; i < rowsMatrix; ++i)
        accessor[i][column] = d;
    }

  template<typename T>
    Vector<T> Matrix<T>::getRow(const unsigned int& row) const{
      if(row >= rowsMatrix)
        throw MathException("Out of bounds");
      Vector<T> v(columnsMatrix);
      for(unsigned int i = 0; i < columnsMatrix; ++i)
        v[i] = accessor[row][i];
      return v;
    }

  template<typename T>
    Vector<T> Matrix<T>::getColumn(const unsigned int& column) const{
      if(column >= columnsMatrix)
        throw MathException("Out of bounds");
      Vector<T> v(rowsMatrix);
      for(unsigned int i = 0; i < rowsMatrix; ++i)
        v[i] = accessor[i][column];
      return v;
    }


  template<typename T>
    void Matrix<T>::swapRows(const unsigned int& row1,
        const unsigned int& row2){
      if(row1 >= rowsMatrix || row2 >= rowsMatrix)
        throw MathException("Out of bounds");

      std::swap_ranges(accessor[row1],accessor[row1]+columnsMatrix,accessor[row2]);
    }


  template<typename T>
    void Matrix<T>::swapColumns(const unsigned int& column1,
        const unsigned int& column2){
      if(column1 >= columnsMatrix || column2 >= columnsMatrix)
        throw MathException("Out of bounds");

      for(unsigned int i = 0; i < rowsMatrix; ++i)
        std::swap(accessor[i][column1],accessor[i][column2]);
    }


  template<typename T>
    void Matrix<T>::setDiagonal(const Vector<T>& v){
      unsigned int min = std::min(columnsMatrix,rowsMatrix);

      if(v.size() != min)
        throw MathException("Vector size is incorrect");

      for(unsigned int i = 0; i < min; ++i)
        accessor[i][i] = v[i];
    }


  template<typename T>
    void Matrix<T>::setDiagonal(const T* a){
      unsigned int min = std::min(columnsMatrix,rowsMatrix);

      for(unsigned int i = 0; i < min; ++i)
        accessor[i][i] = a[i];
    }

  template<typename T>
    void Matrix<T>::setDiagonal(const T& d)
    {
      unsigned int min = std::min(columnsMatrix,rowsMatrix);

      for(unsigned int i = 0; i < min; ++i)
        accessor[i][i] = d;
    }


  template<typename T>
    Vector<T> Matrix<T>::getDiagonal(void) const{
      unsigned int min = std::min(columnsMatrix,rowsMatrix);
      Vector<T> v(min);
      for(unsigned int i = 0; i < min; ++i)
        v[i] = accessor[i][i];
      return v;
    }

  template<typename T>
    T Matrix<T>::trace(void) const{
      T sum = 0;
      unsigned int min = std::min(columnsMatrix,rowsMatrix);
      for(unsigned int i = 0; i < min; ++i)
        sum += accessor[i][i];
      return sum;
    }

  template<typename T>
    void Matrix<T>::cross3x3(const kn::Vector<T> v){
      if(v.size() != 3)
        throw MathException("argument shoud be a 3-vector");

      if(rowsMatrix != 3 && columnsMatrix != 3)
        throw MathException("valid only with 3x3 matrices");

      accessor[0][0] =   0.0;  accessor[0][1] = -v[2];  accessor[0][2] =  v[1];
      accessor[1][0] =  v[2];  accessor[1][1] =   0.0;  accessor[1][2] = -v[0];
      accessor[2][0] = -v[1];  accessor[2][1] =  v[0];  accessor[2][2] =   0.0;
    }


  template<class U> std::ostream& operator<< (std::ostream& stream,
      const Matrix<U>& m){
    for(unsigned int i=0; i<m.rows(); ++i){
      for (unsigned int j=0; j<m.columns(); ++j)
        stream<<m[i][j] << " " ;
      stream<<std::endl;
    }
    return stream;
  }

  template<class U> Vector<U> operator* (const Vector<U>& v,
      const Matrix<U>& m){
    if(m.rows() != v.size())
      throw MathException("Matrix or vector size is incorrect");

    Vector<U> result(m.columns());
    U sum;
    for(unsigned int i=0; i<m.columns(); ++i){
      sum = U(0.0);
      for(unsigned int j=0; j<m.rows(); ++j)
        sum += m[j][i] * v[j];
      result[i] = sum;
    }
    return result;
  }

  template<class U> Matrix<U> operator* (const U& d,
      const Matrix<U>& m){
    return m*d;
  }

  /*
   * Type definition
   */

  typedef Matrix<float>  Matrixf;
  typedef Matrix<double> Matrixd;
  typedef Matrix<int>      Matrixi;


  /*
   * End of Namespace
   */
}


/*
 * End of Anti-doublon
 */
#endif
