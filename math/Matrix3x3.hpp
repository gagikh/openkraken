/*************************************************************************** \
 * Copyright (C) by University Paris-Est - MISS team
 * Matrix3x3.hpp created in 09 2008.
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
\***************************************************************************/

/*
 * Anti-doublon
 */
#ifndef __OPENKN_MATH__MATRIX3X3_HPP__
#define __OPENKN_MATH__MATRIX3X3_HPP__

/*
 * Internal Includes
 */
#include "Matrix.hpp"
#include "Vector3.hpp"

namespace kn
{

        template<class T> class Matrix3x3 : public Matrix <T>
        {

        protected:

        public:

            Matrix3x3(void)
                : Matrix<T>(3)
                {
                    setIdentity();
                }

            Matrix3x3(const Matrix3x3<T>& m)
                : Matrix<T>(3)
                {
                    std::copy(m.begin_, m.end_, this->begin_);
                }

                        Matrix3x3(const Matrix<T>& m);

            Matrix3x3(const Matrix3x3<T>* m);

            Matrix3x3(const T& d)
                : Matrix<T>(3)
                {
                    std::fill(this->begin_, this->end_, d);
                }


            Matrix3x3(const T* d);

            Matrix3x3(const Vector<T>& v,
                      const bool& setasrows = true);

            Matrix3x3(const Vector3<T>* v,
                      const bool& setasrows = true);

            ~Matrix3x3(void)
                {
                }


        public:

            using Matrix<T>::operator*;


            Vector<T> operator*(const Vector<T>& v) const;


            inline size_t rows(void) const
                {
                    return 3;
                }

            inline size_t columns(void) const
                {
                    return 3;
                }

            Matrix3x3<T>& transpose(void);

            inline Matrix3x3<T> getTranspose(void) const
                {
                    return Matrix3x3<T>(*this).transpose();
                }

            void setIdentity(void);

            inline bool isSquare(void) const
                {
                    return true;
                }

           void cross3x3(const kn::Vector3<T> v);

           Matrix3x3<T>& operator=(const Matrix3x3<T>& m);

           Matrix3x3<T> operator+(const Matrix3x3<T>& m) const;

           Matrix3x3<T> operator-(const Matrix3x3<T>& m) const;

           Matrix3x3<T> operator/(const T& d) const;

           Matrix3x3<T> operator*(const Matrix3x3<T>& m) const;

           Matrix3x3<T> operator*(const T& d) const;

           Matrix3x3<T> operator-(void) const;

           Matrix3x3<T>& operator+=(const Matrix3x3<T>& m);

           Matrix3x3<T>& operator-=(const Matrix3x3<T>& m);

           Matrix3x3<T>& operator/=(const T& d);

           Matrix3x3<T>& operator*=(const Matrix3x3<T>& m);

           Matrix3x3<T>& operator*=(const T& d);

           Vector3<T> getRow(const unsigned int& row) const;

           Vector3<T> getColumn(const unsigned int& column) const;

           Vector3<T> getDiagonal(void) const;

           Matrix3x3<T>& power(const unsigned int& p);

           T trace(void) const;

        };



        template<class T>
          Matrix3x3<T>::Matrix3x3(const Matrix3x3<T>* m)
          : Matrix<T>(3)
          {
            if(m==0)
              throw MathException("Pointer is null");
            std::copy(m->begin_, m->end_, this->begin_);
          }


        template<class T>
          Matrix3x3<T>::Matrix3x3(const T* a)
          : Matrix<T>(3)
          {
            if(a==0)
              throw MathException("Pointer is null");
            std::copy(a, a + 9, this->begin_);
          }


        template<class T>
          Matrix3x3<T>::Matrix3x3(const Vector<T>& v,
              const bool& setasrows)
          : Matrix<T>(3)
          {
            if(v.size()!=9)
              throw MathException("Vector size is different from Matrix size");

            if(setasrows)
              for(unsigned int i = 0; i < 9; ++i)
                (this->data)[i] = v[i];
            else
              for(unsigned int j = 0; j < 3; ++j)
                for(unsigned int i = 0; i < 3; ++i)
                  (this->accessor)[i][j] = v[j*3+i];
          }


        template<class T>
          Matrix3x3<T>::Matrix3x3(const Vector3<T>* v,
              const bool& setasrows)
          : Matrix<T>(3)
          {

            if(v==0)
              throw MathException("Pointer is null");

            if(setasrows)
              for(unsigned int i = 0; i < 3; ++i)
                for(unsigned int j = 0; j < 3; ++j)
                  (this->accessor)[i][j] = v[i][j];
            else
              for(unsigned int j = 0; j < 3; ++j)
                for(unsigned int i = 0; i < 3; ++i)
                  (this->accessor)[i][j] = v[j][i];

          }

        template<class T>
          Matrix3x3<T>::Matrix3x3(const Matrix<T>& m)
          : Matrix<T>(3) {
            if ((m.rows() != 3) || (m.columns() != 3))
              throw MathException("Matrix3x3 Constructor","Matrix sizes are different");
            std::copy(m.getMatrixArray(), m.getMatrixArray()+m.rows()*m.columns(), this->begin_);
          }

        template<class T>
          Vector<T> Matrix3x3<T>::operator*(const Vector<T>& v) const
          {
            if(v.size() != 3)
              throw MathException("Matrice rows and vector size are different");

            Vector<T> vtmp(3);
            vtmp[0] = v[0]*(this->accessor)[0][0] + v[1]*(this->accessor)[0][1] + v[2]*(this->accessor)[0][2];
            vtmp[1] = v[0]*(this->accessor)[1][0] + v[1]*(this->accessor)[1][1] + v[2]*(this->accessor)[1][2];
            vtmp[2] = v[0]*(this->accessor)[2][0] + v[1]*(this->accessor)[2][1] + v[2]*(this->accessor)[2][2];
            return vtmp;
          }


        template<class T>
          Matrix3x3<T>& Matrix3x3<T>::transpose(void)
          {

            for(unsigned int i = 0; i < 3; ++i)
              for(unsigned int j = 0; j < i; ++j)
                std::swap((this->accessor)[i][j],(this->accessor)[j][i]);

            return *this;
          }


        template<class T>
          void Matrix3x3<T>::setIdentity(void)
          {
            this->setZero();
            (this->accessor)[0][0] = (this->accessor)[1][1] = (this->accessor)[2][2] = T(1.0);
          }




        template<class U>
          Vector<U> operator* (const Vector3<U>& v,
              const Matrix3x3<U>& m)
          {

            Vector<U> result(3);
            result[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
            result[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
            result[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
            return result;
          }


        template<typename T>
          Matrix3x3<T>& Matrix3x3<T>::operator=(const Matrix3x3<T>& m){
            if(&m == this) return *this;
            std::copy(m.begin_,m.end_,this->begin_);

            return *this;
          }


        template<typename T>
          Matrix3x3<T> Matrix3x3<T>::operator+(const Matrix3x3<T>& m) const{
            Matrix3x3<T> result = *this;
            result += m;
            return result;
          }


        template<typename T>
          Matrix3x3<T> Matrix3x3<T>::operator-(const Matrix3x3<T>& m) const{
            Matrix3x3<T> result = *this;
            result -= m;
            return result;
          }


        template<typename T>
          Matrix3x3<T> Matrix3x3<T>::operator/(const T& d) const{
            Matrix3x3<T> result = *this;
            result /= d;
            return result;
          }


        template<typename T>
          Matrix3x3<T> Matrix3x3<T>::operator*(const Matrix3x3<T>& m) const{

            Matrix3x3<T> result;
            T sum;
            for(unsigned int i = 0; i < 3; ++i)
              for(unsigned int j = 0; j < 3; ++j){
                sum = T(0.0);

                for(unsigned int k = 0; k < 3; ++k){
                  sum += (*this)[i][k]*m[k][j];
                }
                result[i][j] = sum;
              }

            return result;
          }

        template<typename T>
          Matrix3x3<T> Matrix3x3<T>::operator*(const T& d) const{
            Matrix3x3<T> result = *this;
            result *= d;
            return result;
          }


        template<typename T>
          Matrix3x3<T> Matrix3x3<T>::operator-(void) const{
            Matrix3x3<T> result = *this;
            std::transform(this->begin_, this->end_, result.begin_, std::negate<T>());
            return result;
          }


        template<typename T>
          Matrix3x3<T>& Matrix3x3<T>::operator+=(const Matrix3x3<T>& m){
            std::transform(this->begin_, this->end_, m.begin_, this->begin_, std::plus<T>());
            return *this;
          }


        template<typename T>
          Matrix3x3<T>& Matrix3x3<T>::operator-=(const Matrix3x3<T>& m){
            std::transform(this->begin_, this->end_, m.begin_, this->begin_, std::minus<T>());
            return *this;
          }

        template<typename T>
          Matrix3x3<T>& Matrix3x3<T>::operator/=(const T& d){
            if(d==0)
              throw MathException("Value is null");
            //std::transform(this->begin_, this->end_, this->begin_, std::bind2nd(std::divides<T>(),d));
            auto lambda = std::bind(std::divides<T>(), std::placeholders::_1, d);
            std::transform(this->begin_, this->end_, this->begin_, lambda);
            return *this;
          }


        template<typename T>
          Matrix3x3<T>& Matrix3x3<T>::operator*=(const Matrix3x3<T>& m){

            Matrix3x3<T> tmp = *this;
            T sum;
            for(unsigned int i = 0; i < 3; ++i)
              for(unsigned int j = 0; j < 3; ++j){
                sum = T(0.0);

                for(unsigned int k = 0; k < 3; ++k){
                  sum += tmp[i][k]*m[k][j];
                }
                this->accessor[i][j] = sum;
              }
            return *this;
          }

        template<typename T>
          Matrix3x3<T>& Matrix3x3<T>::operator*=(const T& d){
            //std::transform(this->begin_, this->end_, this->begin_, std::bind2nd(std::multiplies<T>(),d));
            auto lambda = std::bind(std::multiplies<T>(), std::placeholders::_1, d);
            std::transform(this->begin_, this->end_, this->begin_, lambda);
            return *this;
          }


        template<typename T>
          Matrix3x3<T>& Matrix3x3<T>::power(const unsigned int& p){

            if(p == 0){
              setIdentity();
              return *this;
            }
            Matrix3x3<T> tmp(*this);

            for(unsigned int i=1; i < p; ++i)
              *this *= tmp;

            return *this;
          }

        template<typename T>
          Vector3<T> Matrix3x3<T>::getDiagonal(void) const{
            Vector3<T> v;
            for(unsigned int i = 0; i < 3; ++i)
              v[i] = this->accessor[i][i];
            return v;
          }

        template<typename T>
          Vector3<T> Matrix3x3<T>::getRow(const unsigned int& row) const{
            if(row >= this->rowsMatrix)
              throw MathException("Out of bounds");
            Vector3<T> v;
            for(unsigned int i = 0; i < 3; ++i)
              v[i] = this->accessor[row][i];
            return v;
          }

        template<typename T>
          Vector3<T> Matrix3x3<T>::getColumn(const unsigned int& column) const{
            if(column >= this->columnsMatrix)
              throw MathException("Out of bounds");
            Vector3<T> v;
            for(unsigned int i = 0; i < 3; ++i)
              v[i] = this->accessor[i][column];
            return v;
          }

        template<typename T>
          void Matrix3x3<T>::cross3x3(const kn::Vector3<T> v){
            this->accessor[0][0] =   0.0;  this->accessor[0][1] = -v[2];  this->accessor[0][2] =  v[1];
            this->accessor[1][0] =  v[2];  this->accessor[1][1] =   0.0;  this->accessor[1][2] = -v[0];
            this->accessor[2][0] = -v[1];  this->accessor[2][1] =  v[0];  this->accessor[2][2] =   0.0;
          }


        template<typename T>
          T Matrix3x3<T>::trace(void) const{
            return (this->accessor[0][0] + this->accessor[1][1] + this->accessor[2][2]);
          }


        template<class U> Matrix3x3<U> operator* (const U& d,
            const Matrix3x3<U>& m){
          return m*d;
        }


        /*
         * Type definition
         */

        typedef Matrix3x3<float>  Matrix3x3f;
        typedef Matrix3x3<double> Matrix3x3d;
        typedef Matrix3x3<int>    Matrix3x3i;


}


#endif
