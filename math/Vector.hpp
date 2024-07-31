/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * Vector.hpp created in 09 2008.
 * Mail : biri@univ-mlv.fr
 *
 * This file is part of the OpenKraken-image.
 *
 * The OpenKraken-image is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * The OpenKraken-image is distributed in the hope that it will be useful,
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
#ifndef __OPENKN_MATH__VECTOR_HPP__
#define __OPENKN_MATH__VECTOR_HPP__

/* Removing some warnings on Microsoft Visual C++ */
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma warning(disable:4996)
#endif

/*
 * External Includes
 */
#include <iostream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>


/*
 * Internal Includes
 */
#include "MathException.hpp"


/*
 * Namespace
 */
namespace kn{

  /*
   * Class definition
   */

  template<class T>
    class Vector {

      /*
       * class variable
       */
      protected :
        size_t sizeVector;

        T *begin_;

        T *end_;

        T *data;


        /*
         * Constructor & destructors
         */
      public:

        Vector();

        explicit Vector(const size_t& size);

        Vector(const Vector<T> & v);

        Vector(Vector<T> * v);

        Vector(const size_t& dataSize,
            const T *a);

        Vector(const size_t& dataSize,
            const T &d);

        /*
         * Protected methods
         */
      protected :
        void allocate(const size_t& d);

        inline void desallocate(void){
          if(sizeVector != 0 && data != 0) delete[] data;
          sizeVector = 0;
        }

        /*
         * Public methods
         */
      public :

        virtual inline ~Vector(void){
          desallocate();
        }

        T& at(const unsigned int& i);

        const T& at(const unsigned int& i) const;

        void resize(const int& size);

        void resizeAndCopy(const Vector<T>& v);

        /*
         * Operators redefinition
         */
        bool operator==(const Vector<T>& v) const;

        inline bool operator!=(const Vector<T>& v) const{
          return !(*this == v);
        }

        Vector<T>& operator=(const Vector<T>& v);

        Vector<T>& operator=(const T & value);

        Vector<T> operator+(const Vector<T>& v) const;


        Vector<T> operator-(const Vector<T>& v) const;

        inline Vector<T> operator/(const T& d) const{
          return Vector(*this) /= d;
        }

        T operator*(const Vector<T>& v) const;

        inline T dot(const Vector<T>& v) const{
                                        return (*this)*v;
                                }

        inline Vector<T> operator*(const T& d) const{
          return Vector(*this) *= d;
        }

        inline Vector<T> operator-(void)const{
          return Vector<T>(*this) * ((T)(-1));
        }

        Vector<T>& operator+=(const Vector<T>& v);

        Vector<T>& operator-=(const Vector<T>& v);

        Vector<T>& operator/=(const T& d);

        Vector<T>& operator*=(const T& d);

        T& operator[](const unsigned int& i);

        const T& operator[](const unsigned int& i) const;

        T& operator()(const unsigned int& i);

        const T& operator()(const unsigned int& i) const;

        Vector<T> operator^(const Vector<T>& v)const;

        inline Vector<T> cross(const Vector<T>& v)const{
                                        return (*this)^v;
                                }



        /*
         * Vector methods
         */

        inline size_t size(void) const{
          return sizeVector;
        }

        inline void fill(const T &d){
          std::fill(begin_,end_,d);
        }

        inline void setZero(void){
          fill(T(0.0));
        }

        void roundZero(const double & d = 1.0e-14);

        void setSubVector(const unsigned int& index, const Vector<T>& v);

        void setSubVector(const unsigned int& index1, const unsigned int& index2,const unsigned int& size, const Vector<T>& v);
        Vector<T> getSubVector(const unsigned int& index, const int& size) const;


        inline virtual T* begin(void) const{
          return begin_;
        }


        inline virtual T* end(void) const{
          return end_;
        }


        double getNorm(void) const;

        Vector<T>& normalize(void);

        Vector<T> getNormalized(void) const;

        void swap(const unsigned int& index1, const unsigned int& index2);

        Vector<T> getHomogeneous(const T &d = T(1.0)) const;


        Vector<T> getUnhomogeneous(const T &zeroValue = T(1.0e-13)) const;


        void setHomogeneousNormalForm(const T &zeroValue = T(1.0e-13)) const;

    };

  template<class T>
    Vector<T>::Vector(){
      sizeVector = 0;
    }

  template<class T>
    Vector<T>::Vector(const size_t& size){
      allocate(size);
      setZero();
    }

  template<class T>
    Vector<T>::Vector(const Vector<T> & v){
      allocate(v.sizeVector);

      std::copy(v.begin_, v.end_, begin_);
    }


  template<class T>
    Vector<T>::Vector(Vector<T> * v){
      if(v == 0) throw MathException("Null pointer exception");

      allocate(v->sizeVector);

      std::copy(v->begin_, v->end_, this->begin_);
    }

  template<class T>
    Vector<T>::Vector(const size_t& dataSize,
        const T *a){
      if(a == 0) throw MathException("Null pointer exception");

      allocate(dataSize);

      std::copy(a,a+dataSize,this->begin_);
    }


  template<class T>
    Vector<T>::Vector(const size_t& dataSize,
        const T &d){
      allocate(dataSize);
      fill(d);
    }

  template<class T>
    void Vector<T>::resize(const int& size){
      if(sizeVector != size) desallocate();
      if(sizeVector == 0) allocate(size);
      std::fill(begin_, end_, T(0));
    }

  template<class T>
    void Vector<T>::resizeAndCopy(const Vector<T>& v){
      if(&v == this) return;
      if(sizeVector != v.sizeVector) desallocate();
      if(sizeVector == 0) allocate(v.sizeVector);
      std::copy(v.begin_, v.end_, begin_);
    }

  template<class T>
    void Vector<T>::allocate(const size_t& d){
      if(d == 0) throw MathException("vector size is 0");

      sizeVector = d;
      data         = new T[d];
      begin_       = data;
      end_         = begin_ + d;
    }

  template<class T>
    T& Vector<T>::at(const unsigned int& i){
      if(i >= sizeVector) throw MathException("Index out of bounds");
      return data[i];
    }

  template<class T>
    const T& Vector<T>::at(const unsigned int& i) const{
      if(i >= sizeVector) throw MathException("Index out of bounds");
      return data[i];
    }

  template<class T>
    bool Vector<T>::operator==(const Vector<T>& v) const{
      if(v.sizeVector != sizeVector) return false;

      return std::equal(begin_, end_, v.begin_);
    }


  template<class T>
    Vector<T>& Vector<T>::operator=(const Vector<T>& v){
      if(&v == this) return *this;

      if(sizeVector == 0) allocate(v.sizeVector);

      if(sizeVector != v.sizeVector) throw MathException("incompatible vector size");

      std::copy(v.begin_, v.end_, begin_);

      return *this;
    }

  template<class T>
    Vector<T>& Vector<T>::operator=(const T & value) {
      std::fill(begin_,end_,T(value));

      return *this;
    }

  template<class T>
    Vector<T> Vector<T>::operator+(const Vector<T>& v) const{
      if(sizeVector != v.sizeVector) throw MathException("incompatible vector size");

      Vector tmp(sizeVector);

      std::transform(begin_, end_, v.begin_, tmp.begin_, std::plus<T>());

      return tmp;
    }


  template<class T>
    Vector<T> Vector<T>::operator-(const Vector<T>& v) const{
      if(sizeVector != v.sizeVector) throw MathException("incompatible vector size");

      Vector tmp(sizeVector);

      std::transform(begin_, end_, v.begin_, tmp.begin_, std::minus<T>());

      return tmp;
    }

  template<class T>
    T Vector<T>::operator*(const Vector<T>& v) const{
      if(sizeVector != v.sizeVector) throw MathException("incompatible vector size");

      return std::inner_product(begin_, end_, v.begin_,T(0.0));
    }

  template<class T>
    Vector<T>& Vector<T>::operator+=(const Vector<T>& v)
    {
      if(sizeVector != v.sizeVector) throw MathException("incompatible vector size");

      std::transform(begin_, end_, v.begin_, begin_, std::plus<T>()); 

      return *this;
    }

  template<class T>
    Vector<T>& Vector<T>::operator-=(const Vector<T>& v){
      if(sizeVector != v.sizeVector) throw MathException("incompatible vector size");

      std::transform(begin_, end_, v.begin_, begin_, std::minus<T>()); 

      return *this;
    }

  template<class T>
    Vector<T>& Vector<T>::operator/=(const T& d){
      //std::transform(begin_, end_, begin_, std::bind2nd(std::divides<T>(), d));
      auto lambda = std::bind(std::divides<T>(), std::placeholders::_1, d);
      std::transform(begin_, end_, begin_, lambda);
      return *this;
    }

  template<class T>
    Vector<T>& Vector<T>::operator*=(const T& d){
      //std::transform(begin_, end_, begin_, std::bind2nd(std::multiplies<T>(), d));
      auto lambda = std::bind(std::multiplies<T>(), std::placeholders::_1, d);
      std::transform(begin_, end_, begin_, lambda);
      return *this;
    }

  template<class T>
    T& Vector<T>::operator[](const unsigned int& i){
      return data[i];
    }

  template<class T>
    const T& Vector<T>::operator[](const unsigned int& i) const{
      return data[i];
    }

  template<class T>
    T& Vector<T>::operator()(const unsigned int& i){
      return data[i];
    }

  template<class T>
    const T& Vector<T>::operator()(const unsigned int& i) const{
      return data[i];
    }

  template<class T>
    Vector<T> Vector<T>::operator^(const Vector<T>& v) const{
      if(sizeVector!= 3) throw MathException("(this) cross product only with vectors of size 3");
      if(v.sizeVector != 3) throw MathException("cross product only with vectors of size 3");

      Vector<T> tmp(3);
      tmp[0] = data[1]*v[2] - data[2]*v[1];
      tmp[1] = data[2]*v[0] - data[0]*v[2];
      tmp[2] = data[0]*v[1] - data[1]*v[0];

      return tmp;
    }



  template<class T>
    void Vector<T>::roundZero(const double & d){
      T *iter = begin_;

      while(iter != end_){
        if(std::fabs((double)(*iter)) < d) (*iter) = T(0.0);
        ++iter;
      }
    }


  template<class T>
    void Vector<T>::setSubVector(const unsigned int& index, const Vector<T>& v){
      if(index >= sizeVector) throw MathException("index is out of bounds");
      if(index + v.sizeVector > sizeVector) throw MathException("index is out of bounds(v.sizeVector is too big)");
      std::copy(v.begin_, v.end_, begin_ + index);
    }

  template<class T>
    void Vector<T>::setSubVector(const unsigned int& index1, const unsigned int& index2,const unsigned int& size, const Vector<T>& v){
      if(index1 >= sizeVector) throw MathException("index1 is out of bounds");
      if(index1 + size > sizeVector) throw MathException("index1 is out of bounds(size is too big)");
      if(index2 >= v.sizeVector) throw MathException("index2 is out of bounds");
      if(index2 + size > v.sizeVector) throw MathException("index2 is out of bounds(size is too big)");
      std::copy(v.begin_+index2, v.begin_+index2+size, begin_ + index1);
  }

  template<class T>
    Vector<T> Vector<T>::getSubVector(const unsigned int& index, const int& size) const{
      if(index >= sizeVector) throw MathException("index is out of bounds");
      if(index + size > sizeVector) throw MathException("index is out of bounds (size is too big)");

      Vector<T> v(size);
      std::copy(begin_ + index, begin_ + index + size , v.begin_);
      return v;
    }

  template<class T>
    double Vector<T>::getNorm(void) const{
      double squareSum = 0.0;
      for(unsigned int i=0;i<sizeVector;i++){
        squareSum += (double)(this->data[i]*this->data[i]);
      }
      return sqrt(squareSum);
    }

  template<class T>
    Vector<T>& Vector<T>::normalize(void){
      double norm = getNorm();
      for(unsigned int i=0;i<sizeVector;i++){
        this->data[i] = T(this->data[i] / norm);
      }
      return *this;
    }

  template<class T>
    Vector<T> Vector<T>::getNormalized(void) const{
      Vector<T> v(*this);
      v.normalize();
      return v;
    }


  template<class T>
    Vector<T> Vector<T>::getHomogeneous(const T &d) const{
      Vector<T> tmp(sizeVector+1);

      std::copy(begin_, end_, tmp.begin_);
      tmp[sizeVector] = d;

      return tmp;
    }


  template<class T>
    Vector<T> Vector<T>::getUnhomogeneous(const T &zeroValue) const{
      if(sizeVector == 1) return *this;

      Vector<T> tmp(sizeVector-1);

      // unhomogenize
      T value = this->data[sizeVector-1];
      if (std::fabs((float)value) > zeroValue && value != T(1.0)) {
          //std::transform(begin_, end_ - 1, tmp.begin_, std::bind2nd(std::divides<T>(), value));
          auto lambda = std::bind(std::divides<T>(), std::placeholders::_1, value);
          std::transform(begin_, end_ - 1, tmp.begin_, lambda);
      }
      else
        std::copy(begin_, end_-1, tmp.begin_);

      return tmp;
    }


  template<class T>
    void Vector<T>::setHomogeneousNormalForm(const T &zeroValue) const{
      T value = this->data[sizeVector-1];
      if(std::fabs((float)value) > (float)zeroValue && value != T(1.0)){
        //std::transform(begin_, end_-1, begin_, std::bind2nd(std::divides<T>(), value));
        auto lambda = std::bind(std::divides<T>(), std::placeholders::_1, value);
        std::transform(begin_, end_ - 1, begin_, lambda);
        this->data[sizeVector-1] = T(1.0);
      }
    }


  template<class T>
    void Vector<T>::swap(const unsigned int& index1,
        const unsigned int& index2){
      if(index1 >= sizeVector || index2 >= sizeVector) throw MathException("Index out of bounds");

      std::swap(data[index1],data[index2]);
    }

  template<class U>
    std::ostream& operator<< (std::ostream& stream,
        const Vector<U>& v){
      stream << "(";
      for (unsigned int i=0; i<v.size()-1; ++i)
        stream << v[i] << " , ";

      stream << v[v.size()-1]<< ")";

      return stream;
    }

  template<class U>
    Vector<U> operator* (const U& d,
        const Vector<U>& v){
      return v*d;
    }

  /*
   * Type definition
   */

  typedef Vector<float>  Vectorf;
  typedef Vector<double> Vectord;
  typedef Vector<int>      Vectori;

  /*
   * End of Namespace
   */
}

/*
 * End of Anti-doublon
 */
#endif
