/*************************************************************************** \
 * Copyright (C) by University Paris-Est - MISS team
 * Vector2.hpp created in 12 2008.
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
#ifndef __OPENKN_MATH__VECTOR2_HPP__
#define __OPENKN_MATH__VECTOR2_HPP__

/*
 * Internal Includes
 */
#include "Vector.hpp"

namespace kn{

  template<class T> class Vector2 : public Vector <T>{

    public:

      Vector2(void) : Vector<T>(2){ }

      Vector2(const Vector2<T>& v) : Vector<T>(2){
        this->data[0] = v[0];
        this->data[1] = v[1];
      }


      Vector2(const Vector2<T>* v);

      Vector2(const T& d1, const T& d2) : Vector<T>(2){
        this->data[0] = d1;
        this->data[1] = d2;
      }

      explicit Vector2(const T& d) : Vector<T>(2, d){}

      Vector2(const T* a);

      Vector2(const Vector<T>& v);

      ~Vector2(void){}


    public:

      /*
       * Accessor
       */

      T& x(void){
        return this->data[0];
      }

      const T& x(void)const{
        return this->data[0];
      }

      T& y(void){
        return this->data[1];
      }

      const T& y(void)const{
        return this->data[1];
      }


      //using Vector<T>::operator*;

      Vector2<T>& normalize(void);

      Vector2<T> getNormalized(void) const;


      T operator*(const Vector2<T>& v) const;

        inline T dot(const Vector2<T>& v) const{
          return (*this)*v;
        }


      Vector2<T>& operator=(const Vector2<T>& v);

      Vector2<T> operator+(const Vector2<T>& v) const;

      Vector2<T> operator-(const Vector2<T>& v) const;

      inline Vector2<T> operator/(const T& d) const{
        return Vector2(*this) /= d;
      }

      inline Vector2<T> operator*(const T& d) const{
        return Vector2(*this) *= d;
      }

      inline Vector2<T> operator-(void)const{
        return Vector2<T>(*this) * ((T)(-1));
      }

      Vector2<T>& operator+=(const Vector2<T>& v);

      Vector2<T>& operator-=(const Vector2<T>& v);

      Vector2<T>& operator/=(const T& d);

      Vector2<T>& operator*=(const T& d);



      inline size_t size(void) const{
        return 2;
      }

      inline void fill(const T &d){
        this->data[0] = this->data[1] = d;
      }

      inline void setZero(void){
        this->data[0] = this->data[1] = T(0.0);
      }
  };


  template<class T>
    Vector2<T>::Vector2(const Vector<T>& v) :Vector<T>(2){
      if(v.size() != 2)
        throw MathException("Incompatible size for Vector2");
      this->data[0] = v[0];
      this->data[1] = v[1];
    }

  template<class T>
    Vector2<T>::Vector2(const Vector2<T>* v) : Vector<T>(2){
      if(v==0) throw MathException("Pointer is null");
      this->data[0] = (*v)[0];
      this->data[1] = (*v)[1];
    }

  template<class T>
    Vector2<T>::Vector2(const T* a) : Vector<T>(2){
      if(a==0) throw MathException("Pointer is null");
      this->data[0] = a[0];
      this->data[1] = a[1];
    }


  template<class T>
    Vector2<T>& Vector2<T>::normalize(void){
      float norm = this->getNorm();
      this->data[0] = T(this->data[0] / norm);
      this->data[1] = T(this->data[1] / norm);
      return *this;
    }

  template<class T>
    Vector2<T> Vector2<T>::getNormalized(void) const{
      Vector2<T> v(*this);
      v.normalize();
      return v;
    }

  template<class T>
    T Vector2<T>::operator*(const Vector2<T>& v) const{
      return v[0]*this->data[0] + v[1]*this->data[1];
    }


  template<class T>
    Vector2<T>& Vector2<T>::operator=(const Vector2<T>& v){
      if(&v == this) return *this;
      this->data[0] = v[0];
      this->data[1] = v[1];
      return *this;
    }


  template<class T>
    Vector2<T> Vector2<T>::operator+(const Vector2<T>& v) const{
      Vector2 tmp;
      tmp[0] = this->data[0] + v[0];
      tmp[1] = this->data[1] + v[1];
      return tmp;
    }


  template<class T>
    Vector2<T> Vector2<T>::operator-(const Vector2<T>& v) const{
      Vector2 tmp;
      tmp[0] = this->data[0] - v[0];
      tmp[1] = this->data[1] - v[1];
      return tmp;
    }

  template<class T>
    Vector2<T>& Vector2<T>::operator+=(const Vector2<T>& v){
      this->data[0] += v[0];
      this->data[1] += v[1];
      return *this;
    }

  template<class T>
    Vector2<T>& Vector2<T>::operator-=(const Vector2<T>& v){
      this->data[0] -= v[0];
      this->data[1] -= v[1];
      return *this;
    }

  template<class T>
    Vector2<T>& Vector2<T>::operator/=(const T& d){
      this->data[0] /= d;
      this->data[1] /= d;
      return *this;
    }

  template<class T>
    Vector2<T>& Vector2<T>::operator*=(const T& d){
      this->data[0] *= d;
      this->data[1] *= d;
      return *this;
    }

  template<class U> Vector2<U> operator*(const U& d,
      const Vector2<U>& v){
    Vector2<U> result(v);
    return result*d;
  }


  /*
   * Type definition
   */

  typedef Vector2<float>  Vector2f;
  typedef Vector2<double> Vector2d;
  typedef Vector2<int>    Vector2i;
  typedef Vector2<unsigned int>    Vector2u;
}

#endif