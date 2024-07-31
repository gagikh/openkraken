/*************************************************************************** \
 * Copyright (C) by University Paris-Est - MISS team
 * Vector4.hpp created in 09 2008.
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
#ifndef __OPENKN_MATH__VECTOR4_HPP__
#define __OPENKN_MATH__VECTOR4_HPP__

/*
 * Internal Includes
 */
#include "Vector3.hpp"

namespace kn{

  template<class T> class Vector4 : public Vector <T>{

    public:

      Vector4(void) : Vector<T>(4){}

      Vector4(const Vector4<T>& v) : Vector<T>(4){
        this->data[0] = v[0];
        this->data[1] = v[1];
        this->data[2] = v[2];
        this->data[3] = v[3];
      }

      Vector4(const Vector3<T>& v,const T& d1) : Vector<T>(4){
        this->data[0] = v[0];
        this->data[1] = v[1];
        this->data[2] = v[2];
        this->data[3] = d1;
      }

      Vector4(const T& d1, const Vector3<T>& v) : Vector<T>(4){
        this->data[0] = d1;
        this->data[1] = v[0];
        this->data[2] = v[1];
        this->data[3] = v[2];
      }

      Vector4(const Vector2<T>& v1,const Vector2<T>& v2) : Vector<T>(4){
        this->data[0] = v1[0];
        this->data[1] = v1[1];
        this->data[2] = v2[0];
        this->data[3] = v2[1];
      }

      Vector4(const Vector2<T>& v1, const T& d1, const T& d2) : Vector<T>(4){
        this->data[0] = v1[0];
        this->data[1] = v1[1];
        this->data[2] = d1;
        this->data[3] = d2;
      }

      Vector4(const T& d1, const T& d2, const Vector2<T>& v1) : Vector<T>(4){
        this->data[0] = d1;
        this->data[1] = d2;
        this->data[2] = v1[0];
        this->data[3] = v1[1];
      }

      Vector4(const T& d1, const Vector2<T>& v1, const T& d2) : Vector<T>(4){
        this->data[0] = d1;
        this->data[1] = v1[0];
        this->data[2] = v1[1];
        this->data[3] = d2;
      }

      Vector4(const Vector4<T>* v);

      Vector4(const T& d1, const T& d2, const T& d3, const T& d4) : Vector<T>(4){
        this->data[0] = d1;
        this->data[1] = d2;
        this->data[2] = d3;
        this->data[3] = d4;
      }

      explicit Vector4(const T& d) : Vector<T>(4, d){}

      Vector4(const T* a);

      Vector4(const Vector<T>& v);

      ~Vector4(void){}


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

      T& z(void){
        return this->data[2];
      }

      const T& z(void)const{
        return this->data[2];
      }

      T& w(void){
        return this->data[3];
      }

      const T& w(void)const{
        return this->data[3];
      }

      Vector4<T>& normalize(void);

      Vector4<T> getNormalized(void) const;

      Vector4<T>& operator=(const Vector4<T>& v);

      T operator*(const Vector4<T>& v) const;

        inline T dot(const Vector4<T>& v) const{
          return (*this)*v;
        }

      Vector4<T> operator+(const Vector4<T>& v) const;

      Vector4<T> operator-(const Vector4<T>& v) const;

      inline Vector4<T> operator/(const T& d) const{
        return Vector4(*this) /= d;
      }

      inline Vector4<T> operator*(const T& d) const{
        return Vector4(*this) *= d;
      }

      inline Vector4<T> operator-(void)const{
        return Vector4<T>(*this) * ((T)(-1));
      }

      Vector4<T>& operator+=(const Vector4<T>& v);

      Vector4<T>& operator-=(const Vector4<T>& v);

      Vector4<T>& operator/=(const T& d);

      Vector4<T>& operator*=(const T& d);


      inline size_t size(void) const{
        return 4;
      }

      inline void fill(const T &d){
        this->data[0] = this->data[1] = this->data[2] = this->data[3] = d;
      }

      inline void setZero(void){
        this->data[0] = this->data[1] = this->data[2] = this->data[3] = T(0.0);
      }

  };


  template<class T>
    Vector4<T>::Vector4(const Vector<T>& v) : Vector<T>(4){
      if(v.size() != 4)
        throw MathException("Incompatible size for Vector4");
      this->data[0] = v[0];
      this->data[1] = v[1];
      this->data[2] = v[2];
      this->data[3] = v[3];
    }

  template<class T>
    Vector4<T>::Vector4(const Vector4<T>* v) : Vector<T>(4){
      if(v==0) throw MathException("Pointer is null");
        this->data[0] = (*v)[0];
        this->data[1] = (*v)[1];
        this->data[2] = (*v)[2];
        this->data[3] = (*v)[3];
    }

  template<class T>
    Vector4<T>::Vector4(const T* a) : Vector<T>(4){
      if(a==0) throw MathException("Pointer is null");
      this->data[0] = a[0];
      this->data[1] = a[1];
      this->data[2] = a[2];
      this->data[3] = a[3];
    }

  template<class T>
    Vector4<T>& Vector4<T>::normalize(void){
      float norm = this->getNorm();
      this->data[0] = T(this->data[0] / norm);
      this->data[1] = T(this->data[1] / norm);
      this->data[2] = T(this->data[2] / norm);
      this->data[3] = T(this->data[3] / norm);
      return *this;
    }

  template<class T>
    Vector4<T> Vector4<T>::getNormalized(void) const{
      Vector4<T> v(*this);
      v.normalize();
      return v;
    }

  template<class T>
    T Vector4<T>::operator*(const Vector4<T>& v) const{
      return v[0]*this->data[0] + v[1]*this->data[1] + v[2]*this->data[2] + v[3]*this->data[3];
    }

  template<class T>
    Vector4<T>& Vector4<T>::operator=(const Vector4<T>& v){
      if(&v == this) return *this;
      this->data[0] = v[0];
      this->data[1] = v[1];
      this->data[2] = v[2];
      this->data[3] = v[3];
      return *this;
    }


  template<class T>
    Vector4<T> Vector4<T>::operator+(const Vector4<T>& v) const{
      Vector4 tmp;
      tmp[0] = this->data[0] + v[0];
      tmp[1] = this->data[1] + v[1];
      tmp[2] = this->data[2] + v[2];
      tmp[3] = this->data[3] + v[3];
      return tmp;
    }


  template<class T>
    Vector4<T> Vector4<T>::operator-(const Vector4<T>& v) const{
      Vector4 tmp;
      tmp[0] = this->data[0] - v[0];
      tmp[1] = this->data[1] - v[1];
      tmp[2] = this->data[2] - v[2];
      tmp[3] = this->data[3] - v[3];
      return tmp;
    }

  template<class T>
    Vector4<T>& Vector4<T>::operator+=(const Vector4<T>& v){
      this->data[0] += v[0];
      this->data[1] += v[1];
      this->data[2] += v[2];
      this->data[3] += v[3];
      return *this;
    }

  template<class T>
    Vector4<T>& Vector4<T>::operator-=(const Vector4<T>& v){
      this->data[0] -= v[0];
      this->data[1] -= v[1];
      this->data[2] -= v[2];
      this->data[3] -= v[3];
      return *this;
    }

  template<class T>
    Vector4<T>& Vector4<T>::operator/=(const T& d){
      this->data[0] /= d;
      this->data[1] /= d;
      this->data[2] /= d;
      this->data[3] /= d;
      return *this;
    }

  template<class T>
    Vector4<T>& Vector4<T>::operator*=(const T& d){
      this->data[0] *= d;
      this->data[1] *= d;
      this->data[2] *= d;
      this->data[3] *= d;
      return *this;
    }

  template<class U> Vector4<U> operator*(const U& d,
      const Vector4<U>& v){
    Vector4<U> result(v);
    return result*d;
  }


  /*
   * Type definition
   */

  typedef Vector4<float>  Vector4f;
  typedef Vector4<double> Vector4d;
  typedef Vector4<int>    Vector4i;
}


#endif