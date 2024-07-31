/*************************************************************************** \
 * Copyright (C) by University Paris-Est - MISS team
 * Vector3.hpp created in 09 2008.
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
#ifndef __OPENKN_MATH__VECTOR3_HPP__
#define __OPENKN_MATH__VECTOR3_HPP__

/*
 * Internal Includes
 */
#include "Vector2.hpp"

namespace kn{


  template<class T> class Vector3 : public Vector <T>{

    public:

      Vector3(void) : Vector<T>(3){ }

      Vector3(const Vector3<T>& v) : Vector<T>(3){
        this->data[0] = v[0];
        this->data[1] = v[1];
        this->data[2] = v[2];
      }


      Vector3(const Vector3<T>* v);


      Vector3(const Vector2<T>& v,const T& d1) : Vector<T>(3){
        this->data[0] = v[0];
        this->data[1] = v[1];
        this->data[2] = d1;
      }

      Vector3(const T& d1, const Vector2<T>& v) : Vector<T>(3){
        this->data[0] = d1;
        this->data[1] = v[0];
        this->data[2] = v[1];
      }


      Vector3(const T& d1, const T& d2, const T& d3) : Vector<T>(3){
        this->data[0] = d1;
        this->data[1] = d2;
        this->data[2] = d3;
      }

      Vector3(const T& d) : Vector<T>(3, d){}

      Vector3(const T* a);

      Vector3(const Vector<T>& v);

      ~Vector3(void){}


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

      //using Vector<T>::operator*;
      Vector3<T>& normalize(void);

      Vector3<T> getNormalized(void) const;

      T operator*(const Vector3<T>& v) const;

      T dot(const Vector3<T>& v) const{
        return (*this)*v;
      }

      Vector3<T>& operator=(const Vector3<T>& v);

      Vector3<T> operator+(const Vector3<T>& v) const;

      Vector3<T> operator-(const Vector3<T>& v) const;

      inline Vector3<T> operator/(const T& d) const{
        return Vector3(*this) /= d;
      }

      inline Vector3<T> operator*(const T& d) const{
        return Vector3(*this) *= d;
      }

      inline Vector3<T> operator-(void)const{
        return Vector3<T>(*this) * ((T)(-1));
      }

      Vector3<T>& operator+=(const Vector3<T>& v);

      Vector3<T>& operator-=(const Vector3<T>& v);

      Vector3<T>& operator/=(const T& d);

      Vector3<T>& operator*=(const T& d);


      Vector3<T> operator^(const Vector3<T>& v)const;

      inline Vector3<T> cross(const Vector3<T>& v)const{
        return (*this)^v;
      }

      inline size_t size(void) const{
        return 3;
      }

      inline void fill(const T &d){
        this->data[0] = this->data[1] = this->data[2] = d;
      }

      inline void setZero(void){
        this->data[0] = this->data[1] = this->data[2] = T(0.0);
      }
  };


  template<class T>
    Vector3<T>::Vector3(const Vector<T>& v) :Vector<T>(3){
      if(v.size() != 3)
        throw MathException("Incompatible size for Vector3");
      this->data[0] = v[0];
      this->data[1] = v[1];
      this->data[2] = v[2];
    }

  template<class T>
    Vector3<T>::Vector3(const Vector3<T>* v) : Vector<T>(3){
      if(v==0) throw MathException("Pointer is null");
      this->data[0] = (*v)[0];
      this->data[1] = (*v)[1];
      this->data[2] = (*v)[2];
    }

  template<class T>
    Vector3<T>::Vector3(const T* a) : Vector<T>(3){
      if(a==0) throw MathException("Pointer is null");
      this->data[0] = a[0];
      this->data[1] = a[1];
      this->data[2] = a[2];
    }

  template<class T>
    Vector3<T>& Vector3<T>::normalize(void){
      float norm = this->getNorm();
      this->data[0] = T(this->data[0] / norm);
      this->data[1] = T(this->data[1] / norm);
      this->data[2] = T(this->data[2] / norm);
      return *this;
    }

  template<class T>
    Vector3<T> Vector3<T>::getNormalized(void) const{
      Vector3<T> v(*this);
      v.normalize();
      return v;
    }

  template<class T>
    T Vector3<T>::operator*(const Vector3<T>& v) const{
      return v[0]*this->data[0] + v[1]*this->data[1] + v[2]*this->data[2];
    }



  template<class T>
    Vector3<T>& Vector3<T>::operator=(const Vector3<T>& v){
      if(&v == this) return *this;
      this->data[0] = v[0];
      this->data[1] = v[1];
      this->data[2] = v[2];
      return *this;
    }


  template<class T>
    Vector3<T> Vector3<T>::operator+(const Vector3<T>& v) const{
      Vector3 tmp;
      tmp[0] = this->data[0] + v[0];
      tmp[1] = this->data[1] + v[1];
      tmp[2] = this->data[2] + v[2];
      return tmp;
    }


  template<class T>
    Vector3<T> Vector3<T>::operator-(const Vector3<T>& v) const{
      Vector3 tmp;
      tmp[0] = this->data[0] - v[0];
      tmp[1] = this->data[1] - v[1];
      tmp[2] = this->data[2] - v[2];
      return tmp;
    }

  template<class T>
    Vector3<T>& Vector3<T>::operator+=(const Vector3<T>& v){
      this->data[0] += v[0];
      this->data[1] += v[1];
      this->data[2] += v[2];
      return *this;
    }

  template<class T>
    Vector3<T>& Vector3<T>::operator-=(const Vector3<T>& v){
      this->data[0] -= v[0];
      this->data[1] -= v[1];
      this->data[2] -= v[2];
      return *this;
    }

  template<class T>
    Vector3<T>& Vector3<T>::operator/=(const T& d){
      if(d==T(0)) throw MathException("variable is null");
      this->data[0] /= d;
      this->data[1] /= d;
      this->data[2] /= d;
      return *this;
    }

  template<class T>
    Vector3<T>& Vector3<T>::operator*=(const T& d){
      this->data[0] *= d;
      this->data[1] *= d;
      this->data[2] *= d;
      return *this;
    }

  template<class T>
    Vector3<T> Vector3<T>::operator^(const Vector3<T>& v)const{
      Vector3<T> tmp;
      tmp[0] = this->data[1]*v[2] - this->data[2]*v[1];
      tmp[1] = this->data[2]*v[0] - this->data[0]*v[2];
      tmp[2] = this->data[0]*v[1] - this->data[1]*v[0];
      return tmp;
    }




  template<class U> Vector3<U> operator*(const U& d,
      const Vector3<U>& v){
    Vector3<U> result(v);
    return result*d;
  }


  /*
   * Type definition
   */

  typedef Vector3<float>  Vector3f;
  typedef Vector3<double> Vector3d;
  typedef Vector3<int>    Vector3i;
}


#endif
