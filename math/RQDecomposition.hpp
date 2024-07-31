/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * RQDecomposition.hpp created in 01 2009.
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
#ifndef __OPENKN_MATH__RQ_DECOMPOSITION_HPP__
#define __OPENKN_MATH__RQ_DECOMPOSITION_HPP__


/*
 * External Includes
 */
#include <cmath>

/*
 * Internal Includes
 */
#include "MathException.hpp"
#include "Matrix.hpp"
#include "Matrix3x3.hpp"


/*
 * Namespace
 */
namespace kn {


  template <class T>
    void qrDecomposition(const Matrix<T> &A, Matrix<T> &R, Matrix<T> &Q){
      if(A.rows() != A.columns()) throw MathException("the input matrix is not a square matrix");
      if(R.rows() != R.columns()) throw MathException("the 'R' matrix is not a square matrix");
      if(Q.rows() != Q.columns()) throw MathException("the 'Q' matrix is not a square matrix");
      if(A.rows() != R.rows() || A.rows() != Q.rows()) throw MathException("input matrix does not have the same size");

      R=A;
      Q.setIdentity();
      Matrix<T> Qtmp(Q.rows());

      //During all this loop we consider Q to be transpose(Q)
      for(unsigned int i=0;i<A.columns()-1;i++){
        //v=first column of sub problem
        Vector<T> v=R.getColumn(i).getSubVector(i,R.rows()-i);
        //v = v +- (||v||,0,0,.....)
        //the sign is the same as Vo
        if(v.at(0)<0.0) {v.at(0) -= v.getNorm();}
        else            {v.at(0) += v.getNorm();}
        //v = v/||v||
        if(v.getNorm()<10e-30){
          throw MathException("Division by zero");
        }
        v.normalize();

        //we resolve the QR decomposition for the sub problem
        Matrix<T> Qn(v.size(),v.size());
        Qn.setIdentity();
        //Qn = Id - 2*tranpose(v)*v;
        for(unsigned int j=0;j<v.size();j++){
          for(unsigned int k=j;k<v.size();k++){
            Qn.at(j,k) -= 2*v.at(j)*v.at(k);
            Qn.at(k,j) = Qn.at(j,k);//symetric matrix
          }
        }
        Qn=-Qn;

        //Qtmp = Id or resolution of the sub problem in the bottom right corner of the matrix
        Qtmp.setIdentity();
        Qtmp.setSubMatrix(i,i,Qn);
        //transpose(Q) = Qn Qn-1 .... Q1
        Q = Qtmp*Q;

        //R = transpose(Q)*A
        R = Q*A;
        R.roundZero(10e-10);

      }
      Q.transpose();
    }


  template <class T>
    void rqDecomposition3x3(const Matrix<T> &A, Matrix<T> &R, Matrix<T> &Q){
      if(A.rows() != 3 || A.columns() != 3) throw MathException("the input matrix is not a 3x3 matrix");
      if(R.rows() != 3 || R.columns() != 3) throw MathException("the 'R' matrix is not a 3x3 matrix");
      if(Q.rows() != 3 || Q.columns() != 3) throw MathException("the 'Q' matrix is not a 3x3 matrix");

      Matrix3x3<T> Qx;
      Matrix3x3<T> Qy;
      Matrix3x3<T> Qz;

      Qx.setIdentity();
      Qy.setIdentity();
      Qz.setIdentity();


      R = A;

      Qx[2][2] = -R[2][2]/sqrt(pow(R[2][1],2)+pow(R[2][2],2));
      Qx[1][1] = Qx[2][2];
      Qx[2][1] = R[2][1]/sqrt(pow(R[2][1],2)+pow(R[2][2],2));
      Qx[1][2] = -Qx[2][1];

      R = R*Qx;

      // Note : there are 2 equivalent solutions for "sqrt(pow(R[][],2)+pow(R[][],2));" correpsonding to the angle
      // "+ alpha" and "- 360° + alpha". The selected solution is "+sqrt(pow(R[][],2)+pow(R[][],2));".

      Qy[2][2] = R[2][2]/sqrt(pow(R[2][0],2)+pow(R[2][2],2));
      Qy[0][0] = Qy[2][2];
      Qy[0][2] = R[2][0]/sqrt(pow(R[2][0],2)+pow(R[2][2],2));
      Qy[2][0] = -Qy[0][2];

      R = R*Qy;

      Qz[0][0] = R[1][1]/sqrt(pow(R[1][1],2)+pow(R[1][0],2));
      Qz[1][1] = Qz[0][0];
      Qz[1][0] = -R[1][0]/sqrt(pow(R[1][1],2)+pow(R[1][0],2));
      Qz[0][1] = -Qz[1][0];

      R = R*Qz;

      Qx.transpose();
      Qy.transpose();
      Qz.transpose();
      Q = Qz * Qy * Qx;
    }


  template <class T>
    void rq3x3MakePositiveDiagonal(Matrix<T> &R, Matrix<T> &Q){
      if(R.rows() != 3 || R.columns() != 3) throw MathException("the 'R' matrix is not a 3x3 matrix");
      if(Q.rows() != 3 || Q.columns() != 3) throw MathException("the 'Q' matrix is not a 3x3 matrix");

      // variations pour garantir R[1][1] >= 0
      if(R[1][1] < 0)
      {
        for(int i=0; i<3; i++) R[i][1] *= -1.0;
        for(int i=0; i<3; i++) Q[1][i] *= -1.0;
      }

      // variations pour garantir R[0][0] >= 0
      if(R[0][0] < 0)
      {
        for(int i=0; i<3; i++) R[i][0] *= -1.0;
        for(int i=0; i<3; i++) Q[0][i] *= -1.0;
      }
    }

  /*
   * End of Namespace
   */
}

/*
 * End of Anti-doublon
 */
#endif






