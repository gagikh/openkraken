/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * determinant.hpp created in 01 2009.
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
#ifndef __OPENKN_MATH__DETERMINANT_HPP__
#define __OPENKN_MATH__DETERMINANT_HPP__

/*
 * External Includes
 */
#include <algorithm>

/*
 * Internal Includes
 */
#include "MathException.hpp"
#include "Matrix.hpp"
#include "Matrix3x3.hpp"


/*
 * Namespace
 */
namespace kn{

  template <class T>
    T determinant3x3(Matrix<T> &mat){
      if(mat.rows() != 3 || mat.columns() != 3)
        throw MathException("the input matrix is not a 3x3 matrix");

      return mat[0][0]*mat[1][1]*mat[2][2] -
        mat[0][0]*mat[2][1]*mat[1][2] -
        mat[0][1]*mat[1][0]*mat[2][2] +
        mat[0][1]*mat[2][0]*mat[1][2] +
        mat[0][2]*mat[1][0]*mat[2][1] -
        mat[0][2]*mat[2][0]*mat[1][1];
    }


  template <class T>
    T determinant3x3(Matrix3x3<T> &mat){

      return mat[0][0]*mat[1][1]*mat[2][2] -
        mat[0][0]*mat[2][1]*mat[1][2] -
        mat[0][1]*mat[1][0]*mat[2][2] +
        mat[0][1]*mat[2][0]*mat[1][2] +
        mat[0][2]*mat[1][0]*mat[2][1] -
        mat[0][2]*mat[2][0]*mat[1][1];
    }


  template <class T>
  T determinant(const Matrix<T> & A_){
    // thresold to decide whether a pivot is zero or not
    const double determinantZero = 1.0e-10;

    // check the size compatibility
    if(A_.columns() != A_.rows())
      throw MathException("The input matrix is not a square matrix" ,"determinant");

    // make a copy
    Matrix<T> A(A_);

    // permuation
    T permutation = T(1);

    // make A triangular
    unsigned int nbIteration = 0;
    while(nbIteration < A.rows()-1){

      // find the pivot
      T maxValue;
      unsigned int maxValueRow,maxValueCol;
      determinantFindPivot(A,nbIteration,maxValueRow,maxValueCol,maxValue);

      // check if the "pivot" is non-null
      if(fabs(maxValue) < determinantZero)
        return T(0); // the matrix is singular (not full rank) => det(A) = 0.

      // if required, swap rows
      if(maxValueRow != nbIteration){
        // swap the non zero part of the rows
        std::swap_ranges(&(A[maxValueRow][nbIteration]),A[maxValueRow]+A.columns(),&(A[nbIteration][nbIteration]));
        // update permutation
        permutation *= T(-1);
      }

      // if required, swap columns
      if(maxValueCol != nbIteration){
        // swap columns
        A.swapColumns(maxValueCol,nbIteration);
        // update permutation
        permutation *= T(-1);
      }

      // row alimination (let's put some 0)
      determinantRowsElimination(A,nbIteration);

      // next iteration
      ++nbIteration;
    }

    // determinant
    T determinant = T(1);
    for(unsigned int i=0; i<A.rows(); ++i)
      determinant = determinant * A[i][i];

    return determinant * permutation;
  }


  template <class T>
  void determinantFindPivot(const Matrix<T> & A, const unsigned int nbIteration,
                            unsigned int & maxValueRow, unsigned int & maxValueCol, T & maxValue){

      maxValueRow = maxValueCol = nbIteration;
      maxValue = A[nbIteration][nbIteration];

      for(unsigned int i = nbIteration ; i < A.rows(); ++i)
        for(unsigned int j = nbIteration ; j < A.columns(); ++j) {
          T tmp = T(fabs((A[i][j])));
          if((tmp > maxValue)){
            maxValueRow = i;
            maxValueCol = j;
            maxValue = tmp;
          }
        }
  }


  template <class T>
  void determinantRowsElimination(Matrix<T> & A, const unsigned int nbIteration){
      // for every remaining rows
      for(unsigned int i = nbIteration+1; i < A.rows(); ++i){
        T tmp = T(A[i][nbIteration] / A[nbIteration][nbIteration]);

        // row_i = row_i – Aik/Akk row_k
        A[i][nbIteration] = T(0);
        for(unsigned int j = nbIteration+1; j < A.columns(); ++j){
            A[i][j] -=  tmp * A[nbIteration][j];
        }
      }
  }


  inline int determinant(const Matrix<int> & A_){
      throw MathException("invalid input matrix type, convert your matrix in float or double" ,"determinant");
  }


  /*
   * End of Namespace
   */
}

/*
 * End of Anti-doublon
 */
#endif






