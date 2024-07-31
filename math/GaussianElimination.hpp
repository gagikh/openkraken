/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * GaussianElimination.hpp created in 10 2009.
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
#ifndef __OPENKN_GAUSSIAN_ELIMINATION_HPP__
#define __OPENKN_GAUSSIAN_ELIMINATION_HPP__

/*
 * External Includes
 */
#include <algorithm>
#include <vector>

/*
 * Internal Includes
 */
#include "MathException.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"


namespace kn{

  template <class T>
  unsigned int gaussianElimination(Matrix<T> & A, Vector<T> &b, Vector<T> &x, const bool total = false){
    if(total) return gaussianEliminationTotal(A,b,x);
    else return gaussianEliminationPartial(A,b,x);
  }


  template <class T>
  unsigned int gaussianEliminationPartial(Matrix<T> & A, Vector<T> &b, Vector<T> &x){
    // thresold to decide whether a pivot is zero or not
    const double gaussianEliminationZero = 1.0e-10;

    // check the size compatibility
    if(A.columns() != A.rows())
      throw MathException("The input matrix is not a square matrix" ,"gaussianEliminationPartial");

    if(A.columns() != b.size())
      throw MathException("The input matrix is not compatible with the 'b' vector" ,"gaussianEliminationPartial");

    if(A.rows() != x.size())
      throw MathException("The input matrix is not compatible with the 'x' vector" ,"gaussianEliminationPartial");

    // make A triangular
    unsigned int nbIteration = 0;
    while(nbIteration < A.columns()-1){

      // find the pivot
      T maxValue;
      unsigned int maxValueRow;
      gaussianEliminationFindPartialPivot(A,nbIteration,maxValueRow,maxValue);

      // check if the "pivot" is non-null
      if(fabs(maxValue) < gaussianEliminationZero)
        return 0; // the matrix is singular (not full rank)

      // if required, swap rows
      if(maxValueRow != nbIteration){
        // swap the non zero part of the rows
        std::swap_ranges(&(A[maxValueRow][nbIteration]),A[maxValueRow]+A.columns(),&(A[nbIteration][nbIteration]));
        // swap the corresponding elements in vector b
        std::swap(b[maxValueRow],b[nbIteration]);
      }

      // check the last pivot
      if(fabs(A[nbIteration][nbIteration]) < gaussianEliminationZero) return 0;

      // row alimination (let's put some 0)
      gaussianRowsElimination(A,b,nbIteration);

      // next iteration
      ++nbIteration;
    }

    // check the last pivot
    if(fabs(A[nbIteration][nbIteration]) < gaussianEliminationZero) return 0;

    // elimination (find the solution) back substitution
    gaussianBackSubstitutionPartial(A,b,x);

    return 1; // full rank matrix
  }


  template <class T>
  unsigned int gaussianEliminationTotal(Matrix<T> & A, Vector<T> &b, Vector<T> &x){
    // thresold to decide whether a pivot is zero or not
    const double gaussianEliminationZero = 1.0e-10;

    // check the size compatibility
    if(A.columns() != A.rows())
      throw MathException("The input matrix is not a square matrix" ,"gaussianEliminationPartial");

    if(A.columns() != b.size())
      throw MathException("The input matrix is not compatible with the 'b' vector" ,"gaussianEliminationPartial");

    if(A.rows() != x.size())
      throw MathException("The input matrix is not compatible with the 'x' vector" ,"gaussianEliminationPartial");

    // sort the permutations
    std::vector<unsigned int> sorter(x.size());
    for(unsigned int i=0; i<sorter.size(); ++i)
      sorter[i] = i;

    // make A triangular
    unsigned int nbIteration = 0;
    while(nbIteration < A.columns()-1){

      // find the pivot
      T maxValue;
      unsigned int maxValueRow,maxValueCol;
      gaussianEliminationFindTotalPivot(A,nbIteration,maxValueRow,maxValueCol,maxValue);

      // check if the "pivot" is non-null
      if(fabs(maxValue) < gaussianEliminationZero)
        return 0; // the matrix is singular (not full rank). Note that rank(A) = nbIteration

      // if required, swap rows
      if(maxValueRow != nbIteration){
        // swap the non zero part of the rows
        std::swap_ranges(&(A[maxValueRow][nbIteration]),A[maxValueRow]+A.columns(),&(A[nbIteration][nbIteration]));
        // swap the corresponding elements in vector b
        std::swap(b[maxValueRow],b[nbIteration]);
      }

      // if required, swap columns
      if(maxValueCol != nbIteration){
        // swap columns
        A.swapColumns(maxValueCol,nbIteration);
        // swap the corresponding elements in vector x
        std::swap(sorter[maxValueCol],sorter[nbIteration]);
      }

      // row alimination (let's put some 0)
      gaussianRowsElimination(A,b,nbIteration);

      // next iteration
      ++nbIteration;
    }

    // check the last pivot
    if(fabs(A[nbIteration][nbIteration]) < gaussianEliminationZero) return 0;

    // elimination (find the solution) back substitution
    gaussianBackSubstitutionTotal(A,b,x,sorter);

    return 1; // full rank matrix
  }


  template <class T>
  unsigned int gaussianEliminationInverseMatrix(const Matrix<T> & A, Matrix<T> & Ainverse, const bool total = false){
    if(total) return gaussianEliminationInverseMatrixTotal(A,Ainverse);
    else return gaussianEliminationInverseMatrixPartial(A,Ainverse);
  }


  template <class T>
  unsigned int gaussianEliminationInverseMatrixPartial(const Matrix<T> & A_, Matrix<T> & Ainverse){
    // thresold to decide whether a pivot is zero or not
    const double gaussianEliminationZero = 1.0e-10;

    // make a copy of A_
    Matrix<T> A(A_);

    // check the size compatibility
    if(A.columns() != A.rows())
      throw MathException("The input matrix is not a square matrix" ,"gaussianEliminationInverseMatrixPartial");

    if(Ainverse.columns() != Ainverse.rows())
      throw MathException("The output matrix is not compatible with the 'b' vector" ,"gaussianEliminationInverseMatrixPartial");

    if(A.rows() != Ainverse.rows())
      throw MathException("The input an d output matrices size are not compatible" ,"gaussianEliminationInverseMatrixPartial");

    // identity matrix
    Matrix<T> Id(A.rows(),A.columns());
    Id.setIdentity();

    // make A triangular
    unsigned int nbIteration = 0;
    while(nbIteration < A.columns()-1){

      // find the pivot
      T maxValue;
      unsigned int maxValueRow;
      gaussianEliminationFindPartialPivot(A,nbIteration,maxValueRow,maxValue);

      // check if the "pivot" is non-null
      if(fabs(maxValue) < gaussianEliminationZero)
        return 0; // the matrix is singular (not full rank), we cannot inverse this matrix

      // if required, swap rows
      if(maxValueRow != nbIteration){
        // swap the non zero part of the rows
        std::swap_ranges(&(A[maxValueRow][nbIteration]),A[maxValueRow]+A.columns(),&(A[nbIteration][nbIteration]));
        // swap the corresponding rows in Id matrix
        Id.swapRows(maxValueRow,nbIteration);
      }

      // row alimination (let's put some 0)
      gaussianRowsElimination(A,Id,nbIteration);

      // next iteration
      ++nbIteration;
    }

    // check the last pivot
    if(fabs(A[nbIteration][nbIteration]) < gaussianEliminationZero) return 0;

    // elimination (find the solution) back substitution
    gaussianBackSubstitutionPartial(A,Id,Ainverse);

    return 1; // full rank matrix
  }


  template <class T>
  unsigned int gaussianEliminationInverseMatrixTotal(const Matrix<T> & A_, Matrix<T> & Ainverse){
    // thresold to decide whether a pivot is zero or not
    const double gaussianEliminationZero = 1.0e-10;

    // make a copy of A_
    Matrix<T> A(A_);

    // check the size compatibility
    if(A.columns() != A.rows())
      throw MathException("The input matrix is not a square matrix" ,"gaussianEliminationInverseMatrixPartial");

    if(Ainverse.columns() != Ainverse.rows())
      throw MathException("The output matrix is not compatible with the 'b' vector" ,"gaussianEliminationInverseMatrixPartial");

    if(A.rows() != Ainverse.rows())
      throw MathException("The input an d output matrices size are not compatible" ,"gaussianEliminationInverseMatrixPartial");

    // identity matrix
    Matrix<T> Id(A.rows(),A.columns());
    Id.setIdentity();

    // sort the permutations
    std::vector<unsigned int> sorter(Id.rows());
    for(unsigned int i=0; i<sorter.size(); ++i)
      sorter[i] = i;

    // make A triangular
    unsigned int nbIteration = 0;
    while(nbIteration < A.columns()-1){

      // find the pivot
      T maxValue;
      unsigned int maxValueRow,maxValueCol;
      gaussianEliminationFindTotalPivot(A,nbIteration,maxValueRow,maxValueCol,maxValue);

      // check if the "pivot" is non-null
      if(fabs(maxValue) < gaussianEliminationZero)
        return 0; // the matrix is singular (not full rank). Note that rank(A) = nbIteration

      // if required, swap rows
      if(maxValueRow != nbIteration){
        // swap the non zero part of the rows
        std::swap_ranges(&(A[maxValueRow][nbIteration]),A[maxValueRow]+A.columns(),&(A[nbIteration][nbIteration]));
        // swap the corresponding rows in Id matrix
        Id.swapRows(maxValueRow,nbIteration);
      }

      // if required, swap columns
      if(maxValueCol != nbIteration){
        // swap columns
        A.swapColumns(maxValueCol,nbIteration);
        // swap the corresponding elements in vector x
        std::swap(sorter[maxValueCol],sorter[nbIteration]);
      }

      // row alimination (let's put some 0)
      gaussianRowsElimination(A,Id,nbIteration);

      // next iteration
      ++nbIteration;
    }

    // check the last pivot
    if(fabs(A[nbIteration][nbIteration]) < gaussianEliminationZero) return 0;

    // elimination (find the solution) back substitution
    gaussianBackSubstitutionTotal(A,Id,Ainverse,sorter);

    return 1; // full rank matrix
  }


  template <class T>
  void gaussianBackSubstitutionPartial(const Matrix<T> & A, const Vector<T> &b, Vector<T> &x){
    for(int i=int(A.rows())-1; i>=0; --i){
      T numerator = T(0);
      for(int j=i+1; j<int(A.columns()); ++j)
        numerator += A[i][j] * x[j];

      x[i] = (b[i] - numerator) / A[i][i];
    }
  }


  template <class T>
  void gaussianBackSubstitutionPartial(const Matrix<T> & A, const Matrix<T> &B, Matrix<T> &X){
    for(unsigned int m=0; m<X.columns(); ++m)
      for(int i=int(A.rows())-1; i>=0; --i){
        T numerator = T(0);
        for(int j=i+1; j<int(A.columns()); ++j)
          numerator += A[i][j] * X[j][m];

        X[i][m] = (B[i][m] - numerator) / A[i][i];
      }
  }


  template <class T>
  void gaussianBackSubstitutionTotal(const Matrix<T> & A, const Vector<T> &b, Vector<T> &x, const std::vector<unsigned int> &sorter){
    for(int i=int(A.rows())-1; i>=0; --i){
      T numerator = T(0);
      for(int j=i+1; j<int(A.columns()); ++j)
        numerator += A[i][j] * x[sorter[j]];

      x[sorter[i]] = (b[i] - numerator) / A[i][i];
    }
  }


  template <class T>
  void gaussianBackSubstitutionTotal(const Matrix<T> & A, const Matrix<T> &B, Matrix<T> &X, const std::vector<unsigned int> &sorter){
    for(unsigned int m=0; m<X.columns(); ++m)
      for(int i=int(A.rows())-1; i>=0; --i){
        T numerator = T(0);
        for(int j=i+1; j<int(A.columns()); ++j)
          numerator += A[i][j] * X[sorter[j]][m]; 

        X[sorter[i]][m] = (B[i][m] - numerator) / A[i][i];
      }
  }


  template <class T>
  void gaussianEliminationFindPartialPivot(const Matrix<T> & A, const unsigned int nbIteration, 
                                           unsigned int & maxValueRow, T & maxValue){

      maxValueRow = nbIteration;
      maxValue = A[nbIteration][nbIteration];

      for(unsigned int i = nbIteration+1 ; i < A.rows(); ++i) {
        T tmp = T(fabs((A[i][nbIteration])));
        if((tmp > maxValue)){
          maxValueRow = i;
          maxValue = tmp;
        }
      }
  }


  template <class T>
  void gaussianEliminationFindTotalPivot(const Matrix<T> & A, const unsigned int nbIteration, 
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
  void gaussianRowsElimination(Matrix<T> & A, Vector<T> & b, const unsigned int nbIteration){
      // for every remaining rows
      for(unsigned int i = nbIteration+1; i < A.rows(); ++i){ 
        T tmp = T(A[i][nbIteration] / A[nbIteration][nbIteration]);

        // row_i = row_i – Aik/Akk row_k
        A[i][nbIteration] = T(0);
        for(unsigned int j = nbIteration+1; j < A.columns(); ++j){ 
            A[i][j] -=  tmp * A[nbIteration][j];
        }

        // the same for b
        b[i] -= tmp * b[nbIteration]; 
      }
  }


  template <class T>
  void gaussianRowsElimination(Matrix<T> & A, Matrix<T> & Id, const unsigned int nbIteration){
      // for every remaining rows
      for(unsigned int i = nbIteration+1; i < A.rows(); ++i){ 
        T tmp = T(A[i][nbIteration] / A[nbIteration][nbIteration]);

        // row_i = row_i – Aik/Akk row_k
        A[i][nbIteration] = T(0);
        for(unsigned int j = nbIteration+1; j < A.columns(); ++j){ 
            A[i][j] -=  tmp * A[nbIteration][j];
        }

        // the same for Id
        for(unsigned int j = 0; j < Id.columns(); ++j){ 
            Id[i][j] -=  tmp * Id[nbIteration][j];
        }
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