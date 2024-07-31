/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * InverseMatrix.cpp created in 11 2008.
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


#include "InverseMatrix.hpp"

namespace kn{

  Matrix<double> inverseMatrixSVD(const Matrix<double>& m){

    if(m.rows() != m.columns())
      throw MathException("The input matrix is not a square matrix" ,"inverseMatrixSVD");


    Matrix<double> u(m);
    Matrix<double> v(m.columns(),m.columns());
    Vector<double> d(m.columns());

    decompositionSVD(u,d,v);

    u.transpose();

    const double delta = 1.0e-12;
    for(unsigned int i=0; i<d.size(); i++){
      if(std::fabs(d[i]) > delta)
        d[i] = 1.0/d[i];
      else
        d[i] = 0.0;
    }

    Matrix<double> dinverse(d.size(),d.size());
    dinverse.setZero();
    dinverse.setDiagonal(d);

    return v*dinverse*u;
  }

  Matrix<double> pseudoInverseMatrixSVD(const Matrix<double>& m){

    Matrix<double> mtrans(m);
    mtrans.transpose();

    Matrix<double> inverse(inverseMatrixSVD(mtrans*m));

    return inverse*mtrans;
  }


  Matrix<double> inverseMatrixGaussianElimination(const Matrix<double>& m, const bool total){
    kn::Matrix<double> inverse(m.rows(), m.columns());
    gaussianEliminationInverseMatrix(m,inverse,total);
    return  inverse;
  }

  Matrix<double> pseudoInverseMatrixGaussianElimination(const Matrix<double>& m, const bool total){
    Matrix<double> mTranspose(m);
    mTranspose.transpose();

    Matrix<double> inverse(inverseMatrixGaussianElimination(mTranspose*m, total));

    return inverse*mTranspose;
  }

}
