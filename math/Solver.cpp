/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * Solver.cpp created in 11 2008.
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

#include <sstream>


#include "Solver.hpp"
#include "GaussianElimination.hpp"

namespace kn{

  void solveSystemSVD(const Matrix<double>& a,
                      Vector<double>& x,
                      const Vector<double>& b){

    Vector<double> w(a.columns());
    Matrix<double> v(a.columns(),a.columns());
    Matrix<double> atmp(a);
    decompositionSVD(atmp,w,v);
    solveSVD(atmp,w,v,b,x);
  }


  void solveNullSystemSVD(const Matrix<double>& a,
                          Vector<double>& x){

    Vector<double> w(a.columns());
    Matrix<double> v(a.columns(),a.columns());
    Matrix<double> atmp(a);
    decompositionSVD(atmp,w,v);
    sortSingularValuesSVD(atmp,w,v);

    x = v.getColumn(v.columns()-1);
    x.setHomogeneousNormalForm();
  }


  void solveSystemGaussianElimination(const Matrix<double>& a,
                                      Vector<double>& x,
                                      const Vector<double>& b, 
                                      const bool total){
    kn::Matrix<double> aInv(a.rows(), a.columns());
    gaussianEliminationInverseMatrix(a,aInv,total);

    x = aInv * b;
  }


}
