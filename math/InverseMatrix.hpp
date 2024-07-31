/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * InverseMatrix.hpp created in 11 2008.
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
#ifndef __OPENKN_MATH__INVERSE_MATRIX_HPP__
#define __OPENKN_MATH__INVERSE_MATRIX_HPP__

#include "SVD.hpp"
#include "GaussianElimination.hpp"

namespace kn{

  Matrix<double> inverseMatrixSVD(const Matrix<double>& m);

  Matrix<double> pseudoInverseMatrixSVD(const Matrix<double>& m);

  Matrix<double> inverseMatrixGaussianElimination(const Matrix<double>& m, const bool total = false);

  Matrix<double> pseudoInverseMatrixGaussianElimination(const Matrix<double>& m, const bool total = false);


}

#endif