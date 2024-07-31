/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * SVD.hpp created in 10 2008.
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
#ifndef __OPENKN_MATH__SVD_HPP__
#define __OPENKN_MATH__SVD_HPP__

#include "Matrix.hpp"

namespace kn{

  void decompositionSVD(Matrix<double>& a,
      Vector<double>& w,
      Matrix<double>& v);


  void sortSingularValuesSVD(Matrix<double>& u,
      Vector<double>& w,
      Matrix<double>& v);

  void solveSVD(const Matrix<double>& u,
      const Vector<double>& w,
      const Matrix<double>& v,
      const Vector<double>& b,
      Vector<double>& x);

}
#endif