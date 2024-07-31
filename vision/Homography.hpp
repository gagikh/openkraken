/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team   
 * Homography.hpp created in 12 2008.
 * Mail : biri@univ-mlv.fr
 *   
 * This file is part of the OpenKraken-vision.
 *
 * The OpenKraken-vision is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * The OpenKraken-vision is distributed in the hope that it will be useful,
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
#ifndef __OPENKN_VISION__HOMOGRAPHY_HPP__
#define __OPENKN_VISION__HOMOGRAPHY_HPP__

/*
 * External Includes
 */
#include <vector>

/*
 * Internal Includes
 */
#include <OpenKraken/math/Matrix3x3.hpp>
#include <OpenKraken/math/Vector3.hpp>

/*
 * Namespace
 */
namespace kn{

            kn::Matrix3x3d computeHomographyNormalized(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d);

            kn::Matrix3x3d computeHomographyUnNormalized(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d);

            kn::Matrix3x3d computeHomography(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d, 
                                             const bool& normalized = true);
  
          /*
           * End of Namespace
           */
}

/*
 * End of Anti-doublon
 */
#endif
