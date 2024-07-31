/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team   
 * Homography.cpp created in 12 2008.
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
 * Internal Includes
 */
#include "Homography.hpp"
#include "VisionException.hpp"
#include <OpenKraken/math/Solver.hpp>
#include <OpenKraken/math/Vector2.hpp>

/*
 * Namespace
 */
namespace kn{

        // normalized DLT (Direct Linera Transformation)
        // cf multiple view geometry 1st edition p92
        static void normalize(std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d,
                              kn::Matrix3x3d& T1,
                              kn::Matrix3x3d& T2,
                              kn::Matrix3x3d& T2Inverse)
        {

                // Compute the average for each set of points
                kn::Vector2d averagePosList1(0.0);
                kn::Vector2d averagePosList2(0.0);
                std::vector< std::pair<kn::Vector3d,kn::Vector3d> >::iterator it = points2d2d.begin();
                while(it != points2d2d.end()){
                        averagePosList1 += it->first.getUnhomogeneous();
                        averagePosList2 += it->second.getUnhomogeneous();
                        ++it;
                }

                averagePosList1 /= double(points2d2d.size());
                averagePosList2 /= double(points2d2d.size());

                // Compute the mean for each set of points
                double meanList1 = 0.0;
                double meanList2 = 0.0;
                kn::Vector2d tmp;
                it = points2d2d.begin();
                while(it != points2d2d.end()){
                        tmp = it->first.getUnhomogeneous() - averagePosList1;
                        meanList1 += tmp.getNorm();
                        tmp = it->second.getUnhomogeneous() - averagePosList2;
                        meanList2 += tmp.getNorm();
                        ++it;
                }
                meanList1 /= double(points2d2d.size());
                meanList2 /= double(points2d2d.size());

                double scaleList1 = meanList1 / 1.4142135623730950488;
                double scaleList2 = meanList2 / 1.4142135623730950488;

                // T matrix
                T1.setIdentity();
                T1[0][0] = 1.0/scaleList1;
                T1[1][1] = 1.0/scaleList1;
                T1[0][2] = -averagePosList1.x()/scaleList1;
                T1[1][2] = -averagePosList1.y()/scaleList1;

                // T' matrix
                T2.setIdentity();
                T2[0][0] = 1.0/scaleList2;
                T2[1][1] = 1.0/scaleList2;
                T2[0][2] = -averagePosList2.x()/scaleList2;
                T2[1][2] = -averagePosList2.y()/scaleList2;

                // T'^-1
                T2Inverse.setIdentity();
                T2Inverse[0][0] = scaleList2;
                T2Inverse[1][1] = scaleList2;
                T2Inverse[0][2] = averagePosList2.x();
                T2Inverse[1][2] = averagePosList2.y();

                // normalize the input data
                it = points2d2d.begin();
                while(it != points2d2d.end()){
                        it->first = T1 * it->first;
                        it->second = T2 * it->second;
                        ++it;
                }
        }


        // cf multiple view geometry 1st edition pp72
        static kn::Matrixd createSystem(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d)
        {
                kn::Matrixd m(points2d2d.size()*2, 9); 
                std::vector< std::pair<kn::Vector3d,kn::Vector3d> >::const_iterator it = points2d2d.begin();
                unsigned int i = 0;
                while(it != points2d2d.end()){
                        // first line
                        m[2*i][0] = 0.0;
                        m[2*i][1] = 0.0;
                        m[2*i][2] = 0.0;
                
                        m[2*i][3] = -it->second.z()*it->first.x();
                        m[2*i][4] = -it->second.z()*it->first.y();  
                        m[2*i][5] = -it->second.z()*it->first.z();
                
                        m[2*i][6] = it->second.y()*it->first.x();
                        m[2*i][7] = it->second.y()*it->first.y();
                        m[2*i][8] = it->second.y()*it->first.z();

                        // second line
                        m[2*i+1][0] = it->second.z()*it->first.x();
                        m[2*i+1][1] = it->second.z()*it->first.y();
                        m[2*i+1][2] = it->second.z()*it->first.z();

                        m[2*i+1][3] = 0.0;
                        m[2*i+1][4] = 0.0;
                        m[2*i+1][5] = 0.0;

                        m[2*i+1][6] = -it->second.x()*it->first.x();
                        m[2*i+1][7] = -it->second.x()*it->first.y(); 
                        m[2*i+1][8] = -it->second.x()*it->first.z();

                        ++it;
                        ++i;
                }
                return m;
        }


        kn::Matrix3x3d computeHomographyNormalized(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d)
        {
                if(points2d2d.size() < 4){
                  throw VisionException("at least 4 pixel correspondance required","computeHomographyNormalized");
                }

                // data normalization
                std::vector< std::pair<kn::Vector3d,kn::Vector3d> > ptmp (points2d2d);
                kn::Matrix3x3d T1(0.0);
                kn::Matrix3x3d T2(0.0);
                kn::Matrix3x3d T2Inverse(0.0);
                normalize(ptmp,T1,T2,T2Inverse);
                
                // H computation
                kn::Matrixd M = createSystem(ptmp);
                kn::Vectord h(9);
                kn::solveNullSystemSVD(M,h);
                kn::Matrix3x3d H(h);

                // data denormalization
                H = T2Inverse * H * T1;
                
                return H;
            }


        kn::Matrix3x3d computeHomographyUnNormalized(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d)
        {
                if(points2d2d.size() < 4){
                    throw VisionException("at least 4 pixel correspondance required","computeHomographyUnNormalized");
                }

                // H computation
                kn::Matrixd A = createSystem(points2d2d);
                kn::Vectord h(9);
                kn::solveNullSystemSVD(A,h);

                return kn::Matrix3x3d(h);
            }



        kn::Matrix3x3d computeHomography(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d, const bool& normalized)
        {
                if(normalized)
                    return computeHomographyNormalized(points2d2d);
                else
                    return computeHomographyUnNormalized(points2d2d);
        }



        /*
         * End of Namespace
         */
}
