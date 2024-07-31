/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team   
 * EpipolarGeometry.cpp created in 12 2008.
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
#include "EpipolarGeometry.hpp"
#include "VisionException.hpp"
#include <OpenKraken/math/Vector.hpp>
#include <OpenKraken/math/Vector2.hpp>
#include <OpenKraken/math/SVD.hpp>
#include <OpenKraken/math/Solver.hpp>
#include <OpenKraken/math/InverseMatrix.hpp>

/*
 * Namespace
 */
namespace kn{

        // normalized DLT (Direct Linera Transformation)
        // cf multiple view geometry 1st edition p92
        static void normalize(std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d,
                              kn::Matrix3x3d& T1,
                              kn::Matrix3x3d& T2){

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

                // normalize the input data
                it = points2d2d.begin();
                while(it != points2d2d.end()){
                        it->first = T1 * it->first;
                        it->second = T2 * it->second;
                        ++it;
                }
        }


        // cf multiple view geometry 1st edition pp265
        static kn::Matrixd createSystem(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d)
        {
                kn::Matrixd A(points2d2d.size(), 9); 
                std::vector< std::pair<kn::Vector3d,kn::Vector3d> >::const_iterator it = points2d2d.begin();
                unsigned int i = 0;
                while(it != points2d2d.end()){
                        A[i][0] = it->second[0]*it->first[0];
                        A[i][1] = it->second[0]*it->first[1];
                        A[i][2] = it->second[0]*it->first[2];
                
                        A[i][3] = it->second[1]*it->first[0];
                        A[i][4] = it->second[1]*it->first[1];  
                        A[i][5] = it->second[1]*it->first[2];
                
                        A[i][6] = it->second[2]*it->first[0];
                        A[i][7] = it->second[2]*it->first[1];
                        A[i][8] = it->second[2]*it->first[2];

                        // next line
                        ++it;
                        ++i;
                }
                return A;
        }


       /*
        * Compute the fundamental matrix F
        */
        kn::Matrix3x3d computeFundamentalMatrix(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d, 
                                                const bool& normalized){
                if(normalized)
                    return computeFundamentalMatrixNormalized(points2d2d);
                else
                    return computeFundamentalMatrixUnNormalized(points2d2d);
            }


       /*
        * Compute the fundamental matrix F
        */
        kn::Matrix3x3d computeFundamentalMatrixNormalized(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d){
                if(points2d2d.size() < 8){
                  throw VisionException("at least 8 pixel correspondance required","computeFundamentalMatrixNormalized");
                }

                // data normalization
                kn::Matrix3x3d T1(0.0);
                kn::Matrix3x3d T2(0.0);
                std::vector< std::pair<kn::Vector3d,kn::Vector3d> > ptmp (points2d2d);
                normalize(ptmp,T1,T2);

                // F computation
                kn::Matrixd A = createSystem(ptmp);
                kn::Vectord f(9);
                kn::solveNullSystemSVD(A,f);
                kn::Matrix3x3d F(f);

                // enforcement : det F = 0
                enforcementFundamentalMatrix(F);

                // data denormalization
                F = T2.getTranspose() * F * T1;
                
                return F;
            }


           /*
            * Compute the fundamental matrix F
            */
            kn::Matrix3x3d computeFundamentalMatrixUnNormalized(const std::vector< std::pair<kn::Vector3d,kn::Vector3d> >& points2d2d){
                if(points2d2d.size() < 8){
                  throw VisionException("at least 8 pixel correspondance required","computeFundamentalMatrixNormalized");
                }

                // F computation
                kn::Matrixd A = createSystem(points2d2d);
                kn::Vectord f(9);
                kn::solveNullSystemSVD(A,f);
                kn::Matrix3x3d F(f);

                // enforcement : det F = 0
                enforcementFundamentalMatrix(F);

                return F;
            }


           /*
            * force the fundamental matrix F to be rank 2
            */
            void enforcementFundamentalMatrix(kn::Matrix3x3d &F){
                kn::Matrixd U(F);
                kn::Matrixd V(3,3);
                kn::Vectord D(3);

                // SVD decompostion
                kn::decompositionSVD(U,D,V);
                kn::sortSingularValuesSVD(U,D,V);

                // V = DV with D[2][2] = 0
                V.transpose();

                V[0][0] *= D[0];
                V[0][1] *= D[0];
                V[0][2] *= D[0];

                V[1][0] *= D[1];
                V[1][1] *= D[1];
                V[1][2] *= D[1];

                V[2][0] = V[2][1] = V[2][2] = 0.0;

                // recompose F
                F = U*V;
            }


           /*
            * compute the fundamental matrix between 2 projective cameras :
            * F  = [e2]xP2.P1+
            * e2 = P2C
            */
            kn::Matrix3x3d computeFundamentalMatrix(const kn::ProjectiveCamera &cam1, const kn::ProjectiveCamera &cam2){
                // e2
                kn::Vector3d e2(cam2.P()*cam1.C());

                // P1+
                kn::Matrix<double> P1pseudoInverse(kn::pseudoInverseMatrixSVD(cam1.P()));

                // [e2]x
                kn::Matrix3x3d eCross;
                eCross.cross3x3(e2);
    
                // F = [e2]xP2.P1+
                return kn::Matrix3x3d(eCross*cam2.P()*P1pseudoInverse);
            }


           /*
            * \brief compute 2 epipoles from a fundamental matrix
            * Solve F.e1=0 and Ft.e2=0
            */
            void computeEpipoles(const kn::Matrix3x3d &F, kn::Vector3d &e1, kn::Vector3d &e2){
                // epipole 1 : F.e1 = 0
                kn::solveNullSystemSVD(F,e1);
                e1.setHomogeneousNormalForm();

                // epipole 2 : Ft.e2 = 0
                kn::Matrix3x3d Ftranspose(F.getTranspose());
                kn::solveNullSystemSVD(Ftranspose,e2);          
                e2.setHomogeneousNormalForm();
            }



        /*
         * End of Namespace
         */
}
