/*************************************************************************** \
 * Copyright (C) by University Paris-Est - MISS team
 * ProjectiveCamera.hpp created in 01 2009.
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
#ifndef __OPENKN_VISION__PROJECTIVE_CAMERA_HPP__
#define __OPENKN_VISION__PROJECTIVE_CAMERA_HPP__

/*
 * External Includes
 */
#include <vector>


/*
 * Internal Includes
 */
#include <OpenKraken/math/Matrix.hpp>
#include <OpenKraken/math/Matrix3x3.hpp>
#include <OpenKraken/math/Matrix4x4.hpp>
#include <OpenKraken/math/Vector.hpp>
#include <OpenKraken/math/Vector3.hpp>

#ifdef INTERFACEDLL_EXPORTS
#define UserApp_API __declspec(dllexport)
#else
#define UserApp_API __declspec(dllimport)
#endif

/*
 * Namespace
 */
namespace kn{

        /*
         * Class definition
         */

        class UserApp_API ProjectiveCamera {


                /*
                 * class variable
                 */
        protected :
                kn::Matrix<double> matP;

                kn::Matrix3x3<double> matK;

                kn::Matrix3x3<double> matKinverse;

                kn::Matrix3x3<double> matR;

                kn::Vector4<double> vecC;


        public:

                ProjectiveCamera(const ProjectiveCamera &myProjectiveCamera);

                ProjectiveCamera(const kn::Matrix<double>& myP);

                ProjectiveCamera(const kn::Matrix3x3<double> &myK,
                                 const kn::Matrix3x3<double> &myR,
                                 const kn::Vector4<double> &myC);

                ProjectiveCamera(const std::vector<std::pair<kn::Vector4<double>,kn::Vector3<double> > > & myList3d2d);

                ProjectiveCamera(const kn::Matrix3x3<double> &myK);

                ProjectiveCamera(const std::vector<kn::Matrix3x3d> &Hlist);

                ~ProjectiveCamera();


        protected :
                void decompose();

                kn::Matrix<double> compose(const kn::Matrix3x3<double> &K, 
                                           const kn::Matrix3x3<double> &R, 
                                           const kn::Vector4<double> &C) const;

                void extractCentre();

                void computeP(const std::vector< std::pair<kn::Vector4<double>,kn::Vector3<double> > > & myList3d2d);

                void normalizationP(std::vector< std::pair<kn::Vector4<double>,kn::Vector3<double> > > &myList, 
                                    kn::Matrix3x3<double> &T, 
                                    kn::Matrix3x3<double> &Tinverse, 
                                    kn::Matrix4x4<double> &U) const; 

                void makeSystemDLT(kn::Matrix<double> &A, 
                                   const std::vector<std::pair<kn::Vector4<double>,kn::Vector3<double> > > &myList) const;

                void rotationEnforce(kn::Matrix<double> &R)const;

                kn::Matrix3x3d computeInternalParameters(const std::vector<kn::Matrix3x3d> &Hlist);

                kn::Matrix3x3<double> composeInternalParameterZhang(const kn::Vector<double> &b);

                kn::Vector<double> vZhang(const kn::Matrix<double> &V,
                                          const int i,
                                          const int j)const;

                void makeSystemZhang(kn::Matrix<double> &V,
                                     const std::vector<kn::Matrix3x3d> &Hlist,
                                     const bool zeroSkew = true) const;


        public :
                void computeExternalParameters(const kn::Matrix3x3<double> &H);

                inline kn::Matrix<double> P() const{
                        return matP;
                };

                inline kn::Matrix3x3<double> K() const{
                        return matK;
                };

                inline kn::Matrix3x3<double> R() const{
                        return matR;
                };

                inline kn::Vector4<double> C() const{
                        return vecC;
                };

                kn::Matrix<double> M() const;


                inline kn::Vector<double> project(const kn::Vector4d &X) const{
                       kn::Vector<double> x(matP*X);
                       x.setHomogeneousNormalForm();
                       return x; 
                };

                inline void updateRC(const kn::Matrix3x3d &R, const kn::Vector4d &C){
                        matR = R;
                        vecC = C;
                        vecC.setHomogeneousNormalForm();
                        matP = compose(matK, matR, vecC);
                };

                 void updateM(const kn::Matrixd &M);

                void updateP(const kn::Matrixd& P);

                void getGLProjectionMatrix(const size_t& cameraWidth, 
                                           const size_t& cameraHeight,
                                           const float& glnear,
                                           const float& glfar,
                                           float GLPMat[16]) const;


                void getGLModelviewMatrix(float GLMMat[16]) const;


                kn::Vector3d principalRay();


                static void weakCalibration(kn::ProjectiveCamera &cam1, kn::ProjectiveCamera &cam2, const kn::Matrix3x3d &F);

        };

        /*
         * End of Namespace
         */
}

/*
 * End of Anti-doublon
 */
#endif
