/*************************************************************************** \
 * Copyright (C) by University Paris-Est - MISS team
 * ProjectiveCamera.cpp created in 01 2009.
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
 * External Includes
 */
#include <cstring>



/*
 * Internal Includes
 */
#include "ProjectiveCamera.hpp"
#include "VisionException.hpp"
#include "Homography.hpp"
#include "EpipolarGeometry.hpp"
#include <OpenKraken/math/SVD.hpp>
#include <OpenKraken/math/Solver.hpp>
#include <OpenKraken/math/InverseMatrix.hpp>
#include <OpenKraken/math/Determinant.hpp>
#include <OpenKraken/math/RQDecomposition.hpp>


/*
 * Namespace
 */
namespace kn{

        /* \brief Copy Constructor
        * \param myProjectiveCamera : Source projective camera
        */
        ProjectiveCamera::ProjectiveCamera(const ProjectiveCamera &myProjectiveCamera)
                :matP(myProjectiveCamera.matP), 
                 matK(myProjectiveCamera.matK), 
                 matKinverse(myProjectiveCamera.matKinverse),
                 matR(myProjectiveCamera.matR), 
                 vecC(myProjectiveCamera.vecC){
        }


        /*
        * \brief constructor from a projection matrix P. 
        * This function performs the following decomposition : P = K[R|-RC].
        * \param myP : a 3x4 projection matrix
        * \throw VisionException argument P matrix : invalid size
        */
        ProjectiveCamera::ProjectiveCamera(const kn::Matrix<double>& myP)
                :matP(myP){
                // check P size
                if(matP.rows() != 3 && matP.columns()!= 4)
                        throw VisionException("invalid size : the input P matrix is not a 3x4 matrix");

                // extract K, R, C from P
                decompose();

                // inverse K 
                matKinverse = kn::inverseMatrixSVD(matK);
        }


        /*
        * \brief constructor from all the required matrices : P = K[R|-RC].
        * \param myK : a 3x3 internal parameters matrix.
        * \param myR : a 3x3 rotation matrix.
        * \param myC : a 4 homogeneous vector (camera centre coordinates).
        * \throw VisionException argument : invalid size
        */
        ProjectiveCamera::ProjectiveCamera(const kn::Matrix3x3<double> &myK,
                                           const kn::Matrix3x3<double> &myR,
                                           const kn::Vector4<double>   &myC)
                :matP(3,4), matK(myK), matR(myR), vecC(myC){

                // homogeneous normal form 
                vecC.setHomogeneousNormalForm();

                // P
                matP = compose(matK, matR, vecC);

                // inverse K 
                matKinverse = kn::inverseMatrixSVD(matK);
        }


        /*
        * \brief constructor from a set of 2d-3d correspondances. 
        * Performs the following decomposition: P = K[R|-RC].
        * \param myList : a list of 3d-2d correspondances, in homogenous coordinates.
        */
        ProjectiveCamera::ProjectiveCamera(const std::vector< std::pair<kn::Vector4<double>,kn::Vector3<double> > > & myList3d2d)
                :matP(3,4){
                // find P matrix 
                computeP(myList3d2d);

                // extract K, R, C from P = K[R|-RC]
                decompose();

                // inverse K 
                matKinverse = kn::inverseMatrixSVD(matK);
        }

        /*
        * \brief constructor from an internal parameter matrix. The projection matrix will be P=K[Id|0]
        * \param myK : a 3x3 internal parameters matrix.
        * \throw VisionException argument : invalid size
        */
        ProjectiveCamera::ProjectiveCamera(const kn::Matrix3x3<double> &myK)
                :matP(3,4),matK(myK){

                // K inverse
                matKinverse = kn::inverseMatrixSVD(matK);

                // rotation
                matR.setIdentity();

                // translation
                vecC.setZero();
                vecC[3] = 1.0;

                // P
                matP = compose(matK, matR, vecC);
        }


        /*
        * \brief constructor from a list of homography : only internal parameters (K matrix) will be computed with Zhang method. The other parameters are set to identity.
        * \param Hlist : a list of 3x3 homography matrices (at least 3).
        * \throw VisionException not enought homographies
        */
        ProjectiveCamera::ProjectiveCamera(const std::vector<kn::Matrix3x3d> &Hlist)
                :matP(3,4){

                // K matrix with Zhang method
                matK = computeInternalParameters(Hlist);

                // K inverse
                matKinverse = kn::inverseMatrixSVD(matK);

                // Update external parameters with the last homography
                computeExternalParameters(Hlist[Hlist.size()-1]);
        }


        /*
        * \brief destructor
        */
        ProjectiveCamera::~ProjectiveCamera(){
        }


        /*
        * \brief compute P from K, R and C : P = K[R|-RC].
        */
        kn::Matrix<double> ProjectiveCamera::compose(const kn::Matrix3x3<double> &K, 
                                                     const kn::Matrix3x3<double> &R, 
                                                     const kn::Vector4<double>   &C) const{
                // P = K[R|-RC]
                kn::Matrix<double> P(3,4);

                // R contribution
                for(unsigned int i=0; i<3; i++)
                        for(unsigned int j=0; j<3; j++)
                                P[i][j] = R[i][j];

                // -RC contribution
                P.setColumn(3,-R*C.getUnhomogeneous());
                        
                // K contribution
                P = K * P;

                return P;
        }


        /*
        * \brief extract K, R and C from P :
        * P = [M|-MC] et M = KR from Multiple View Geometry, Hartley Zisserman, camera model
        */
        void ProjectiveCamera::decompose(){
                // M
                kn::Matrix<double> M(3,3);
                for(unsigned int i=0; i<3; i++)
                        for(unsigned int j=0; j<3; j++)
                                M[i][j] = matP[i][j];

                // K and R
                kn::rqDecomposition3x3(M,matK,matR); 
                kn::rq3x3MakePositiveDiagonal(matK,matR); 

                // we set K[3][3] to 1 (cf Multiple view geometry page 143 "finite camera")
                matK = matK / matK[2][2]; 
                matK.roundZero(1.0e-10);

                // camera centre
                extractCentre();
        }

        /*
        * \brief extract C from P
        */
        void ProjectiveCamera::extractCentre(){
                // matrix used for every component extraction
                kn::Matrix<double> det(3,3);

                // x 
                det[0][0]=matP[0][1]; det[0][1]=matP[0][2]; det[0][2]=matP[0][3];
                det[1][0]=matP[1][1]; det[1][1]=matP[1][2]; det[1][2]=matP[1][3];
                det[2][0]=matP[2][1]; det[2][1]=matP[2][2]; det[2][2]=matP[2][3];
                vecC[0] = kn::determinant3x3(det);

                // y
                det[0][0]=matP[0][0]; det[0][1]=matP[0][2]; det[0][2]=matP[0][3];
                det[1][0]=matP[1][0]; det[1][1]=matP[1][2]; det[1][2]=matP[1][3];
                det[2][0]=matP[2][0]; det[2][1]=matP[2][2]; det[2][2]=matP[2][3];
                vecC[1] = -kn::determinant3x3(det);

                // z
                det[0][0]=matP[0][0]; det[0][1]=matP[0][1]; det[0][2]=matP[0][3];
                det[1][0]=matP[1][0]; det[1][1]=matP[1][1]; det[1][2]=matP[1][3];
                det[2][0]=matP[2][0]; det[2][1]=matP[2][1]; det[2][2]=matP[2][3];
                vecC[2] = kn::determinant3x3(det);

                // w
                det[0][0]=matP[0][0]; det[0][1]=matP[0][1]; det[0][2]=matP[0][2];
                det[1][0]=matP[1][0]; det[1][1]=matP[1][1]; det[1][2]=matP[1][2];
                det[2][0]=matP[2][0]; det[2][1]=matP[2][1]; det[2][2]=matP[2][2];
                vecC[3] = -kn::determinant3x3(det);

                // homogeneous normal form 
                vecC.setHomogeneousNormalForm();
        }


        /*
        * \brief compute the P matrix from a set of 3d-2d correspondances. 
        * \param myList : a list of 3d-2d correspondances, in homogenous coordinates.
        */
        void ProjectiveCamera::computeP(const std::vector< std::pair<kn::Vector4<double>,kn::Vector3<double> > > & myList3d2d){
                // number of correspondances
                int nbElements = myList3d2d.size();
                if(nbElements < 6) 
                        throw VisionException("camera calibration : at least 6 correspondances 3d-2d are required","computeP");

                kn::Matrix<double> A(2*nbElements,12);
                kn::Vector<double> p(12);

                // data normalization
                kn::Matrix3x3<double> T;
                kn::Matrix3x3<double> Tinverse;
                kn::Matrix4x4<double> U;
                std::vector< std::pair<kn::Vector4<double>,kn::Vector3<double> > > normalizedList(myList3d2d);
                normalizationP(normalizedList,T,Tinverse,U);

                // build the linear system to solve
                makeSystemDLT(A,normalizedList);

                // solve the system
                kn::solveNullSystemSVD(A, p);

                // build P matrix from p vecrtor
                matP = kn::Matrix<double>(3,4,p,true);

                // data denormalization
                matP = Tinverse * matP * U;
        }


        /*
        * \brief data normalization when computing the P matrix from a set of 3d-2d correspondances.
        * cf. Multiple View Geometry, Hartley Zisserman, Computation of the camera matrix P
        * \param myList : a list of 3d-2d correspondances, in homogenous coordinates, to be normalized.
        * \param T : transformation matrix (will be set by the function)
        * \param Tinverse : inverse of T, to be set by the function
        * \param U : transformation matrix (will be set by the function)s
        */
        void ProjectiveCamera::normalizationP(std::vector< std::pair<kn::Vector4<double>,kn::Vector3<double> > > &myList,                                                             kn::Matrix3x3<double> &T, 
                                              kn::Matrix3x3<double> &Tinverse, 
                                              kn::Matrix4x4<double> &U) const {

                double nbPoints = (double)(myList.size());

                // list : set normal form
                std::vector< std::pair<kn::Vector4<double>,kn::Vector3<double> > >::iterator it = myList.begin();
                while(it != myList.end()){
                        it->first.setHomogeneousNormalForm();
                        it->second.setHomogeneousNormalForm();
                        ++it;
                }                               

                // Compute average for each set of points
                kn::Vector3<double> averagePosList4d(0.0);
                kn::Vector2<double> averagePosList3d(0.0);
                it = myList.begin();
                while(it != myList.end()){
                        averagePosList4d += it->first.getUnhomogeneous();
                        averagePosList3d += it->second.getUnhomogeneous();
                        ++it;
                }

                averagePosList4d /= nbPoints;
                averagePosList3d /= nbPoints;

                // compute mean
                double meanList4d = 0.0;
                double meanList3d = 0.0;
                kn::Vector3d tmp;
                it = myList.begin();
                while(it != myList.end()){
                        meanList4d += (it->first.getUnhomogeneous() - averagePosList4d).getNorm();
                        meanList3d += (it->second.getUnhomogeneous() - averagePosList3d).getNorm();
                        ++it;
                }

                meanList4d /= nbPoints;
                meanList3d /= nbPoints;

                // scale to perform
                double scaleList4d = meanList4d / 1.73205081;            // sqrt(3)
                double scaleList3d = meanList3d / 1.4142135623730950488; // sqrt(2)

                // T matrix
                T.setIdentity();
                T[0][0] = 1.0/scaleList3d;
                T[1][1] = 1.0/scaleList3d;
                T[0][2] = -averagePosList3d[0]/scaleList3d;
                T[1][2] = -averagePosList3d[1]/scaleList3d;

                // T inverse matrix
                Tinverse.setIdentity();
                Tinverse[0][0] = scaleList3d;
                Tinverse[1][1] = scaleList3d;
                Tinverse[0][2] = averagePosList3d[0];
                Tinverse[1][2] = averagePosList3d[1];

                // U matrix
                U.setIdentity();
                U[0][0] = 1.0/scaleList4d;
                U[1][1] = 1.0/scaleList4d;
                U[2][2] = 1.0/scaleList4d;
                U[0][3] = -averagePosList4d[0]/scaleList4d;
                U[1][3] = -averagePosList4d[1]/scaleList4d;
                U[2][3] = -averagePosList4d[2]/scaleList4d;

                // list normalization
                it = myList.begin();
                while(it != myList.end()){
                        it->first = U * it->first;
                        it->first.setHomogeneousNormalForm();
                        it->second = T * it->second;
                        it->second.setHomogeneousNormalForm();
                        ++it;
                }
        }


        /*
        * \brief Build the linear system to solve for the projection matrix computation. DLT : Direct Linear Transform.
        * \param A : the main matrix of the linear system. This matrix is filled in this method.
        * \param myList : a list of 3d-2d correspondances, in homogenous coordinates.
        */
        void ProjectiveCamera::makeSystemDLT(kn::Matrix<double> &A, 
                                             const std::vector< std::pair<kn::Vector4<double>,kn::Vector3<double> > > &myList) const{
                int i = 0;
                std::vector< std::pair<kn::Vector4<double>,kn::Vector3<double> > >::const_iterator iter = myList.begin();

                while(iter != myList.end()){
                        // 1st line
                        A[2*i][0]=0.0;  A[2*i][1]=0.0; 
                        A[2*i][2]=0.0;  A[2*i][3]=0.0;

                        A[2*i][4] = -iter->second[2] * iter->first[0];
                        A[2*i][5] = -iter->second[2] * iter->first[1];    // -pixel.wi * point.XiT
                        A[2*i][6] = -iter->second[2] * iter->first[2];
                        A[2*i][7] = -iter->second[2] * iter->first[3];

                        A[2*i][8] =  iter->second[1] * iter->first[0];
                        A[2*i][9] =  iter->second[1] * iter->first[1];    // pixel.yi * point.XiT
                        A[2*i][10]=  iter->second[1] * iter->first[2];
                        A[2*i][11]=  iter->second[1] * iter->first[3];

                        // 2nd line
                        A[2*i+1][0] = iter->second[2] * iter->first[0];
                        A[2*i+1][1] = iter->second[2] * iter->first[1];   // pixel.wi * point.XiT
                        A[2*i+1][2] = iter->second[2] * iter->first[2];
                        A[2*i+1][3] = iter->second[2] * iter->first[3];

                        A[2*i+1][4]=0.0;  A[2*i+1][5]=0.0;
                        A[2*i+1][6]=0.0;  A[2*i+1][7]=0.0;

                        A[2*i+1][8] = -iter->second[0] * iter->first[0];
                        A[2*i+1][9] = -iter->second[0] * iter->first[1]; // pixel.xi*point.XiT
                        A[2*i+1][10]= -iter->second[0] * iter->first[2];
                        A[2*i+1][11]= -iter->second[0] * iter->first[3];

                        // next
                        ++iter;
                        ++i;
                }
        }


        /*
        * \brief transform the matrix R into the "nearest" orthogonal matrix
        * \param R : A nearly rotation matrix.
        */
        void ProjectiveCamera::rotationEnforce(kn::Matrix<double> &R)const{
                kn::Matrix<double> V(R.rows(), R.columns());
                kn::Vector<double> D(R.columns());

                kn::decompositionSVD(R,D,V); // R = R * D * V
                V.transpose();
                R = R * V;      // we change D to identity
        }


        /*
        * \brief compute the intrinsic parameters of the camera from a list of homographies (Zhang method). In this function, we use Zhang notation.
        * \param Hlist : a list of 3x3 homography matrices (at least 3).
        * \throw VisionException not enought homographies
        */
        kn::Matrix3x3d ProjectiveCamera::computeInternalParameters(const std::vector<kn::Matrix3x3d> &Hlist){
                // check the number of homographies
                if(Hlist.size() < 3)
                        throw VisionException("not enough homographies to calibrate the camera");

                // calibration step (cf Zhang articles)
                // one additional line to handle zero skew
                kn::Matrix<double> V(2*Hlist.size()+1, 6);
                V.setZero();
                makeSystemZhang(V,Hlist,true);

                kn::Vector<double> b(6);
                kn::solveNullSystemSVD(V, b);

                // extract parameters (K matrix)
                return composeInternalParameterZhang(b);
        }


        /*
        * \brief Build the system to be solved to compute the intrinsic parameter matrix K with Zhang method.
        * \param V : the matrix to be built
        * \param Hlist : a list of homographies (at least 3)
        * \param zeroSkew : true if we want to force zero skew
        */
        void ProjectiveCamera::makeSystemZhang(kn::Matrix<double> &V,
                                               const std::vector<kn::Matrix3x3d> &Hlist,
                                               const bool zeroSkew) const {

                std::vector<kn::Matrix3x3d>::const_iterator iter = Hlist.begin();
                int i = 0;

                while(iter != Hlist.end())
                {
                        V.setRow(i,vZhang(*iter,0,1));
                        i++;
                        V.setRow(i,vZhang(*iter,0,0)-vZhang(*iter,1,1));
                        i++;
                        iter++;
                }

                // zeroSkew
                if(zeroSkew) V[V.rows()-1][1] = 1.0;
        }


        /*
        * \brief returns a Vector3d : cf Zhang article and Zhang notation
        * \param h : a matrix
        * \param i : matrix index
        * \param j : matrix index
        */
        kn::Vector<double> ProjectiveCamera::vZhang(const kn::Matrix<double> &h,
                                                    const int i,
                                                    const int j) const {
                kn::Vector<double> v(6);
                v[0] = h[0][i] * h[0][j];
                v[1] = h[0][i] * h[1][j] + h[1][i] * h[0][j];
                v[2] = h[1][i] * h[1][j];
                v[3] = h[2][i] * h[0][j] + h[0][i] * h[2][j];
                v[4] = h[2][i] * h[1][j] + h[1][i] * h[2][j];
                v[5] = h[2][i] * h[2][j];

                return v;
        }


        /*
        * \brief build the intrinsic parameter matrix K from a vector b (Zhang notation).
        * \param b : a vector to be converted into an internal parameter 3x3 matrix
        */
        kn::Matrix3x3<double> ProjectiveCamera::composeInternalParameterZhang(const kn::Vector<double> &b){
                // B = A-T.A-1
                kn::Matrix<double> B(3,3);
                B[0][0] = b[0];
                B[0][1] = B[1][0] = b[1];
                B[1][1] = b[2];
                B[0][2] = B[2][0] = b[3];
                B[1][2] = B[2][1] = b[4];
                B[2][2] = b[5];

                // parameters extraction (cf Zhang article appendix)
                kn::Matrix3x3d K;
                K.setIdentity();

                K[1][2]  = (B[1][0]*B[2][0] - B[0][0]*B[2][1]) / (B[0][0]*B[1][1] - B[1][0]*B[1][0]);
                double l = B[2][2] - (B[2][0]*B[2][0] + K[1][2]*(B[1][0]*B[2][0] - B[0][0]*B[2][1]))/B[0][0];
                K[0][0]  = l/B[0][0] > 0 ? sqrt(l/B[0][0]) : sqrt(-l/B[0][0]);
                K[1][1]  = l*B[0][0] / (B[0][0]*B[1][1] - B[1][0]*B[1][0]);
                K[1][1]  = K[1][1] > 0 ? sqrt(K[1][1]) : sqrt(-K[1][1]);
                K[0][1]  = -B[1][0] * K[0][0]*K[0][0] * K[1][1] / l;
                K[0][2]  = (K[0][1]*K[1][2] / K[1][1]) - (B[2][0] *K[0][0]*K[0][0] / l);

                return K;
        }

        /*
        * \brief compute the new position C and orientation R of the camera (Zhang method). 
        * \param H : A 3x3 homography matrix.
        * \throw VisionException argument : invalid size
        */
        void ProjectiveCamera::computeExternalParameters(const kn::Matrix3x3<double> &H){

                // scale factor and variables
                kn::Vector3d h1(H.getColumn(0));
                kn::Vector3d h2(H.getColumn(1));
                kn::Vector3d h3(H.getColumn(2));
                double l = 0.5 * (1.0 / (matKinverse*h1).getNorm()) + 0.5 * (1.0 / (matKinverse*h2).getNorm());

                // rotation matrix
                matR.setColumn(0, matKinverse * h1 * l);
                matR.setColumn(1, matKinverse * h2 * l);
                matR.setColumn(2, (matR.getColumn(0))^(matR.getColumn(1)));
                rotationEnforce(matR);

                // translation vector
                vecC = (matKinverse * h3 * l).getHomogeneous(1.0);

                // R contribution
                for(unsigned int i=0; i<3; i++)
                        for(unsigned int j=0; j<3; j++)
                                matP[i][j] = matR[i][j];

                // t contribution
                matP.setColumn(3,vecC.getUnhomogeneous());
                        
                // K contribution
                matP = matK * matP;

                // Compute C thanks to t = -RC
                vecC = (-matR.getTranspose()*vecC.getUnhomogeneous()).getHomogeneous(1.0);
        }


        /*
         * update the camera external parameters
         * M : 4x3 matrix : M=[R|-RC] (C=camera position)
         */
        void ProjectiveCamera::updateM(const kn::Matrixd &M){
                if(M.rows()!=3 && M.columns()!=4) 
                  throw VisionException("updateM : M should be a 3x4 matrix");

                matP = matK * M;
                decompose();
        }


        /*
         * update the camera projection matrix parameters
         * P : 4x3 matrix : P=K[R|-RC] (C=camera position)
         */
        void ProjectiveCamera::updateP(const kn::Matrixd &P){
                // check P size
                if(matP.rows() != 3 && matP.columns()!= 4)
                        throw VisionException("invalid size : the input P matrix is not a 3x4 matrix");

                matP = P;
                // extract K, R, C from P
                decompose();

                // inverse K 
                matKinverse = kn::inverseMatrixSVD(matK);
         }


        /*
         * \brief return M Matrix (M=[R|t]) (read-only)
         */
        kn::Matrix<double> ProjectiveCamera::M() const{

                // P = K[R|-RC] = KM
                kn::Matrix<double> M(3,4);

                // R contribution
                for(unsigned int i=0; i<3; i++)
                        for(unsigned int j=0; j<3; j++)
                                M[i][j] = matR[i][j];

                // -RC contribution
                M.setColumn(3,-matR*vecC.getUnhomogeneous());

                // finish
                return M;
        }


        /*
         * \brief compute the GL projection matrix
         * The way to use this matrix in an OpenGL program is :
         * glMatrixMode(GL_PROJECTION);
         * glLoadMatrixf(GLPMat);
         */
        void ProjectiveCamera::getGLProjectionMatrix(const size_t& camerawidth, 
                                                     const size_t& cameraheight,
                                                     const float& glnear,
                                                     const float& glfar,
                                                     float GLPMat[16]) const{
                // Reset the GL matrix
                memset(GLPMat,0,16*sizeof(float));

                // Original internal parameters matrix is 
                // fx skew Cx
                // 0  fy   Cy
                // 0  0    1
                //
                // OpenGL projection Matrix is
                // 2n/(r-l) 0            (r+l)/(r-l)   0
                // 0            2n/(t-b) (t+b)/(t-b)   0
                // 0            0                 -(f+n)/(f-n) -2fn/(f-n)
                // 0            0                 -1               0
                // Keep in mind that GL matrix are transposed and r,l,b and t are in (-1,1) range
                // Y = height-y we also have to negate y parameters (elements 5 and 9)

                GLPMat[0] = float(matK[0][0])*2.0f / float(camerawidth);    // fx -> (-1,1)
                GLPMat[4] = float(matK[0][1])*2.0f / float(camerawidth);    // Skew
                GLPMat[5] = -float(matK[1][1])*2.0f / float(cameraheight);  // fy -> (-1,1)
                GLPMat[8] = 1.0f - 2.0f*float(matK[0][2]) / float(camerawidth);
                GLPMat[9] = 2.0f*float(matK[1][2]) / float(cameraheight)-1.0f;
                GLPMat[10] = -(glnear+glfar)/(glfar-glnear);
                GLPMat[11] = -1.0; 
                GLPMat[14] = -2.0*(glnear*glfar)/(glfar-glnear);

        }


        /*
         * \brief compute the GL modelview matrix
         * The way to use this matrix in an OpenGL program is :
         * glMatrixMode(GL_MODELVIEW);
         * glLoadMatrixf(GLMMat);
         */
        void ProjectiveCamera::getGLModelviewMatrix(float GLMMat[16]) const{
                
                kn::Vector3<double> tmp = -matR*vecC.getUnhomogeneous();

                // Original internal parameters matrix is 
                // R1x R1y R1z Tx
                // R2x R2y R2z Ty
                // R3x R3y R3z Tz
                //
                // OpenGL modelview Matrix (transposed version) is
                // 
                // R1x R2x R3x 0
                // R1y R2y R3y 0
                // R1z R2z R3z 0
                // Tx  Ty  Tz  1
                //
                // Since original referential is not direct we have to invert the third column of the GL matrix
                //
                // R2x  R1x -R3x 0
                // R2y  R1y -R3y 0
                // R2z  R1z -R3z 0
                // Ty   Tx  -Tz  1
                
                GLMMat[0] = matR[0][0];
                GLMMat[1] = matR[1][0];
                GLMMat[2] = -matR[2][0];
                GLMMat[3] = 0.0;
                GLMMat[4] = matR[0][1];
                GLMMat[5] = matR[1][1];
                GLMMat[6] = -matR[2][1];
                GLMMat[7] = 0.0;  
                GLMMat[8] = matR[0][2];
                GLMMat[9] = matR[1][2];
                GLMMat[10] = -matR[2][2];
                GLMMat[11] = 0.0;
                GLMMat[12] = tmp[0];
                GLMMat[13] = tmp[1];
                GLMMat[14] = -tmp[2];
                GLMMat[15] = 1.0;
        }



      /*
        * \brief compute the principal ray (axis) of the camera corresponding to the ray passing through the camera center C with direction vector m3T (3rd row of M=KR), this principal ray is represented by the axis vector v = det(M)m3T that is directed towards the front of the camera.
        * \author Vincent
        */
        kn::Vector3d ProjectiveCamera::principalRay(){
          kn::Matrix3x3d M(matK*matR);
          return kn::Vector3d(kn::determinant3x3(M) * M.getColumn(2));
        }



       /*
        * compute a weak calibration for 2 matrices related by an fundamental matrix.
        * P1=[Id|0] and P2=[e2xF|e2]
        */
         void ProjectiveCamera::weakCalibration(kn::ProjectiveCamera &cam1, 
                                                kn::ProjectiveCamera &cam2, 
                                                const kn::Matrix3x3d &F){
                // P1 = [Id|0]
                cam1.matP.setIdentity();
                cam1.decompose();
                cam1.matKinverse = kn::inverseMatrixSVD(cam1.matK);

                // e2 and e2x
                kn::Vector3d e1;
                kn::Vector3d e2;
                computeEpipoles(F,e1,e2);
                kn::Matrix3x3d e2cross;
                e2cross.cross3x3(e2);

                // P2 = [e2xF|e2]
                kn::Matrix3x3d M = e2cross * F;
                cam2.matP.setSubMatrix(0,0,M);
                cam2.matP.setColumn(3,e2);
                cam1.decompose();
                cam1.matKinverse = kn::inverseMatrixSVD(cam1.matK);
        }


        /*
        * End of Namespace
        */
}
