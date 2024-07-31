/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * SVD.cpp created in 10 2008.
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

#include "SVD.hpp"
#include "MathTools.hpp"

namespace kn{


  // decompositionSVD
  void decompositionSVD(Matrix<double>& a,
                        Vector<double>& w,
                        Matrix<double>& v){

    if(w.size()!=a.columns()){
      std::ostringstream out;
      out << a.columns();
      throw MathException("size of the second parameter is not correct (should be " + out.str() + ")" ,"decompositionSVD");
    }

    if(v.columns()!=a.columns() || v.rows()!=a.columns()){
      std::ostringstream out;
      out << a.columns();
      throw MathException("size of the third parameter is not correct (should be " + out.str() + "x" + out.str() + ")" ,"decompositionSVD");
    }

    bool flag;
    int i,its,j,jj,k,l=0,nm=0;
    double anorm, c,f=0.0,g,h,s,scale,x,y,z=0.0;
    int m = a.rows();
    int n = a.columns();
    Vector<double> rv1(n);

    g = scale = anorm = 0.0;

    for(i=1; i<=n; ++i)
      {
        l = i+1;
        rv1[i-1] = scale*g;
        g = s = scale = 0.0;
        if(i <= m)
          {
            for(k=i; k<=m; ++k)
              scale += std::fabs(a[k-1][i-1]);

            if(scale)
              {
                for(k=i; k<=m; ++k)
                  {
                    a[k-1][i-1] /= scale;
                    s += a[k-1][i-1]*a[k-1][i-1];
                  }

                f = a[i-1][i-1];
                g = -setSign(sqrt(s),f);
                h = f*g - s;
                a[i-1][i-1] = f - g;

                for(j=l; j<=n; ++j)
                  {
                    for(s=0.0,k=i; k<=m; ++k)
                      s += a[k-1][i-1]*a[k-1][j-1];

                    f = s/h;

                    for(k=i; k<=m; ++k)
                      a[k-1][j-1] += f*a[k-1][i-1];
                  }

                for(k=i; k<=m; ++k)
                  a[k-1][i-1] *= scale;
              }
          }
        w[i-1] = scale * g;
        g = s = scale = 0.0;

        if(i <= m && i != n)
          {
            for(k=l; k<=n; ++k)
              scale += std::fabs(a[i-1][k-1]);

            if(scale)
              {
                for(k=l; k<=n; ++k)
                  {
                    a[i-1][k-1] /= scale;
                    s += a[i-1][k-1]*a[i-1][k-1];
                  }

                f = a[i-1][l-1];
                g = -setSign(sqrt(s),f);
                h = f*g - s;
                a[i-1][l-1] = f-g;

                for(k=l; k<=n; ++k)
                  rv1[k-1] = a[i-1][k-1] / h;

                for(j=l; j<=m; ++j)
                  {
                    for(s=0.0,k=l; k<=n; ++k)
                      s += a[j-1][k-1] * a[i-1][k-1];

                    for(k=l; k<=n; ++k)
                      a[j-1][k-1] += s * rv1[k-1];
                  }

                for(k=l; k<=n; ++k)
                  a[i-1][k-1] *= scale;
              }
          }

        anorm = std::max(anorm,(std::fabs(w[i-1])+ std::fabs(rv1[i-1])));
      }

    for(i = n; i>=1; i--)
      {
        if(i < n)
          {
            if(g)
              {
                for(j=l; j<=n; ++j)
                  v[j-1][i-1] = (a[i-1][j-1]/a[i-1][l-1])/g;

                for(j=l; j<=n; ++j)
                  {
                    for(s=0.0,k=l; k<=n; ++k)
                      s += a[i-1][k-1] * v[k-1][j-1];

                    for(k=l; k<=n; ++k)
                      v[k-1][j-1] += s * v[k-1][i-1];
                  }
              }
            for(j=l; j<=n; ++j)
              v[i-1][j-1] = v[j-1][i-1] = 0.0;
          }

        v[i-1][i-1] = 1.0;
        g = rv1[i-1];
        l = i;
      }

    for(i=(int)std::min(m,n); i>=1; i--)
      {
        l = i+1;
        g = w[i-1];

        for(j=l; j<=n; ++j)
          a[i-1][j-1] = 0.0;

        if(g)
          {
            g = 1.0/g;

            for(j=l; j<=n; ++j)
              {
                for(s=0.0,k=l; k<=m; ++k)
                  s += a[k-1][i-1] * a[k-1][j-1];

                f = (s/a[i-1][i-1])*g;

                for(k=i;k<=m;++k)
                  a[k-1][j-1] += f*a[k-1][i-1];
              }

            for(j=i;j<=m;++j)
              a[j-1][i-1] *= g;
          }
        else for(j=i; j<=m; ++j)
          a[j-1][i-1] = 0.0;

        ++a[i-1][i-1];
      }

    for(k=n; k>=1; k--)
      {
        for(its=1; its<=30; ++its)
          {
            flag = true;

            for(l=k; l>=1; l--)
              {
                nm = l-1;
                if((double)(std::abs(rv1[l-1])+anorm) == anorm)
                  {
                    flag = false;
                    break;
                  }

                if((double)(std::abs(w[nm-1])+anorm) == anorm)
                  break;
              }
            if(flag)
              {
                c = 0.0;
                s = 1.0;

                for(i=l;i<=k;++i)
                  {
                    f=s*rv1[i-1];
                    rv1[i-1] = c*rv1[i-1];
                    if((double)(std::abs(f)+anorm) == anorm)
                      break;
                    g=w[i-1];
                    h = pythag(f,g);
                    w[i-1]=h;
                    h=1.0/h;
                    c=g*h;
                    s=-f*h;

                    for(j=1;j<=m;++j)
                      {
                        y=a[j-1][nm-1];
                        z=a[j-1][i-1];
                        a[j-1][nm-1]=y*c+z*s;
                        a[j-1][i-1]=z*c-y*s;
                      }
                  }
              }

            z=w[k-1];
            if(l==k)
              {
                if(z < 0.0)
                  {
                    w[k-1] = -z;
                    for(j=1;j<=n;++j)
                      v[j-1][k-1] = -v[j-1][k-1];
                  }
                break;
              }
            if(its == 30)  return;
            //nrerror("No convergence in 30 SVDCMP iterations");
            x=w[l-1];
            nm=k-1;
            y=w[nm-1];
            g=rv1[nm-1];
            h=rv1[k-1];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g= pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f + setSign(g,f)))-h))/x;
            c=s=1.0;
            for (j=l;j<=nm;j++)
              {
                i=j+1;
                g=rv1[i-1];
                y=w[i-1];
                h=s*g;
                g=c*g;
                z = pythag(f,h);
                rv1[j-1]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y*=c;
                for (jj=1;jj<=n;jj++)
                  {
                    x=v[jj-1][j-1];
                    z=v[jj-1][i-1];
                    v[jj-1][j-1]=x*c+z*s;
                    v[jj-1][i-1]=z*c-x*s;
                  }
                z = pythag(f,h);
                w[j-1]=z;
                if (z)
                  {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                  }
                f=(c*g)+(s*y);
                x=(c*y)-(s*g);
                for (jj=1;jj<=m;jj++)
                  {
                    y=a[jj-1][j-1];
                    z=a[jj-1][i-1];
                    a[jj-1][j-1]=y*c+z*s;
                    a[jj-1][i-1]=z*c-y*s;
                  }
              }
            rv1[l-1]=0.0;
            rv1[k-1]=f;
            w[k-1]=x;
          }
      }
  }


  // sortSingularValuesSVD
  void sortSingularValuesSVD(Matrix<double>& u,
                             Vector<double>& w,
                             Matrix<double>& v){

    if(w.size()!=u.columns()){
      std::ostringstream out;
      out << u.columns();
      throw MathException("size of the second parameter is not correct (should be " + out.str() + ")" ,"sortSingularValuesSVD");
    }

    if(v.columns()!=u.columns() || v.rows()!=u.columns()){
      std::ostringstream out;
      out << u.columns();
      throw MathException("size of the third parameter is not correct (should be " + out.str() + "x" + out.str() + ")" ,"sortSingularValuesSVD");
    }

    for(unsigned int i=0; i<w.size(); ++i){
      for(unsigned int j=0; j<w.size()-i-1; ++j){
        if(w[j] < w[j+1]){
          w.swap(j,j+1);
          v.swapColumns(j,j+1);
          u.swapColumns(j,j+1);
        }
      }
    }
  }


  // solveSVD
  void solveSVD(const Matrix<double>& u,
                const Vector<double>& w,
                const Matrix<double>& v,
                const Vector<double>& b,
                Vector<double>& x){

    if(w.size()!=u.columns()){
      std::ostringstream out;
      out << u.columns();
      throw MathException("size of the second parameter is not correct (should be " + out.str() + ")" ,"solveSVD");
    }

    if(v.columns()!=u.columns() || v.rows()!=u.columns()){
      std::ostringstream out;
      out << u.columns();
      throw MathException("size of the third parameter is not correct (should be " + out.str() + "x" + out.str() + ")" ,"solveSVD");
    }

    if(x.size()!=u.columns()){
      std::ostringstream out;
      out << u.columns();
      throw MathException("size of the fourth parameter is not correct (should be " + out.str() + ")" ,"solveSVD");
    }

    if(b.size()!=u.rows()){
      std::ostringstream out;
      out << u.rows();
      throw MathException("size of the fifth parameter is not correct (should be " + out.str() + ")" ,"solveSVD");
    }


    int jj, j, i;
    double s;

    int m = u.rows();
    int n = u.columns();

    Vector<double> tmp(n);

    // calculate UT.b and W-1
    for(j=0; j<n;++j){
      s = 0.0;
      if(w[j] != 0.0){
        for(i=0; i<m; ++i)
          s += u[i][j]*b[i];
        s /= w[j];
      }
      tmp[j] = s;
    }

    for(j=0; j<n;++j){
      s = 0.0;
      for(jj=0; jj<n; ++jj)
        s += v[j][jj]*tmp[jj];
      x[j] = s;
    }

  }

}