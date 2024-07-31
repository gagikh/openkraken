/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team
 * MathException.hpp created in 09 2008.
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
#ifndef __OPENKN_MATH__MATHEXCEPTION_HPP__
#define __OPENKN_MATH__MATHEXCEPTION_HPP__

/*
 * External Includes
 */
#include <iostream>
#include <string>
#include <exception>

#ifdef INTERFACEDLL_EXPORTS
#define UserApp_API __declspec(dllexport)
#else
#define UserApp_API __declspec(dllimport)
#endif

/*
 * Internal Includes
 */

/*
 * Namespace
 */
namespace kn{


  /*
   * Class definition
   */

  class UserApp_API MathException : public std::exception {

    /*
     * Constructor & destructors
     */
    public:
      MathException(const std::string& err="", const std::string& funcname="");
      MathException(const MathException & e);
      ~MathException() throw() {}

    private :
      std::string str;
    public :
      inline std::string errorString() const {return str;}
      virtual const char* what() const throw() {return str.c_str();}

  };

  inline std::ostream& operator <<(std::ostream& stream, const MathException & err){
    return stream << err.errorString();
  }


  /*
   * End of Namespace
   */
}

/*
 * End of Anti-doublon
 */
#endif
