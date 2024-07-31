/***************************************************************************\
 * Copyright (C) by University Paris-Est - MISS team   
 * VisionException.hpp created in 12 2008.
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
\***************************************************************************/

/*
 * Anti-doublon
 */
#ifndef __OPENKN_VISION__VISIONEXCEPTION_HPP__
#define __OPENKN_VISION__VISIONEXCEPTION_HPP__

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

                class UserApp_API VisionException : std::exception{

            /*
             * Constructor & destructors
             */      
        public:
            VisionException(const std::string& err="", const std::string& funcname="");
            VisionException(const VisionException & e);
                ~VisionException() throw() {}
      
        private :
            std::string str;
        public :
            inline std::string errorString() const {return str;}
                virtual const char* what() const throw() {return str.c_str();}
        };

        inline std::ostream& operator <<(std::ostream& stream, const VisionException & err){
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
