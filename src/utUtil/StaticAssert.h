/*
 * Ubitrack - Library for Ubiquitous Tracking
 * Copyright 2006, Technische Universitaet Muenchen, and individual
 * contributors as indicated by the @authors tag. See the 
 * copyright.txt in the distribution for a full listing of individual
 * contributors.
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.
 */


/**
 * @file
 * a more readable static assert
 * @author Florian Echtler <echtler@in.tum.de>
 */


#ifndef STATICASSERT_H
#define STATICASSERT_H

namespace Ubitrack { namespace Util {

	/**
	 * The default implementation of CompileTimeCheck
	 * has a constructor that accepts every argument and
	 * does nothing else.
	 */
	template< bool b > class CompileTimeCheck
	{
		public:
			CompileTimeCheck(...) {}
	};

	/**
	 * The specialization for false of CompileTimeCheck
	 * can't be instanciated as it has no constructor,
	 * leading to a compile time error.
	 */
	template<> class CompileTimeCheck<false>
	{
	};

} }


/**
 * This macro tries to create an object of class CompileTimeCheck.
 * When the test evaluates to false, this fails. The error message
 * is passed as a parameter and is visible in the compiler errors.
 * The tmp variable and sizeof assignment are there to avoid warnings
 * about unused variables and will likely not generate any unneeded
 * code (at most one single variable assignment).
 */

#define UBITRACK_STATIC_ASSERT(test, errormsg)                        \
  do {                                                                \
    struct ERROR_##errormsg {};                                       \
    typedef Ubitrack::Util::CompileTimeCheck< (test) != 0 > tmplimpl; \
    tmplimpl aTemp = tmplimpl(ERROR_##errormsg());                    \
    int tmp;                                                          \
    tmp = sizeof(aTemp);                                              \
  } while (0)

#endif // STATICASSERT_H

