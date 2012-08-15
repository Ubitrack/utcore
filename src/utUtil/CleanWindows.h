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
 * Provides a wrapper around the standard <Windows.h> include that resolves 
 * some problems, especially with VS2003
 *
 * @author Daniel Pustka <daniel.pustka@in.tum.de>
 */ 

#ifndef __UBITRACK_UTIL_CLEANWINDOWS_H_INCLUDED__
#define __UBITRACK_UTIL_CLEANWINDOWS_H_INCLUDED__

// target XP
#ifndef WINVER
	#define WINVER 0x0501
#endif

#ifndef _WIN32_WINNT
	#define _WIN32_WINNT WINVER
#endif

#ifndef WIN32_LEAN_AND_MEAN
	#define WIN32_LEAN_AND_MEAN
	#define UBITRACK_UNDEF_WIN32_LEAN_AND_MEAN
#endif

#ifndef NOMINMAX
	#define NOMINMAX
	#define UBITRACK_UNDEF_NOMINMNAX
#endif

#include <Windows.h>

#ifdef UBITRACK_UNDEF_NOMINMNAX
	#undef NOMINMAX
#endif

#ifdef UBITRACK_UNDEF_WIN32_LEAN_AND_MEAN
	#undef WIN32_LEAN_AND_MEAN
#endif

#endif
