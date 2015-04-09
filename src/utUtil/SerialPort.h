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
 * Serial port class.
 *
 * @author Manuel Huber <huberma@in.tum.de>
 */ 


#ifndef __UBITRACK_UTIL_SERIALPORT_H_INCLUDED__
#define __UBITRACK_UTIL_SERIALPORT_H_INCLUDED__



#include <string>
#include <utCore.h>

#ifdef _WIN32

	#include <utUtil/CleanWindows.h>

#else // *nix

	#ifndef __APPLE__
		#include <termio.h>
	#else
		#include <sys/ioctl.h>
		#include <termios.h>
	#endif

#endif

// parity defines
#define N 0
#define O 1
#define E 2

namespace Ubitrack { namespace Util {

class UBITRACK_EXPORT SerialPort
{
public:
	SerialPort( std::string port, unsigned long baudRate, int bits = 8, int parity = 0, int stop = 1 )
		: m_portName( port )
		, m_baudRate( baudRate )
		, m_portOpen( false )
		, m_bits( bits )
		, m_parity( parity )
		, m_stop( stop )
	{}

	~SerialPort( )
	{
		close();
	}


	void open( int vtime = 5, int vmin = 0 );
	void close();

	unsigned long send( const unsigned char* buffer, unsigned long size );
	unsigned long read( unsigned char* buffer, unsigned long size );

	unsigned long bytesOnRead();
	void sendBreak();
	void flush();

	unsigned long getBaudRate()
	{
		return m_baudRate;
	}

protected:

	std::string m_portName;
	unsigned long m_baudRate;
	bool m_portOpen;
	int m_bits, m_parity, m_stop;

#ifdef _WIN32

	HANDLE        m_hSerialPort;

#else // *nix

	int m_fileDescriptor;
	struct termios m_termiosCurrent;
	struct termios m_termiosOriginal;

#endif

};

} } // namespace Ubitrack::Util

#endif


