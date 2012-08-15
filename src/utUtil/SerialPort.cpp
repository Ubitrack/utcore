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
 * Implementation of serial port driver.
 *
 * @author Manuel Huber <huberma@in.tum.de>
 */ 

 
#include "SerialPort.h"

#ifdef _WIN32


#include <stdio.h>
#include <string.h>
#include <time.h>


namespace Ubitrack { namespace Util {

/** @fn void SerialPort::open()
	* @brief Opens the serial port and configures it
	*/
void SerialPort::open( int vtime, int vmin )
{
	DCB   dcb;

	if (m_portOpen == true) return;

	// Open the serial port

	m_hSerialPort = CreateFile( m_portName.c_str(),            // Filename.
	                            GENERIC_READ | GENERIC_WRITE,  // Desired access.
	                            0,                             // Device isn't shared.
	                            NULL,                          // Get default security.
	                            OPEN_EXISTING,                 // Specify action to take.
	                            0,                             // Flags and Atrributes.
	                            NULL);                         // Template file.

	if (m_hSerialPort == NULL) UBITRACK_THROW("Failed to open port.");

	// Setup the COMMTIMEOUTS
	{
		COMMTIMEOUTS  timeOuts;

		timeOuts.ReadIntervalTimeout          = 0xFFFFFFFF;
		timeOuts.ReadTotalTimeoutConstant     = vtime*100;
		timeOuts.ReadTotalTimeoutMultiplier   = 0;
		timeOuts.WriteTotalTimeoutConstant    = vtime*100;
		timeOuts.WriteTotalTimeoutMultiplier  = 0;

		if (SetCommTimeouts(m_hSerialPort, &timeOuts) == 0) UBITRACK_THROW("Failed to set timeouts.");
	}

	// Configure the serial port
	{
		dcb.DCBlength = sizeof(DCB);
		GetCommState( m_hSerialPort, &dcb );
		dcb.BaudRate    = m_baudRate;
		dcb.ByteSize    = m_bits;

		if (m_parity == 0) dcb.Parity   =   NOPARITY;
		if (m_parity == 1) dcb.Parity   =  ODDPARITY;
		if (m_parity == 2) dcb.Parity   = EVENPARITY;

		if (m_stop == 1)   dcb.StopBits = ONESTOPBIT;
		if (m_stop == 2)   dcb.StopBits = TWOSTOPBITS;

		dcb.fRtsControl = RTS_CONTROL_DISABLE;

		if(!SetCommState(m_hSerialPort, &dcb) ||
			 !SetupComm(m_hSerialPort, 10000, 10000))
		{
			CloseHandle(m_hSerialPort);
			UBITRACK_THROW("Failed to set comm parameters.");
		}
	}

	// The serial port was opened successfully!
	m_portOpen = true;
}


/** @fn void SerialPort::close()
	* @brief Closes the serial port.
	*/
void SerialPort::close()
{
	if (m_hSerialPort != NULL)
	{
		CloseHandle(m_hSerialPort);
	}

	// Set member values to benign state.
	m_portOpen = false;
}


/** @fn unsigned long SerialPort::read(unsigned char* buffer, unsigned long size)
	* @brief Reads a stream of data from the serial port.
	*
	* @param buffer Copies data from serial port into this buffer.
	* @param size The size of buffer in bytes.
	* @return Number of bytes read.
	*/
unsigned long SerialPort::read( unsigned char* buffer, unsigned long size )
{
	BOOL    readStatus;
	DWORD   bytesRead;
	DWORD   errorFlags;
	COMSTAT comStat;


	if((m_portOpen != true) || (m_hSerialPort == NULL))
		return -1;

	ClearCommError(m_hSerialPort, &errorFlags, &comStat);
	if(!comStat.cbInQue)
		return 0;

	if(comStat.cbInQue < size)
	{
		bytesRead = comStat.cbInQue;
	}
	else
	{
		bytesRead = size;
	}

	readStatus = ReadFile(m_hSerialPort, buffer, bytesRead, &bytesRead, NULL);

	return bytesRead;
}


/** @fn unsigned long SerialPort::send( const unsigned char* buffer, unsigned long size )
	* @brief Sends data to the serial port.
	* 
	* @param buffer Buffer to send to serial port.
	* @param size Size of buffer to send.
	* @return Number of bytes actually sent to serial port.
	*/
unsigned long SerialPort::send( const unsigned char* buffer, unsigned long size )
{
	DWORD numBytesWritten;

	// Check to see if the port is open
	if((m_portOpen != true) || (m_hSerialPort == NULL))
		return(0);

	WriteFile(m_hSerialPort, buffer, size, &numBytesWritten, NULL);
//  printf("sent: ");
//  for (int i=0; i<numBytesWritten; ++i) {
//	  printf("%02X ", sendBuf[i]);
//  }
//  printf("\n");

	return (unsigned long)numBytesWritten;
}


/** @fn unsigned long SerialPort::bytesOnRead()
	* @brief Checks the COM port read FIFO for received characters.
	*
	* @return Number of characters in the read FIFO.
	*/
unsigned long SerialPort::bytesOnRead()
{
	DWORD   errorFlags;
	COMSTAT comStat;

	// Check to see if the port is open
	if((m_portOpen != true) || (m_hSerialPort == NULL))
		return(0);

	ClearCommError(m_hSerialPort, &errorFlags, &comStat);

	return comStat.cbInQue;
}


/** @fn void SerialPort::flush()
	* @brief Clears the serial port FIFO.  
	*
	* @param void
	* @return none
	*/
void SerialPort::flush()
{
	unsigned char buf[1000];
	// Read serial port data
	while(bytesOnRead() != 0)
		read(buf, 1000);
}


void SerialPort::sendBreak()
{
	SetCommBreak(m_hSerialPort);
	Sleep(500);
	ClearCommBreak(m_hSerialPort);
}


} } // namespace Ubitrack::Util

#else // *nix


#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

#include <iostream>

#include <utUtil/Exception.h>


namespace Ubitrack { namespace Util {

void SerialPort::open( int vtime, int vmin )
{
	if ( m_portOpen )
	{
		return;
	}

	long baud;
	switch ( m_baudRate )
	{
	case 9600:
		baud = B9600;
		break;
	case 19200:
		baud = B19200;
		break;
	case 38400:
		baud = B38400;
		break;
	case 57600:
		baud = B57600;
		break;
	case 115200:
		baud = B115200;
		break;
	case 230400:
		baud = B230400;
		break;
	default:
		UBITRACK_THROW( "Unsupported baud rate" );
	}

	int size = CS8;
	switch ( m_bits ) {
		case 5: size = CS5; break;
		case 6: size = CS6; break;
		case 7: size = CS7; break;
	}

	int parity = 0;
	if (m_parity != 0) {
		parity = PARENB;
		if (m_parity == 1)
			parity |= PARODD;
	}

	int stop = (m_stop == 2 ? CSTOPB : 0);

	m_fileDescriptor = ::open( m_portName.c_str(), O_RDWR | O_NOCTTY );
	if ( m_fileDescriptor < 0 )
		UBITRACK_THROW( "Cannot open special file" );

	if ( tcgetattr( m_fileDescriptor, &m_termiosCurrent ) )
		UBITRACK_THROW( "Cannot get port parameter" );

	m_termiosOriginal = m_termiosCurrent;

	m_termiosCurrent.c_cflag = baud | size | parity | stop | CLOCAL | CREAD;
	m_termiosCurrent.c_oflag = 0;
	m_termiosCurrent.c_iflag = IGNBRK | IGNPAR;
	m_termiosCurrent.c_lflag = 0;

	cfmakeraw( &m_termiosCurrent );

	m_termiosCurrent.c_cc[VTIME] = vtime;
	m_termiosCurrent.c_cc[VMIN]  = vmin;

	cfsetospeed( &m_termiosCurrent, baud );
	cfsetispeed( &m_termiosCurrent, baud );

	if ( tcsetattr( m_fileDescriptor, TCSANOW, &m_termiosCurrent ) < 0 )
		UBITRACK_THROW( "Cannot set port parameter" );

	m_portOpen = true;
}

void SerialPort::close()
{
	if ( !m_portOpen )
	{
		return;
	}

	if ( tcsetattr( m_fileDescriptor, TCSANOW, &m_termiosOriginal ) < 0 )
	{
		// fail silently
	}

	::close( m_fileDescriptor );
}

unsigned long SerialPort::read( unsigned char* buffer, unsigned long size )
{
	long readBytes;

	if ( !m_portOpen )
		UBITRACK_THROW( "Port is not open" );

	readBytes = ::read( m_fileDescriptor, (void*) buffer, size );

	if ( readBytes < 0 )
		UBITRACK_THROW( "Error reading bytes" );

	return readBytes;
}

unsigned long SerialPort::send( const unsigned char* buffer, unsigned long size )
{
	long writtenBytes;

	if ( !m_portOpen )
		UBITRACK_THROW( "Port is not open" );

	writtenBytes = write ( m_fileDescriptor, buffer, size );

	return writtenBytes;
}

unsigned long SerialPort::bytesOnRead()
{
	if ( !m_portOpen )
		UBITRACK_THROW( "Port is not open" );

	int onread;
	ioctl( m_fileDescriptor, FIONREAD, (char*) &onread );

	return onread;
}

void SerialPort::sendBreak()
{
	tcsendbreak( m_fileDescriptor, 0 );
}

void SerialPort::flush()
{
	tcflush( m_fileDescriptor, TCIOFLUSH );
}

} } // namespace Ubitrack::Util

#endif

