/*
 * AndroidLogAppender.cpp
 *
 */

#include "PortabilityImpl.hh"
#ifdef LOG4CPP_HAVE_UNISTD_H
#    include <unistd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <android/log.h>

#include <log4cpp/AndroidLogAppender.hh>



namespace log4cpp {

    AndroidLogAppender::AndroidLogAppender(const std::string& name) : 
        LayoutAppender(name) {
    }
    
    AndroidLogAppender::~AndroidLogAppender() {
    }

    void AndroidLogAppender::close() {
        // empty
    }

    void AndroidLogAppender::_append(const LoggingEvent& event) {
	std::string message(_getLayout().format(event));
	switch(event.priority)
	{
//TODO do the rest
		case log4cpp::Priority::EMERG:
		case log4cpp::Priority::ALERT:
		case log4cpp::Priority::CRIT:
		case log4cpp::Priority::ERROR:{
		__android_log_print(ANDROID_LOG_ERROR,"ubitrack",message.c_str());
		break;
		}
		case log4cpp::Priority::WARN:{
		__android_log_print(ANDROID_LOG_WARN,"ubitrack",message.c_str());
		break;
		}
		case log4cpp::Priority::NOTICE:{
		__android_log_print(ANDROID_LOG_INFO,"ubitrack",message.c_str());
		break;
		}
		case log4cpp::Priority::INFO:{
		__android_log_print(ANDROID_LOG_INFO,"ubitrack",message.c_str());
		break;
		}
		case log4cpp::Priority::DEBUG:{
		__android_log_print(ANDROID_LOG_DEBUG,"ubitrack",message.c_str());
		break;
		}
		default:{
		__android_log_print(ANDROID_LOG_VERBOSE,"ubitrack",message.c_str());
		break;
		}
	}

        
    }

    bool AndroidLogAppender::reopen() {
        return true;
    }      
}
