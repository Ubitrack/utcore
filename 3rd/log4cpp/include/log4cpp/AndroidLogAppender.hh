/*
 * OstreamAppender.hh
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef _LOG4CPP_ANDROIDLOGAPPENDER_HH
#define _LOG4CPP_ANDROIDLOGAPPENDER_HH

#include <log4cpp/Portability.hh>
#include <string>
#include <log4cpp/LayoutAppender.hh>


namespace log4cpp {

    /**
     * AndroidLogAppender appends LoggingEvents to standard android logging.
     **/
    class LOG4CPP_EXPORT AndroidLogAppender : public LayoutAppender {
        public:
        AndroidLogAppender(const std::string& name);
        virtual ~AndroidLogAppender();
        
        virtual bool reopen();
        virtual void close();

        protected:
        virtual void _append(const LoggingEvent& event);

    };
}

#endif
