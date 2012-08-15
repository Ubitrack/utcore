/*
 * Priority.cpp
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#include "PortabilityImpl.hh"
#include <log4cpp/Priority.hh>
#include <cstdlib>

namespace log4cpp {

    namespace {
		
        struct { std::string name; Priority::Value level; } names[12] = {
			{ "TRACE", Priority::TRACE },
			{ "DEBUG", Priority::DEBUG },
			{ "INFO", Priority::INFO },
            { "NOTICE", Priority::NOTICE },
			{ "WARN",  Priority::WARN },
			{ "ERROR", Priority::ERROR },
			{ "CRIT",  Priority::CRIT },
			{ "ALERT", Priority::ALERT },
            { "FATAL", Priority::FATAL },
            { "EMERG", Priority::EMERG },
			{ "NOTSET", Priority::NOTSET },
			{ "UNKNOWN", -1 }
        };
    }

    const std::string& Priority::getPriorityName(int priority) throw() {
        int i = 0;
        for (; names[ i ].level >= 0 && names[ i ].level != priority; i++);
        return names[ i ].name;
    }

    Priority::Value Priority::getPriorityValue(const std::string& priorityName) 
    throw(std::invalid_argument) {
		for (unsigned int i = 0; names[i].level >= 0; i++)
		    if (priorityName == names[i].name)
	            return names[i].level;

		char* endPointer;
		Priority::Value value = std::strtoul(priorityName.c_str(), &endPointer, 10);
		if (*endPointer != 0) {
			throw std::invalid_argument(
				std::string("unknown priority name: '") + priorityName + "'"
			);
		}
		
		return value;
    }
}
