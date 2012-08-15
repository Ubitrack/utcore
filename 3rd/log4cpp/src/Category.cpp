/*
 * Category.cpp
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#include "PortabilityImpl.hh"

#ifdef LOG4CPP_HAVE_UNISTD_H
#    include <unistd.h>
#endif

#include <log4cpp/Category.hh>
#include <log4cpp/HierarchyMaintainer.hh>
#include <log4cpp/NDC.hh>
#include "StringUtil.hh"

namespace log4cpp {

    Category& Category::getRoot() {
        return getInstance("");
    }

    void Category::setRootPriority(Priority::Value priority) {
        getRoot().setPriority(priority);
    }

    Priority::Value Category::getRootPriority() throw() {
        return getRoot().getPriority();
    }

    Category& Category::getInstance(const std::string& name) {
        return HierarchyMaintainer::getDefaultMaintainer().getInstance(name);
    }

    Category* Category::exists(const std::string& name) {
        return HierarchyMaintainer::getDefaultMaintainer().getExistingInstance(name);
    }

    std::vector<Category*>* Category::getCurrentCategories() {
        return HierarchyMaintainer::getDefaultMaintainer().
            getCurrentCategories();
    }

    void Category::shutdown() {
        HierarchyMaintainer::getDefaultMaintainer().shutdown();
    }

    Category::Category(const std::string& name, Category* parent, Priority::Value priority) : 
        _name(name),
        _parent(parent),
        _priority(priority),
		_chainedPriority(priority),
        _isAdditive(true) {
		if (parent != NULL) {
			if (_chainedPriority == Priority::NOTSET)
				_chainedPriority = parent->getChainedPriority();
			parent->_addChild( this );
		}
    }

    Category::~Category() {
        removeAllAppenders();
    }

    void Category::setPriority(Priority::Value priority)
    throw(std::invalid_argument) {
        if (priority < Priority::NOTSET) 
            _priority = _chainedPriority = priority;
		else
			if (getParent() != NULL) {
				_priority = priority;
				_chainedPriority = getParent()->getChainedPriority();
	        } else {
	            /* caller tried to set NOTSET priority to root Category. 
	               Bad caller!
	            */
	            throw std::invalid_argument("cannot set priority NOTSET on Root Category");
	        }
		/* update children */
		threading::ScopedLock lock(_childrenListMutex);
		for (ChildrenList::iterator it=_children.begin(); it!=_children.end(); it++)
			(*it)->_updateChainedPriority(_chainedPriority);
    }

	void Category::_updateChainedPriority(Priority::Value priority) throw() {
		if(getPriority() == Priority::NOTSET) {
			_chainedPriority = priority;
			/* update children */
			threading::ScopedLock lock(_childrenListMutex);
			for (ChildrenList::iterator it=_children.begin(); it!=_children.end(); it++)
				(*it)->_updateChainedPriority(_chainedPriority);
		}
	}
    
    void Category::addAppender(Appender* appender) 
    throw(std::invalid_argument) {
        if (appender) {
            threading::ScopedLock lock(_appenderSetMutex);
            {
                AppenderSet::iterator i = _appender.find(appender);
                if (_appender.end() == i) {
                    // not found
                    _appender.insert(appender);
                    _ownsAppender[appender] = true;
                }
            }
        } else {
            throw std::invalid_argument("NULL appender");
        }
    }
    
    void Category::addAppender(Appender& appender) {
        threading::ScopedLock lock(_appenderSetMutex);
        {
            AppenderSet::iterator i = _appender.find(&appender);
            if (_appender.end() == i) {
                _appender.insert(&appender);
                _ownsAppender[&appender] = false;
            }
        }
    }
    
    Appender* Category::getAppender() const {
        threading::ScopedLock lock(_appenderSetMutex);
        {
            AppenderSet::const_iterator i = _appender.begin();
            return (_appender.end() == i) ? NULL : *i;
        }
    }

    Appender* Category::getAppender(const std::string& name) const {
        threading::ScopedLock lock(_appenderSetMutex);
        {
            AppenderSet::const_iterator i = _appender.begin();
            if (_appender.end() != i) {
                // found
                return((*i)->getAppender(name));
            }
            else {
                return(NULL);
            }
        }
    }

    AppenderSet Category::getAllAppenders() const {
        threading::ScopedLock lock(_appenderSetMutex);
        {
            return _appender;
        }
    }

    void Category::removeAllAppenders() {
        threading::ScopedLock lock(_appenderSetMutex);
        {
            for (AppenderSet::iterator i = _appender.begin();
                 i != _appender.end(); i++) {
                // found
                OwnsAppenderMap::iterator i2;
                if (ownsAppender(*i, i2)) {
                    delete (*i);
                }
            }

            _ownsAppender.clear();
            _appender.clear();           
        }
    }

    void Category::removeAppender(Appender* appender) {
        threading::ScopedLock lock(_appenderSetMutex);
        {
            AppenderSet::iterator i = _appender.find(appender);
            if (_appender.end() != i) {            
                OwnsAppenderMap::iterator i2;
                if (ownsAppender(*i, i2)) {
                    _ownsAppender.erase(i2);
                    delete (*i);
                }
                _appender.erase(i);
            } else {
                // appender not found 
            }
        }
    }

    bool Category::ownsAppender(Appender* appender) const throw() {
        bool owned = false;

        threading::ScopedLock lock(_appenderSetMutex);
        {
            if (NULL != appender) {            
                OwnsAppenderMap::const_iterator i =
                    _ownsAppender.find(appender);
                if (_ownsAppender.end() != i) {
                    owned = (*i).second;
                }
            }
        }

        return owned;
    }

    /* assume lock is held */
    bool Category::ownsAppender(Appender* appender, 
                                Category::OwnsAppenderMap::iterator& i2) throw() {
        bool owned = false;

        if (NULL != appender) {
            OwnsAppenderMap::iterator i = _ownsAppender.find(appender);
            if (_ownsAppender.end() != i) {
                owned = (*i).second;
                if (owned) {
                    i2 = i;
                }
            }
        }

        return owned;
    }

    void Category::callAppenders(const LoggingEvent& event) throw() {
        threading::ScopedLock lock(_appenderSetMutex);
        {
            if (!_appender.empty()) {
                for(AppenderSet::const_iterator i = _appender.begin();
                    i != _appender.end(); i++) {
                    (*i)->doAppend(event);
                }
            }
        }
        if (getAdditivity() && (getParent() != NULL)) {
            getParent()->callAppenders(event);
        }
    }

    void Category::setAdditivity(bool additivity) {
        _isAdditive = additivity;
    }

    bool Category::getAdditivity() const throw() {
        return _isAdditive; 
    }

	void Category::_addChild(Category* child) throw() {
		threading::ScopedLock lock(_childrenListMutex);
		_children.push_back(child);
	}
	
    void Category::_logUnconditionally(Priority::Value priority, 
                                       const char* format, 
                                       va_list arguments) throw() {
        _logUnconditionally2(priority, StringUtil::vform(format, arguments));
    }
    
    void Category::_logUnconditionally2(Priority::Value priority, 
                                        const std::string& message, 
										const char* file, 
										unsigned line) throw() {
        LoggingEvent event(getName(), message, NDC::get(), priority, file, line);
        callAppenders(event);
    }
    
    void Category::log(Priority::Value priority, 
                       const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(priority)) {
            va_list va;
            va_start(va, stringFormat);
            _logUnconditionally(priority, stringFormat, va);
            va_end(va);
        }
    }

    void Category::log(Priority::Value priority, 
                       const std::string& message, const char* file, unsigned line) throw() { 
        if (isPriorityEnabled(priority))
            _logUnconditionally2(priority, message, file, line);
    }
    
    void Category::logva(Priority::Value priority, 
                         const char* stringFormat,
                         va_list va) throw() { 
        if (isPriorityEnabled(priority)) {
            _logUnconditionally(priority, stringFormat, va);
        }
    }

    void Category::trace(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::TRACE)) {
            va_list va;
            va_start(va,stringFormat);
            _logUnconditionally(Priority::TRACE, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::trace(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::TRACE))
            _logUnconditionally2(Priority::TRACE, message);
    }
    
    void Category::debug(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::DEBUG)) {
            va_list va;
            va_start(va,stringFormat);
            _logUnconditionally(Priority::DEBUG, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::debug(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::DEBUG))
            _logUnconditionally2(Priority::DEBUG, message);
    }
    
    void Category::info(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::INFO)) {
            va_list va;
            va_start(va,stringFormat);
            _logUnconditionally(Priority::INFO, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::info(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::INFO))
            _logUnconditionally2(Priority::INFO, message);
    }
    
    void Category::notice(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::NOTICE)) {
            va_list va;
            va_start(va,stringFormat);
            _logUnconditionally(Priority::NOTICE, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::notice(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::NOTICE))
            _logUnconditionally2(Priority::NOTICE, message);
    }
    
    void Category::warn(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::WARN)) {
            va_list va;
            va_start(va,stringFormat);
            _logUnconditionally(Priority::WARN, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::warn(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::WARN))
            _logUnconditionally2(Priority::WARN, message);
    }
    
    void Category::error(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::ERROR)) {
            va_list va;
            va_start(va,stringFormat);
                       _logUnconditionally(Priority::ERROR, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::error(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::ERROR))
            _logUnconditionally2(Priority::ERROR, message);
    }

    void Category::crit(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::CRIT)) {
            va_list va;
            va_start(va,stringFormat);
            _logUnconditionally(Priority::CRIT, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::crit(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::CRIT))
            _logUnconditionally2(Priority::CRIT, message);
    }

    void Category::alert(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::ALERT)) {
            va_list va;
            va_start(va,stringFormat);
            _logUnconditionally(Priority::ALERT, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::alert(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::ALERT))
            _logUnconditionally2(Priority::ALERT, message);
    }

    void Category::emerg(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::EMERG)) {
            va_list va;
            va_start(va,stringFormat);
            _logUnconditionally(Priority::EMERG, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::emerg(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::EMERG))
            _logUnconditionally2(Priority::EMERG, message);
    }

    void Category::fatal(const char* stringFormat, ...) throw() { 
        if (isPriorityEnabled(Priority::FATAL)) {
            va_list va;
            va_start(va,stringFormat);
            _logUnconditionally(Priority::FATAL, stringFormat, va);
            va_end(va);
        }
    }
    
    void Category::fatal(const std::string& message) throw() { 
        if (isPriorityEnabled(Priority::FATAL))
            _logUnconditionally2(Priority::FATAL, message);
    }

    CategoryStream Category::getStream(Priority::Value priority) {
        return CategoryStream(*this, isPriorityEnabled(priority) ?
                              priority : Priority::NOTSET);
    }

    CategoryStream Category::operator<<(Priority::Value priority) {
        return getStream(priority);
    }
} 

