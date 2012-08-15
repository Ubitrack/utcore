This version of log4cpp is based on log4cpp 0.3.5rc3, but some improvements 
were made for the Ubitrack library:

* Less virtual functions in Category
* Category::getChainedPriority no longer walks through hierarchy, priority 
  changes instead are propagated when setPriority is called.
* Support for file name and line number: use %f (without path), %F (with path) and %l in PatternLayout
* Added LOG4CPP_DEBUG, ... macros to Category.hh for performance improvement
  and file/line numbers
* New log level: TRACE (like log4j)