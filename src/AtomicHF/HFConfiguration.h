#ifndef OHMMS_MAIN_HFCONFIGURATION_H
#define OHMMS_MAIN_HFCONFIGURATION_H
#include <iostream>
// #if !defined(LOGMSG)
// #define LOGMSG(msg) {std::cout << "HF: " << msg << std::endl;}
// #endif

// #if !defined(DEBUGMSG)
// #define DEBUGMSG(msg) {std::cout << "DEBUG: " << msg << std::endl;}
// #endif

// #if !defined(ERRORMSG)
// #define ERRORMSG(msg) {std::cout << "ERROR: " << msg << std::endl;}
// #endif
// int Z = 0;

#if !defined(LOGMSG)
#define LOGMSG(msg) {std::cout<< "QMC " << msg << std::endl;}
#endif

#if !defined(WARNMSG)
#define WARNMSG(msg) {std::cout<< "WARNING " << msg << std::endl;}
#endif

#if !defined(DEBUGMSG)
#if defined(PRINT_DEBUG)
#define DEBUGMSG(msg) {std::cout<< "DEBUG " << msg << std::endl;}
#else
#define DEBUGMSG(msg)
#endif
#endif

#if !defined(ERRORMSG)
#define ERRORMSG(msg) {std::cout<< "ERROR " << msg << std::endl;}
#endif

#if !defined(XMLReport)
//#if defined(PRINT_DEBUG)
#define XMLReport(msg) {std::cout<< "XML " << msg << std::endl;}
//#else
//#define XMLReport(msg)
//#endif
#endif

#endif

