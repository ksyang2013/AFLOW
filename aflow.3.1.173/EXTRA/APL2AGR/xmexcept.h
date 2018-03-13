#ifndef XMLIB_EXCEPTIONS_H
#define XMLIB_EXCEPTIONS_H

#include <string>
#include <stdexcept>

namespace xmlib
{

class XMLIBRuntimeError : public std::runtime_error
{
    public:
      XMLIBRuntimeError(const std::string& s) : std::runtime_error(s) { }
};

class XMLIBLogicError : public std::logic_error
{
    public:
      XMLIBLogicError(const std::string& s) : std::logic_error(s) { }
};

}

#endif
