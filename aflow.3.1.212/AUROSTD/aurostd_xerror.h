//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                  Marco Esters - Duke University 2018                    *
// *                                                                         *
//****************************************************************************
#ifndef _AUROSTD_XERROR_H_
#define _AUROSTD_XERROR_H_

namespace aurostd {
// BEGIN CONSTANTS
// Definitions for the named error code constant (see README, section 2.1.).

#define _GENERIC_ERROR_        1
#define _ILLEGAL_CODE_         2
#define _INPUT_ERROR_         10
#define _INPUT_UNKNOWN_       11
#define _INPUT_MISSING_       12
#define _INPUT_AMBIGUOUS_     13
#define _INPUT_ILLEGAL_       14
#define _INPUT_NUMBER_        15
#define _FILE_ERROR_          20
#define _FILE_NOT_FOUND_      21
#define _FILE_WRONG_FORMAT_   22
#define _FILE_CORRUPT_        23
#define _VALUE_ERROR_         30
#define _VALUE_ILLEGAL_       31
#define _VALUE_RANGE_         32
#define _INDEX_ERROR_         40
#define _INDEX_ILLEGAL_       41
#define _INDEX_BOUNDS_        42
#define _INDEX_MISMATCH_      43
#define _RUNTIME_ERROR_       50
#define _RUNTIME_INIT_        51
#define _ALLOC_ERROR_         60
#define _ALLOC_ALLOCATE_      61
#define _ALLOC_INSUFFICIENT_  62

// END CONSTANTS

class xerror {
  public:
    xerror(const std::string&, const std::string&, int = 1);
    xerror(const std::string&, const std::stringstream&, int = 1);
    xerror(const std::stringstream&, const std::stringstream&, int = 1);
    xerror(const std::stringstream&, const std::string&, int = 1);
    ~xerror() throw() {};
    int error_code;
    std::string what();
    std::string where();
    std::string error_message;
  private:
    int error_type, error_number;
    std::string f_name, message;

    void buildException(const std::string&, const std::string&, const int&);
    bool codeValid();
    std::string buildMessageString();
    std::string error_string();
};
} // namespace aurostd

#endif
//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *                  Marco Esters - Duke University 2018                    *
// *                                                                         *
//****************************************************************************
