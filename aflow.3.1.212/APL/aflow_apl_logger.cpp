#include "aflow_apl.h"

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

Logger::Logger(ofstream& os, const _aflags& flags) {
  _os = &os;
  _ss.clear();
  _barCode = "00000";
  _typeofMessage = "MESSAGE";
  _moduleName = "";
  _aflowFlags = flags;
  _isQuiet = false;
}

// ///////////////////////////////////////////////////////////////////////////

Logger::Logger(Logger& l) {
  *this = l;
}

// ///////////////////////////////////////////////////////////////////////////

Logger::~Logger() {
  _ss.clear();
}

// ///////////////////////////////////////////////////////////////////////////

ofstream& Logger::getOutputStream() {
  return *_os;
}

// ///////////////////////////////////////////////////////////////////////////

void Logger::setTypeOfMessage(const string& s) {
  _typeofMessage = s;
}

// ///////////////////////////////////////////////////////////////////////////

void Logger::setQuietMode(bool b) {
  _isQuiet = b;
}

// ///////////////////////////////////////////////////////////////////////////

void Logger::setBarCode(const string& s) {
  _barCode = s;
}

void Logger::setBarCode(const char& c) {
  _barCode.clear();
  for (int i = 0; i < 5; i++) {
    _barCode.push_back(c);
  }
}

// ///////////////////////////////////////////////////////////////////////////

void Logger::setModuleName(const string& s) {
  _moduleName = s;
}

// ///////////////////////////////////////////////////////////////////////////

Logger& Logger::operator=(Logger& l) {
  if (this == &l) return *this;

  _os = l._os;
  _barCode = l._barCode;
  _typeofMessage = l._typeofMessage;
  _moduleName = l._moduleName;
  _aflowFlags = l._aflowFlags;
  _isQuiet = l._isQuiet;

  // Ok, stringstream does not have any kind of copy constructur, but...
  _ss << l._ss.str();

  return *this;
}

// ///////////////////////////////////////////////////////////////////////////

Logger& Logger::operator<<(const string& s) {
  _ss << s;
  return *this;
}

Logger& Logger::operator<<(const int& n) {
  _ss << n;
  return *this;
}

Logger& Logger::operator<<(const uint& n) {
  _ss << n;
  return *this;
}

Logger& Logger::operator<<(const double& d) {
  _ss << d;
  return *this;
}

// ///////////////////////////////////////////////////////////////////////////

Logger& Logger::operator<<(MyStreamManipulator manip) {
  return manip(*this);
}

Logger& Logger::operator<<(StandardEndLine manip) {
  manip(_ss);
  return *this;
}

// ///////////////////////////////////////////////////////////////////////////

Logger& endl(Logger& l) {
  l.endl();
  return l;
}

Logger& message(Logger& l) {
  l.setBarCode('0');
  l.setTypeOfMessage("MESSAGE");
  return l;
}

Logger& notice(Logger& l) {
  cursor_fore_green();
  l.setBarCode('0');
  l.setTypeOfMessage("NOTICE ");
  return l;
}

Logger& warning(Logger& l) {
  cursor_fore_yellow();
  l.setBarCode('W');
  l.setTypeOfMessage("WARNING");
  return l;
}

Logger& error(Logger& l) {
  cursor_fore_red();
  l.setBarCode('E');
  l.setTypeOfMessage("ERROR  ");
  return l;
}

Logger& quiet(Logger& l) {
  l.setQuietMode(true);
  return l;
}

// ///////////////////////////////////////////////////////////////////////////

void Logger::endl() {
  // Get user's message
  _ss.flush();

  // Put it to screen, if no quiet....
  if (!_isQuiet) {
    cout << _barCode;
    if (!_typeofMessage.empty())
      cout << "  " << _typeofMessage;
    if (!_moduleName.empty())
      cout << " " << _moduleName;
    if (!_ss.str().empty())
      cout << " " << _ss.str();
    cout << " " << Message(_aflowFlags, "user,host,time") << std::endl;
    // Default
    cursor_attr_none();
  }
  cout.flush();

  // Put it to log file in any case...
  *_os << _barCode;
  if (!_moduleName.empty())
    *_os << "  " << _typeofMessage;
  if (!_moduleName.empty())
    *_os << " " << _moduleName;
  if (!_ss.str().empty())
    *_os << " " << _ss.str();
  *_os << " " << Message(_aflowFlags, "user,host,time") << std::endl;
  _os->flush();

  // Reuse our _ss
  // http://www.velocityreviews.com/forums/t278533-std-stringstream-reset.html
  _ss.str(std::string());  // reconstruct std::stringstream, but some architectures can be buggy and it does not work...

  // Set default
  if (_barCode[0] != '0') {
    setBarCode('0');
    setTypeOfMessage("MESSAGE");
  }
}

// ///////////////////////////////////////////////////////////////////////////

Logger& setwidth(Logger& l, int n) {
  l._ss.width(n);
  return l;
}

LMANIP<int> sw(int n) {
  return LMANIP<int>(setwidth, n);
}

// ///////////////////////////////////////////////////////////////////////////

// Based on http://www.horstmann.com/cpp/iostreams.html
Logger& setformat(Logger& l, const char* fmt) {
  int i = 0;
  while (fmt[i] != 0) {
    if (fmt[i] != '%') {
      l._ss << fmt[i];
      i++;
    } else {
      i++;
      if (fmt[i] == '%') {
        l._ss << fmt[i];
        i++;
      } else {
        bool ok = TRUE;
        int istart = i;
        bool more = TRUE;
        int width = 0;
        int precision = 6;
        l._ss.unsetf(ios::showpos);
        ios::fmtflags flags = l._ss.flags();
        char fill = ' ';
        bool alternate = FALSE;
        while (more) {
          switch (fmt[i]) {
            case '+':
              flags |= ios::showpos;
              break;
            case '-':
              flags |= ios::left;
              break;
            case '0':
              flags |= ios::internal;
              fill = '0';
              break;
            case '#':
              alternate = TRUE;
              break;
            case ' ':
              break;
            default:
              more = FALSE;
              break;
          }
          if (more) i++;
        }
        if (isdigit(fmt[i])) {
          width = atoi(fmt + i);
          do
            i++;
          while (isdigit(fmt[i]));
        }
        if (fmt[i] == '.') {
          i++;
          precision = atoi(fmt + i);
          while (isdigit(fmt[i])) i++;
        }
        switch (fmt[i]) {
          case 'd':
            flags |= ios::dec;
            break;
          case 'x':
            flags |= ios::hex;
            if (alternate) flags |= ios::showbase;
            break;
          case 'X':
            flags |= ios::hex | ios::uppercase;
            if (alternate) flags |= ios::showbase;
            break;
          case 'o':
            flags |= ios::hex;
            if (alternate) flags |= ios::showbase;
            break;
          case 'f':
            flags |= ios::fixed;
            if (alternate) flags |= ios::showpoint;
            break;
          case 'e':
            flags |= ios::scientific;
            if (alternate) flags |= ios::showpoint;
            break;
          case 'E':
            flags |= ios::scientific | ios::uppercase;
            if (alternate) flags |= ios::showpoint;
            break;
          case 'g':
            if (alternate) flags |= ios::showpoint;
            break;
          case 'G':
            flags |= ios::uppercase;
            if (alternate) flags |= ios::showpoint;
            break;
          default:
            ok = FALSE;
            break;
        }
        i++;
        if (fmt[i] != 0) ok = FALSE;
        if (ok) {
          l._ss.unsetf(ios::adjustfield | ios::basefield | ios::floatfield);
          l._ss.setf(flags);
          l._ss.width(width);
          l._ss.precision(precision);
          l._ss.fill(fill);
        } else {
          i = istart;
        }
      }
    }
  }
  return l;
}

LMANIP<const char*> sf(const char* s) {
  return LMANIP<const char*>(setformat, s);
}

// ProgressBar ///////////////////////////////////////////////////////////////

void Logger::initProgressBar(const string& s) {
  initProgressBar(s.c_str());
}

void Logger::initProgressBar(const char* s) {
  string title = string(s);

  int titlepos = 50 - ((title.size()) / 2);

  int i = 1;
  for (; i < titlepos; i++) {
    if ((i == 1) || (i % 10 == 0))
      cout << "|";
    else
      cout << ".";
  }
  cout << " " << title << " ";
  for (i += title.size() + 2; i <= 100; i++) {
    if ((i == 1) || (i % 10 == 0))
      cout << "|";
    else
      cout << ".";
  }
  cout << std::endl;
  cout.flush();
  _progressBarLastPercent = 0;
  _progressBarPercent = 0.0;  // ME 180831
}

void Logger::updateProgressBar(int i, int n) {
  double percent = 100.0 * i / (n - 1);
  if ((percent - _progressBarLastPercent) >= 1.0) {
    cout << "|";
    cout.flush();
    _progressBarLastPercent += 1;
  }
}

// ME180831 - Update the progress bar to be updated by a set increment
void Logger::updateProgressBar(double p) {
  _progressBarPercent += 100 * p;
  double snapshot = _progressBarPercent - 100 * p;
  int bars = (int) _progressBarPercent - (int) snapshot;
  for (int i = 0; i < bars; i++) {
    cout << "|";
    cout.flush();
    _progressBarLastPercent += 1;
  }
}

void Logger::finishProgressBar() {
  // ME 180831 -- Due to numerical errors, _progressBarPercent may end up
  // being less than 100%.
  if (_progressBarLastPercent == 99) {
    cout << "|";
  }
  cout << std::endl;
  cout.flush();
  _progressBarLastPercent = 0;
  _progressBarPercent = 0.0;
}

// ///////////////////////////////////////////////////////////////////////////

}  // namespace apl
