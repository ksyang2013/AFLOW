#include "aflow_apl.h"

using namespace std;

namespace apl {

//////////////////////////////////////////////////////////////////////////////

void tokenize(const string& str, vector<string>& tokens, string del) {
  string delimiters = del;

  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // Find first "non-delimiter".
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // Skip delimiters. Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);

    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

//////////////////////////////////////////////////////////////////////////////

string getVASPVersionString(const string& binfile) {
  // Get the full path to this binary
  //CO - START
  //string tmpfile = "/tmp/apl." + stringify(rand()) + ".tmp";
  //aurostd::execute( string("which ") + binfile + " > " + tmpfile );
  //ifstream infile(tmpfile.c_str(),ios_base::in);
  //if( !infile.is_open() )
  //    throw apl::APLRuntimeError("apl::getVASPVersionString(); Cannot open temporary file.");
  string fullPathBinaryName = XHOST.command(binfile);
  //infile >> fullPathBinaryName;
  //infile.close();
  //infile.clear();
  if (fullPathBinaryName.empty())
    throw apl::APLLogicError("apl::getVASPVersionString(); Cannot estimate the location of the binary file.");
  //aurostd::execute(string("rm ") + tmpfile);

  // Open our binary
  //infile.open(fullPathBinaryName.c_str(),ios::in|ios::binary);
  ifstream infile(fullPathBinaryName.c_str(), ios::in | ios::binary);
  //CO - END
  if (!infile.is_open())
    throw apl::APLLogicError("apl::getVASPVersionString(); Cannot open binary file.");
  // Read bytes...
  int bufferSize = 1024;
  char buffer[bufferSize];
  string versionString;
  while (true) {
    if (!infile.read(buffer, bufferSize))
      bufferSize = infile.gcount();

    for (int i = 0; i < bufferSize; i++) {
      if ((buffer[i] == 'v') &&
          (buffer[i + 1] == 'a') &&
          (buffer[i + 2] == 's') &&
          (buffer[i + 3] == 'p') &&
          (buffer[i + 4] == '.')) {
        int j = i + 5;
        while (buffer[j] != ' ')
          versionString.push_back(buffer[j++]);
        break;
      }
    }
    if (!versionString.empty()) break;
    if (infile.eof()) break;

    // Shift get reader to avoid the case the "vasp." string is on the boundary of two buffers...
    infile.seekg(-20, ios::cur);
  }

  infile.close();
  infile.clear();

  return versionString;
}

//////////////////////////////////////////////////////////////////////////////

//http://en.wikipedia.org/wiki/Fletcher%27s_checksum
unsigned int getFletcher32(unsigned short* data, size_t len) {
  unsigned int sum1 = 0xffff, sum2 = 0xffff;

  while (len) {
    unsigned tlen = len > 360 ? 360 : len;
    len -= tlen;
    do {
      sum1 += *data++;
      sum2 += sum1;
    } while (--tlen);
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
  }

  // Second reduction step to reduce sums to 16 bits
  sum1 = (sum1 & 0xffff) + (sum1 >> 16);
  sum2 = (sum2 & 0xffff) + (sum2 >> 16);

  return sum2 << 16 | sum1;
}

//////////////////////////////////////////////////////////////////////////////

unsigned int getFileCheckSum(const string& filename) {
  ifstream infile(filename.c_str(), ios::in | ios::binary);
  if (!infile.is_open())
    throw apl::APLLogicError("apl::getFileCheckSum(); Cannot open file.");

  // Get file length
  infile.seekg(0, ios::end);
  unsigned long length = infile.tellg();
  infile.seekg(0, ios::beg);

  // Setup read buffer (for whole file)
  if (length % 2 != 0)
    length++;
  char* buffer = new char[length];
  buffer[length - 1] = 0x00;

  // Read it in!
  infile.read(buffer, length);
  infile.close();

  // Get checksum
  unsigned int checksum = getFletcher32((unsigned short*)buffer, length >> 1);
  delete[] buffer;

  // Return value
  return checksum;
}

//////////////////////////////////////////////////////////////////////////////

void printXVector(const xvector<double>& v, bool newline) {
  cout << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  cout << setprecision(4);
  for (int m = v.lrows; m <= v.rows; m++)
    cout << setw(10) << v(m) << "   ";
  if (newline)
    cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void printXVector(const xvector<xcomplex<double> >& v) {
  cout << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  cout << setprecision(4);
  for (int m = v.lrows; m <= v.rows; m++) {
    // TODO: operator() doest not work for xvector< xcomplex<double> > ? Why?
    cout << setw(10) << "(" << v[m].re << "," << v[m].im << ")   ";
  }
  cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void printXMatrix(const xmatrix<double>& M) {
  cout << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  cout << setprecision(4);
  for (int m = M.lrows; m <= M.rows; m++) {
    for (int n = M.lcols; n <= M.cols; n++)
      cout << setw(10) << M(m, n) << "   ";
    cout << std::endl;
  }
  cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void printXMatrix(const xmatrix<xcomplex<double> >& M) {
  cout << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  cout << setprecision(4);
  for (int m = M.lrows; m <= M.rows; m++) {
    for (int n = M.lcols; n <= M.cols; n++)
      cout << setw(10) << "(" << M(m, n).re << "," << M(m, n).im << ")   ";
    cout << std::endl;
  }
  cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void printXMatrix2(const xmatrix<xcomplex<double> >& M) {
  cout << setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  cout << setprecision(4);
  cout << "Real part:" << std::endl;
  for (int m = M.lrows; m <= M.rows; m++) {
    for (int n = M.lcols; n <= M.cols; n++)
      cout << setw(10) << M(m, n).re << "  ";
    cout << std::endl;
  }
  cout << "Imag part:" << std::endl;
  for (int m = M.lrows; m <= M.rows; m++) {
    for (int n = M.lcols; n <= M.cols; n++)
      cout << setw(10) << M(m, n).im << "  ";
    cout << std::endl;
  }
  cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

}  // namespace APL
