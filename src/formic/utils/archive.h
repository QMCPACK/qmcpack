///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/archive.h
///
/// \brief   header file for the archiving class
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FORMIC_ARCHIVE_HEADER
#define FORMIC_ARCHIVE_HEADER

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <complex>

#include <formic/utils/exception.h>
#include <formic/utils/matrix.h>

namespace formic {

  class Archive {

    protected:

      typedef formic::Matrix<                 int  >   i_mat;
      typedef formic::Matrix<              double  >   d_mat;
      typedef formic::Matrix< std::complex<double> >   c_mat;

      typedef formic::ConstMatrix<                 int  >   i_const_mat;
      typedef formic::ConstMatrix<              double  >   d_const_mat;
      typedef formic::ConstMatrix< std::complex<double> >   c_const_mat;

      /// \brief  returns true if the internal stream is in a bad state (for example if the last io operation failed) and false otherwise
      virtual bool fail() const = 0;

    public:

      virtual ~Archive() {} //{ std::cout << "destroyed Archive" << std::endl; }

      /// \brief  returns true if the internal stream is in a bad state (for example if the last io operation failed) and false otherwise
      bool operator!() const { return this->fail(); }

      /// \brief  returns null (i.e. false) if the internal stream's state is bad and a pointer to the object itself otherwise
      operator void*() const { return ( this->fail() ? 0 : (void*)this ); }

      virtual int getline(std::string & str) = 0;
      virtual int getline(std::string & str, char delim) = 0;

      virtual Archive & operator>>(                              bool& val ) = 0;
      virtual Archive & operator>>(                             short& val ) = 0;
      virtual Archive & operator>>(                    unsigned short& val ) = 0;
      virtual Archive & operator>>(                               int& val ) = 0;
      virtual Archive & operator>>(                      unsigned int& val ) = 0;
      virtual Archive & operator>>(                              long& val ) = 0;
      virtual Archive & operator>>(                     unsigned long& val ) = 0;
      virtual Archive & operator>>(                             float& val ) = 0;
      virtual Archive & operator>>(                            double& val ) = 0;
      virtual Archive & operator>>(                       long double& val ) = 0;
      virtual Archive & operator>>(              std::complex<double>& val ) = 0;
      virtual Archive & operator>>(                             void*& val ) = 0;
      virtual Archive & operator>>(                              char& val ) = 0;
      virtual Archive & operator>>(                       signed char& val ) = 0;
      virtual Archive & operator>>(                     unsigned char& val ) = 0;

      //virtual Archive & operator>>(                    std::streambuf*  sb ) = 0;

      //virtual Archive & operator>>(                              char* str ) = 0;
      //virtual Archive & operator>>(                       signed char* str ) = 0;
      //virtual Archive & operator>>(                     unsigned char* str ) = 0;
      virtual Archive & operator>>(                      std::string & str ) = 0;

      virtual Archive & operator>>(   std::istream& ( *pf )(std::istream&) ) = 0;
      virtual Archive & operator>>(           std::ios& ( *pf )(std::ios&) ) = 0;
      virtual Archive & operator>>( std::ios_base& ( *pf )(std::ios_base&) ) = 0;

      virtual Archive & operator>>(                std::vector<short>& val ) = 0;
      virtual Archive & operator>>(       std::vector<unsigned short>& val ) = 0;
      virtual Archive & operator>>(                  std::vector<int>& val ) = 0;
      virtual Archive & operator>>(         std::vector<unsigned int>& val ) = 0;
      virtual Archive & operator>>(                 std::vector<long>& val ) = 0;
      virtual Archive & operator>>(        std::vector<unsigned long>& val ) = 0;
      virtual Archive & operator>>(                std::vector<float>& val ) = 0;
      virtual Archive & operator>>(               std::vector<double>& val ) = 0;
      virtual Archive & operator>>(          std::vector<long double>& val ) = 0;
      virtual Archive & operator>>(                 std::vector<char>& val ) = 0;
      virtual Archive & operator>>(          std::vector<signed char>& val ) = 0;
      virtual Archive & operator>>(        std::vector<unsigned char>& val ) = 0;
      virtual Archive & operator>>(          std::vector<std::string>& val ) = 0;
      virtual Archive & operator>>(std::vector<std::complex<double> >& val ) = 0;

      virtual Archive & operator>>(                   std::set<short>& val ) = 0;
      virtual Archive & operator>>(          std::set<unsigned short>& val ) = 0;
      virtual Archive & operator>>(                     std::set<int>& val ) = 0;
      virtual Archive & operator>>(            std::set<unsigned int>& val ) = 0;
      virtual Archive & operator>>(                    std::set<long>& val ) = 0;
      virtual Archive & operator>>(           std::set<unsigned long>& val ) = 0;
      virtual Archive & operator>>(                    std::set<char>& val ) = 0;
      virtual Archive & operator>>(             std::set<signed char>& val ) = 0;
      virtual Archive & operator>>(           std::set<unsigned char>& val ) = 0;
      virtual Archive & operator>>(             std::set<std::string>& val ) = 0;

      virtual Archive & operator>>(                             i_mat& val ) = 0;
      virtual Archive & operator>>(                             d_mat& val ) = 0;
      virtual Archive & operator>>(                             c_mat& val ) = 0;

      virtual Archive & operator<<(                         const bool val ) = 0;
      virtual Archive & operator<<(                        const short val ) = 0;
      virtual Archive & operator<<(               const unsigned short val ) = 0;
      virtual Archive & operator<<(                          const int val ) = 0;
      virtual Archive & operator<<(                 const unsigned int val ) = 0;
      virtual Archive & operator<<(                         const long val ) = 0;
      virtual Archive & operator<<(                const unsigned long val ) = 0;
      virtual Archive & operator<<(                        const float val ) = 0;
      virtual Archive & operator<<(                       const double val ) = 0;
      virtual Archive & operator<<(                  const long double val ) = 0;
      virtual Archive & operator<<(        const std::complex<double>& val ) = 0;
      virtual Archive & operator<<(                        const void* val ) = 0;
      virtual Archive & operator<<(                         const char val ) = 0;
      virtual Archive & operator<<(                  const signed char val ) = 0;
      virtual Archive & operator<<(                const unsigned char val ) = 0;

      //virtual Archive & operator<<(                    std::streambuf*  sb ) = 0;

      virtual Archive & operator<<(                        const char* str ) = 0;
      //virtual Archive & operator<<(                 const signed char* str ) = 0;
      //virtual Archive & operator<<(               const unsigned char* str ) = 0;
      virtual Archive & operator<<(                const std::string & str ) = 0;

      virtual Archive & operator<<(   std::ostream& ( *pf )(std::ostream&) ) = 0;
      virtual Archive & operator<<(           std::ios& ( *pf )(std::ios&) ) = 0;
      virtual Archive & operator<<( std::ios_base& ( *pf )(std::ios_base&) ) = 0;

      virtual Archive & operator<<(          const std::vector<short>& val ) = 0;
      virtual Archive & operator<<( const std::vector<unsigned short>& val ) = 0;
      virtual Archive & operator<<(            const std::vector<int>& val ) = 0;
      virtual Archive & operator<<(   const std::vector<unsigned int>& val ) = 0;
      virtual Archive & operator<<(           const std::vector<long>& val ) = 0;
      virtual Archive & operator<<(  const std::vector<unsigned long>& val ) = 0;
      virtual Archive & operator<<(          const std::vector<float>& val ) = 0;
      virtual Archive & operator<<(         const std::vector<double>& val ) = 0;
      virtual Archive & operator<<(    const std::vector<long double>& val ) = 0;
      virtual Archive & operator<<(           const std::vector<char>& val ) = 0;
      virtual Archive & operator<<(    const std::vector<signed char>& val ) = 0;
      virtual Archive & operator<<(  const std::vector<unsigned char>& val ) = 0;
      virtual Archive & operator<<(    const std::vector<std::string>& val ) = 0;
      virtual Archive & operator<<( const std::vector<std::complex<double> >& val ) = 0;

      virtual Archive & operator<<(             const std::set<short>& val ) = 0;
      virtual Archive & operator<<(    const std::set<unsigned short>& val ) = 0;
      virtual Archive & operator<<(               const std::set<int>& val ) = 0;
      virtual Archive & operator<<(      const std::set<unsigned int>& val ) = 0;
      virtual Archive & operator<<(              const std::set<long>& val ) = 0;
      virtual Archive & operator<<(     const std::set<unsigned long>& val ) = 0;
      virtual Archive & operator<<(              const std::set<char>& val ) = 0;
      virtual Archive & operator<<(       const std::set<signed char>& val ) = 0;
      virtual Archive & operator<<(     const std::set<unsigned char>& val ) = 0;
      virtual Archive & operator<<(       const std::set<std::string>& val ) = 0;

      virtual Archive & operator<<(               const   i_const_mat& val ) = 0;
      virtual Archive & operator<<(               const   d_const_mat& val ) = 0;
      virtual Archive & operator<<(               const   c_const_mat& val ) = 0;

  };

  template<class T> class TextArchiveTemplate : public Archive {

    protected:

      boost::shared_ptr<T> _s;

      template<class S> Archive & write_fmat(const formic::ConstMatrix<S> & val) {
        // we print the matrix in c-indexing, but store it it internally in fortran-indexing, so we need to print out its transpose
        const formic::ConstMatrix<S> tp = val.t();
        return *this << std::endl << val.rows() << " " << val.cols() << " " << std::vector<S>(tp.begin(), tp.end());
      }

      template<class S> Archive & read_fmat(formic::Matrix<S> & val) {
        size_t n, m;
        std::vector<S> temp;
        *this >> n >> m >> temp;
        // we print the matrix in c-indexing, but store it it internally in fortran-indexing, so we need to transpose what we read in
        val = ( formic::Matrix<S>(m,n) <<= temp ).tip();
        return *this;
      }

      template<class U> Archive & write_vector(const std::vector<U> & val) {
        *this << std::endl;
        *this << '[' << std::endl;
        for (size_t i = 0; i < val.size(); i++)
          *this << " " << val[i] << std::endl;
        *this << ']' << std::endl;
        return *this;
      }

      template<class U> Archive & read_vector(std::vector<U> & val) {

        // read in the beginning-of-vector parentheses
        {
          char c;
          *this >> c;
          if (c != '[')
            throw formic::Exception("expected the \"[\" character at the start of a vector, found \"%s\" instead") % c;
        }

        // read in the body of the vector, which is terminated by a ']' preceeded by whitespace, e.g. " ]"
        std::string vec_text;
        std::getline(*_s, vec_text, ']');
        if ( _s->eof() || _s->fail() )
          throw formic::Exception("expected a \"]\" character preceeded by whitespace at the end of the vector");

        // count the number of elements and resize the vector
        size_t n = 0;
        {
          std::stringstream vec_stream(vec_text);
          for (U x; vec_stream >> x; n++) {}
        }
        val.resize(n);

        // populate the vector
        {
          std::stringstream vec_stream(vec_text);
          for (size_t i = 0; i < n; i++)
            vec_stream >> val[i];
        }

        // return the archive
        return *this;

      }

      template<class U> Archive & write_set(const std::set<U> & val) {
        *this << std::endl;
        *this << '[' << std::endl;
        for (typename std::set<U>::const_iterator it = val.begin(); it != val.end(); it++)
          *this << " " << *it << std::endl;
        *this << ']' << std::endl;
        return *this;
      }

      template<class U> Archive & read_set(std::set<U> & val) {
        val.clear();
        std::vector<U> vec;
        *this >> vec;
        val.insert(vec.begin(), vec.end());
        return *this;
      }

      /// \brief  returns true if the internal stream is in a bad state (for example if the last io operation failed) and false otherwise
      bool fail() const {
        if ( *_s )
          return false;
        return true;
      }

    public:

      ~TextArchiveTemplate() {} //{ std::cout << "destroyed TextArchiveTemplate" << std::endl; }

      int getline(std::string & str) { std::getline(*_s, str); return 0; }
      int getline(std::string & str, char delim) { std::getline(*_s, str, delim); return 0; }

      Archive & operator>>(                              bool& val ) { *_s >> std::boolalpha >> val; return *this; }
      Archive & operator>>(                             short& val ) { *_s >> val; return *this; }
      Archive & operator>>(                    unsigned short& val ) { *_s >> val; return *this; }
      Archive & operator>>(                               int& val ) { *_s >> val; return *this; }
      Archive & operator>>(                      unsigned int& val ) { *_s >> val; return *this; }
      Archive & operator>>(                              long& val ) { *_s >> val; return *this; }
      Archive & operator>>(                     unsigned long& val ) { *_s >> val; return *this; }
      Archive & operator>>(                             float& val ) { *_s >> val; return *this; }
      Archive & operator>>(                            double& val ) { *_s >> val; return *this; }
      Archive & operator>>(                       long double& val ) { *_s >> val; return *this; }
      Archive & operator>>(             std::complex<double> & val ) { *_s >> val; return *this; }
      Archive & operator>>(                             void*& val ) { *_s >> val; return *this; }
      Archive & operator>>(                              char& val ) { *_s >> val; return *this; }
      Archive & operator>>(                       signed char& val ) { *_s >> val; return *this; }
      Archive & operator>>(                     unsigned char& val ) { *_s >> val; return *this; }

      //Archive & operator>>(                    std::streambuf*  sb ) { *_s >>  sb; return *this; }

      //Archive & operator>>(                              char* str ) { *_s >> str; return *this; }
      //Archive & operator>>(                       signed char* str ) { *_s >> str; return *this; }
      //Archive & operator>>(                     unsigned char* str ) { *_s >> str; return *this; }
      Archive & operator>>(                      std::string & str ) { *_s >> str; return *this; }

      Archive & operator>>(   std::istream& ( *pf )(std::istream&) ) { *_s >>  pf; return *this; }
      Archive & operator>>(           std::ios& ( *pf )(std::ios&) ) { *_s >>  pf; return *this; }
      Archive & operator>>( std::ios_base& ( *pf )(std::ios_base&) ) { *_s >>  pf; return *this; }

      Archive & operator>>(                std::vector<short>& val ) { return this->read_vector(val); }
      Archive & operator>>(       std::vector<unsigned short>& val ) { return this->read_vector(val); }
      Archive & operator>>(                  std::vector<int>& val ) { return this->read_vector(val); }
      Archive & operator>>(         std::vector<unsigned int>& val ) { return this->read_vector(val); }
      Archive & operator>>(                 std::vector<long>& val ) { return this->read_vector(val); }
      Archive & operator>>(        std::vector<unsigned long>& val ) { return this->read_vector(val); }
      Archive & operator>>(                std::vector<float>& val ) { return this->read_vector(val); }
      Archive & operator>>(               std::vector<double>& val ) { return this->read_vector(val); }
      Archive & operator>>(          std::vector<long double>& val ) { return this->read_vector(val); }
      Archive & operator>>(                 std::vector<char>& val ) { return this->read_vector(val); }
      Archive & operator>>(          std::vector<signed char>& val ) { return this->read_vector(val); }
      Archive & operator>>(        std::vector<unsigned char>& val ) { return this->read_vector(val); }
      Archive & operator>>(          std::vector<std::string>& val ) { return this->read_vector(val); }
      Archive & operator>>(std::vector<std::complex<double> >& val ) { return this->read_vector(val); }

      Archive & operator>>(                   std::set<short>& val ) { return this->read_set(val); }
      Archive & operator>>(          std::set<unsigned short>& val ) { return this->read_set(val); }
      Archive & operator>>(                     std::set<int>& val ) { return this->read_set(val); }
      Archive & operator>>(            std::set<unsigned int>& val ) { return this->read_set(val); }
      Archive & operator>>(                    std::set<long>& val ) { return this->read_set(val); }
      Archive & operator>>(           std::set<unsigned long>& val ) { return this->read_set(val); }
      Archive & operator>>(                    std::set<char>& val ) { return this->read_set(val); }
      Archive & operator>>(             std::set<signed char>& val ) { return this->read_set(val); }
      Archive & operator>>(           std::set<unsigned char>& val ) { return this->read_set(val); }
      Archive & operator>>(             std::set<std::string>& val ) { return this->read_set(val); }

      Archive & operator>>(                             i_mat& val ) { return this->read_fmat(val); }
      Archive & operator>>(                             d_mat& val ) { return this->read_fmat(val); }
      Archive & operator>>(                             c_mat& val ) { return this->read_fmat(val); }

      Archive & operator<<(                         const bool val ) { *_s << " " << std::boolalpha << val; return *this; }
      Archive & operator<<(                        const short val ) { *_s << " " << val; return *this; }
      Archive & operator<<(               const unsigned short val ) { *_s << " " << val; return *this; }
      Archive & operator<<(                          const int val ) { *_s << " " << val; return *this; }
      Archive & operator<<(                 const unsigned int val ) { *_s << " " << val; return *this; }
      Archive & operator<<(                         const long val ) { *_s << " " << val; return *this; }
      Archive & operator<<(                const unsigned long val ) { *_s << " " << val; return *this; }
      Archive & operator<<(                        const float val ) {
        *_s << " " << std::scientific << std::setprecision(15) << val;
        return *this;
      }
      Archive & operator<<(                       const double val ) {
        *_s << " " << std::scientific << std::setprecision(15) << val;
        return *this;
      }
      Archive & operator<<(                  const long double val ) {
        *_s << " " << std::scientific << std::setprecision(15) << val;
        return *this;
      }
      Archive & operator<<(        const std::complex<double>& val ) {
        *_s << " " << std::scientific << std::setprecision(15) << val;
        return *this;
      }
      Archive & operator<<(                        const void* val ) { *_s << " " << val; return *this; }
      Archive & operator<<(                         const char val ) { *_s << " " << val; return *this; }
      Archive & operator<<(                  const signed char val ) { *_s << " " << val; return *this; }
      Archive & operator<<(                const unsigned char val ) { *_s << " " << val; return *this; }

      //Archive & operator<<(                    std::streambuf*  sb ) { *_s << " " <<  sb; return *this; }

      Archive & operator<<(                        const char* str ) { *_s << " " << str; return *this; }
      //Archive & operator<<(                 const signed char* str ) { *_s << " " << str; return *this; }
      //Archive & operator<<(               const unsigned char* str ) { *_s << " " << str; return *this; }
      Archive & operator<<(                const std::string & str ) { *_s << " " << str; return *this; }

      Archive & operator<<(   std::ostream& ( *pf )(std::ostream&) ) { *_s << pf; return *this; }
      Archive & operator<<(           std::ios& ( *pf )(std::ios&) ) { *_s << pf; return *this; }
      Archive & operator<<( std::ios_base& ( *pf )(std::ios_base&) ) { *_s << pf; return *this; }

      Archive & operator<<(          const std::vector<short>& val ) { return this->write_vector(val); }
      Archive & operator<<( const std::vector<unsigned short>& val ) { return this->write_vector(val); }
      Archive & operator<<(            const std::vector<int>& val ) { return this->write_vector(val); }
      Archive & operator<<(   const std::vector<unsigned int>& val ) { return this->write_vector(val); }
      Archive & operator<<(           const std::vector<long>& val ) { return this->write_vector(val); }
      Archive & operator<<(  const std::vector<unsigned long>& val ) { return this->write_vector(val); }
      Archive & operator<<(          const std::vector<float>& val ) { return this->write_vector(val); }
      Archive & operator<<(         const std::vector<double>& val ) { return this->write_vector(val); }
      Archive & operator<<(    const std::vector<long double>& val ) { return this->write_vector(val); }
      Archive & operator<<(           const std::vector<char>& val ) { return this->write_vector(val); }
      Archive & operator<<(    const std::vector<signed char>& val ) { return this->write_vector(val); }
      Archive & operator<<(  const std::vector<unsigned char>& val ) { return this->write_vector(val); }
      Archive & operator<<(    const std::vector<std::string>& val ) { return this->write_vector(val); }
      Archive & operator<<( const std::vector<std::complex<double> >& val ) { return this->write_vector(val); }

      Archive & operator<<(             const std::set<short>& val ) { return this->write_set(val); }
      Archive & operator<<(    const std::set<unsigned short>& val ) { return this->write_set(val); }
      Archive & operator<<(               const std::set<int>& val ) { return this->write_set(val); }
      Archive & operator<<(      const std::set<unsigned int>& val ) { return this->write_set(val); }
      Archive & operator<<(              const std::set<long>& val ) { return this->write_set(val); }
      Archive & operator<<(     const std::set<unsigned long>& val ) { return this->write_set(val); }
      Archive & operator<<(              const std::set<char>& val ) { return this->write_set(val); }
      Archive & operator<<(       const std::set<signed char>& val ) { return this->write_set(val); }
      Archive & operator<<(     const std::set<unsigned char>& val ) { return this->write_set(val); }
      Archive & operator<<(       const std::set<std::string>& val ) { return this->write_set(val); }

      Archive & operator<<(               const   i_const_mat& val ) { return this->write_fmat(val); }
      Archive & operator<<(               const   d_const_mat& val ) { return this->write_fmat(val); }
      Archive & operator<<(               const   c_const_mat& val ) { return this->write_fmat(val); }

  };

  class TextArchive : public TextArchiveTemplate<std::stringstream> {

    private:

      /// \brief  the copy constructor is disabled
      TextArchive(const TextArchive &);

      /// \brief  the assignment operator is disabled
      TextArchive & operator=(const TextArchive &);

    public:

      TextArchive() {
        _s = boost::shared_ptr<std::stringstream>( new std::stringstream() );
      }

      ~TextArchive() {} //{ std::cout << "destroyed TextArchive" << std::endl; }

      std::string str() { return _s->str(); }

  };

  class TextFileArchive : public TextArchiveTemplate<std::fstream> {

    private:

      /// \brief  the copy constructor is disabled
      TextFileArchive(const TextFileArchive &);

      /// \brief  the assignment operator is disabled
      TextFileArchive & operator=(const TextFileArchive &);

    public:

      TextFileArchive(const std::string & filename, const std::string & mode = "rwt") {

        // get the mode in which to open the file
        std::ios_base::openmode which;
        if      (mode == "r")
          which = std::ios_base::in;
        else if (mode == "w")
          which = std::ios_base::out;
        else if (mode == "wt")
          which = std::ios_base::out | std::ios_base::trunc;
        else if (mode == "rw")
          which = std::ios_base::in | std::ios_base::out;
        else if (mode == "rwt")
          which = std::ios_base::in | std::ios_base::out | std::ios_base::trunc;
        else
          throw formic::Exception("mode for TextFileArchive must be \"r\", \"w\", \"wt\", \"rw\", or \"rwt\"");

        // open the file
        _s = boost::shared_ptr<std::fstream>( new std::fstream(filename.c_str(), which) );

        // verify file is open
        if (!_s->is_open())
          throw formic::Exception("Failed to open file \"%s\" in TextFileArchive constructor.") % filename;

      }

      ~TextFileArchive() {
        _s->close();
        //std::cout << "destroyed TextFileArchive" << std::endl;
      }

  };

  template<class T> class BinaryArchiveTemplate : public Archive {

    protected:

      boost::shared_ptr<T> _s;

      template<class U> void write(const U * ptr, const size_t n) {
        _s->write((const char *)ptr, n*sizeof(U));
      }

      template<class U> void read(U * ptr, const size_t n) {
        _s->read((char *)ptr, n*sizeof(U));
      }

      template<class U> Archive & write(const U & val) {
        this->write(&val, 1);
        return *this;
      }

      template<class U> Archive & read(U & val) {
        this->read(&val, 1);
        return *this;
      }

      template<class S> Archive & write(const formic::ConstMatrix<S> & val) {
        const size_t n = val.rows();
        const size_t m = val.cols();
        this->write(n);
        this->write(m);
        if ( val.size() > 0 )
          this->write(val.begin(), val.size());
        return *this;
      }

      template<class S> Archive & read(formic::Matrix<S> & val) {
        size_t n, m;
        this->read(n);
        this->read(m);
        if ( *this ) {
          val.reset(n,m);
          if (val.size() > 0)
            this->read(val.begin(), val.size());
        } else {
          throw formic::Exception("failure in stream object when attempting to read formic::Matrix<S> data in formic::BinaryArchiveTemplate::read");
        }
        return *this;
      }

      template<class U> Archive & write(const std::vector<U> & val) {
        const size_t n = val.size();
        this->write(n);
        if (n > 0)
          this->write(&val[0], n);
        return *this;
      }

      template<class U> Archive & read(std::vector<U> & val) {
        size_t n;
        if ( this->read(n) ) {
          val.resize(n);
          if (n > 0)
            this->read(&val[0], n);
        }
        return *this;
      }

      template<class U> Archive & write(const std::set<U> & val) {
        const size_t n = val.size();
        this->write(n);
        for (typename std::set<U>::const_iterator it = val.begin(); it != val.end(); it++)
          this->write(U(*it));
        return *this;
      }

      template<class U> Archive & read(std::set<U> & val) {
        size_t n;
        if ( !this->read(n) )
          return *this;
        std::set<U> temp;
        for (size_t i = 0; i < n; i++) {
          U u;
          if ( !this->read(u) )
            return *this;
          temp.insert(u);
        }
        val.clear();
        val.insert(temp.begin(), temp.end());
        return *this;
      }

      /// \brief  returns true if the internal stream is in a bad state (for example if the last io operation failed) and false otherwise
      bool fail() const {
        if ( *_s )
          return false;
        return true;
      }

    public:

      ~BinaryArchiveTemplate() {} //{ std::cout << "destroyed BinaryArchiveTemplate" << std::endl; }

      int getline(std::string & str) { std::getline(*_s, str); return 0; }
      int getline(std::string & str, char delim) { std::getline(*_s, str, delim); return 0; }

      // define input operators
      Archive & operator>>(                              bool& val ) { return this->read(val); }
      Archive & operator>>(                             short& val ) { return this->read(val); }
      Archive & operator>>(                    unsigned short& val ) { return this->read(val); }
      Archive & operator>>(                               int& val ) { return this->read(val); }
      Archive & operator>>(                      unsigned int& val ) { return this->read(val); }
      Archive & operator>>(                              long& val ) { return this->read(val); }
      Archive & operator>>(                     unsigned long& val ) { return this->read(val); }
      Archive & operator>>(                             float& val ) { return this->read(val); }
      Archive & operator>>(                            double& val ) { return this->read(val); }
      Archive & operator>>(                       long double& val ) { return this->read(val); }
      Archive & operator>>(             std::complex<double> & val ) { return this->read(val); }
      Archive & operator>>(                             void*& val ) { return this->read(val); }
      Archive & operator>>(                              char& val ) { return this->read(val); }
      Archive & operator>>(                       signed char& val ) { return this->read(val); }
      Archive & operator>>(                     unsigned char& val ) { return this->read(val); }

      //Archive & operator>>(                    std::streambuf*  sb ) { *_s >>  sb; return *this; }

      //Archive & operator>>(                              char* str ) { return this->rea }
      //Archive & operator>>(                       signed char* str ) { return this->rea }
      //Archive & operator>>(                     unsigned char* str ) { return this->rea }
      Archive & operator>>(                      std::string & val ) {
        std::vector<char> vec;
        *this >> vec;
        val.assign(vec.begin(), vec.end());
        return *this;
      }

      Archive & operator>>(   std::istream& ( *pf )(std::istream&) ) { return *this; }
      Archive & operator>>(           std::ios& ( *pf )(std::ios&) ) { return *this; }
      Archive & operator>>( std::ios_base& ( *pf )(std::ios_base&) ) { return *this; }

      Archive & operator>>(                std::vector<short>& val ) { return this->read(val); }
      Archive & operator>>(       std::vector<unsigned short>& val ) { return this->read(val); }
      Archive & operator>>(                  std::vector<int>& val ) { return this->read(val); }
      Archive & operator>>(         std::vector<unsigned int>& val ) { return this->read(val); }
      Archive & operator>>(                 std::vector<long>& val ) { return this->read(val); }
      Archive & operator>>(        std::vector<unsigned long>& val ) { return this->read(val); }
      Archive & operator>>(                std::vector<float>& val ) { return this->read(val); }
      Archive & operator>>(               std::vector<double>& val ) { return this->read(val); }
      Archive & operator>>(          std::vector<long double>& val ) { return this->read(val); }
      Archive & operator>>(                 std::vector<char>& val ) { return this->read(val); }
      Archive & operator>>(          std::vector<signed char>& val ) { return this->read(val); }
      Archive & operator>>(        std::vector<unsigned char>& val ) { return this->read(val); }
      Archive & operator>>(          std::vector<std::string>& val ) { return this->read(val); }
      Archive & operator>>(std::vector<std::complex<double> >& val ) { return this->read(val); }

      Archive & operator>>(                   std::set<short>& val ) { return this->read(val); }
      Archive & operator>>(          std::set<unsigned short>& val ) { return this->read(val); }
      Archive & operator>>(                     std::set<int>& val ) { return this->read(val); }
      Archive & operator>>(            std::set<unsigned int>& val ) { return this->read(val); }
      Archive & operator>>(                    std::set<long>& val ) { return this->read(val); }
      Archive & operator>>(           std::set<unsigned long>& val ) { return this->read(val); }
      Archive & operator>>(                    std::set<char>& val ) { return this->read(val); }
      Archive & operator>>(             std::set<signed char>& val ) { return this->read(val); }
      Archive & operator>>(           std::set<unsigned char>& val ) { return this->read(val); }
      Archive & operator>>(             std::set<std::string>& val ) { return this->read(val); }

      Archive & operator>>(                             i_mat& val ) { return this->read(val); }
      Archive & operator>>(                             d_mat& val ) { return this->read(val); }
      Archive & operator>>(                             c_mat& val ) { return this->read(val); }

      Archive & operator<<(                         const bool val ) { return this->write(val); }
      Archive & operator<<(                        const short val ) { return this->write(val); }
      Archive & operator<<(               const unsigned short val ) { return this->write(val); }
      Archive & operator<<(                          const int val ) { return this->write(val); }
      Archive & operator<<(                 const unsigned int val ) { return this->write(val); }
      Archive & operator<<(                         const long val ) { return this->write(val); }
      Archive & operator<<(                const unsigned long val ) { return this->write(val); }
      Archive & operator<<(                        const float val ) { return this->write(val); }
      Archive & operator<<(                       const double val ) { return this->write(val); }
      Archive & operator<<(                  const long double val ) { return this->write(val); }
      Archive & operator<<(        const std::complex<double>& val ) { return this->write(val); }
      Archive & operator<<(                        const void* val ) { return this->write(val); }
      Archive & operator<<(                         const char val ) { return this->write(val); }
      Archive & operator<<(                  const signed char val ) { return this->write(val); }
      Archive & operator<<(                const unsigned char val ) { return this->write(val); }

      //Archive & operator<<(                    std::streambuf*  sb ) { *_s << " " <<  sb; return *this; }

      Archive & operator<<(                        const char* str ) { return *this << std::string(str); }
      //Archive & operator<<(                 const signed char* str ) { *_s << " " << str; return *this; }
      //Archive & operator<<(               const unsigned char* str ) { *_s << " " << str; return *this; }
      Archive & operator<<(                const std::string & val ) {
        return *this << std::vector<char>(val.begin(), val.end());
      }

      Archive & operator<<(   std::ostream& ( *pf )(std::ostream&) ) { return *this; }
      Archive & operator<<(           std::ios& ( *pf )(std::ios&) ) { return *this; }
      Archive & operator<<( std::ios_base& ( *pf )(std::ios_base&) ) { return *this; }

      Archive & operator<<(          const std::vector<short>& val ) { return this->write(val); }
      Archive & operator<<( const std::vector<unsigned short>& val ) { return this->write(val); }
      Archive & operator<<(            const std::vector<int>& val ) { return this->write(val); }
      Archive & operator<<(   const std::vector<unsigned int>& val ) { return this->write(val); }
      Archive & operator<<(           const std::vector<long>& val ) { return this->write(val); }
      Archive & operator<<(  const std::vector<unsigned long>& val ) { return this->write(val); }
      Archive & operator<<(          const std::vector<float>& val ) { return this->write(val); }
      Archive & operator<<(         const std::vector<double>& val ) { return this->write(val); }
      Archive & operator<<(    const std::vector<long double>& val ) { return this->write(val); }
      Archive & operator<<(           const std::vector<char>& val ) { return this->write(val); }
      Archive & operator<<(    const std::vector<signed char>& val ) { return this->write(val); }
      Archive & operator<<(  const std::vector<unsigned char>& val ) { return this->write(val); }
      Archive & operator<<(    const std::vector<std::string>& val ) { return this->write(val); }
      Archive & operator<<( const std::vector<std::complex<double> >& val ) { return this->write(val); }

      Archive & operator<<(             const std::set<short>& val ) { return this->write(val); }
      Archive & operator<<(    const std::set<unsigned short>& val ) { return this->write(val); }
      Archive & operator<<(               const std::set<int>& val ) { return this->write(val); }
      Archive & operator<<(      const std::set<unsigned int>& val ) { return this->write(val); }
      Archive & operator<<(              const std::set<long>& val ) { return this->write(val); }
      Archive & operator<<(     const std::set<unsigned long>& val ) { return this->write(val); }
      Archive & operator<<(              const std::set<char>& val ) { return this->write(val); }
      Archive & operator<<(       const std::set<signed char>& val ) { return this->write(val); }
      Archive & operator<<(     const std::set<unsigned char>& val ) { return this->write(val); }
      Archive & operator<<(       const std::set<std::string>& val ) { return this->write(val); }

      Archive & operator<<(               const   i_const_mat& val ) { return this->write(val); }
      Archive & operator<<(               const   d_const_mat& val ) { return this->write(val); }
      Archive & operator<<(               const   c_const_mat& val ) { return this->write(val); }

  };

  class BinaryArchive : public BinaryArchiveTemplate<std::stringstream> {

    private:

      /// \brief  the copy constructor is disabled
      BinaryArchive(const BinaryArchive &);

      /// \brief  the assignment operator is disabled
      BinaryArchive & operator=(const BinaryArchive &);

    public:

      BinaryArchive() {
        _s = boost::shared_ptr<std::stringstream>( new std::stringstream() );
      }

      ~BinaryArchive() {} //{ std::cout << "destroyed BinaryArchive" << std::endl; }

  };

  class BinaryFileArchive : public BinaryArchiveTemplate<std::fstream> {

    private:

      /// \brief  the copy constructor is disabled
      BinaryFileArchive(const BinaryFileArchive &);

      /// \brief  the assignment operator is disabled
      BinaryFileArchive & operator=(const BinaryFileArchive &);

    public:

      BinaryFileArchive() {

        _s = boost::shared_ptr<std::fstream>( new std::fstream() );

      }

      BinaryFileArchive(const std::string & filename, const std::string & mode = "rwt") {

        // get the mode in which to open the file
        std::ios_base::openmode which;
        if      (mode == "r")
          which = std::ios_base::binary | std::ios_base::in;
        else if (mode == "w")
          which = std::ios_base::binary | std::ios_base::out;
        else if (mode == "wt")
          which = std::ios_base::binary | std::ios_base::out | std::ios_base::trunc;
        else if (mode == "rw")
          which = std::ios_base::binary | std::ios_base::in | std::ios_base::out;
        else if (mode == "rwt")
          which = std::ios_base::binary | std::ios_base::in | std::ios_base::out | std::ios_base::trunc;
        else
          throw formic::Exception("mode for BinaryFileArchive must be \"r\", \"w\", \"wt\", \"rw\", or \"rwt\"");

        // open the file
        _s = boost::shared_ptr<std::fstream>( new std::fstream(filename.c_str(), which) );

        // verify file is open
        if (!_s->is_open())
          throw formic::Exception("Failed to open file \"%s\" in BinaryFileArchive constructor.") % filename;

      }

      ~BinaryFileArchive() {
        if (_s->is_open())
          _s->close();
        //std::cout << "destroyed BinaryFileArchive" << std::endl;
      }

  };

}

#endif
