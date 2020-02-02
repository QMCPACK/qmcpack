

#include<cmath>
#include<cassert>
#include<iostream>
#include<cstdlib>
#include<string>
#include <bitset>
#include <sys/time.h>
#include <ctime>
#if defined(USE_MPI)
  #include<mpi.h>
#endif
#include"Utilities/SQCUtils.h"
#include "sprng.h"

void SQCAbort(std::string text, int id) {
  std::cerr<<text <<std::endl;
  std::cerr.flush();
  std::cout.flush();
  if(id==0) id=1;
#if defined(USE_MPI)
  MPI_Abort(MPI_COMM_WORLD,id);
#else
  exit(id);
#endif
}

unsigned long int toUL(unsigned char *arr, int offset)
{
  unsigned long int res;
  memcpy(&res,&(arr[offset]),sizeof(unsigned long int));
  return res;
}

/*
#if defined(_LINUX_)
  #include "sys/sysinfo.h"
  inline size_t freemem()
  {
    struct sysinfo si;
    sysinfo(&si);
    si.freeram+=si.bufferram;
    return si.freeram>>20;
  }
#else
  inline size_t freemem()
  {
     return 0;
  }
#endif
*/

myRNG::myRNG(int np, int nt, int rank, int &seed):nth(nt),npr(np),rk(rank) {
   strm = new int*[nth];
   if ( seed == 0 )
     seed = make_sprng_seed();
   for(int i=0; i<nth; i++)
     strm[i] = init_sprng(nth*rk+i,npr*nth,seed,SPRNG_DEFAULT);
}

/*
int myRNG::Irand(int n) {
  return isprng(strm[n]);
}

double myRNG::Drand(int n) {
  return sprng(strm[n]);
}
*/

double myTimer::getTime() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
}

myTimer::myTimer()
{
   for(int i=0; i<100; i++) {
     timeInt[i]=0.0;
     counter[i]=0;
     timeI[i]=0;
   }
   tag=0;
}

int myTimer::getTag() {
  int n=tag;
  tag++;
  return n;
}
/*
void myTimer::start(int n)  {
// assert(n<100)
  timeI[n] = getTime();
}

void myTimer::stop(int n)  {
// assert(n<100)
  double tm=getTime();
  timeInt[n]+=(tm-timeI[n]);
  counter[n]++;
}

double myTimer::average(int n)  {
  return timeInt[n]/double(counter[n]);
}
*/

void myTimer::reset(int n)  {
  timeInt[n]=timeI[n]=0.0;
  counter[n]=0;
}

int parseLine(std::ifstream &in, vector<std::string> &text2)
{
   int len, num, pos1, pos2, pos3;
   std::string text;

   getline(in,text);

   len = text.length();

   pos1 = text.find_first_not_of (" \t\n=,");
   if(pos1 ==  std::string::npos || text[pos1] == '#')
      return 10;

   text2.clear();
   pos2 = 0;
   while(pos1 != std::string::npos)
   {
      pos2 = text.find_first_of(" \t\n=,",pos1);
      text2.push_back(text.substr(pos1, pos2-pos1));
      pos1 = text.find_first_not_of(" \t\n=,",pos2);
   }
   return 0;
}

std::string itostr(int value, int base)
{
        enum { kMaxDigits = 35 };

        std::string buf;

        buf.reserve( kMaxDigits );

        if (base < 2 || base > 16) return buf;

        int quotient = value;

        do {

                buf += "0123456789abcdef"[ std::abs( quotient % base ) ];

                quotient /= base;

        } while ( quotient );

        if ( value < 0 && base == 10) buf += '-';

        reverse( buf.begin(), buf.end() );

        return buf;
}

std::string mystring_lower(std::string &str) {
 std::string res=str;
 std::transform(str.begin(),str.end(),res.begin(),::tolower);
 return res;
}


