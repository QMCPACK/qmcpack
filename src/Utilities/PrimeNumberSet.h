//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file PrimeNumberSet.h
 * @brief define a class to generate prime numbers for random number generators
 */
#ifndef QMCPLUSPLUS_PRIME_NUMBER_SET_H
#define QMCPLUSPLUS_PRIME_NUMBER_SET_H
#include <vector>
#include <limits>

///dummy declaration
template<typename UIntType>
struct PrimeConstants {};

///specialization for uint32_t
template<>
struct PrimeConstants<uint32_t>
{
  enum {max_prime=49979687, min_num_primes=4096, max_prime_offset=779156};
};

///specialization for uint64_t
template<>
struct PrimeConstants<uint64_t>
{
  enum {max_prime=3037000501U, min_num_primes=55108, max_prime_offset=146138719U};
};

/** class to generate prime numbers
 */
template<typename UIntType>
struct PrimeNumberSet: public PrimeConstants<UIntType>
{
  typedef UIntType result_type;
  std::vector<UIntType> primes;
  using PrimeConstants<UIntType>::max_prime;
  using PrimeConstants<UIntType>::min_num_primes;
  using PrimeConstants<UIntType>::max_prime_offset;

  /** default constructor
   *
   * Reserve space for primes and construct min_prime-lowest  prime numbers
   * excluding 2
   */
  inline PrimeNumberSet()
  {
    primes.reserve(2*min_num_primes);
    primes.push_back(3);
    result_type largest=3;
    int n=min_num_primes;
    while(n)
    {
      largest+=2;
      bool is_prime=true;
      for(int j=0; j<primes.size(); j++)
      {
        if(largest%primes[j]==0)
        {
          is_prime=false;
          break;
        }
        else
          if(primes[j]*primes[j]>largest)
          {
            break;
          }
      }
      if(is_prime)
      {
        primes.push_back(largest);
        n--;
      }
    }
  }

  inline result_type operator[](size_t i) const
  {
    return primes[i];
  }

  inline size_t size() const
  {
    return primes.size();
  }

  /** add n new primes starting with an offset
   * @param offset offset of a set
   * @param n number of primes requested
   * @param primes_add contains n prime numbers if successful
   * @return true, if successful
   *
   * For i=[0,n), primes_add[i]=primes[offset+i]
   */
  inline bool get(UIntType offset, int n, std::vector<UIntType>& primes_add)
  {
    offset=offset%max_prime_offset; //roll back
    //have enough prime numbers, assign them in an array
    if(n+offset+1<primes.size())
    {
      primes_add.insert(primes_add.end(), primes.begin()+offset,primes.begin()+offset+n);
      return true;
    }
    UIntType largest=primes.back();
    int n2add=offset+n-primes.size();
    while(n2add>0 && largest<max_prime)
    {
      largest+=2;
      bool is_prime=true;
      for(int j=0; j<primes.size(); j++)
      {
        if(largest%primes[j]==0)
        {
          is_prime=false;
          break;
        }
        else
          if(primes[j]*primes[j]>largest)
          {
            break;
          }
      }
      if(is_prime)
      {
        primes.push_back(largest);
        n2add--;
      }
    }
    if(n2add)
    {
      std::ostringstream o;
      o << "  PrimeNumberSet::get Failed to generate " << n2add << " prime numbers among "
        << n << " requested.";
      APP_ABORT(o.str());
      return false; //to make compiler happy
    }
    primes_add.insert(primes_add.end(),primes.begin()+offset, primes.begin()+offset+n);
    return true;
  }

};
#endif
