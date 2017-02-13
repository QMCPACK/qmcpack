#ifndef AFQMC_READHEADER_H
#define AFQMC_READHEADER_H

#include<cstdlib>
#include<iostream>
#include<fstream>
#include<vector> 
#include<string> 
#include<ctype.h>

#include "Utilities/SimpleParser.h"

#include "AFQMC/config.h"

namespace qmcplusplus
{

  bool readHeader( std::ifstream& in,
     int& NMAX, int& NMO, int& NETOT, int& NAEA, int& NAEB, int& NCA, int& NCB, int& MS2, bool& spinRestricted, int& ISYM, std::vector<IndexType>& occup_alpha, std::vector<IndexType>& occup_beta, std::vector<IndexType>& orbSymm, std::vector<IndexType>& occupPerSymm_alpha, std::vector<IndexType>& occupPerSymm_beta, bool& orderStates, bool factorizedHam)
  {

     factorizedHam=false;
     // Read header, but do not overwrite variables that are >= 0 (initialized in xml input). 
     std::vector<std::string> words;
     getwords(words,in);       
     do {
       if(words.size() == 0) { 
         app_error()<<"Format error in ASCII integral file. End of file in header. \n";
         return false;
       } 
       for(std::vector<std::string>::iterator it=words.begin(); it!=words.end(); it++) {
         if(*it == "&FCI") {
           // do nothing   
         } else if(*it == "NORB" || *it == "NMO") {
           if( it+1 == words.end() ) {
             app_error()<<"Format error in ASCII integral file. NORB \n";
             return false;
           }
           if(NMO < 0) NMO = atoi((++it)->c_str());   
           else it++;
         } else if(*it == "NMAX") {
           if( it+1 == words.end() ) {
             app_error()<<"Format error in ASCII integral file. NMAX \n";
             return false;
           }
           if(NMAX < 0) NMAX = atoi((++it)->c_str());   
           else it++;
         } else if(*it == "NAEA") {
           if( it+1 == words.end() )  {
             app_error()<<"Format error in ASCII integral file. NAEA \n";
             return false;
           }
           if(NAEA < 0) NAEA = atoi((++it)->c_str());   
           else it++;
         } else if(*it == "NAEB") {
           if( it+1 == words.end() )  {
             app_error()<<"Format error in ASCII integral file. NAEB \n";
             return false;
           }
           if(NAEB < 0) NAEB = atoi((++it)->c_str());   
           else it++;
         } else if(*it == "NCB") {
           if( it+1 == words.end() )  {
             app_error()<<"Format error in ASCII integral file. NAEB \n";
             return false;
           }
           if(NCB <= 0) NCB = atoi((++it)->c_str());  
           else it++;
         } else if(*it == "NCA") {
           if( it+1 == words.end() )  {
             app_error()<<"Format error in ASCII integral file. NAEB \n";
             return false;
           }
           if(NCA <= 0) NCA = atoi((++it)->c_str());  
           else it++;
         } else if(*it == "NELEC") {
           if( it+1 == words.end() )  {
             app_error()<<"Format error in ASCII integral file. NETOT \n";
             return false;
           }
           if(NETOT < 0) NETOT = atoi((++it)->c_str());   
           else it++;
         } else if(*it == "MS2") {
           if( it+1 == words.end() ) { 
             app_error()<<"Format error in ASCII integral file. MS2 \n";
             return false;
           }
           if(MS2 < -50) MS2 = atoi((++it)->c_str());   
           else it++;
         } else if(*it == "ORBSYM") {
           if( NMO < 0 ) { 
             app_error()<<"NMO (NORB) must be defined before ORBSYM in ASCII integral file.\n"; 
             return false;
           }
           orbSymm.clear();
           orbSymm.reserve(2*NMO);
           int n = NMO; //spinRestricted?NMO:2*NMO; 
           while(orbSymm.size() < n) {
             it++;
             if(it==words.end()) {
               getwords(words,in);
               if(words.size() == 0)
                 app_error()<<"Format error in ASCII integral file. End of file in header. \n";                
               it=words.begin();
             } 
             bool isNumber = true;
             for(std::string::const_iterator k = it->begin(); k != it->end(); ++k)
               isNumber = (isNumber&&isdigit(*k));
             if(isNumber) {
               orbSymm.push_back( atoi(it->c_str()) );
             } else {
               app_error()<<"  Format error in section ORBSYM" <<std::endl; 
               app_error()<<"  Expecting an integer, found: " <<*it <<std::endl;    
               app_error()<<"  Number of terms found so far: " <<orbSymm.size() <<std::endl;
               return false;
             }
           };
         } else if(*it == "ISYM") {
           if( it+1 == words.end() ) { 
             app_error()<<"Format error in ASCII integral file. ISYM \n";
             return false;
           }
           if(ISYM < 0) ISYM = atoi((++it)->c_str());  
           else it++;
         } else if(*it == "UHF" || *it == "IUHF") {
           if( it+1 == words.end() ) { 
             app_error()<<"Format error in ASCII integral file. UHF \n";
             return false;
           }
           int uhf = atoi((++it)->c_str());
           spinRestricted = (uhf==0); 
         } else if(*it == "OCCUP_ALPHA") {
           if( NAEA < 0 ) { 
             app_error()<<"NCA and NAEA  must be defined before OCCUP_ALPHA in ASCII integral file.\n"; 
             return false;
           }
           if( it+(NAEA) == words.end() ) {
             app_error()<<"Format error in ASCII integral file. OCCUP_ALPHA \n";
             return false;
           }
           occup_alpha.resize(NAEA);  
           it++;
           for(int i=0; i<NAEA; i++,it++) occup_alpha[i] = atoi(it->c_str())-1;  
           std::sort(occup_alpha.begin(),occup_alpha.end());
         } else if(*it == "OCCUP_BETA") {
           if( NAEB < 0 ) {
             app_error()<<"NCB and NAEB  must be defined before OCCUP_ALPHA in ASCII integral file.\n";
             return false;
           }
           if( it+(NAEB) == words.end() ) {
             app_error()<<"Format error in ASCII integral file. OCCUP_BETA \n";
             return false;
           }
           occup_beta.resize(NAEB);
           it++;
           for(int i=0; i<NAEB; i++,it++) occup_beta[i] = atoi(it->c_str())-1+NMO-NCA;
           std::sort(occup_beta.begin(),occup_beta.end());
         } else if(*it == "OCCUP") {
           if( NAEB < 0 || NAEA < 0 || NAEA != NAEB || NCA!=NCB ) {
             app_error()<<"OCCUP std::string in ASCII integral file requires NCA=NCB,NAEA=NAEB,NAEA>0,NAEB>0. \n" << std::endl;
             return false;
           }
           if( words.size() < NAEA+1 ) {
             app_error()<<"Format error in ASCII integral file. OCCUP \n" <<std::endl;
             return false;
           }
           occup_alpha.resize(NAEA);
           occup_beta.resize(NAEB);
           for(int i=0; i<NAEA; i++) occup_alpha[i] = atoi((++it)->c_str())-1;
           for(int i=0; i<NAEB; i++) occup_beta[i] = occup_alpha[i]+NMO-NCA;
           std::sort(occup_beta.begin(),occup_beta.end());
           std::sort(occup_alpha.begin(),occup_alpha.end());

         } else if(*it == "OCCUPSYMM_ALPHA") {
           if( NAEA < 0 ) {
             app_error()<<"NCA and NAEA  must be defined before OCCUPSYMM_ALPHA in ASCII integral file.\n";
             return false;
           }
           occupPerSymm_alpha.clear();
           for(;(it+1)!=words.end();it++)
             occupPerSymm_alpha.push_back(atoi((it+1)->c_str()));
           int cnt=0;
           for(int i=0; i<occupPerSymm_alpha.size(); i++)
             cnt+=occupPerSymm_alpha[i]; 
           if(cnt != NCA+NAEA) {
             app_error()<<" Problems with OCCUPSYMM_ALPHA. Number of orbitals does not add to NCA+NAEA. \n"; 
             return false;
           }
           if(!orderStates) orderStates=true;
         } else if(*it == "OCCUPSYMM_BETA") {
           if( NAEB < 0 ) {
             app_error()<<"NCB and NAEB  must be defined before OCCUPSYMM_BETA in ASCII integral file.\n";   
             return false;
           }
           occupPerSymm_beta.clear();
           for(;(it+1)!=words.end();it++)
             occupPerSymm_beta.push_back(atoi((it+1)->c_str()));
           int cnt=0;
           for(int i=0; i<occupPerSymm_beta.size(); i++)
             cnt+=occupPerSymm_beta[i];
           if(cnt != NCB+NAEB) {
             app_error()<<" Problems with OCCUPSYMM_BETA. Number of orbitals does not add to NCB+NAEB. \n";
             return false;
           }
           if(!orderStates) orderStates=true;
         } else if(*it == "OCCUPSYMM") {

           if( NAEB < 0 || NAEA < 0 || NAEA != NAEB || NCA!=NCB ) {
             app_error()<<"OCCUPSYMM std::string in ASCII integral file requires NCA=NCB,NAEA=NAEB,NAEA>0,NAEB>0. \n" << std::endl;
             return false;
           }
           occupPerSymm_alpha.clear();
           occupPerSymm_beta.clear();
           for(;(it+1)!=words.end();it++) {
             occupPerSymm_alpha.push_back(atoi((it+1)->c_str()));
             occupPerSymm_beta.push_back(atoi((it+1)->c_str()));
           }
           int cnt=0;
           for(int i=0; i<occupPerSymm_alpha.size(); i++)
             cnt+=occupPerSymm_alpha[i];
           if(cnt != NCA+NAEA) {
             app_error()<<" Problems with OCCUPSYMM_ALPHA. Number of orbitals does not add to NCA+NAEA. \n";
             return false;
           }
           cnt=0;
           for(int i=0; i<occupPerSymm_beta.size(); i++)
             cnt+=occupPerSymm_beta[i];
           if(cnt != NCB+NAEB) {
             app_error()<<" Problems with OCCUPSYMM_BETA. Number of orbitals does not add to NCB+NAEB. \n";	  
             return false;
           }
           if(!orderStates) orderStates=true;
         } else if(*it == "NPROP" || *it == "PROPBITLEN") {
           break; // ignore the rest of the line 
         } else {
           app_log()<<"Ignoring unknown tag in ASCII integral file: " <<*it <<std::endl;
         }
       }
       getwords(words,in);       
       if(in.eof() && words.size() == 0)
         app_error()<<"Format error in ASCII integral file. End of file in header. \n";
       while(!in.eof() && words.size() == 0) {
         if(in.eof())
           app_error()<<"Format error in ASCII integral file. End of file in header. \n";
         getwords(words,in);       
       }
//       if(words.size() == 0) 
//         app_error()<<"Format error in ASCII integral file. End of file in header. \n";
     } while( (words[0].find(std::string("/"))==std::string::npos && words[0].find(std::string("&END"))==std::string::npos));  

  if(NMAX < 0) NMAX = NMO;

  return true;
}
 
}

#endif
