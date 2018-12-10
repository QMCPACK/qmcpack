
#ifndef AFQMC_MYTIMER_H 
#define AFQMC_MYTIMER_H 


#include<tuple>
#include<vector>
#include<string>
#include<map>
#include<ctime>
#include <sys/time.h>
#include<cstdlib>
#include<ctype.h>
#include<algorithm>
#include<iostream>
#include<ostream>

class myTimer {

   // TimeData<0>: number of intervals accumulated 
   // TimeData<1>: time of last call to start 
   // TimeData<2>: cumulative sum of intervals 
   typedef std::tuple<int,double,double> TimeData;
   typedef std::tuple<int,double,double>* TimeDataPtr;

   private: 
     
     std::vector<TimeData > timer; 
     std::map<std::string, int> id2pos; 

     double getTime() {
       struct timeval tv;
       gettimeofday(&tv, NULL);
       return double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
     }
 
     // You can either register a timer (with call to add)
     // or you can just call start or reset (if timer doesn't exist, it is created there).
     // Calls to any other function with non-existent timers will do nothing 
     int getPos(const std::string& str)  {
       std::map<std::string, int>::iterator it = id2pos.find(str);
       if(it!=id2pos.end()) 
         return it->second; 
       else
         return -1;
     }

   public: 

     myTimer()
     {
       timer.reserve(100);
     }

     // if a std::string is already associated with a timer, it does nothing
     void add(const std::string& str) {
       int n = getPos(str); 
       if(n < 0) {
         timer.push_back(std::make_tuple(0,0.0,0.0)); 
         id2pos[str] = timer.size()-1; 
       }
     } 

     void start(const std::string& str) {
       int n = getPos(str);          
       if(n < 0) {
         timer.push_back(std::make_tuple(0,0.0,0.0));
         id2pos[str] = timer.size()-1;
         n = timer.size()-1;
       }  
       std::get<1>(timer[n]) = getTime();
     }

     void stop(const std::string& str) {
       double tm=getTime();
       int n = getPos(str);
       if(n >= 0) {
         std::get<0>(timer[n])++;  
         std::get<2>(timer[n]) += (tm-std::get<1>(timer[n])); 
       }
     }

     // add a time interval to the total time without incrementing the counter
     void add_time(const std::string& str) {
       double tm=getTime();
       int n = getPos(str);
       if(n >= 0) {
         std::get<2>(timer[n]) += (tm-std::get<1>(timer[n]));
       }
     }

     double elapsed(const std::string& str) {
       int n = getPos(str);
       if(n >= 0) 
         return getTime()-std::get<1>(timer[n]);
       return -1;
     }

     double average(const std::string& str)  {
       int n = getPos(str);
       if(n >= 0)
         return (std::get<0>(timer[n])==0)?(0.0):(std::get<2>(timer[n])/static_cast<double>(std::get<0>(timer[n])));
       return -1;
     }

     double total(const std::string& str )  {
       int n = getPos(str);
       if(n >= 0)
         return std::get<2>(timer[n]); 
       return -1;
     }

     void reset(const std::string& str) {
       int n = getPos(str);
       if(n < 0) {
         timer.push_back(std::make_tuple(0,0.0,0.0));
         id2pos[str] = timer.size()-1;
         return;
       }
       std::get<0>(timer[n])=0;
       std::get<1>(timer[n])=std::get<2>(timer[n])=0.0;
     }

     void reset_all() {
       for( TimeData& t: timer) { 
         std::get<0>(t)=0;
         std::get<1>(t)=std::get<2>(t)=0.0; 
       }
     }

     void print_elapsed(const std::string& str, std::ostream& out)
     {
       int n = getPos(str);
       if(n >= 0) 
         out<<" Elapsed time in " <<str <<": " <<getTime()-std::get<1>(timer[n]) <<"\n";    
       else
         out<<" Elapsed time in " <<str <<": Undefined Timer" <<"\n";  
     } 

     void print_average(const std::string& str, std::ostream& out)
     {
       int n = getPos(str);
       if(n >= 0) 
         out<<" Average time in " <<str <<": " <<((std::get<0>(timer[n])==0)?(0.0):(std::get<2>(timer[n])/static_cast<double>(std::get<0>(timer[n])))) <<"\n";  
       else
         out<<" Average time in " <<str <<": Undefined Timer" <<"\n";  
     }

     void print_total(const std::string& str, std::ostream& out)
     {
       int n = getPos(str);
       if(n >= 0)
         out<<" Total time in " <<str <<": " <<std::get<2>(timer[n]) <<"\n";
       else
         out<<" Total time in " <<str <<": Undefined Timer" <<"\n";  
 
     }

     void print_average_all(std::ostream& out)
     {
       for(std::map<std::string,int>::iterator it=id2pos.begin() ; it!=id2pos.end(); it++) {
         int n = it->second;
         out<<" Average time in " <<it->first <<": " <<((std::get<0>(timer[n])==0)?(0.0):(std::get<2>(timer[n])/static_cast<double>(std::get<0>(timer[n])))) <<"\n";
       }
     }

};

#endif  // myTimer 
