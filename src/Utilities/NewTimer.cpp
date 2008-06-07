#include "Utilities/NewTimer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"

namespace qmcplusplus  {
  TimerManagerClass TimerManager;

  void TimerManagerClass::reset()
  {
    for (int i=0; i<TimerList.size(); i++) 
      TimerList[i]->reset();
  }

  void
  TimerManagerClass::print(Communicate* comm) 
  {
    // Sort all the timers by name, and summing the ones
    // with the same name into lists.
    TimerComparator comp;
    std::sort(TimerList.begin(), TimerList.end(), comp);
    std::vector<std::string> nameList;
    std::vector<double> timeList;
    std::vector<long>   callList;
    std::string lastName = "";
    int numDistinct = 0;

    for (int i=0; i<TimerList.size(); i++) {
      NewTimer &timer = *TimerList[i];
      if (timer.get_name() == lastName && lastName != "") {
	timeList[numDistinct-1]  += timer.get_total();
	callList[numDistinct-1] += timer.get_num_calls();
      }
      else {
	nameList.push_back(timer.get_name());
	timeList.push_back(timer.get_total());
	callList.push_back(timer.get_num_calls());
	lastName = timer.get_name();
	numDistinct++;
      }
    }
    // Now, we collect date from all nodes in the communicator, and
    // add it up.
    int numToSend = numDistinct;
    comm->allreduce(timeList);
    comm->allreduce(callList);
    //for (int i=0; i<numToSend; i++) {
    //  std::string myName = nameList[i];
    //  comm->bcast(myName);
    //  double myTime = 0.0;
    //  long myCalls = 0;
    //  for (int j=0; j<nameList.size(); j++) 
    //    if (nameList[j] == myName) {
    //      myTime  += timeList[j];
    //      myCalls += callList[j];
    //    }
    //  comm->allreduce(myTime);
    //  comm->allreduce(myCalls);
    //  if (comm->rank() == 0) {
    //    timeList[i] = myTime;
    //    callList[i] = myCalls;
    //  }
    //}
    
    bool omp_rank0 = true;
#ifdef ENABLE_OPENMP
    if (omp_get_thread_num() != 0)
      omp_rank0 = false;
#endif
    if (omp_rank0 && comm->rank() == 0) {
      fprintf (stderr, "Routine name                             Total time"
	       "      Num Calls     Time per call\n");
      for (int i=0; i<numDistinct; i++) {
	fprintf (stderr, "%-40s  %9.4f  %13ld  %16.9f\n", nameList[i].c_str(),
		 timeList[i], callList[i], timeList[i]/(double)callList[i]);
      }
    }
  }
  
}
