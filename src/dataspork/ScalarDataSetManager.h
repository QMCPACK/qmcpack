#ifndef QMCPLUSPLUS_SCALARDATAFILEMAGANER_H
#define QMCPLUSPLUS_SCALARDATAFILEMAGANER_H
#include <map>
///forward declation of ScalarObservable
class ScalarDataSet;

class ScalarDataSetManager
{

  std::vector<std::string> Observables;
  std::vector<std::string> Collectables;


  std::map<std::string,ScalarDataSet*> Files;

public:

  ScalarDataSetManager();
  ~ScalarDataSetManager();

  void addDataSet(const std::vector<std::string>& flist,
                  const std::string& gname, int first, int last);
  void addDataSet(const std::string& fname, int first, int last);

  void registerObservables(std::map<std::string, int>& olist);
  void registerCollectables(std::map<std::string, int>& olist);

  void write(std::ostream& os);

};
#endif
