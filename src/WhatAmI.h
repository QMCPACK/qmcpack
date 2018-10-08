
#ifndef QMCPLUSPLUS_WHATAMI_H
#define QMCPLUSPLUS_WHATAMI_H

class WhatAmI
{
public:
  WhatAmI() : created_type_("nothing") {}
  WhatAmI(const std::string& created_type) : created_type_(created_type) {}
  virtual void setWhatAmI(const std::string& wai) { created_type_ = wai; }
  virtual const std::string& ask() const { return created_type_; }
private:
  std::string created_type_;
};

#endif
