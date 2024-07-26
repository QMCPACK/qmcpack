#ifndef QMCPLUSPLUS_VALIDINPUTHELP_H
#define QMCPLUSPLUS_VALIDINPUTHELP_H

#define TEST_INPUT_ACCESSORS(ENUM_CLASS)  std::string_view operator[] (ENUM_CLASS val) const { \
    return xml[static_cast<std::size_t>(val)];                                            \
  }                                                                                       \
                                                                                          \
  auto begin() { return xml.begin(); }                                                    \
  auto end() {return xml.end(); }

#endif
