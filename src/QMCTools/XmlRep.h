#ifndef XML_REP_H
#define XML_REP_H
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <utility>

// perhaps a project for the future would be to create a global xmlTree structure
// where everything here would be members of it.  It could handle the istream
// using RAII and could take a fair amount of complexity out of xmlNode

enum class tagType { opening, closing, selfClosing, pod };

// helper class, stores a few things about each xmlElement, whether
// it be a tag or POD
class xmlElement {
public:
  std::streampos startLoc;
  std::streampos endLoc;
  std::string name;
  tagType type;
  xmlElement() { };
  xmlElement(const xmlElement& e) {
    startLoc = e.startLoc;
    endLoc = e.endLoc;
    name = e.name;
    type = e.type;
  }
};    

// stream that will build up a sequential list of where all
// tags and POD are in the filestream that contains the xml file
// this is used to create an xmlNode structure
class xmlStream {
private:
  int addNextTag();
  std::string getTagName(const std::string& tag, tagType type) const;
  int checkForPOD(const xmlElement& before, const xmlElement& after) const;
  int startComment(long position, long length) const;
  int endComment(long position, long length) const;
  
public:
  xmlStream(std::istream* is);
  std::string getTag(const xmlElement& e) const;
  std::string getTag(int i) const;
  std::string getStreamSection(const std::streampos& start, const std::streampos& end) const;
  void findChildElements(int start, std::vector<int>& childIndices, int& podIndex) const;
  void listAll() const;
  
private:
  std::istream* stream;

public:
  std::vector<xmlElement> elements;
};


// structure that holds information about an xmlNode
// knows its name, attributes and holds its children xmlNodes
class xmlNode {
public:
  using attrpair = std::pair<std::string, std::string>;
private:
  // private information about how we store our state
  bool isSelfClosing;
  std::string name;
  std::vector<attrpair> attributes;
  bool valInline; // true for inline, false for putting the value on the next line
  std::string value;
  bool valueDeferred; // if true, need to look into the stream to get the value,
                      // if false it is stored in value
  std::streampos podStart;
  std::streampos podEnd;
  std::vector<xmlNode> children;
  std::istream* stream;

private:
  // private helper functions
  void readToString(std::string& s) const; // make s = to this xmlNode's POD
  void createFromStream(const xmlStream& stream, int start = 0, bool deferValue = true); 
  std::string getInStr(int is) const; // get whitespace to indent by is characters
  std::string trimWhitespace(const std::string& str) const;
  size_t getPosNextLiveChar(const std::string& str, char c) const;
  void getNextKeyVal(std::string& contentstr, std::string& key, std::string& val) const;
  void handleTagString(const xmlStream& xs, int loc); // populate name and attributes from loc of xmlStream

public:
  // public functions to interact with the xml data

  // deal with the name
  std::string getName() const { return name; }
  // be careful here, it is illegal to start a name with a number,
  // but this would allow you to screw that up...
  template<typename T> void setName(const T& val);

  // deal with attributes
  int getAttributeIndex(const std::string& name, int strict = 0) const;
  // note this is somewhat dangerous.
  // If the conversion to the specified type using a stringstream is not good, you will get garbage
  // You are protected in that if an attribute with that name is not found, the program will die
  template<typename T> void getAttribute(const std::string& name, T& result) const;
  template<typename T> void getAttribute(const char* name, T& result) const;
  std::string getAttribute(const std::string& name) const;
  std::string getAttribute(const char* name) const;
  std::string getAttribute(int index) const;
  std::string getAttributeName(int index) const;
  int getNumAttributes() const { return attributes.size(); }
  template<typename T, typename TT> void addAttribute(const T& n, const TT& v);  
  template<typename T> void addYesNoAttribute(const T& n, int i);

  // deal with the value (POD)
  template<typename T> void getValue(T& result) const;
  template<typename T> void getValue(std::vector<T>& result) const;
  std::string getValue() const;
  int getValueSize() const; // gives the number of characters in the value
  template<typename T> void setValue(const T& v);

  // deal with children nodes
  int getNumChildren() const { return children.size(); }
  int getChildIndex(const std::string& name, int strict = 0) const;
  xmlNode& getChild(const std::string& name);
  const xmlNode& getChild(const std::string& name) const;
  const xmlNode& getChild(int i) const;
  xmlNode& getChild(int i);
  xmlNode& addChild(const xmlNode& nd); // add nd to this's xmlNodes' list of children
  xmlNode& addChild(); // add a blank xmlNode to this xmlNode's list of children
  template<typename T, typename TT> void addParameterChild(const T& n, const TT& v);
  template<typename T> void addYesNoParameterChild(const T& n, int v);

  // constructors
  xmlNode() : stream(NULL) { isSelfClosing = false; valueDeferred = false; }
  xmlNode(std::istream* stream, int start = 0, bool deferValue = false);
  xmlNode(const xmlStream& xstream, std::istream* const stream, int start, bool deferValue = false);
  xmlNode(const xmlNode& c);

  // methods to write out the whole hierarchy
  void write(std::ostream& os, int indentLevel = 0) const;
  std::string getString(int indentLevel = 0) const;
};






// implemenation of templates here
template<typename T>
void xmlNode::setName(const T& val) {
  std::stringstream ss;
  ss << val;
  name = ss.str();
}

template<typename T>
void xmlNode::getAttribute(const std::string& name, T& result) const {
  std::string value = attributes[getAttributeIndex(name, 1)].second;
  std::stringstream ss(value);
  ss >> result;
}

template<typename T>
void xmlNode::getAttribute(const char* name, T& result) const {
  std::string sname(name);
  return getAttribute(sname, result);
}

template<typename T>
void xmlNode::getValue(T& result) const {
  if (valueDeferred) {
    std::string val;
    readToString(val);
    std::stringstream ss(val);
    ss >> result;
  } else {
    std::stringstream ss(value);
    ss >> result;
  }
}
template<typename T>
void xmlNode::getValue(std::vector<T>& result) const {
  std::stringstream ss;
  if (valueDeferred) {
    std::string val;
    readToString(val);
    ss << val;
  } else {
    ss << value;
  }
  T temp;
  while (ss >> temp) {
    result.push_back(temp);
  }
}

template<typename T>
void xmlNode::setValue(const T& v) {
  valueDeferred = false;
  std::stringstream ss;
  ss << v;
  value = ss.str();
}

template<typename T, typename TT>
void xmlNode::addAttribute(const T& n, const TT& v) {
  std::stringstream ss1;
  ss1 << n;
  std::stringstream ss2;
  ss2 << v;
  attrpair ap(ss1.str(), ss2.str());
  attributes.push_back(ap);
}

template<typename T>
void xmlNode::addYesNoAttribute(const T& n, int i) {
  std::stringstream ss;
  ss << n;
  attrpair ap;
  ap.first = ss.str();
  if (i == 0) {
    ap.second = "no";
  } else {
    ap.second = "yes";
  }
  attributes.push_back(ap);
}

template<typename T, typename TT>
void xmlNode::addParameterChild(const T& n, const TT& v) {
  xmlNode nd;
  nd.setName("parameter");
  nd.addAttribute("name", n);
  nd.valInline = true;
  nd.setValue(v);
    addChild(nd);
}

template<typename T>
void xmlNode::addYesNoParameterChild(const T& n, int v) {
  std::string vstr;
  if (v == 0) {
    vstr = "no";
  } else {
    vstr = "yes";
  }
  addParameterChild(n, vstr);
}


#endif
