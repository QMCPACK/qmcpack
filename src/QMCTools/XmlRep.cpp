#include "XmlRep.h"
#include <ctype.h>
#include <stdio.h>
#include <string.h>
using namespace std;

string xmlStream::getTagName(const string& tagstr, tagType type) const {
  string tagName;
  
  size_t spaceLoc = tagstr.find(" ");
  size_t closeLoc = tagstr.find('>');
  size_t slashLoc = string::npos;

  if (type == tagType::selfClosing) {
    slashLoc = tagstr.find("/");
  }
  
  int endChar = tagstr.size();
  if (spaceLoc != string::npos) {
    endChar = spaceLoc;
  }
  if (closeLoc < endChar && closeLoc != string::npos) {
    endChar = closeLoc;
  }
  if (slashLoc < endChar && slashLoc != string::npos) {
    endChar = slashLoc;
  }

  if (type == tagType::closing) {
    endChar -= 2;
    tagName = tagstr.substr(2, endChar);  
  } else {
    endChar -= 1;
    tagName = tagstr.substr(1, endChar);  
  }
  return tagName;
}

int xmlStream::startComment(long position, long length) const {
  int isCommentStart = 0;
  std::streampos curLoc = stream->tellg();
  if ((length - position) > 4) {
    char buf[4];
    stream->read(buf, 4);
    if (strncmp(buf, "<!--", 4) == 0) {
      isCommentStart = 1;
    }
  }
  stream->seekg(curLoc);
  return isCommentStart;
}

int xmlStream::endComment(long position, long length) const {
  int isCommentEnd = 0;
  std::streampos curLoc = stream->tellg();
  if ((length - position) > 3) {
    char buf[3];
    stream->read(buf, 3);
    if (strncmp(buf, "-->", 3) == 0) {
      isCommentEnd = 1;
    }
  }
  stream->seekg(curLoc);
  return isCommentEnd;
}

int xmlStream::checkForPOD(const xmlElement& before, const xmlElement& after) const {
  int foundPOD = 0;
  std::streampos curLoc = stream->tellg();

  stream->seekg(before.endLoc);

  long length = after.startLoc - before.endLoc;
  long position = 0;
  char c;
  int inComment = 0;
  while (foundPOD == 0 && position < length) {
    c = stream->get();
    if (!isspace(c)) {
    //if (!isspace(c) && c != '\n') {
      stream->unget();
      // if we're in a comment, check that we are not ending it
      if (inComment == 1) {
	if (endComment(position, length) == 1) {
	  // position ourselves after the end comment tag
	  stream->seekg(3,std::ios_base::cur);
	  position+=3;
	  inComment = 0;
	}
      } else {
	// if we're not in a comment, check that we aren't starting one
	if (endComment(position,length) == 1) {
	  // position ourselves after the comment tag
	  stream->seekg(4,std::ios_base::cur);
	  position+=4;
	  inComment = 1;
	} else {
	  // we've found POD !!
	  foundPOD = 1;
	}
      }
    }
    position++;
  }
  
  stream->seekg(curLoc);
  return foundPOD;
}



int xmlStream::addNextTag() {
  std::streampos start;
  std::streampos end;
  char twoprev = '0';
  char prev = '0';
  char current = '0';
  
  int isProcessingInstruction = 0;
  int isComment = 0;
  int isClosingTag = 0;
  int isSelfClosing = 0;
  
  int openCaretFound = 0;
  int closeCaretFound = 0;
  
  int numSingleQuotes = 0;
  int numDoubleQuotes = 0;
  while ((current = stream->get()) && current != EOF && closeCaretFound == 0) {
    if (current == '<') {
      if (prev != '\\') {
	// we found a live start string
	stream->unget();
	start = stream->tellg();
	stream->get();
	openCaretFound = 1;
      }
      char one = stream->get();
      char two = stream->get();
      char three = stream->get();
      if (one == '\?') {
	isProcessingInstruction = 1;
      }
      if (one == '!' && two == '-' && three == '-') {
	isComment = 1;
      }
      if (one == '/') {
	isClosingTag = 1;
      }
      stream->unget(); stream->unget(); stream->unget();
    }
    if (openCaretFound == 1) {
      if (current == '\'') {
	numSingleQuotes++;
      } else if (current == '\"') {
	numDoubleQuotes++;
      } else if (current == '>') {
	// check that we aren't currently in a quoted section
	if (numSingleQuotes % 2 == 0 && numDoubleQuotes % 2 == 0) {
	  // check that this close caret isn't escaped
	  if (prev != '\\') {
	    if (isComment == 1) {
	      if (prev == '-' && twoprev == '-') {
		closeCaretFound = 1;
		end = stream->tellg();
	      }
	    } else {
	      closeCaretFound = 1;
	      end = stream->tellg();
	      if (prev == '/') {
		isSelfClosing = 1;
	      }
	    }
	  }
	}
      }
    }
    twoprev = prev;
    prev = current;
  }
  
  if (current == EOF) {
    stream->clear();
    stream->seekg(0);
    return -1;
  }
  
  if (isProcessingInstruction == 0 && isComment == 0) {
    xmlElement elem;
    elem.startLoc=start;
    elem.endLoc=end;

    if (isSelfClosing == 1) {
      elem.type = tagType::selfClosing;
    } else if (isClosingTag == 1) {
      elem.type = tagType::closing;
    } else {
      elem.type = tagType::opening;
    }
    elem.name = getTagName(getTag(elem), elem.type);
    elements.push_back(elem);
    return 1;
  }
  if (isProcessingInstruction == 1 || isComment == 1) {
    return 0;
  }
  return -1;
}

string xmlStream::getStreamSection(const std::streampos& start, const std::streampos& end) const {
  // save current place in the stream
  std::streampos curLoc = stream->tellg();

  char* buffer = new char[end - start];
  stream->seekg(start);
  stream->read(buffer,end-start);
  string result(buffer, end-start);
  delete[] buffer;
  
  // go back to current place in the stream
  stream->seekg(curLoc);
  return result;  
}

void xmlStream::findChildElements(int start, vector<int>& childIndices, int& podIndex) const {
  int index = start+1;
  int level = 1;

  while (index < elements.size() && level > 0) {
    const tagType tt = elements[index].type;
    if (tt == tagType::opening) {
      if (level == 1) {
	childIndices.push_back(index);
      }
      level++;
    } else if (tt == tagType::closing) {
      level--;
    } else if (tt == tagType::selfClosing) {
      if (level == 1) {
	childIndices.push_back(index);
      }
    } else if (tt == tagType::pod) {
      if (level == 1) {
	podIndex = index;
      }
    } else {
      //cout << "There is an error, in findChildElements and the xmlElement at index: " << index << " did not have a recognized tag type" << endl;
      exit(1);
    }
    index++;
  }
}
  
void xmlStream::listAll() const {
  for (int i = 0; i < elements.size(); i++) {
    std::cout << i << "   ";
    tagType tt = elements[i].type;
    if (tt == tagType::opening) {
      std::cout << "opening, name = " << elements[i].name << std::endl;
    } else if (tt == tagType::closing) {
      std::cout << "closing, name = " << elements[i].name << std::endl;
    } else if (tt == tagType::selfClosing) {
      std::cout << "selfClosing, name = " << elements[i].name << std::endl;
    } else if (tt == tagType::pod) {
      std::cout << "POD" << std::endl;
    }
  }
}

xmlStream::xmlStream(std::istream* is) : stream(is) {
  // on first pass, try to find all of the tags and encode their names
  while(addNextTag() != -1) { }
  
  //now go through and look at space between live tags and see if there is POD there
  std::vector<xmlElement> tags;
  elements.swap(tags);
  
  elements.push_back(tags[0]);
  for (int i = 1; i < tags.size(); i++) {
    if (checkForPOD(tags[i-1],tags[i]) == 1) {
      // add POD element
      xmlElement podElem;
      podElem.startLoc = tags[i-1].endLoc;
      podElem.endLoc = tags[i].startLoc;
      podElem.type = tagType::pod;
      elements.push_back(podElem);
    }
    elements.push_back(tags[i]);
  } 
}
  
  
string xmlStream::getTag(const xmlElement& e) const {
  return getStreamSection(e.startLoc, e.endLoc);
}

string xmlStream::getTag(int i) const {
  if (i >= elements.size()) {
    cout << "requested a tag index past the end of the vector" << endl;
    exit(1);
  }
  return getTag(elements[i]);
}

void xmlNode::readToString(std::string& s) const {
  if (valueDeferred) {
    std::streampos curLoc = stream->tellg();
    
    char* buffer = new char[podEnd - podStart];
    stream->seekg(podStart);
    stream->read(buffer,podEnd-podStart);
    
    s.assign(buffer,podEnd-podStart);
    delete[] buffer;

    // go back to current place in the stream
    stream->seekg(curLoc);
    
    // now strip out any xml comments
    std::stringstream ss;
    int inComment = 0;
    int commentStart = -1;
    int commentEnd = 0;

    // if we find a comment, will put everything but the comments in ss
    while (s.find("<!--",commentEnd) != string::npos) {
      commentStart = s.find("<!--",commentEnd);
      ss << s.substr(commentEnd, commentStart-commentEnd);
      commentEnd = s.find("-->",commentStart);
    }
    // this means we didn't find a comment, so the string s is OK to return
    if (commentStart == -1) {
      return;
    } else {
      // we found comments and we're putting everything after the last comment in ss
      ss << s.substr(commentEnd+3, s.size()-commentEnd-3);
    }
    s = ss.str();
  } else {
    s = value;
  }
}

// note, will return only the first child index that matches!
int xmlNode::getChildIndex(const string& childName, int strict) const {
  int index = -1;
  for (int i = 0; i < children.size(); i++) {
    if (children[i].getName() == childName) {
      index = i;
    }
  }
  if (strict != 0 && index < 0) {
    cout << "In xmlNode with name: " << name << ", could not find index for child with name: " << childName << endl;
    exit(1);
  }
  return index;
}

xmlNode& xmlNode::getChild(const string& name) {
  return getChild(getChildIndex(name, 1));
}		  

const xmlNode& xmlNode::getChild(const string& name) const {
  return getChild(getChildIndex(name, 1));
}		  

xmlNode& xmlNode::getChild(int i) {
  if (i < 0 || i >= children.size()) {
    cout << "Asked to get child node: " << i << ", but there are only " << getNumChildren() << "nodes" << endl;
    exit(1);
  }
  return children[i];
}

const xmlNode& xmlNode::getChild(int i) const {
  if (i < 0 || i >= children.size()) {
    cout << "Asked to get child node: " << i << ", but there are only " << getNumChildren() << "nodes" << endl;
    exit(1);
  }
  return children[i];
}

int xmlNode::getAttributeIndex(const string& attrName, int strict) const {
  int index = -1;
  for (int i = 0; i < attributes.size(); i++) {
    if (attributes[i].first == attrName) {
      index = i;
    }
  }
  if (strict != 0 && index < 0) {
    cout << "In xmlNode with name: " << name << ", could not find index for attribute with name: " << attrName << endl;
    exit(1);
  }
  return index;
}

string xmlNode::getAttributeName(int index) const {
  return attributes[index].first;
}

string xmlNode::getAttribute(int index) const {
  if (index < 0 || index >= attributes.size()) {
    cout << "in xmlNode with name: " << name << ", requested attribute with index " << index << ", but this index is not present." << endl;
    exit(1);
  } 
  return attributes[index].second;
}

string xmlNode::getAttribute(const string& name) const {
  return attributes[getAttributeIndex(name, 1)].second;
}

string xmlNode::getAttribute(const char* name) const {
  string sname(name);
  return getAttribute(sname);
}

std::string xmlNode::getValue() const {
  if (valueDeferred) {
    std::string val;
    readToString(val);
    return val;
  }
  return value;
}

int xmlNode::getValueSize() const {
  if (valueDeferred) {
    return (podEnd - podStart);
  }
  return value.size();
}

xmlNode& xmlNode::addChild(const xmlNode& nd) {
  children.push_back(nd);
  return children.back();
}

xmlNode& xmlNode::addChild() {
  xmlNode nd;
  children.push_back(nd);
  return children.back();
}

// removes whitespace at the beginning and end of a string
string xmlNode::trimWhitespace(const string& str) const {
  size_t first = str.find_first_not_of(" \n\t\v\f\r");
  if (string::npos == first) {
    return string("");
  }
  size_t last = str.find_last_not_of(" \n\t\v\f\r");
  return str.substr(first, (last - first + 1));
}

void xmlNode::handleTagString(const xmlStream& xs, int loc) {
  const xmlElement& elem = xs.elements[loc];
  name = elem.name;
  string tagstr = xs.getTag(elem);
  if (elem.type == tagType::selfClosing) {
    isSelfClosing = true;
  } else {
    isSelfClosing = false;
  }
  
  // take everything after the name
  string decliningstring = tagstr.substr(tagstr.find(name)+name.size()+1);
  
  // remove trailing > or /> if appropriate and trim whitespace from left and right ends
  int numToRemove = 1;
  if (isSelfClosing) {
    numToRemove++;
  }
  decliningstring = decliningstring.substr(0,decliningstring.size()-numToRemove);
  decliningstring = trimWhitespace(decliningstring);
  
  while(decliningstring.size() > 1) {
    attrpair att;
    getNextKeyVal(decliningstring, att.first, att.second);
    attributes.push_back(att);
  }
}

void xmlNode::getNextKeyVal(string& contentstr, string& key, string& val) const {
  size_t breakone = getPosNextLiveChar(contentstr, '=');

  key = contentstr.substr(0,breakone);
  key = trimWhitespace(key);
  //cout << "in getNextKeyVal, key = \'" << key << "\'" << endl;
  contentstr = contentstr.substr(breakone+1);
 
  size_t firstquote = getPosNextLiveChar(contentstr, '\"');
  if (firstquote == string::npos) {
    firstquote = getPosNextLiveChar(contentstr, '\'');
  }
  contentstr = contentstr.substr(firstquote+1);
  size_t secondquote = getPosNextLiveChar(contentstr, '\"');
  if (secondquote == string::npos) {
    secondquote = getPosNextLiveChar(contentstr, '\'');
  }
  val = contentstr.substr(0,secondquote);
  val = trimWhitespace(val);

  contentstr = contentstr.substr(secondquote+1);
}

size_t xmlNode::getPosNextLiveChar(const std::string& str, char c) const {
  size_t index = string::npos;
  if (str[0] == c) {
    index = 0;
  } else {
    for(int i = 1; i < str.size(); i++) {
      if (str[i] == c && str[i-1] != '\\') {
	index = i;
	break;
      }
    }
  }
  return index;
}

std::string xmlNode::getInStr(int is) const {
  std::stringstream ss;
  for (int i = 0; i < is; i++) {
    ss << " ";
  }
  return ss.str(); 
}

void xmlNode::write(ostream& os, int indentLevel) const {
  string str = getString(indentLevel);
  os << str;
}

string xmlNode::getString(int indentLevel) const {
  stringstream ss;
  ss << getInStr(indentLevel);
  ss << "<" << name;
  
  for (int i = 0; i < attributes.size(); i++) {
    ss << " " << attributes[i].first << "=\"" << attributes[i].second << "\"";
  }
  if (isSelfClosing == 1) {
    ss << "/>" << endl;
    return ss.str();
  } else {
    ss << ">";
  } 
  if (valInline) {
    ss << getValue() << "</" << name << ">" << endl;
  } else {
    ss << endl;
    for (int i = 0; i < children.size(); i++) {
      ss << children[i].getString(indentLevel+2);
    }
    
    if(getValueSize() > 0) ss << getInStr(indentLevel+2) << getValue() << endl;
    ss << getInStr(indentLevel) << "</" << name << ">" << endl;
  }
  return ss.str();
}

xmlNode::xmlNode(const xmlNode& c) : stream(c.stream) {
  isSelfClosing = c.isSelfClosing;
  name = c.name;
  attributes = c.attributes;
  valInline = c.valInline;
  value = c.value;
  valueDeferred = c.valueDeferred;
  podStart = c.podStart;
  podEnd = c.podEnd;
  children = c.children;
}

xmlNode::xmlNode(istream* _stream, int start, bool deferValue) : stream(_stream) {
  xmlStream xl(_stream);
  cout << "Finished creating stream object" << endl;
  //xl.listAll();
  createFromStream(xl, start, deferValue);
}

xmlNode::xmlNode(const xmlStream& xstream, std::istream* const _stream, int start, bool deferValue) : stream(_stream) {
  createFromStream(xstream, start, deferValue);
}

void xmlNode::createFromStream(const xmlStream& xstream, int start, bool deferValue) {
  valueDeferred = deferValue;

  // this will populate the name and attributes
  handleTagString(xstream, start);
  tagType tt = xstream.elements[start].type;
  if (tt == tagType::selfClosing) {
    // if self closing, then there is not POD and we are at the end
    isSelfClosing = true;
  } else {
    // otherwise need to look for POD and subtags
    isSelfClosing = false;

    vector<int> childIndices;
    int podIndex = -1;
    xstream.findChildElements(start, childIndices, podIndex);

    // if no children, try to put and POD inline
    if (childIndices.size() == 0) {
      valInline = true;
    } else {
      valInline = false;
    }

    // deal with POD if it exists
    if (podIndex > 0) {
      podStart = xstream.elements[podIndex].startLoc;
      podEnd = xstream.elements[podIndex].endLoc;
      if (valueDeferred == false) {
	valueDeferred = true;
	readToString(value);
	valueDeferred = false;
      }
    }

    // now sequentially create xmlNodes from subelements and add them to children vector
    for (int i = 0; i < childIndices.size(); i++) {
      xmlNode child(xstream, this->stream, childIndices[i], valueDeferred);
      children.push_back(child);
    }
  }
  
}
