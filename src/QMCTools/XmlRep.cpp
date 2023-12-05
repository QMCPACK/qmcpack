/////////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois / NCSA Open Source License.
// See LICENSE file in top directory for details .
//
// Copyright ( c ) 2018 QMCPACK developers
//
// File developed by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
/////////////////////////////////////////////////////////////////////////////////////////

#include "XmlRep.h"
#include <ctype.h>
#include <stdio.h>
#include <string.h>
using namespace std;

string XmlStream::getTagName(const string& tagstr, tagType type) const
{
  string tagName;

  size_t spaceLoc = tagstr.find(" ");
  size_t closeLoc = tagstr.find('>');
  size_t slashLoc = string::npos;

  if (type == tagType::selfClosing)
  {
    slashLoc = tagstr.find("/");
  }

  int endChar = tagstr.size();
  if (spaceLoc != string::npos)
  {
    endChar = spaceLoc;
  }
  if (closeLoc < endChar && closeLoc != string::npos)
  {
    endChar = closeLoc;
  }
  if (slashLoc < endChar && slashLoc != string::npos)
  {
    endChar = slashLoc;
  }

  if (type == tagType::closing)
  {
    endChar -= 2;
    tagName = tagstr.substr(2, endChar);
  }
  else
  {
    endChar -= 1;
    tagName = tagstr.substr(1, endChar);
  }
  return tagName;
}

int XmlStream::startComment(long position, long length) const
{
  int isCommentStart    = 0;
  std::streampos curLoc = stream_->tellg();
  if ((length - position) > 4)
  {
    char buf[4];
    stream_->read(buf, 4);
    if (strncmp(buf, "<!--", 4) == 0)
    {
      isCommentStart = 1;
    }
  }
  stream_->seekg(curLoc);
  return isCommentStart;
}

int XmlStream::endComment(long position, long length) const
{
  int isCommentEnd      = 0;
  std::streampos curLoc = stream_->tellg();
  if ((length - position) > 3)
  {
    char buf[3];
    stream_->read(buf, 3);
    if (strncmp(buf, "-->", 3) == 0)
    {
      isCommentEnd = 1;
    }
  }
  stream_->seekg(curLoc);
  return isCommentEnd;
}

int XmlStream::checkForPOD(const XmlElement& before, const XmlElement& after) const
{
  int foundPOD          = 0;
  std::streampos curLoc = stream_->tellg();

  stream_->seekg(before.endLoc);

  long length   = after.startLoc - before.endLoc;
  long position = 0;
  int c;
  int inComment = 0;
  while (foundPOD == 0 && position < length)
  {
    c = stream_->get();
    if (!isspace(c))
    {
      //if (!isspace(c) && c != '\n') {
      stream_->unget();
      // if we're in a comment, check that we are not ending it
      if (inComment == 1)
      {
        if (endComment(position, length) == 1)
        {
          // position ourselves after the end comment tag
          stream_->seekg(3, std::ios_base::cur);
          position += 3;
          inComment = 0;
        }
      }
      else
      {
        // if we're not in a comment, check that we aren't starting one
        if (endComment(position, length) == 1)
        {
          // position ourselves after the comment tag
          stream_->seekg(4, std::ios_base::cur);
          position += 4;
          inComment = 1;
        }
        else
        {
          // we've found POD !!
          foundPOD = 1;
        }
      }
    }
    position++;
  }

  stream_->seekg(curLoc);
  return foundPOD;
}


int XmlStream::addNextTag()
{
  std::streampos start;
  std::streampos end;
  int twoprev = 0;
  int prev    = 0;
  int current = 0;

  int isProcessingInstruction = 0;
  int isComment               = 0;
  int isClosingTag            = 0;
  int isSelfClosing           = 0;

  int openCaretFound  = 0;
  int closeCaretFound = 0;

  int numSingleQuotes = 0;
  int numDoubleQuotes = 0;
  while ((current = stream_->get()) && current != EOF && closeCaretFound == 0)
  {
    if (current == '<')
    {
      if (prev != '\\')
      {
        // we found a live start string
        stream_->unget();
        start = stream_->tellg();
        stream_->get();
        openCaretFound = 1;
      }
      int one   = stream_->get();
      int two   = stream_->get();
      int three = stream_->get();
      if (one == '\?')
      {
        isProcessingInstruction = 1;
      }
      if (one == '!' && two == '-' && three == '-')
      {
        isComment = 1;
      }
      if (one == '/')
      {
        isClosingTag = 1;
      }
      stream_->unget();
      stream_->unget();
      stream_->unget();
    }
    if (openCaretFound == 1)
    {
      if (current == '\'')
      {
        numSingleQuotes++;
      }
      else if (current == '\"')
      {
        numDoubleQuotes++;
      }
      else if (current == '>')
      {
        // check that we aren't currently in a quoted section
        if (numSingleQuotes % 2 == 0 && numDoubleQuotes % 2 == 0)
        {
          // check that this close caret isn't escaped
          if (prev != '\\')
          {
            if (isComment == 1)
            {
              if (prev == '-' && twoprev == '-')
              {
                closeCaretFound = 1;
                end             = stream_->tellg();
              }
            }
            else
            {
              closeCaretFound = 1;
              end             = stream_->tellg();
              if (prev == '/')
              {
                isSelfClosing = 1;
              }
            }
          }
        }
      }
    }
    twoprev = prev;
    prev    = current;
  }

  if (current == EOF)
  {
    stream_->clear();
    stream_->seekg(0);
    return -1;
  }

  if (isProcessingInstruction == 0 && isComment == 0)
  {
    XmlElement elem;
    elem.startLoc = start;
    elem.endLoc   = end;

    if (isSelfClosing == 1)
    {
      elem.type = tagType::selfClosing;
    }
    else if (isClosingTag == 1)
    {
      elem.type = tagType::closing;
    }
    else
    {
      elem.type = tagType::opening;
      stream_->unget();
    }
    elem.name = getTagName(getTag(elem), elem.type);
    elements.push_back(elem);
    return 1;
  }
  if (isProcessingInstruction == 1 || isComment == 1)
  {
    stream_->unget();
    return 0;
  }
  return -1;
}

string XmlStream::getStreamSection(const std::streampos& start, const std::streampos& end) const
{
  // save current place in the stream
  std::streampos curLoc = stream_->tellg();

  std::string result(end - start, '\0');
  stream_->seekg(start);
  stream_->read(result.data(), end - start);

  // go back to current place in the stream
  stream_->seekg(curLoc);
  return result;
}

void XmlStream::findChildElements(int start, vector<int>& childIndices, int& podIndex) const
{
  int index = start + 1;
  int level = 1;

  while (index < elements.size() && level > 0)
  {
    const tagType tt = elements[index].type;
    if (tt == tagType::opening)
    {
      if (level == 1)
      {
        childIndices.push_back(index);
      }
      level++;
    }
    else if (tt == tagType::closing)
    {
      level--;
    }
    else if (tt == tagType::selfClosing)
    {
      if (level == 1)
      {
        childIndices.push_back(index);
      }
    }
    else if (tt == tagType::pod)
    {
      if (level == 1)
      {
        podIndex = index;
      }
    }
    else
    {
      //cout << "There is an error, in findChildElements and the XmlElement at index: " << index << " did not have a recognized tag type" << endl;
      exit(1);
    }
    index++;
  }
}

void XmlStream::listAll() const
{
  for (int i = 0; i < elements.size(); i++)
  {
    std::cout << i << "   ";
    tagType tt = elements[i].type;
    if (tt == tagType::opening)
    {
      std::cout << "opening, name = " << elements[i].name << std::endl;
    }
    else if (tt == tagType::closing)
    {
      std::cout << "closing, name = " << elements[i].name << std::endl;
    }
    else if (tt == tagType::selfClosing)
    {
      std::cout << "selfClosing, name = " << elements[i].name << std::endl;
    }
    else if (tt == tagType::pod)
    {
      std::cout << "POD" << std::endl;
    }
  }
}

XmlStream::XmlStream(std::istream* is) : stream_(is)
{
  // on first pass, try to find all of the tags and encode their names
  while (addNextTag() != -1) {}

  //now go through and look at space between live tags and see if there is POD there
  std::vector<XmlElement> tags;
  elements.swap(tags);

  elements.push_back(tags[0]);
  for (int i = 1; i < tags.size(); i++)
  {
    if (checkForPOD(tags[i - 1], tags[i]) == 1)
    {
      // add POD element
      XmlElement podElem;
      podElem.startLoc = tags[i - 1].endLoc;
      podElem.endLoc   = tags[i].startLoc;
      podElem.type     = tagType::pod;
      elements.push_back(podElem);
    }
    elements.push_back(tags[i]);
  }
}


string XmlStream::getTag(const XmlElement& e) const { return getStreamSection(e.startLoc, e.endLoc); }

string XmlStream::getTag(int i) const
{
  if (i >= elements.size())
  {
    cerr << "requested a tag index past the end of the vector" << endl;
    exit(1);
  }
  return getTag(elements[i]);
}

void XmlNode::readToString(std::string& s) const
{
  if (valueDeferred_)
  {
    std::streampos curLoc = stream_->tellg();

    s.resize(podEnd_ - podStart_);
    stream_->seekg(podStart_);
    stream_->read(s.data(), podEnd_ - podStart_);

    // go back to current place in the stream
    stream_->seekg(curLoc);

    // now strip out any xml comments
    std::stringstream ss;
    int commentStart = -1;
    int commentEnd   = 0;

    // if we find a comment, will put everything but the comments in ss
    while (s.find("<!--", commentEnd) != string::npos)
    {
      commentStart = s.find("<!--", commentEnd);
      ss << s.substr(commentEnd, commentStart - commentEnd);
      commentEnd = s.find("-->", commentStart);
    }
    // this means we didn't find a comment, so the string s is OK to return
    if (commentStart == -1)
    {
      return;
    }
    else
    {
      // we found comments and we're putting everything after the last comment in ss
      ss << s.substr(commentEnd + 3, s.size() - commentEnd - 3);
    }
    s = ss.str();
  }
  else
  {
    s = value_;
  }
}

// note, will return only the first child index that matches!
int XmlNode::getChildIndex(const string& childName, int strict) const
{
  int index = -1;
  for (int i = 0; i < children_.size(); i++)
  {
    if (children_[i].getName() == childName)
    {
      index = i;
    }
  }
  if (strict != 0 && index < 0)
  {
    cerr << "In XmlNode with name: " << name_ << ", could not find index for child with name: " << childName << endl;
    exit(1);
  }
  return index;
}

XmlNode& XmlNode::getChild(const string& name) { return getChild(getChildIndex(name, 1)); }

const XmlNode& XmlNode::getChild(const string& name) const { return getChild(getChildIndex(name, 1)); }

XmlNode& XmlNode::getChild(int i)
{
  if (i < 0 || i >= children_.size())
  {
    cerr << "Asked to get child node: " << i << ", but there are only " << getNumChildren() << "nodes" << endl;
    exit(1);
  }
  return children_[i];
}

const XmlNode& XmlNode::getChild(int i) const
{
  if (i < 0 || i >= children_.size())
  {
    cerr << "Asked to get child node: " << i << ", but there are only " << getNumChildren() << "nodes" << endl;
    exit(1);
  }
  return children_[i];
}

int XmlNode::getAttributeIndex(const string& attrName, int strict) const
{
  int index = -1;
  for (int i = 0; i < attributes_.size(); i++)
  {
    if (attributes_[i].first == attrName)
    {
      index = i;
    }
  }
  if (strict != 0 && index < 0)
  {
    cerr << "In XmlNode with name: " << name_ << ", could not find index for attribute with name: " << attrName << endl;
    exit(1);
  }
  return index;
}

string XmlNode::getAttributeName(int index) const { return attributes_[index].first; }

string XmlNode::getAttribute(int index) const
{
  if (index < 0 || index >= attributes_.size())
  {
    cerr << "in XmlNode with name: " << name_ << ", requested attribute with index " << index
         << ", but this index is not present." << endl;
    exit(1);
  }
  return attributes_[index].second;
}

string XmlNode::getAttribute(const string& name) const { return attributes_[getAttributeIndex(name, 1)].second; }

string XmlNode::getAttribute(const char* name) const
{
  string sname(name);
  return getAttribute(sname);
}

std::string XmlNode::getValue() const
{
  if (valueDeferred_)
  {
    std::string val;
    readToString(val);
    return val;
  }
  return value_;
}

int XmlNode::getValueSize() const
{
  if (valueDeferred_)
  {
    return (podEnd_ - podStart_);
  }
  return value_.size();
}

XmlNode& XmlNode::addChild(const XmlNode& nd)
{
  children_.push_back(nd);
  return children_.back();
}

XmlNode& XmlNode::addChild()
{
  XmlNode nd;
  children_.push_back(nd);
  return children_.back();
}

// removes whitespace at the beginning and end of a string
string XmlNode::trimWhitespace(const string& str) const
{
  size_t first = str.find_first_not_of(" \n\t\v\f\r");
  if (string::npos == first)
  {
    return string("");
  }
  size_t last = str.find_last_not_of(" \n\t\v\f\r");
  return str.substr(first, (last - first + 1));
}

void XmlNode::handleTagString(const XmlStream& xs, int loc)
{
  const XmlElement& elem = xs.elements[loc];
  name_                  = elem.name;
  string tagstr          = xs.getTag(elem);
  if (elem.type == tagType::selfClosing)
  {
    isSelfClosing_ = true;
  }
  else
  {
    isSelfClosing_ = false;
  }

  // take everything after the name
  string decliningstring = tagstr.substr(tagstr.find(name_) + name_.size() + 1);

  // remove trailing > or /> if appropriate and trim whitespace from left and right ends
  int numToRemove = 1;
  if (isSelfClosing_)
  {
    numToRemove++;
  }
  decliningstring = decliningstring.substr(0, decliningstring.size() - numToRemove);
  decliningstring = trimWhitespace(decliningstring);

  while (decliningstring.size() > 1)
  {
    attrpair att;
    getNextKeyVal(decliningstring, att.first, att.second);
    attributes_.push_back(att);
  }
}

void XmlNode::getNextKeyVal(string& contentstr, string& key, string& val) const
{
  size_t breakone = getPosNextLiveChar(contentstr, '=');

  key = contentstr.substr(0, breakone);
  key = trimWhitespace(key);
  //cout << "in getNextKeyVal, key = \'" << key << "\'" << endl;
  contentstr = contentstr.substr(breakone + 1);

  size_t firstquote = getPosNextLiveChar(contentstr, '\"');
  if (firstquote == string::npos)
  {
    firstquote = getPosNextLiveChar(contentstr, '\'');
  }
  contentstr         = contentstr.substr(firstquote + 1);
  size_t secondquote = getPosNextLiveChar(contentstr, '\"');
  if (secondquote == string::npos)
  {
    secondquote = getPosNextLiveChar(contentstr, '\'');
  }
  val = contentstr.substr(0, secondquote);
  val = trimWhitespace(val);

  contentstr = contentstr.substr(secondquote + 1);
}

size_t XmlNode::getPosNextLiveChar(const std::string& str, char c) const
{
  size_t index = string::npos;
  if (str[0] == c)
  {
    index = 0;
  }
  else
  {
    for (int i = 1; i < str.size(); i++)
    {
      if (str[i] == c && str[i - 1] != '\\')
      {
        index = i;
        break;
      }
    }
  }
  return index;
}

std::string XmlNode::getInStr(int is) const
{
  std::stringstream ss;
  for (int i = 0; i < is; i++)
  {
    ss << " ";
  }
  return ss.str();
}

void XmlNode::write(ostream& os, int indentLevel) const
{
  string str = getString(indentLevel);
  os << str;
}

string XmlNode::getString(int indentLevel) const
{
  stringstream ss;
  ss << getInStr(indentLevel);
  ss << "<" << name_;

  for (int i = 0; i < attributes_.size(); i++)
  {
    ss << " " << attributes_[i].first << "=\"" << attributes_[i].second << "\"";
  }
  if (isSelfClosing_ == 1)
  {
    ss << "/>" << endl;
    return ss.str();
  }
  else
  {
    ss << ">";
  }
  if (valInline_)
  {
    ss << getValue() << "</" << name_ << ">" << endl;
  }
  else
  {
    ss << endl;
    for (int i = 0; i < children_.size(); i++)
    {
      ss << children_[i].getString(indentLevel + 2);
    }

    if (getValueSize() > 0)
      ss << getInStr(indentLevel + 2) << getValue() << endl;
    ss << getInStr(indentLevel) << "</" << name_ << ">" << endl;
  }
  return ss.str();
}

XmlNode::XmlNode(const XmlNode& c) : stream_(c.stream_)
{
  isSelfClosing_ = c.isSelfClosing_;
  name_          = c.name_;
  attributes_    = c.attributes_;
  valInline_     = c.valInline_;
  value_         = c.value_;
  valueDeferred_ = c.valueDeferred_;
  podStart_      = c.podStart_;
  podEnd_        = c.podEnd_;
  children_      = c.children_;
}

XmlNode::XmlNode(istream* _stream, int start, bool deferValue) : stream_(_stream)
{
  XmlStream xl(_stream);
  //cout << "Finished creating stream object" << endl;
  //xl.listAll();
  createFromStream(xl, start, deferValue);
}

XmlNode::XmlNode(const XmlStream& xstream, std::istream* const _stream, int start, bool deferValue) : stream_(_stream)
{
  createFromStream(xstream, start, deferValue);
}

void XmlNode::createFromStream(const XmlStream& xstream, int start, bool deferValue)
{
  valueDeferred_ = deferValue;
  valInline_     = false;

  // this will populate the name and attributes
  handleTagString(xstream, start);
  tagType tt = xstream.elements[start].type;
  if (tt == tagType::selfClosing)
  {
    // if self closing, then there is not POD and we are at the end
    isSelfClosing_ = true;
  }
  else
  {
    // otherwise need to look for POD and subtags
    isSelfClosing_ = false;

    vector<int> childIndices;
    int podIndex = -1;
    xstream.findChildElements(start, childIndices, podIndex);

    // if no children, try to put and POD inline
    if (childIndices.size() == 0)
    {
      valInline_ = true;
    }
    else
    {
      valInline_ = false;
    }

    // deal with POD if it exists
    if (podIndex > 0)
    {
      podStart_ = xstream.elements[podIndex].startLoc;
      podEnd_   = xstream.elements[podIndex].endLoc;
      if (valueDeferred_ == false)
      {
        valueDeferred_ = true;
        readToString(value_);
        valueDeferred_ = false;
      }
    }

    // now sequentially create XmlNodes from subelements and add them to children vector
    for (int i = 0; i < childIndices.size(); i++)
    {
      XmlNode child(xstream, this->stream_, childIndices[i], valueDeferred_);
      children_.push_back(child);
    }
  }
}
