#ifndef XML_WRITER_CLASS2
#define XML_WRITER_CLASS2

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>


class XMLAttribute
{
  std::string Name, Content;
public:
  void Write (std::ostream &out);
  void Write ( std::string  &out);
  
  XMLAttribute ( std::string name, std::string content);
  
  XMLAttribute ( std::string name, double val);
  XMLAttribute ( std::string name, int val);
  XMLAttribute (const XMLAttribute &attr);  
};

class XMLElement
{
  std::vector<XMLAttribute*> Attributes;
  std::vector<XMLElement*> Children;
  std::string Name, Content;
  int Level;
  void Indent(std::ostream &out);
public:
  void Write (std::ostream &out);

  inline int GetLevel() { return Level; }

  void AddElement (XMLElement *elem)
  { Children.push_back(elem); }

  void AddAttribute (XMLAttribute *attr) 
  { Attributes.push_back(attr); }

  void AddContent ( std::string content) 
  { Content += content; }

  void AddContent (std::vector<double> &data);

  XMLElement ( std::string name, int level=0)
  { Name = name; Level = level; }

  XMLElement (const XMLElement &elem);
};


class XMLWriterClass
{
private:
  std::vector<XMLElement*> Elements;
  void Write ();
  std::ofstream Out;
public:
  bool StartDocument( std::string fname, std::string version="",
		     std::string encoding="", std::string standalone="");
  bool EndDocument();
  bool StartElement ( std::string name);
  bool EndElement();
  bool FullEndElement();
  bool WriteAttribute ( std::string name, std::string content);
  bool WriteAttribute ( std::string name, double val, bool scientific=false);
  bool WriteAttribute ( std::string name, int val);
  bool WriteData(std::vector<double> data);
  bool WriteData ( std::string content);
  bool WriteElement( std::string name, std::vector<double> data);
  template<typename T> bool WriteData(T data);
};

template<typename T> bool
XMLWriterClass::WriteData(T data)
{
  std::stringstream content;
  content << " " << data << " ";
  Elements.back()->AddContent (content.str());
  return true;
}


#endif
