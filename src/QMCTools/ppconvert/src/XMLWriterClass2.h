#ifndef XML_WRITER_CLASS2
#define XML_WRITER_CLASS2

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

class XMLAttribute
{
  string Name, Content;
public:
  void Write (ostream &out);
  void Write (string  &out);
  
  XMLAttribute (string name, string content);
  
  XMLAttribute (string name, double val);
  XMLAttribute (string name, int val);
  XMLAttribute (const XMLAttribute &attr);  
};

class XMLElement
{
  vector<XMLAttribute*> Attributes;
  vector<XMLElement*> Children;
  string Name, Content;
  int Level;
  void Indent(ostream &out);
public:
  void Write (ostream &out);

  inline int GetLevel() { return Level; }

  void AddElement (XMLElement *elem)
  { Children.push_back(elem); }

  void AddAttribute (XMLAttribute *attr) 
  { Attributes.push_back(attr); }

  void AddContent (string content) 
  { Content += content; }

  void AddContent (vector<double> &data);

  XMLElement (string name, int level=0)
  { Name = name; Level = level; }

  XMLElement (const XMLElement &elem);
};


class XMLWriterClass
{
private:
  vector<XMLElement*> Elements;
  void Write ();
  ofstream Out;
public:
  bool StartDocument(string fname, string version="",
		     string encoding="", string standalone="");
  bool EndDocument();
  bool StartElement (string name);
  bool EndElement();
  bool FullEndElement();
  bool WriteAttribute (string name, string content);
  bool WriteAttribute (string name, double val, bool scientific=false);
  bool WriteAttribute (string name, int val);
  bool WriteData(vector<double> data);
  bool WriteData (string content);
  bool WriteElement(string name, vector<double> data);
  template<typename T> bool WriteData(T data);
};

template<typename T> bool
XMLWriterClass::WriteData(T data)
{
  stringstream content;
  content << " " << data << " ";
  Elements.back()->AddContent (content.str());
  return true;
}


#endif
