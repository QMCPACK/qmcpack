#include "XMLWriterClass2.h"
#include <sstream>
#include <iomanip>


XMLAttribute::XMLAttribute(string name, string content)
{
  Name = name;
  Content = content;
}

XMLAttribute::XMLAttribute (const XMLAttribute &attr)
{ 
  Name = attr.Name; 
  Content = attr.Content; 
}


void
XMLAttribute::Write(string &out)
{
  stringstream str;
  str << Name << "=\"" << Content << "\"";
  out = str.str();
}


XMLElement::XMLElement(const XMLElement &elem)
{
  Attributes = elem.Attributes;
  Children   = elem.Children;
  Name       = elem.Name;
  Content    = elem.Content;
  Level      = elem.Level;
}

void 
XMLElement::Indent(ostream &out)
{
  for (int i=0; i<Level; i++)
    out << "  ";
}


void 
XMLElement::Write(ostream &out)
{
  Indent(out);
  out << "<" << Name;
  int lineCount = Name.size()+1;
  for (int i=0; i<Attributes.size(); i++) {
    out << " ";
    string str;
    Attributes[i]->Write(str);
    if (lineCount + str.size() + 2*Level > 73) {
      out << "\n";
      Indent (out);
      out << " ";
      lineCount = 1;
    }
    lineCount += str.size();
    out << str;
  }
  int numLines = 0;
  for (int i=0; i<Content.size(); i++)
    numLines += (Content[i] == '\n');

  if (Content == "" && Children.size()==0)
    out << "/>\n";
  else {
    out << ">";
    if (numLines > 0 || Children.size() > 0)
      out << "\n";
    out << Content;
    for (int i=0; i<Children.size(); i++)
      Children[i]->Write(out);
    if (numLines > 0 || Children.size()>0)
      Indent(out);
    out << "</" << Name << ">\n";
  }
}
    

void 
XMLElement::AddContent(vector<double> &data)
{
  stringstream str;
  str.setf(ios_base::scientific, ios_base::floatfield);
  str << setprecision(14);
  for (int i=0; i<data.size(); i++) {
    if ((i%3 == 0)) 
      Indent(str);
    str.width(22);
    str << data[i];
    if ((i%3) == 2)
      str << "\n";
  }
  if ((data.size()%3) != 0)
   str << "\n";
  Content += str.str();
}



bool
XMLWriterClass::StartDocument(string fname, string version,
			      string encoding, string standalone)
{
  Out.open (fname.c_str(),ofstream::out | ofstream::trunc);
  if (!Out.is_open())
    return false;
  if (version == "")
    version = "1.0";
  if (encoding == "")
    encoding = "UTF-8";
  
  Out << "<?xml version=\"" << version 
      << "\" encoding=\""   << encoding << "\"?>\n";

  return true;
}

bool
XMLWriterClass::EndDocument()
{
  if (Elements.size() <=0)
    return false;
  Elements[0]->Write(Out);
  return true;
}

bool
XMLWriterClass::StartElement(string name)
{
  int level = Elements.size();
  
  XMLElement* elem = new XMLElement(name, level);
  if (level > 0)
    Elements.back()->AddElement (elem);
  Elements.push_back(elem);
  return true;
}

bool
XMLWriterClass::EndElement()
{
  if (Elements.size() > 1)
    Elements.pop_back();
  return true;
}

bool
XMLWriterClass::FullEndElement()
{
  Elements.pop_back();
  return true;
}

bool
XMLWriterClass::WriteAttribute (string name, string content)
{
  XMLAttribute *attr = new XMLAttribute (name, content);
  Elements.back()->AddAttribute (attr);
  return true;
}

bool
XMLWriterClass::WriteAttribute (string name, double val, bool scientific)
{
  stringstream content;
  if (scientific) {
    content.setf(ios_base::scientific, ios_base::floatfield);
    content << setprecision(14);
  }
  content << val;
  XMLAttribute *attr = new XMLAttribute (name, content.str());
  Elements.back()->AddAttribute (attr);
  return true;
}

bool
XMLWriterClass::WriteAttribute (string name, int val)
{
  stringstream content;
  content << val;
  XMLAttribute *attr = new XMLAttribute (name, content.str());
  Elements.back()->AddAttribute (attr);
  return true;
}

bool
XMLWriterClass::WriteData(vector<double> data)
{
  Elements.back()->AddContent (data);
  return true;
}

bool
XMLWriterClass::WriteData(string data)
{
  Elements.back()->AddContent (data);
  return true;
}


bool
XMLWriterClass::WriteElement (string name, vector<double> data)
{
  int level = Elements.size();
  XMLElement *elem = new XMLElement (name, level);
  elem->AddContent (data);
  if (level > 0)
    Elements.back()->AddElement (elem);
  return true;
}

