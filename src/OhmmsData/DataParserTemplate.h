//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#if defined(ENABLE_LIBXML2)
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#ifndef PUTDOMTODATA
#define PUTDOMTODATA bool put(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur);
#endif

struct DataParserBase {
  virtual bool put(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) = 0;
};

template<typename AType>
class DataParser : public DataParserBase
{
  AType& ref_;
public:
  DataParser(AType& a): ref_(a){ }
  bool put(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) {
    cur = cur->xmlChildrenNode;
    istringstream 
      stream((const char*)(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1)));
    stream >> ref_;
    return true;
  }
};

template<>
class DataParser<ProjectData>: public DataParserBase {
  ProjectData& ref_;
public:
  DataParser<ProjectData>(ProjectData& proj):ref_(proj){ }
  bool put(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) {
    ref_.m_title = (const char*)(xmlGetProp(cur, (const xmlChar *) "ID"));
    ref_.m_series = atoi((const char*)(xmlGetProp(cur, (const xmlChar *) "series")));
    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
      if ((!xmlStrcmp(cur->name, (const xmlChar *)"userid")) &&
	  (cur->ns == ns)){
	ref_.m_userid = (const char*)(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
      }
      if ((!xmlStrcmp(cur->name, (const xmlChar *)"host")) &&
	  (cur->ns == ns)){
	ref_.m_host = (const char*)(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
      }
      if ((!xmlStrcmp(cur->name, (const xmlChar *)"date")) &&
	  (cur->ns == ns)){
	ref_.m_date = (const char*)(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1));
      }
      cur = cur->next;
    }
    ref_.put();
    return true;
  }
};


template<>
class DataParser<MDRunData>: public DataParserBase {

  MDRunData& ref_;
public:
  DataParser<MDRunData>(MDRunData& proj):ref_(proj){ }

  bool put(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) {
    cur = cur->xmlChildrenNode;
    while (cur != NULL) {

      if ((!xmlStrcmp(cur->name, (const xmlChar *)"Loops")) && (cur->ns == ns)){
	
	ref_.blocks = atoi((const char*)(xmlGetProp(cur, (const xmlChar *) "blocks")));
	ref_.steps = atoi((const char*)(xmlGetProp(cur, (const xmlChar *) "steps")));

	xmlNodePtr subcur = cur->xmlChildrenNode;
	while(subcur!=NULL) {	
	  if ((!xmlStrcmp(subcur->name, (const xmlChar *)"timestep")) &&
	      (cur->ns == ns)){
	    ref_.dt = atof((const char*)(xmlNodeListGetString(doc, subcur->xmlChildrenNode, 1)));
	  }

	  ///\todo Will create PropertyRecord classes
	  if ((!xmlStrcmp(subcur->name, (const xmlChar *)"record")) &&
	      (cur->ns == ns)){
	    ref_.stride = atoi((const char*)(xmlGetProp(subcur, (const xmlChar *) "stride")));
	    ref_.outlevel= atoi((const char*)(xmlGetProp(subcur, (const xmlChar *) "level")));
	  }
	  subcur = subcur->next;
	}
      }
      if ((!xmlStrcmp(cur->name, (const xmlChar *)"ensemble")) &&
	  (cur->ns == ns)){
	ref_.mdtype = (const char*)(xmlGetProp(cur, (const xmlChar *) "name"));
	xmlNodePtr subcur = cur->xmlChildrenNode;
	while(subcur != NULL) {
	  if ((!xmlStrcmp(subcur->name, (const xmlChar *)"temperature")) &&
	      (subcur->ns == ns)){
	    ref_.temperature = atof((const char*)
				    (xmlNodeListGetString(doc, subcur->xmlChildrenNode, 1)));
	  }
	  if ((!xmlStrcmp(subcur->name, (const xmlChar *)"pressure")) &&
	      (subcur->ns == ns)){
	    ref_.pressure = 
	      atof((const char*)(xmlNodeListGetString(doc, subcur->xmlChildrenNode, 1)));
	  }
	  subcur = subcur->next;
	}
      }
      cur = cur->next;
    }
    return true;
  }
};

template<>
class DataParser<ParticleSets>: public DataParserBase {
  ParticleSets& ref_;
public:
  DataParser<ParticleSets>(ParticleSets& pts):ref_(pts){ }
  bool put(xmlDocPtr doc, xmlNsPtr ns, xmlNodePtr cur) {
    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
      if ((!xmlStrcmp(cur->name, (const xmlChar *)"ParticleSet")) &&
	  (cur->ns == ns)){
	//Particle_t* anewptcl = new Particle_t;
 	//DataParser<Particle_t> pparser(*anewptcl);
 	//pparser.put(doc,ns,cur);
 	ref_.add(OHMMS::createParticle(doc,ns,cur));
      }
      cur = cur->next;
    }
    ref_.put();
    return true;
  }
};


template<typename T>
class OhmmsElement: public OhmmsElementBase {

  T& ref_;

public:

  OhmmsElement(T& item): ref_(item){ }

  virtual void put(const string& s) const {
    istringstream stream(s);
    stream >> ref_;
  }

  bool get(ostream& ) const { return true;}
  bool put(istream& )  { return true;}
  bool put()  { return true;}
  bool action()  { return true;}
};

template<typename T>
void AddElement(map<string, 
		const OhmmsElementBase*>& emap, 
		T& aref, const string& apath) 
{

  emap.insert(make_pair(apath,new OhmmsElement<T>(aref)));
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
