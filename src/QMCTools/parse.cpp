#include "QMCTools/CasinoParser.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"

int main(int argc, char **argv) {

  if(argc<2) {
    std::cout << "Usage: casinon <filename> [-gridtype log|log0|linear -first ri -last rf -size npts]" << std::endl;
    std::cout << "Defaults : -gridtype log -first 1e-6 -last 100 -size 1001" << std::endl;
    return 1;
  }

  OHMMS::Controller->initialize(argc,argv);

  OhmmsInfo welcome(argc,argv,OHMMS::Controller->mycontext());

  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right,std::ios::adjustfield);
  std::cout.precision(12);

  QMCGaussianParserBase::init();

  CasinoParser parser(argc,argv);
  parser.parse(argv[1]);

  //parser.search(fin, "GEOMETRY");
  //parser.getGeometry(fin);
  //parser.search(fin, "BASIS");
  //parser.getGaussianCenters(fin);
  //parser.search(fin, "MULTIDETERMINANT");
  //parser.search(fin, "EIGENVECTOR");

  xmlDocPtr doc = xmlNewDoc((const xmlChar*)"1.0");
  xmlNodePtr wf_root = xmlNewNode(NULL, BAD_CAST "determinantset");
  xmlNewProp(wf_root,(const xmlChar*)"type",(const xmlChar*)"MolecularOrbital");
  xmlNewProp(wf_root,(const xmlChar*)"usegrid",(const xmlChar*)"yes");

  xmlNodePtr bset = parser.createBasisSet();
  xmlAddChild(wf_root,bset);

  xmlNodePtr slaterdet = parser.createDeterminantSet();
  xmlAddChild(wf_root,slaterdet);

  xmlDocSetRootElement(doc, wf_root);

  xmlXPathContextPtr m_context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr result
    = xmlXPathEvalExpression((const xmlChar*)"//atomicBasisSet",m_context);
  if(!xmlXPathNodeSetIsEmpty(result->nodesetval)) {
    for(int ic=0; ic<result->nodesetval->nodeNr; ic++) {
      xmlNodePtr cur = result->nodesetval->nodeTab[ic];
      parser.map2GridFunctors(cur);
    }
  }
  xmlXPathFreeObject(result);

  std::string fname = parser.Title+"."+parser.basisName+".xml";
  xmlSaveFormatFile(fname.c_str(),doc,1);
  xmlFreeDoc(doc);
  return 0;
}

/*
int main(int argc, char **argv) {

  char buffer[200];
  std::string _txt;
  std::string _data;
  double _t;
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right,std::ios::adjustfield);
  std::cout.precision(12);
  while(std::cin.getline( buffer, sizeof ( buffer ) ) ){
    std::istringstream stream(buffer);
    if(isdigit(buffer[1])) {
      while(stream>>_t) {
        std::cout << std::setw(21) << _t ;
      }
      std::cout << std::endl;
    } else {
      if(stream>>_t) {
        std::cout << "probably numbers " << buffer << std::endl;
      } else {
        std::cout << "probably statement " << buffer << std::endl;
      }
    }
  }
}
*/
