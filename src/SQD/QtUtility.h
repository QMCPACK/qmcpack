//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_SIMPLETREEBUILDER_H
#define OHMMS_SIMPLETREEBUILDER_H

// template<class NP>
// void
// SimpleTreeBuilder(QListViewItem *parentItem, NP parentNode){}
void SimpleTreeBuilder(QListViewItem *parentItem, xmlNodePtr parentNode)
{
  QListViewItem *thisItem = 0;
  xmlNodePtr cur=parentNode->children;
  while(cur != NULL)
  {
    QString cname((const char*)(cur->name));
    if(cname != "comment")
    {
      if(parentItem != 0)
      {
        if(cname == "text")
        {
          parentItem->setText(1,(const char*)(cur->content));
        }
        else
          if(cname == "PARAMETER")
          {
            thisItem = new QListViewItem(parentItem,thisItem);
            ///process attributes: type or format
            xmlAttrPtr att = cur->properties;
            while(att != NULL)
            {
              QString aname((const char*)(att->name));
              if(aname == "name")
              {
                thisItem->setText(0,(const char*)(att->children->content));
              }
              else
                if(aname == "InUnit")
                {
                  thisItem->setText(2,(const char*)(att->children->content));
                }
              att = att->next;
            }
            thisItem->setText(1,(const char*)(cur->children->content));
          }
          else
          {
            thisItem = new QListViewItem(parentItem,thisItem);
            thisItem->setText(0,cname);
            SimpleTreeBuilder(thisItem,cur);
            ///process attributes: type or format
            xmlAttrPtr att = cur->properties;
            while(att != NULL)
            {
              QString aname((const char*)(att->name));
              if(aname == "name")
              {
                thisItem->setText(1,(const char*)(att->children->content));
              }
              else
              {
                new QListViewItem(thisItem,aname,(const char*)(att->children->content));
              }
              att = att->next;
            }
          }
      }
    }
    cur = cur->next;
  }
}
/*
void SimpleTreeBuilder(QListViewItem *parentItem, const QDomElement& parentNode) {

  QListViewItem *thisItem = 0;
  QDomNode node = parentNode.firstChild();
  while ( !node.isNull() ) {
    if(node.isElement() && node.nodeName() != "comment") {
      if(parentItem != 0) {
 	if(node.isText()) {
 	  parentItem->setText(1,node.nodeValue());
	} else {
 	  if(node.nodeName() == "PARAMETER"){
	    thisItem = new QListViewItem(parentItem,thisItem);
	    QDomNamedNodeMap att_map=node.attributes();
	    for(int i=0; i<att_map.count(); i++) {
	      QDomNode att=att_map.item(i);
	      if(att.nodeName() == "name")
		thisItem->setText(0,att.nodeValue());
	      else if(att.nodeName() == "Unit")
		thisItem->setText(2,att.nodeValue());
	    }
	    if(!node.firstChild().isNull())
	      thisItem->setText(1,node.firstChild().nodeValue());
	  } else {
	    thisItem = new QListViewItem(parentItem,thisItem);
	    thisItem->setText(0,node.nodeName());
	    SimpleTreeBuilder(thisItem,node.toElement());
	    QDomNamedNodeMap att_map=node.attributes();
	    for(int i=0; i<att_map.count(); i++) {
	      QDomAttr att=att_map.item(i).toAttr(); //QDomNode att=att_map.item(i);
	      thisItem->setText(1,att.nodeValue());

	    }
	  }
	}
      }
    }
    node = node.nextSibling();
  }
  if(thisItem) thisItem->setOpen(true);
}
*/
#endif
