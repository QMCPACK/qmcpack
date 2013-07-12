#!/bin/env python
from xml.dom.minidom import *
import sys
from pylab import *
from numpy import *
doc = parse(sys.argv[1])

#add cuspCorrection attribute
detset=doc.getElementsByTagName("determinantset");
detset[0].setAttribute("cuspCorrection","yes")

#add cuspInfo attribute
dets=doc.getElementsByTagName("determinant");
a=dets[0].attributes["id"]
dets[0].setAttribute("cuspInfo",a.value+".cuspInfo.xml")

a=dets[1].attributes["id"]
dets[1].setAttribute("cuspInfo",a.value+".cuspInfo.xml")

print doc.toxml()
