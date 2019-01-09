##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  vasp_analyzer.py                                                  #
#    Supports analysis of vasp output data including OUTCAR and      #
#    xml data (vasprun.xml).  All data (bands, forces, stress, etc.) #
#    are available in a structured format.                           #
#                                                                    #
#  Content summary:                                                  #
#    VaspAnalyzer                                                    #
#      SimulationAnalyzer class for VASP.                            #
#                                                                    #
#    OutcarData                                                      #
#      Reads/processes data from OUTCAR file.                        #
#                                                                    #
#    VaspLines                                                       #
#      Convenience class for line advancing read of OUTCAR file.     #
#                                                                    #
#    read_vxml                                                       #
#      Function to read vasprun.xml.  Returns nested object          #
#      representation of the file.                                   #
#                                                                    #
#    VXML                                                            #
#      Represents an xml node of vasprun.xml.                        #
#                                                                    #
#====================================================================#


import os
from numpy import array,zeros
from generic import obj
from simulation import Simulation,SimulationAnalyzer
from vasp_input import VaspInput,Incar
from developer import DevBase
from debug import *


# vasp xml reader classes/functions

class VXML(DevBase):
    basic_types = set('i v dimension field set time'.split())

    data_types = obj(int=int,string=str,float=float)
    
    def __init__(self,tag,attr=None):
        self._tag   = tag
        self._lines = []
        self._attr  = None
        self._value = None
        if attr!=None and len(attr)>0:
            self._attr = obj()
            tokens = attr.split('"')
            next=0
            for token in tokens:
                next+=1
                if len(token)>0 and token[-1]=='=':
                    name = token.strip()[:-1]
                    self._attr[name]=tokens[next].replace(' ','_').replace('-','_').lower()
                #end if
            #end for
        #end if
    #end def __init__

    def _is_empty(self):
        return len(self)-4==0
    #end def _is_empty

    def _add(self,new):
        tag = new._tag
        if not tag in self:
            self[tag] = new
        else:
            cur = self[tag]
            if 0 in cur:
                cur._append(new)
            else:
                coll = VXMLcoll(tag)
                coll._append(cur)
                coll._append(new)
                self[tag] = coll
            #end if
        #end if
    #end def _add

    def _parse(self):
        # rename sub objects if name is present
        for name in list(self.keys()):
            value = self[name]
            if isinstance(value,VXML):
                if value._attr!=None and 'name' in value._attr:
                    del self[name]
                    self[value._attr.name] = value
                    del value._attr.name
                elif isinstance(value,VXMLcoll):
                    for n in list(value.keys()):
                        v = value[n]
                        if isinstance(v,VXML):
                            if v._attr!=None and 'name' in v._attr:
                                del value[n]
                                self[v._attr.name] = v
                                del v._attr.name
                            #end if
                        #end if
                    #end for
                    if len(value)==0:
                        del self[name]
                    else:
                        value._reorder()
                    #end if
                #end if
            #end if
        #end for
        # have all sub objects parse (read fields, set _value)
        for v in self:
            if isinstance(v,VXML):
                v._parse()
            #end if
        #end for
        lines = self._lines
        if len(lines)>0:
            #fail = False
            #try:
            if self._tag=='array':
                self._parse_array(lines)
            else:
                self._parse_values(lines)
            #end if
            #except Exception,e:
            #    print '_parse failed!'
            #    print e
            #    print self
            #    print tag
            #    print attr
            #    print lines[0:min(10,len(lines))]
            #    print
            #    print
            #    fail = True
            ##end try
            #if fail:
            #    self.error('parse failed please see information above','read_vxml')
            ##end if
        #end if

        # if sub-objects resolve to a value, replace with that value
        for name in list(self.keys()):
            value = self[name]
            if isinstance(value,VXML) and value._value!=None: 
                self[name] = value._value
            #end if
        #end for

        # assign attributes
        if self._attr!=None:
            if 'type' in self._attr:
                del self._attr.type
            #end if
            self.transfer_from(self._attr)
        #end if

        return
    #end def _parse


    def _parse_values(self,lines):
        if len(lines)==1 and not '<' in lines[0]:
            self._value = readval(lines[0])
        else:
            arr = None
            for line in lines:
                start = line.startswith('<') and not line.startswith('</')
                end   = line.endswith('>')
                if start:
                    fline = line
                else:
                    fline += line
                #end if
                if end:
                    tokens = fline.replace('<','|').replace('>','|').split('|')
                    tn = tokens[1]
                    val = readval(tokens[2].strip())
                    if 'name=' in tn:
                        name = tn.split('name="',1)[1].split('"')[0].lower()
                        self[name] = val
                    elif arr is None:
                        arr = [val]
                    else:
                        arr.append(val)
                    #end if
                #end if
            #end for
            if arr!=None:
                self._value = array(arr)
            #end if
        #end if
    #end def parse_values


    def _parse_array(self,lines):
        #print 'parsing array'
        dims       = obj()
        fields     = obj()
        dim_counts = None
        field_list = []
        level      = -1
        set_dims   = False
        for line in lines:
            if line.startswith('<'):
                if line.startswith('<dimension'):
                    tokens = line.replace('<','|').replace('>','|').split('|')
                    tn = tokens[1]
                    dname = tokens[2].lower().replace(' ','_').replace('-','_')
                    if 'dim=' in tn:
                        d = int(tn.split('dim="',1)[1].split('"')[0])
                        dims[d] = dname
                    else:
                        dims.append(dname)
                    #end if
                elif line.startswith('<field'):
                    tokens = line.replace('<','|').replace('>','|').split('|')
                    tn = tokens[1]
                    fname = tokens[2].lower().replace(' ','_').replace('-','_')
                    if 'type=' in tn:
                        t = tn.split('type="',1)[1].split('"')[0]
                        if t in VXML.data_types:
                            dtype = VXML.data_types[t]
                        else:
                            self.error('field type {0} is unrecognized: {1}'.format(t,line))
                        #end if
                    else:
                        dtype = float
                    #end if
                    fields.append(obj(name=fname,dtype=dtype))
                elif line.startswith('<set'):
                    if not set_dims:
                        dims   = dims.list()
                        dims.reverse()
                        dims = tuple(dims)
                        dim_counts = zeros((len(dims),),dtype=int)
                        set_dims = True
                    #end if
                    level += 1
                    dim_counts[level]=0
                elif line.startswith('</set'):
                    level -= 1
                    if level!=-1:
                        dim_counts[level]+=1
                    #end if
                else:
                    self.error('array parsing failed\n unrecognized xml encountered: {0}'.format(line),'read_vxml')
                #end if
            else:
                dim_counts[level]+=1
                field_list.append(line.split())
            #end if
        #end for
        self.dims   = dims
        for findex,field in fields.iteritems():
            lst = []
            for field_vals in field_list:
                lst.append(field_vals[findex])
            #end for
            arr = array(lst,dtype=field.dtype).ravel()
            arr.shape = tuple(dim_counts)
            self[field.name] = arr
        #end for
        #print '  done'
    #end def _parse_array


    def _remove_empty(self):
        for n in list(self.keys()):
            v = self[n]
            if isinstance(v,VXML):
                v._remove_empty()
                if isinstance(v,VXMLcoll) and v._is_empty():
                    del self[n]
                #end if
            #end if
        #end for 
    #end def _remove_empty


    def _remove_hidden(self):
        del self._tag
        del self._attr
        del self._lines
        del self._value
        for v in self:
            if isinstance(v,VXML):
                v._remove_hidden()
            #end if
        #end for
    #end def _remove_hidden()
#end class VXML


class VXMLcoll(VXML):
    def _append(self,new):
        index = len(self)-4
        self[index]=new
    #end def _append

    def _reorder(self):
        n=0
        for key in sorted(self.keys()):
            value = self[key]
            if isinstance(value,VXML):
                del self[key]
                self[n]=value
                n+=1
            #end if
        #end for
    #end def _reorder
#end class VXMLcoll


booldict = dict(T=True,F=False)

def readval(val):
    fail = False
    split = False
    if isinstance(val,str):
        split = ' ' in val
    #end if
    if isinstance(val,list) or split:
        if split:
            val = val.split()
        #end if
        try:
            v = array(val,dtype=int)
        except:
            try:
                v = array(val,dtype=float)
            except:
                try:
                    v = array(val,dtype=str)
                except:
                    fail = True
                #end try
            #end try
        #end try
    elif val in booldict:
        v = booldict[val]
    else:
        try:
            v = int(val)
        except:
            try:
                v = float(val)
            except:
                v = val
            #end try
        #end try
    #end if
    if fail:
        VXML.class_error('failed to read value: "{0}"'.format(val),'read_vxml')
    #end if
    return v
#end def readval
            


def read_vxml(filepath):
    if not os.path.exists(filepath):
        VXML.class_error('file {0} does not exist'.format(filepath),'read_vxml')
    #end if
    #print 'read'
    contents = open(filepath,'r').read()
    #print 'replace'
    contents = contents.replace('<rc>',' ').replace('</rc>',' ')
    contents = contents.replace('<r>' ,' ').replace('</r>' ,' ')
    contents = contents.replace('<c>' ,' ').replace('</c>' ,' ')
    #print 'split lines'
    lines    = contents.splitlines()

    #print 'process lines'
    root = VXML('vasprun')
    stack = [root]
    cur = stack[0]
    for line in lines:
        ls = line.strip()
        if ls.startswith('</'):
            tag = ls[2:-1]
            if tag==cur._tag:
                stack.pop()
                cur = stack[-1]
                #print len(stack)*'  '+'end '+tag
            else:
                cur._lines.append(ls)
            #end if
        elif ls.startswith('<?'):
            None
        elif ls.startswith('<'):
            ta,rest = ls[1:].split('>',1)
            tokens = ta.split(' ',1)
            tag = tokens[0]
            if not tag in VXML.basic_types:
                if len(tokens)==1:
                    attr = None
                else:
                    attr = tokens[1].strip()
                #end if
                if ls.endswith('</{0}>'.format(tag)):
                    new = VXML(tag,attr)
                    new._lines.append(ls.replace('<','|').replace('>','|').split('|')[2])
                    cur._add(new)
                    #print len(stack)*'  '+'end '+tag
                else:
                    #print len(stack)*'  '+tag
                    new = VXML(tag,attr)
                    cur._add(new)
                    cur = new
                    stack.append(cur)
                #end if
            else:
                cur._lines.append(ls)
            #end if
        else:
            cur._lines.append(ls)
        #end if
    #end for

    if len(stack)!=1:
        VXML.class_error('read failed\nxml tree did not seem to close')
    #end if

    #print 'parse'
    root._parse()
    root._remove_empty()
    root._remove_hidden()

    #print 'done'
    return root
#end def read_vxml



# vasp outcar functions

class VaspLines(DevBase):
    def __init__(self,lines):
        self.pointer = 0
        self.lines   = lines
    #end def __init__

    def advance_line(self,amount):
        self.pointer += amount
        return self.lines[self.pointer]
    #end def advance_line

    def advance_token(self,token):
        psave = self.pointer
        for line in self.lines[self.pointer:]:
            if token in line:
                return line
            #end if
            self.pointer += 1
        #end while
        self.pointer = psave
        return None
    #end def advance

    def advance(self,amount):
        self.pointer += amount
    #end def advance

    def remainder(self):
        return self.lines[self.pointer:]
    #end def remainder

    def rewind(self,point=0):
        self.pointer = point
    #end def rewind

    def get_line(self,point=None):
        if point is None:
            point = self.pointer
        #end if
        return self.lines[point]
    #end def get_line

    def get_line_ahead(self,nahead):
        return self.lines[self.pointer+nahead]
    #end def get_line_ahead
#end class VaspLines


def read_outcar_header_values(vlines,odata):
    line = vlines.advance_token('TOTEN')
    odata.total_energy = float(line.split()[4])
    vlines.advance_token('energy without entropy')
#end def read_outcar_header_values


def read_outcar_core_potentials(vlines,odata):
    line = vlines.advance_token('the test charge radii are')
    odata.core_potential_radii = array(line.split()[5:],dtype=float)
    vlines.advance(2)
    n = 0
    cpots = []
    for line in vlines.remainder():
        ls = line.strip()
        n+=1
        if len(ls)==0:
            break
        #end if
        tokens = line.replace('-',' -').split()
        cpots.extend(tokens[1::2])
    #end for
    odata.core_potentials = array(cpots,dtype=float)
    vlines.advance(n)
#end def read_outcar_core_potentials


def read_outcar_fermi_energy(vlines,odata):
    line = vlines.advance_token('E-fermi')
    odata.Efermi = float(line.split()[2])
#end def read_outcar_fermi_energy


def read_outcar_bands(vlines,odata):
    bands = obj()
    line = vlines.advance_token('spin component')
    if line!=None:
        last_empty = True
        n  = 0
        for line in vlines.remainder():
            if len(line)>2:
                if line[1]=='s':
                    ns = int(line.split()[2])
                    spin = obj()
                    bands[ns] = spin
                elif line[1]=='k':
                    tokens = line.split()
                    nk = int(tokens[1])
                    kp = array(tokens[3:],dtype=float)
                    kpoint = obj(kpoint=kp,energies=[],occupations=[])
                    spin[nk]=kpoint
                elif line[2]=='b':
                    None
                else:
                    bnum,energy,occ = line.split()
                    kpoint.energies.append(float(energy))
                    kpoint.occupations.append(float(occ))
                #end if
                last_empty = False
            else:
                if last_empty:
                    break
                #end if
                last_empty = True
            #end if
            n+=1
        #end for
        vlines.advance(n)
    #end if
    for ns,spin in bands.iteritems():
        for nk,kpoint in spin.iteritems():
            kpoint.energies    = array(kpoint.energies,dtype=float)
            kpoint.occupations = array(kpoint.occupations,dtype=float)
        #end for
    #end for
    odata.bands = bands
#end def read_outcar_bands


def read_outcar_charge_mag(vlines,odata,token):
    ion   = obj(s=[],p=[],d=[],tot=[])
    total = obj()
    vlines.advance_token(token)
    vlines.advance(4)
    prev_end = False
    n=0
    for line in vlines.remainder():
        n+=1
        if prev_end:
            break
        #end if
        if line[0]=='-':
            prev_end = True
        else:
            vals = array(line.split()[1:],dtype=float)
            ion.s.append(vals[0])
            ion.p.append(vals[1])
            ion.d.append(vals[2])
            ion.tot.append(vals[3])
        #end if
    #end for
    for channel,vals in ion.iteritems():
        ion[channel] = array(vals,dtype=float)
    #end for
    vlines.advance(n)
    vals = array(line.split()[1:],dtype=float)
    total.s   = vals[0]
    total.p   = vals[1]
    total.d   = vals[2]
    total.tot = vals[3]
    return ion,total
#end def read_outcar_charge_mag


def read_outcar_total_charge(vlines,odata):
    ion,total = read_outcar_charge_mag(vlines,odata,'total charge ') # trailing space is important
    odata.ion_charge   = ion
    odata.total_charge = total
#end def read_outcar_total_charge


def read_outcar_magnetization(vlines,odata):
    ion,total = read_outcar_charge_mag(vlines,odata,'magnetization')
    odata.ion_magnetization   = ion
    odata.total_magnetization = total
#end def read_outcar_magnetization


def read_outcar_stress(vlines,odata):
    vlines.advance_token('FORCE on cell')
    line = vlines.advance_line(1)
    dirs = line.split()[1:]
    st       = array(vlines.advance_token('Total').split()[1:],dtype=float)
    st_kb    = array(vlines.advance_line(1).split()[2:],dtype=float)
    pressure = float(vlines.advance_line(1).split()[3])
    stress = obj()
    stress_kb = obj()
    for i in range(len(dirs)):
        d = dirs[i].lower()
        stress[d]    = st[i]
        stress_kb[d] = st_kb[i]
    #end for
    odata.stress    = stress
    odata.stress_kb = stress_kb
    odata.pressure  = pressure
#end def read_outcar_stress


def read_outcar_cell(vlines,odata):
    vlines.advance_token('VOLUME and BASIS')
    volume = float(vlines.advance_line(3).split()[-1])
    a1 = vlines.advance_line(2).split()[0:3]
    a2 = vlines.advance_line(1).split()[0:3]
    a3 = vlines.advance_line(1).split()[0:3]
    lattice_vectors = array([a1,a2,a3],dtype=float)
    odata.volume = volume
    odata.lattice_vectors = lattice_vectors
#end def read_outcar_cell


def read_outcar_position_force(vlines,odata):
    position    = []
    force       = []
    vlines.advance_token('POSITION')
    vlines.advance(2)
    prev_end = False
    for line in vlines.remainder():
        if prev_end:
            break
        #end if
        if line[1]=='-':
            prev_end = True
        else:
            tokens = line.split()
            position.append(tokens[0:3])
            force.append(tokens[3:6])
        #end if
    #end for
    total_drift = line.split()[2:5]
    odata.position = array(position,dtype=float)
    odata.force    = array(force,dtype=float)
    odata.total_drift = array(total_drift,dtype=float)
#end def read_outcar_position_force


def read_outcar_accounting(vlines,odata):
    time = obj()
    memory = obj()
    vlines.advance_token('General timing and accounting')
    vlines.advance(2)
    time.cpu     = float(vlines.advance_line(1).split()[-1])
    time.user    = float(vlines.advance_line(1).split()[-1])
    time.system  = float(vlines.advance_line(1).split()[-1])
    time.elapsed = float(vlines.advance_line(1).split()[-1])
    vlines.advance(1)
    memory.maximum = float(vlines.advance_line(1).split()[-1])
    memory.average = float(vlines.advance_line(1).split()[-1])
    odata.time   = time
    odata.memory = memory
#end def read_outcar_accounting




class OutcarData(DevBase):
    any_functions = [
        ('header_values'   , read_outcar_header_values  ),
        ]
    elast_functions = [
        ('core_potentials' , read_outcar_core_potentials),
        ('fermi_energy'    , read_outcar_fermi_energy   ),
        ('bands'           , read_outcar_bands          ),
        ('total_charge'    , read_outcar_total_charge   ),
        ('magnetization'   , read_outcar_magnetization  ),
        ('stress'          , read_outcar_stress         ),
        ('cell'            , read_outcar_cell           ),
        ('position_force'  , read_outcar_position_force ),
        ]
    ilast_functions = [
        ('accounting'      , read_outcar_accounting     ),
        ]

    read_outcar_functions = any_functions + elast_functions + ilast_functions

    def __init__(self,filepath=None,lines=None):
        if filepath!=None:
            if not os.path.exists(filepath):
                self.error('file {0} does not exist'.format(filepath))
            #end if
            f = open(filepath,'r')
            lines = f.read().splitlines()
            f.close()
        #end if
        self.vlines   = VaspLines(lines)
    #end def __init__


    def read(self,ilast=False,elast=False,all=True):
        ilast |= all
        elast |= all
        vlines = self.vlines
        del self.vlines
        read_functions = []
        read_functions.extend(self.any_functions)
        if elast:
            read_functions.extend(self.elast_functions)
            if ilast:
                read_functions.extend(self.ilast_functions)
            #end if
        #end if
        for quantity,read_function in read_functions:
            try:
                read_function(vlines,self)
            except:
                None
            #end try
        #end for
    #end def read
#end class OutcarData




# main analyzer class

class VaspAnalyzer(SimulationAnalyzer):
    def __init__(self,arg0=None,xml=False,analyze=False):
        path    = None
        prefix  = None
        incar   = None
        outcar  = None
        xml_file = None
        if isinstance(arg0,Simulation):
            sim = arg0
            file = sim.infile
            path   = sim.locdir
        elif arg0!=None:
            path,file = os.path.split(arg0)
        else:
            file = ''
            path = ''
            xml  = False
        #end if
        if len(file)>0:
            if file.endswith('INCAR'):
                incar  = file
                prefix = file.replace('INCAR','').strip()
            elif file.endswith('OUTCAR'):
                prefix = file.replace('OUTCAR','').strip()
            else:
                self.error('please provide the path to an INCAR or OUTCAR file')
            #end if
            incar   = prefix+'INCAR'
            outcar  = prefix+'OUTCAR'
            xml_file = prefix+'vasprun.xml'
            if prefix=='':
                prefix=None
            #end if
        #end if
        incar_file  = incar
        outcar_file = outcar
        incar = None
        incar_path = os.path.join(path,incar_file)
        neb = False
        if incar_file!=None and os.path.exists(incar_path):
            incar = Incar(incar_path)
            if 'images' in incar:
                neb = True
            #end if
        #end if
        self.info = obj(
            path         = path,
            prefix       = prefix,
            incar_file   = incar_file,
            outcar_file  = outcar_file,
            incar        = incar,
            xml_file     = xml_file,
            xml          = xml,
            neb          = neb
            )
        if analyze:
            self.analyze()
        #end if
    #end def __init__


    def analyze(self,outcar=None):
        if self.info.neb:
            self.neb_analyzers = obj()
            for i in range(self.info.incar.images):
                n = i+1
                self.neb_analyzers[n] = VaspAnalyzer(
                    arg0    =  os.path.join(self.info.path,str(n).zfill(2),'OUTCAR'),
                    xml     = self.info.xml,
                    analyze = True
                    )
            #end for
            return
        #end if
            
        if outcar is None and self.info.outcar_file!=None:
            outcar = os.path.join(self.info.path,self.info.outcar_file)
        #ned if
        if self.info.xml:
            self.xmldata = read_vxml(os.path.join(self.info.path,self.info.xml_file))
        #end if
        if outcar!=None:
            self.analyze_outcar(outcar)
        #end if
    #end def analyze


    def analyze_outcar(self,outcar):
        if not os.path.exists(outcar):
            self.error('outcar file {0} does not exist'.format(outcar))
        #end if
        oc = open(outcar,'r')
        lines = oc.read().splitlines()
        oc.close()
        del oc
        # gather initialization lines
        init = []
        n = 0
        for line in lines:
            if len(line)>0 and line[0]=='-' and 'Iteration' in line:
                break
            #end if
            init.append(line)
            n+=1
        #end for
        # gather lines for each iteration
        ion_steps = obj()
        for line in lines[n:]:
            if len(line)>0 and line[0]=='-' and 'Iteration' in line:
                iteration = []
                inum,enum = line.strip(' -Iteration)').split('(')
                inum = int(inum)
                enum = int(enum)
                if not inum in ion_steps:
                    ion_steps[inum] = obj()
                #end if
                ion_steps[inum][enum] = OutcarData(lines=iteration)
            #end if
            iteration.append(line)
        #end for
        del lines
        del n
        # read data from each iteration
        if len(ion_steps)>0:
            imax = array(ion_steps.keys(),dtype=int).max()
            for inum,ion_step in ion_steps.iteritems():
                ilast = inum==imax
                if len(ion_step)>0:
                    emax = array(ion_step.keys(),dtype=int).max()
                    for enum,elec_step in ion_step.iteritems():
                        elast = enum==emax
                        elec_step.read(ilast,elast,all=False)
                        if ilast and elast:
                            self.transfer_from(elec_step)
                        #end if
                    #end for
                #end if
            #end for
        #end if
        self.ion_steps = ion_steps
    #end def analyze_outcar
#end class VaspAnalyzer


