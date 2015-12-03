##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  superstring                                                       #
#    Small collection of additional string utility functions.        #
#    Used infrequently.                                              #
#                                                                    #
#====================================================================#


'''
  superstring
    a collection of functions to manipulate strings

  general purpose
    next_visible_character
    remove_whitespace
    shrink_whitespace
    var2string
    string2array
    is_string
    stringmap
    stringbreak
    find_matching_brace
    remove_comment_lines
    contains_any
    contains_all

  C/C++
    find_preprocessor_end
    find_commend_block_end
    find_matching_cbrace

'''

from numpy import array
import sys
import string

#/////////////////////////////////////////////////
#///////        general purpose            ///////
#/////////////////////////////////////////////////


def next_visible_character(string,start,end):
    i = start
    character_visible = False
    while not character_visible and i<end:
        c = string[i]
        character_visible = c!=' ' and c!='\t' and c!='\n'
        i+=1
    #end while

    if character_visible:
        vis_char = c
        vis_loc  = i-1
    else:
        vis_char = ''
        vis_loc  = -1
    #end if

    return (vis_char,vis_loc)
#end def next_visible_character


def remove_whitespace(s):
    sr = s.replace('\n','').replace('\t','').replace(' ','')
    return
#end def remove_whitespace


def shrink_whitespace(si):
    sw = si.strip().replace('\n','')
    lst = sw.split(' ')
    s = ''
    for t in lst:
        if(t!=''):
            s += t+' '
        #end if
    #end for
    return s
#end def shrink_whitespace


def var2string(v):

    vt = type(v)

    nt = type(None)
    st = type(str(1))
    it = type(1)
    rt = type(1.0)
    at = type(array([[1]]))

    simple_set = set([nt,st,it,rt]) 

    s = ''
    if(vt == at):
        (n,m) = v.shape
        for i in range(n):
            for j in range(m):
                s += str(v[i,j]) + ' '
            #end for
            s += '\n'
        #end for
    elif(vt in simple_set):
        s = str(v)
    else:
        print 'ERROR: in var2string'
        print '   type '+str(vt)+' not implemented'
        sys.exit()
    #end if

    return s
#end def var2string

#string2val = lambda x: x.isalpha() and x \
#    or x.isdigit() and int(x) \
#    or x.isalnum() and x \
#    or len(set(string.punctuation).intersection(x)) == 1 and  x.count('.') == 1 and float(x) \
#    or x

def sbool(var):
    if var=='True':
        return True
    elif var=='False':
        return False
    else:
        return var
    #end if
#end def sbool

def is_bool(var):
    return var==True or var==False or var in ['True','False']
#end def is_bool

def is_int(var):
    try:
        int(var)
        return True
    except ValueError:
        return False
#end def is_float

def is_float(var):
    try:
        float(var)
        return True
    except ValueError:
        return False
#end def is_float

def is_array(var,type,delim=None):
    try:
        if isinstance(var,str):
            array(var.split(delim),type)
        else:
            array(var,type)
        #end if
        return True
    except ValueError:
        return False
#end def is_float_array

def string2val(s,delim=None):
    if is_bool(s):
        val = sbool(s)
    elif is_int(s):
        val = int(s)
    elif is_float(s):
        val = float(s)
    elif is_array(s,int,delim):
        val = array(s.split(delim),int)
    elif is_array(s,float,delim):
        val = array(s.split(delim),float)
    else:
        val = s
    #end if
    return val
#end def string2val



def string2array(string):
    ilst = string.strip().split(' ')
    lst = []
    for l in ilst:
        if(l.strip()!=''):
            lst.append(float(l))
        #end if
    #end for            
    return array(lst)
#end def string2array

def is_string(var):
    return type(var)==type("s")
#end def is_string


def stringmap(s):
    smap=[]
    quotes=set(['"',"'"])
    altquote={'"':"'","'":'"'}
    instr=False
    depth=0
    for i in range(len(s)):
        c=s[i]
        if not instr and c in quotes:
            instr=True
            lastquote=c
            depth=1
            direction=1
        elif instr and c in quotes:
            if c!=altquote[lastquote]:
                direction=-1
            #end if
            lastquote=c
            depth+=direction
        #end if
        smap+=[instr]
        if depth==0:
            instr=False
        #end if
    #end for
    return smap
#end def stringmap


def stringbreak(s,delimiter):
    strings=[]
    blocks=''
    strstart=s.startswith('"') or s.startswith("'")
    nblocks=0
    smap=[]
    quotes=set(['"',"'"])
    altquote={'"':"'","'":'"'}
    instr=False
    bstart=0
    depth=0
    for i in range(len(s)):
        c=s[i]
        if not instr and c in quotes:
            instr=True
            lastquote=c
            depth=1
            direction=1
            sstart=i
            bend=i
            if bend>0:
                blocks+=s[bstart:bend]+delimiter
            #end if
        elif instr and c in quotes:
            if c!=altquote[lastquote]:
                direction=-1
            #end if
            lastquote=c
            depth+=direction
        #end if
        #smap+=[instr]
        if depth==0 and instr:
            send=i+1
            strings+=[s[sstart:send]]
            instr=False
            bstart=send
        #end if
    #end for
    if not instr:
        bend=len(s)
        blocks+=s[bstart:bend]+delimiter
    #end if
    return strings,blocks,strstart
#end def stringbreak


def find_matching_brace(string,start,end):
    brace_dict = dict( [ ('(',')'), ('[',']'), ('{','}'), ('<','>') ] )
    left_brace  = string[start]
    right_brace = brace_dict[left_brace]
    found_match = False
    i = start + 1
    left_scope  = 0
    right_scope = 0
    while not found_match and i<end:
        if string[i]==left_brace:
            right_scope+=1
        elif string[i]==right_brace:
            found_match = right_scope==left_scope
            right_scope-=1
        #end if
        i+=1
    #end while
    if found_match:
        brace_loc = i-1
    else:
        brace_loc = -1
    #end if
    return brace_loc
#end def find_matching_brace


def find_matching_pair(s,pair,start=0,end=-1):
    if end==-1:
        end=len(s)
    #end if

    left  = pair[0]
    right = pair[1]

    llen=len(left)
    rlen=len(right)

    ileft  = s.find(left,start,end)
    iright = -1
    if ileft==-1:
        return ileft,iright
    else:
        i=ileft+llen
        left_scope  = 0
        right_scope = 0
        found_match = False
        failed = False
        while not found_match and i<end:
            nleft  = s.find(left,i,end)
            nright = s.find(right,i,end)
            if nleft!=-1 and nleft<nright:
                right_scope+=1
                i=nleft+llen
            elif nright!=-1:
                found_match = right_scope==left_scope
                right_scope-=1
                i=nright+rlen
            elif nright==-1:
                failed=True
                break
            #end if
        #end while
        if found_match:
            iright = i
        #end if
        if failed:
            ileft,iright=-1,-1
        #end if
    #end if
    return ileft,iright
#end def find_matching_pair


def remove_pair_sections(s,pair):
    sc=s
    ir=0
    n=0
    while ir!=-1 and n<10:
        il,ir = find_matching_pair(sc,pair)
        sc=sc.replace(sc[il:ir],'')
    #end while
    return sc
#end def


def remove_comment_lines(comment_char,s_in):
        lines = s_in.splitlines()
        s_out=''
        for l in lines:
            if not l.strip().startswith(comment_char):
                s_out=s_out+l+'\n'
            #end if
        #end if
        return s_out
#def remove_comment_lines


def remove_empty_lines(s):
    sr=''
    lines = s.splitlines()
    for l in lines:
        if l.strip()!='':
            sr+=l + '\n'
        #end if
    #end for
    return sr
#end def remove_empty_lines


def contains_any(str, set):
    for c in set:
        if c in str: return 1;
    return 0;
#end def contains_any

def contains_all(str, set):
    for c in set:
        if c not in str: return 0;
    return 1;
#end def contains_all


invalid_variable_name_chars=set('!"#$%&\'()*+,-./:;<=>?@[\\]^`{|}-\n\t ')
def valid_variable_name(s):
    return not contains_any(s,invalid_variable_name_chars)
#end def valid_variable_name


def split_delims(s,delims=['.','-','_']):
    sp = s
    for d in delims:
        sp = sp.replace(d,' ')
    #end for
    return sp.split()
#end def split_delims



#/////////////////////////////////////////////////
#///////             C/C++                 ///////
#/////////////////////////////////////////////////
def find_preprocessor_end(string,start,end):
    newline_loc = string.find('\n',start,end)
    prep_end = newline_loc
    line_continues = string[start:prep_end+1].rstrip(' \t\n').endswith('\\')
    continued_preprocessor = line_continues
    while line_continues:
        newline_loc = string.find('\n',prep_end+1,end)
        prep_end = newline_loc
        line_continues = string[start:prep_end+1].rstrip(' \t\n').endswith('\\')                        
    #end while

    return prep_end
#end def find_preprocessor_end


def find_comment_block_end(string,start,end):
    loc = string.find('*/',start,end)
    if loc!=-1:
        loc +=1
        #print 'fcbe',string[loc-1],string[loc]
    #end if
    return loc
#end def find_comment_block_end


def find_matching_cbrace(string,start,end,verbose=True):

    brace_dict = dict( [ ('(',')'), ('[',']'), ('{','}'), ('<','>') ] )

    left_brace  = string[start]
    right_brace = brace_dict[left_brace]

    found_match = False
    i = start + 1
    left_scope  = 0
    right_scope = 0

    in_comment_line  = False
    in_comment_block = False
    in_preprocessor  = False

    comment_block = False
    while not found_match and i<end:
##         if comment_block:
##             print 'fmb2',string[i],string[i+1]
##         #end if
        comment_block = False
        if string[i]=='#':            
            preprocessor_end = find_preprocessor_end(string,i,end)
            if preprocessor_end!=-1:
                i = preprocessor_end
            else:
                if verbose:
                    print 'ERROR: in find_matching_brace'
                    print '       end of preprocessor statement not found'
                #end if
                brace_loc = -1
            #end if
        elif string[i]=='/':
            comment_end = -1
            if string[i+1]=='/':
                comment_end = find_line_end(string,i,end)                
            elif string[i+1]=='*':
                comment_block = True
                comment_end = find_comment_block_end(string,i,end)
            else:
                comment_end = i #this is in the case of regular division
            #end if            
            if comment_end != -1:
                i = comment_end
            else:
                if verbose:
                    print 'ERROR: in find_matching_brace'
                    print '       comment mis-processed'
                #end if

                print string[i:end]
                print string[end+325]
                
                brace_loc = -1
            #end if
        elif string[i]==left_brace:
            right_scope+=1
        elif string[i]==right_brace:
            found_match = right_scope==left_scope
            right_scope-=1
        #end if

##         if comment_block:
##             print 'fmb1',string[i],string[i+1]
##         #end if
        
        i+=1
    #end while

    if found_match:
        brace_loc = i-1
    else:
        brace_loc = -1
    #end if

    return brace_loc
#end def find_matching_cbrace
