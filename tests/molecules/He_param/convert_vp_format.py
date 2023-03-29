#!/usr/bin/env python3

import sys
import h5py
import numpy as np
import argparse

# Converts variational parameter files from HDF to text and back.
# The suffix of the input determines the conversion direction

# Sample text format
#version 1.0.0
#timestamp 2021-12-02 10:48:05 CST<
#jud_0 0.72053
# Format for complex parameters
#jud_0 0.72053 0.0

#Currently assumes that
# - No parameters are named "version" or "timestamp"
# - Parameter names have no spaces in the name


# Version 2 text format
#version 2.0.0
#timestamp 2021-12-02 10:48:05 CST<
#group jud
#  group name_value_lists
#    jud_0 0.72053
#  end_group
#end_group

# Sample HDF:
#HDF5 "he_opt.s009.vp.h5" {
#GROUP "/" {
#  GROUP "jud" {  // Version 2
#   GROUP "name_value_lists" {
#      DATASET "names" {
#         DATATYPE  H5T_STRING {
#            STRSIZE H5T_VARIABLE;
#            STRPAD H5T_STR_NULLTERM;
#            CSET H5T_CSET_ASCII;
#            CTYPE H5T_C_S1;
#         }
#         DATASPACE  SIMPLE { ( 4 ) / ( 4 ) }
#         DATA {
#         (0): "jud_0", "jud_1", "jud_2", "jud_3"
#         }
#      }
#      DATASET "values" {
#         DATATYPE  H5T_IEEE_F64LE
#         DATASPACE  SIMPLE { ( 4 ) / ( 4 ) }
#         DATA {
#         (0): 0.716782, 0.148293, -0.645633, -0.214129
#         }
#      }
#    } // Version 2
#   }
#   DATASET "timestamp" {
#      DATATYPE  H5T_STRING {
#         STRSIZE 23;
#         STRPAD H5T_STR_NULLTERM;
#         CSET H5T_CSET_ASCII;
#         CTYPE H5T_C_S1;
#      }
#      DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
#      DATA {
#      (0): "2021-12-02 10:28:39 CST"
#      }
#   }
#   DATASET "version" {
#      DATATYPE  H5T_STD_I32LE
#      DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
#      DATA {
#      (0): 1, 0, 0
#      }
#   }


class Group:
    def __init__(self, name):
        self.name = name
        self.groups = []
        # list of tuples (for complex values)
        self.name_value_pairs = []

class VP:
    def __init__(self):
        self.version = ""
        self.timestamp = ""
        # list of tuples (for complex values)
        self.name_value_pairs = []
        self.groups = []

    def version_as_array(self):
        array1 = self.version.split(".")
        return [int(a) for a in array1]

    def set_version_from_array(self, a):
        ar = [str(a1) for a1 in a]
        self.version = ".".join(ar)

    def major_version(self):
        return self.version_as_array()[0]

# Read text file

def parse_text_group(f, group_name):
    group = Group(group_name)
    while True:
        line = f.readline()
        # End of file
        if len(line) == 0 :
            print("Error: reached end of file inside group")
            break
        line = line.strip()
        if len(line) == 0 or line.startswith("#"):
            continue
        if line == "end_group":
            break
        elems = line.split(' ',1)
        name = elems[0]
        value = elems[1].strip()
        if name == 'group':
            new_group = parse_text_group(f, value)
            group.groups.append(new_group)

        if group_name == "name_value_lists":
            vals = value.split()
            if len(vals) == 1:
                val = (float(value), 0.0)
            if len(vals) == 2:
                r = float(vals[0])
                i = float(vals[1])
                val = (r,i)
            group.name_value_pairs.append( (name, val) )
    return group



def read_from_text(fname_in):
    vp = VP()
    with open(fname_in,'r') as f:
        while True:
            line = f.readline()
            # End of file
            if len(line) == 0 :
                break
            line = line.strip()
            elems = line.split(' ',1)
            name = elems[0]
            value = elems[1].strip()
            if len(line) == 0 or line.startswith("#"):
                continue
            if name == 'version':
                vp.version = value
                continue
            if name == 'timestamp':
                vp.timestamp = value
                continue
            if vp.major_version() == 1:
                vals = value.split()
                if len(vals) == 1:
                    val = (float(value), 0.0)
                if len(vals) == 2:
                    r = float(vals[0])
                    i = float(vals[1])
                    val = (r,i)
                vp.name_value_pairs.append( (name, val) )
            else:
                if name == 'group':
                    vp.groups.append(parse_text_group(f, value))

    return vp

# Write text file

def write_text_group(f, group, output_complex, level):
    pad = "  "*level
    if group.name == "name_value_lists":
        f.write(pad + "group " + group.name + "\n")
        for n,v in group.name_value_pairs:
            if output_complex:
                v_str = str(v[0]) + " " + str(v[1])
            else:
                v_str = str(v[0])
            line = n + " " + v_str + "\n"
            f.write(pad + "  " + line)
        f.write(pad + "end_group\n")
    else:
        f.write(pad + "group " + group.name + "\n")
        for g in group.groups:
            write_text_group(f, g, output_complex, level+1)
        f.write(pad + "end_group\n")


def write_to_text(vp, fname_out, output_complex=False):
    with open(fname_out, 'w') as f:
        f.write("version " + vp.version+"\n")
        f.write("timestamp " + vp.timestamp+"\n")
        if vp.major_version() == 1:
            for n,v in vp.name_value_pairs:
                if output_complex:
                    v_str = str(v[0]) + " " + str(v[1])
                else:
                    v_str = str(v[0])

                line = n + " " + v_str + "\n"
                f.write(line)
        else:
            for g in vp.groups:
                write_text_group(f, g, output_complex, 0)

#  Read HDF file

def read_hdf_group(name, hgroup):
    group = Group(name)
    for key,item in hgroup.items():
        if isinstance(item, h5py.Group):
            if key == "name_value_lists":
                nvl_group = Group("name_value_lists")
                g = hgroup["name_value_lists"]
                names = g["parameter_names"]
                values = g["parameter_values"]
                for n,v in zip(names, values):
                    name = n.decode("utf-8")
                    try:
                        val = (v[0], v[1])
                    except (TypeError,IndexError):
                        val = (v, 0.0)
                    nvl_group.name_value_pairs.append( (name, val) )
                group.groups.append(nvl_group)
            else:
                group = read_hdf_group(key,item)
                vp.groups.append(group)
    return group


def read_from_hdf(fname_in):
    f = h5py.File(fname_in,"r")
    vp = VP()
    vp.set_version_from_array(f["version"])

    vp.timestamp = f["timestamp"][0].decode("utf-8")

    if vp.major_version() == 1:
        g = f["name_value_lists"]
        names = g["parameter_names"]
        values = g["parameter_values"]
        for n,v in zip(names, values):
            name = n.decode("utf-8")
            try:
                val = (v[0], v[1])
            except (TypeError,IndexError):
                val = (v, 0.0)
            vp.name_value_pairs.append( (name, val) )
    else:
        for key,item in f.items():
            if isinstance(item, h5py.Group):
                group = read_hdf_group(key,item)
                vp.groups.append(group)

    return vp


# Write HDF file

def write_hdf_name_value_lists(g, group, output_complex):
    names = []
    values = []
    for n,v in group.name_value_pairs:
        names.append(n)
        if output_complex:
            values.append(v)
        else:
            values.append(v[0])

    g.create_dataset("parameter_names",data=np.array(names,dtype=object),dtype=h5py.string_dtype('ascii'))
    g.create_dataset("parameter_values",data=values)

def write_hdf_group(g, group, output_complex):
    g2 = g.create_group(group.name)
    if group.name == "name_value_lists":
        write_hdf_name_value_lists(g2, group, output_complex)
    for group2 in group.groups:
        write_hdf_group(g2, group2, output_complex)

def write_to_hdf(vp, fname_out, output_complex):


    f = h5py.File(fname_out,"w")
    f.create_dataset("timestamp",data=np.array([vp.timestamp],dtype=object),dtype=h5py.string_dtype('ascii'))
    f.create_dataset("version",data=vp.version_as_array())

    if vp.major_version() == 1:
        names = []
        values = []
        for n,v in vp.name_value_pairs:
            names.append(n)
            if output_complex:
                values.append(v)
            else:
                values.append(v[0])
        g = f.create_group("name_value_lists")
        g.create_dataset("parameter_names",data=np.array(names,dtype=object),dtype=h5py.string_dtype('ascii'))
        g.create_dataset("parameter_values",data=values)
    else:
        for g in vp.groups:
            write_hdf_group(f, g, output_complex)



def convert_from_text_to_hdf(fname_in, fname_out=None, output_complex=False):
    if not fname_out:
        fname_out = fname_in.replace(".txt",".h5")

    if fname_in == fname_out:
        print("Filenames identical, skipping h5 output")
        print("in = ",fname_in," out = ",fname_out)

    vp = read_from_text(fname_in)
    write_to_hdf(vp, fname_out, output_complex)


def convert_from_hdf_to_text(fname_in, fname_out=None, output_complex=False):
    if not fname_out:
        fname_out = fname_in.replace(".h5",".txt")

    if fname_in == fname_out:
        print("Filenames identical, skipping text output")
        print("in = ",fname_in," out = ",fname_out)

    vp = read_from_hdf(fname_in)
    write_to_text(vp, fname_out, output_complex)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert format of QMCPACK Variational Parameter files")
    parser.add_argument('input_file',help="Input file (HDF or text)")
    parser.add_argument('-o','--output',help="Output file name (default is input file name with suffix changed)")
    parser.add_argument('--complex',action='store_true',help="Output complex values")

    args = parser.parse_args()
    fname_in = args.input_file

    fname_out = None
    if args.output: fname_out = args.output

    if fname_in.endswith(".h5"):
        convert_from_hdf_to_text(fname_in, fname_out, args.complex)

    if fname_in.endswith(".txt"):
        convert_from_text_to_hdf(fname_in, fname_out, args.complex)

    if not fname_in.endswith(".h5") and not fname_in.endswith(".txt"):
        print("Expecting .h5 or .txt file suffix")
