import h5py

# Convert contents of an HDF file to C++ data structures.
# Small amounts of data can be stored in text form rather than binary.


def convert_to_std_vector(name, a):
    size = a.shape[0]

    #out = "  std::vector<double> " + name + " = {\n"
    out = "  " + name + " = {\n"
    for i in range(a.shape[0]):
        out += "    %12.10g, "%a[i]
        #if i%10 == 9:
        #    out += "\n  "
    out += "};\n";

    return out


# Example of desired output
#  deriv_records.resize(1,1);
#  double tmp[] = {1.1};
#  std::copy(tmp, tmp+1, deriv_records.begin());

def convert_to_ohmms_matrix(name, a):

    typename = "double"

    size = a.shape
    tmp_name = "tmp_%s"%name.replace(".","_")

    #out = "  Matrix<%s> %s(%d,%d);\n"%(typename, name, size[0], size[1])
    #out += "  %s %s[] = {\n"%(typename,tmp_name);
    out = "  %s.resize(%d, %d);\n"%(name, size[0], size[1])
    out += "  %s %s[] = {\n"%(typename,tmp_name);
    for i in range(a.shape[0]):
        out += "  // Row %d\n"%i
        for j in range(a.shape[1]):
            out += "%15.10g, "%a[i,j]
            if j%10 == 9:
                out += "\n  "
        out += "\n"
    out += "  };\n"

    out += "  std::copy(%s, %s+%d, %s.begin());\n"%(tmp_name, tmp_name, size[0]*size[1], name)

    return out

def convert_to_int_scalar(name, a):
    out = "  %s = %d;\n"%(name, a)
    return out

def convert_to_scalar(name, a):
    out = "  %s = %15.10g;\n"%(name,a[()])
    return out


inputs = h5py.File("qmc_short.s000.matrix_inputs.h5","r")
outputs = h5py.File("qmc_short.s000.linear_matrices.h5","r")

out_func = """
// Generated with convert_hdf_to_cpp.py
// clang-format off

void get_diamond_fill_data(qmcplusplus::FillData &fd)
{
"""


a = inputs["deriv_records"]

out_func += convert_to_int_scalar("fd.numSamples",a.shape[0])
out_func += convert_to_int_scalar("fd.numParam",a.shape[1])
out_func += convert_to_scalar("fd.sum_wgt",inputs["weight"])
out_func += convert_to_scalar("fd.sum_e_wgt",inputs["e_weight"])
out_func += convert_to_scalar("fd.sum_esq_wgt",inputs["e_sq_weight"])
out_func += convert_to_std_vector("fd.reweight",inputs["reweight"])
out_func += convert_to_std_vector("fd.energy_new",inputs["energy_new"])
out_func += convert_to_ohmms_matrix("fd.derivRecords", inputs["deriv_records"])
out_func += convert_to_ohmms_matrix("fd.HDerivRecords", inputs["H_deriv_records"])
out_func += convert_to_ohmms_matrix("fd.ovlp_gold", outputs["overlap"])
out_func += convert_to_ohmms_matrix("fd.ham_gold", outputs["Hamiltonian"])

out_func += "}\n"

print(out_func)

