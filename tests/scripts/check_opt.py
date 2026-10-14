#! /usr/bin/env python3


import os
import xml.etree.ElementTree as ET
from optparse import OptionParser


def exit_fail(msg=None):
    if msg != None:
        print(msg)
    # end if
    print("Test status: fail")
    exit(1)


# end def exit_fail


def exit_pass(msg=None):
    if msg != None:
        print(msg)
    # end if
    print("Test status: pass")
    exit(0)


# end def exit_pass


# Open the XML file and return coefficient values
def get_opt_coeff(file):
    tree = ET.parse(file)
    root = tree.getroot()
    j_coeff = []
    bf_coeff = []
    for elem in root.findall("wavefunction/jastrow/correlation/"):
        for i in elem.text.split():
            j_coeff.append(float(i))

    for elem in root.findall(
        "wavefunction/determinantset/backflow/transformation/correlation/"
    ):
        for i in elem.text.split():
            bf_coeff.append(float(i))

    return j_coeff, bf_coeff


# end def get_opt_coeff


passfail = {True: "pass", False: "fail"}


def run_opt_test(options):

    prefix_file = options.prefix + ".s" + str(options.series).zfill(3) + ".opt.xml"

    # Set default for the reference
    ref_file = (
        f"./qmc-ref/{prefix_file}" if options.ref is None else options.ref
    )

    if not os.path.exists(prefix_file):
        exit_fail("Test not found:" + prefix_file)

    if not os.path.exists(ref_file):
        exit_fail("Reference not found:" + ref_file)

    j_output, bf_output = get_opt_coeff(prefix_file)
    j_reference, bf_reference = get_opt_coeff(ref_file)

    if len(j_output) != len(j_reference):
        exit_fail(
            f"Number of coefficient in test({len(j_output)}) does not match with the reference({len(j_reference)})"
        )

    success = True
    j_tolerance = 1e-06
    bf_tolerance = 1e-05
    j_deviation = []
    bf_deviation = []

    for i in range(len(j_output)):
        j_deviation.append(abs(float(j_output[i]) - float(j_reference[i])))
        quant_success = j_deviation[i] <= j_tolerance
        if quant_success is False:
            success &= quant_success
        # end if
    # end for

    msg = f"\n  Testing Series: {options.series}\n"
    msg += f"   reference Jastrow coefficients   : {j_reference}\n"
    msg += f"   computed  Jastrow coefficients   : {j_output}\n"
    msg += f"   pass tolerance           : {j_tolerance: 12.6f}\n"
    msg += f"   deviation from reference : {j_deviation}\n"
    msg += f"   status of this test      :   {passfail[success]}\n"

    if bf_output or bf_reference:
        if len(bf_output) != len(bf_reference):
            exit_fail(
                f"Number of coefficient in test({len(bf_output)}) does not match with the reference({len(bf_reference)})"
            )

        for i in range(len(bf_output)):
            bf_deviation.append(abs(float(bf_output[i]) - float(bf_reference[i])))
            quant_success = bf_deviation[i] <= bf_tolerance
            if quant_success is False:
                success &= quant_success
            # end if
        # end for

        msg += f"\n  Testing Series: {options.series}\n"
        msg += f"   reference Backflow coefficients   : {bf_reference}\n"
        msg += f"   computed  Backflow coefficients   : {bf_output}\n"
        msg += f"   pass tolerance           : {bf_tolerance: 12.6f}\n"
        msg += f"   deviation from reference : {bf_deviation}\n"
        msg += f"   status of this test      :   {passfail[success]}\n"

    return success, msg


# end def run_opt_test


if __name__ == "__main__":
    parser = OptionParser(
        usage="usage: %prog [options]",
        add_help_option=False,
    )
    parser.add_option(
        "-h",
        "--help",
        dest="help",
        action="store_true",
        default=False,
        help="Print help information and exit (default=%default).",
    )
    parser.add_option(
        "-p",
        "--prefix",
        dest="prefix",
        default="qmc",
        help="Prefix for output files (default=%default).",
    )
    parser.add_option(
        "-s",
        "--series",
        dest="series",
        default="0",
        help="Output series to analyze (default=%default).",
    )
    parser.add_option(
        "-r",
        "--ref",
        dest="ref",
        help="Reference to check output files (default=./qmc-ref/$PREFIX).",
    )

    options, files_in = parser.parse_args()

    if options.help:
        print("\n" + parser.format_help().strip())
        exit()
    # end if

    success, msg = run_opt_test(options)

    if success:
        exit_pass(msg)
    else:
        exit_fail(msg)

# end if
