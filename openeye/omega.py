#!/usr/bin/env python
# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

import sys
from openeye import oechem
from openeye import oeomega


def main(argv=[__name__]):
    if len(argv) != 3:
        oechem.OEThrow.Usage("%s <infile> <outfile>" % argv[0])

    ifs = oechem.oemolistream()
    if not ifs.open(argv[1]):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % argv[1])

    ofs = oechem.oemolostream()
    if not ofs.open(argv[2]):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % argv[2])

    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.SetMaxConfs(1000)
    omega = oeomega.OEOmega(omegaOpts)

    for mol in ifs.GetOEMols():
        oechem.OEThrow.Info("Title: %s" % mol.GetTitle())
        if omega(mol):
            oechem.OEWriteMolecule(ofs, mol)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
