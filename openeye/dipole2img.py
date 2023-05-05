#!/usr/bin/env python3
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

#############################################################################
# Calculates and visualizes the dipole moment of a molecule using property
# map
#############################################################################

import sys
from openeye import oechem
from openeye import oedepict
from openeye import oegrapheme


def main(argv=[__name__]):

    itf = oechem.OEInterface()
    oechem.OEConfigure(itf, InterfaceData)
    oedepict.OEConfigureImageWidth(itf, 400.0)
    oedepict.OEConfigureImageHeight(itf, 400.0)
    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    iname = itf.GetString("-in")
    oname = itf.GetString("-out")

    ifs = oechem.oemolistream()
    if not ifs.open(iname):
        oechem.OEThrow.Fatal("Cannot open %s input file!" % iname)

    ext = oechem.OEGetFileExtension(oname)
    if not oedepict.OEIsRegisteredImageFile(ext):
        oechem.OEThrow.Fatal("Unknown image type!")

    ofs = oechem.oeofstream()
    if not ofs.open(oname):
        oechem.EThrow.Fatal("Cannot open output file!")

    mol = oechem.OEMol()
    if not oechem.OEReadMolecule(ifs, mol):
        oechem.OEThrow.Fatal("Cannot read molecule from %s input file!" % iname)

    if mol.GetDimension() != 3:
        oechem.OEThrow.Fatal("3D coordinates are requires!")

    width, height = oedepict.OEGetImageWidth(itf), oedepict.OEGetImageHeight(itf)
    image = oedepict.OEImage(width, height)

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)

    stag = "dipole moment"
    itag = oechem.OEGetTag(stag)
    if calculate_dipole_moment(mol, itag):
        depict_molecule_with_dipole(image, mol, opts, stag)
    else:
        depict_molecule(image, mol, opts)

    oedepict.OEWriteImage(ofs, ext, image)


def depict_molecule_with_dipole(image, mol, opts, stag):
    """
    Depicts the molecule property map.

    :type image: oedepict.OEImageBase
    :type mol: oechem.OEMolBase
    :type opts: oedepict.OE2DMolDisplayOptions
    :type stag: string
    """

    oegrapheme.OEPrepareDepictionFrom3D(mol, True)

    opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
    opts.SetScale(oegrapheme.OEGetMoleculeSurfaceScale(mol, opts))

    disp = oedepict.OE2DMolDisplay(mol, opts)

    propmap = oegrapheme.OE2DPropMap(opts.GetBackgroundColor())
    propmap.SetLegendLocation(oegrapheme.OELegendLocation_Left)
    propmap.Render(disp, stag)

    oedepict.OERenderMolecule(image, disp)


def depict_molecule(image, mol, opts):

    oegrapheme.OEPrepareDepictionFrom3D(mol, True)
    disp = oedepict.OE2DMolDisplay(mol, opts)
    oedepict.OERenderMolecule(image, disp)


def calculate_dipole_moment(mol, itag):
    """
    Calculates dipole moment for each atom.

    :type mol: oechem.OEMolBase
    :type: itag: int
    """

    oechem.OEMMFFAtomTypes(mol)
    oechem.OEMMFF94PartialCharges(mol)

    ncrg = 0.0
    pcrg = 0.0
    ncen = [0.0, 0.0, 0.0]
    pcen = [0.0, 0.0, 0.0]
    for atom in mol.GetAtoms():
        charge = atom.GetPartialCharge()
        x, y, z = mol.GetCoords(atom)
        if charge < 0.0:
            ncen[0] -= charge * x
            ncen[1] -= charge * y
            ncen[2] -= charge * z
            ncrg -= charge
        elif charge > 0.0:
            pcen[0] += charge * x
            pcen[1] += charge * y
            pcen[2] += charge * z
            pcrg += charge

    for idx in range(0, 3):
        ncen[idx] = ncen[idx] / ncrg
        pcen[idx] = pcen[idx] / pcrg
    if pcrg > ncrg:
        pcrg = ncrg

    dpmag = 0.0
    for i in range(0, 3):
        dpmag += (ncen[i] - pcen[i]) * (ncen[i] - pcen[i])

    if abs(dpmag) < 0.001:
        # no dipole moment
        return False

    dpmag *= 4.80324 * pcrg

    dipdir = [0.0, 0.0, 0.0]
    for idx in range(0, 3):
        dipdir[idx] = 4.80324 * pcrg * ((pcen[idx] - ncen[idx]) / dpmag)

    dipcen = [0.0, 0.0, 0.0]
    for idx in range(0, 3):
        dipcen[idx] = 0.5 * (pcen[idx] + ncen[idx])

    dpip = oechem.OEFloatArray(mol.GetMaxAtomIdx())

    for atom in mol.GetAtoms():
        x, y, z = mol.GetCoords(atom)
        dpip[atom.GetIdx()] = sum([(x - dipcen[0]) * dipdir[0],
                                   (y - dipcen[1]) * dipdir[1],
                                   (z - dipcen[2]) * dipdir[2]])

    max_dpip = max(dpip)
    min_dpip = min(dpip)

    dpip = [-i / min_dpip if i < 0.0 else i for i in dpip]
    dpip = [+i / max_dpip if i > 0.0 else i for i in dpip]

    for atom, dipole in zip(mol.GetAtoms(), dpip):
        atom.SetData(itag, dipole)

    for bond in mol.GetBonds():
        avg = (bond.GetBgn().GetData(itag) + bond.GetEnd().GetData(itag)) / 2.0
        bond.SetData(itag, avg)

    return True


InterfaceData = '''
!CATEGORY "input/output options" 1

    !PARAMETER -in 1
      !ALIAS -i
      !TYPE string
      !REQUIRED true
      !KEYLESS 1
      !VISIBILITY simple
      !BRIEF Input molecule filename
    !END

    !PARAMETER -out 2
      !ALIAS -o
      !TYPE string
      !REQUIRED true
      !KEYLESS 2
      !VISIBILITY simple
      !BRIEF Output filename of the generated image
    !END

!END
'''


if __name__ == "__main__":
    sys.exit(main(sys.argv))
