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

#############################################################################
# Depicts shape and color overlap between a 3D reference structure and
# a sets of 3D fit molecules. The molecules have to be pre-aligned.
# The first molecule is expected be the reference.
#############################################################################

import sys
from openeye import oechem
from openeye import oedepict
from openeye import oegrapheme
from openeye import oeshape
from openeye import oegraphsim


def main(argv=[__name__]):

    itf = oechem.OEInterface(InterfaceData)

    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    iname = itf.GetString("-in")
    oname = itf.GetString("-out")
    maxhits = itf.GetInt("-maxhits")
    depictsim = itf.GetBool("-depictsim")

    pagebypage = itf.GetBool("-pagebypage")

    # check input/output files

    ifs = oechem.oemolistream()
    if not ifs.open(iname):
        oechem.OEThrow.Fatal("Cannot open input molecule file!")

    ext = oechem.OEGetFileExtension(oname)
    if not pagebypage and not oedepict.OEIsRegisteredMultiPageImageFile(ext):
        oechem.OEThrow.Warning("Report will be generated into separate pages!")
        pagebypage = True

    refmol = oechem.OEMol()
    if not oechem.OEReadMolecule(ifs, refmol):
        oechem.OEThrow.Fatal("Cannot read reference molecule!")

    ropts = oedepict.OEReportOptions(3, 1)
    ropts.SetHeaderHeight(40.0)
    ropts.SetFooterHeight(20.0)
    report = oedepict.OEReport(ropts)

    cff = oeshape.OEColorForceField()
    cff.Init(oeshape.OEColorFFType_ImplicitMillsDean)
    cffdisplay = oegrapheme.OEColorForceFieldDisplay(cff)

    qopts = get_shape_query_display_options(depictsim)

    refmol_displays = dict()
    init_multi_query_displays(refmol_displays, refmol, cff, qopts)
    print("Shape overlaps will be generated for the reference with {} conformations.".format(len(refmol_displays)))

    fitmols = list()
    for fitidx, fitmol in enumerate(ifs.GetOEGraphMols()):
        if maxhits > 0 and fitidx >= maxhits:
            break
        fitmols.append(oechem.OEGraphMol(fitmol))
    maxhits = len(fitmols)

    depict_shape_color_graphsim_overlaps(report, refmol, refmol_displays, fitmols, depictsim)

    cffopts = oegrapheme.OEColorForceFieldLegendDisplayOptions(1, 6)
    for header in report.GetHeaders():
        oegrapheme.OEDrawColorForceFieldLegend(header, cffdisplay, cffopts)
        oedepict.OEDrawCurvedBorder(header, oedepict.OELightGreyPen, 10.0)

    font = oedepict.OEFont(oedepict.OEFontFamily_Default, oedepict.OEFontStyle_Default, 12,
                           oedepict.OEAlignment_Center, oechem.OEBlack)
    for idx, footer in enumerate(report.GetFooters()):
        oedepict.OEDrawTextToCenter(footer, "-" + str(idx + 1) + "-", font)

    if pagebypage:
        oedepict.OEWriteReportPageByPage(oname, report)
    else:
        oedepict.OEWriteReport(oname, report)

    return 0


def depict_shape_color_graphsim_overlaps(report, refmol, refmol_displays,
                                         fitmols, depictsim):
    """
    Depict shape, color, and 2D similarities.

    :type report: oedepict.OEReport
    :type refmol: oechem.OEMol
    :type refmol_displays: dict[string, oegrapheme.OEShapeQueryDisplay]
    :type fitmols: list[oechem.OEMol]
    :type depictsim: boolean
    """

    fptag, fptype, refmolfp, bondglyph = None, None, None, None
    if depictsim:
        fptag = oechem.OEGetTag("fpoverlap")
        fptype = oegraphsim.OEGetFPType(oegraphsim.OEFPType_Tree)
        print("Using fingerprint type %s" % fptype.GetFPTypeString())
        refmolfp = oegraphsim.OEFingerPrint()
        oegraphsim.OEMakeFP(refmolfp, refmol, fptype)
        fpcolorg = get_fingerprint_colorgradient(get_max_bond_self_similarity_score(refmol, fptype))
        bondglyph = ColorBondByOverlapScore(fpcolorg, fptag)

    sopts = get_shape_overlap_display_options()
    copts = get_color_overlap_display_options()

    ftableopts = get_fit_table_options(depictsim)
    rtableopts = get_ref_table_options()

    scorefont = oedepict.OEFont(oedepict.OEFontFamily_Default, oedepict.OEFontStyle_Bold, 9,
                                oedepict.OEAlignment_Center, oechem.OEBlack)

    tracer = oechem.OEConsoleProgressTracer()
    tracer.SetTask("Generating overlays")

    for fitidx, fitmol in enumerate(fitmols):

        tracer.SetProgress(fitidx, len(fitmols))

        if not oechem.OEHasSDData(fitmol, "ROCS_ShapeQuery"):
            warning = "Shape query reference is not available for molecule '{}'"
            oechem.OEThrow.Warning(warning.format(fitmol.GetTitle()))
            continue

        reftitle = oechem.OEGetSDData(fitmol, "ROCS_ShapeQuery")

        if reftitle not in refmol_displays:
            warning = "Shape query reference '{}'' is not valid for molecule '{}'"
            oechem.OEThrow.Warning(warning.format(reftitle, fitmol.GetTitle()))
            continue
        refdisp = refmol_displays[reftitle]

        cell = report.NewCell()
        fittable = oedepict.OEImageTable(cell, ftableopts)

        # title + score graph + query
        maintitle = "Hit: {}".format(fitmol.GetTitle())
        fittable.DrawText(fittable.GetCell(1, 1), maintitle)

        reftable = oedepict.OEImageTable(fittable.GetCell(2, 1), rtableopts)

        reftable.DrawText(reftable.GetCell(1, 1), "Rank: {}".format(fitidx+1))
        render_score(reftable.GetCell(2, 1), fitmol, "ROCS_TanimotoCombo", "Tanimoto Combo", scorefont)

        simscore = None if not depictsim else calc_fingerprint_similarity(refmol, refmolfp, fitmol, fptype, fptag)
        render_score_radial(reftable.GetCell(3, 1), fitmol, simscore)

        oegrapheme.OERenderShapeQuery(reftable.GetCell(4, 1), refdisp)
        reftable.DrawText(reftable.GetCell(5, 1), "query : {}".format(reftitle))

        odisp = oegrapheme.OEShapeOverlapDisplay(refdisp, fitmol, sopts, copts)

        # shape overlap
        render_score(fittable.GetHeaderCell(1), fitmol, "ROCS_ShapeTanimoto", "Shape Tanimoto", scorefont)
        oegrapheme.OERenderShapeOverlap(fittable.GetCell(2, 2), odisp)

        # color overlap
        render_score(fittable.GetHeaderCell(2), fitmol, "ROCS_ColorTanimoto", "Color Tanimoto", scorefont)
        oegrapheme.OERenderColorOverlap(fittable.GetCell(2, 3), odisp)

        # 2D similarity
        if depictsim:
            simtitle = "2D Graph Tanimoto = {:4.3f}".format(simscore)
            oedepict.OEDrawTextToCenter(fittable.GetHeaderCell(3), simtitle, scorefont)
            depict_molecule_similarity(fittable.GetCell(2, 4), fitmol, refdisp, bondglyph, fptag)

    tracer.Pop()


def init_multi_query_displays(query_map, refmol, cff, qopts):
    """
    Generated multiple shape display objects for each conformation
    of the reference molecule molecule.

    :type query_map: dict[string, oegrapheme.OEShapeQueryDisplay]
    :type refmol: oechem.OEMol
    :type cff: oeshape.OEColorForceField()
    :type qopts: oegrapheme.OEShapeQueryDisplayOptions()
    """

    for conf in refmol.GetConfs():
        title = conf.GetTitle()
        query_map[title] = oegrapheme.OEShapeQueryDisplay(conf, cff, qopts)

    return len(query_map)


def add_common_display_options(opts):

    opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)
    opts.SetAtomLabelFontScale(1.5)
    pen = oedepict.OEPen(oechem.OEBlack, oechem.OEBlack, oedepict.OEFill_Off, 1.5)
    opts.SetDefaultBondPen(pen)


def get_shape_query_display_options(depictsim):

    qopts = oegrapheme.OEShapeQueryDisplayOptions()
    add_common_display_options(qopts)
    arcpen = oedepict.OEPen(oedepict.OELightGreyPen)
    qopts.SetSurfaceArcFxn(oegrapheme.OEDefaultArcFxn(arcpen))
    if depictsim:
        qopts.SetDepictOrientation(oedepict.OEDepictOrientation_Vertical)
    else:
        qopts.SetDepictOrientation(oedepict.OEDepictOrientation_Square)
    qopts.SetBackgroundColor(oechem.OETransparentColor)
    return qopts


def get_shape_overlap_display_options():

    sopts = oegrapheme.OEShapeOverlapDisplayOptions()
    add_common_display_options(sopts)
    arcpen = oedepict.OEPen(oechem.OEGrey, oechem.OEGrey, oedepict.OEFill_Off, 1.0, 0x1111)
    sopts.SetQuerySurfaceArcFxn(oegrapheme.OEDefaultArcFxn(arcpen))
    sopts.SetOverlapColor(oechem.OEColor(110, 110, 190))
    sopts.SetOverlapDisplayStyle(oegrapheme.OEShapeOverlapDisplayStyle_PropertyCloud)
    sopts.SetBackgroundColor(oechem.OETransparentColor)
    return sopts


def get_color_overlap_display_options():

    copts = oegrapheme.OEColorOverlapDisplayOptions()
    add_common_display_options(copts)
    arcpen = oedepict.OEPen(oechem.OEGrey, oechem.OEGrey, oedepict.OEFill_Off, 1.0, 0x1111)
    copts.SetQuerySurfaceArcFxn(oegrapheme.OEDefaultArcFxn(arcpen))
    copts.SetBackgroundColor(oechem.OETransparentColor)
    return copts


def get_fit_table_options(depictsim):

    rows = 2
    cols = 4 if depictsim else 3
    tableopts = oedepict.OEImageTableOptions(rows, cols, oedepict.OEImageTableStyle_LightGrey)
    tableopts.SetHeader(True)
    tableopts.SetStubColumn(True)
    tableopts.SetRowHeights([10, 90])
    if depictsim:
        tableopts.SetColumnWidths([20, 30, 30, 30])
    else:
        tableopts.SetColumnWidths([18, 30, 30])
    tableopts.SetCellColor(oechem.OEWhite, False)
    font = oedepict.OEFont(oedepict.OEFontFamily_Default, oedepict.OEFontStyle_Bold, 9,
                           oedepict.OEAlignment_Center, oechem.OEBlack)
    tableopts.SetHeaderFont(font)
    return tableopts


def get_ref_table_options():

    tableopts = oedepict.OEImageTableOptions(5, 1, oedepict.OEImageTableStyle_NoStyle)
    tableopts.SetHeader(False)
    tableopts.SetStubColumn(False)
    tableopts.SetRowHeights([6, 6, 25, 60, 6])
    tableopts.SetMargins(0)
    font = oedepict.OEFont(oedepict.OEFontFamily_Default, oedepict.OEFontStyle_Default, 8,
                           oedepict.OEAlignment_Center, oechem.OEBlack)
    tableopts.SetCellFont(font)
    return tableopts


def get_score(mol, sdtag):

    if oechem.OEHasSDData(mol, sdtag):
        return float(oechem.OEGetSDData(mol, sdtag))
    return 0.0


def render_score_radial(image, mol, fpscore=None):

    sscore = max(min(get_score(mol, "ROCS_ShapeTanimoto"), 1.0), 0.0)
    cscore = max(min(get_score(mol, "ROCS_ColorTanimoto"), 1.0), 0.0)
    if sscore > 0.0 or cscore > 0.0:
        if fpscore is None:
            scores = oechem.OEDoubleVector([sscore, cscore])
        else:
            scores = oechem.OEDoubleVector([sscore, cscore, fpscore])
        oegrapheme.OEDrawROCSScores(image, scores)


def render_score(image, mol, sdtag, label, scorefont):

    score = get_score(mol, sdtag)
    if score == 0.0:
        return
    oedepict.OEDrawTextToCenter(image, label + " = " + str(score), scorefont)


#############################################################################
# 2D similarity depiction
#############################################################################

def get_max_bond_self_similarity_score(mol, fptype):

    rbonds = oechem.OEUIntArray(mol.GetMaxBondIdx())
    for match in oegraphsim.OEGetFPOverlap(mol, mol, fptype):
        for bond in match.GetPatternBonds():
            rbonds[bond.GetIdx()] += 1
    return max(rbonds)


def calc_fingerprint_similarity(refmol, fprefmole, fitmol, fptype, tag):

    fbonds = oechem.OEUIntArray(fitmol.GetMaxBondIdx())

    for match in oegraphsim.OEGetFPOverlap(refmol, fitmol, fptype):
        for bond in match.GetTargetBonds():
            fbonds[bond.GetIdx()] += 1

    for bond in fitmol.GetBonds():
        bond.SetData(tag, fbonds[bond.GetIdx()])

    # calculate 2D similarity score
    fpfitmol = oegraphsim.OEFingerPrint()
    oegraphsim.OEMakeFP(fpfitmol, fitmol, fptype)
    return oegraphsim.OETanimoto(fprefmole, fpfitmol)


def render_similarity_score(image, refmol, fitmol, fptype, label):

    score = oegraphsim.OETanimoto(refmol, fitmol)
    font = oedepict.OEFont(oedepict.OEFontFamily_Default, oedepict.OEFontStyle_Default, 9,
                           oedepict.OEAlignment_Center, oechem.OEBlack)
    oedepict.OEDrawTextToCenter(image, "{:s} = {:6.3f}".format(label, score), font)


def get_fingerprint_colorgradient(selfscore):
    colorg = oechem.OELinearColorGradient()
    colorg.AddStop(oechem.OEColorStop(0.0, oechem.OEPinkTint))
    colorg.AddStop(oechem.OEColorStop(1.0, oechem.OEYellow))
    colorg.AddStop(oechem.OEColorStop(selfscore, oechem.OEDarkGreen))
    return colorg


def depict_molecule_similarity(cell, mol3D, refdisp, bondglyph, fptag):
    """
    Aligns the molecule to the reference and depicts molecule similarity.

    :type mol3D: oechem.OEMol
    :type refdisp: oegrapheme.OEShapeQueryDisplay
    :type bondglyph: ColorBondByOverlapScore class
    :type fptag: integer
    """

    mol2D = oechem.OEGraphMol(mol3D)
    oegrapheme.OEPrepareAlignedDepictionFrom3D(mol2D, mol3D, refdisp)

    width, height, scale = cell.GetWidth(), cell.GetHeight(), oedepict.OEScale_AutoScale
    opts = oedepict.OE2DMolDisplayOptions(width, height, scale)
    opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
    opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)
    opts.SetScale(oegrapheme.OEGetMoleculeSurfaceScale(mol2D, opts))

    disp = oedepict.OE2DMolDisplay(mol2D, opts)
    oegrapheme.OEAddGlyph(disp, bondglyph, oechem.IsTrueBond())
    oedepict.OERenderMolecule(cell, disp, False)


class ColorBondByOverlapScore(oegrapheme.OEBondGlyphBase):
    def __init__(self, cg, tag):
        oegrapheme.OEBondGlyphBase.__init__(self)
        self.colorg = cg
        self.tag = tag

    def RenderGlyph(self, disp, bond):

        bdisp = disp.GetBondDisplay(bond)
        if bdisp is None or not bdisp.IsVisible():
            return False

        if not bond.HasData(self.tag):
            return False

        linewidth = disp.GetScale() / 3.0
        color = self.colorg.GetColorAt(bond.GetData(self.tag))
        pen = oedepict.OEPen(color, color, oedepict.OEFill_Off, linewidth)

        adispB = disp.GetAtomDisplay(bond.GetBgn())
        adispE = disp.GetAtomDisplay(bond.GetEnd())

        layer = disp.GetLayer(oedepict.OELayerPosition_Below)
        layer.DrawLine(adispB.GetCoords(), adispE.GetCoords(), pen)

        return True

    def ColorBondByOverlapScore(self):
        return ColorBondByOverlapScore(self.colorg, self.tag).__disown__()


#############################################################################
# INTERFACE
#############################################################################

InterfaceData = '''
!BRIEF [-in] <input> [-out] <output pdf>

!CATEGORY "input/output options:" 0

  !PARAMETER -in
    !ALIAS -i
    !TYPE string
    !REQUIRED true
    !KEYLESS 1
    !VISIBILITY simple
    !BRIEF Input molecule filename
    !DETAIL
         The first molecule in the file is expected to be the reference
         molecule
  !END

  !PARAMETER -out
    !ALIAS -o
    !TYPE string
    !REQUIRED true
    !KEYLESS 2
    !VISIBILITY simple
    !BRIEF Output image filename
  !END

!END

!CATEGORY "general options:" 1

  !PARAMETER -maxhits
    !ALIAS -mhits
    !TYPE int
    !REQUIRED false
    !DEFAULT 0
    !LEGAL_RANGE 0 500
    !VISIBILITY simple
    !BRIEF Maximum number of hits depicted
    !DETAIL
            The default of 0 means there is no limit.
  !END

  !PARAMETER -depictsim
    !ALIAS -sim
    !TYPE bool
    !REQUIRED false
    !DEFAULT false
    !VISIBILITY simple
    !BRIEF Calculate and depict 2D molecule similarity
  !END

!END

!CATEGORY "report options" 2

  !PARAMETER -pagebypage
    !ALIAS -p
    !TYPE bool
    !REQUIRED false
    !DEFAULT false
    !VISIBILITY simple
    !BRIEF Write individual numbered separate pages
  !END

!END

'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
