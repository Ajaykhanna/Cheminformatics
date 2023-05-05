import sys
from openeye import oechem
from openeye import oedepict


smiles = ["NC1=CC=C(C2=NON=C12)[N+](=O)[O-]", "CCN(CC)c1ccc2c(c1)oc-3cc(=O)c4ccccc4c3n2",
	  "COC1=CC2=C(C=C1)C(=CC(=O)O2)CC(=O)O"]


multi = oedepict.OEMultiPageImageFile(oedepict.OEPageOrientation_Landscape,
                                      oedepict.OEPageSize_US_Letter)
image = multi.NewPage()

opts = oedepict.OE2DMolDisplayOptions()

rows, cols = 2, 2
grid = oedepict.OEImageGrid(image, rows, cols)
grid.SetCellGap(20)
grid.SetMargins(20)
citer = grid.GetCells()

for smi in smiles:
    if not citer.IsValid():
        # go to next page
        image = multi.NewPage()
        grid = oedepict.OEImageGrid(image, rows, cols)
        grid.SetCellGap(20)
        grid.SetMargins(20)
        citer = grid.GetCells()

    cell = citer.Target()
    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, smi)
    oedepict.OEPrepareDepiction(mol)
    opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), oedepict.OEScale_AutoScale)
    disp = oedepict.OE2DMolDisplay(mol, opts)
    oedepict.OERenderMolecule(cell, disp)
    oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OERedPen))
    citer.Next()

oedepict.OEWriteMultiPageImage("MultiPage.png", multi)
