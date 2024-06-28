same resid as(within 5 of resid 1) # Will result in all molecules "without" broken bonds in a readion of 5 of resid 1
within 5 of resid 1 				# Will result in all molecules "with" broken bonds in a readion of 5 of resid 1
not within 4 of resid 1	  	# Selection of Solvent molecules away 5 angstrom from resid 1
resid 1 										# Selecting particular residue id
resid 1 to 50 							# Selecting bunch of resids
#----------------------------------------------------------------------------------------
name CA
resid 35
name CA and resname ALA
backbone
not protein
protein (backbone or name H)
name 'A 1'
name 'A *'
name "C.*"
mass < 5
numbonds = 2
abs(charge) > 1
x < 6 and x > 3
sqr(x-5)+sqr(y+4)+sqr(z) > sqr(5)
within 5 of name FE
protein within 5 of nucleic
same resname as (protein within 5 of nucleic)
protein sequence "C..C"
name eq $atomname 
protein and @myselection
index 0 to 95 # Selecting Atoms from 0 to 95
# There are two types of selection modes. The first is a keyword followed by a list of either values or a range of values. For example,
		name CA
# selects all atoms with the name CA (which could be a C  or a calcium);
		resname ALA PHE ASP
# selects all atoms in either alanine, phenylalanine, or asparagine;
		index 5
# selects the 6th atom (in the internal VMD numbering scheme).
# VMD can also do range selections, similar to X-PLOR's `:' notation:
		mass 5 to 11.5
# selects atoms with mass between 5 and 11.5 inclusive,
		resname ALA to CYS TYR
# selects atoms in alanine, arginine, asparagine, aspartic acid, cystine, and also tyrosine.
# The keyword selection works by checking each term on the list following the keyword. The term is either a single word (eg, name CA) or a range (eg resid 35 to 90).

# The method for determining the range checking is determined from the keyword data type; numeric comparisons are different than string comparisons.
The comparison should work as expected so that ``8'' is between ``1'' and ``11'' in a numeric context but not in a string one.
This may lead to some peculiar problems. Some keywords, such as segname, can take on string values but can also be used by some people as a number field.
Suppose someone labeled the segname field with the numbers 1 through 12 on the assumption that they are numbers. That person would be rather confused to find that segname 1 to 11 only returns two segments.
Also, strings will be converted (via atof()) to a number so if the string isn't a number, it will be given the value of 0. It is possible to force a search to be done in either a string or numeric context using the relational operator

# Selections can be combined with the boolean operators and and or, collected inside of parenthesis, and modified by not, as in
		(name CA or name CB) and mass 12 to 17
# which selects all atoms name CA or CB and have masses between 12 and 17 amu (this could be used to distinguish a C-alpha from a calcium). VMD has operator precedence similar to C so leaving the parentheis out of the previous expression, as in:
		name CA or name CB and mass 12 to 17
# actually selects all atoms named CA or those that are named CB and have the appropriate mass.
