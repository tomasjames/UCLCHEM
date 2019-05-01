#enter two reaction files to see all reactions that are one file but not the other

import math
import os
import string
import struct
import sys
import time
import fileinput
import itertools
from Functions import *


reactions1='inputFiles/umist12-ucledit.csv'
reactions2='inputFiles/kida.uva.2014.csv'
speciesFile = 'inputFiles/uclspeciesbasic.csv'

#differences are only relevant insofar as the missing reactions contain your species
nSpecies, speciesList = read_species_file(speciesFile)

make_capitals(reactions1)
make_capitals(reactions2)
make_capitals(speciesFile)


print '\nReading reactions'
nReactions1, reactions1 = read_reaction_file(reactions1, speciesList, 'UMIST')
nReactions2, reactions2 = read_reaction_file(reactions2,speciesList,'KIDA')

print "Reactions from UMIST not in KIDA"
for reaction1 in reactions1:
	match=False
	for reaction2 in reactions2:
		if set(reaction2.reactants)==set(reaction1.reactants):
			if set(reaction2.products)==set(reaction2.products):
				match=True
	if not match:
		if ("CH3OH" in reaction1.reactants) or ("CH3OH" in reaction1.products):
			print reaction1.reactants, "-->", reaction1.products, "alpha:", reaction1.alpha, "beta:", reaction1.beta, "gamma:", reaction1.gamma, "delta T:", reaction1.templow, "-", reaction1.temphigh

print "Reactions from KIDA not in UMIST"
for reaction1 in reactions2:
	match=False
	for reaction2 in reactions1:
		if set(reaction2.reactants)==set(reaction1.reactants):
			if set(reaction2.products)==set(reaction2.products):
				match=True
	if not match:
		if ("CH3OH" in reaction1.reactants) or ("CH3OH" in reaction1.products):
			print reaction1.reactants, "-->", reaction1.products, "alpha:", reaction1.alpha, "beta:", reaction1.beta, "gamma:", reaction1.gamma, "delta T:", reaction1.templow, "-", reaction1.temphigh
