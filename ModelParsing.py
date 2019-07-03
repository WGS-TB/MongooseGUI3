# This file contains functions for parsing the Excel files containing the networks
# Created by: Leonid Chindelevitch
# Last modified: November 30, 2016

from fractions import Fraction
import re
import xlrd
import libsbml
from ClassDefinitions import *
from Utilities import *
from FeaturePreparation import convertFormula
zero, one = Fraction(0), Fraction(1)
geneRE = re.compile('[0-9]*\.[0-9-]*\.[0-9-]*\.[0-9-]')
compRE = re.compile('\([a-zA-Z]\)|\[[a-zA-Z]\]')
coeffRE = re.compile('[0-9]+.?[0-9]*')
formRE = re.compile('FORMULA: [a-zA-Z0-9]*</p>')
chargeRE = re.compile('CHARGE: [-a-zA-Z0-9]*</p>')

def parseSBML(filename, specialEnding = '', Cobra = False):
    # Parses an SBML file containing the network description
    document = libsbml.SBMLReader().readSBMLFromFile(filename)
    model = document.getModel()
    version, level = model.getVersion(), model.getLevel()
    modelName, modelId = model.getName(), model.getId()
    modelDescription = {'Source File': filename, 'ID': modelId, 'Name': modelName}
    allSpecies = model.getListOfSpecies()
    if (level == 2 and version in [2, 3, 4]):
        allSpeciesTypes = model.getListOfSpeciesTypes()
    else:
        allSpeciesTypes = None
    allCompartments = [comp.getId() for comp in model.getListOfCompartments()]
    extraCompartment = findExactlyOne(allCompartments, ['extra', 'outside'])
    if extraCompartment is not None:
        extraCompartment = allCompartments[extraCompartment]
    (species, metabolites) = createSpeciesList(allSpecies, allSpeciesTypes, extraCompartment, Cobra, specialEnding)
    allReactions = model.getListOfReactions()
    network = createNetwork(allReactions, species, metabolites, Cobra)
    network.description = modelDescription
    network.createMatrices()
    return network

def createSpeciesList(speciesObjects, speciesTypes, extr, Cobra = False, specialEnding = ''):
    # Creates a list of species from given species, types and extracellular compartment
    if not speciesTypes:
        if Cobra:
            detectedDiff = False
        allSpecies = [[]]*len(speciesObjects)
        allMetabs  = [[]]*len(speciesObjects)
        for ind, species in enumerate(speciesObjects):
            curName = species.getName()
            curID   = species.getId()
            curComp = species.getCompartment()
            curExtr = bool(extr and curComp == extr)
            if specialEnding and curName.endswith(specialEnding):
                curName = curName[:-len(specialEnding)]
                if extr is not None:
                    curComp = extr
                else:
                    curComp = 'e'
                curExtr = True
            if Cobra:
                newExtr = bool(species.getBoundaryCondition() or curID.endswith('_b'))
                if newExtr != curExtr:
                    detectedDiff = True
                    curExtr = newExtr
            curNote = species.getNotesString()
            curFormula = species.getPlugin('fbc').getChemicalFormula()
            if not curFormula.strip():
                curFormula = getFormulaAndCharge(curNote)
            else:
                curFormula= convertFormula(curFormula)
            allSpecies[ind] = Species(curID, ind, {'name': curName, 'notes': curNote, 'formula': curFormula})
            allMetabs[ind]  = Metabolite(allSpecies[ind], curComp, ind, curExtr)
        if Cobra and detectedDiff:
            print("At least one metabolite constraint changes in COBRA!")
        return (allSpecies, allMetabs)
    else:
        print("Warning: I may not be able to handle this type of file yet!")
        return

def getFormulaAndCharge(noteString):
    # Extracts a formula, including a possible charge, from a note field for a species
    Formula = {}
    formula = re.search(formRE, noteString)
    if formula is not None:
        Formula = formula.group().replace('FORMULA: ', '').replace('</p>', '')
        Formula = convertFormula(Formula)
    charge = re.search(chargeRE, noteString)
    if charge is not None:
        Charge = charge.group().replace('CHARGE: ', '').replace('</p>', '')
        Charge = int(Charge)
        if 'e' in Formula:
            print('Error: there should never be any electrons in a formula!')
            print(noteString)
        else:
            Formula['e'] = -Charge
    return Formula

def createNetwork(reactionObjects, species, metabolites, Cobra = False):
    # Creates a network from a list of reactions, species and metabolites
    n = len(reactionObjects)
    allReactions = [[]]*n
    allMetabNames = [x.species.name for x in metabolites]
    if len(species) != len(metabolites):
        allComparts   = [x.compartment  for x in metabolites]
    allNames = [[]]*n
    biomassCandidates = []
    biomassVector = [zero] * n
    if Cobra:
        detectedDiff = False
    for ind, reaction in enumerate(reactionObjects):
        curName  = reaction.getName()
        allNames[ind] = curName
        if 'biomass' in curName.lower():
            biomassCandidates.append(curName)
        curID    = reaction.getId()
        curReags = reaction.getListOfReactants()
        curProds = reaction.getListOfProducts()
        curMods  = reaction.getListOfModifiers()
        curNote  = {'notes': reaction.getNotesString()}
        curRev   = reaction.getReversible()
        mult = 1
        if Cobra:
            KL = reaction.getKineticLaw()
            if KL:
                curLower, curUpper = None, None
                allParams = KL.getListOfParameters()
                for param in allParams:
                    if param.getId() == "LOWER_BOUND":
                        curLower = param.value
                    elif param.getId() == "UPPER_BOUND":
                        curUpper = param.value
                    elif param.getId() == 'OBJECTIVE_COEFFICIENT':
                        biomassVector[ind] = param.value
                if curLower is not None and curUpper is not None:
                    if curLower >= 0:
                        newRev = False
                    elif curUpper <= 0:
                        newRev = False
                        mult = -1
                    else:
                        newRev = True
                    if newRev != curRev:
                        detectedDiff = True
                        curRev = newRev
        allPairs = []
        for item in curReags:
            allPairs.append(processItem(item, allMetabNames, -mult))
        for item in curProds:
            allPairs.append(processItem(item, allMetabNames, mult))
        if curMods:
            geneCombo = CNF([[processItem(item, allMetabNames, 0)] for item in curMods])
        else:
            geneCombo = None
        allReactions[ind] = Reaction(curName, allPairs, ind, curRev, description = curNote, geneCombination = geneCombo)
    if Cobra and detectedDiff:
        print("At least one reversibility constraint changes in COBRA!")
    network = Network(species, metabolites, allReactions)
    network.createMatrices()
    if not Cobra:
        biomassIndex = findExactlyOne(biomassCandidates, 'biomass')
        if biomassIndex is not None:
            biomassName = biomassCandidates[biomassIndex]
            biomassVector[allNames.index(biomassName)] = one
    else:
        if sum(biomassVector) == 1:
            if len(biomassCandidates) == 0 or not allNames[biomassVector.index(one)] in biomassCandidates:
                print("The biomass reaction chosen by COBRA is different from what MONGOOSE would find!")
            elif len(biomassCandidates) >= 1:
                print(("The biomass reaction chosen by COBRA is among the " + str(len(biomassCandidates)) + " candidates found by MONGOOSE"))
        elif len(biomassCandidates) >= 1:
            print("COBRA does not detect any biomass reaction while MONGOOSE would")
        else:
            print("Neither MONGOOSE nor COBRA detect any biomass reaction")
    network.changeObjectiveFunction(biomassVector)
    return network

def processItem(item, metabNames, mult = 1):
    # Processes a single metabolite; the multiplier multiplies the coefficient
    curIndex = metabNames.index(item.getSpecies())
    if mult:
        curCoeff = convertToFraction(item.getStoichiometry())
        curCoeff = mult * curCoeff
        return [curIndex, curCoeff]
    else:
        return curIndex

def parseExcel(filename, biomassFilename = None, specialEnding = 'xt'):
    # Parses an Excel file containing the network description
    if biomassFilename is None:
        biomassFilename = filename
    workbook = xlrd.open_workbook(filename)
    sheetNames = workbook.sheet_names()
    reactSheet = findExactlyOne(sheetNames, 'reaction')
    if reactSheet is not None:
        reactSheet = workbook.sheet_by_index(reactSheet)
        network = parseReactionSheet(reactSheet, specialEnding)
        network.description = {'Source File': filename}
        metabSheet = findExactlyOne(sheetNames, ['metabolite', 'compound'])
        if metabSheet is not None:
            metabSheet = workbook.sheet_by_index(metabSheet)
            speciesInfo = parseMetaboliteSheet(metabSheet)
            mergeSpecies(network.species, speciesInfo)
        if not any(network.biomassCoefficients):
            biomassWorkbook = xlrd.open_workbook(biomassFilename)
            biomassSheetNames = biomassWorkbook.sheet_names()
            biomassSheet = findExactlyOne(biomassSheetNames, 'biomass')
            if biomassSheet is not None:
                biomassSheet = biomassWorkbook.sheet_by_index(biomassSheet)
                biomassReaction = parseBiomassSheet(biomassSheet, specialEnding)
                biomassPairs = convertPairs(biomassReaction, network.metabolites)
                network.addReaction('biomass', biomassPairs, False, description = {'Name': 'biomass'}, biomass = True)
        network.createMatrices()
        return network

def findSpecialSymbols(reactionList):
    # Finds special symbols used in representing the reactions in a given list of reactions
    keys = ['separator', 'coefficientOpen', 'coefficientClose', 'reversible', 'irreversible', 'compartmentOpen', 'compartmentClose', 'compartmentSeparate']
    symbols = dict.fromkeys(keys)
    allPieces = sum([x.split() for x in reactionList],[])
    if '+' in allPieces:
        symbols['separator'] = ' +'
    else:
        print('Error: no separator found!')
    candidates = [x for x in allPieces if not x.isalnum()]
    coeffs = [(x[0], x[1:-1], x[-1]) for x in candidates if len(x) > 2 and not x[0].isalnum() and not x[-1].isalnum() and re.match(coeffRE, x[1:-1])]
    coeffOpenCandidates = list(set([x[0] for x in coeffs]))
    symbols['coefficientOpen']  = processCandidates(coeffOpenCandidates, 'opening coefficients')
    coeffCloseCandidates = list(set([x[2] for x in coeffs]))
    symbols['coefficientClose'] = processCandidates(coeffCloseCandidates, 'closing coefficients')
    minCount = len(reactionList)/20 # assumption is that 5% or more of reactions are reversible/irreversible
    candSet = list(set(candidates).difference('+'))
    revs = [x for x in candSet if candidates.count(x) >= minCount]
    symbols['reversible'] = processCandidates(revs, 'reversible reactions')
    irrevCandidates = [x for x in revs if x != symbols['reversible']]
    symbols['irreversible'] = processCandidates(irrevCandidates, 'irreversible reactions')
    comps = [(x[0], x[1], x[2], x[3:]) for x in candidates if len(x) in [3,4] and not x[0].isalnum() and x[1].isalpha() and not x[2].isalnum()]
    compOpenCandidates = list(set([x[0] for x in comps]))
    symbols['compartmentOpen'] = processCandidates(compOpenCandidates, 'opening compartments')
    compCloseCandidates = list(set([x[2] for x in comps]))
    symbols['compartmentClose'] = processCandidates(compCloseCandidates, 'closing compartments')
    compSepCandidates = list(set([x[3] for x in comps]))
    if ':' in candidates:
        compSepCandidates.append(':')
    symbols['compartmentSeparate'] = processCandidates(compSepCandidates, 'separating compartments')
    return symbols

def extractMetabolites(reactionList, specialSymbols, special, defaultCompartment = 'c'):
    # Extracts the reactions, metabolite, species and compartment names and reversibility
    # information from a specified list of reactions and a dictionary of special symbols
    # The special string describes an ending corresponding to the external compartment
    # A default compartment for unspecified metabolites can also be optionally specified
    metabolites = set([])
    species = set([])
    compartments = set([])
    reversible = []
    allReacts = []
    allPairs = []
    compOpen = specialSymbols['compartmentOpen']
    compClose = specialSymbols['compartmentClose']
    compSeparate = specialSymbols['compartmentSeparate']
    revSymbol = specialSymbols['reversible']
    irrevSymbol = specialSymbols['irreversible']
    separator = specialSymbols['separator']
    coeffOpen = specialSymbols['coefficientOpen']
    coeffClose = specialSymbols['coefficientClose']
    for (number, reaction) in enumerate(reactionList):
        reaction = reaction.strip()
        if compOpen and reaction.startswith(compOpen):
            reaction = reaction[len(compOpen):]
            if compClose and compClose in reaction:
                index = reaction.find(compClose)
                compartment = reaction[:index]
                reaction = reaction[(index + len(compClose)):]
                reaction = reaction.strip()
                if compSeparate and reaction.startswith(compSeparate):
                    reaction = reaction[len(compSeparate):]
            else:
                print('Error: compartment opens, does not close!')
                print((str(number) + ') ' + reaction))
        else:
            compartment = defaultCompartment
        if revSymbol and revSymbol in reaction:
            reversible.append(number)
            index = reaction.find(revSymbol)
            LHS = reaction[:index]
            RHS = reaction[(index + len(revSymbol)):]
        elif irrevSymbol and irrevSymbol in reaction:
            index = reaction.find(irrevSymbol)
            LHS = reaction[:index]
            RHS = reaction[(index + len(irrevSymbol)):]
        else:
            print("Error: the reaction's reversibility is unclear!")
            print((str(number) + ') ' + reaction))
        leftMetabs = LHS.split(separator)
        leftMetabs = [processMetabolite(x, coeffOpen, coeffClose, compOpen, compClose, compartment, -1, special) for x in leftMetabs]
        rightMetabs = RHS.split(separator)
        rightMetabs = [processMetabolite(x, coeffOpen, coeffClose, compOpen, compClose, compartment, 1, special) for x in rightMetabs]
        allMetabs = [_f for _f in leftMetabs + rightMetabs if _f]
        metabolites.update([(x[1], x[2]) for x in allMetabs])
        species.update([x[1] for x in allMetabs])
        compartments.update([x[2] for x in allMetabs])
        allReacts.append(allMetabs)
    metabolites = sorted(list(metabolites))
    species = sorted(list(species))
    compartments = sorted(list(compartments))
    allPairs = [sorted([[metabolites.index((met[1],met[2])), met[0]] for met in react]) for react in allReacts]
    return(allPairs, metabolites, species, compartments, reversible)

def processMetabolite(metabolite, coeffOpen, coeffClose, compOpen, compClose, baseComp, mult, specialEnding, nullMetab = 'Nothing'):
    # Processes a metabolite using symbols for opening/closing coefficients/compartments, a base compartment and a
    # multiplier for the coefficient; additionally, a special ending treated as compartment indicator can be given.
    metabolite = metabolite.strip()
    coeff = one
    if coeffOpen and metabolite.startswith(coeffOpen):
        if coeffClose in metabolite:
            index = metabolite.find(coeffClose)
            try:
                coeff = convertToFraction(metabolite[len(coeffOpen):index])
                metabolite = metabolite[(index + len(coeffClose)):]
            except:
                print(('Error: cannot convert the coefficient in ' + metabolite + '; using 1 by default'))
        else:
            print('Error: coefficient opens, does not close!')
            print(metabolite)
        metabolite = metabolite.strip()
    elif ' ' in metabolite:
        index = metabolite.find(' ')
        try:
            coeff = convertToFraction(metabolite[:index])
            metabolite = metabolite[(index + len(' ')):]
        except:
            print(('Error: cannot convert the coefficient in ' + metabolite + '; using 1 by default'))
    if re.match(compRE, metabolite[-3:]):
        comp = metabolite[-2]
        metabolite = metabolite[:-3]
    elif metabolite.endswith(specialEnding): # special case!
        comp = specialEnding
        metabolite = metabolite[:-len(specialEnding)]
    else:
        comp = baseComp
    # special case: underscores are not meaningful symbols!
    metabolite = metabolite.strip()
    metabolite = metabolite.strip('_')
    if metabolite and metabolite != nullMetab:
        return (mult * coeff, metabolite, comp)
    else:
        return

def parseReactionSheet(sheet, specialEnding):
    # Parses a reaction sheet and creates a Network object out of it
    headerIndex = findHeaderIndex(sheet)
    startIndex = headerIndex + 1
    headerFields = [x.strip() for x in sheet.row_values(headerIndex)]
    reactColumn = findExactlyOne(headerFields, ['reaction', 'equation'])
    if reactColumn is not None:
        reactions = sheet.col_values(reactColumn, startIndex)
        usedRows = [x for x in range(len(reactions)) if reactions[x].strip()]
        usedRowsOffset = [x + startIndex for x in usedRows]
        allReactions = [[]] * len(usedRows)
        reactions = [reactions[x] for x in usedRows]
        specialSymbols = findSpecialSymbols(reactions)
        (allPairs, metabolites, species, compartments, reversible) = extractMetabolites(reactions, specialSymbols, specialEnding)
        geneColumn = findExactlyOne(headerFields, 'gene')
        if geneColumn is not None:
            genes = sheet.col_values(geneColumn, startIndex)
            genes = [genes[x].strip() for x in usedRows]
            formatSymbols = findFormatSymbols(genes)
            (allGenes, geneList) = parseGenes(genes, formatSymbols)
        enzymeColumn = findExactlyOne(headerFields, ['EC', 'E.C.'])
        if enzymeColumn is not None:
            enzymes = sheet.col_values(enzymeColumn, startIndex)
            enzymes = [enzymes[x].strip() for x in usedRows]
            (allEnzymes, enzymeList) = parseEnzymes(enzymes)
        usedColumns = [x for x in [reactColumn, geneColumn, enzymeColumn] if x is not None]
        extraColumns = [x for x in range(sheet.ncols) if x not in usedColumns]
        extraColumnNames = [headerFields[x] for x in extraColumns]
        for ind in range(len(usedRows)):
            currentPairs = allPairs[ind]
            currentRev = bool(ind in reversible)
            if geneColumn:
                currentGene = CNF(allGenes[ind])
            else:
                currentGene = None
            if enzymeColumn:
                currentEnzyme = allEnzymes[ind]
            else:
                currentEnzyme = None
            currentRow = sheet.row_values(usedRowsOffset[ind])
            extraFields = [currentRow[x] for x in extraColumns]
            description = dict(list(zip(extraColumnNames, extraFields)))
            allReactions[ind] = Reaction(currentEnzyme, currentPairs, ind, currentRev, currentGene, description = description)
        allSpecies, allMetabolites = createMetaboliteObjects(metabolites)
        network = Network(allSpecies, allMetabolites, allReactions)
        network.createMatrices()
        # try to find the biomass reaction
        foundBiomass = False
        if geneColumn is not None:
            biomassCandidates = [x for x in genes if 'biomass' in x]
            chosen = processCandidates(biomassCandidates, 'biomass reaction')
            if chosen is not None:
                biomassIndex = genes.index(chosen)
                network.biomassCoefficients[biomassIndex] = one
                foundBiomass = True
        if not foundBiomass:
            allSpecies = [x.name for x in network.species]
            biomassMetabs = [x for x in allSpecies if 'biomass' in x.lower()]
            chosen = processCandidates(biomassMetabs, 'biomass')
            if chosen is not None:
                biomassSpecies = network.species[allSpecies.index(chosen)]
                biomassMetabs = [x for x in network.metabolites if x.species == biomassSpecies and not x.external]
                biomassComparts = [x.compartment for x in biomassMetabs]
                chosenCompartment = processCandidates(biomassComparts, 'compartments for biomass')
                if chosenCompartment is not None:
                    biomassMetab = biomassMetabs[biomassComparts.index(chosenCompartment)].index
                    # find a reaction containing this metabolite as a product
                    fullMatrix = network.fullMatrix
                    reactCandidates = [x for x in range(len(fullMatrix[0])) if fullMatrix[biomassMetab][x] > 0]
                    if len(reactCandidates) == 1:
                        network.biomassCoefficients[reactCandidates[0]] = one
                        foundBiomass = True
                    else:
                        biomassCandidates = [network.printReactionFormula(x) for x in reactCandidates]
                        chosen = processCandidates(biomassCandidates, 'biomass reaction')
                        if chosen is not None:
                            biomassIndex = reactCandidates[biomassCandidates.index(chosen)]
                            network.biomassCoefficients[biomassIndex] = one
                            foundBiomass = True
        if not foundBiomass:
            print('Warning: biomass reaction not found; you may need to edit biomassCoefficients manually!')
        return network
    else:
        print('Error: cannot create the network!')
        return

def createMetaboliteObjects(metabolites, external = ['e','x']):
    # Creates species and metabolite objects from an ordered list of species-compartment pairs
    allSpecies = sorted(list(set([x[0] for x in metabolites])))
    speciesList = [Species(name, index, {}) for (index, name) in enumerate(allSpecies)]
    metaboliteList = [[]]*len(metabolites)
    for ind, metabolite in enumerate(metabolites):
        curSpecies = speciesList[allSpecies.index(metabolite[0])]
        curCompart = metabolite[1]
        metaboliteList[ind] = Metabolite(curSpecies, curCompart, ind, any([(x in curCompart.lower()) for x in external]))
    return speciesList, metaboliteList

def findFormatSymbols(geneList):
    # Finds the symbols used to format genes in a given list of gene combinations
    keys = ['andSymbol', 'orSymbol', 'groupOpen', 'groupClose']
    symbols = dict.fromkeys(keys)
    andSymbols = [' and ', ' AND ', '+', ';']
    orSymbols  = [' or ', ' OR ', '|', ',', '/', ' ']
    openSymbols = ['(', '[']
    closeSymbols = [')', ']']
    andSymbolCandidates = [x for x in andSymbols if any([x in y for y in geneList])]
    symbols['andSymbol'] = processCandidates(andSymbolCandidates, 'AND symbol')
    orSymbolCandidates  = [x for x in orSymbols if any([x in y for y in geneList])]
    symbols['orSymbol'] = processCandidates(orSymbolCandidates, 'OR symbol')
    if symbols['andSymbol'] is not None and symbols['orSymbol'] is not None:
        if symbols['andSymbol'].startswith(symbols['orSymbol']):
            # the or symbol is a space, the and symbol starts with a space => correct it!
            symbols['orSymbol'] = None
    openSymbolCandidates = [x for x in openSymbols if any([x in y for y in geneList])]
    symbols['groupOpen'] = processCandidates(openSymbolCandidates, 'open group symbol')
    closeSymbolCandidates = [x for x in closeSymbols if any([x in y for y in geneList])]
    symbols['groupClose'] = processCandidates(closeSymbolCandidates, 'close group symbol')
    return symbols
    
def parseGenes(geneList, formatSymbols):
    # Parses a given list of gene combinations given a dictionary of format symbols
    # NOTE: the gene format has to be either CNF or DNF; the returned value is in DNF
    andSymbol = formatSymbols['andSymbol']
    orSymbol = formatSymbols['orSymbol']
    openSymbol = formatSymbols['groupOpen']
    closeSymbol = formatSymbols['groupClose']
    allGenes = [[]]*len(geneList)
    geneSet = set([])
    for ind, gene in enumerate(geneList):
        geneGroups = [[]]
        DNF = False
        gene = gene.strip()
        if openSymbol and gene.startswith(openSymbol) and closeSymbol and gene.endswith(closeSymbol):
            gene = gene[len(openSymbol):(-len(closeSymbol))]
            enclosed = True
        else:
            enclosed = False
        if (orSymbol is None or not orSymbol in gene) and (andSymbol is None or not andSymbol in gene):
            geneGroups = [[gene]]
        elif (orSymbol is not None and orSymbol in gene) and (andSymbol is None or not andSymbol in gene):
            geneGroups = [[x.strip() for x in gene.split(orSymbol)]]
        elif (andSymbol is not None and andSymbol in gene) and (orSymbol is None or not orSymbol in gene):
            geneGroups = [[x.strip()] for x in gene.split(andSymbol)]
        elif (orSymbol is not None and orSymbol in gene) and (andSymbol is not None and andSymbol in gene):
            if enclosed:
                gene = openSymbol + gene + closeSymbol
            firstOpen, firstClose = -1, -1
            if openSymbol is not None:
                firstOpen = gene.find(openSymbol)
            if closeSymbol is not None:
                firstClose = gene.find(closeSymbol)
            if (firstOpen == -1 and firstClose != -1) or (firstOpen != -1 and firstClose == -1):
                print(('Error: unbalanced grouping symbols in ' + gene))
            else:
                if firstOpen == -1 and firstClose == -1:
                    print(('Warning: the grouping of ' + gene + ' is ambiguous! Assuming DNF form!'))
                    inner, outer = andSymbol, orSymbol
                    DNF = True
                else:
                    firstGroup = gene[firstOpen:firstClose]
                    if andSymbol in firstGroup and not orSymbol in firstGroup:
                        inner, outer = andSymbol, orSymbol
                        DNF = True
                    elif orSymbol in firstGroup and not andSymbol in firstGroup:
                        inner, outer = orSymbol, andSymbol
                    elif orSymbol in firstGroup and andSymbol in firstGroup:
                        if firstGroup.index(orSymbol) < firstGroup.index(andSymbol):
                            inner, outer = andSymbol, orSymbol
                            DNF = True
                            print(('Warning: the grouping of ' + gene + ' is possibly ambiguous! Assuming DNF form!'))
                        else:
                            inner, outer = orSymbol, andSymbol
                            print(('Warning: the grouping of ' + gene + ' is possibly ambiguous! Assuming CNF form!'))
                    else:
                        print(('Error: neither and nor or are present in the first group of ' + gene))
                groups = [x.strip() for x in gene.split(outer)]
                geneGroups = []
                for group in groups:
                    if openSymbol and group.startswith(openSymbol) and closeSymbol and group.endswith(closeSymbol):
                        group = group[len(openSymbol):-len(closeSymbol)]
                    geneGroups.append([x.strip() for x in group.split(inner)])
        geneSet.update(sum(geneGroups,[]))
        if DNF:
            geneGroups = DtoC(geneGroups)
        allGenes[ind] = geneGroups
    geneSet = sorted(list(geneSet))
    return(allGenes, geneSet)

def parseEnzymes(enzymes):
    # Parses a list of enzyme combinations, using a regular expression for EC numbers
    allEnzymes = [[]]*len(enzymes)
    enzymeSet = set([])
    for ind, enzyme in enumerate(enzymes):
        if enzyme:
            curEnzymes = re.findall(geneRE, enzyme)
            if not curEnzymes:
                print(('Warning: no enzymes found in ' + enzyme))
            allEnzymes[ind] = curEnzymes
            enzymeSet.update(curEnzymes)
    enzymeSet = sorted(list(enzymeSet))
    return(allEnzymes, enzymeSet)

def parseMetaboliteSheet(sheet):
    # This function parses a metabolite sheet into a dictionary of dictionaries (one per metabolite)
    headerIndex = findHeaderIndex(sheet)
    startIndex = headerIndex + 1
    headerFields = [x.strip() for x in sheet.row_values(headerIndex)]
    metabColumn = findExactlyOne(headerFields, 'abbrev')
    nameColumn  = findExactlyOne(headerFields, 'name')
    usedColumns = [x for x in [metabColumn, nameColumn] if x is not None]
    extraColumns = [x for x in range(sheet.ncols) if x not in usedColumns]
    extraColumnNames = [headerFields[x] for x in extraColumns]
    if (metabColumn is not None and nameColumn is None) or (metabColumn is None and nameColumn is not None):
        newColumn = findExactlyOne(extraColumnNames, 'metabolite')
        newColumn = headerFields.index(extraColumnNames[newColumn])
        extraColumns.remove(newColumn)
        extraColumnNames = [headerFields[x] for x in extraColumns]
        usedColumns.append(newColumn)
        if metabColumn is None:
            metabColumn = newColumn
        else:
            nameColumn = newColumn
    if 'compartment' in [x.lower() for x in extraColumnNames]:
        badIndex = [x.lower() for x in extraColumnNames].index('compartment')
        extraColumns.pop(badIndex)
        extraColumnNames.pop(badIndex)
    metabInfo = {}
    if metabColumn is not None:
        metabs = sheet.col_values(metabColumn, startIndex)
        usedRows = [x for x in range(len(metabs)) if metabs[x].strip()]
        usedRowsOffset = [x + startIndex for x in usedRows]
        if nameColumn is not None:
            names = sheet.col_values(nameColumn, startIndex)
            names = [names[x].strip() for x in usedRows]
        for ind in range(len(usedRows)):
            metab = metabs[ind].strip()
            if re.match(compRE, metab[-3:]): # remove compartment information
                metab = metab[:-3]
            currentRow = sheet.row_values(usedRowsOffset[ind])       
            extraFields = [currentRow[x] for x in extraColumns]
            curInfo = dict(list(zip(extraColumnNames, extraFields)))
            if nameColumn is not None:
                curInfo['name'] = names[ind]
            if metab in metabInfo:
                if metabInfo[metab] != curInfo:
                    for key in metabInfo[metab]:
                        if key in curInfo and curInfo[key] != metabInfo[metab][key]:
                            pass
            else:
                metabInfo[metab] = curInfo
    return metabInfo

def mergeSpecies(speciesList, infoDict):
    # Merges information on species into existing descriptions of species
    # Prints out a list of species with no extra information as well as
    # a list of species which have extra information but are not known.
    found = []
    for species in speciesList:
        if species.name in infoDict:
            if species.description is None:
                species.description = infoDict[species.name]
            else:
                species.description.update(infoDict[species.name])
            found.append(species.name)
        else:
            print(('No extra information found on ' + species.name))
    notFound = set(infoDict.keys()).difference(found)
    if len(notFound):
        print('Error: some species remain unused in the database!')
        print(('; '.join(sorted(notFound))))
    return

def parseBiomassSheet(sheet, specialEnding):
    # Processes a biomass sheet to extract a biomass reaction (as a list of pairs)
    headerIndex = findHeaderIndex(sheet)
    if headerIndex is not None:
        startIndex = headerIndex + 1
        headerFields = sheet.row_values(headerIndex)
        metabColumn = findExactlyOne(headerFields, 'component')
        coeffColumn = findExactlyOne(headerFields, 'coefficient')
        if metabColumn is not None and coeffColumn is not None:
            metabs = sheet.col_values(metabColumn, startIndex)
            coeffs = sheet.col_values(coeffColumn, startIndex)
            biomassReaction = [[metabs[ind], convertToFraction(str(coeffs[ind]))] for ind in range(len(metabs))]
            return biomassReaction
        else:
            print('Error: no component column or coefficient column found')
            return
    else: # assume there is only one row and one column, containing the reaction!
        reactions = sheet.row_values(0)
        specialSymbols = findSpecialSymbols(reactions)
        (allPairs, metabolites, species, compartments, reversible) = extractMetabolites(reactions, specialSymbols, specialEnding)
        biomassReaction = [[metabolites[x[0]][0], convertToFraction(str(x[1]))] for x in allPairs[0]]
        return biomassReaction

def convertPairs(reaction, metabolites, compartment = 'c'):
    # Converts the names in a list of [name, coeff] pairs into indices
    # the compartment specifies which compartment metabolites are in!
    speciesNames = ['']*len(metabolites)
    allNames = [metabolite.species.name for metabolite in metabolites]
    nameSet = set(allNames)
    allComps = [metabolite.compartment  for metabolite in metabolites]
    goodInds = [x for x in range(len(metabolites)) if allComps[x] == compartment]
    for ind in goodInds:
        curName = allNames[ind]
        speciesNames[ind] = curName
        nameSet.difference_update([curName])
    if len(nameSet):
        for name in nameSet:
            if allNames.count(name) == 1:
                speciesNames[allNames.index(name)] = name
                print(('Rescued the species ' + name + ' unambiguously'))
            else:
                print(('Warning: the species ' + name + ' is not in ' + compartment))
    newReaction = [[speciesNames.index(x[0]), x[1]] for x in reaction]
    return newReaction

def findHeaderIndex(sheet):
    # Finds the index of the header row in a given Excel sheet, or None if none is found
    ind = 0
    maxInd = sheet.nrows
    while(len([_f for _f in sheet.row_values(ind) if _f]) <= 1):
        ind += 1
        if ind == maxInd:
            print('Error: no header found!')
            return None
    return ind
