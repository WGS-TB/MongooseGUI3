# This file contains functions for defining the network formats
# Created by: Leonid Chindelevitch
# Last modified: January 30, 2017

from OutputProcessing import CreateSMatrix
from ModelProcessing import *
from Utilities import *
from Unrelated import *
from fractions import Fraction
from functools import reduce
from free_lunch_witnesses import *
from free_lunch_analysis import *
zero, one = Fraction(0), Fraction(1)
UNBLOCKED, TOPO, IRREV, STOICH, SEMI, SUBSET = list(range(6))
SURVIVING, DEADEND, REDUNDANT = 0, 1, 6

class Network:
    def __init__(self, species, metabolites, reactions, description = None):
        self.species = species
        self.metabolites = metabolites
        self.reactions = reactions
        self.compartments = sorted([x.compartment for x in self.metabolites])
        self.Matrix = None
        self.fullMatrix = None
        self.reducedMatrix = None
        self.reactionSubsets = None
        self.biomassCoefficients = [zero] * len(self.reactions)
        self.description = description
    def createMatrices(self):
        selfReacts = [reaction.pairs for reaction in self.reactions]
        numReacts = len(self.reactions)
        numMetabs = len(self.metabolites)
        self.fullMatrix = CreateSMatrix(selfReacts, numMetabs)
        intMetabs = [x for x in range(numMetabs) if not self.metabolites[x].external]
        self.Matrix = [self.fullMatrix[x] for x in intMetabs]
    def checkElementalBalance(self, excludeExchange = True):
        allDescriptions = [x.species.description for x in self.metabolites]
        allFormulas = [(y['formula'] if 'formula' in y.keys() else {}) for y in allDescriptions]
        if excludeExchange:
            ER = self.findExchangeReactions()
            checkAtomicBalance(allFormulas, [x for i,x in enumerate(self.reactions) if i not in ER])
        else:
            checkAtomicBalance(allFormulas, self.reactions)
    def updateReduction(self):
        if self.reducedMatrix is not None:
            print('Rerunning network reduction to update the reduced network.')
            self.reduceNetwork('NewReduction.txt')
    def reduceNetwork(self, filename = 'Reduction.txt'):
        N = self.Matrix
        Irr = self.findIrreversibleReactions()
        reduction = reduceMatrix(N, Irr, filename)
        self.applyReduction(reduction)
    def applyReduction(self, reducedMatrix):
        self.reducedMatrix = reducedMatrix[0]
        irreversibles = reducedMatrix[1]
        metabStatuses = reducedMatrix[2]
        reactStatuses = reducedMatrix[3]
        forwardOnly  = reducedMatrix[4]
        backwardOnly = reducedMatrix[5]
        subsets = reducedMatrix[6]
        internalMetabs = [x for x in self.metabolites if not x.external]
        for i, metab in enumerate(internalMetabs):
            metab.reductionStatus = metabStatuses[i]
        self.reactionSubsets = []
        anchor_to_subset = {}
        for i, subset in enumerate(subsets):
            anchor = subset[0][0]
            anchor_to_subset[anchor] = subset
        anchors = list(anchor_to_subset.keys())
        index = 0
        for i, reaction in enumerate(self.reactions):
            reaction.reductionStatus = reactStatuses[i]
            if reactStatuses[i] in [UNBLOCKED, SEMI]:
                geneCombo = self.reactions[reaction.index].geneCombination
                mult = one
                if reactStatuses[i] == SEMI:
                    if i in backwardOnly:
                        mult = -mult
                    elif i not in forwardOnly:
                        print('Error: this should never happen!')
                new_subset = Reaction_Subset([[reaction.index, mult]], index, not (index in irreversibles), geneCombo)
                index += 1
                self.reactionSubsets.append(new_subset)
            elif reactStatuses[i] == SUBSET:
                if i in anchors:
                    subset = anchor_to_subset[i]
                    geneCombo = reduce(lambda u, v: u.logicalAnd(v), [_f for _f in [self.reactions[x].geneCombination for x in [y[0] for y in subset]] if _f], CNF([]))
                    transformedSubset = [[y[0], y[1] * (1 - 2 * int(y[0] in backwardOnly))] for y in subset]
                    self.reactionSubsets.append(Reaction_Subset(transformedSubset, index, not (index in irreversibles), geneCombo))
                    index += 1
    def checkReduced(self):
        if self.reactionSubsets is None:
            print('Error: the network has not been reduced yet!')
            return False
        else:
            return True
    def findInternalMetabolites(self):
        Int = [i for i, m in enumerate(self.metabolites) if not m.external]
        return Int
    def convertToFullNumbering(self, metabs):
        Int = self.findInternalMetabolites()
        fullMetabs = [Int[x] for x in metabs]
        return fullMetabs
    def findIrreversibleReactions(self):
        Irr = [i for i, r in enumerate(self.reactions) if not r.reversible]
        return Irr
    def translateToReduced(self, reactionList):
        if self.checkReduced():
            badReactions = [r for r in reactionList if self.reactions[r].reductionStatus in [TOPO, STOICH, IRREV]]
            if len(badReactions):
                print('Warning: not all reactions on the list have survived the reduction!')
                reactionList = [r for r in reactionList if r not in badReactions] 
            L = len(reactionList)
            translatedReactions = [-1] * L
            count = 0
            for ind, subset in enumerate(self.reactionSubsets):
                if count == L:
                    break
                for pair in subset.pairs:
                    if pair[0] in reactionList:
                        translatedReactions[reactionList.index(pair[0])] = ind
                        count += 1
                        if count == L:
                            break
            return sorted(list(set(translatedReactions)))
    def findExchangeReactions(self):
        exchange = []
        externalInds = [i for i, m in enumerate(self.metabolites) if m.external]
        for ind, react in enumerate(self.reactions):
            for pair in react.pairs:
                if pair[0] in externalInds:
                    exchange.append(ind)
                    break
        return exchange
    def findExchangeReactionsReduced(self):
        if self.checkReduced():
            exchange = self.findExchangeReactions()
            exchangeR = self.translateToReduced(exchange)
            return exchangeR
    def findTopoBlockedMetabolites(self):
        return self.findMetabolites(Type = DEADEND)
    def findMetabolites(self, Type):
        if self.checkReduced():
            found = [i for i, m in enumerate(self.metabolites) if m.reductionStatus == Type]
            return found
    def findTopoBlockedReactions(self):
        return self.findReactions(Type = TOPO)
    def findIrrevBlockedReactions(self):
        return self.findReactions(Type = IRREV)
    def findStoichBlockedReactions(self):
        return self.findReactions(Type = STOICH)
    def findSemiBlockedReactions(self):
        return self.findReactions(Type = SEMI)
    def findReactions(self, Type):
        if self.checkReduced():
            found = [i for i, r in enumerate(self.reactions) if r.reductionStatus == Type]
            return found
    def findBiomassReaction(self):
        if sum(self.biomassCoefficients) == one and one in self.biomassCoefficients:
            return self.biomassCoefficients.index(one)
        else:
            return -1
    def findBiomassReactionReduced(self):
        biomass = self.findBiomassReaction()
        if biomass == -1:
            print('Cannot do anything as there is no biomass reaction!')
            return
        Target = self.translateToReduced([biomass])
        if not Target:
            print('Cannot do anything as the biomass reaction is blocked!')
            return
        return Target[0]
    def unblockBiomassReaction(self):
        if self.checkReduced():
            N = self.Matrix
            Irr = self.findIrreversibleReactions()
            growth = self.findBiomassReaction()
            status = self.reactions[growth].reductionStatus
            if status in [UNBLOCKED, SEMI, SUBSET]:
                print('The biomass reaction is not blocked; nothing to do here.')
            elif status in [TOPO, STOICH]:
                print('Returning a list of metabolites that can be unconstrained')
                return self.convertToFullNumbering(minimalUnblock(N, Irr, growth))
            else: # status = IRREV
                print('Returning a list of reactions that can be allowed to be reversible')
                return minimalThermoUnblock(N, Irr, growth)
    def changeObjectiveFunction(self, coefficients):
        assert len(coefficients) == len(self.reactions)
        assert min(coefficients) >= 0
        self.biomassCoefficients = coefficients
    def deleteReactions(self, reactionIndices):
        self.reactions = [r for i, r in enumerate(self.reactions) if i not in reactionIndices]
        self.createMatrices()
        self.updateReduction()
    def addReaction(self, reactionName, reactionPairs, reversible = False, description = {}, biomass = False):
        nextIndex = len(self.reactions)
        if any([x[0] < 0 for x in reactionPairs]):
            print("Error: cannot add reaction because some metabolites have negative indices")
            return
        if any([x[0] >= len(self.metabolites) for x in reactionPairs]):
            print("Error: cannot add reaction because a metabolite exceeds the maximum index, " + str(len(self.metabolites)))
            return
        if not all(x[1] for x in reactionPairs):
            print("Error: cannot add reaction because some metabolites have zero coefficients")
            return
        self.reactions.append(Reaction(reactionName, reactionPairs, nextIndex, reversible, description = description))
        if biomass:
            self.biomassCoefficients = [zero for x in range(len(self.reactions))]
            self.biomassCoefficients[-1] =  one
        print('Added the specified reaction with number ' + str(nextIndex))
        self.createMatrices()
        self.updateReduction()
    def deleteMetabolites(self, metaboliteIndices):
        assert min(metaboliteIndices) >= 0
        assert max(metaboliteIndices) < len(self.metabolites)
        self.metabolites = [m for i, m in enumerate(self.metabolites) if i not in metaboliteIndices]
        self.createMatrices()
        self.updateReduction()
    def addMetabolite(self, metaboliteName, compartment = None, external = False, description = {}):
        nextIndex = len(self.metabolites)
        curSpecies = [x.species.name for x in self.metabolites]
        if metaboliteName in curSpecies:
            metaboliteSpecies = self.metabolites[curSpecies.index(metaboliteName)].species
            print('Warning: a metabolite with the same species name already exists in the model.')
        else:
            metaboliteSpecies = Species(name = metaboliteName, index = len(curSpecies), description = {})
        self.metabolites.append(Metabolite(metaboliteSpecies, compartment, nextIndex, external, description = description))
        print('Added the specified metabolite with number ' + str(nextIndex))
        self.createMatrices()
        self.updateReduction()
    def findEssentialReactions(self):
        if self.checkReduced():
            Network = self.reducedMatrix
            growth = self.findBiomassReactionReduced()
            if growth is not None:
                Exchange = list(range(len(self.reactionSubsets)))
                allowed = Exchange
                Irr = [i for i, subset in enumerate(self.reactionSubsets) if not subset.reversible]
                Essential = findEssential(Network, growth, Exchange, allowed, rec = False, I = Irr)
                EssentialSubsets = [self.reactionSubsets[y] for y in Essential]
                EssentialFull = sum([[x[0] for x in y.pairs] for y in EssentialSubsets], [])
                numBasicEssential = len(self.reactionSubsets[growth].pairs)
                print(('There are ' + str(numBasicEssential - 1) + ' essential reactions in a subset with the growth reaction.'))
                print('Returning the remaining essential reactions')
                return sorted(EssentialFull)
    def findSyntheticLethalPairs(self,parallel=0):
        if self.checkReduced():
            Network = self.reducedMatrix
            Target = self.findBiomassReactionReduced()
            if Target is not None:
                Irr = [i for i, subset in enumerate(self.reactionSubsets) if not subset.reversible]
                Essential, Lethal = findEssentialLethal(Network, Target, rec = False, I = Irr, verbose = True,parallel=parallel)
                allPairs = []
                for item in Lethal:
                    first, second = item[0], item[1]
                    firstSubset = [x[0] for x in self.reactionSubsets[first].pairs]
                    secondSubset = [x[0] for x in self.reactionSubsets[second].pairs]
                    allPairs += [sorted([x,y]) for x in firstSubset for y in secondSubset]
                return sorted(allPairs)
    def findMinimalMedia(self):
        N = self.Matrix
        biomassIndex = self.findBiomassReaction()
        Irr = self.findIrreversibleReactions()
        Ext = set(range(len(self.metabolites)))
        Ext.difference_update(self.findInternalMetabolites())
        Ext = list(Ext)
        Exch = classifyExchange(self.fullMatrix, Ext, Irr, extra = True)
        Exch = sum([Exch[x] for x in range(6) if x != 2], []) # omit irreversible export reactions
        return findMinimalMedia(N, biomassIndex, Exch, rec = False, I = Irr)
    def printReactionFormula(self, reactionNumber):
        assert reactionNumber > 0
        assert reactionNumber < len(self.reactions)
        curReact = self.reactions[reactionNumber]
        curPairs = curReact.pairs
        curMetabs = self.metabolites
        curSpecies = [x.species.name for x in curMetabs]
        curComparts = [x.compartment for x in curMetabs]
        curRev = curReact.reversible
        reagents = [x for x in curPairs if x[1] < 0]
        products = [x for x in curPairs if x[1] > 0]
        LHS = ' + '.join([str(-x[1]) + ' ' + curSpecies[x[0]] + '[' + curComparts[x[0]] + ']' for x in reagents])
        RHS = ' + '.join([str( x[1]) + ' ' + curSpecies[x[0]] + '[' + curComparts[x[0]] + ']' for x in products])
        if curRev:
            sep = ' <-> '
        else:
            sep = ' -> '
        return (LHS + sep + RHS)
    def adjustReversibility(self):
        if all([reaction.reversible for reaction in self.reactions]):
            print('All the reactions are currently reversible; adjusting using upper and lower bounds.')
            descriptors = list(self.reactions[0].description.keys())
            UB = findExactlyOne(descriptors, 'upper')
            LB = findExactlyOne(descriptors, 'lower')
            if UB is not None and LB is not None:
                UB, LB = descriptors[UB], descriptors[LB]
                for ind, reaction in enumerate(self.reactions):
                    curDescription = reaction.description
                    reaction.reversible = bool(float(curDescription[UB]) > 0 and float(curDescription[LB]) < 0)
            else:
                print('Error: cannot find either upper bounds or lower bounds in the reaction descriptions.')
    def adjustReactionNames(self, oldName = 'EC numbers'):
        print('Adjusting reaction names')
        descriptors = list(self.reactions[0].description.keys())
        chosen = findExactlyOne(descriptors, 'name')
        if chosen is not None:
            goodDescriptor = descriptors[chosen]
            for ind, reaction in enumerate(self.reactions):
                if reaction.description is not None:
                    if reaction.name is not None:
                        reaction.description[oldName] = reaction.name
                    if goodDescriptor in reaction.description:
                        reaction.name = reaction.description[goodDescriptor]
                        del reaction.description[goodDescriptor]
    def adjustExternal(self, compart):
        print('Adjusting external metabolites')
        for metabolite in self.metabolites:
            if metabolite.compartment == compart:
                metabolite.external = True
        self.createMatrices()
    def adjustCompartments(self, compart = 'E', startPos = 1):
        if len(self.metabolites) == len(self.species):
            print('Adjusting compartments')
            speciesNames = [x.name for x in self.species]
            Map = {}
            for name in speciesNames:
                curPos = len(name) + startPos if startPos < 0 else startPos
                if name[curPos:(curPos + len(compart))] == compart:
                    reducedName = name[:curPos] + name[(curPos + len(compart)):]
                    if reducedName not in speciesNames:
                        print(('Appending a new species: ' + reducedName)) 
                        self.species.append(Species(name = reducedName, index = len(speciesNames), description = {}))
                        speciesNames.append(reducedName)
                    Map[speciesNames.index(name)] = speciesNames.index(reducedName)
            for (key, value) in Map.items():
                self.metabolites[key].species = self.species[value]
                self.metabolites[key].compartment = compart
                self.metabolites[key].external = True
            self.species = [x for i,x in enumerate(self.species) if i not in Map]
            self.createMatrices()
    def findFreeLunchWitnesses (self,y_lower_bound=None,max_num_atoms=1e3):
        return free_lunch_witness_finder(self,y_lower_bound,max_num_atoms)
    def findFreeLunchMetabolites(self):
        return find_free_lunch_metabolites(self)
    def findMinimalFreeLunch(self,free_lunch_meta_list=None,print_mess=True):
        return minimal_free_lunch_finder(self,free_lunch_meta_list,print_mess)
    def findIndividuallyProducibleFLM(self):
        return find_ind_flm(self)
    def ElementalBalanceCheckWithChemicalFormula(self):
        return totalEBCHeck(self)

            
class NetworkDescription():
    def __init__(self, genus, species, iden = None, authors = None, pmid = None, date = None):
        self.genus = genus
        self.species = species
        self.iden = iden
        self.authors = authors
        self.pmid = pmid
        self.date = date
    def search(self, s):
        #return a score corresponding to how close the match is
        matchVal = 0
        for t in s.split():
            t = t.lower()
            if t is self.genus.lower():
                matchVal += 15
            if t is self.species.lower():
                matchVal += 15
            if t == str(self.iden).lower():
                matchVal += 100
            if t is self.authors[0].lower():
                matchVal += 20
            if t is self.authors[-1].lower():
                matchVal += 10
            if t in [a.lower() for a in self.authors[1:-1]]:
                matchVal += 5
            if t[-2:] == self.date[-2:]:
                matchVal += 5
            if t == self.pmid:
                matchVal += 100
        return matchVal
    def write(self):
        print(self.genus,self.species)
        print(self.iden)
        print(self.authors)
        print(self.pmid)
        print(self.date)

class ExchangeStatus:
    IMPORT = 0
    EXPORT = 1
    BOTH = 2
    NONE = 3

class Reaction:
    def __init__(self, name, pairs, index = None, reversible = True, geneCombination = None, reductionStatus = None, exchangeStatus = None, description = None):
        self.name = name
        # pairs = [(metabolite_index,coefficient)]
        # omitted metabolites have coefficient 0
        self.pairs = pairs
        self.index = index
        self.reversible = reversible
        self.geneCombination = geneCombination
        self.reductionStatus = reductionStatus
        self.exchangeStatus = exchangeStatus
        self.description = description
    def show(self):
        print(self.name, self.pairs, self.index, self.reversible, self.geneCombination, self.reductionStatus, self.exchangeStatus, self.description)

class Metabolite:
    def __init__(self, species, compartment = None, index = None, external = False, reductionStatus = None, description = None):
        self.species = species
        self.compartment = compartment
        self.index = index
        self.external = external
        self.reductionStatus = reductionStatus
        self.description = description

class Species:
    def __init__(self, name = None, index = None, description = None):
        self.name = name
        self.index = index
        self.description = description

class Reaction_Subset:
    def __init__(self, pairs, index = None, reversible = True, geneCombination = None):
        self.pairs = pairs
        self.index = index
        self.reversible = reversible
        self.geneCombination = geneCombination

class CNF:
    def __init__(self, clauses):
        self.clauses = clauses
    def logicalAnd(self, other):
        new = CNF(clauses = self.clauses)
        new.clauses = new.clauses + other.clauses
        new.logicalSimplify()
        return new
    def logicalOr(self,other):
        new = CNF(clauses = self.clauses)
        full = DtoC(new.clauses) + DtoC(other.clauses)
        new.clauses = DtoC(full)
        new.logicalSimplify()
        return new
    def logicalSimplify(self):
        curClauses = self.clauses
        curClauses = [list(set(x)) for x in curClauses]
        curClauses = extractUnique(curClauses)
        self.clauses = curClauses
