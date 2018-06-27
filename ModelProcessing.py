# This file contains functions for processing metabolic models with Mongoose
# Created by: Leonid Chindelevitch
# Last modified: January 30, 2017

import os, copy, itertools, random, subprocess, shelve

from functools import reduce
from decimal import Decimal
from math import gcd
from fractions import Fraction
from Utilities import *
import multiprocessing
import qsoptex

zero, one = Fraction(0), Fraction(1)

#Docker ESOLVER_PATH:
ESOLVER_PATH = "/qsopt-ex/build/esolver/.libs/esolver"
#ESOLVER_PATH = "/Users/christopherle/Documents/Leonid/qsopt-ex/build/esolver/.libs/esolver"
#ESOLVER_PATH = "/Users/Admin/Downloads/DownloadedSoftware/qsopt-ex/build/esolver/.libs/esolver"


def reduceMatrix(N, Irr, Filename = 'Reduction.txt'):
    # This function computes the reduced form of a given stoichiometric matrix
    # assuming that the specified list of reactions is irreversible.
    # The reduction proceeds in several steps, each of which is verified using
    # exact (rational arithmetic-based) linear programming whenever possible.
    # The sequence of steps is:
    # 1) Find and delete topology-blocked reactions.
    # 2) Find and delete stoichiometry-blocked reactions.
    # 3) Find and delete irreversibility-blocked reactions.
    # 4) Find and process the semi-blocked reactions.
    # 5) Find and group reaction subsets.
    # 6) Find and delete redundant constraints.

    # rows[i], cols[i] contains the iteration at which the i-th row, column was deleted
    # print("<> Here!")
    m, n = getSize(N)
    rows, cols = [0]*m, [0]*n
    # currows[i], curcols[i] is the index of the current i-th row, column in original matrix
    currows, curcols = list(range(m)), list(range(n))
    # currirrev is True for the current irreversible reactions, False otherwise
    curirrev = [bool(x in Irr) for x in range(n)]
    Iter = 1
    print((str(Iter) + ') Finding and deleting topology-blocked reactions.'))
    (deadRows, deadCols, current) = processDeadEnds(N, False)
    updateRecord(rows, currows, deadRows, Iter)
    updateRecord(cols, curcols, deadCols, Iter)
    currows  = filterOut(currows,  deadRows)
    curcols  = filterOut(curcols,  deadCols)
    curirrev = filterOut(curirrev, deadCols)
    Iter = 2
    print((str(Iter) + ') Finding and deleting irreversibility-blocked reactions.'))
    (TBlocked, current) = processTBlocked(current, curirrev)
    thermoBlocked = mapList(TBlocked, curcols)
    updateRecord(cols, curcols, TBlocked, Iter)
    curcols  = filterOut(curcols,  TBlocked)
    curirrev = filterOut(curirrev, TBlocked)
    Iter = 3
    print((str(Iter) + ') Finding and deleting stoichiometry-blocked reactions'))
    (SBlocked, current, curB) = processSBlocked(current)
    stoichioBlocked = mapList(SBlocked, curcols)
    updateRecord(cols, curcols, SBlocked, Iter)
    curcols  = filterOut(curcols,  SBlocked)
    curirrev = filterOut(curirrev, SBlocked)
    Iter = 4
    print((str(Iter) + ') Finding and processing the semi-blocked reactions'))
    (onlyPos, onlyNeg, current) = processUnidirectional(current, curirrev, option = 'null')
    onlyPositive = mapList(onlyPos, curcols)
    onlyNegative = mapList(onlyNeg, curcols)
    updateRecord(cols, curcols, onlyPos + onlyNeg, Iter)
    for i in onlyPos + onlyNeg:
        curirrev[i] = True
    # adjust the nullspace by flipping the reactions that can carry only negative fluxes
    for i in onlyNeg:
        curB[i] = [-curB[i][k] for k in range(len(curB[i]))]
    Iter = 5
    print((str(Iter) + ') Finding and grouping reaction subsets.'))
    (Enzymes, lumpedReacts, subsetReacts, current) = processSubsets(current, curB)
    Enzymes = [list(zip(mapList(subset[0], curcols), subset[1])) for subset in Enzymes]
    updateRecord(cols, curcols, subsetReacts, Iter)
    curcols  = filterOut(curcols,  lumpedReacts)
    curirrev = filterOut(curirrev, lumpedReacts)
    Iter = 6
    print((str(Iter) + ') Finding and deleting redundant constraints'))
    (redRows, current) = processRedundant(current)
    updateRecord(rows, currows, redRows, Iter)
    currows = filterOut(currows, redRows)
    reductionRecord = (rows, cols, onlyPositive, onlyNegative, Enzymes)
    describeReduction(reductionRecord, Filename, EBA = False)
    Irrev = findTrueIndices(curirrev)
    extraFilename = Filename[:-4] + 'Full' + Filename[-4:]
    writeResults(current, Irrev, reductionRecord, extraFilename)
    return (current, Irrev) + reductionRecord


def energyBalanceReduce(Matrix, Irrev, External, Filename = "EnergyReduction.txt", signs = False):
    m, n = getSize(Matrix)
    curirrev   = [bool(x in Irrev) for x in range(n)]
    # create a list of external reactions
    Internal = [x for x in range(n) if x not in External]
    # remove the external reactions from Matrix
    current = [filterOut(Matrix[x], External) for x in range(len(Matrix))]
    # rows[i], cols[i] contains the iteration at which the i-th row, column was deleted
    m, n = getSize(current)
    rows, cols = [0]*m, [0]*n
    # currirrev is True for the current irreversible reactions, False otherwise
    # remove the external reactions from curirrev
    curirrev = filterOut(curirrev, External)
    # currows[i], curcols[i] is the index of the current i-th row, column in original matrix
    currows, curcols = list(range(m)), list(range(n))
    # Step 1) find and remove any zero loops that are present in the network
    Iter = 1
    (loopInds, current) = processZeroLoops(current)
    loops = mapList(loopInds, curcols)
    Loops = mapList(loops, Internal)
    updateRecord(cols, curcols, loopInds, Iter)
    curcols  = filterOut(curcols,  loopInds)
    curirrev = filterOut(curirrev, loopInds)
    # Step 2) find and remove any energy-blocked irreversible reactions
    Iter = 2
    (EBlocked, current) = processEBlocked(current, curirrev)
    energyBlocked = mapList(EBlocked, curcols)
    EnergyBlocked = mapList(energyBlocked, Internal)
    updateRecord(cols, curcols, EBlocked, Iter)
    curcols  = filterOut(curcols,  EBlocked)
    curirrev = filterOut(curirrev, EBlocked)
    # Step 3) find and process the unidirectional reactions
    Iter = 3
    (onlyPos, onlyNeg, current) = processUnidirectional(current, curirrev, option = 'row')
    onlyPositive = mapList(onlyPos, curcols)
    OnlyPositive = mapList(onlyPositive, Internal)
    onlyNegative = mapList(onlyNeg, curcols)
    OnlyNegative = mapList(onlyNegative, Internal)
    updateRecord(cols, curcols, onlyPos + onlyNeg, Iter)
    for i in onlyPos + onlyNeg:
        curirrev[i] = True
    # Step 4) find and group together all isozyme subsets
    Iter = 4
    (Isozymes, deletedReacts, current, curirrev) = processIsozymes(current, curirrev, True)
    Isozymes = [mapList(subset, curcols) for subset in Isozymes]
    updateRecord(cols,  curcols, deletedReacts, Iter)
    curcols  = filterOut(curcols,  deletedReacts)
    curirrev = filterOut(curirrev, deletedReacts)
    # Step 5) find all possible sign combinations
    allSigns = []
    if signs:
        curIrrev = findTrueIndices (curirrev)
        curRev   = findFalseIndices(curirrev)
        r = len(curRev)
        if (r > 0 and r <= 20): # NOTE: Do NOT use if r > 20, as it might take too long!
            allSigns = checkAllSigns(current, curRev)
        else:
            print('No or too many internal reversible reactions!')
    # Step 6) put back the original exchange reactions
    Iter = 5
    origMat = transpose(Matrix)
    for y in OnlyPositive:     # multiply OnlyPositive by -1 (see the paper for the explanation)
            origMat[y] = [-origMat[y][x] for x in range(m)]
    allInds = sorted([Internal[x] for x in curcols] + External)
    fullMat = [origMat[x] for x in allInds]
    current = transpose(fullMat)
    curirrev = [bool(x in Irrev + OnlyPositive + OnlyNegative) for x in allInds]
    cols = [cols[Internal.index(x)] if x in Internal else 0 for x in range(len(origMat))]
    # Step 7) remove any redundant constraints
    Iter = 6
    (redRows, current) = processRedundant(current)
    updateRecord(rows, currows, redRows, Iter)
    currows = filterOut(currows, redRows)
    reductionRecord = (rows, cols, OnlyPositive, OnlyNegative, Isozymes)
    describeReduction(reductionRecord, Filename, EBA = True)
    Irrev = findTrueIndices(curirrev)
    newExternal = [i for i,x in enumerate(allInds) if x in External]
    extraFilename = Filename[:-4] + 'Full' + Filename[-4:]
    writeResults(current, Irrev, tuple([newExternal] + [allSigns]) + reductionRecord, extraFilename)
    return (current, Irrev, newExternal, allSigns) + reductionRecord


def fullIterativeReduce(Matrix, Irrev, External, Filename = "IterativeReduction.txt"):
    # Iteratively reduces a system by applying flux-balance and energy-balance constraints
    Iter = 1
    (mF,nF) = getSize(Matrix)
    while(True):
        print(('Performing iteration ' + str(Iter)))
        curFilename = Filename[:-4] + 'Energy' + str(Iter) + Filename[-4:]
        result = energyBalanceReduce(Matrix, Irrev, External, curFilename, False)
        Matrix, Irrev, External, Record = result[0], result[1], result[2], result[4:]
        (mE,nE) = getSize(Matrix)
        print(("The current size of the system is " + str(mE) + " by " + str(nE)))
        if mE == mF and nE == nF and not (Record[2] + Record[3]):
            break
        curFilename = Filename[:-4] + 'Flux' + str(Iter) + Filename[-4:]
        result = reduceMatrix(Matrix, Irrev, curFilename)
        Matrix, Irrev, Record = result[0], result[1], result[2:]
        External = findExternal(External, Record)
        (mF,nF) = getSize(Matrix)
        print(("The current size of the system is " + str(mF) + " by " + str(nF)))
        if mF == mE and nF == nE and not (Record[2] + Record[3]):
            break
        Iter += 1
    if len(Matrix) > 0:
        allRev = [x for x in range(len(Matrix[0])) if x not in Irrev + External]
        if len(allRev) < 20:
            print("Extracting the signs")
            Iter += 1
            curFilename = Filename[:-4] + 'Energy' + str(Iter) + Filename[-4:]
            result = energyBalanceReduce(Matrix, Irrev, External, curFilename, True)
            Matrix, Irrev, External = result[0], result[1], result[2]
    return (Matrix, Irrev, External)


def findExternal(initExternal, reductionRecord):
    External = []
    (rows, cols, onlyPositive, onlyNegative, Enzymes) = reductionRecord
    anchor_to_subset = {}
    for i, subset in enumerate(Enzymes):
        anchor = subset[0][0]
        anchor_to_subset[anchor] = subset
    anchors = list(anchor_to_subset.keys())
    index = 0
    for i, status in enumerate(cols):
        if status in [0,4]:
            if i in initExternal:
                External.append(index)
            index += 1
        if status == 5 and i in anchors:
            subset = anchor_to_subset[i]
            subsetReacts = [x[0] for x in subset]
            if any([j in initExternal for j in subsetReacts]):
                External.append(index)
            index += 1
    return External


def writeResults(Matrix, Irrev, Record, Filename, sep = '\t'):
    f = open(Filename, 'w')
    for row in Matrix:
        f.write(sep.join([str(x) for x in row]) + '\n')
    f.write('\n')
    f.write(sep.join([str(x) for x in Irrev]) + '\n')
    f.write('\n')
    for record in Record:
        f.write(sep.join([str(x) for x in record]) + '\n')
        f.write('\n')
    f.close()
    return


def describeReduction(reductionRecord, Filename, EBA = False):
    # This function creates a description of the reduction (based on a given record) in the specified filename
    # (rows, cols, onlyPositive, onlyNegative, Enzymes, loops, Isozymes) = reductionRecord
    (rows, cols, onlyPositive, onlyNegative, Enzymes) = reductionRecord
    f = open(Filename,'w')
    if EBA:
        f.write('Deleted ' + str(cols.count(1)) + ' zero loops\n')
        f.write('Removed ' + str(cols.count(2)) + ' energy blocked reactions!\n')
    else:
        f.write('Deleted ' + str(rows.count(1)) + ' topologically blocked rows!\n')
        f.write('Deleted ' + str(cols.count(1)) + ' topologically blocked columns!\n')
        f.write('Removed ' + str(cols.count(2)) + ' thermodynamically blocked reactions!\n')
        f.write('Removed ' + str(cols.count(3)) + ' stoichiometrically blocked reactions!\n')
    f.write('Identified ' + str(len(onlyPositive)) + ' reactions that can only go forward!\n')
    f.write('Identified and flipped ' + str(len(onlyNegative)) + ' reactions that can only go backward!\n')
    f.write('Grouped ' + str(sum([len(x) for x in Enzymes])) + ' reactions into ' + str(len(Enzymes)) + ' subsets!\n')
    f.write('Removed ' + str(rows.count(6)) + ' redundant constraints!\n')
    f.close()
    return


def matrixToGraph(Matrix):
    # This function creates a graph structure for the bipartite graph represented by Matrix.
    # The graph structure is a dictionary of incidence lists; 'R' for rows, 'C' for columns.
    D = {}
    for ind0, row in enumerate(Matrix):
        for ind1, elt in enumerate(row):
            if elt:
                str0 = 'R' + str(ind0)
                str1 = 'C' + str(ind1)
                myAdd(D, str0, str1)
                myAdd(D, str1, str0)
    return D

def reduceByDegree(D, dmin, specialStart, Causes = None):
    # This function removes all nodes of a graph structure with degree not exceeding dmin.
    # An additional argument specialStart specifies the prefix of the nodes to be examined.
    # The optional Causes dictionary contains the causes of deletion for each deleted node.
    for x in list(D.keys()):
        if x.startswith(specialStart):
            curList = D[x]
            if len(curList) <= dmin:
                for y in curList:
                    for z in D[y]:
                        D[z].remove(y)
                        if Causes is not None and len(D[z]) <= dmin and z not in Causes:
                            Causes[z] = y
                    if Causes is not None and y not in Causes:
                        Causes[y] = x
                    del D[y]
                del D[x]
    return

def traceBack(item, record):
    # This function traces back the item through a dictionary of records
    trace = []
    while item in record:
        if item in trace:
            break
        else:
            trace.append(item)
            item = record[item]
    trace = list(reversed(trace))
    return trace

def fullReduce(D, dmin, specialStart, details):
    # This function iteratively prunes a graph structure until all degrees exceed dmin.
    # An additional argument specialStart specifies the prefix of the nodes to be examined.
    # If details = True, additionally returns the causes of deletion for each of the nodes.
    numNodes = len(D)
    prevNodes = numNodes + 1
    if details:
        Causes = {}
    else:
        Causes = None
    while (numNodes < prevNodes):
        reduceByDegree(D, dmin, specialStart, Causes)
        prevNodes = numNodes
        numNodes = len(D)
    return (D, Causes)

def processDeadEnds(Matrix, details):
    # This function finds and returns the dead ends (rows and columns) of a given
    # stoichiometric matrix, as well as a reduced version of the stoichiometric
    # matrix; if the details flag is True, also returns the cause of each "death"
    m, n = getSize(Matrix)
    D = matrixToGraph(Matrix)
    (reducedD, Causes) = fullReduce(D, 1, 'R', details)
    goodRows = sorted([int(x[1:]) for x in reducedD if x.startswith('R')])
    goodCols = sorted([int(x[1:]) for x in reducedD if x.startswith('C')])
    badRows = [x for x in range(m) if x not in goodRows]
    badCols = [x for x in range(n) if x not in goodCols]
    reducedMatrix = [[Matrix[x][y] for y in goodCols] for x in goodRows]
    if details:
        return (badRows, badCols, reducedMatrix, Causes)
    else:
        return (badRows, badCols, reducedMatrix)

def findIsozymes(N, Full = False):
    # This function determines all the sets of isozymes in a given stoichiometric matrix.
    # Isozymes are defined as reactions that are identical (not reverses of each other).
    # If Full = True, groups together all reactions that differ by a constant multiple.
    ISubsets = {}
    if Full:
        ID = findSubsets(transpose(N), True)
        for x in ID:
            ISubsets[x] = [y[0] for y in ID[x]]
    else:
        ID = groupIdentical(transpose(N))
        for x in ID:
            if len(x) > 1:
                ISubsets[x[0]] = x[1:]
    return ISubsets

def processIsozymes(N, irrev, Full = False):
    # This function processes all the sets of isozymes in a given stoichiometric matrix.
    # It returns an isozyme list and a set of duplicated columns, as well as the proper
    # reversibility information for each isozyme subset and a new stoichiometric matrix.
    ISubsets = findIsozymes(N, Full = Full)
    badCols = sum(list(ISubsets.values()), [])
    Isozymes = [[]]*len(ISubsets)
    for (ind, key) in enumerate(ISubsets.keys()):
        Isozymes[ind] = [key] + ISubsets[key]
        irrev[key] = all([irrev[x] for x in Isozymes[ind]])
    newIrrev = filterOut(irrev, badCols)
    (m, n) = getSize(N)
    newN = [filterOut(N[x], badCols) for x in range(m)]
    return (Isozymes, badCols, newN, newIrrev)

def processZeroLoops(N):
    # This function finds the indices of zero loops and removes them from a given matrix.
    (m, n) = getSize(N)
    loopInds = [y for y in range(n) if not any([N[x][y] for x in range(m)])]
    if len(loopInds) == n:
        loopInds = list(range(1,n))
    newN = [filterOut(N[x], loopInds) for x in range(m)]
    return (loopInds, newN)


def findRedundant(N):
    # This function returns a list of non-redundant rows in a given matrix.
    # These non-redundant rows form a basis for the rowspace of the matrix.
    Nt = transpose(N)
    (A, active) = GaussJordan(Nt, Gauss = True)
    return active


def processRedundant(N):
    # This function returns a list of redundant rows and a matrix without them.
    m, n = getSize(N)
    zeroRows = [x for x in range(m) if not any(N[x])]
    goodRows = filterOut(list(range(m)), zeroRows)
    if goodRows:
        newN = filterOut(N, zeroRows)
        active = findRedundant(newN)
        inactive = filterOut(list(range(len(newN))), active)
        inactive = mapList(inactive, goodRows)
        redRows = zeroRows + inactive
        newN = [newN[x] for x in active]
    else:
        redRows = list(range(1,m))
        if n > 0:
            newN = [N[0]]
        else:
            newN = N
    return (redRows, newN)

def findTBlocked(N, Irrev, basename = 'TBlocked.lp', restricted = True, option = 'row', negated = False, rev = False):
    # This function finds all the thermodynamically blocked reactions in a metabolic network
    # given by its stoichiometric matrix and a list of irreversible reactions. See paper for
    # a detailed justification of the algorithm. Note: returned reactions are irreversible!
    # For the meaning of the "restricted" parameter, see the function findPosSupport below.
    Iter = 0
    index = basename.find('.lp')
    weight = [-1 if negated else 1]*len(Irrev)
    found = set()
    Min = -1 if rev else 0
    while len(found) < len(Irrev):
        curName = basename[:index] + str(Iter) + basename[index:]
        # print("FINE")
        (val, vec) = findPosSupport(N, Irrev, weight, curName, Min = Min, restricted = restricted, option = option)
        # print("STILL FINE")
        if val > 0:
            Iter = Iter + 1
            if rev:
                newlyFound = set([int(y[1:]) for y,v in vec.items() if y.startswith('Y') and (-1 if negated else 1) * v > 0])
            else:
                newlyFound = set([int(y[1:]) for y in list(vec.keys()) if y.startswith('Y')])
            found.update(newlyFound.intersection(Irrev))
            # Recall that the solution contains only the variables whose values are > 0 (< 0 if negated and rev are True)
            for y in found:
                x = Irrev.index(y) # translate into the original indices!
                weight[x] = 0
        else:
            break
    found = sorted(list(set(found)))
    return found


def processTBlocked(N, irrev):
    # This function returns a list of thermodynamically blocked reactions and a matrix without them.
    Irrev = findTrueIndices(irrev)
    TBlocked = findTBlocked(N, Irrev)
    newN = [filterOut(N[x], TBlocked) for x in range(len(N))]
    return (TBlocked, newN)


def processEBlocked(N, irrev):
    # This function returns a list of energy blocked reactions and a matrix without them.
    Irrev = findTrueIndices(irrev)
    EUnblocked = findTBlocked(N, Irrev, basename = 'EBlocked.lp', restricted = False)
    EBlocked = [x for  x in Irrev if x not in EUnblocked]
    newN = [filterOut(N[x], EBlocked) for x in range(len(N))]
    return (EBlocked, newN)


def prepareForCplex(Matrix):
    Matrix, Mults = makeIntegral(transpose(Matrix), mults = True, decim = True) # integralize column-wise
    Matrix = transpose(Matrix)
    return [[Decimal(x.numerator)/Decimal(x.denominator) for x in y] for y in Matrix], Mults


# def findPosSupportNew(N, support, weight = [1], Filename = 'trial.lp', Min = 0, restricted = True, option = 'row'):
#     # This function finds the vector optimizing a given weight in the row/nullspace of N whose
#     # support is restricted to a given set of entries; those entries must be non-negative!
#     # Note: if the weight vector has a single component, it is automatically taken to be 1!
#     # If a nonzero Min value is specified, the components in the support are at least Min.
#     # If restricted = False, same but support is NOT restricted to the given set of entries.
#     print("NEW fps!")
#     m, n = getSize(N)
#     # intro = 'Problem\n'
#     p = qsoptex.ExactProblem()
#     # opt = 'Maximize\n'
#     p.set_objective_sense(qsoptex.ObjectiveSense.MAXIMIZE)
#
#
#
#     # if Min:
#     #     if Min > 0:
#     #         bounds += ''.join(['Y' + str(j) + ' >= ' + str(Min) + '\n' for j in support])
#     #     else:
#     #         bounds += ''.join([str(Min) + ' <= Y' + str(j) + ' <= ' + str(-Min) + '\n' for j in support])
#     # else:
#     #     bounds += ''.join(['Y' + str(j) + ' <= 1\n' for j in support])
#     # if Cplex:
#     #     const = const.replace('+ -', '-')
#
#     if Min:
#         if Min > 0:
#             curLower, curUpper = Min, None
#         else:
#             curLower, curUpper = Min, -Min
#     else:
#         curLower, curUpper = 0, 1
#
#     # Cplex is false
#     # if Cplex:
#     #     print()
#
#     # if len(weight) == len(support):
#     #     opt += ' + '.join([str(weight[i]) + ' Y' + str(support[i]) for i in range(len(support)) if weight[i]])
#     # elif len(weight) == 1:   # equal weights, take all of them equal to 1
#     #     opt += ' + '.join(['Y' + str(i) for i in support])
#     # else:
#     #     print('Error: the weight vector is not of the right length!')
#     # opt += '\n'
#     # const = 'Subject' + ' To' * int(Cplex) + '\n'
#     variables_to_pass = []
#
#
#     if len(weight) == len(support):
#         for ind, item in enumerate(support):
#             p.add_variable(name = 'Y' + str(item), objective = weight[ind], lower = curLower, upper = curUpper)
#             variables_to_pass.append('Y' + str(item))
#     elif len(weight) == 1:
#         for ind, item in enumerate(support):
#             p.add_variable(name = 'Y' + str(item), objective = 1, lower = curLower, upper = curUpper)
#             variables_to_pass.append('Y' + str(item))
#     else:
#         print('Error: the weight vector is not of the right length!')
#
#     # if restricted:
#     #     bounds += ''.join(['Y' + str(j) + ' = 0\n' for j in range(n) if j not in support])
#     for ind in range(n):
#         if ind not in support:
#             p.add_variable(name = 'Y' + str(ind), objective = 0, lower = 0, upper = (0 if restricted else None))
#             variables_to_pass.append('Y' + str(ind))
#
#     # if option == 'row':
#     #     const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(i) for i in range(m) if N[i][j]]) + ' + -1 Y' + str(j) + ' = 0\n' for j in range(n)])
#     # else:
#     #     const += ''.join([' + '.join([str(N[i][j]) + ' Y' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
#     # bounds = 'Bounds\n'
#     # if option == 'row':
#     #     bounds += ''.join(['X' + str(i) + ' free\n' for i in range(m) if [_f for _f in N[i] if _f]])
#
#     if option == 'row':
#         for ind in range(m):
#             if [_f for _f in N[ind] if _f]:
#                 p.add_variable(name = 'X' + str(ind), objective = 0, lower = None, upper = None)
#                 variables_to_pass.append('Y' + str(ind))
#         for j in range(n):
#             curDict = {'X'+str(i) : N[i][j] for i in range(m)}
#             curDict.update({'Y'+str(j) : -1})
#             p.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, curDict, rhs = 0)
#     else: # option == 'null'
#         for i in range(m):
#             print("NEW fps! (17)")
#             p.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, {'Y'+str(j) : N[i][j] for j in range(n)}, rhs = 0)
#
#     # opt += 'obj: '
#     # f = open(Filename, 'w')
#     # if not Cplex:
#     #     f.write(intro)
#     # f.write(opt)
#     # f.write(const)
#     # f.write(bounds)
#     # f.write('End\n')
#     # f.close()
#     # return processFile(Filename, True)
#
#     return processProblem(p, variables_to_pass)

def findPosSupport(N, support, weight = [1], Filename = 'trial.lp', Min = 0, restricted = True, Cplex = False, option = 'row'):
    # This function finds the vector optimizing a given weight in the row/nullspace of N whose
    # support is restricted to a given set of entries; those entries must be non-negative!
    # Note: if the weight vector has a single component, it is automatically taken to be 1!
    # If a nonzero Min value is specified, the components in the support are at least Min.
    # If restricted = False, same but support is NOT restricted to the given set of entries.
    # print("OLD fps!")
    m, n = getSize(N)
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Maximize\n'
    opt += 'obj: '
    if len(weight) == len(support):
        opt += ' + '.join([str(weight[i]) + ' Y' + str(support[i]) for i in range(len(support)) if weight[i]])
    elif len(weight) == 1:   # equal weights, take all of them equal to 1
        opt += ' + '.join(['Y' + str(i) for i in support])
    else:
        print('Error: the weight vector is not of the right length!')
    opt += '\n'
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    if option == 'row':
        const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(i) for i in range(m) if N[i][j]]) + ' + -1 Y' + str(j) + ' = 0\n' for j in range(n)])
    else:
        const += ''.join([' + '.join([str(N[i][j]) + ' Y' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
    bounds = 'Bounds\n'
    if option == 'row':
        bounds += ''.join(['X' + str(i) + ' free\n' for i in range(m) if [_f for _f in N[i] if _f]])
    if Min:
        if Min > 0:
            bounds += ''.join(['Y' + str(j) + ' >= ' + str(Min) + '\n' for j in support])
        else:
            bounds += ''.join([str(Min) + ' <= Y' + str(j) + ' <= ' + str(-Min) + '\n' for j in support])
    else:
        bounds += ''.join(['Y' + str(j) + ' <= 1\n' for j in support])
    if restricted:
        bounds += ''.join(['Y' + str(j) + ' = 0\n' for j in range(n) if j not in support])
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    # exit()
    # f = open(Filename, 'r')
    # print("content:\n" + str(f.readlines()))
    # f.close
    return processFile(Filename, True)

def findSBlocked(N, NB = False):
    # This function finds all the stoichiometrically blocked reactions in a metabolic network
    # given by its stoichiometric matrix. NB is True if the nullspace basis is to be returned.
    SBlocked = []
    (A,R) = GaussJordan(N)
    B = NullspaceBasis((A,R), True)
    SBlocked = [j for j in range(len(B)) if not any(B[j])]
    if NB:
        return (SBlocked, B)
    else:
        return SBlocked

def processSBlocked(N):
    # This function returns a list of stoichiometrically blocked reactions and a matrix without them.
    (SBlocked, B) = findSBlocked(N, True)
    newN = [filterOut(N[x], SBlocked) for x in range(len(N))]
    newB = filterOut(B, SBlocked)
    return (SBlocked, newN, newB)

def checkAllSigns(N, Rev):
    # This function produces a list of all possible sign combinations for the reversible reactions
    # in the row combinations of N, assuming that the irreversible reactions must be non-negative.
    r = len(Rev)
    allSigns = []
    counter = 0
    maxCounter = 2**r
    print(('There are ' + str(maxCounter) + ' combinations to process'))
    for signs in itertools.product('+-', repeat = r):
        counter += 1
        if (counter % 100 == 0):
            print(('Processed ' + str(counter) + ' combinations so far'))
        result = checkSigns(N, Rev, signs, Filename = 'signs' + str(counter) + '.lp')
        if result:
            allSigns.append(flipSigns(signs)) # need to flip them because entropy increases!
    return allSigns

def flipSigns(signs):
    # This function flips - to + and + to - in a given tuple (first converted to string)
    signs = ''.join(signs)
    signsF = signs.replace('-','p').replace('+','-').replace('p','+')
    return signsF

def checkSigns(N, Rev, signs, Filename = 'signs.lp', Cplex = False):
    # Thus function checks whether a particular sign pattern on the reversible reactions gives
    # a feasible vector in the rowspace of the input matrix; signs is a string of '+' and '-'.
    m, n = getSize(N)
    negIndices = [i for i, x in enumerate(signs) if x == '-']
    negatives = mapList(negIndices, Rev)
    posIndices = [i for i, x in enumerate(signs) if x == '+']
    positives = mapList(posIndices, Rev)
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Minimize\n'
    if Cplex:
        opt += ('Y' + str(0) + '\n')
    else:
        opt += '\n'
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(i) for i in range(m) if N[i][j]]) + ' + -1 Y' + str(j) + ' = 0\n' for j in range(n)])
    bounds = 'Bounds\n'
    bounds += ''.join(['X' + str(i) + ' free\n' for i in range(m) if [_f for _f in N[i] if _f]])
    bounds += ''.join([   '1 <= Y' + str(j) + '<= inf\n' for j in positives])
    bounds += ''.join(['-inf <= Y' + str(j) + '<= -1 \n' for j in negatives])
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    val = processFile(Filename, False)
    if (type(val) == type([]) and len(val) == 0): # infeasible
        return False
    else:
        return True

def vectorInSpan(N, vec, Filename = 'trial.lp', Cplex = False):
    # This function determines whether a given vector vec is in the row span of a matrix.
    # This is achieved by solving a linear program using QSOpt_ex and checking its value.
    m, n = getSize(N)
    n1 = len(vec)
    if n != n1:
        print('Problem: the vector is not compatible with the matrix!')
        return False
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Maximize\n'
    opt += 'obj: '
    opt += ('Y' + '\n')
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(i) for i in range(m) if N[i][j]]) + (' + ' + str(vec[j]) + ' Y')*int(bool(vec[j])) + '= 0\n' for j in range(n)])
    # Note that this slightly strange syntax ensures that no zero coefficients appear in the LP file
    bounds = 'Bounds\n'
    bounds += ''.join(['-1 <= X' + str(i) + ' <= 1\n' for i in range(m) if [_f for _f in N[i] if _f]])
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    val = processFile(Filename, False)
    if val:
        return True
    else:
        return False


def computeDistance(N, vec, norm = 'inf', Irrev = [], Filename = 'Distance.lp', Cplex = False):
    # This function computes the distance from a given vector vec to the row span of a matrix.
    # The possible options for norm are 'one', 'two' and 'inf'; 'two' is currently unavailable.
    # The optional list Irrev specifies the set of row coefficients required to be nonnegative.
    m, n = getSize(N)
    n1 = len(vec)
    if n != n1:
        print('Problem: the vector is not compatible with the matrix!')
        return
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Minimize\n'
    opt += 'obj: '
    opt += ('Y' + '\n')
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(i) for i in range(m) if N[i][j]]) + ' - V' + str(j) + ' = 0\n' for j in range(n)])
    const += ''.join(['V' + str(j) + ' - T' + str(j) + ' = ' + str(vec[j]) + '\n' for j in range(n)])
    const += ''.join(['Y' + str(j) + ' + T' + str(j) + '>= 0\n' for j in range(n)])
    const += ''.join(['Y' + str(j) + ' - T' + str(j) + '>= 0\n' for j in range(n)])
    if norm == 'one':
        const += ' + '.join(['Y' + str(j) for j in range(n)]) + '- Y = 0\n'
    elif norm == 'inf':
        const += ''.join(['Y' + str(j) + ' - Y <= 0\n' for j in range(n)])
    else:
        print(('Error: the option ' + norm + ' is not a valid option for the norm!'))
        return
    bounds = 'Bounds\n'
    bounds += ''.join(['X' + str(i) + ' free\n' for i in range(m) if i not in Irrev and [_f for _f in N[i] if _f]])
    bounds += ''.join(['V' + str(j) + ' free\n' for j in range(n)])
    bounds += ''.join(['T' + str(j) + ' free\n' for j in range(n)])
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    return processFile(Filename, True)

def findDistance(N, special, Irrev, norm = 'inf'):
    # This function determines the smallest change in biomass coefficients (reaction indexed
    # by special) needed in a stoichiometric matrix N to enable a unit of biomass production
    # while respecting irreversibility constraints specified by the vector of indices Irrev.
    Nt = transpose(N)
    negBiomass = [-x for x in Nt[special]]
    Nt = filterOut(Nt, [special])
    # recompute the indices of the irreversible reactions for Nt
    redIrrev = [x for x in Irrev if x < special] + [x - 1 for x in Irrev if x > special]
    val, vec = computeDistance(Nt, negBiomass, norm = norm, Irrev = redIrrev)
    redVec  = {}
    for x in list(vec.keys()):
        if x.startswith('V'):
            redVec[x] = vec[x]
    return val, redVec

####parallelization
def parallelizedNegFindFeasible(N, Irrev, pos  ,index, option,subList,out_q_neg):
    count = len(subList) * index
    outdict = {}
    for ind, react in enumerate(subList): # changed from subOnlyNeg to subList!
        (val1, vec1) = findFeasible(N, react, Irrev, pos, ('sub+' + str(index) + 'V' + str(ind) + 'sets.lp'), option=option)
        if (type(val1) == type([]) and len(val1) == 0):  # infeasible
            outdict[count+ind] = react #puts infeasible reactions into dictionary
    out_q_neg.put(outdict) #places dictionary into queue asynchronously

def parallelizedPosFindFeasible(N, Irrev, pos  ,index, option,subList,out_q_pos):
    outdict = {}
    count = len(subList) * index
    for ind, react in enumerate(subList):# changed from subOnlyPos to subList!
        (val0, vec0) = findFeasible(N, react, Irrev, pos, ('sub-' + str(index) + 'V' + str(ind) + 'sets.lp'), option=option)
        if (type(val0) == type([]) and len(vec0) == 0):  # infeasible
            outdict[count+ind] = react #puts infeasible reactions into dictionary
    out_q_pos.put(outdict) #places dictionary into queue asynchronously


##end functions for parallelization

def findUnidirectional(N, Irrev, option = 'null', verbose = False, parallel = 0):
    # This function finds all unidirectional (effectively irreversible) reactions.
    # NOTE: It assumes that all the reactions in the network can have nonzero flux;
    # otherwise may incorrectly classify blocked reactions as only negative.
    # The option can be 'null' for nullspace (default) or 'row' for rowspace.
    # If parallel > 0, this specifies the number of cores available for processing.
    m, n = getSize(N)
    onlyPos, onlyNeg = [], []
    allRev = [i for i in range(n) if not i in Irrev]
    canBePositive = findTBlocked(N, allRev, basename = 'canBePositive.lp', restricted = False, option = option, rev = True)
    canBeNegative = findTBlocked(N, allRev, basename = 'canBeNegative.lp', restricted = False, option = option, rev = True, negated = True)
    onlyPosCandidates = [i for i in allRev if i not in canBeNegative]
    onlyNegCandidates = [i for i in allRev if i not in canBePositive]
    parallel = 0 #testing

    if len(onlyNegCandidates) <= 1:
        onlyNeg = onlyNegCandidates
    else:
        if parallel>0: # this is to be distributed between the threads
            splitNegPairs = [onlyNegCandidates[i::parallel] for i in range(parallel)]
            # Each process will get  a queue to put its outdict into
            out_q_neg = multiprocessing.Queue()
            procs = []
            for index, subList in enumerate(splitNegPairs):
                p = multiprocessing.Process(
                        target=parallelizedNegFindFeasible,
                        args=(N, Irrev, True, index, option,subList,out_q_neg))
                procs.append(p)
                p.daemon = True
                p.start()

            # Collect all results into a single result dict. We know how many dicts
            # with results to expect.
            resultdict = {}
            for i in range(parallel):
                resultdict.update(out_q_neg.get())

            # Wait for all worker processes to finish
            for p in procs:
                p.join()
            onlyNeg = list(resultdict.values())

        else:
            for ind, react in enumerate(onlyNegCandidates):
                (val1, vec1) = findFeasible(N, react, Irrev, True,  'sub+' + str(ind) + 'sets.lp', option = option)
                if (type(val1) == type([]) and len(val1) == 0): # infeasible
                    onlyNeg.append(react)
    onlyPosCandidates = [x for x in onlyPosCandidates if x not in onlyNeg]
    if len(onlyPosCandidates) <= 1:
        onlyPos = onlyPosCandidates
    else:
        if parallel>0: # this is to be distributed between the threads
            splitPosPairs = [onlyPosCandidates[i::parallel] for i in range(parallel)]

            # Each process will get queue to put its out dict into
            out_q_pos = multiprocessing.Queue()
            procs = []
            for index, subList in enumerate(splitPosPairs):
                p = multiprocessing.Process(
                        target=parallelizedPosFindFeasible,
                        args=(N, Irrev, False, index,option, subList,out_q_pos))
                procs.append(p)
                p.daemon = True
                p.start()

            # Collect all results into a single result dict. We know how many dicts
            # with results to expect.
            resultdict = {}
            for i in range(parallel):
                resultdict.update(out_q_pos.get())

            # Wait for all worker processes to finish
            for p in procs:
                p.join()
            onlyPos = list(resultdict.values())

        else:
            for ind, react in enumerate(onlyPosCandidates):
                (val0, vec0) = findFeasible(N, react, Irrev, False,  'sub-' + str(ind) + 'sets.lp', option = option)
                if (type(val0) == type([]) and len(val0) == 0): # infeasible
                    onlyPos.append(react)
    if verbose:
        print('This required ' + str(len(onlyNegCandidates) + len(onlyPosCandidates)) + ' linear programs')

    return (onlyPos, onlyNeg)

def processUnidirectional(N, irrev, option = 'null'):
    # This function returns a list of unidirectional reactions and a matrix without them.
    # The option can be 'null' for nullspace (default) or 'row' for rowspace.
    m, n = getSize(N)
    Irrev = findTrueIndices(irrev)
    (onlyPos, onlyNeg) = findUnidirectional(N, Irrev, option = option)
    newN = [[N[x][y] for y in range(n)] for x in range(m)]
    for x in range(m):
        for y in onlyNeg:
            newN[x][y] = -newN[x][y]
    return (onlyPos, onlyNeg, newN)

def findSubsets(N, NB = False):
    # This function returns all enzyme subsets as well as the multipliers for each reaction.
    # If NB is True, the input is assumed to be the nullspace basis; if not, it is computed.
    ESubsets = {}
    if NB:
        Matrix = N
    else:
        Matrix = NullspaceBasis(N)
    n = len(Matrix)
    redMatrix = [[] for i in range(n)]
    Mults = [0 for i in range(n)]
    for i in range(n):
        nzRow = [_f for _f in Matrix[i] if _f]
        if not nzRow:
            print('Error: zero row found!')
            return ESubsets
        else:
            # save the multiplier, make first non-zero entry 1
            Mults[i] = convertToFraction(nzRow[0])
            redMatrix[i] = [x / nzRow[0] for x in Matrix[i]]
    groups = groupIdentical(redMatrix)
    for group in groups:
        if len(group) > 1:
            mults = [Mults[x] for x in group]
            ESubsets[group[0]] = [[group[x], mults[x] / mults[0]] for x in range(1, len(group))]
    return ESubsets

def processSubsets(N, curB):
    # This function returns a list of reactions in enzyme subsets and a matrix without them.
    # It requires a basis for the nullspace of the matrix N as well as the matrix N itself!
    m, n = getSize(N)
    ESubsets = findSubsets(curB, True)
    newN = [[N[x][y] for y in range(n)] for x in range(m)]
    lumpedReacts = []
    Enzymes = [[]]*len(ESubsets)
    anchors = list(ESubsets.keys())
    for (ind, key) in enumerate(anchors):
        value = ESubsets[key]
        (reacts, ratios) = list(zip(*value))
        Enzymes[ind] = (tuple([key]) + reacts, tuple([one]) + ratios)
        lumpedReacts += list(reacts)
        for num, react in enumerate(reacts):
            ratio = ratios[num]
            for x in range(m):
                newN[x][key] += (ratio * newN[x][react])
    newN = [filterOut(newN[x], lumpedReacts) for x in range(m)]
    subsetReacts = anchors + lumpedReacts
    return (Enzymes, lumpedReacts, subsetReacts, newN)

def reconfigureNetwork(N, irrev):
    # This function returns a reconfigured version of the given stoichiometric matrix
    # along with a list of reversible reaction indices that have been split into two.
    m, n = getSize(N)
    rev = findFalseIndices(irrev)
    newN = [[N[x][y] for y in range(n)] for x in range(m)]
    for x in range(m):
        newN[x] += [-newN[x][y] for y in rev]
    m, n0 = getSize(newN)
    newIrrev = [True] * n0
    # find columns duplicated by the reconfiguration and remove them!
    (Isozymes, badCols, newN, newIrrev) = processIsozymes(newN, newIrrev)
    if any([len(subset) != 2 for subset in Isozymes]):
        print('Error: all isozyme subsets should have size 2 during reconfiguration!')
    badInds = [subset[1] - n for subset in Isozymes]
    rev = filterOut(rev, badInds)
    return (newN, rev)


def findFeasible(N, special, Irrev = [], pos = True, Filename = 'trial.lp', disable = [], negative = [], option = 'null', Cplex = False):
    # This function finds a feasible vector in the row/nullspace of N whose set of irreversible
    # reactions is given. The entry corresponding to special is 1 if pos is True, -1 otherwise.
    # Additional features: it is now possible to specify a subset of reactions to be disabled
    # as well as a subset of reactions which are constrained to have only negative flux.
    # The option can be 'null' for nullspace (default) or 'row' for rowspace.
    m, n = getSize(N)
    Rev = [x for x in range(n) if x not in Irrev]
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Maximize\n'
    if Cplex:
        opt += ('V' + str(special) + '\n')
    else:
        opt += '\n'
    # note: we are only looking for a feasible vector, hence no objective function required!
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    if option == 'row':
        const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(i) for i in range(m) if N[i][j]]) + ' + -1 V' + str(j) + ' = 0\n' for j in range(n)])
    else: # assumes option = 'null'
        const += ''.join([' + '.join([str(N[i][j]) + ' V' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
    if pos:
        const += 'V' + str(special) + ' = 1\n'
    else:
        const += 'V' + str(special) + ' = -1\n'
    if disable:
        const += ''.join(['V' + str(i) + ' = 0\n' for i in disable])
    bounds = 'Bounds\n'
    if negative:
        bounds += ''.join(['-inf <= V' + str(i) + ' <= 0\n' for i in negative])
    if option == 'row':
        bounds += ''.join(['X' + str(i) + ' free\n' for i in range(m) if [_f for _f in N[i] if _f]])
    bounds += ''.join(['V' + str(i) + ' free\n' for i in Rev if i not in negative and [_f for _f in [N[k][i] for k in range(m)] if _f]])
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    if not Cplex:
        return processFile(Filename, True)


def findRatio(N, react1, react2, Irrev, Max = True, ratio = 0, Filename = 'trial.lp', Cplex = False):
    # This function finds the minimum or the maximum ratio of two given entries in the nullspace
    # of N whose set of irreversible reactions is given. Maximize if Max is True, minimize if not.
    # If a ratio is specified, checks whether the difference react1 - ratio * react2 can equal 1.
    m, n = getSize(N)
    Rev = [x for x in range(n) if x not in Irrev]
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    if Max:
        opt = 'Maximize\n'
    else:
        opt = 'Minimize\n'
    if ratio != 0:
        num, den = ratio.numerator, ratio.denominator
        if Cplex:
            if num > 0:
                opt += str(den) + ' V' + str(react1) + ' - ' + str(num)  + ' V' + str(react2) + '\n'
            else:
                opt += str(den) + ' V' + str(react1) + ' + ' + str(-num) + ' V' + str(react2) + '\n'
        else:
            opt += '\n'
    else:
        opt += 'obj: '
        opt += 'V' + str(react2) + '\n'
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    if ratio != 0:
        if num > 0:
            const += str(den) + ' V' + str(react1) + ' - ' + str(num)  + ' V' + str(react2) + ' = 1\n'
        else:
            const += str(den) + ' V' + str(react1) + ' + ' + str(-num) + ' V' + str(react2) + ' = 1\n'
    else:
        const += 'V' + str(react1) + ' = 1\n'
    const += ''.join([' + '.join([str(N[i][j]) + ' V' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
    bounds = 'Bounds\n'
    bounds += ''.join(['V' + str(i) + ' free\n' for i in Rev if [_f for _f in [N[k][i] for k in range(m)] if _f]])
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    if not Cplex:
        return processFile(Filename)

def processFile(Filename, opt = False, destroyIn = True, destroyOut = True, suppressOutput = True):
    # This function processes the linear programming problem described in the specified file.
    # It submits it to the exact rational solver, processes the output and returns the solution
    outFile = Filename.replace('.lp', '.sol')
    if suppressOutput:
        fnull = open(os.devnull, 'w')
        subprocess.call([ESOLVER_PATH, "-O", outFile, Filename], stdout = fnull, stderr = fnull)
        fnull.close()
    else:
        subprocess.call([ESOLVER_PATH, "-O", outFile, Filename])
    if destroyIn:
        subprocess.call(["rm", Filename])
    result = parseOutput(outFile, opt)
    if destroyOut:
        subprocess.call(["rm", outFile])
    # print(result)
    # print("the answer is ")
    # print(result[0])
    # exit()
    return result


def parseOutput(Filename, opt = False, verbose = False):
    # This function parses a solution file in the QSOpt_ex format
    # and returns a list containing the value of the objective function ([] if
    # the problem is infeasible, [float('Inf')] if it is unbounded, else [num, den]).
    # If opt = TRUE, also returns a dictionary with values of the nonzero variables.
    f = open(Filename, 'r')
    g = f.readlines()
    g = [x.strip() for x in g]
    dico = {}
    if g[0] == 'status = OPTIMAL':
        valueLine = g[2]
        index = valueLine.index('=')
        value = convertToFraction(valueLine[index+2:])
    elif g[0] == 'status = UNBOUNDED':
        if verbose:
            print('Note: problem is unbounded!')
        value = [float('Inf')]
    elif g[0] == 'status = INFEASIBLE':
        if verbose:
            print('Note: problem is infeasible!')
        value = []
    else:
        if verbose:
            print('Problem: A solution was not found!')
        return
    if opt:
        if type(value) == type(zero):
            # the optimal value is finite
            pos = 4
            while "=" in g[pos]:
                index = g[pos].index('=')
                varName = g[pos][:index-1]
                varValue = g[pos][index+2:]
                dico[varName] = convertToFraction(varValue)
                pos += 1
        return (value, dico)
    else:
        return value

def mapBackCutset(cutsets, cols, prevRev, Isozymes):
    # Reverses the steps performed in reducing a stoichiometric matrix
    # in order to create cutsets corresponding to the given ones.
    m = len(cutsets)
    L = len(prevRev)
    survivors = [x for x in range(len(cols)) if cols[x] == 0]
    nPrev = len(survivors)
    n = nPrev + L
    reCutsets = [0]*m
    # Step 1: undo the reconfiguration
    for i in range(m):
        reCutsets[i] = [survivors[x] for x in cutsets[i] if x < nPrev]
        reCutsets[i] += [survivors[prevRev[x - nPrev]] for x in cutsets[i] if x >= nPrev]
    # Step 2: reconstruct the isozyme subsets
    K = len(Isozymes)
    for i in range(m):
        for k in range(K):
            if ([x for x in reCutsets[i] if x in Isozymes[k]]):
                reCutsets[i] += Isozymes[k]
        # As a last step, eliminate duplicates
        reCutsets[i] = sorted(list(set(reCutsets[i])))
    return reCutsets

def mapBackFlux(fluxes, cols, prevRev, Subsets):
    # Reverses the steps performed in reducing a stoichiometric matrix
    # in order to create flux vectors corresponding to the given ones.
    m = len(fluxes)
    n = len(fluxes[0])
    L = len(prevRev)
    nPrev = n - L
    n0 = len(cols)
    survivors = [x for x in range(n0) if cols[x] == 0]
    if len(survivors) != nPrev:
        print('Error: incompatible dimensions!')
    # Step 1: undo the reconfiguration
    reFluxes = [0]*m
    for i in range(m):
        reFluxes[i] = [0]*n0
        for j in range(nPrev):
            reFluxes[i][survivors[j]] = fluxes[i][j]
        for k in range(L):
            reFluxes[i][survivors[prevRev[k]]] -= fluxes[i][nPrev + k]
    # Step 2: reconstruct the enzyme subsets
    K = len(Subsets)
    for j in range(K):
        head = Subsets[j][0][0]
        rest = Subsets[j][1:]
        for i in range(m):
            base = reFluxes[i][head]
            for k in range(len(rest)):
                reFluxes[i][rest[k][0]] = base * rest[k][1]
    return reFluxes

def convertVector(Dict, length):
    # Converts a dictionary representation of a vector into a list representation
    vector = [zero]*length
    for (key, value) in Dict.items():
        vector[int(key[1:])] = value
    return vector

def convertBackVector(List, variable = 'V'):
    # Converts a list representation of a vector into a dictionary representation
    dict = {}
    for x in range(len(List)):
        if List[x]:
             dict[variable + str(x)] = List[x]
    return dict

def findEFMs(N, prevRev, rec = True, I = []):
    # This function finds a large set of EFMs in a stoichiometric network N
    # excluding the "trivial" ones of the type e_i + e_j where i, j come from
    # the reconfiguration of a single initially reversible reaction.
    # Note: assumes N has been reconfigured. Otherwise, set rec to be False
    # and specify the set of irreversible reactions in I (ignoring prevRev!)
    m, n = getSize(N)
    allEFMs = [0]*n
    if rec:
        L = len(prevRev)
        # the last L reactions come from reconfiguration
        for i in range(n):
            if i in prevRev:
                j = prevRev.index(i)
                z = (n - L) + j
                allEFMs[i] = findEFM(N, i, [z], 'EFM' + str(i) + '.lp')
            elif i >= (n - L):
                j = i - (n - L)
                z = prevRev[j]
                allEFMs[i] = findEFM(N, i, [z], 'EFM' + str(i) + '.lp')
            else:
                allEFMs[i] = findEFM(N, i, [], 'EFM' + str(i) + '.lp')
    else:
        for i in range(n):
            allEFMs[i] = findEFM(N, i, [], 'EFM' + str(i) + '.lp', rec = False, I = I)
    return allEFMs

def findEFM(N, special, zeros = [], Filename = 'trial.lp', rec = True, I = []):
    # This function finds a set of EFMs in a stoichiometric network N that
    # contain a special reaction with coefficient 1 and do not contain any
    # of the reactions listed in zeros. Note: assume N has been reconfigured.
    # Otherwise, set rec to be False and specify irreversible reactions in I.
    myEFMs = []
    if len(N) > 0:
        n = len(N[0])
        weight = [1]*n
        # contains the actual vectors
        EFMpos = []
        # contains the boolean lists of positions
        new = True
        Iter = 0
        while (new):
            (val, vec) = findMin1Norm(N, special, weight, zeros, EFMpos, 1e-5, Filename[:-3] + 'I' + str(Iter) + '.lp', 'null', rec, I)
            if not rec: # make sure only the relevant variables are kept!
                for x in list(vec.keys()):
                    if x[0] == 'T':
                        del vec[x]
            peeledVec = peelOff(vec, myEFMs)
            if peeledVec:
                positions = [False]*n
                weight = [1]*n
                for key in peeledVec:
                    pos = int(key[1:])
                    positions[pos] = True
                    weight[pos] = 0
                if positions in EFMpos or not positions[special]: # NOTE THE MORE STRINGENT REQUIREMENT HERE!
                    new = False
                else:
                    EFMpos.append(positions)
                    myEFMs.append(peeledVec)
            else:
                new = False
            Iter += 1
    return myEFMs

def peelOff(Vector, otherVectors):
    # This function removes multiples of given vectors from a given one,
    # so as to make the result as sparse as possible. The vectors are
    # assumed to be nonnegative and to have rational entries [num, den].
    # NOTE: In general, the order in which this is done matters; however,
    # if the vector is an FM and the other vectors are EFMs, it does not!
    L = len(otherVectors)
    givenVector = {}
    for x in Vector:
        givenVector[x] = Vector[x]
    for i in list(reversed(list(range(L)))):
        # try to subtract a multiple of the i-th vector from given one
        otherVector = {}
        for x in otherVectors[i]:
            otherVector[x] = otherVectors[i][x]
        if not ([x for x in otherVector if x not in givenVector]):
            # the support of otherVector is a subset of that of Vector
            ratios = [givenVector[x] / otherVector[x] for x in otherVector]
            myMin = min(ratios)
            for x in otherVector:
                givenVector[x] -= (myMin * otherVector[x])
                if givenVector[x] == zero:
                    del givenVector[x]
    return givenVector

def findUniqueEFMs(EFMs):
    # This function finds a subset of unique EFMs from a given set of EFMs
    # The defining criteria for equality is the same set of nonzero entries
    # Note that this function applies without modification to MCSs as well!
    L = len(EFMs)
    if L > 0:
        aux = [(sorted([x for x in EFMs[i]]),i) for i in range(L)]
        aux.sort(key = lambda x:x[0])
        unique = [aux[0][1]] + [aux[i][1] for i in range(1, L) if aux[i-1][0] != aux[i][0]]
        return [EFMs[x] for x in unique]
    else:
        return EFMs

def findMin1Norm(N, special, weight = [1], zeros = [], exclude = [], eps = 1e-5, Filename = 'trial.lp', option = 'null', rec = True, I = [], Cplex = False):
    # This function finds the vector of smallest overall weight in the nullspace of N whose
    # reactions are assumed to be all irreversible unless rec is specified to be False.
    # In that case, the reactions considered to be irreversible should be specified in I.
    # Note that rec is ignored if option = 'col'! The entry corresponding to special is 1.
    # If a zeros constraint is present, the corresponding entries are forced to equal 0.
    # If exclusion constraints are present, ensures the support differs from given supports.
    # Note: if the weight vector has a single component, it is automatically taken to be 1!
    # Available option values currently are 'null' for nullspace and 'col' for columnspace.
    m, n = getSize(N)
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Minimize\n'
    opt += 'obj: '
    if rec:
        var = 'V'
    else:
        var = 'T'
    if option == 'null':
        if len(weight) == n:
            opt += ' + '.join([str(weight[i]) + ' ' + var + str(i) for i in range(n) if weight[i]])
        else:   # equal weights, take all of them equal to 1
            opt += ' + '.join([var + str(i) for i in range(n)])
    elif option == 'col':
        if len(weight) == m:
            opt += ' + '.join([str(weight[i]) + ' T' + str(i) for i in range(m) if weight[i]])
        else:   # equal weights, take all of them equal to 1
            opt += ' + '.join(['T' + str(i) for i in range(m)])
    else:
        print('Error: unrecognized option!')
        return
    opt += '\n'
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    if option == 'null':
        const += ''.join([' + '.join([str(N[i][j]) + ' V' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
        if not rec:
            const += ''.join(['T' + str(j) + ' + V' + str(j) + '>= 0\n' for j in range(n)])
            const += ''.join(['T' + str(j) + ' - V' + str(j) + '>= 0\n' for j in range(n)])
    elif option == 'col':
        const += ''.join([' + '.join([str(N[i][j]) + ' V' + str(j) for j in range(n) if N[i][j]]) + ' - W' + str(i) + ' = 0\n' for i in range(m)])
        const += ''.join(['T' + str(j) + ' + W' + str(j) + '>= 0\n' for j in range(m)])
        const += ''.join(['T' + str(j) + ' - W' + str(j) + '>= 0\n' for j in range(m)])
    const += 'V' + str(special) + ' = 1\n'
    if zeros:
        const += ''.join(['V' + str(x) + ' = 0\n' for x in zeros])
    if exclude:
        k = len(exclude)
        s = len(exclude[0])
        if s == n:
            const += ''.join([' + '.join([var + str(j) for j in range(n) if not exclude[i][j]]) + ' >' + str(eps) + '\n' for i in range(k)])
        else:
            print("Error: the length of the excluded vectors is incorrect!")
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    if option == 'col':
        bounds = 'Bounds\n'
        bounds += ''.join(['W' + str(j) + ' free\n' for j in range(m) if [_f for _f in N[j] if _f]])
        f.write(bounds)
    elif not rec:
        bounds = 'Bounds\n'
        bounds += ''.join(['V' + str(j) + ' free\n' for j in range(n) if j not in I and [_f for _f in [N[i][j] for i in range(m)] if _f]])
        f.write(bounds)
    f.write('End\n')
    f.close()
    return processFile(Filename, True)

def LCM(a,b):
    # This function computes the least common multiple of two integers
    return (a//gcd(a,b))*b

def integralize(Vector, mult = False, decim = False):
    # This function takes a rational vector and returns the smallest integral multiple
    # Note: this function works for a dictionary as well as a list representation(!)
    # If you need the mutliplier to be returned as well, change mult to a True value.
    # If you only need the scaled vector to be an exact decimal, set decim to True.
    if type(Vector) == type({}):
        keys = list(Vector.keys())
        values = list(Vector.values())
        intVector = {}
    elif type(Vector) == type([]):
        keys = list(range(len(Vector)))
        values = Vector
        intVector = [[] for i in keys]
    lcm = reduce(LCM, [convertToFraction(z).denominator for z in values])
    if decim:
        while lcm % 2 == 0:
            lcm /= 2
        while lcm % 5 == 0:
            lcm /= 5
    for x in keys:
        curValue = lcm * convertToFraction(Vector[x])
        if decim:
            intVector[x] = curValue
        else:
            intVector[x] = curValue.numerator
    if not mult:
        return intVector
    else:
        return (intVector, lcm)

def makeIntegral(Matrix, mults = False, decim = False):
    # This function makes a matrix into an integral one by multiplying rows by integers.
    # If you need the mutlipliers to be returned as well, change mults to a True value.
    intMatrix = [integralize(Matrix[i], mults, decim) for i in range(len(Matrix))]
    if not mults:
        return intMatrix
    else:
        return ([x[0] for x in intMatrix], [x[1] for x in intMatrix])

def findMCSs(N):
    # This function finds a large set of MCSs in a stoichiometric network N,
    # several for each reaction. Note: assume N has been reconfigured.
    m, n = getSize(N)
    allMCSs = [0]*n
    for i in range(n):
            allMCSs[i] = findMCS(N, i, 'MCS' + str(i) + '.lp')
    return allMCSs

def findMCS(N, special, Filename = 'trial.lp'):
    # This function finds a set of MCSs in a stoichiometric network N that
    # block a specified target reaction. Note: assume N has been reconfigured.
    MCSpos = []
    if len(N) > 0:
        n = len(N[0])
        weight = [1]*n
        # contains the boolean lists of positions
        new = True
        Iter = 0
        while (new):
            (val, vec) = findMinAdded(N, special, weight, MCSpos, Filename = Filename[:-3] + 'I' + str(Iter) + '.lp')
            cutset = [int(x[1:]) for x in list(vec.keys()) if x.startswith('T')]
            curMCS = extractMinimal(cutset, testCutSet, [N, special, Filename[:-3] + 'I' + str(Iter) + '.lp'])
            if curMCS:
                positions = [False]*n
                weight = [1]*n
                for pos in curMCS:
                    positions[pos] = True
                    weight[pos] = 0
                if positions in MCSpos:
                    new = False
                else:
                    MCSpos.append(positions)
            else:
                new = False
            Iter += 1
    MCSs = [[y for y in range(len(x)) if x[y]] for x in MCSpos]
    return MCSs

def testCutSet(Cutset, N, Target, Filename = 'trial.lp', rec = True, I = [], Cplex = False):
    # This function determines whether a given subset of reactions represents a cutset
    # for a given target reaction in a network which is assumed to be irreversible,
    # unless rec is specified to be False. In that case, the reactions considered to be
    # irreversible should be specified in I.
    m, n = getSize(N)
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Minimize\n'
    if Cplex:
        opt += 'V' + str(Target) + '\n'
    else:
        opt += '\n'
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    const += 'V' + str(Target) + ' = 1\n'
    const += ''.join([' + '.join([str(N[i][j]) + ' V' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
    for react in Cutset:
        const += 'V' + str(react) + ' = 0\n'
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    if not rec:
        bounds = 'Bounds\n'
        bounds += ''.join(['V' + str(j) + ' free\n' for j in range(n) if j not in I and [_f for _f in [N[i][j] for i in range(m)] if _f]])
        f.write(bounds)
    f.write('End\n')
    f.close()
    val = processFile(Filename)
    return (type(val) == type([]) and len(val) == 0) # TRUE IFF THE PROBLEM IS INFEASIBLE

def testSet(Set, Function, Args):
    # General function for testing a monotone property and identifying a local minimum for it.
    # This function determines whether a given set represents a MINIMAL set that returns TRUE
    # when the Function is called with this Set as the first argument and the other given Args.
    # If the Function fails, returns False; otherwise, returns a subset that is minimal for it.
    if not Function(Set, *Args):
        return False
    else:
        return extractMinimal(Set, Function, Args)

def extractMinimal(Set, Function, Args):
    # This function returns a minimal subset of a given set that returns TRUE when
    # the Function is called with it as the first argument and the other given Args
    # Note that the processing is performed in the reverse order for simplicity.
    # NOTE: This function assumes that the Set itself returns TRUE for the Function.
    for i in list(reversed(list(range(len(Set))))):
        minusI = Set[:i] + Set[i+1:]
        if Function(minusI, *Args):
            Set = minusI
    return Set

def findSmallCutsets(Network, Target, skip = [], Kmax = 3):
    # This function systematically finds all cutsets of size at most 3 in a REVERSIBLE network
    # It exploits the property that X is a cutset for i if and only if X-j+i is a cutset for j
    # skip is a subset of indices that are omitted from the triplet testing; for debugging only!
    m,n = getSize(Network)
    Subsets = [[] for i in range(min(Kmax,3))]
    # Step 0: Reduce the network by removing the Target reaction
    curcols = [x for x in range(n) if x != Target]
    n = len(curcols)
    redNetwork = [[Network[x][y] for y in curcols] for x in range(m)]
    # Step 1: Identify the cutsets of size 1 (NOTE: if the matrix is reduced, there will be none)
    (SBlocked, B) = findSBlocked(redNetwork, True)
    Subsets[0] = [curcols[x] for x in SBlocked]
    kept = [x for x in range(n) if x not in SBlocked]
    curcols = [curcols[x] for x in kept]
    n = len(curcols)
    redNetwork = [[redNetwork[x][y] for y in kept] for x in range(m)]
    curB = [B[x] for x in kept]
    # Step 2: identify subsets of size 2 (enzyme subsets in the reduced network)!
    ESubsets = findSubsets(curB, True)
    Pairs = []
    for (key, value) in ESubsets.items(): # converting to original indices
        Pairs.append([[curcols[key]] + [curcols[x[0]] for x in value]])
    Pairs = sum(Pairs, [])
    Subsets[1] = sum([[[x,y] for x in pair for y in pair if x < y] for pair in Pairs],[])
    if Kmax == 3:
        # Step 3: identify subsets of size 3 (enzyme subsets with one reaction removed)!
        for i in range(n):
            print(('Processing reaction ' + str(i)))
            if i not in skip:
                cur = curcols[i]
                redNetworkT = [[redNetwork[x][y] for y in range(n) if y != i] for x in range(m)]
                ESubsetsT = findSubsets(redNetworkT)
                PairsT = []
                for (key, value) in ESubsetsT.items(): # converting to original indices
                    PairsT.append([[curcols[key + int(key >= i)]] + [curcols[x[0] + int(x[0] >= i)] for x in value]])
                    # NOTE: We are shifting the indices for subsets after the i-th position
                PairsT = sum(PairsT, [])
                Subsets[2] += sum([[sorted([cur,x,y]) for x in pair for y in pair if x != y] for pair in PairsT],[])
    return Subsets
####parallelization
def parallelizedTestCutSet(Network, Target,  rec, I,subList,index,out_q_lethal,out_q_iter):
    Filename = 'lethal.lp'
    count = len(subList) * index
    outdict_lethal = {}
    outdict_iter = {}
    subIter = 0
    for ind, pair in enumerate(subList):
        if testCutSet(pair, Network, Target, (Filename[:-3] + str(index) + 'V' + str(ind) + Filename[-3:]), rec, I):
            outdict_lethal[count+ind] = pair #place pair into dictionary
        subIter += 1

    #place items into dictionary then into queue
    outdict_iter[index] = subIter
    out_q_lethal.put(outdict_lethal)
    out_q_iter.put(outdict_iter)
###end parallelization

def findEssentialLethal(Network, Target, Filename = 'lethal.lp', rec = True, I = [], verbose = False, parallel = 0):
    # This function identifies all essential and synthetic lethal pairs in a given
    # network for a specified target reaction.
    # It returns a tuple, containing a list of singletons and a list of pairs.
    # The network is assumed to be irreversible unless rec is specified to be False.
    # In that case, the reactions considered to be irreversible should be specified in I.
    # If parallel > 0, it specifies the number of cores - the checks are split into as many parts.
    m, n = getSize(Network)
    Essential, Lethal = [], []
    if n == 1:
        return (Essential, Lethal)
    (val, vec) = findMin1Norm(Network, Target, [1]*n, [], [], 1e-5, Filename, 'null', rec, I)
    firstEntries = [int(y[1:]) for y in list(vec.keys()) if y.startswith('V')]
    L = len(firstEntries)
    # start by checking for essentiality; if a reaction is not essential, add the feasible vector to the collection!
    Collection = [firstEntries]
    Iter = 0
    for Iter, entry in enumerate(firstEntries):
        if entry != Target:
            (valN, vecN) = findMin1Norm(Network, Target, [1]*n, [entry], [], 1e-5, Filename[:-3] + 'I' + str(Iter) + '.lp', 'null', rec, I)
            if valN: # feasible
                Collection.append([int(y[1:]) for y in list(vecN.keys()) if y.startswith('V')])
            else:
                Essential.append(entry)
    # eliminate Target reaction and essentials from the collection
    Collection = [set(vector).difference(Essential + [Target]) for vector in Collection]
    Collection = [_f for _f in Collection if _f]
    # now look for synthetic lethal candidates among the collection
    allCandidates = sorted(list(reduce(lambda x,y:x.union(y), Collection, set([]))))
    CandidatePairs = [[x,y] for x in allCandidates for y in allCandidates if x < y]
    if verbose:
        print(('There are ' + str(len(CandidatePairs)) + ' pairs to be processed'))
    #parallel = 6 # testing
    if parallel > 0:
        #chunkify
        CandidatePairs = [pair for pair in CandidatePairs if all([(pair[0] in z or pair[1] in z) for z in Collection])]
        splitPairs = [CandidatePairs[i::parallel] for i in range(parallel)]


        out_q_lethal = multiprocessing.Queue() #each worker process puts a result dictionary into it
        out_q_iter = multiprocessing.Queue() #each worker process puts a result dictionary into it
        procs = []
        for index, subList in enumerate(splitPairs):
            #creates process, target function and parameters
            p = multiprocessing.Process(
                    target=parallelizedTestCutSet,
                    args=(Network, Target, rec, I,subList,index,out_q_lethal,out_q_iter ))
            procs.append(p)
            p.daemon = True # runs as a background process
            p.start()

        # Collect all results into a single result dict
        resultdict_lethal = {}
        for i in range(parallel):
            resultdict_lethal.update(out_q_lethal.get())

        resultdict_iter = {}
        for i in range(parallel):
            resultdict_iter.update(out_q_iter.get())

        # Wait for all worker processes to finish
        for p in procs:
            p.join()
        Lethal = list(resultdict_lethal.values())
        Iter += sum(resultdict_iter.values())

    else:
        for ind, pair in enumerate(CandidatePairs):
            if verbose and (ind + 1) % 1000 == 0:
                print(('Processed ' + str(ind + 1) + ' pairs so far'))
            if all([(pair[0] in z or pair[1] in z) for z in Collection]):
                if testCutSet(pair, Network, Target, Filename[:-3] + str(Iter) + Filename[-3:], rec, I):
                    Lethal.append(pair)
                Iter += 1
    if verbose:
        print(("This required a total of " + str(Iter) + " linear programs"))

    return (Essential, Lethal)

def testSubsets(Network, Target, fluxes = [], Kmax = 3, Filename = 'trial.lp', rec = True, I = [], startInd = 0, startSubsets = []):
    # This function systematically tests all subsets of size at most Kmax to see
    # if they are cutsets of a given network for a specified target reaction.
    # It returns a list of all minimal such subsets in order of increasing size.
    # The network is assumed to be irreversible unless rec is specified to be False.
    # In that case, the reactions considered to be irreversible should be specified in I.
    # If the processing had already been partly done, change startInd to reflect that!
    # In that case, the optional argument startSubsets should contain the initial set.
    m,n = getSize(Network)
    Subsets = [[] for i in range(Kmax)]
    # Use initialization if necessary!
    if startSubsets:
        for i in range(len(startSubsets)):
            Subsets[i] = startSubsets[i]
    Iter = 0
    if Kmax >= 4:
        print('This function is not implemented for Kmax >= 4 as it is likely to take too long!')
        return Subsets
    if not fluxes:
        fluxes = findEFM(Network, Target, [], Filename, rec, I)
    indices = [[int(y[1:]) for y in list(x.keys()) if y.startswith('V')] for x in fluxes]
    relevant = sorted(list(set(sum(indices, [])).difference([Target])))
    additional = [x for x in range(n) if x not in relevant and x != Target]
    L = len(fluxes)
    G = len(relevant)
    H = len(additional)
    print(('There are ' + str(L) + ' fluxes containing ' + str(G) + ' reactions out of a total ' + str(n)))
    # Construct matrix M with M[i,j] = number of co-occurrences of i and j
    M = {}
    for x in relevant:
        for y in relevant:
            M[x,y] = 0
    for i in range(L):
        indices[i].remove(Target)
        for x in indices[i]:
            for y in indices[i]:
                M[x,y] += 1
    if Kmax >= 1:
        print('Processing subsets of size 1')
        for index, subset in enumerate(relevant):
            print(('Processing element number ' + str(index)))
            if M[subset,subset] == L:
                if Iter >= startInd:
                    result = testCutSet([subset], Network, Target, Filename[:-3] + str(Iter) + Filename[-3:], rec, I)
                    if result:
                        Subsets[0].append([subset])
                    else:
                        # see if any extension of size 1 or 2 would work, if necessary
                        if Kmax >= 2:
                            simpleExtensions = []
                            print(('Checking ' + str(H) + ' possible extensions of size 1'))
                            for ind, item in enumerate(additional):
                                if (ind % 25 == 0):
                                    print(('Checked ' + str(ind) + ' options so far'))
                                result = testCutSet([subset]+[item], Network, Target, Filename[:-3] + str(Iter) + Filename[-3:], rec, I)
                                if result:
                                    Subsets[1].append(sorted([subset] + [item]))
                                    simpleExtensions.append(item)
                                Iter += 1
                        if Kmax >= 3:
                            candidates = [x for x in additional if x not in simpleExtensions]
                            candidatePairs = [[x,y] for x in candidates for y in candidates if x < y]
                            print(('Checking ' + str(len(candidatePairs)) + ' possible extensions of size 2'))
                            for ind, pair in enumerate(candidatePairs):
                                if (ind % 25 == 0):
                                    print(('Checked ' + str(ind) + ' options so far'))
                                result = testCutSet([subset] + pair, Network, Target, Filename[:-3] + str(Iter) + Filename[-3:], rec, I)
                                if result:
                                    Subsets[2].append(sorted([subset] + pair))
                                Iter += 1
                    Iter += 1
    if Kmax >= 2:
        print('Processing subsets of size 2')
        candidates = [x for x in relevant if [x] not in Subsets[0]]
        candidatePairs = [[x,y] for x in candidates for y in candidates if x < y]
        print(('There are ' + str(len(candidatePairs)) + ' subsets to process'))
        for index, subset in enumerate(candidatePairs):
            print(('Processing pair number ' + str(index)))
            p, q = subset[0], subset[1]
            if M[p,p] + M[q,q] - M[p,q] == L:
                if Iter >= startInd:
                    result = testCutSet(subset, Network, Target, Filename[:-3] + str(Iter) + Filename[-3:], rec, I)
                    if result:
                        Subsets[1].append(subset)
                    else:
                        # see if any extension of size 1 would work
                        if Kmax >= 3:
                            print(('Checking ' + str(H) + ' possible extensions of size 1'))
                            for ind, item in enumerate(additional):
                                if (ind % 25 == 0):
                                    print(('Checked ' + str(ind) + ' options so far'))
                                result = testCutSet(subset + [item], Network, Target, Filename[:-3] + str(Iter) + Filename[-3:], rec, I)
                                if result:
                                    Subsets[2].append(sorted(subset + [item]))
                                Iter += 1
                Iter += 1
    if Kmax >= 3:
        print('Processing subsets of size 3')
        candidates = [x for x in relevant if [x] not in Subsets[0]]
        candidateTriples = [[x,y,z] for x in candidates for y in candidates for z in candidates if x < y and y < z]
        for index, subset in enumerate(candidateTriples):
            print(('Processing triplet number ' + str(index)))
            p, q, r = subset[0], subset[1], subset[2]
            if [p,q] in Subsets[1] or [p,r] in Subsets[1] or [q,r] in Subsets[1]:
                continue
            if M[p,p] + M[q,q] + M[r,r] - M[p,q] - M[p,r] - M[q,r] +  min([M[p,q], M[p,r], M[q,r]]) >= L:
                if Iter >= startInd:
                    result = testCutSet(subset, Network, Target, Filename[:-3] + str(Iter) + Filename[-3:], rec, I)
                    if result:
                        Subsets[2].append(subset)
                Iter += 1
    print(("This required a total of " + str(Iter) + " linear programs"))
    return Subsets

def generateSubsets(n, k):
    # This function generates all k-subsets of range(n) in lexicographic order
    # http://snipplr.com/view/17859/generate-all-ksubsets-of-an-nset-sequentially/
    m, h = 0, k
    a = list(range(k))
    yield a
    while True:
        if m < n-h:
            h = 1
        else:
            h += 1
        m = a[k-h] + 1
        for j in range(h):
            a[k+j-h] = m + j
        yield a
        if a[0] == n - k:
            break

def findMinAdded(N, special, weight = [1], exclude = [], eps = 1e-5, Filename = 'trial.lp', extra = {}, option = 'row', Cplex = False):
    # This function finds the nonnegative vector of smallest weight to be added to a vector
    # in the rowspace of N (nullspace if option is 'null') to make the resulting vector >=0.
    # The entry corresponding to special is 1. Extra lower bounds may be supplied as extra.
    # Note: if the weight vector has a single component, it is automatically taken to be 1!
    m, n = getSize(N)
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Minimize\n'
    opt += 'obj: '
    if len(weight) == n:
        opt += ' + '.join([str(weight[i]) + ' T' + str(i) for i in range(n) if weight[i]])
    else:   # equal weights, take all of them equal to 1
        opt += ' + '.join(['T' + str(i) for i in range(n)])
    opt += '\n'
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    if option == 'row':
        const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(i) for i in range(m) if N[i][j]]) + ' + -1 Y' + str(j) + ' = 0\n' for j in range(n)])
        const += ''.join(['Y' + str(j) + ' + T' + str(j) + '>= 0\n' for j in range(n)])
        const += 'Y' + str(special) + ' = 1\n'
    elif option == 'null':
        const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
        const += ''.join(['X' + str(j) + ' + T' + str(j) + '>= 0\n' for j in range(n)])
        const += 'X' + str(special) + ' = 1\n'
    if exclude:
        k = len(exclude)
        s = len(exclude[0])
        if s == n:
            const += ''.join([' + '.join(['T' + str(j) for j in range(n) if not exclude[i][j]]) + ' >' + str(eps) + '\n' for i in range(k)])
        else:
            print("Error: the length of the excluded vectors is incorrect!")
    bounds = 'Bounds\n'
    if option == 'row':
        bounds += ''.join(['X' + str(i) + ' free\n' for i in range(m) if [_f for _f in N[i] if _f]])
        bounds += ''.join(['Y' + str(j) + ' free\n' for j in range(n) if j not in extra])
        ## AL; changed order of bounds to make QSopt_ex not complain
        # bounds += ''.join(['Y'+str(j)+'<='+str(round(extra[j],5)) + '\n' for j in extra])
        print(''.join([str(round(extra[j],5))+'<=' + 'Y' + str(j)+'<= 0'+ '\n' for j in extra]))
        bounds += ''.join([str(round(extra[j],5))+'<=' + 'Y' + str(j)+'<= 0\n' for j in extra])
        print([round(extra[j],5) for j in extra])
        #note that this puts the automatic upper bound that Yj<=0 as per QSopt_ex defaults
        #const += ''.join(['Y' + str(j) + ' >= ' + str(round(extra[j],5)) + '\n' for j in extra])
    elif option == 'null':
        bounds += ''.join(['X' + str(j) + ' free\n' for j in range(n) if [_f for _f in [N[i][j] for i in range(m)] if _f]])
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    return processFile(Filename, True,destroyIn=False,destroyOut=False)

def checkUnblocked(RemoveConst, N, Irr, growth, Filename = 'Unblock.lp'):
    # This function checks whether removing a given subset of constraints
    # produces a metabolic network that allows for positive growth.
    # It returns TRUE if that is the case, and FALSE otherwise.
    Nred = [N[x] for x in range(len(N)) if x not in RemoveConst]
    (val, vec) = findFeasible(Nred, growth, Irr, True, Filename)
    return not (type(val) == type([]) and len(val) == 0) # TRUE IFF THE PROBLEM IS FEASIBLE

def minimalUnblock(N, Irr, growth, algo = 'norm', Filename = 'Unblock.lp'):
    # This function uses a heuristic to compute the smallest set of metabolites
    # that need to be unbalanced in order to produce growth in a blocked model.
    # N is the stoichiometric matrix, Irr are the irreversible reactions, growth
    # is the reaction to be unblocked. The returned set is only inclusion-minimal!
    # Algorithms currently implemented (specified by option algo):
    # greedy (removes metabolites in order of decreasing degree in the network);
    # norm (minimizes the 1-norm of the metabolite balance vector defined as Nv);
    # sense (same as norm, but with iterative weight updates until convergence);
    # redund (return a set of "redundant" metabolites plus an extra relevant one);
    # full (exhaustive enumeration of all possible reaction subsets of size <= 5).
    m, n = getSize(N)
    if algo == 'greedy':
        Iter = 0
        degrees = getDegrees(N)
        order = sorted(list(range(m)), key = lambda x:degrees[x])
        # NOTE: The sorting is in increasing order, but removal is from the end!
        return extractMinimal(order, checkUnblocked, [N, Irr, growth, Filename])
    elif algo == 'norm' or algo == 'sense':
        Ntemp = [[] for i in range(m)]
        for i in range(m):
            Ntemp[i] = [0]*n
            for j in range(n):
                Ntemp[i][j] = N[i][j]
        Rev = [x for x in range(n) if x not in Irr]
        n += len(Rev)
        for x in range(m):
            Ntemp[x] += [-Ntemp[x][y] for y in Rev]
        if growth in Rev:
            j = Rev.index(growth)
            z = [n - len(Rev) + j]
        else:
            z = []
        (val, vec) = findMin1Norm(Ntemp, growth, [1], z, [], Filename = Filename, option = 'col')
        for x in list(vec.keys()):
            if not x.startswith('T'):
                del vec[x]
        if algo == 'sense':
            L = len(vec)
            eps = min(list(vec.values()))
            Iter = 0
            while(True):
                weight = [one / eps] * m
                for x in vec:
                    if x.startswith('T'):
                        weight[int(x[1:])] = one / (eps + vec[x])
                (val, vec) = findMin1Norm(Ntemp, growth, weight, [], [], Filename = Filename[:-3] + 'I' + str(Iter) + '.lp', option = 'col')
                for x in list(vec.keys()):
                    if not x.startswith('T'):
                        del vec[x]
                if len(vec) == L:
                    break
                else:
                    L = len(vec)
                    Iter += 1
        return extractMinimal([int(x[1:]) for x in vec], checkUnblocked, [N, Irr, growth, Filename])
    elif algo == 'redund':
        # NOTE: this option ignores the thermodynamic constraints, strictly speaking!
        nonRed = findRedundant(N)
        Nred = [N[x] for x in nonRed]
        Vec = [0]*n
        Vec[growth] = 1
        (val, vec) = vectorInSpan(Nred, Vec, Filename)
        if len(val):
            opts = [int(x[1:]) for x in list(vec.keys()) if x.startswith('X')]
            return [x for x in range(m) if x not in nonRed and [_f for _f in N[x] if _f]] + [opts[0]]
        else:
            print("Error: the given reaction wasn't stoichiometrically blocked in the first place")
            return []
    elif algo == 'full':
        found = False
        active = [x for x in range(m) if [_f for _f in N[x] if _f]]
        q = len(active)
        for k in range(1,5):
            if not found:
                Iter = 0
                for subset in generateSubsets(q,k):
                    Subset = [active[x] for x in subset]
                    Ntemp = [N[x] for x in range(m) if x not in Subset]
                    if len(findFeasible(Ntemp, growth, Irr, True, Filename[:-3] + str(k) + 'I' + str(Iter) + '.lp')[0]):
                        found = True
                        break
                    else:
                        Iter += 1
            else:
                break
        return Subset
    else:
        print('Error: unrecognized option!')
        return []

def checkThermoUnblocked(RemoveConst, N, Irr, growth, Filename = 'Unblock.lp'):
    # This function checks whether removing a given subset of irreveresibility
    # constraints produces a metabolic network that allows for positive growth.
    # It returns TRUE if that is the case, and FALSE otherwise.
    IrrRed = [x for x in Irr if x not in RemoveConst]
    (val, vec) = findFeasible(N, growth, IrrRed, True, Filename)
    return not (type(val) == type([]) and len(val) == 0) # TRUE IFF THE PROBLEM IS FEASIBLE

def minimalThermoUnblock(N, Irr, growth, Filename = 'Unblock.lp'):
    # This function uses a heuristic to compute the smallest set of reactions
    # that need to be made reversible in order to allow growth in a given model.
    # N is the stoichiometric matrix, Irr are the irreversible reactions, growth
    # is the reaction to be unblocked. The returned set is only inclusion-minimal!
    # NOTE: If the network is not thermodynamically unblockable, returns None!
    n = len(N[0])
    weight = [int(i in Irr) for i in range(n)]
    (val, vec) = findMinAdded(N, growth, weight, option = 'null')
    if (type(val) == type([]) and len(val) == 0):
        print('The growth reaction is in fact stoichiometrically blocked!')
        return
    else:
        vec = [int(x[1:]) for x in vec if x.startswith('T')]
        used = [x for x in vec if x in Irr]
        return extractMinimal(used, checkThermoUnblocked, [N, Irr, growth, Filename])

def getDegrees(N, opt = 'metabs', direct = 0):
    # This function returns the set of degrees associated with the given network.
    # If opt = 'metabs', degrees of metabolites; if opt = 'reacts', of reactions.
    # If direct is set to a positive number count only the positive coefficients;
    # if to a negative number only the negative ones; otherwise count everything.
    if opt == 'metabs':
        if direct < 0:
            return [len([x for x in x if x < 0]) for x in N]
        elif direct > 0:
            return [len([x for x in x if x > 0]) for x in N]
        else: # direct = 0
            return [len([_f for _f in x if _f]) for x in N]
    elif opt == 'reacts':
        if direct < 0:
            return [len([x for x in [N[x][y] for x in range(len(N))] if x < 0]) for y in range(len(N[0]))]
        elif direct > 0:
            return [len([x for x in [N[x][y] for x in range(len(N))] if x > 0]) for y in range(len(N[0]))]
        else:
            return [len([_f for _f in [N[x][y] for x in range(len(N))] if _f]) for y in range(len(N[0]))]
    else:
        print('Error: unrecognized option')
        return

def checkEnergyBalance(N, EFMs, Internal):
    # Returns a boolean array corresponding to given EFMs (True if feasible, False if not)
    # See the paper for a justification of the linear program used for checking purposes.
    m, n = getSize(N)
    L = len(EFMs)
    checks = [False]*L
    Nred = [[N[i][j] for j in Internal] for i in range(m)]
    for i in range(L):
        support = [x for x in Internal if EFMs[i][x]]
        (val, vec) = findPosSupport(Nred, support, Filename = 'EBA' + str(i) + '.lp', Min = 1, restricted = False, option = 'row')
        if len(val):
            checks[i] = True
    return checks

def findMinimalMedia(N, growth, Exchange, rec = False, I = [], opt = True):
    # This function attempts to find minimal several minimal media allowing growth
    # for a given stoichiometric matrix and a given subset of "exchange" reactions.
    n = len(N[0])
    MediaPos = []
    MediaVec = []
    weight = [0]*n
    new = True
    Iter = 0
    for i in Exchange:
        weight[i] = 1
    while (new):
        (val, vec) = findMin1Norm(N, growth, weight, [], MediaPos, 1e-5, 'Media' + str(Iter) + '.lp', 'null', rec, I)
        vec = peelOff(vec, MediaVec)
        if vec:
            MediaVec.append(vec)
            active = [int(x[1:]) for x in vec if x.startswith('V')]
            relevant = [False]*n
            for x in Exchange:
                if x in active:
                    relevant[x] = True
            if opt:
                y = findTrueIndices(relevant)
                minMedium = extractMinimal(y, checkMedium, [N, growth, Exchange, rec, I])
                relevant = [(True if z in minMedium else False) for z in range(n)]
            if relevant not in MediaPos:
                MediaPos.append(relevant)
                Iter += 1
            else:
                new = False
        else:
            new = False
    Media = [[y for y in range(len(x)) if x[y]] for x in MediaPos]
    return Media

def checkMedium(Set, N, growth, Exchange, rec = False, I = []):
    # returns True iff the given set of exchange reactions allows growth
    m, n = getSize(N)
    if rec:
        I = list(range(n))
    block = [x for x in Exchange if x not in Set]
    (val, vec) = findFeasible(N, growth, Irrev = I, pos = True, Filename = 'trial.lp', disable = block, negative = [])
    return not (type(val) == type([]) and len(val) == 0)

def GaussJordan(N, pivoting = True, Gauss = False, verbose = False):
    # This function returns the (pivoted) Row-Reduced Echelon Form of N.
    # It was adapted from http://adorio-research.org/wordpress/?p=192.
    # pivoting = True for partial row pivoting, False if no pivoting.
    # Gauss = True if a Row Echelon Form of N is needed, False otherwise.
    # Returns (A, R); A: transformed input matrix; R: indices of pivots.
    nrow, ncol = getSize(N)
    A = [[convertToFraction(N[i][j]) for j in range(ncol)] for i in range(nrow)]
    R = []
    # Triangularization
    currow = 0
    # This holds the index of the row we are currently searching
    if verbose:
        print(('There are ' + str(ncol) + ' columns to be processed'))
    for pivot in range(ncol):
        if verbose and pivot % 10 == 0:
            print(('Processing column ' + str(pivot)))
        if pivoting:
            absVals = [abs(A[x][pivot]) for x in range(currow, nrow)]
            bestPivot = max(absVals)
            best = currow + absVals.index(bestPivot)
            # exchange pivot row with best row.
            if currow != best:
                A[best], A[currow] = A[currow], A[best]
        m = A[currow][pivot]
        if m == zero:
            # In solving a system, this would mean the algorithm broke down;
            # however, in our case, this simply means we found a zero column!
            continue
        else:
            R.append(pivot)
            for row in range(0, nrow):
                # NOTE: For Gaussian elimination, only eliminate the rows below the current one!
                if ((Gauss and row > currow) or (not Gauss and row != currow)) and (A[row][pivot] != zero):
                    kmul = A[row][pivot] / m
                    # Apply rectangular rule (Reduction)
                    for col in range(pivot + 1, ncol):
                        A[row][col] -= (kmul * A[currow][col])
                    A[row][pivot] = zero
            currow += 1
            if currow == nrow:
                # we are all the way down!
                break
    if not Gauss:
        # NOTE: For Gaussian elimination, we also don't need to normalize the rows at the end!
        nonPivots = [x for x in range(ncol) if x not in R]
        for row in range(len(R)):
            pivot = R[row]
            m = A[row][pivot]
            A[row][pivot] = one
            for col in nonPivots:
                if A[row][col]:
                    A[row][col] /= m
    return (A, R)

def NullspaceBasis(N, GJ = False):
    # This function computes the nullspace basis using a Gauss-Jordan elimination
    if not GJ:
        (A, R) = GaussJordan(N)
    else:
        (A, R) = N
    m = len(R)
    if len(A) > 0:
        n = len(A[0])
    else:
        n = 0
    Rcomp = [x for x in range(n) if x not in R]
    if m > n:
        print ('Error: this cannot happen!')
    elif m == n:
        return []
    NullBasis = [[] for i in range(n)]
    for i in range(n):
        if i in R:
            ind = R.index(i)
            NullBasis[i] = [-A[ind][j] for j in Rcomp]
        elif i in Rcomp:
            ind = Rcomp.index(i)
            NullBasis[i] = [zero]*(n - m)
            NullBasis[i][ind] = one
    return NullBasis

def transpose(Matrix):
    # This function creates the transpose of an input matrix
    m, n = getSize(Matrix)
    return [[Matrix[i][j] for i in range(m)] for j in range(n)]

def ReadMatrix(Filename):
    # reads a sparse matrix file and converts contents into fractions
    f = open(Filename, 'r')
    g = f.readlines()
    f.close()
    I = [int(x) for x in g[0].strip()[1:-2].split(',')]
    J = [int(x) for x in g[1].strip()[1:-2].split(',')]
    S = [float(x) for x in g[2].strip()[1:-2].split(',')]
    m = max(I)
    n = max(J)
    Mat = [[] for i in range(m)]
    for i in range(m):
        Mat[i] = [0]*n
    for k in range(len(I)):
        # NOTE: Matlab numbering starts at 1!
        Mat[I[k]-1][J[k]-1] = S[k]
    return Mat

def WriteASCIIMatrix(Matrix, Filename = 'Matrix.txt'):
    # Creates an ASCII representation of an input matrix in a file with the specified name.
    # WARNING: The matrix MUST contain ONLY integers; otherwise the result will be wrong...
    m, n = getSize(Matrix)
    output = '\n'.join(['\t'.join([str(int(Matrix[x][y])) for y in range(n)]) for x in range(m)])
    f = open(Filename, 'w')
    f.write(output)
    f.close()
    return

def ReadASCIIMatrix(Filename, Rows = False, sep = '\t', Blocks = False):
    # Reads an ASCII matrix from a file with the specified name.
    # If Rows is True, looks for row identifiers at the end of each row.
    # In that case, the rows are reordered before the matrix is returned!
    # To change the separator which is used to split the entries use sep.
    # If Blocks is True, the parser reads and returns a block structure!
    f = open(Filename, 'r')
    g = f.readlines()
    f.close()
    L = len(g)
    if Blocks:
        # Start by getting the blocks
        endBlocks = [-1] + [x for x in range(L) if g[x].startswith('---')] + [L]
        B = len(endBlocks) - 1
        allBlocks = [[] for i in range(B)]
        allRows = [[] for i in range(B)]
        # Then get the row indices and delete them
        rInds = [g[x].rfind('R:') for x in range(L)]
        rowInds = [int(g[x][rInds[x] + len('R:'):]) - 1 for x in range(L) if x not in endBlocks]
        # NOTE: The numbering starts from 0 in Python, but from 1 in the BlocDiag program!
        g = [g[x][:rInds[x]] for x in range(L)]
        g = [x.split('|\t')[:-1] for x in g]
        # Finally, split the matrix into the relevant blocks and return
        for i in range(B):
            allRows[i] = rowInds[(endBlocks[i]+1-i):(endBlocks[i+1]-i)]
            curBlock = g[(endBlocks[i]+1):endBlocks[i+1]]
            if i == B-1 and (not curBlock or len(curBlock[0]) <= i): # the last block is empty!
                curBlock = []
            else:
                curBlock = [x[i].split(sep)[:-1] for x in curBlock]
                curBlock = [[int(y+'0')/10 for y in x] for x in curBlock]
                # Note that in this way, the empty string is converted to a 0
            allBlocks[i] = curBlock
        return(allBlocks, allRows)
    else:
        g = [x.strip().split(sep) for x in g]
        if Rows:
            rowInds = [x[-1] for x in g]
            g = [x[:-1] for x in g]
        # TEMPORARY: Accommodating integers written as a float (MATLAB notation)
        Matrix = [[float(x) for x in y] for y in g]
        # Remove any empty rows at the end!
        while not Matrix[-1]:
            Matrix = Matrix[:-1]
        if Rows:
            rowInds = [int(x[x.find(':')+1:]) for x in rowInds]
            order = [x[0] for x in sorted([(i,rowInds[i]) for i in range(len(rowInds))], key=lambda x:x[1])]
            Matrix = [Matrix[order[i]] for i in range(len(order))]
        return Matrix

def findReaction(cols, Subsets, origReaction):
    # Finds the index of the column corresponding to the origReaction reaction
    survivors = [x for x in range(len(cols)) if cols[x] == 0]
    if origReaction in survivors:
        return survivors.index(origReaction)
    elif cols[origReaction] == 5: # enzyme subset
        subReacts = [[z[0] for z in x] for x in Subsets]
        subset = [x for x in range(len(subReacts)) if origReaction in subReacts[x]]
        if len(subset) != 1:
            print('Error: incorrectly identified subset!')
            return []
        else:
            return survivors.index(subReacts[subset[0]][0])
    else:
        print("Sorry, couldn't find the desired reaction")
        return []

def SubspaceIntersection(ListOfBases, NullB = False):
    # Returns the basis of the intersection of a set of vector spaces
    # given by a list of their basis vectors (assumed to be row vectors)
    # If NullB is TRUE, returns a bais for the intersection's nullspace.
    L = len(ListOfBases)
    Basis = []
    if L == 0:
        print('Error: no bases supplied!')
        return
    else:
        n = len(ListOfBases[0][0])
    for i in range(L):
        curBasis = ListOfBases[i]
        if len(curBasis[0]) != n:
            print(('Error: the basis number ' + str(i) + ' has incorrect dimensions!'))
            return
        else:
            curNull = NullspaceBasis(ListOfBases[i])
            if curNull:
                Basis += transpose(curNull)
    if NullB:
        return Basis
    else:
        return transpose(NullspaceBasis(Basis))

def findEssential(Network, growth, Exchange, allowed, Filename = 'trial.lp', rec = True, I = [], forbidden = [], GeneDeletes = [], J = []):
    # This procedure finds all the reactions essential for growth in a medium defined
    # by the set of exchange reactions and the subset of allowed exchange reactions.
    # The computations are performed using files derived from the specified Filename.
    # The network is assumed to be irreversible unless rec is specified to be False.
    # In that case, the reactions considered irreversible should be specified in I.
    # If any additional reactions are forbidden, they should be specified at the end.
    # The last argument can be dictionary of reactions disabled by each gene deletion,
    # and in that case the procedure finds all the genes essential for growth instead.
    # Latest addition: there is now an option J for reactions with negative-only flux.
    m = len(Network)
    n = len(Network[0])
    Forbidden = [x for x in (Exchange + forbidden) if x not in allowed]
    if not rec:
        Irrev = I
    else:
        Irrev = list(range(n))
    Essential = []
    (val0, vec0) = findFeasible(Network, growth, Irrev, True, Filename, Forbidden, J)
    if type(val0) == type([]) and val0 == []:
        print("Error: the organism cannot grow under the given condition!")
        return
    else:
        active = [int(x[1:]) for x in vec0 if x.startswith('V')]
        active.remove(growth)
    if GeneDeletes:
        i = 0
        for gene in GeneDeletes:
            (val, vec) = findFeasible(Network, growth, Irrev, True, Filename[:-3] + str(i) + '.lp', Forbidden + GeneDeletes[gene], J)
            i = i + 1
            if type(val) == type([]) and val == []:
                Essential.append(gene)
    else:
        for i in active:
            (val, vec) = findFeasible(Network, growth, Irrev, True, Filename[:-3] + str(i) + '.lp', Forbidden + [i], J)
            if type(val) == type([]) and val == []:
                Essential.append(i)
    return Essential

def classifyExchange(FullNetwork, externalMetabs, Irrev, extra = False):
    # This function finds all the exchange reactions in a given full network
    # (with the external metabolites included) and classifies them into 3 types:
    # input only, output only, and mixed, and additionally by irreversibility.
    # If extra is True, also classifies reactions with a single internal metabolite.
    m, n = getSize(FullNetwork)
    Exchange = [y for y in range(n) if any([FullNetwork[x][y] for x in externalMetabs])]
    InputIrr, InputRev, OutputIrr, OutputRev, MixedIrr, MixedRev = [], [], [], [], [], []
    for react in Exchange:
        In, Out = False, False
        for x in externalMetabs:
            if FullNetwork[x][react] < 0:
                In = True
            elif FullNetwork[x][react] > 0:
                Out = True
        if In and Out:
            if react in Irrev:
                MixedIrr.append(react)
            else:
                MixedRev.append(react)
        elif In and (not Out):
            if react in Irrev:
                InputIrr.append(react)
            else:
                InputRev.append(react)
        elif (not In) and Out:
            if react in Irrev:
                OutputIrr.append(react)
            else:
                OutputRev.append(react)
        else:
            print('This should never happen!')
    if extra: # identify non-exchange reactions with a single metabolite (must be internal)
        other = [y for y in set(range(n)).difference(Exchange) if len(list(filter(None, [FullNetwork[x][y] for x in range(m)]))) == 1]
        for react in other:
            coeff = list(filter(None, [FullNetwork[x][react] for x in range(m)]))[0]
            if coeff > 0:
                if react in Irrev:
                    OutputIrr.append(react)
                else:
                    OutputRev.append(react)
            else:
                if react in Irrev:
                    InputIrr.append(react)
                else:
                    InputRev.append(react)
    return (InputIrr, InputRev, OutputIrr, OutputRev, MixedIrr, MixedRev)

def findFreeLunch(N, Irrev, weight = [1], freeMetabs = [], Filename = 'trial.lp', Cplex = False):
    # This function finds the vector optimizing a given weight in the column space of N with
    # nonnegative entries and whose flux vector satisfies the irreversibility conditions.
    # The free metabolites are ones that can be considered given (i.e. they can be consumed).
    # NOTE: To get meaningful results, the input matrix should contain external metabolites!
    m, n = getSize(N)
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Maximize\n'
    opt += 'obj: '
    if len(weight) == m:
        opt += ' + '.join([str(weight[i]) + ' Y' + str(i) for i in range(m) if weight[i]])
    elif len(weight) == 1:   # equal weights, take all of them equal to 1
        opt += ' + '.join(['Y' + str(i) for i in range(m)])
    else:
        print('Error: the weight vector is not of the right length!')
    opt += '\n'
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(j) for j in range(n) if N[i][j]]) + ' + -1 Y' + str(i) + ' = 0\n' for i in range(m)])
    bounds = 'Bounds\n'
    bounds += ''.join(['X' + str(j) + ' free\n' for j in range(n) if j not in Irrev])
    bounds += ''.join(['-1 <= Y' + str(i) + ' <= 1\n' for i in freeMetabs])
    bounds += ''.join(['Y' + str(i) + ' <= 1\n' for i in range(m) if i not in freeMetabs])
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    return processFile(Filename, True)

def FBA(N, growth, Exchange, allowed, limits = [1], Filename = 'trial.lp', rec = True, I = [], forbidden = [], Negative = [], Cplex = False):
    # This procedure finds the maximal growth rate of an organism in a medium defined
    # by the set of exchange reactions and the subset of allowed exchange reactions.
    # The bounds on the flux of each allowed reaction are given by the vector limits.
    # If limits is a scalar, all the limits are assumed to be 1. 'inf' is acceptable!
    # The network is assumed to be irreversible unless rec is specified to be False.
    # In that case, the reactions considered irreversible should be specified in I.
    # If any additional reactions are forbidden, they should be specified at the end.
    # Latest addition: there is now a Negative for reactions with negative-only flux.
    m, n = getSize(N)
    Forbidden = [x for x in (Exchange + forbidden) if x not in allowed]
    if not rec:
        Irrev = I
    else:
        Irrev = list(range(n))
    Rev = [x for x in range(n) if x not in Irrev]
    intro = 'Problem\n'
    intro += (Filename.replace('.lp', '') + '\n')
    opt = 'Maximize' + '\n' +  'V' + str(growth) + '\n'
    const = 'Subject' + ' To' * int(Cplex) + '\n'
    const += ''.join([' + '.join([str(N[i][j]) + ' V' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
    if Forbidden:
        const += ''.join(['V' + str(i) + ' = 0\n' for i in Forbidden])
    if len(limits) == 1:
        const += ''.join(['V' + str(i) + ' <= 1\n' for i in allowed])
    elif len(limits) == len(allowed):
        const += ''.join(['V' + str(allowed[i]) + ' <= ' + str(allowed[i]) + '\n' for i in range(len(allowed))])
    else:
        print('Error: incompatible dimension of limits!')
    bounds = 'Bounds\n'
    if Negative:
        bounds += ''.join(['-inf <= V' + str(i) + ' <= 0\n' for i in Negative])
    bounds += ''.join(['V' + str(i) + ' free\n' for i in Rev if i not in Negative and [_f for _f in [N[k][i] for k in range(m)] if _f]])
    if Cplex:
        const = const.replace('+ -', '-')
    f = open(Filename, 'w')
    if not Cplex:
        f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    return processFile(Filename)
