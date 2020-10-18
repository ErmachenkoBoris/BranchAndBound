import multiprocessing
import urllib.request
import numpy as np
import time
import cplex

timeLimitValue = 1
delta = 0.00001
local = True

paths = [
      'http://iridia.ulb.ac.be/~fmascia/files/DIMACS/C125.9.clq',
    # 'http://iridia.ulb.ac.be/~fmascia/files/DIMACS/brock200_2.clq',
    #  'http://iridia.ulb.ac.be/~fmascia/files/DIMACS/keller4.clq'
]

localPath = [
    'graphs/c-fat200-1.clq',
    'graphs/c-fat200-2.clq',
    'graphs/c-fat200-5.clq',
    'graphs/c-fat500-1.clq',
    'graphs/c-fat500-10.clq',
    'graphs/c-fat500-2.clq',
    'graphs/c-fat500-5.clq',
    'graphs/MANN_a9.clq',
    'graphs/hamming6-2.clq',
    'graphs/hamming6-4.clq',
    'graphs/gen200_p0.9_44.clq',
    'graphs/gen200_p0.9_55.clq',
    'graphs/san200_0.7_1.clq',
    'graphs/san200_0.9_1.clq',
    'graphs/san200_0.9_2.clq',
    'graphs/san200_0.9_3.clq',
    'graphs/sanr200_0.7.clq',
    'graphs/C125.9.clq',
    'graphs/keller4.clq',
    'graphs/brock200_1.clq',
    'graphs/brock200_2.clq',
    'graphs/brock200_3.clq',
    'graphs/brock200_4.clq',
    'graphs/p_hat300-1.clq',
    'graphs/p_hat300-2.clq'
]

# -----------------------------------------------TEST-------------------------------
# matrixTest = np.zeros((4, 4))
# matrixTest = np.array([
#            [1,1,0,0],
#           [1,1,1,1],
#           [0,1,1,1],
#           [0,1,1,1]])
matrixTest = np.array([
            [1,1,1,1],
           [1,1,1,1],
          [1,1,1,1],
    [1,1,1,1]])
# -----------------------------------------------TEST-------------------------------



def openGraph(filePath):
    n = -1
    m = -1
    global local
    if local==True:
        file = open(filePath)
    else:
        file = urllib.request.urlopen(filePath)
    for line in file:
        if local == False:
            line = line.decode('ascii')
            line = line.strip('\n')
        if line.startswith('p'):
            n = int(line.split(' ')[2])
            m = int(line.split(' ')[3])
            break
    graphMatrix = np.zeros((n, n))
    for line in file:
        if local == False:
            line = line.decode('ascii')
            line = line.strip('\n')
        if line.startswith('e'):
            i = int(line.split(' ')[1]) - 1
            j = int(line.split(' ')[2]) - 1
            graphMatrix[i, j] = 1
            graphMatrix[j, i] = 1
    return n, m, graphMatrix

def colorGreedy(matrix, edges):
    V = [i for i in range(edges)]
    colorGroups = [[]]
    coloredV = [-1 for p in range(edges)]
    k = 0
    for i in range(edges):
        if i not in V:
            continue
        colorGroups[k].append(i)
        V.remove(i)
        while len(matrix[i].nonzero()[0]) != edges: # пока есть ненулевые
            for j in range(i, edges): # ?
                if matrix[i, j] == 0 and j in V:
                    break
            if j == edges:
                break
            if j == edges-1 and matrix[i, j] != 0 or j not in V:
                break

            colorGroups[k].append(j)
            V.remove(j)
            matrix[i] = matrix[i] + matrix[j]

        k = k+1
        colorGroups.insert(k, [])
    for i in range(k):
        for j in range(len(colorGroups[i])):
            coloredV[colorGroups[i][j]] = i
    return coloredV

def getWithMaxColorNumber(matrix, n, coloredEdges):
    maxColorCount = 0
    maxColorCountEdges = []
    for i in range(n):
        colorTmpCounter = []
        for j in range(n):
            if matrix[i,j] == 1 and i!=j and coloredEdges[j] not in colorTmpCounter:
                colorTmpCounter.append(coloredEdges[j])
        if len(colorTmpCounter) == maxColorCount:
            maxColorCountEdges.append(i)
        if len(colorTmpCounter) > maxColorCount:
            maxColorCount = len(colorTmpCounter)
            maxColorCountEdges = [i]
    return maxColorCountEdges

def findEvristicClique(maxColorCandidat, matrix, n, coloredEd):
    clickEvr = [maxColorCandidat]
    clickCandidat = []
    for i in range(n):
        if matrix[maxColorCandidat, i] == 1 and i != maxColorCandidat:
            clickCandidat.append(i)

    def findClickEvr(clickEvrF, clickCandidatF, matrix):
       maxColorLocalValue = -1
       maxColorLocal = clickCandidatF[0]

       # ищем максимального соседа
       for z in clickCandidatF:
           if coloredEd[z] >=maxColorLocalValue and z not in clickEvrF:
               maxColorLocalValue = coloredEd[z]
               maxColorLocal = z
       clickEvrF.append(maxColorLocal) # добавили в клику

       clickLocalCandidat = []  # ищем соседей новых
       for i in range(n):
           if matrix[maxColorLocal, i] == 1 and i != maxColorLocal:
               clickLocalCandidat.append(i)

        # находим пересечение со старыми соседями
       newCandidats = list(set(clickCandidatF) & set(clickLocalCandidat))
       if len(newCandidats) > 0:
           return findClickEvr(clickEvrF, newCandidats, matrix)
       return clickEvrF

    return findClickEvr(clickEvr, clickCandidat, matrix)



def evristic(path):

    print('---evristic')
    n, m, confusion_matrix = openGraph(path)

    # ----------------TEST
    # n = 4
    # confusion_matrix = matrixTest
    # ----------------TEST

    start_evr_time = time.time()

    coloredEd = colorGreedy(confusion_matrix.copy(), n)

    maxColor = getWithMaxColorNumber(confusion_matrix.copy(), n, coloredEd)

    bestEvrValue = -1
    bestEvr = []

    for i in range(min(100, len(maxColor))):
        clickEvristic = findEvristicClique(maxColor[i], confusion_matrix.copy(), n, coloredEd)
        if len(clickEvristic) > bestEvrValue:
            bestEvrValue = len(clickEvristic)
            bestEvr = clickEvristic.copy()

    clickValue = [0 for i in range(n)]
    for i in bestEvr:
        clickValue[i] = 1

    print('')
    print('')
    print('--------------- grath N ', path)
    print("--- %s seconds EVRISTIC---" % (time.time() - start_evr_time))
    print('coloredEd ', coloredEd)
    print('maxColor ', maxColor)
    print('bestEvr ', bestEvr)
    print('clickEvristicPower ', bestEvrValue)
    print(clickValue)
    print('--------------- grath N ', path)

    return bestEvrValue, clickValue, confusion_matrix, n

def initalClickCPLEX(constrains,
            constrainsNames,
            constrainsTypes,
            constrainsRightParts,
            maxCliqueModel,
            matrix,
            n):
    maxCliqueModel.variables.add(names=["y" + str(i) for i in range(n)],
                                   types=[maxCliqueModel.variables.type.continuous for i in range(n)])

    for i in range(n):
        maxCliqueModel.variables.set_lower_bounds(i, 0.0)

    for i in range(matrix.shape[0]):
        for j in range(i + 1, matrix.shape[1]):
            if matrix[i][j] == 0.0:
                constrains.append([["y" + str(i), "y" + str(j)], [1, 1]])
                constrainsNames.append("constraint_" + str(i) + "_" + str(j))
                constrainsTypes.append('L')
                constrainsRightParts.append(1.0)

    maxCliqueModel.linear_constraints.add(
        lin_expr=constrains,
        rhs=constrainsRightParts,
        names=constrainsNames,
        senses=constrainsTypes)

    maxCliqueModel.set_log_stream(None)
    maxCliqueModel.set_warning_stream(None)
    maxCliqueModel.set_results_stream(None)

    for i in range(n):
        maxCliqueModel.objective.set_linear("y" + str(i), 1)

    maxCliqueModel.objective.set_sense(maxCliqueModel.objective.sense.maximize)

    maxCliqueModel.solve()
    values = maxCliqueModel.solution.get_values()
    result = 0
    for v in values:
        result = result + v
    return values

def bnbContainer(evristicPower, evristicValues, matrix, n, graph, return_dict):
    start_BNB_time = time.time()
    constrains = []
    constrainsNames = []
    constrainsTypes = []
    constrainsRightParts = []

    global bestDecision

    maxCliqueModel = cplex.Cplex()

    initalClickCPLEX(
        constrains,
        constrainsNames,
        constrainsTypes,
        constrainsRightParts,
        maxCliqueModel,
        matrix,
        n, )

    bestDecision = evristicPower

    return_dict['power'] = evristicPower
    return_dict['values'] = evristicValues

    result, resultValues = BNB(evristicValues, maxCliqueModel, return_dict)
    print('grapth NAME ', graph)
    print("--- %s seconds BNB ---" % (time.time() - start_BNB_time))
    print('!!!!! result ', result)
    print('!!!!! resultValues ', resultValues)

def bnbStartEngine(graphs):
    for i in range(len(graphs)):
        evristicPower, evristicValues, matrix, n = evristic(graphs[i])
        if __name__ == '__main__':
            global timeLimitValue
            manager = multiprocessing.Manager()
            return_dict = manager.dict()
            return_dict['power'] = 0
            return_dict['values'] = []
            print('start bnb for ', graphs[i])
            p = multiprocessing.Process(target=bnbContainer, args=(evristicPower, evristicValues, matrix, n, graphs[i], return_dict))
            p.start()
            p.join(timeLimitValue)
            if p.is_alive():
                p.terminate()
                print("--- %s seconds LIMIT BNB ---" % timeLimitValue)
                print('!!!!! TIMEOUT result ', return_dict['power'] )
                print('!!!!! TIMEOUT resultValues ', return_dict['values'])

def numberWithDelta(number, maxValue, minVale, eps):
    if number + eps >= maxValue:
        return maxValue
    if number - eps <= minVale:
        return minVale
    return number

def addConstrain(i,value,maxCliqueModel):
    maxCliqueModel.linear_constraints.add(
        lin_expr=[[["y" + str(i)], [1]]],
        rhs=[value],
        names=["constraint_" + str(i)],
        senses=['E']
    )

def removeConstrain(i,maxCliqueModel):
    maxCliqueModel.linear_constraints.delete("constraint_" + str(i))

def solveWithCPLX(maxCliqueModel):
    maxCliqueModel.solve()
    return maxCliqueModel.solution.get_values()

def BNB(bestDecisionValue, maxCliqueModel, return_dict):
    global bestDecision
    currentDecisionValue = solveWithCPLX(maxCliqueModel)
    N = len(currentDecisionValue)

    currentDecision = 0

    for i in range(N):
        currentDecision = currentDecision + currentDecisionValue[i]

    global delta
    if currentDecision - delta <= bestDecision or currentDecision + delta <= bestDecision:
        return bestDecision, bestDecisionValue

    flag = False
    for index in range(N):
        currentDecisionValue[index] = numberWithDelta(currentDecisionValue[index], 1, 0, delta)
        if currentDecisionValue[index] !=0 and currentDecisionValue[index] !=1:
            flag = True
            break
    if index == N - 1 and flag == False:
        if currentDecision > bestDecision:

            print('new bestDecision ', bestDecision, currentDecision)

            bestDecision = currentDecision
            bestDecisionValue = currentDecisionValue

            return_dict['power'] = bestDecision
            return_dict['values'] = bestDecisionValue

        return bestDecision, bestDecisionValue


    addConstrain(index, 1, maxCliqueModel)

    BNB(bestDecisionValue, maxCliqueModel, return_dict)

    removeConstrain(index, maxCliqueModel)

    addConstrain(index, 0, maxCliqueModel)

    BNB(bestDecisionValue, maxCliqueModel, return_dict)

    removeConstrain(index, maxCliqueModel)

bestDecision = 0
# MAIN
if __name__ == '__main__':
    bnbStartEngine(localPath)
