"""
EvoDock Version 0.1

created and developed by Maximilian Edich at Universitaet Bielefeld.
"""

import random
import TestScore
import argparse
import matplotlib.pyplot as plt

"""
Alanin	        Ala	A
Arginin	        Arg	R
Asparagin	    Asn	N
Asparaginsaeure	Asp	D
Cystein     	Cys	C
Glutamin	    Gln	Q
Glutaminsaeure	Glu	E
Glycin	        Gly	G
Histidin	    His	H
Isoleucin	    Ile	I
Leucin	        Leu	L
Lysin	        Lys	K
Methionin	    Met	M
Phenylalanin	Phe	F
Prolin	        Pro	P
Serin	        Ser	S
Threonin	    Thr	T
Tryptophan	    Trp	W
Tyrosin	        Tyr	Y
Valin	        Val	V
"""

# strings for commands
MUTATE = "mutate"
MUTATE_PLUS = "mutate+"
CROSSOVER = "crossover"
CROSSOVER_CLASSIC = "crossoverClassic"
SELECT = "select"
LOOP = "loop"

# indices for lists
ITERATION_COUNT_MUTATION_INDEX = 0
ITERATION_COUNT_CROSSOVER_INDEX = 1

# global parameters
global MAX_REC_DEPTH_GET_RANDOM_INDIVIDUAL
MAX_REC_DEPTH_GET_RANDOM_INDIVIDUAL = 100
global MAX_REC_DEPTH_GET_NEW_RANDOM_MUTANT
MAX_REC_DEPTH_GET_NEW_RANDOM_MUTANT = 100


# ## read out input arguments
parser = argparse.ArgumentParser(description="EvoDock - An Evolutionary Algorithm for Protein optimization")
parser.add_argument("-s", "--settings", type=str, required=True,
                    help="Path to the initial settings file.\nThe initial settings file must be a text file"
                         " containing only one single command per line.\nThe order of the commands is not important. "
                         "Spaces should only be used to seperate a command and its parameter. These are the specified"
                         " commands:\n'startpopulation n', where n is an int >= 1, specifying the number of random"
                         "generated individuals + the initial individual given by the initial protein.")
parser.add_argument("-r", "--routine", type=str, required=True, help="Path to the routine file.")
args = parser.parse_args()
# TODO check if path leads to existing file!
# specify routine file
# specify optional technical settings


# ## Specification of mutable amino acids, possible mutations and target score
# TODO: check if number inputs are real numbers
startPopulationSize = 0
targetScore = 0
original = [[], 0]
numberOfMutableAA = 0
allowedMutations = []
definedTargetScore = False

# try to load initial settings file
initialContent = ""
try:
    initialFile = open(args.settings, 'r')
    initialContent = initialFile.readlines()
    initialFile.close()
except FileNotFoundError:
    exit("Error: Initial settings file does not exist!")
except PermissionError:
    exit("Error: Initial settings file can not be opened, permission denied!")


lineIndex = 0
for line in initialContent:
    lineIndex += 1
    # prepare line
    line = line.strip()
    splitText = line.split(' ')
    if splitText[0] == "startpopulation":
        # set value of initial population size, only int > 0
        try:
            startPopulationSize = int(splitText[1])
        except ValueError:
            exit("Error, startpopulation value in initial settings file is not a valid number.")
        if startPopulationSize < 1:
            exit("Error, illegal value in line " + str(lineIndex) + " of the initial settings file."
                                                                    "Must be >= 1.")
    elif splitText[0] == "targetscore":
        definedTargetScore = True
        # set the target score, any value possible
        targetScore = int(splitText[1])
    elif splitText[0][0:len("initial>")] == "initial>":
        # original individual
        initialGene = splitText[0][len("initial>"):].split(',')
        original[0] = initialGene
        numberOfMutableAA = len(initialGene)
    elif splitText[0][0] == ">":
        # read allowed substitutions for mutations per position
        if splitText[0][1:] == "":
            # empty, allow all mutations
            allowedMutations.append(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                                     'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'O'])
        else:
            # not empty, allow only specified mutations
            allowed = splitText[0][1:].split(',')
            allowedMutations.append(allowed)
    else:
        exit("Error, undefined keywords in line " + str(lineIndex) + " of the initial settings file.")


# catch undefined values
if not original[0]:
    exit("Error in initial settings file: You have to specify the mutable amino acids of the original protein!")
if numberOfMutableAA != len(allowedMutations):
    exit("Error in initial settings file: Number of substitutions per position does not match with the number"
         "of mutable amino acids in the original protein!")
if startPopulationSize == 0:
    exit("Error in initial settings file: You have to specify the \"startpopulation n\" with n > 0!")
if not definedTargetScore:
    exit("Error in initial settings file: You have to specify the \"targetscore n\"!")


# TODO load optional history as look up table for score to avoid long re-calculations
lookUpScores = []


# ## Score Function
def get_and_write_score(targetIndividual):
    """
    Calculates the score of an individual and writes it into its score value.
    :param targetIndividual: An individual in the form of [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.
    :return: the calculated score as the individuals fitness, that was written into the score value of the individual.
    """
    # calculate the score (by using external software)
    score = TestScore.get_score(targetIndividual)
    # write the score into the individual and return score
    targetIndividual[1] = score
    return score


# ## Create Population - Functions
def get_random_genes():
    """
    Get a random combination of possible genes (in terms of genetic algorithms).
    :return: A list, containing random amino acids chosen by the specified pattern.
    """

    genes = []
    for aa in range(numberOfMutableAA):
        # for each amino acid, add a random amino acid from specified lists
        genes.append(random.choice(allowedMutations[aa]))

    return genes


def get_random_individual(history, recDepth):
    """
    Generate a new individual with random gene combination.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history within
    this function.
    :param recDepth: The recursion depth. If a generated individual is not new, the
    functions is recalled with increased recDepth as long as recDepth is smaller
    than the maximal recursion depth.
    :return: A new and random individual in the form [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid. If the maximal recursion depth is reached, before a new individual was found,
    None is returned.
    """
    if recDepth > MAX_REC_DEPTH_GET_RANDOM_INDIVIDUAL:
        # return None, if maximum recursion depth is reached
        return

    newIndividual = [get_random_genes(), 0]
    # check if the new individual appears in history
    for indiv in history:
        if indiv[0] == newIndividual[0]:
            # individual not new, repeat with new random individual
            return get_random_individual(history, recDepth + 1)
    # individual is new, add to history
    history.append(newIndividual)

    return newIndividual


def generate_initial_population(numberOfInitialIndividuals, history, originalIndividual):
    """
    Generate the initial population with random and new individuals.
    :param numberOfInitialIndividuals: The number ( > 0) of total individuals in the generated population.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param originalIndividual: An individual with the initial amino acids of interest. The original
    individual is added as first one to the new population.
    :return: A new population. A population is a list of individuals.
    """
    # initialize population
    newPopulation = [originalIndividual]
    print("Create initial Population of size " + str(numberOfInitialIndividuals) + " ...")
    for n in range(numberOfInitialIndividuals - 1):
        # generate new individual
        newIndividual = get_random_individual(history, 0)
        if newIndividual is None:
            break
        # add to population, if new individual is available
        newPopulation.append(newIndividual)

    # calculate score of initial population
    for indiv in newPopulation:
        get_and_write_score(indiv)
        print(indiv)

    return newPopulation


# ## Mutation and Crossing Over Functions
def get_new_random_mutant(parent, history, recDepth, statsDataList):
    """
    Creates a copy of the given individual (parent) and chooses by random choice a mutable position.
    For this position a random character of the allowed characters on this position is chosen.
    :param parent: The parents gene combination is copied, before a mutation is applied.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history within
    this function.
    :param recDepth: The recursion depth. If a generated individual is not new, the
    functions is recalled with increased recDepth as long as recDepth is smaller
    than the maximal recursion depth.
    :param statsDataList: Used to gather data for several stats.
    :return: An individual, if the generated individual is new.
    Otherwise a new recursion is performed, until the max recursion depth is reached, then None is
    returned. The returned mutant follows the format [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.
    """
    if recDepth > MAX_REC_DEPTH_GET_NEW_RANDOM_MUTANT:
        # return None, if maximum recursion depth is reached
        return

    # create copy of individual
    newMutant = [parent[0].copy(), 0]
    # choose random position and mutate it until it is not equal to its parent
    while newMutant[0] == parent[0]:
        r = random.choice(range(numberOfMutableAA))
        newMutant[0][r] = random.choice(allowedMutations[r])
    # check the history
    for indiv in history:
        if indiv[0] == newMutant[0]:
            # repeat with new random mutation
            return get_new_random_mutant(parent, history, recDepth + 1, statsDataList)

    history.append(newMutant)
    statsDataList[ITERATION_COUNT_MUTATION_INDEX] += 1
    return newMutant


def get_random_mutants(parent, numberOfNewMutants, history, statsDataList):
    """
    Get from a given parent individual several new mutants.
    :param parent: The parents gene combination is copied, before mutations are applied.
    :param numberOfNewMutants: The maximum number of new mutants generated by this function.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param statsDataList: Used to gather data for several stats.
    :return: A list of new individuals. Since the mutation function can fail, this list can be
    empty. If all mutations were successful, the list contains 'numberOfNewMutants' items.
    """
    mutants = []
    # generate several mutants
    for mutant in range(numberOfNewMutants):
        # generate new individual, if available. Start at recursion depth 0.
        newMutant = get_new_random_mutant(parent, history, 0, statsDataList)
        if newMutant is None:
            break
        # only keep new individuals
        mutants.append(newMutant)

    return mutants


def mutate_population(inputPopulation, numberOfNewMutants, history, statsDataList):
    """
    Iterate through the whole population and create several mutants for each. All the original individuals from
    the inputPopulation and all new are transferred into the returned new population.
    :param inputPopulation: The population from which each individual is passed trough the mutation process.
    :param numberOfNewMutants: The maximum number of new mutants generated by this function.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param statsDataList: Used to gather data for several stats.
    :return: A new population containing the whole inputPopulation and all new mutants.
    """
    print("\ncreate Mutants from individuals")
    newPopulation = []
    # for each individual in population
    for parent in inputPopulation:
        # keep the individual in population
        newPopulation.append(parent)
        # generate new mutants
        mutants = get_random_mutants(parent, numberOfNewMutants, history, statsDataList)
        # add new mutants to population
        for mutant in mutants:
            newPopulation.append(mutant)

    return newPopulation


def mutate_and_keep_improvements(inputPopulation, numberOfNewMutants, history, statsDataList):
    """
    Iterate through the whole population and create several mutants for each. Only the original individuals from
    the inputPopulation and mutants with an improved score are transferred into the returned new population.
    :param inputPopulation: The population from which each individual is passed trough the mutation process.
    :param numberOfNewMutants: The maximum number of new mutants generated by this function.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param statsDataList: Used to gather data for several stats.
    :return: A new population containing the whole inputPopulation and improved mutants.
    """
    print("\ncreate Mutants from individuals and only keep improvements")
    newPopulation = []
    # for each individual in population
    for parent in inputPopulation:
        # keep the individual in population
        newPopulation.append(parent)
        # generate new mutants
        mutants = get_random_mutants(parent, numberOfNewMutants, history, statsDataList)
        # check, if mutants are better than parent
        for mutant in mutants:
            get_and_write_score(mutant)
            if get_individuals_score_relative_to_targetScore(mutant) < get_individuals_score_relative_to_targetScore(parent):
                # if the mutant is an improvement, keep it in new population
                newPopulation.append(mutant)

    return newPopulation


def get_random_mating_partner(inputPopulation, matingPartnerOne):
    """
    Get a random mating partner for cross over operations.
    :param inputPopulation: The population to pick from.
    :param matingPartnerOne: The individual you search a mating partner for.
    :return: If the population contains only one individual, the given matingPartnerOne is returned.
    Otherwise, a random picked other individual is returned.
    """
    # prevent infinite recursion
    if len(inputPopulation) <= 1:
        return matingPartnerOne

    matingPartnerTwo = random.choice(inputPopulation)
    if matingPartnerOne == matingPartnerTwo:
        return get_random_mating_partner(inputPopulation, matingPartnerOne)

    return matingPartnerTwo


def perform_crossing_over_classic(inputPopulation, repetitions, history, statsDataList):
    """
    Perform a classic cross over in terms of Genetic Algorithms with the whole population. All recombination
    are kept.
    :param inputPopulation: The population from which each individual is passed trough the crossing over process.
    :param repetitions: Number, how many mating partners are chosen for each individual in the population.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history within
    this function.
    :param statsDataList: Used to gather data for several stats.
    :return: A new population containing the whole inputPopulation and all new individuals from recombination, if it
    did not appear in the history.
    """
    newPopulation = []
    # add whole input population
    newPopulation.extend(inputPopulation)
    index = 0
    for matingPartnerOne in inputPopulation:
        index += 1
        # perform crossover with each following individual
        for matingPartnerNumber in range(repetitions):
            matingPartnerTwo = get_random_mating_partner(inputPopulation, matingPartnerOne)
            # set random crossover point, included in tail part
            crossPoint = random.choice(range(numberOfMutableAA - 1)) + 1
            # create new gene containing parts of the parents genes
            newGene = matingPartnerOne[0][0:crossPoint].copy()
            newGene.extend(matingPartnerTwo[0][crossPoint:].copy())
            # check history
            isNew = True
            for historyEntry in history:
                if historyEntry[0] == newGene:
                    isNew = False
                    break
            # if new, create individual and calculate score
            if isNew:
                newIndividual = [newGene, 0]
                get_and_write_score(newIndividual)
                # add to history and population
                history.append(newIndividual)
                statsDataList[ITERATION_COUNT_CROSSOVER_INDEX] += 1
                newPopulation.append(newIndividual)

    return newPopulation


def get_random_bit_mask(length):
    """
    Generate a random bit mask with the given length. The mask is used to perform a uniform crossing over.
    :param length: The length must match the number of mutable amino acids, if used for a uniform crossing over.
    :return: A list with the specified length, containing a zero or a one randomly chosen for each position.
    """
    mask = []
    for n in range(length):
        mask.append(random.getrandbits(1))
    # check if mask is not only 1 or only 0
    justZero = True
    justOne = True
    for n in range(length):
        if mask[n] == 1:
            justZero = False
        else:
            justOne = False
    if justZero or justOne:
        return get_random_bit_mask(length)

    return mask


def perform_crossing_over_uniform(inputPopulation, repetitions, history, statsDataList):
    """
    Perform a uniform cross over in terms of Genetic Algorithms with the whole population. All recombination
    are kept. For each position, a coin flip chooses the parent.
    :param inputPopulation: The population from which each individual is passed trough the crossing over process.
    :param repetitions: Number, how many mating partners are chosen for each individual in the population.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history within
    this function.
    :param statsDataList: Used to gather data for several stats.
    :return: A new population containing the whole inputPopulation and all new individuals from recombination, if it
    did not appear in the history.
    """
    newPopulation = []
    # add whole input population
    newPopulation.extend(inputPopulation)
    for matingPartnerOne in inputPopulation:
        # perform crossover with 'repetition' random individuals
        for matingPartnerNumber in range(repetitions):
            matingPartnerTwo = get_random_mating_partner(inputPopulation, matingPartnerOne)
            # create mask for choosing parent for each gene
            mask = get_random_bit_mask(numberOfMutableAA)
            newGene = []
            for n in range(numberOfMutableAA):
                if mask[n]:
                    newGene.append(matingPartnerOne[0][n])
                else:
                    newGene.append(matingPartnerTwo[0][n])
            print(matingPartnerOne[0])
            print(matingPartnerTwo[0])
            print(mask)
            print(newGene)
            print("\n")
            # check history
            isNew = True
            for historyEntry in history:
                if historyEntry[0] == newGene:
                    isNew = False
                    break
            # if new, create individual and calculate score
            if isNew:
                newIndividual = [newGene, 0]
                get_and_write_score(newIndividual)
                # add to history and population
                history.append(newIndividual)
                statsDataList[ITERATION_COUNT_CROSSOVER_INDEX] += 1
                newPopulation.append(newIndividual)

    return newPopulation


# ## Selection Functions
def get_individuals_score_relative_to_targetScore(inputIndividual):
    """
    Get the score value of an individual relative to the target score.
    :param inputIndividual: The individual of interest.
    :return: The difference between the individuals score value and the target score.
    If the score value was not previously calculated, it will return the default score.
    A difference of 0 is interpreted as best result.
    """
    return abs(targetScore - inputIndividual[1])


def ger_average_score(inputPopulation):
    """
    Get the average absolute score of the population.
    :param inputPopulation: The population of interest.
    :return: A float number, representing the average score without consideration of the target score.
    """
    score = 0
    for indiv in inputPopulation:
        score += indiv[1]

    return score / len(inputPopulation)


def select_fittest_by_number(selectionNumber, numberOfRandomPicks, inputPopulation):
    """
    Select a given number of individuals from the population, after it was sorted by fitness.
    :param selectionNumber: The number of individuals you want to select.
    :param numberOfRandomPicks: The number of randomly picked individuals, which do not fall into the selected ones
    with the best score.
    :param inputPopulation: The population you want to select from.
    :return: A list of the fittest individuals from the input population. If the selectionNumber
    exceeds the size of the population, the whole population is returned.
    In any case the list of individuals is sorted by fitness, the fittest individuals being at the beginning
    of the list.
    """
    # check if selection number exceeds population size
    if selectionNumber > len(inputPopulation):
        selectionNumber = len(inputPopulation)

    # check if number of random picked weak individuals exceeds the selection number
    if numberOfRandomPicks > selectionNumber:
        numberOfRandomPicks = selectionNumber

    keepPopulation = []
    # sort by score
    inputPopulation.sort(reverse=False, key=get_individuals_score_relative_to_targetScore)
    # select first ones, representing the fittest individuals of the current population
    for index in range(int(selectionNumber - numberOfRandomPicks)):
        selectedIndividual = inputPopulation.pop(0)
        print(selectedIndividual)
        keepPopulation.append(selectedIndividual)
    # select random individuals not included in keepPopulation yet
    for index in range(int(numberOfRandomPicks)):
        selectedIndividual = inputPopulation.pop(random.choice(range(len(inputPopulation) - 1)))
        print(selectedIndividual)
        keepPopulation.append(selectedIndividual)

    return keepPopulation


def select_fittest_by_fraction(fractionPercent, randomPicksPercent, inputPopulation):
    """
    Select a given proportion of individuals from the inputPopulation, after it was sorted by fitness.
    :param fractionPercent: A fraction in percent, written as float value in range [0, 1].
    Determines number of individuals you want to select from the inputPopulation.
    :param randomPicksPercent: A fraction in percent, written as float value in range [0, 1].
    Determines fraction of the selected individuals, that are randomly picked and not by the best score.
    :param inputPopulation: The population you want to select from.
    :return: A list of the fittest individuals from the input population.
    The list of individuals is sorted by fitness, the fittest individuals being at the beginning
    of the list. If the input value is out of range,
    the inputPopulation is returned.
    """
    # return inputPopulation, if input values are not valid
    if fractionPercent < 0 or fractionPercent > 1:
        print("Warning: 'fractionPercent' value in 'SelectFittestByFraction' is out of range!")
        return inputPopulation
    if randomPicksPercent < 0 or randomPicksPercent > 1:
        print("Warning: 'randomPicksPercent' value in 'SelectFittestByFraction' is out of range!")
        return inputPopulation

    # calculate number of individuals being selected
    keepInt = int(fractionPercent * len(inputPopulation))
    randomPicks = int(randomPicksPercent * keepInt)

    return select_fittest_by_number(keepInt, randomPicks, inputPopulation)


# ## Perform Evolution
def check_routine():
    """
    Iterates through the specified routine file and checks, if all commands
    are valid. The commands are not executed.
    :return: A list of errors, represented as strings.
    """
    routineFile = open("routine.txt", 'r')
    routines = routineFile.readlines()
    routineFile.close()

    routineErrors = []
    index = 0
    for routine in routines:
        error = None
        index += 1
        # strip line and separate it by spaces
        splitCommand = routine.strip().split(' ')
        undefined = False
        if len(splitCommand) == 1:
            if not routine == "\n":
                undefined = True
        else:
            if splitCommand[0] == MUTATE:
                mutationNumber = int(splitCommand[1])
                if mutationNumber < 0:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        splitCommand[1])
            elif splitCommand[0] == MUTATE_PLUS:
                mutationNumber = int(splitCommand[1])
                if mutationNumber < 0:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        splitCommand[1])
            elif splitCommand[0] == CROSSOVER:
                repetitionNumber = int(splitCommand[1])
                if repetitionNumber < 0:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        splitCommand[1])
            elif splitCommand[0] == CROSSOVER_CLASSIC:
                repetitionNumber = int(splitCommand[1])
                if repetitionNumber < 0:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        splitCommand[1])
            elif splitCommand[0] == SELECT:
                selectionParam = float(splitCommand[1])
                if selectionParam < 0:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        splitCommand[1])
                if len(splitCommand) > 2:
                    selectionParam = float(splitCommand[2])
                    if selectionParam < 0:
                        error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                            splitCommand[1])

            elif splitCommand[0] == LOOP:
                repeatNumber = int(splitCommand[1])
                if repeatNumber < 1:
                    error = "Error in routine file, illegal value in line " \
                            + str(index) + ": " + str(splitCommand[1] + ". Must be an int, >= 1")
            else:
                undefined = True
        if undefined:
            error = "Error in routine file, undefined keyword in line " + str(index) + ": " + str(splitCommand[0])

        if error is not None:
            print(error)
            routineErrors.append(error)

    return routineErrors


def perform_routine(inputPopulation, history, statsDataList):
    """
    Iterates through the specified routine file and executes each command in the given order.
    :param inputPopulation: The population on which the evolution will be performed.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param statsDataList: Used to gather data for several stats.
    :return: A list with the final population after the execution of the whole routine, and some statistical values.
    First item: The final population.
    Second item: The best score over time.
    Third item: The average score over time.
    """
    # load file content
    routineFile = open("routine.txt", 'r')
    routines = routineFile.readlines()
    routineFile.close()

    # init stat variables
    bestScoresOverTime = [inputPopulation[0][1]]
    averageScoresOverTime = [inputPopulation[0][1]]
    inputPopulation.sort(reverse=False, key=get_individuals_score_relative_to_targetScore)
    bestScoresOverTime.append(inputPopulation[0][1])
    averageScoresOverTime.append(ger_average_score(inputPopulation))

    # run routine
    routineStep = 0
    loopNumber = 0
    loopJump = 0
    while routineStep < len(routines):
        routine = routines[routineStep]
        print("Population Size: " + str(len(inputPopulation)))
        # strip line and separate it by spaces
        splitCommand = routine.strip().split(' ')
        if len(splitCommand) >= 2:
            if splitCommand[0] == MUTATE:
                # perform mutation
                mutationNumber = int(splitCommand[1])
                inputPopulation = mutate_population(inputPopulation, mutationNumber, history, statsDataList)
            elif splitCommand[0] == MUTATE_PLUS:
                # perform mutation
                mutationNumber = int(splitCommand[1])
                inputPopulation = mutate_and_keep_improvements(inputPopulation, mutationNumber, history, statsDataList)
            elif splitCommand[0] == CROSSOVER:
                # perform uniform crossing over
                repetitionNumber = int(splitCommand[1])
                inputPopulation = perform_crossing_over_uniform(inputPopulation, repetitionNumber, history, statsDataList)
            elif splitCommand[0] == CROSSOVER_CLASSIC:
                # perform uniform crossing over
                repetitionNumber = int(splitCommand[1])
                inputPopulation = perform_crossing_over_classic(inputPopulation, repetitionNumber, history, statsDataList)
            elif splitCommand[0] == SELECT:
                # select number or of fraction of mutants
                selectionParam1 = float(splitCommand[1])
                if len(splitCommand) > 2:
                    selectionParam2 = float(splitCommand[2])
                else:
                    selectionParam2 = 0

                # check selection method
                if selectionParam1 > 1:
                    # select by number
                    inputPopulation = select_fittest_by_number(int(selectionParam1), int(selectionParam2), inputPopulation)
                elif selectionParam1 >= 0:
                    # select by fraction
                    inputPopulation = select_fittest_by_fraction(selectionParam1, selectionParam2, inputPopulation)
                bestScoresOverTime.append(inputPopulation[0][1])
                averageScoresOverTime.append(ger_average_score(inputPopulation))
            elif splitCommand[0] == LOOP:
                # set repeat
                loopNumber = int(splitCommand[1]) - 1
                loopJump = routineStep + 1
        elif routine == "\n":
            if loopNumber > 0:
                loopNumber -= 1
                routineStep = loopJump - 1
        routineStep += 1
        if routineStep >= len(routines):
            if loopNumber > 0:
                loopNumber -= 1
                routineStep = loopJump

    print("Population Size: " + str(len(inputPopulation)))
    return [inputPopulation, bestScoresOverTime, averageScoresOverTime]


# ## Perform Evolution
# check errors in routine file
errors = check_routine()
if len(errors) > 0:
    exit("Exit algorithm due to errors in routine file.")
# no errors, proceed with evolution


# init variables for evolution
initialScore = get_and_write_score(original)
print("Initial amino acids: " + str(original[0]))
print("Initial Score: " + str(initialScore) + "\n")
totalHistory = []
iterationCounts = [0, 0]

# generate random population
population = generate_initial_population(startPopulationSize, totalHistory, original)
# perform the whole evolution routine
routineResults = perform_routine(population, totalHistory, iterationCounts)
population = routineResults[0]
bestScores = routineResults[1]
averageScores = routineResults[2]

# analyze results
# count max target
count = 0
bestScore = population[0][1]
print("Best score: " + str(bestScore) + " | with difference to target score: "
      + str(get_individuals_score_relative_to_targetScore(population[0])))
for individual in population:
    if individual[1] == bestScore:
        count += 1
print("Number of mutants sharing best score: " + str(count))
print("Number of accepted mutations: " + str(iterationCounts[ITERATION_COUNT_MUTATION_INDEX]))
print("Number of accepted crossing overs: " + str(iterationCounts[ITERATION_COUNT_CROSSOVER_INDEX]))

# ## Save Output
# output results
resultsFileContent = "Population Size: " + str(len(population)) + "\n" \
                     + str("Best score: " + str(bestScore) + " | with difference to target score: "
                           + str(get_individuals_score_relative_to_targetScore(population[0]))) + "\n"
resultsFileContent += str("Number of mutants sharing best score: " + str(count)) + "\n"
resultsFileContent += "Number of accepted mutations: " + str(iterationCounts[ITERATION_COUNT_MUTATION_INDEX]) + "\n"
resultsFileContent += "Number of accepted crossing overs: " + str(iterationCounts[ITERATION_COUNT_CROSSOVER_INDEX])\
                      + "\n"
resultsFileContent += "Best Score over Time: " + str(bestScores) + "\n"
resultsFileContent += "Average Score over Time: " + str(averageScores) + "\n"
for individual in population:
    resultsFileContent += str(individual[0]) + ", " + str(individual[1]) + "\n"

resultsFile = open("EvoDock_results.txt", 'w')
resultsFile.write(resultsFileContent)
resultsFile.close()

# output history
totalHistory.sort(reverse=False, key=get_individuals_score_relative_to_targetScore)
historyFileContent = "Entries in history: " + str(len(totalHistory)) + "\n"
for individual in totalHistory:
    historyFileContent += str(individual[0]) + ", " + str(individual[1]) + "\n"

historyFile = open("EvoDock_history.txt", 'w')
historyFile.write(historyFileContent)
historyFile.close()

# Plot results
plt.plot(bestScores)
plt.plot(bestScores, "ob")
plt.plot(averageScores)
plt.plot(averageScores, "or")
plt.show()
