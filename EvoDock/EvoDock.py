"""
EvoDock Version 0.1

created and developed by Maximilian Edich at Universitaet Bielefeld.
"""

# handle imports
try:
    import random
except ImportError as e:
    random = None
    exit("Error: Import of \"random\" failed. Make sure to provide this module, since it is essential. "
         "Error Message: " + str(e))
try:
    import argparse
except ImportError as e:
    argparse = None
    exit("Error: Import of \"argparse\" failed. Make sure to provide this module, since it is essential. "
         "Error Message: " + str(e))
try:
    import datetime
except ImportError as e:
    datetime = None
    exit("Error: Import of \"datetime\" failed. Make sure to provide this module, since it is essential. "
         "Error Message: " + str(e))
try:
    from pathlib import Path
except ImportError as e:
    Path = None
    exit("Error: Import of \"Path\" from \"pathlib\" failed. Make sure to provide this module, since it is essential. "
         "Error Message: " + str(e))
try:
    import MutateDockScoreModule
except ImportError as e:
    MutateDockScoreModule = None
    exit("Error: Import of \"MutateDockScoreModule\" failed. Make sure to provide this module, since it is essential. "
         "This module is part of EvoDock, so probably something went wrong with the installation. "
         "Make sure, that all EvoDock modules are in the same folder. Check the documentation for more information. "
         "Error Message: " + str(e))

PLOT = False
if PLOT:
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        plt = None
        exit("Error: Import of \"matplotlib.pyplot\" failed. Make sure to provide this module, since it is essential.")

print("\nIMPORTS DONE!")

# initial settings file values
START_POPULATION = "startpopulation"
TARGET_SCORE = "targetscore"
MODE = "mode"
MODE_LIGAND_BINDING = "LIGAND-BINDING"
MODE_PROTEIN_BINDING = "PROTEIN-BINDING"
MODE_THERMO_STABILITY = "THERMO-STABILITY"

# strings for commands
MUTATE = "mutate"
MUTATE_PLUS = "mutate+"
RECOMBINATION = "recombine"
RECOMBINATION_CLASSIC = "recombine_classic"
SELECT = "select"
LOOP = "loop"

# indices for lists
ITERATION_COUNT_MUTATION_INDEX = 0
ITERATION_COUNT_RECOMBINATION_INDEX = 1

# technical fixed parameters
MAX_REC_DEPTH_GET_RANDOM_INDIVIDUAL = 500
MAX_REC_DEPTH_GET_NEW_RANDOM_MUTANT = 500

PRINT_OUT = False

# ## read out input arguments
parser = argparse.ArgumentParser(description="EvoDock - An Evolutionary Algorithm for Protein Optimization")
parser.add_argument("-s", "--settings", type=str, required=True,
                    help="Path to the initial settings file.\nThe initial settings file must be a text file"
                         " containing only one single command per line.\nThe order of the commands is not important. "
                         "Spaces should only be used to seperate a command and its parameter. See the documentation"
                         " for the full list of commands.")
parser.add_argument("-r", "--routine", type=str, required=True, help="Path to the routine file. See the documentation"
                                                                     "for a list of all commands.")
parser.add_argument("-pre", "--prefix", type=str, required=False, help="Prefix of the output folder.")
args = parser.parse_args()
# TODO specify optional technical settings


# ## Specification of mutable amino acids, possible mutations and target score
# TODO: check if number inputs are real numbers
mds = MutateDockScoreModule.MutateDockScore()
start_population_size = 0
target_score = 0
original = [[], 0]
protein_name = ""
mode = ""
amino_acid_paths = []
number_of_mutable_aa = 0
allowed_mutations = []
defined_target_score = False

# load files
# # try to load routine settings file
routine_file_content = ""
try:
    routine_file = open(args.routine, 'r')
    routine_file_content = routine_file.readlines()
    routine_file.close()
except FileNotFoundError:
    exit("Error: Routine file does not exist!")
except PermissionError:
    exit("Error: Routine file can not be opened, permission denied!")

# # try to load initial settings file
initial_settings_file_content = ""
try:
    initial_file = open(args.settings, 'r')
    initial_settings_file_content = initial_file.readlines()
    initial_file.close()
except FileNotFoundError:
    exit("Error: Initial settings file does not exist!")
except PermissionError:
    exit("Error: Initial settings file can not be opened, permission denied!")

# # check content of initial settings file and load data
line_index = 0
for line in initial_settings_file_content:
    line_index += 1
    # prepare line
    line = line.strip()
    split_text = line.split(' ')
    if split_text[0] == START_POPULATION:
        # set value of initial population size, only int > 0
        try:
            start_population_size = int(split_text[1])
        except ValueError:
            exit("Error, " + START_POPULATION + " value in initial settings file is not a valid number, "
                 "must be int: \""
                 + split_text[1] + "\" in line " + str(line_index) + " of the initial settings file.")
        if start_population_size < 1:
            exit("Error, illegal value in line " + str(line_index) + " of the initial settings file. Must be >= 1.")

    elif split_text[0] == TARGET_SCORE:
        defined_target_score = True
        # set the target score, any value possible
        try:
            target_score = float(split_text[1])
        except ValueError:
            exit("Error, " + TARGET_SCORE + " value in initial settings file is not a valid number, must be float: \""
                 + split_text[1] + "\" in line " + str(line_index) + " of the initial settings file.")

    elif split_text[0][0:len("initial>")] == "initial>":
        # original individual
        initial_genes = split_text[0][len("initial>"):].split(',')
        original[0] = initial_genes
        number_of_mutable_aa = len(initial_genes)
    elif split_text[0][0] == ">":
        # read allowed substitutions for mutations per position
        if split_text[0][1:] == "":
            # empty, allow all mutations
            allowed_mutations.append(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                                      'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'O'])
        else:
            # not empty, allow only specified mutations
            allowed = split_text[0][1:].split(',')
            allowed_mutations.append(allowed)
    elif split_text[0] == "aa":
        # amino acid path
        amino_acid_paths.append(str(split_text[1]))
    elif split_text[0] == "prot":
        # protein pdb code
        protein_name = str(split_text[1])
    elif split_text[0] == MODE:
        # modification mode
        mode = str(split_text[1])
        if mode != MODE_LIGAND_BINDING:
            exit("Error in line " + str(line_index) + " of the initial settings file. Unknown mode. Use one of these: "
                 + MODE_LIGAND_BINDING)
    else:
        exit("Error, undefined keywords in line " + str(line_index) + " of the initial settings file: " + split_text[0])

# # catch undefined values
if not original[0]:
    exit("Error in initial settings file: You have to specify the mutable amino acids of the original protein!")
if number_of_mutable_aa != len(allowed_mutations):
    exit("Error in initial settings file: Number of substitutions per position does not match with the number"
         "of mutable amino acids in the original protein!")
if start_population_size == 0:
    exit("Error in initial settings file: You have to specify the \"" + START_POPULATION + " n\" with n > 0!")
if not defined_target_score:
    exit("Error in initial settings file: You have to specify the \"" + TARGET_SCORE + " f\"!")

# TODO read initial genes from pdb
# do this via new class

print("FILE LOADINGS DONE!\n")
print("original individual + amino acid paths:")
print(original[0])
print(amino_acid_paths)

# TODO load optional history as look up table for score to avoid long re-calculations
look_up_scores = []


# ## Score Function
def get_and_write_score(target_individual):
    """
    Calculates the score of an individual and writes it into its score value.
    :param target_individual: An individual in the form of [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.
    :return: The calculated score as the individuals fitness, that was written into the score value of the individual.
    """
    # calculate the score (by using external software)
    score = MutateDockScoreModule.get_score(target_individual, mds)
    # write the score into the individual and return score
    target_individual[1] = score
    return score


# ## Create Population - Functions
def get_random_genes():
    """
    Get a random combination of possible genes (in terms of genetic algorithms).
    :return: A list, containing random amino acids chosen by the specified pattern.
    """

    genes = []
    for aa in range(number_of_mutable_aa):
        # for each amino acid, add a random amino acid from specified lists
        genes.append(random.choice(allowed_mutations[aa]))

    return genes


def get_random_individual(history, rec_depth):
    """
    Generate a new individual with random gene combination.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history within
    this function.
    :param rec_depth: The recursion depth. If a generated individual is not new, the
    functions is recalled with increased recDepth as long as recDepth is smaller
    than the maximal recursion depth.
    :return: A new and random individual in the form [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid. If the maximal recursion depth is reached, before a new individual was found,
    None is returned.
    """
    if rec_depth > MAX_REC_DEPTH_GET_RANDOM_INDIVIDUAL:
        # return None, if maximum recursion depth is reached
        return

    new_individual = [get_random_genes(), 0]
    # check if the new individual appears in history
    for individual in history:
        if individual[0] == new_individual[0]:
            # individual not new, repeat with new random individual
            try:
                return get_random_individual(history, rec_depth + 1)
            except RecursionError as error:
                print("Please lower the new-random-individual-recursion-limit in the technical settings!"
                      "Error Message: " + str(error))
                return
    # individual is new, add to history
    history.append(new_individual)

    return new_individual


def generate_initial_population(number_of_initial_individuals, history, original_individual):
    """
    Generate the initial population with random and new individuals.
    :param number_of_initial_individuals: The number ( > 0) of total individuals in the generated population.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param original_individual: An individual with the initial amino acids of interest. The original
    individual is added as first one to the new population.
    :return: A new population. A population is a list of individuals.
    """
    # initialize population
    new_population = [original_individual]
    history.append(original_individual)
    if PRINT_OUT:
        print("Create initial Population of size " + str(number_of_initial_individuals) + " ...")
    for n in range(number_of_initial_individuals - 1):
        # generate new individual
        new_individual = get_random_individual(history, 0)
        if new_individual is None:
            break
        # add to population, if new individual is available
        new_population.append(new_individual)

    # calculate score of initial population
    for individual in new_population:
        get_and_write_score(individual)
        if PRINT_OUT:
            print(individual)

    return new_population


# ## Mutation and Recombination Functions
def get_new_random_mutant(parent, history, rec_depth, stats_data_list):
    """
    Creates a copy of the given individual (parent) and chooses by random choice a mutable position.
    For this position a random character of the allowed characters on this position is chosen.
    :param parent: The parents gene combination is copied, before a mutation is applied.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history within
    this function.
    :param rec_depth: The recursion depth. If a generated individual is not new, the
    functions is recalled with increased recDepth as long as recDepth is smaller
    than the maximal recursion depth.
    :param stats_data_list: Used to gather data for several stats.
    :return: An individual, if the generated individual is new.
    Otherwise a new recursion is performed, until the max recursion depth is reached, then None is
    returned. The returned mutant follows the format [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.
    """
    if rec_depth > MAX_REC_DEPTH_GET_NEW_RANDOM_MUTANT:
        # return None, if maximum recursion depth is reached
        return

    # create copy of individual
    new_mutant = [parent[0].copy(), 0]
    # choose random position and mutate it until it is not equal to its parent
    while new_mutant[0] == parent[0]:
        r = random.choice(range(number_of_mutable_aa))
        new_mutant[0][r] = random.choice(allowed_mutations[r])
    # check the history
    for indiv in history:
        if indiv[0] == new_mutant[0]:
            # repeat with new random mutation
            try:
                return get_new_random_mutant(parent, history, rec_depth + 1, stats_data_list)
            except RecursionError as error:
                print("Please lower the new-mutant-recursion-limit in the technical settings!"
                      "Error Message: " + str(error))
                return

    history.append(new_mutant)
    stats_data_list[ITERATION_COUNT_MUTATION_INDEX] += 1
    return new_mutant


def get_random_mutants(parent, number_of_new_mutants, history, stats_data_list):
    """
    Get from a given parent individual several new mutants.
    :param parent: The parents gene combination is copied, before mutations are applied.
    :param number_of_new_mutants: The maximum number of new mutants generated by this function.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param stats_data_list: Used to gather data for several stats.
    :return: A list of new individuals. Since the mutation function can fail, this list can be
    empty. If all mutations were successful, the list contains 'numberOfNewMutants' items.
    """
    mutants = []
    # generate several mutants
    for mutant in range(number_of_new_mutants):
        # generate new individual, if available. Start at recursion depth 0.
        new_mutant = get_new_random_mutant(parent, history, 0, stats_data_list)
        if new_mutant is None:
            break
        # only keep new individuals
        mutants.append(new_mutant)

    return mutants


def mutate_population(input_population, number_of_new_mutants, history, stats_data_list):
    """
    Iterate through the whole population and create several mutants for each. All the original individuals from
    the inputPopulation and all new are transferred into the returned new population.
    :param input_population: The population from which each individual is passed trough the mutation process.
    :param number_of_new_mutants: The maximum number of new mutants generated by this function.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param stats_data_list: Used to gather data for several stats.
    :return: A new population containing the whole inputPopulation and all new mutants.
    """
    if PRINT_OUT:
        print("\ncreate Mutants from individuals")
    new_population = []
    # for each individual in population
    for parent in input_population:
        # keep the individual in population
        new_population.append(parent)
        # generate new mutants
        mutants = get_random_mutants(parent, number_of_new_mutants, history, stats_data_list)
        # add new mutants to population
        for mutant in mutants:
            new_population.append(mutant)

    return new_population


def mutate_and_keep_improvements(input_population, number_of_new_mutants, history, stats_data_list):
    """
    Iterate through the whole population and create several mutants for each. Only the original individuals from
    the inputPopulation and mutants with an improved score are transferred into the returned new population.
    :param input_population: The population from which each individual is passed trough the mutation process.
    :param number_of_new_mutants: The maximum number of new mutants generated by this function.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param stats_data_list: Used to gather data for several stats.
    :return: A new population containing the whole inputPopulation and improved mutants.
    """
    if PRINT_OUT:
        print("\ncreate Mutants from individuals and only keep improvements")
    new_population = []
    # for each individual in population
    for parent in input_population:
        # keep the individual in population
        new_population.append(parent)
        # generate new mutants
        mutants = get_random_mutants(parent, number_of_new_mutants, history, stats_data_list)
        # check, if mutants are better than parent
        for mutant in mutants:
            get_and_write_score(mutant)
            if get_individuals_score_relative_to_targetScore(mutant) < \
                    get_individuals_score_relative_to_targetScore(parent):
                # if the mutant is an improvement, keep it in new population
                new_population.append(mutant)

    return new_population


def get_random_mating_partner(input_population, mating_partner_one):
    """
    Get a random mating partner for recombination operations.
    :param input_population: The population to pick from.
    :param mating_partner_one: The individual you search a mating partner for.
    :return: If the population contains only one individual, the given matingPartnerOne is returned.
    Otherwise, a random picked other individual is returned.
    """
    # prevent infinite recursion
    if len(input_population) <= 1:
        return mating_partner_one

    mating_partner_two = random.choice(input_population)
    if mating_partner_one == mating_partner_two:
        return get_random_mating_partner(input_population, mating_partner_one)

    return mating_partner_two


def perform_recombination_classic(input_population, repetitions, history, stats_data_list):
    """
    Perform a classic cross over in terms of Genetic Algorithms with the whole population. All recombination
    are kept.
    :param input_population: The population from which each individual is passed trough the recombination process.
    :param repetitions: Number, how many mating partners are chosen for each individual in the population.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history within
    this function.
    :param stats_data_list: Used to gather data for several stats.
    :return: A new population containing the whole inputPopulation and all new individuals from recombination, if it
    did not appear in the history.
    """
    new_population = []
    # add whole input population
    new_population.extend(input_population)
    index = 0
    for mating_partner_one in input_population:
        index += 1
        # perform recombination with each following individual
        for mating_partner_number in range(repetitions):
            mating_partner_two = get_random_mating_partner(input_population, mating_partner_one)
            # set random crossover point, included in tail part
            cross_point = random.choice(range(number_of_mutable_aa - 1)) + 1
            # create new genes containing parts of the parents genes
            new_genes = mating_partner_one[0][0:cross_point].copy()
            new_genes.extend(mating_partner_two[0][cross_point:].copy())
            # check history
            is_new = True
            for history_entry in history:
                if history_entry[0] == new_genes:
                    is_new = False
                    break
            # if new, create individual and calculate score
            if is_new:
                new_individual = [new_genes, 0]
                get_and_write_score(new_individual)
                # add to history and population
                history.append(new_individual)
                stats_data_list[ITERATION_COUNT_RECOMBINATION_INDEX] += 1
                new_population.append(new_individual)

    return new_population


def get_random_bit_mask(length):
    """
    Generate a random bit mask with the given length. The mask is used to perform a uniform recombination.
    :param length: The length must match the number of mutable amino acids, if used for a uniform recombination.
    :return: A list with the specified length, containing a zero or a one randomly chosen for each position.
    """
    mask = []
    for n in range(length):
        mask.append(random.getrandbits(1))
    # check if mask is not only 1 or only 0
    just_zeros = True
    just_ones = True
    for n in range(length):
        if mask[n] == 1:
            just_zeros = False
        else:
            just_ones = False
    if just_zeros or just_ones:
        return get_random_bit_mask(length)

    return mask


def perform_uniform_recombination(input_population, repetitions, history, stats_data_list):
    """
    Perform a uniform recombination in terms of Genetic Algorithms with the whole population. All recombination
    are kept. For each position, a coin flip chooses the parent.
    :param input_population: The population from which each individual is passed trough the recombination process.
    :param repetitions: Number, how many mating partners are chosen for each individual in the population.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history within
    this function.
    :param stats_data_list: Used to gather data for several stats.
    :return: A new population containing the whole inputPopulation and all new individuals from recombination, if it
    did not appear in the history.
    """
    new_population = []
    # add whole input population
    new_population.extend(input_population)
    for mating_partner_one in input_population:
        # perform recombination with 'repetition' random individuals
        for mating_partner_number in range(repetitions):
            mating_partner_two = get_random_mating_partner(input_population, mating_partner_one)
            # create mask for choosing parent for each gene
            mask = get_random_bit_mask(number_of_mutable_aa)
            new_genes = []
            for n in range(number_of_mutable_aa):
                if mask[n]:
                    new_genes.append(mating_partner_one[0][n])
                else:
                    new_genes.append(mating_partner_two[0][n])
            if PRINT_OUT:
                print(mating_partner_one[0])
                print(mating_partner_two[0])
                print(mask)
                print(new_genes)
                print("\n")
            # check history
            is_new = True
            for history_entry in history:
                if history_entry[0] == new_genes:
                    is_new = False
                    break
            # if new, create individual and calculate score
            if is_new:
                new_individual = [new_genes, 0]
                get_and_write_score(new_individual)
                # add to history and population
                history.append(new_individual)
                stats_data_list[ITERATION_COUNT_RECOMBINATION_INDEX] += 1
                new_population.append(new_individual)

    return new_population


# ## Selection Functions
def get_individuals_score_relative_to_targetScore(input_individual):
    """
    Get the score value of an individual relative to the target score.
    :param input_individual: The individual of interest.
    :return: The difference between the individuals score value and the target score.
    If the score value was not previously calculated, it will return the default score.
    A difference of 0 is interpreted as best result.
    """
    return abs(target_score - input_individual[1])


def ger_average_score(input_population):
    """
    Get the average absolute score of the population.
    :param input_population: The population of interest.
    :return: A float number, representing the average score without consideration of the target score.
    """
    score = 0
    for indiv in input_population:
        score += indiv[1]

    return score / len(input_population)


def select_fittest_by_number(selectionNumber, number_of_random_picks, input_population):
    """
    Select a given number of individuals from the population, after it was sorted by fitness.
    :param selectionNumber: The number of individuals you want to select.
    :param number_of_random_picks: The number of randomly picked individuals, which do not fall into the selected ones
    with the best score.
    :param input_population: The population you want to select from.
    :return: A list of the fittest individuals from the input population. If the selectionNumber
    exceeds the size of the population, the whole population is returned.
    In any case the list of individuals is sorted by fitness, the fittest individuals being at the beginning
    of the list.
    """
    # check if selection number exceeds population size
    if selectionNumber > len(input_population):
        selectionNumber = len(input_population)
    if selectionNumber < 1:
        selectionNumber = 1

    # check if number of random picked weak individuals exceeds the selection number
    if number_of_random_picks > selectionNumber:
        number_of_random_picks = selectionNumber

    keep_population = []
    # sort by score
    input_population.sort(reverse=False, key=get_individuals_score_relative_to_targetScore)
    # select first ones, representing the fittest individuals of the current population
    for index in range(int(selectionNumber - number_of_random_picks)):
        selected_individual = input_population.pop(0)
        if PRINT_OUT:
            print(selected_individual)
        keep_population.append(selected_individual)
    # select random individuals not included in keepPopulation yet
    for index in range(int(number_of_random_picks)):
        selected_individual = input_population.pop(random.choice(range(len(input_population) - 1)))
        if PRINT_OUT:
            print(selected_individual)
        keep_population.append(selected_individual)

    return keep_population


def select_fittest_by_fraction(fraction_percent, random_picks_percent, input_population):
    """
    Select a given proportion of individuals from the inputPopulation, after it was sorted by fitness.
    :param fraction_percent: A fraction in percent, written as float value in range [0, 1].
    Determines number of individuals you want to select from the inputPopulation.
    :param random_picks_percent: A fraction in percent, written as float value in range [0, 1].
    Determines fraction of the selected individuals, that are randomly picked and not by the best score.
    :param input_population: The population you want to select from.
    :return: A list of the fittest individuals from the input population.
    The list of individuals is sorted by fitness, the fittest individuals being at the beginning
    of the list. If the input value is out of range,
    the inputPopulation is returned.
    """
    # return inputPopulation, if input values are not valid
    if fraction_percent < 0 or fraction_percent > 1:
        if PRINT_OUT:
            print("Warning: 'fractionPercent' value in 'SelectFittestByFraction' is out of range!")
        return input_population
    if random_picks_percent < 0 or random_picks_percent > 1:
        if PRINT_OUT:
            print("Warning: 'randomPicksPercent' value in 'SelectFittestByFraction' is out of range!")
        return input_population

    # calculate number of individuals being selected
    keep_int = round(fraction_percent * len(input_population))
    random_picks = int(random_picks_percent * keep_int)
    print(str(keep_int) + " - " + str(random_picks))
    if keep_int < 1:
        keep_int = 1

    return select_fittest_by_number(keep_int, random_picks, input_population)


# ## Perform Evolution
def check_routine():
    """
    Iterates through the specified routine file and checks, if all commands
    are valid. The commands are not executed.
    :return: A list of errors, represented as strings.
    """
    # init error list and index
    routine_errors = []
    index = 0
    # for each line in file
    for routine in routine_file_content:
        error = None
        index += 1
        # strip line and separate it by spaces
        split_command = routine.strip().split(' ')
        undefined = False
        if len(split_command) == 1:
            if not routine == "\n":
                undefined = True
        else:
            if split_command[0] == MUTATE:
                try:
                    mutation_number = int(split_command[1])
                    if mutation_number < 0:
                        error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                            split_command[1]) + ". Must be >= 0"
                except ValueError:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        split_command[1]) + ". Must be an integer"
            elif split_command[0] == MUTATE_PLUS:
                try:
                    mutation_number = int(split_command[1])
                    if mutation_number < 0:
                        error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                            split_command[1]) + ". Must be >= 0"
                except ValueError:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        split_command[1]) + ". Must be an integer"
            elif split_command[0] == RECOMBINATION:
                try:
                    repetition_number = int(split_command[1])
                    if repetition_number < 0:
                        error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                            split_command[1]) + ". Must be >= 0"
                except ValueError:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        split_command[1]) + ". Must be an integer"
            elif split_command[0] == RECOMBINATION_CLASSIC:
                try:
                    repetition_number = int(split_command[1])
                    if repetition_number < 0:
                        error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                            split_command[1]) + ". Must be >= 0"
                except ValueError:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        split_command[1]) + ". Must be an integer"
            elif split_command[0] == SELECT:
                try:
                    selection_param = float(split_command[1])
                    if selection_param < 0:
                        error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                            split_command[1]) + ". Must be >= 0"
                    if len(split_command) > 2:
                        try:
                            selection_param = float(split_command[2])
                            if selection_param < 0:
                                error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                                    split_command[1]) + ". Must be >= 0"
                        except ValueError:
                            error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                                split_command[1]) + ". Must be float"
                except ValueError:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        split_command[1]) + ". Must be float or integer"

            elif split_command[0] == LOOP:
                try:
                    repeat_number = int(split_command[1])
                    if repeat_number < 1:
                        error = "Error in routine file, illegal value in line " \
                                + str(index) + ": " + str(split_command[1] + ". Must be an int, >= 1")
                except ValueError:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        split_command[1]) + ". Must bean integer"
            else:
                undefined = True
        if undefined:
            error = "Error in routine file, undefined keyword in line " + str(index) + ": " + str(split_command[0])

        if error is not None:
            print(error)
            routine_errors.append(error)

    return routine_errors


def save_output(input_population, best_scores_over_time, average_scores_over_time):
    """

    :return:
    """

    best_score = input_population[0][1]
    individual_count = 0
    for individual in population:
        if individual[1] == population[0][1]:
            individual_count += 1

    # output results
    results_file_content = "Final Population Size: " + str(len(input_population)) + "\n" \
                           + str("Best score: " + str(best_score) + " | with difference to target score: "
                                 + str(get_individuals_score_relative_to_targetScore(input_population[0]))) + "\n"
    results_file_content += str("Number of mutants sharing best score: " + str(individual_count)) + "\n"
    results_file_content += "Number of accepted mutations: " + str(
        iterationCounts[ITERATION_COUNT_MUTATION_INDEX]) + "\n"
    results_file_content += "Number of accepted recombination products: " \
                            + str(iterationCounts[ITERATION_COUNT_RECOMBINATION_INDEX]) + "\n"
    results_file_content += "Best Score over Time: " + str(best_scores_over_time) + "\n"
    results_file_content += "Average Score over Time: " + str(average_scores_over_time) + "\n"
    for individual in input_population:
        results_file_content += str(individual[0]) + ", " + str(individual[1]) + "\n"

    results_file = open(out_path + "/EvoDock_results.txt", 'w')
    results_file.write(results_file_content)
    results_file.close()

    # output history
    total_history.sort(reverse=False, key=get_individuals_score_relative_to_targetScore)
    history_file_content = "Entries in history: " + str(len(total_history)) + "\n"
    for individual in total_history:
        history_file_content += str(individual[0]) + ", " + str(individual[1]) + "\n"

    history_file = open(out_path + "/EvoDock_history.txt", 'w')
    history_file.write(history_file_content)
    history_file.close()

    # output run info TODO finish external module params integration
    # # basic run info
    current_time = datetime.datetime.now()
    run_info_file_content = "Run information of run " + run_string + "\n"
    run_info_file_content += "Run time: " + str(current_time - start_time) + "\n"
    # # Initial Settings File - Content
    run_info_file_content += "\n\n\nInitial Settings File - Content:\n"
    for file_line in initial_settings_file_content:
        run_info_file_content += file_line
    # # Routine File - Content
    run_info_file_content += "\n\n\nRoutine File - Content:\n"
    for file_line in routine_file_content:
        run_info_file_content += file_line
    # # External Settings
    run_info_file_content += "\n\n\nExternal Tool Settings File Contents:\n"

    run_info_file = open(out_path + "/EvoDock_run_information.txt", 'w')
    run_info_file.write(run_info_file_content)
    run_info_file.close()

    return


def perform_routine(input_population, history, stats_data_list):
    """
    Iterates through the specified routine file and executes each command in the given order.
    :param input_population: The population on which the evolution will be performed.
    :param history: A list containing all previous individuals. It is necessary to determine if
    a generated individual is new. The new individual will be added to the history during the process of this
    function.
    :param stats_data_list: Used to gather data for several stats.
    :return: A list with the final population after the execution of the whole routine, and some statistical values.
    First item: The final population.
    Second item: The best score over time.
    Third item: The average score over time.
    """

    # init stat variables
    best_scores_over_time = [input_population[0][1]]
    average_scores_over_time = [input_population[0][1]]
    input_population.sort(reverse=False, key=get_individuals_score_relative_to_targetScore)
    best_scores_over_time.append(input_population[0][1])
    average_scores_over_time.append(ger_average_score(input_population))

    # run routine
    routine_step = 0
    loop_number = 0
    loop_jump = 0
    while routine_step < len(routine_file_content):
        routine = routine_file_content[routine_step]
        if True:
            print("Population Size: " + str(len(input_population)))
            print("Do " + routine.strip() + "...")
        # strip line and separate it by spaces
        split_command = routine.strip().split(' ')
        if len(split_command) >= 2:
            if split_command[0] == MUTATE:
                # perform mutation
                mutation_number = int(split_command[1])
                input_population = mutate_population(input_population, mutation_number, history, stats_data_list)
            elif split_command[0] == MUTATE_PLUS:
                # perform mutation
                mutation_number = int(split_command[1])
                input_population = mutate_and_keep_improvements(input_population, mutation_number, history,
                                                                stats_data_list)
            elif split_command[0] == RECOMBINATION:
                # perform uniform recombination
                repetition_number = int(split_command[1])
                input_population = perform_uniform_recombination(input_population, repetition_number, history,
                                                                 stats_data_list)
            elif split_command[0] == RECOMBINATION_CLASSIC:
                # perform uniform recombination
                repetition_number = int(split_command[1])
                input_population = perform_recombination_classic(input_population, repetition_number, history,
                                                                 stats_data_list)
            elif split_command[0] == SELECT:
                # select number or of fraction of mutants
                selection_param1 = float(split_command[1])
                if len(split_command) > 2:
                    selection_param2 = float(split_command[2])
                else:
                    selection_param2 = 0

                # check selection method
                if selection_param1 > 1:
                    # select by number
                    input_population = select_fittest_by_number(int(selection_param1), int(selection_param2),
                                                                input_population)
                elif selection_param1 >= 0:
                    # select by fraction
                    input_population = select_fittest_by_fraction(selection_param1, selection_param2, input_population)

                best_scores_over_time.append(input_population[0][1])
                average_scores_over_time.append(ger_average_score(input_population))
                save_output(input_population, best_scores_over_time, average_scores_over_time)
            elif split_command[0] == LOOP:
                # set repeat
                loop_number = int(split_command[1]) - 1
                loop_jump = routine_step + 1
        elif routine == "\n":
            if loop_number > 0:
                loop_number -= 1
                routine_step = loop_jump - 1
        routine_step += 1
        if routine_step >= len(routine_file_content):
            if loop_number > 0:
                loop_number -= 1
                routine_step = loop_jump

    if PRINT_OUT:
        print("Population Size: " + str(len(input_population)))
    return [input_population, best_scores_over_time, average_scores_over_time]


# ## Perform Evolution
# check errors in routine file
errors = check_routine()
if len(errors) > 0:
    exit("Exit algorithm due to errors in routine file.")
# no errors, proceed with evolution

# create output folder
# get start time of current run in string format
start_time = datetime.datetime.now()
month = str(start_time.month)
if len(month) == 1:
    month = "0" + month
day = str(start_time.day)
if len(day) == 1:
    day = "0" + day
hour = str(start_time.hour)
if len(hour) == 1:
    hour = "0" + hour
minute = str(start_time.minute)
if len(minute) == 1:
    minute = "0" + minute
# build string containing start time
run_string = str(start_time.year) + "-" + month + "-" + day + "_" + hour + "-" + minute
# build output path and create folder
if args.prefix is None:
    out_path = "EvoDock_output_run_" + run_string
else:
    out_path = str(args.prefix) + "_EvoDock_output_run_" + run_string
try:
    Path(out_path).mkdir(parents=True, exist_ok=False)
except FileExistsError:
    suffix = 2
    folder_created = False
    while not folder_created:
        try:
            Path(out_path + "_" + str(suffix)).mkdir(parents=True, exist_ok=False)
            folder_created = True
        except FileExistsError:
            suffix += 1
    out_path = out_path + "_" + str(suffix)

# prepare and run evolution
# init variables for evolution
total_history = []
iterationCounts = [0, 0]
mds.set_values(original, out_path, protein_name, amino_acid_paths)

# generate random population
print("Generate initial Population...")
population = generate_initial_population(start_population_size, total_history, original)

# perform the whole evolution routine
routine_results = perform_routine(population, total_history, iterationCounts)
population = routine_results[0]
best_scores = routine_results[1]
average_scores = routine_results[2]

# analyze results
# count max target
count = 0
print("Best score: " + str(population[0][1]) + " | with difference to target score: "
      + str(get_individuals_score_relative_to_targetScore(population[0])))
for count_individual in population:
    if count_individual[1] == population[0][1]:
        count += 1
print("Number of mutants sharing best score: " + str(count))
print("Number of accepted mutations: " + str(iterationCounts[ITERATION_COUNT_MUTATION_INDEX]))
print("Number of accepted recombination products: " + str(iterationCounts[ITERATION_COUNT_RECOMBINATION_INDEX]))

# ## Save final Output
save_output(population, best_scores, average_scores)
end_time = datetime.datetime.now()
print("Run time: " + str(end_time - start_time))

# Plot results
if PLOT:
    plt.plot(best_scores)
    plt.plot(best_scores, "ob")
    plt.plot(average_scores)
    plt.plot(average_scores, "or")
    plt.show()
