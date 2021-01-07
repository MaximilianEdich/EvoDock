"""
EvoDock Version 0.1

created and developed by Maximilian Edich at Universitaet Bielefeld.
"""

VERSION = "0.21_01_07"

# region Imports
print("Import core modules...")

import random
from itertools import product
import argparse
import datetime
import ast
from pathlib import Path
import multiprocessing as multi_p

try:
    import MutateApplyScoreModule
except ImportError as e:
    MutateApplyScoreModule = None
    exit("Error: Import of \"MutateApplyScoreModule\" failed. Make sure to provide this module, since it is essential. "
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

# endregion

# start run time counter
start_time = datetime.datetime.now()
print("start timer")

# region String Vars and Indices
# initial settings file values
# task and termination
TARGET_SCORE = "-targetscore"
INFINITE = "INFINITE"
NEGATIVE_INFINITE = "NEGATIVE-INFINITE"
TERMINATE_TIME = "-terminate-by-time"
TERMINATE_SCORE = "-terminate-by-score-range"
TASK_TYPE = "-task"
TASK_LIGAND_BINDING = "LIGAND-BINDING"
TASK_PROTEIN_BINDING = "PROTEIN-BINDING"
TASK_THERMO_STABILITY = "THERMO-STABILITY"
TASKS = []
# TODO add here more tasks
TASKS.append(TASK_LIGAND_BINDING)

# init-pop-run-modes
INITIAL_POPULATION = "-init-pop-run-mode"
INITIAL_POPULATION_NEW_STOP = "create-and-quit"
INITIAL_POPULATION_LOAD = "load-and-evolve"
INITIAL_POPULATION_NEW_FULL_RUN = "create-and-evolve"
# init-pop-create-mode fold/mutate
INITIAL_POPULATION_CREATE_MODE = "-init-pop-create-mode"
CREATE_VIA_FOLD = "fold"
CREATE_VIA_MUTATE = "mutate"

# protein specifications and mutation settings
PROTEIN_PATH = "-prot-path"
RES_ID = "-res-id"
SUBSTITUTIONS = "-substitutions"
SYMMETRY = "-symmetry"
LOOK_UP_TABLE = "-look-up-table-path"
SUBST_MATRIX = "-substitution-matrix"  # TODO
BRUTE_FORCE = "-brute-force"
INCLUDE_MUTANT = "-include-mutant"

# optional technical settings
CPU_CORE_NUMBER = "-cpu"
MAX = "MAX"
SEED = "-seed"

# pipeline modules
MODULE_NAME_MUTATE = "-mutate"
MODULE_NAME_APPLY = "-apply"
MODULE_NAME_SCORE = "-score"
MODULE_NAME_FOLD = "-fold"
USE_SPECIFIC_MUTATE_OUT = "-use-specific-mutate-out"
USE_EXISTENT_MUTATE_OUT_PATH = "-use-existent-mutate-out-path"
SKIP_APPLICATION = "-skip-application"

# strings for routine commands
MUTATE = "mutate"
MUTATE_PLUS = "mutate+"
RECOMBINATION = "recombine"
RECOMBINATION_CLASSIC = "recombine-classic"
SELECT = "select"
LOOP = "loop"

# others
TRUE = "TRUE"
FALSE = "FALSE"

# indices for lists
ITERATION_COUNT_MUTATION_INDEX = 0
ITERATION_COUNT_RECOMBINATION_INDEX = 1

# technical fixed parameters
MAX_REC_DEPTH_GET_RANDOM_INDIVIDUAL = 200
MAX_REC_DEPTH_GET_NEW_RANDOM_MUTANT = 200
BREAK_LOOP_AFTER_MAX_REC_DEPTH = True

PRINT_OUT = False
# endregion

# region ArgParse
parser = argparse.ArgumentParser(description="EvoDock - An Evolutionary Algorithm for Protein Optimization")
parser.add_argument("-s", "--settings", type=str, required=True,
                    help="Path to the initial settings file.\nThe initial settings file must be a text file"
                         " containing only one single command per line.\nThe order of the commands is not important. "
                         "Spaces should only be used to separate a command and its parameter. See the documentation"
                         " for the full list of commands.")
parser.add_argument("-r", "--routine", type=str, required=True, help="Path to the routine file. See the documentation"
                                                                     "for a list of all commands.")
parser.add_argument("-o", "--out", type=str, required=False, help="Path to desired output destination.")
parser.add_argument("-doc", "--documentation", action='store_true', required=False,
                    help="Show documention of all commands.")
parser.add_argument("-p", "--protein_path", type=str, required=False, help="Path to input protein.")
parser.add_argument("-fn", "--folder_name", type=str, required=False, help="Set alternative folder name.")
args = parser.parse_args()
# TODO specify optional technical settings
# endregion

# region Documentation print
# TODO args.documentation print out all commands + commands of specified modules
if args.documentation:
    print("Documentation of EvoDock commands:")
    print("Initial settings file flags:")
    print("required flags:")
    print(INITIAL_POPULATION + " <mode> <arg>")
    print("\tSpecify how to get the initial population and how to proceed, where mode is either '"
          "" + INITIAL_POPULATION_NEW_STOP + "', '" + INITIAL_POPULATION_LOAD + "', or '"
          + INITIAL_POPULATION_NEW_FULL_RUN + "', with arg as the number of individuals in the initial "
                                              "population or the path to an initial population in case of loading it.")
    print("\n\texamples:")
    print("\t" + INITIAL_POPULATION + " " + INITIAL_POPULATION_NEW_STOP + " 100")
    print("\t" + INITIAL_POPULATION + " " + INITIAL_POPULATION_NEW_FULL_RUN + " 100")
    print("\t" + INITIAL_POPULATION + " " + INITIAL_POPULATION_LOAD + " load_population.txt")
    print(TASK_TYPE + " <mode>")
    print("\tSpecify the type of optimization problem. Valid inputs are:")
    for task_mode in TASKS:
        print("\t" + str(task_mode))
    print(TARGET_SCORE + " <n>")
    print("\tSpecify the target score, which defines the best fitness. The meaning of the score depends on the "
          "specified task. <n> >must be either a float "
          "value or '" + INFINITE + "' or '" + NEGATIVE_INFINITE + "'.")
    print(PROTEIN_PATH + " <path>")
    print("\tSpecify the path to the input PDB.")
    print(RES_ID + " <chain> <id>")
    print("\tEnter the residue chain and id of a residue of interest, that can be modified. "
          "Use this specification multiple times to specify multiple residues.")
    print(SUBSTITUTIONS + " <list>")
    print("\tSpecify allowed mutations, where <list> consists of single letter AA codes separated by a comma. "
          "Giving no argument allows all 20 AAs. Use this specification multiple times, where the order is respective "
          "to the order of '" + RES_ID + "'-commands.")
    print("\n\texample:")
    print("\t" + SUBSTITUTIONS + " D,L,A")
    print(MODULE_NAME_MUTATE + " <path>")
    print("\tSpecify the path to the used mutagenesis module.")
    print(MODULE_NAME_APPLY + " <path>")
    print("\tSpecify the path to the used application module.")
    print(MODULE_NAME_SCORE + " <path>")
    print("\tSpecify the path to the used evaluation module.")
    print("\noptional flags:")
    print(TERMINATE_SCORE + " <n>")
    print("\t")
    print(TERMINATE_TIME + " <time>")
    print("\t")
    print(CPU_CORE_NUMBER + " <n>")
    print("\tSpecify the number of cores used for parallel computation. If <n> is "
          "'" + MAX + "', all available cores are used.")
    print(SEED + " <seed>")
    print("\tSet the seed of the random generator for the run. Affects only random choices of the evolution process, "
          "does not affect pipeline modules.")
    print(LOOK_UP_TABLE + " <file-path>")
    print("\tSpecify path to a text file with list of individuals, like a history")
    print()
    print("\t")
    print()
    print("\t")
    print()
    print("\t")

    print("\nRoutine file commands:")
    print("TODO")
    # TODO load modules and print doc

    exit()

# endregion

# region Specification of initial essential and optional values, set default values
target_score = None
terminate_time = None
terminate_score_range = None
task = ""
initial_population_run_mode = ""
initial_population_create_mode = ""
start_population_size = 0
load_population_path = ""

reference_protein = [[], '']
protein_path = ""
amino_acid_paths = []
allowed_mutations = []
number_of_mutable_aa = 0

mds = MutateApplyScoreModule.MutateApplyScore()
module_name_mutate = ""
module_name_apply = ""
module_name_score = ""
module_name_fold = ""
use_specific_mutate_out = None
use_existent_mutate_out_path = None
skip_application = False
keep_improvements = False # TODO access
tabu_search = True

usable_cpu = multi_p.cpu_count()
look_up_table_path = ""
symmetry = None
brute_force_init_pop = True
include_mutants = []

# endregion

# region load files
# # try to load routine settings file
print("load initial settings file content and routine file content...")
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


# endregion

# region Error Message prints and exits
def error_msg_initial_population_options():
    """
    Prints out the error message for missing or wrong arguments for the initial population setting.
    Then it exits.
    """
    print(INITIAL_POPULATION + " " + INITIAL_POPULATION_NEW_STOP + " <initial-population-size>")
    print(INITIAL_POPULATION + " " + INITIAL_POPULATION_LOAD + " <path-to-population>")
    print(INITIAL_POPULATION + " " + INITIAL_POPULATION_NEW_FULL_RUN + " <initial-population-size>")
    exit()


def error_msg_create_init_pop_options():
    """
    Prints out the error message for missing or wrong arguments for the create mode of initial population setting.
    Then it exits.
    """
    print(CREATE_VIA_MUTATE + ":  Each pdb of a random mutant of the initial population is created via mutagenesis.")
    print(CREATE_VIA_FOLD + ":   Each pdb of a random mutant of the initial population is created via folding."
                            "folding.")
    print("See the documentation for more details.")
    exit()


# endregion

# region check initial settings file content and validate input format
print("check initial settings file content flags and validate input format...")
res_id_is_pdb_code = False

line_index = 0
module_params = []
# TODO check exceptions (IndexError)
for line in initial_settings_file_content:
    line_index += 1
    # prepare line
    line = line.strip()
    split_text = line.split(' ')
    if split_text[0] == INITIAL_POPULATION:
        # handle initial population settings
        if len(split_text) < 3:
            print("Error, " + INITIAL_POPULATION + " in initial settings file in line " + str(line_index) +
                  "requires two arguments!")
            print("These are the 3 possible commands:")
            error_msg_initial_population_options()
        else:
            if (split_text[1] == INITIAL_POPULATION_NEW_STOP or split_text[1] == INITIAL_POPULATION_LOAD
                    or split_text[1] == INITIAL_POPULATION_NEW_FULL_RUN):
                initial_population_run_mode = split_text[1]
            else:
                print("Error in line " + str(line_index) + ": Unknown mode, choose one of these options:")
                error_msg_initial_population_options()
            # set arguments and check if they are correct
            if (initial_population_run_mode == INITIAL_POPULATION_NEW_STOP
                    or initial_population_run_mode == INITIAL_POPULATION_NEW_FULL_RUN):
                try:
                    start_population_size = int(split_text[2])
                except ValueError:
                    exit("Error, <initial-population-size> value in initial settings file is not a valid number, "
                         "must be int: \""
                         + split_text[2] + "\" in line " + str(line_index) + " of the initial settings file.")
                if start_population_size < 1:
                    exit("Error, illegal value in line " + str(
                        line_index) + " of the initial settings file. Must be >= 1.")
            elif initial_population_run_mode == INITIAL_POPULATION_LOAD:
                load_population_path = str(split_text[2])
                initial_population_create_mode = INITIAL_POPULATION_LOAD

    elif split_text[0] == INITIAL_POPULATION_CREATE_MODE:
        # specify how the initial population is created, via mutation or via folding
        if initial_population_run_mode != INITIAL_POPULATION_LOAD:
            # not relevant, if loaded from file
            try:
                initial_population_create_mode = split_text[1]
            except IndexError:
                print("Error in line " + str(line_index) + ": Missing mode, choose one of these options:")
                error_msg_create_init_pop_options()
            if not (initial_population_create_mode == CREATE_VIA_FOLD
                    or initial_population_create_mode == CREATE_VIA_MUTATE):
                # no specified mode
                print("Error in line " + str(line_index) + ": Unknown mode, choose one of these options:")
                error_msg_create_init_pop_options()

    elif split_text[0] == TARGET_SCORE:
        # set the target score, any value possible
        try:
            if split_text[1] == INFINITE or split_text[1] == NEGATIVE_INFINITE:
                # set it to max number
                try:
                    import sys

                    if split_text[1] == INFINITE:
                        target_score = sys.maxsize
                    else:
                        target_score = sys.maxsize * -1
                except ImportError as e:
                    print("ERROR in EvoDock: Cannot set target score to max value!")
                    exit(e)
            else:
                target_score = float(split_text[1])
        except ValueError:
            exit("Error, " + TARGET_SCORE + " value in initial settings file is not a valid number, must be float: \""
                 + split_text[1] + "\" in line " + str(line_index) + " of the initial settings file.")
    elif split_text[0] == TERMINATE_SCORE:
        try:
            terminate_score_range = float(split_text[1])
        except ValueError as e:
            print("ERROR in EvoDock: " + TERMINATE_SCORE + " expects a float argument!")
            exit(e)
        except IndexError as e:
            print("ERROR in EvoDock: " + TERMINATE_SCORE + " expects a float argument!")
            exit(e)
    elif split_text[0] == TERMINATE_TIME:
        try:
            time_text = split_text[1]
            if time_text.find(':') < 0:
                print("ERROR in EvoDock: " + TERMINATE_TIME + " expects an argument of this format:")
                exit("'hhh:mm' or 'hh:mm' or 'h:mm', where 'h' and 'm' stand for hours and minutes respectively.")
            time_list = time_text.split(':')
            hours = time_list[0]
            mins = time_list[1]
            terminate_time = (hours * 60 * 60) + (mins * 60)
        except IndexError as e:
            print("ERROR in EvoDock: " + TERMINATE_TIME + " expects an argument of this format:")
            print("'hhh:mm' or 'hh:mm' or 'h:mm', where 'h' and 'm' stand for hours and minutes respectively.")
            exit(e)

    elif split_text[0] == USE_SPECIFIC_MUTATE_OUT:
        # specify if the specific mutate output should be used as input for the application module. May cause errors
        # if modules are not compatible!
        try:
            setting = str(split_text[1])
            if setting == TRUE:
                use_specific_mutate_out = True
            elif setting == FALSE:
                use_specific_mutate_out = False
            else:
                exit("Error in line " + str(line_index) + " of the initial settings file. Illegal argument!"
                                                          "\n" + setting + ". Expected are: 'TRUE' or 'FALSE'")
        except IndexError:
            exit("Error in line " + str(line_index) + " of the initial settings file. Argument missing!")
    elif split_text[0] == USE_EXISTENT_MUTATE_OUT_PATH:
        # specify path to already existing mutate out path, from were results are used, if available
        try:
            use_existent_mutate_out_path = str(split_text[1])
            if use_existent_mutate_out_path[-1:] != '/':
                use_existent_mutate_out_path += "/"
                print("WARNING: path of '"
                      "" + USE_EXISTENT_MUTATE_OUT_PATH + "' was extended to folder: " + use_existent_mutate_out_path)
            # TODO check if it is a folder
        except IndexError:
            exit("Error in line " + str(line_index) + " of the initial settings file. Argument missing!")
    elif split_text[0] == SKIP_APPLICATION:
        # set setting for skipping the application step and do right after mutagenesis the evaluation
        try:
            if split_text[1] == TRUE:
                skip_application = True
        except IndexError:
            exit("Error in line " + str(line_index) + " of the initial settings file. Argument missing!")

    elif split_text[0] == MODULE_NAME_MUTATE:
        # specified mutation module
        try:
            module_name_mutate = str(split_text[1])
        except IndexError:
            exit("Error in line " + str(line_index) + " of the initial settings file. Argument missing!")
    elif split_text[0] == MODULE_NAME_APPLY:
        # specified mutation module
        try:
            module_name_apply = str(split_text[1])
        except IndexError:
            exit("Error in line " + str(line_index) + " of the initial settings file. Argument missing!")
    elif split_text[0] == MODULE_NAME_SCORE:
        # specified mutation module
        try:
            module_name_score = str(split_text[1])
        except IndexError:
            exit("Error in line " + str(line_index) + " of the initial settings file. Argument missing!")
    elif split_text[0] == MODULE_NAME_FOLD:
        # specified mutation module
        try:
            module_name_fold = str(split_text[1])
        except IndexError:
            exit("Error in line " + str(line_index) + " of the initial settings file. Argument missing!")

    elif split_text[0] == SUBSTITUTIONS:
        # read allowed substitutions for mutations per position
        if len(split_text) == 1:
            # empty, allow all mutations TODO fix list on top and make it adjustable
            allowed_mutations.append(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                                      'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
        else:
            # not empty, allow only specified mutations
            allowed = split_text[1].split(',')
            allowed_mutations.append(allowed)
    elif split_text[0] == RES_ID:
        # residue id and/or path, depends on the used tools TODO argument check
        if len(split_text) == 3:
            # use chain + id
            res_id_is_pdb_code = True
            amino_acid_paths.append([split_text[1], int(split_text[2])])
        else:
            # use only number of residue, which is not the same as id (So 1 = 1st AA, even if its id is 16)
            amino_acid_paths.append(int(split_text[1]))
    elif split_text[0] == SYMMETRY:
        # if PDB consists of symmetric chains
        if len(split_text) == 3:
            symmetry = [split_text[1], split_text[2]]  # TODO use it
        else:
            print("ERROR in EvoDock: '" + SYMMETRY + " <a> <b>' requires two arguments:")
            exit("<a>: Chain that is given in res-id\n<b>: Chain id that replaces <a>.")
    elif split_text[0] == BRUTE_FORCE:
        # set if brute force is used to find new solutions in initial population creation.
        param = split_text[1]
        if param == TRUE:
            brute_force_init_pop = True
        else:
            brute_force_init_pop = False
    elif split_text[0] == INCLUDE_MUTANT:
        # set if brute force is used to find new solutions in initial population creation.
        if len(split_text) == 1:
            exit("ERROR in EvoDock: '" + INCLUDE_MUTANT + "' requires list as second argument, where items are"
                                                          "separated only with a comma!")
        else:
            mutant_AAs_list = split_text[1].split(',')
            include_mutants.append(mutant_AAs_list)

    elif split_text[0] == PROTEIN_PATH:
        # protein pdb code
        try:
            protein_path = str(split_text[1])
        except IndexError:
            exit("Error in line " + str(line_index) + " of the initial settings file. Argument missing!")
    elif split_text[0] == TASK_TYPE:
        # modification mode
        task = str(split_text[1])
        if not task in TASKS:
            exit("Error in line " + str(line_index) + " of the initial settings file. Unknown mode. Use one of these: "
                 + str(TASKS))
    elif split_text[0] == CPU_CORE_NUMBER:
        # number of usable cpu cores
        try:
            if split_text[1] == MAX:
                usable_cpu = multi_p.cpu_count()
            else:
                usable_cpu = int(split_text[1])
                if usable_cpu < 1:
                    usable_cpu = 1
                if usable_cpu > multi_p.cpu_count():
                    usable_cpu = multi_p.cpu_count()
        except IndexError as e:
            exit("Error in line " + str(line_index) + ": Missing argument, the number of usable CPU cores.")
        except ValueError as e:
            exit("Error in line " + str(line_index) + ": Argument must be a positive integer or 'MAX'!")
    elif split_text[0] == SEED:
        # set seed for random number generator
        try:
            random.seed(split_text[1])
        except IndexError as e:
            exit("Error in line " + str(line_index) + ": Missing argument, the seed.")
    elif split_text[0] == LOOK_UP_TABLE:
        try:
            look_up_table_path = split_text[1]
        except IndexError as e:
            exit("Error in line " + str(line_index) + ": Missing argument, path to look-up-table.")
    elif split_text[0] == "":
        # allow blank lines
        pass
    elif split_text[0][0] == "#":
        # allow comments
        pass

    elif (split_text[0] == MutateApplyScoreModule.MODULE_PARAM_MUTATE
          or split_text[0] == MutateApplyScoreModule.MODULE_PARAM_APPLY
          or split_text[0] == MutateApplyScoreModule.MODULE_PARAM_SCORE
          or split_text[0] == MutateApplyScoreModule.MODULE_PARAM_FOLD):
        module_params.append(split_text)

    else:
        exit("Error, undefined keywords in line " + str(line_index) + " of the initial settings file: " + split_text[0])

    # if input protein is given via command line
    if args.protein_path is not None:
        protein_path = args.protein_path

# endregion

# region catch undefined essential values
number_of_mutable_aa = len(amino_acid_paths)
if res_id_is_pdb_code:
    for entry in amino_acid_paths:
        if isinstance(entry, int):
            exit("ERROR in EvoDock: Unconsistent res-id format. Use either consecutive numbering order ids or "
                 "PDB format with chain + residue number. Use only one format for all 'res-id' occurrences!")
elif symmetry is not None:
    exit("ERROR in EvoDock: Symmetry can only be used, if residues use chain information!")
if len(amino_acid_paths) != len(allowed_mutations):
    print("Error in initial settings file: Number of residue ids and number of lists with allowed substitutions "
          "are not equal!")
    exit("use equal amount of '" + RES_ID + "' and '" + SUBSTITUTIONS + "'!")
if initial_population_run_mode == "":
    print("Error in initial settings file: You have to specify \"" + INITIAL_POPULATION + "\"! Use these options:")
    error_msg_initial_population_options()
if initial_population_create_mode == "":
    initial_population_create_mode = CREATE_VIA_MUTATE
if target_score is None:
    exit("Error in initial settings file: You have to specify the \"" + TARGET_SCORE + " f\"!")
if protein_path != "":
    try:
        test_file = open(protein_path)
        test_file.close()
    except FileNotFoundError:
        exit("Error in initial settings file: Specified protein file not found!")
    except PermissionError:
        exit("Error in initial settings file: No permission to open file!")
else:
    exit("Error in initial settings file: You have to specify the \"" + PROTEIN_PATH + " leading to your .pdb\"!")
if module_name_mutate == "":
    exit("Error in initial settings file: You have to specify the mutation-module via \"" + MODULE_NAME_MUTATE + "\"!")
if module_name_apply == "":
    exit(
        "Error in initial settings file: You have to specify the application-module via \"" + MODULE_NAME_APPLY + "\"!")
if module_name_score == "":
    exit("Error in initial settings file: You have to specify the scoring-module via \"" + MODULE_NAME_SCORE + "\"!")
if module_name_fold == "" and initial_population_create_mode == CREATE_VIA_FOLD:
    exit("Error in initial settings file: You have to specify the folding-module "
         "via \"" + MODULE_NAME_FOLD + "\" or create the initial population via mutagenesis instead of folding!")
if not (use_specific_mutate_out == True or use_specific_mutate_out == False or use_specific_mutate_out is None):
    exit("Error in initial settings file: You have to specify which mutagenesis output is used via"
         " \"" + USE_SPECIFIC_MUTATE_OUT + " <b>\", with <b> = 'TRUE' or 'FALSE'!")
# if list "include_mutants" contains elements
if include_mutants:
    # check if all included mutants have same number of mutable AAs as reference protein
    for item in include_mutants:
        if len(item) != number_of_mutable_aa:
            exit("ERROR in EvoDock: One included mutant has not expected format."
                 " Expected were " + str(number_of_mutable_aa) + " AAs, input was: " + str(len(item)))

# endregion

# region create output folder
# get start time of current run in string format
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
# take path to destination
if args.out is None:
    out_destination = ""
else:
    out_destination = str(args.out)

# build output path and create folder
out_path = out_destination + "EvoDock_output_run_" + run_string
if args.folder_name is not None:
    out_path = out_destination + args.folder_name
try:
    Path(out_path).mkdir(parents=True, exist_ok=False)
except FileExistsError:
    # repeat process with different number suffixes until successful
    suffix = 2
    folder_created = False
    while not folder_created:
        try:
            Path(out_path + "_" + str(suffix)).mkdir(parents=True, exist_ok=False)
            folder_created = True
        except FileExistsError:
            suffix += 1
    out_path = out_path + "_" + str(suffix)
except FileNotFoundError as e:
    exit("Directory not found. " + str(e))
except PermissionError as e:
    exit("No permission for directory. " + str(e))
# endregion

# region MASM and amino acid paths set up
# init MASM with modules
print("initialize MutateApplyScoreModule: import specified modules...")
# import modules
MutateApplyScoreModule.init(module_name_mutate, module_name_apply, module_name_score, module_name_fold)
# check if each module is assigned correctly
MutateApplyScoreModule.check_imported_modules()

# apply module parameters
print("\nset up specified modules with specific module parameters...")
for module_param_list in module_params:
    returned_value = MutateApplyScoreModule.handle_module_params(module_param_list)
    if returned_value is not None:
        # TODO maybe work with return values, catch warnings
        pass

# validate and prepare tools and maybe update protein path
MutateApplyScoreModule.validate_module_data(protein_path, out_path, skip_application)
if use_specific_mutate_out is None or use_specific_mutate_out == True:
    use_specific_mutate_out = MutateApplyScoreModule.check_compatibility_mutate_apply()
MutateApplyScoreModule.prepare_tool(protein_path, out_path, skip_application)
protein_path = MutateApplyScoreModule.preparation_result_path(protein_path, out_path)
print(protein_path)

# reformat res-ids if required
if res_id_is_pdb_code:
    # update symmetry information
    if symmetry is not None:
        symmetry_list = []
        for aap in amino_acid_paths:
            if aap[0] == symmetry[0]:
                symmetry_list.append([symmetry[1], aap[1]])
        symmetry = MutateApplyScoreModule.get_reformatted_amino_acids(protein_path, symmetry_list)
        amino_acid_paths = MutateApplyScoreModule.get_reformatted_amino_acids(protein_path, amino_acid_paths)
        # check if AAs identical
        sym_AAs = MutateApplyScoreModule.get_ref_protein_amino_acids(protein_path, symmetry)
        path_AAs = MutateApplyScoreModule.get_ref_protein_amino_acids(protein_path, amino_acid_paths)
        if sym_AAs == path_AAs:
            print("Symmetry validated, AAs on both sides identical.")
        else:
            exit("ERROR in EvoDock: Symmetry was set, but AAs are not identical!" + str(sym_AAs) + str(path_AAs))
    else:
        amino_acid_paths = MutateApplyScoreModule.get_reformatted_amino_acids(protein_path, amino_acid_paths)


# get AAs of reference protein
reference_protein[0] = MutateApplyScoreModule.get_ref_protein_amino_acids(protein_path, amino_acid_paths)
# quick reference for several non-changing params
mds.set_values(reference_protein, out_path, protein_path, amino_acid_paths, use_specific_mutate_out,
               use_existent_mutate_out_path, skip_application, keep_improvements, symmetry)

# endregion


print("IMPORTS AND FILE LOADINGS DONE!\n")
# init reference protein individual and print loaded info summary
print("reference protein's AAs + amino acid paths:")
print(reference_protein[0])
print(amino_acid_paths)
print("\nNumber of cores, detected for multi-processing: " + str(multi_p.cpu_count()))
print("Use of multi-processing is restricted to: " + str(usable_cpu) + "\n")

# region look up table
look_up_scores = []
if look_up_table_path != "":
    print("Load look up table")
    content = []
    try:
        look_up_file = open(look_up_table_path, 'r')
        content = look_up_file.readlines()
        look_up_file.close()

    except FileNotFoundError:
        exit("File for look up table not found! " + str(look_up_table_path))
    except PermissionError:
        exit("Error: File for look up table can not be opened, permission denied!")
    # TODO check format
    for line in content[1:]:
        load_individual = [[], 0]
        content = line.strip().split(';')
        load_individual[1] = float(content[1])
        load_individual[0] = ast.literal_eval(content[0])
        look_up_scores.append(load_individual)


# endregion

# region Score Functions
def get_and_write_score(target_individual):
    """
    Calculates the score of an individual and writes it into its score value.
    :param target_individual: An individual in the form of [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid.
    :return: The individual with its new calculated score as the individuals fitness.
    """

    for look_up in look_up_scores:
        if target_individual[0] == look_up[0]:
            target_individual[1] = look_up[1]
            return target_individual

    # calculate the score (by using external software) TODO check parameter and not initial setting
    if initial_population_create_mode == CREATE_VIA_MUTATE:
        score = MutateApplyScoreModule.get_fitness_score(target_individual, mds, False)
    else:
        score = MutateApplyScoreModule.get_fitness_score(target_individual, mds, True)
    # write the score into the individual and return score
    target_individual[1] = score

    return target_individual


def update_history_scores(input_population):
    """
    Updates the scores of the history. Necessary, if multiprocessing was used. Otherwise the data updates automatically
    in the correct object. Call this right after scoring and joining multi processes.
    :param input_population: Actual population right after a scoring process.
    :return: None.
    """
    # check all individuals from population, since these are new
    for individual in input_population:
        for history_entry in total_history:
            if history_entry[0] == individual[0]:
                history_entry[1] = individual[1]
    return


def score_on_multi_core(score_population):
    """
    Takes a new gathered population and scores each individual via multi processing, if enabled.
    :param score_population: New population with individuals without score.
    :return: The input population but with scored individuals.
    """
    if usable_cpu > 1:
        # calculate score of initial population in parallel processes
        # create pool of processes and map tasks
        pool = multi_p.Pool(usable_cpu)
        result = pool.map(get_and_write_score, score_population)
        # close processes and wait for end
        pool.close()
        pool.join()
        update_history_scores(result)

        return result
    else:
        # calculate score in single run
        for individual in score_population:
            get_and_write_score(individual)

        return score_population


# endregion

# region Create Population - Functions
def get_random_genes():
    """
    Get a random combination of possible genes (in terms of genetic algorithms).
    :return: A list, containing random amino acids chosen by the specified pattern.
    """
    genes = []
    for aa in range(number_of_mutable_aa):
        # for each AA position of interest, add a random amino acid from specified lists
        genes.append(random.choice(allowed_mutations[aa]))

    return genes


def get_brute_forced_genes():
    """
    Create all possible combinations of gene sets and test them on total history.
    :return: A new unique set if one was found, otherwise None.
    """
    # check total history on all combinations
    input_list = []
    for aa in range(number_of_mutable_aa):
        input_list.append(allowed_mutations[aa])
    print("brute force!")
    print(input_list)
    all_combinations = list(product(*[]))
    for comb in all_combinations:
        genes = list(comb)
        # test history
        is_old = False
        for history_entry in total_history:
            if history_entry[0] == genes:
                is_old = True
                break
        if not is_old:
            return genes

    return


def get_random_individual(rec_depth=0):
    """
    Generate a new individual with random gene combination and a score of 0.
    :param rec_depth: The recursion depth. If a generated individual is not new, the
    functions is recalled with increased recDepth as long as recDepth is smaller
    than the maximal recursion depth. Use the default value of 0.
    :return: A new and random individual in the form [g, s], where s is the score
    value and g a list of genes (in terms of genetic algorithms) in form of [g1, g2, ..., g_n], where each
    gene represents an amino acid. If the maximal recursion depth is reached, before a new individual was found,
    None is returned.
    """
    if rec_depth > MAX_REC_DEPTH_GET_RANDOM_INDIVIDUAL:
        # return None, if maximum recursion depth is reached
        return

    new_individual = [get_random_genes(), '']
    # check if the new individual appears in history
    for history_entry in total_history:
        if history_entry[0] == new_individual[0]:
            # individual not new, repeat with new random individual
            try:
                return get_random_individual(rec_depth + 1)
            except RecursionError as error:
                print("Please lower the new-random-individual-recursion-limit in the technical settings!"
                      "Error Message: " + str(error))
                # return None, if maximum recursion depth is reached
                # if brute-force is on, try finding possible solutions
                if brute_force_init_pop:
                    brute_forced_genes = get_brute_forced_genes()
                    if brute_forced_genes is None:
                        # no further solutions, quit full initial population generation
                        return
                    else:
                        new_individual = [brute_forced_genes, '']
                        return new_individual

                return

    # individual is new, add to history
    total_history.append(new_individual)

    return new_individual


def generate_initial_population(number_of_initial_individuals, ref_prot):
    """
    Generate the initial population with random and new individuals.
    :param number_of_initial_individuals: The number ( > 0) of total individuals in the generated population.
    :param ref_prot: An individual with the initial amino acids of interest. The reference protein.
    individual is added as first one to the new population.
    :return: A new population. A population is a list of individuals.
    """
    # initialize population with reference protein and included mutants
    new_population = [ref_prot]
    for item in include_mutants:
        new_population.append([item, ''])
    for item in new_population:
        total_history.append(item)
    print("Population with reference protein and included mutants:")
    print(new_population)
    # generate new individuals with genes only, without scoring yet
    for n in range(number_of_initial_individuals - len(new_population)):
        # generate new individual
        new_individual = get_random_individual()
        if new_individual is not None:
            # add to population, if new individual is available
            new_population.append(new_individual)

    # calculate scores and return result
    return score_on_multi_core(new_population)


# endregion

# region Mutation Functions
def get_new_random_mutant(parent, rec_depth=0):
    """
    Creates a copy of the given individual (parent) and chooses by random choice a mutable position.
    For this position a random character of the allowed characters on this position is chosen.
    :param parent: The parents gene combination is copied, before a mutation is applied.
    :param rec_depth: The recursion depth. If a generated individual is not new, the
    functions is recalled with increased recDepth as long as recDepth is smaller
    than the maximal recursion depth.
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
    new_mutant = [parent[0].copy(), '']
    # choose random position and mutate it until it is not equal to its parent
    while new_mutant[0] == parent[0]:
        r = random.choice(range(number_of_mutable_aa))
        new_mutant[0][r] = random.choice(allowed_mutations[r])
    # check the history
    for individual in total_history:
        if individual[0] == new_mutant[0]:
            # repeat with new random mutation
            try:
                return get_new_random_mutant(parent, rec_depth + 1)
            except RecursionError as error:
                print("Please lower the new-mutant-recursion-limit in the technical settings!"
                      "Error Message: " + str(error))
                return

    total_history.append(new_mutant)
    evolution_track.append([parent, new_mutant])
    iterationCounts[ITERATION_COUNT_MUTATION_INDEX] += 1
    return new_mutant


def get_random_mutants(parent, number_of_new_mutants):
    """
    Get from a given parent individual several new mutants.
    :param parent: The parents gene combination is copied, before mutations are applied.
    :param number_of_new_mutants: The maximum number of new mutants generated by this function.
    :return: A list of new individuals. Since the mutation function can fail, this list can be
    empty. If all mutations were successful, the list contains 'numberOfNewMutants' items.
    """
    mutants = []
    # generate several mutants
    for m in range(number_of_new_mutants):
        # generate new individual, if available
        new_mutant = get_new_random_mutant(parent)
        if new_mutant is not None:
            # only keep new individuals
            mutants.append(new_mutant)
        elif BREAK_LOOP_AFTER_MAX_REC_DEPTH:
            break

    return mutants


def mutate_population(input_population, number_of_new_mutants):
    """
    Iterate through the whole population and create several mutants for each. All the original individuals from
    the inputPopulation and all new are transferred into the returned new population.
    :param input_population: The population from which each individual is passed trough the mutation process.
    :param number_of_new_mutants: The maximum number of new mutants generated by this function.
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
        mutants = get_random_mutants(parent, number_of_new_mutants)
        # add new mutants to population
        for mutant in mutants:
            new_population.append(mutant)

    # calculate scores of new individuals
    return score_on_multi_core(new_population)


def mutate_and_keep_improvements(input_population, number_of_new_mutants):
    """
    Iterate through the whole population and create several mutants for each. Only the original individuals from
    the inputPopulation and mutants with an improved score are transferred into the returned new population.
    :param input_population: The population from which each individual is passed trough the mutation process.
    :param number_of_new_mutants: The maximum number of new mutants generated by this function.
    :return: A new population containing the whole inputPopulation and improved mutants.
    """
    if PRINT_OUT:
        print("\ncreate Mutants from individuals and only keep improvements")
    new_population = []
    offsprings = []
    offspring_population = []
    # for each individual in population
    for parent in input_population:
        # keep the individual in population
        new_population.append(parent)
        # generate new mutants
        mutants = get_random_mutants(parent, number_of_new_mutants)
        offsprings.append(mutants)
        for mutant in mutants:
            offspring_population.append(mutant)

    # score offsprings
    if usable_cpu > 1:
        # calculate score in parallel processes
        # create pool of processes and map tasks
        pool = multi_p.Pool(usable_cpu)
        offspring_population = pool.map(get_and_write_score, offspring_population)
        # close processes and wait for end
        pool.close()
        pool.join()
        update_history_scores(offspring_population)
    else:
        # calculate score in single run
        for individual in offspring_population:
            get_and_write_score(individual)
    # compare offsprings to parent and only keep improvements
    index = 0
    for parent in input_population:
        for offspring in offsprings[index]:
            if get_individuals_score_relative_to_targetScore(offspring) < \
                    get_individuals_score_relative_to_targetScore(parent):
                # if the mutant is an improvement, keep it in new population
                new_population.append(offspring)
        index += 1

    return new_population


# endregion

# region Recombination Functions
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


def perform_recombination_classic(input_population, repetitions):
    """
    Perform a classic cross over in terms of Genetic Algorithms with the whole population. All recombination
    are kept.
    :param input_population: The population from which each individual is passed trough the recombination process.
    :param repetitions: Number, how many mating partners are chosen for each individual in the population.
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
            for history_entry in total_history:
                if history_entry[0] == new_genes:
                    is_new = False
                    break
            # if new, create individual and calculate score
            if is_new:
                new_individual = [new_genes, '']
                # add to history and population
                total_history.append(new_individual)
                iterationCounts[ITERATION_COUNT_RECOMBINATION_INDEX] += 1
                new_population.append(new_individual)

    # calculate scores of new individuals
    return score_on_multi_core(new_population)


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


def perform_uniform_recombination(input_population, repetitions):
    """
    Perform a uniform recombination in terms of Genetic Algorithms with the whole population. All recombination
    are kept. For each position, a coin flip chooses the parent.
    :param input_population: The population from which each individual is passed trough the recombination process.
    :param repetitions: Number, how many mating partners are chosen for each individual in the population.
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
                print("")
            # check history
            is_new = True
            for history_entry in total_history:
                if history_entry[0] == new_genes:
                    is_new = False
                    break
            # if new, create individual and calculate score
            if is_new:
                new_individual = [new_genes, '']
                # add to history and population
                total_history.append(new_individual)
                iterationCounts[ITERATION_COUNT_RECOMBINATION_INDEX] += 1
                new_population.append(new_individual)

    # calculate scores of new individuals
    return score_on_multi_core(new_population)


# endregion

# region Selection Functions
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
    for individual in input_population:
        score += individual[1]

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
            print("Warning: 'fractionPercent' value in 'SelectFittestByFraction' is out of range!"
                  "Whole population was kept.")
        return input_population
    if random_picks_percent < 0 or random_picks_percent > 1:
        if PRINT_OUT:
            print("Warning: 'randomPicksPercent' value in 'SelectFittestByFraction' is out of range!"
                  "Whole population was kept.")
        return input_population

    # calculate number of individuals being selected
    keep_int = round(fraction_percent * len(input_population))
    random_picks = int(random_picks_percent * keep_int)
    if keep_int < 1:
        keep_int = 1

    return select_fittest_by_number(keep_int, random_picks, input_population)


# endregion

# region Save Functions
def save_population_list(individuals, name):
    """
    Take a list containing only individuals and save it with the given name included in file name and header.
    The header contains also the number of entries and consists only of the top line.
    :param individuals: List with individuals, where each individual is a list of a genes list and a score.
    :param name: name is included in the file name and in the header.
    :return: None.
    """
    individuals.sort(reverse=False, key=get_individuals_score_relative_to_targetScore)
    file_content = "Entries in " + str(name) + ": " + str(len(individuals)) + "\n"
    for individual in individuals:
        file_content += str(individual[0]) + ";" + str(individual[1]) + "\n"

    out_file = open(out_path + "/EvoDock_" + str(name) + ".txt", 'w')
    out_file.write(file_content)
    out_file.close()
    return


def save_output(input_population, best_scores_over_time, average_scores_over_time, extra_info=""):
    """
    Save the current results, history, and run information, like input files and settings.
    :return: None.
    """

    best_score = input_population[0][1]
    individual_count = 0
    for individual in population:
        if individual[1] == population[0][1]:
            individual_count += 1

    # output results
    results_file_content = ""
    if extra_info != "":
        results_file_content += extra_info + "\n"
    results_file_content += "Final Population Size: " + str(len(input_population)) + "\n" \
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
    save_population_list(total_history, "history")

    # output run info
    # # basic run info
    current_time = datetime.datetime.now()
    run_info_file_content = "Run information of run " + run_string + "\n"
    run_info_file_content += "Run time: " + str(current_time - start_time) + "\n"
    # # version
    run_info_file_content += "\n\nEvoDock Version: " + VERSION + "\n"
    # # reference protein AAs
    run_info_file_content += "\n\n\nReference Protein amino acids:\n"
    for aa in reference_protein:
        run_info_file_content += str(aa) + "\n"
    # # Initial Settings File - Content
    run_info_file_content += "\n\nInitial Settings File - Content:\n"
    for file_line in initial_settings_file_content:
        run_info_file_content += file_line
    # # Routine File - Content
    run_info_file_content += "\n\n\nRoutine File - Content:\n"
    for file_line in routine_file_content:
        run_info_file_content += file_line
    # # External Settings
    run_info_file_content += "\n\n\nExternal Tool Settings File Contents:\n"
    # TODO get version numbers of MutateApplyScore Modules used

    run_info_file = open(out_path + "/EvoDock_run_information.txt", 'w')
    run_info_file.write(run_info_file_content)
    run_info_file.close()

    return


# endregion

# region Routine Check and Performing
def check_routine():
    """
    Iterates through the specified routine file and checks, if all commands
    are valid. The commands are not executed.
    :return: A list of errors, represented as strings.
    """
    # init error list and index
    routine_errors = []
    index = 0
    in_loop = False
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
                in_loop = False
        else:
            if split_command[0][0] == '#':
                # comment line
                pass
            elif split_command[0] == MUTATE:
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
                    selection_param1 = float(split_command[1])
                    # catch negative values
                    if selection_param1 < 0:
                        error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                            selection_param1) + ". Must be >= 0"
                    # check second parameter if existent
                    if len(split_command) > 2:
                        try:
                            selection_param2 = float(split_command[2])
                            # catch negative values
                            if selection_param2 < 0:
                                error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                                    selection_param2) + ". Must be >= 0"
                            # check if second parameter is float for fraction selection
                            if selection_param1 is float:
                                if selection_param2 > 1:
                                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                                        selection_param2) + ". Must be >= 0 and <= 1"
                        except ValueError:
                            error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                                split_command[2]) + ". Must be float or int"
                except ValueError:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        split_command[1]) + ". Must be float or integer"

            elif split_command[0] == LOOP:
                if in_loop:
                    # nested loop detected
                    # TODO warning collection
                    print("WARNING: Nested loop in routine file. Only the inner loop will be repeated, since nested "
                          "loops do not work!")
                in_loop = True
                try:
                    repeat_number = int(split_command[1])
                    if repeat_number < 1:
                        error = "Error in routine file, illegal value in line " \
                                + str(index) + ": " + str(split_command[1] + ". Must be an int, >= 1")
                except ValueError:
                    error = "Error in routine file, illegal value in line " + str(index) + ": " + str(
                        split_command[1]) + ". Must be a positive integer"
            else:
                undefined = True
        if undefined:
            error = "Error in routine file, undefined keyword in line " + str(index) + ": " + str(split_command[0])

        if error is not None:
            print(error)
            routine_errors.append(error)

    return routine_errors


def terminate_by_time():
    """
    Check if the runtime exceeds the termination time.
    :return: True, if runtime exceeds specified termination time. Otherwise, and in case of no specified termination
    time, return False.
    """
    if terminate_time is None:
        return False
    # TODO check time format and test
    # if (datetime.datetime.now() - start_time) > int(terminate_time):
    # print("EvoDock: Termination caused by specified runtime.")
    # return True

    return False


def terminate_by_score(input_population):
    """
    Check if at least one individual has a score, whose difference to the target score is less then the termination
    by score range.
    :param input_population: Population, which has to be checked for the condition.
    :return: True, if at least one indidivual's score is in termination range to the target score. Otherwise, and
    if termination score range is not specified, return False.
    """
    if terminate_score_range is None:
        return False
    for individual in input_population:
        if abs(individual[1] - target_score) <= terminate_score_range:
            print("EvoDock: Termination caused by score.")
            print("One individual's score is in specified range of target score.")
            return True
    return False


def perform_routine(input_population):
    """
    Iterates through the specified routine file and executes each command in the given order.
    :param input_population: The population on which the evolution will be performed.
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
    selection_step = 0
    loop_number = 0
    loop_jump = 0
    while routine_step < len(routine_file_content):
        routine = routine_file_content[routine_step]
        # strip line and separate it by spaces
        split_command = routine.strip().split(' ')
        if split_command[0] != '' and not split_command[0][0] == "#":
            print("Population Size: " + str(len(input_population)))
            print("Do " + routine.strip() + "...")
        if len(split_command) >= 2:
            if split_command[0] == MUTATE:
                # perform mutation
                mutation_number = int(split_command[1])
                input_population = mutate_population(input_population, mutation_number)
            elif split_command[0] == MUTATE_PLUS:
                # perform mutation
                mutation_number = int(split_command[1])
                input_population = mutate_and_keep_improvements(input_population, mutation_number)
            elif split_command[0] == RECOMBINATION:
                # perform uniform recombination
                repetition_number = int(split_command[1])
                input_population = perform_uniform_recombination(input_population, repetition_number)
            elif split_command[0] == RECOMBINATION_CLASSIC:
                # perform uniform recombination
                repetition_number = int(split_command[1])
                input_population = perform_recombination_classic(input_population, repetition_number)
            elif split_command[0] == SELECT:
                # select number or fraction of mutants
                selection_param1 = float(split_command[1])
                if len(split_command) > 2:
                    selection_param2 = float(split_command[2])
                else:
                    selection_param2 = 0

                # check selection method
                if selection_param1 > 1:
                    # select by number
                    selection_param2 = int(selection_param2)
                    input_population = select_fittest_by_number(int(selection_param1), int(selection_param2),
                                                                input_population)
                elif selection_param1 >= 0:
                    # select by fraction
                    input_population = select_fittest_by_fraction(selection_param1, selection_param2, input_population)

                best_scores_over_time.append(input_population[0][1])
                average_scores_over_time.append(ger_average_score(input_population))
                selection_step += 1
                save_output(input_population, best_scores_over_time, average_scores_over_time,
                            "-intermediate result- step " + str(selection_step))
                # if improvements should be saved as pdb, check population for all improvements compared to ref protein
                if keep_improvements:
                    improvements = []
                    ref_score = get_individuals_score_relative_to_targetScore(reference_protein)
                    for indiv in input_population:
                        if get_individuals_score_relative_to_targetScore(indiv) < ref_score:
                            improvements.append(indiv)
                    MutateApplyScoreModule.save_improvements(improvements, mds)

            elif split_command[0] == LOOP:
                # set repeat number and point
                loop_number = int(split_command[1]) - 1
                loop_jump = routine_step + 1
        elif routine == "\n":
            if loop_number > 0:
                loop_number -= 1
                # lower by one, because it increases after this line
                routine_step = loop_jump - 1
        routine_step += 1
        if routine_step >= len(routine_file_content):
            if loop_number > 0:
                loop_number -= 1
                routine_step = loop_jump
        # check termination rules
        if terminate_by_time():
            break
        if terminate_by_score(input_population):
            break

    if PRINT_OUT:
        print("Population Size: " + str(len(input_population)))
    return [input_population, best_scores_over_time, average_scores_over_time]


# endregion

# ## EvoDock Core
# region check errors in routine file
errors = check_routine()
if len(errors) > 0:
    exit("Exit algorithm due to errors in routine file.")
# no errors, proceed with evolution
# endregion

# region init variables and prepare evolution
total_history = []
evolutionary_trees = []
evolution_track = []
iterationCounts = [0, 0]

if look_up_scores:
    save_population_list(look_up_scores, "used_look_up_table")

# endregion

# region generate random or load initial population
population = []
if initial_population_run_mode == INITIAL_POPULATION_NEW_STOP:
    print("Generate initial Population...")
    population = generate_initial_population(start_population_size, reference_protein)
    # create input lists for saving function
    best_scores = [population[0][1]]
    average_scores = [population[0][1]]
    population.sort(reverse=False, key=get_individuals_score_relative_to_targetScore)
    best_scores.append(population[0][1])
    average_scores.append(ger_average_score(population))
    # save
    save_output(population, best_scores, average_scores, "initial population only")
    save_population_list(population, "init_pop_via_" + initial_population_create_mode)

    # save improvements PDBs
    print("new...")
    if keep_improvements:
        print("true")
        improvements = []
        ref_score = get_individuals_score_relative_to_targetScore(reference_protein)
        for indiv in population:
            if get_individuals_score_relative_to_targetScore(indiv) < ref_score:
                improvements.append(indiv)
        print(improvements)
        MutateApplyScoreModule.save_improvements(improvements, mds)

    end_time = datetime.datetime.now()
    # exit
    print("Run time: " + str(end_time - start_time))
    exit("Initial population was successfully created. EvoDock stops here, since that was specified in "
         "the initial settings.")
elif initial_population_run_mode == INITIAL_POPULATION_LOAD:
    # load file
    population_file_content = ""
    try:
        population_file = open(load_population_path, 'r')
        population_file_content = population_file.readlines()
        population_file.close()
    except FileNotFoundError:
        exit("Error: Initial population file does not exist!")
    except PermissionError:
        exit("Error: Initial population file can not be opened, permission denied!")
    # translate into population TODO check format
    for entry in population_file_content[1:]:
        load_individual = [[], 0]
        content = entry.strip().split(';')
        load_individual[1] = float(content[1])
        load_individual[0] = ast.literal_eval(content[0])
        population.append(load_individual)
    # save copy in output folder
    save_population_list(population, "init_pop_via_" + initial_population_create_mode)
elif initial_population_run_mode == INITIAL_POPULATION_NEW_FULL_RUN:
    # create population and save initial population
    print("Generate initial Population...")
    population = generate_initial_population(start_population_size, reference_protein)
    save_population_list(population, "init_pop_via_" + initial_population_create_mode)

# endregion

# region perform evolution routine
routine_results = perform_routine(population)
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
if look_up_scores != []:
    extended_histroy = []
    for entry in total_history:
        extended_histroy.append(entry)
    for entry in look_up_scores:
        contained = False
        for added_entry in extended_histroy:
            if entry[0] == added_entry[0]:
                contained = True
                break
        if not contained:
            extended_histroy.append(entry)
    save_population_list(extended_histroy, "history_plus_look_up_table")
end_time = datetime.datetime.now()
print("Run time: " + str(end_time - start_time))

# build evolution tree
evolutionary_trees = []
for individual in population:
    for entry in evolution_track:
        if entry[1][0] == individual[0]:
            evolutionary_trees.append([individual, entry[0]])
run = True
while run:
    run = False
    for tree in evolutionary_trees:
        for entry in evolution_track:
            if tree[len(tree) - 1][0] == entry[1][0]:
                run = True
                tree.append(entry[0])
for tree in evolutionary_trees:
    print(len(tree))
    print(tree)
    pass

# Plot results
if PLOT:
    plt.plot(best_scores)
    plt.plot(best_scores, "ob")
    plt.plot(average_scores)
    plt.plot(average_scores, "or")
    plt.show()

# endregion

if __name__ == "__main__":
    print("MAIN")
