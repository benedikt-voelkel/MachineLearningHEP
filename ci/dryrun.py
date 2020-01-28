#############################################################################
##  © Copyright CERN 2018. All rights not expressly granted are reserved.  ##
##                 Author: Gian.Michele.Innocenti@cern.ch                  ##
## This program is free software: you can redistribute it and/or modify it ##
##  under the terms of the GNU General Public License as published by the  ##
## Free Software Foundation, either version 3 of the License, or (at your  ##
## option) any later version. This program is distributed in the hope that ##
##  it will be useful, but WITHOUT ANY WARRANTY; without even the implied  ##
##     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    ##
##           See the GNU General Public License for more details.          ##
##    You should have received a copy of the GNU General Public License    ##
##   along with this program. if not, see <https://www.gnu.org/licenses/>. ##
#############################################################################

import argparse

# Collect all processers
from machine_learning_hep.multiprocesser import MultiProcesser
from machine_learning_hep.processer import Processer
from machine_learning_hep.processerDhadrons import ProcesserDhadrons
from machine_learning_hep.processerDhadrons_mult import ProcesserDhadrons_mult
from machine_learning_hep.processerDhadrons_jet import ProcesserDhadrons_jet

# Collect all optimisers
from machine_learning_hep.optimiser import Optimiser

# Collect all analyzers
from machine_learning_hep.analysis.analyzer_manager import AnalyzerManager
from machine_learning_hep.analysis.analyzer_Dhadrons import AnalyzerDhadrons
from machine_learning_hep.analysis.analyzer_Dhadrons_mult import AnalyzerDhadrons_mult
from machine_learning_hep.analysis.analyzer_jet import AnalyzerJet

# Collect all systematics
from machine_learning_hep.analysis.systematics import Systematics


from machine_learning_hep.io import parse_yaml

def create_dirs(common_name, n_dirs, parent_dir="./"):
    
    dir_list = [join(paren_dir, f"{common_name}_{i}") for i in range(n_dirs)]
    for d in dir_list:
        os.makedirs(d)
    return dir_list

def replace_all_paths(db_param, typean, parent_dir="./"):

    for datamc in ["data", "mc"]:
        for state in ["unmerged_tree_dir", "pkl", "pkl_skimmed", "pkl_skimmed_merge_for_ml"]:
            db_param["multi"][datamc][state] = \
                    create_dirs(state, len(db_param["multi"][datamc][state]), parent_dir)

        if datamc == "mc":
            state = "mcreweights"
            db_param["multi"][datamc][state] = \
                    create_dirs(state, len(db_param["multi"][datamc][state]), parent_dir)

        for state in ["pkl_skimmed_merge_for_ml_all", "pkl_eventcounter_all"]:
            db_param["multi"][datamc][state] = create_dirs(state, 1, parent_dir)[0]



        for state in ["pkl_skimmed_dec", "pkl_skimmed_decmerged"]:
            db_param["mlapplication"][datamc][state] = \
                    create_dirs(state, len(db_param["mlapplication"][datamc][state]), parent_dir)

        state = "results"
        db_param["analysis"][typean][datamc][state] = \
                create_dirs(state, len(db_param["analysis"][typean][datamc][state]), parent_dir)

        state = "resultsallp"
        db_param["analysis"][typean][datamc][state] = \
                create_dirs(state, 1, parent_dir)[0]

    db_param["ml"]["mlout"] = create_dirs("mlout", 1, parent_dir)[0]
    db_param["ml"]["mlplot"] = create_dirs("mlplot", 1, parent_dir)[0]

    db_param["analysis"]["inputfonllpred"] = create_dirs("inputfonllpred", 1, parent_dir)[0]
    db_param["analysis"]["dir_general_plots"] = create_dirs("dir_general_plots", 1, parent_dir)[0]


def dryrun(case, db_run_config, db_param, db_ml, db_ml_grid, db_run_list):

    typean = db_run_config["analysis"]["type"]

    # Change and construct all paths which might be picked up by any of the following objects
    # Jsut to make sure directories are there in case any of the objects relies on their existence
    replace_all_paths(db_param, typean, "OUTPUT")

    proc_type = db_param["analysis"][typean].get("proc_type", None)

    if not proc_type:
        print(f"ERROR: Required field \"proc_type\" not defined in analysis section {typean} for case {case}")
        exit(1)

    proc_class = None
    anal_class = None
    if proc_type == "Dhadrons":
        print("Using new feature for Dhadrons")
        proc_class = ProcesserDhadrons
        ana_class = AnalyzerDhadrons
    elif proc_type == "Dhadrons_mult":
        print("Using new feature for Dhadrons_mult")
        proc_class = ProcesserDhadrons_mult
        ana_class = AnalyzerDhadrons_mult
    elif proc_type == "Dhadrons_jet":
        print("Using new feature for Dhadrons_jet")
        proc_class = ProcesserDhadrons_jet
        ana_class = AnalyzerJet

    if not proc_class:
        print(f"INFO: Processing type {proc_type} not supported")
        exit(0)

    # Test creating processer
    multi_proc_mc = MultiProcesser(case, proc_class, db_param[case], typean, db_run_list, "mc")
    mymultiprocessdata = MultiProcesser(case, proc_class, db_param[case], typean, db_run_list, "data")

    # Test creating analyzer
    ana_mgr = AnalyzerManager(ana_class, db_param[case], case, typean, True)

    # Test creating systematic
    # For now don't test systematics
    #syst_mgr = AnalyzerManager(Systematics, db_param[case], case, typean, True, db_run_list)

    # Test creating optimiser
    binminarray = db_param[case]["ml"]["binmin"]
    binmaxarray = db_param[case]["ml"]["binmax"]
    mltype = db_param[case]["ml"]["mltype"]
    for i, (binmin, binmax) in enumerate(zip(binminarray, binmaxarray)):
        myopt = Optimiser(data_param[case], case, typean,
                          data_model[mltype], grid_param, binmin, binmax,
                          raahp[i], training_vars[i])


parser = argparse.ArgumentParser()
parser.add_argument("--run-config", "-r", dest="run_config",
                    help="the run configuration to be used", required=True)
parser.add_argument("--database-analysis", "-d", dest="database_analysis",
                    help="analysis database to be used", required=True)
parser.add_argument("--database-ml-models", dest="database_ml_models",
                    help="ml model database to be used", required=True)
parser.add_argument("--database-ml-gridsearch", dest="database_ml_gridsearch",
                    help="ml gridsearch database to be used", required=True)
parser.add_argument("--database-run-list", dest="database_run_list",
                    help="run list database to be used", required=True)
parser.add_argument("--analysis", "-a", dest="type_ana",
                    help="choose type of analysis", required=True)

args = parser.parse_args()

configure_logger(args.debug, args.log_file)

# Extract which database and run config to be used
pkg_data = "machine_learning_hep.data"
pkg_data_run_config = "machine_learning_hep.submission"
run_config = load_config(args.run_config, (pkg_data_run_config, "default_complete.yaml"))

# Load all config in dictionaries

# Since all are required arguments we are save to do that
db_run_config = parse_yaml(args.run_config)
db_param = parse_yaml(args.database_analysis)
db_ml = parse_yaml(args.database_ml_models)
db_ml_grid = parse_yaml(args.database_ml_gridsearch)
db_run_list = parse_yaml(args.database_run_list)

run_config["analysis"]["type"] = args.type_ana

# Extract the case
case = None
for k in data_param.keys():
    case = k
    break

dryrun(case, db_run_config, db_param, db_ml, db_ml_grid, db_run_list)
