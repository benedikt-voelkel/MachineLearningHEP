#import ROOT
import argparse
import os
import glob
from machine_learning_hep.logger import get_logger
from machine_learning_hep.io import parse_yaml
import uproot
import matplotlib.pyplot as plt
import numpy as np

def process(top_dir, file_signature, recursive, n_files, plot_config):
    top_dir = os.path.expanduser(top_dir)
    logger = get_logger()
    if os.path.isdir(top_dir) is not True:
        logger_string = f"Directory {top_dir} does not exist."
        logger.fatal(logger_string)
    
    file_signature = os.path.join(top_dir, "**/", file_signature)
    logger_string = f"Globbing files with matching {file_signature}"
    logger.info(logger_string)

    root_files = glob.glob(file_signature, recursive=recursive)
    if not root_files:
        logger_fatal("No ROOT files selected")

    if len(root_files) < n_files:
        logger_string = f"Although {n_files} were requested, could only glob {len(root_files)}"
        logger.warning(logger_string)
    elif len(root_files) > n_files and n_files > 0:
        # Skim list of ROOT files to requested number
        root_files = root_files[:n_files]

    for name, prop in plot_config.items():
        if "observable" not in prop:
            logger_string = f"Cannot find key \"observable\" for plot {name}"
            logger.fatal(logger_string)
        # ...
        to_plot = None
        hist, edges = np.histogram([], prop["bins"], prop["range"])
        for f in root_files:
            values_selected = None
            if prop["observable"]["tree"] != prop["cut"]["tree"]:
                match_values = uproot.open(f)[prop["cut"]["tree"]].pandas.df(branches=[*prop["cut"]["required_branches"], prop["match_branch"]]).query(prop["cut"]["expression"])[prop["match_branch"]]
                df_observable = uproot.open(f)[prop["observable"]["tree"]].pandas.df(branches=[prop["observable"]["name"], prop["match_branch"]])
                values_selected = df_observable.loc[df_observable[prop["match_branch"]].isin(match_values)][prop["observable"]["name"]]
            else:
                values_selected = uproot.open(f)[prop["cut"]["tree"]].pandas.df(branches=[ *prop["cut"]["required_branches"], prop["observable"]["name"] ]).query(prop["cut"]["expression"])[prop["observable"]["name"]]
            hist, edges = np.histogram([], prop["bins"], prop["range"])
            hist += np.histogram(values_selected, prop["bins"], prop["range"])[0]

        bin_centers = [ (edges[i+1] + edges[i]) / 2 for i in range(len(edges) - 1) ]
        plt.bar(bin_centers, hist)
        plt.show()


    exit(0)

    RDF = ROOT.ROOT.RDataFrame
    RCH = ROOT.TChain
    
    vector_root_files = ROOT.vector('string')()
    vector_branches = ROOT.vector('string')()
    for f in root_files:
        vector_root_files.push_back(f)
    branches = ["pt_cand", "evt_id", "run_number"]
    for b in branches:
        vector_branches.push_back(b)


    event_ids = RDF("PWGHF_TreeCreator/tree_Lc2V0bachelor", vector_root_files, vector_branches).Take("evt_id") #Filter("pt_cand > 1. && pt_cand < 3.").GetValue().Take("evt_id")
    print(event_ids)
    exit(0)

    chain_cand = RCH("PWGHF_TreeCreator/tree_Lc2V0bachelor")
    chain_event = RCH("PWGHF_TreeCreator/tree_event_char")
    for f in root_files:
        chain_cand.Add(f)
        chain_event.Add(f)

    chain_event.AddFriend(chain_cand, "friend")

    df = RDF(chain_event)

    #histo_low_pt = df.Filter("friend.pt_cand > 1. && friend.pt_cand < 3.").Histo1D("n_tracklets")
    #histo_high_pt = df.Filter("friend.pt_cand > 4. && friend.pt_cand < 6.").Histo1D("n_tracklets").GetValue()
    count = df.Filter("friend.pt_cand > 4. && friend.pt_cand < 6.").Count()

    print("%s entries passed all filters" %count.GetValue())

    #print(f"Entries: {count.GetValue()}")

    #histo_low_pt.Draw("hist")


parser = argparse.ArgumentParser()
parser.add_argument("top", help="directory where it should be looked for ROOT files")
parser.add_argument("--file-signature", dest="file_signature", help="regular expression matching filenames", default=r"*.root")
parser.add_argument("--recursive", action="store_true", help="recursive search for files")
parser.add_argument("--n-files", dest="n_files", help="number of files to be processed", type=int, default=-1)
parser.add_argument("--plot-config", dest="plot_config", help="plots to be produced", required=True)


args = parser.parse_args()

plot_config = parse_yaml(args.plot_config)

process(args.top, args.file_signature, args.recursive, args.n_files, plot_config)

