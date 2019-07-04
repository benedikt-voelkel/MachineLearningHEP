#import ROOT
import argparse
import os
import glob
from machine_learning_hep.logger import configure_logger, get_logger
from machine_learning_hep.io import parse_yaml
import uproot
import matplotlib.pyplot as plt
import numpy as np

class Histo1D:
    def __init__(self, name, bins, axis_range, label):
        self.name = name
        self.bins = bins
        self.axis_range = axis_range
        self.hist, self.edges = np.histogram([], bins, axis_range)
        self.label = label

    def add_values(self, values, weights=None):
        self.hist += np.histogram(values, self.bins, self.axis_range, None, weights)[0]


def plot1D(histograms, plot_name, output_dir="./"):

        logger = get_logger()
        try:
            it = iter(histograms)
        except TypeError as te:
            histograms = [histograms]

        output_dir = os.path.expanduser(output_dir)
        fig = plt.figure()
        ax = fig.add_subplot()
        # Derive bin centers from edges of first histogram
        common_edges = histograms[0].edges
        bin_widths = [ 0.9 * (common_edges[i+1] - common_edges[i] ) for i in range(len(common_edges) - 1) ]
        bin_centers = [ (common_edges[i+1] + common_edges[i]) / 2 for i in range(len(common_edges) - 1) ]
        for h in histograms:
            if not np.array_equal(h.edges,common_edges):
                logger_string = f"Incompatible edges found in histogram {h.name}"
                logger.fatal(logger_string)
            ax.bar(bin_centers, h.hist, width=bin_widths, alpha=0.7, label=h.label)
        output_png = os.path.join(output_dir, plot_name + ".png")
        #output_eps = os.path.join(output_dir, plot_name + ".eps")
        logger_string = f"Save histogram {plot_name} in directory {output_dir}"
        logger.info(logger_string)
        plt.savefig(output_png)
        #plt.savefig(output_eps)


def process(top_dir, file_signature, recursive, n_files, plot_config, output_dir):
    top_dir = os.path.expanduser(top_dir)
    logger = get_logger()
    if not os.path.isdir(top_dir):
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

    available_cuts = plot_config["cuts"]
    available_cuts[None] = None
    plots = plot_config["plots"]
    tree_matches = plot_config["tree_matches"]

    """
    required_cuts = []
    # Collect all required cuts
    for p in plots:
        for o in p["observables"]:
            o["cuts"] = o.get("cuts",[None])
            if o["cuts"][0] is None and [None, None] not in required_cuts:
                required_cuts.append([None,None])
                continue
            for c in o["cuts"]:
                if c not in required_cuts:
                    required_cuts.append(c)
    # Check if chosen cuts are available
    # To match more than one match branch name to one cut name
    required_cuts_dict = {}
    for c in required_cuts:
        if c[0] not in available_cuts:
            logger_string = f"Cut {c[0]} is not available"
            logger.fatal(logger_string)
        if c[0] not in required_cuts_dict:
            required_cuts_dict[c[0]] = []
        if c[0] is None:
            required_cuts_dict[None] = [None]
        elif c[1] is not None:
            required_cuts_dict[c[0]].append(c[1])

    print(required_cuts_dict)
    """

    histograms = {}
    # Now, go though all files
    for i_file, f in enumerate(root_files):
        logger_string = f"Process {i_file}. file at {f}"
        logger.debug(logger_string)
        # Brute force for now
        for plot in plots:
            for o in plot["observables"]:
                for cut in o["cuts"]:
                    if cut is None:
                        values = uproot.open(f)[plot["tree"]].pandas.df(branches=[o["name"]])[o["name"]]
                        histo_name = o["name"] + "_no_cut"
                        if o["name"] not in histograms:
                            histograms[o["name"]] = {}
                        if histo_name not in histograms[o["name"]]:
                                histograms[o["name"]][histo_name] = Histo1D(histo_name, o["bins"], o["range"], "no_cut")
                        histograms[o["name"]][histo_name].add_values(values)
                    elif plot["tree"] == available_cuts[cut]["tree"]:
                        values = uproot.open(f)[plot["tree"]].pandas.df(branches=[o["name"], *available_cuts[cut]["required_branches"]]).query(available_cuts[cut]["expression"])[o["name"]]
                        histo_name = o["name"] + "_" + cut
                        if o["name"] not in histograms:
                            histograms[o["name"]] = {}
                        if histo_name not in histograms[o["name"]]:
                                histograms[o["name"]][histo_name] = Histo1D(histo_name, o["bins"], o["range"], cut)
                        histograms[o["name"]][histo_name].add_values(values)
                    else:
                        match_branch = None
                        for tree_match in tree_matches:
                            if plot["tree"] in tree_match["trees"] and available_cuts[cut]["tree"] in tree_match["trees"]:
                                match_branch = tree_match["common_branch"]
                                break
                        if match_branch is None:
                            logger_srting = f"Could not find match_branch for trees {plot['tree']} and {available_cuts[cut]['tree']}"
                            logger.fatal(logger_string)
                        match_values = uproot.open(f)[available_cuts[cut]["tree"]].pandas.df(branches=[match_branch, *available_cuts[cut]["required_branches"]]).query(available_cuts[cut]["expression"])[match_branch]
                        df_observable = uproot.open(f)[plot["tree"]].pandas.df(branches=[o["name"], match_branch])
                        #values = df_observable.loc[df_observable[prop["match_branch"]].isin(match_values)][prop["observable"]["name"]]
                        values = df_observable.loc[df_observable[match_branch].isin(match_values)][o["name"]]
                        histo_name = o["name"] + "_" + cut
                        if o["name"] not in histograms:
                            histograms[o["name"]] = {}
                        if histo_name not in histograms[o["name"]]:
                                histograms[o["name"]][histo_name] = Histo1D(histo_name, o["bins"], o["range"], cut)
                        histograms[o["name"]][histo_name].add_values(values)

    for h_name, histos in histograms.items():
        histo_list = [ v for v in histos.values() ]
        plot1D(histo_list, h_name)
    exit(0)





"""

    # Now, go though all files
    for f in [root_files[0]]:
        # Now got through each cut and fill histograms
        for cut_name, match_branches in required_cuts_dict.items():
            if cut_name is None:
                for plot in plots:
                    branches = []
                    for o in plot["observables"]:
                        if [None,None] in o["cuts"]:
                            branches.append(o["name"])
                    df = uproot.open(f)[plot["tree"]].pandas.df(branches=[branches])
                    for o in plot["observables"]:
                        if [None, None] not in o["cuts"]:
                            continue
                        histo_name = o["name"] + "_no_cut"
                        if histo_name not in histograms:
                            histograms[histo_name] = Histo1D(histo_name, o["bins"], o["range"])
                        histograms[histo_name].add_values(df[o["name"]])
                # Cut name None end
                continue
            # Collect plot indices
            plot_indices = []
            # Extract further branches if plot has same source tree as cut tree
            cut_required_branches = available_cuts[cut_name]["required_branches"]
            match_branches = []
            variables_in_cut_tree_branches = []
            # Now get all histograms requiring this cut and collect match branch names if any
            for plot in plots:
                if plot["tree"] == available_cuts[cut_name]["tree"]:
                    for o in plot["observables"]:
                        variables_in_cut_tree_branches.append(o["name"])
                elif plot["tree"] != available_cuts[cut_name]["tree"]:
                    for o in plot["observables"]:
                        for c in o["cuts"]:
                            if c[0] == cut_name and c[1] not in branches:
                                match_branches.append(c[1])
            branches = cut_required_branches + match_branches + variables_in_cut_tree_branches
            df_cut = uproot.open(f)[available_cuts[cut_name]["tree"]].pandas.df(branches=branches).query(available_cuts[cut_name]["expression"])
            for plot in plots:
                if plot["tree"] == available_cuts[cut_name]["tree"]:
                    for o in plot["observables"]:
                        histo_name = o["name"] + "_" + cut_name
                        for c in o["cuts"]:
                            if c[0] == cut_name:
                                if histo_name not in histograms:
                                    histograms[histo_name] = Histo1D(histo_name, o["bins"], o["range"])
                                histograms[histo_name].add_values(df_cut[o["name"]])
                else:
                    df_other = 
                    for o in plot["observables"]:
                        histo_name = o["name"] + "_" + cut_name
                        for c in o["cuts"]:
                            if c[0] == cut_name:
                                if histo_name not in histograms:
                                    histograms[histo_name] = Histo1D(histo_name, o["bins"], o["range"])
                                histograms[histo_name].add_values(df_cut[o["name"]])


                                    
                    



            print(cut_name)
            print(plot_indices)
            print(branches_in_cut_tree)
            print("---")
            


                # We can skip this if there is either
                # 1) No branch we need to match with, so same tree
                if "cut_match_variables" not in p or "cuts" not in p or c not in p["cuts"]:
                    continue
                if c in p["cuts"]:
                    # If this histogram requires the cut check if a match_variable needs to be extracted
                    index = p["cuts"].index(c)
                    cut_match_variables.append(p["cut_match_variables"][index])
            # Now there are all variables to be extracted form the cut tree


                





    # Prepare required histograms
    histograms = {}
    for p in plots:
        cuts = p.get




    for name, prop in plot_config.items():
        if "observable" not in prop:
            logger_string = f"Cannot find key \"observable\" for plot {name}"
            logger.fatal(logger_string)
        # ...
        to_plot = None
        histograms = {}
        edges = {}
        for o in prop["observables"]:
            
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
"""




parser = argparse.ArgumentParser()
parser.add_argument("top", help="directory where it should be looked for ROOT files")
parser.add_argument("--file-signature", dest="file_signature", help="regular expression matching filenames", default=r"*.root")
parser.add_argument("--recursive", action="store_true", help="recursive search for files")
parser.add_argument("--debug", action="store_true", help="activate debug log level")
parser.add_argument("--n-files", dest="n_files", help="number of files to be processed", type=int, default=-1)
parser.add_argument("--plot-config", dest="plot_config", help="plots to be produced", required=True)
parser.add_argument("--output-dir", dest="output_dir", help="directory where plots are dumped", default="./plots_output/")


args = parser.parse_args()

configure_logger(args.debug)
plot_config = parse_yaml(args.plot_config)

process(args.top, args.file_signature, args.recursive, args.n_files, plot_config, args.output_dir)

