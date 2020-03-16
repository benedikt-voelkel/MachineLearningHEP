from os.path import join as osjoin
from os.path import exists as osexists
from os import makedirs
import argparse
from glob import glob
import pickle
import numpy as np
import pandas as pd
import lz4.frame
import matplotlib.pyplot as plt
from machine_learning_hep.utilities import openfile
from machine_learning_hep.io import parse_yaml
from machine_learning_hep.logger import get_logger
import multiprocessing as mp
import sys
from seaborn import heatmap
from math import sqrt


LOGGER = get_logger()

def callback(e):
    print(e)



def parallelizer(function, argument_list, kw_argument_list, maxperchunk):
        chunks_args = [argument_list[x:x+maxperchunk] \
                  for x in range(0, len(argument_list), maxperchunk)]
        if not kw_argument_list:
            kw_argument_list = [{} for _ in argument_list]
        chunks_kwargs = [kw_argument_list[x:x+maxperchunk] \
                  for x in range(0, len(kw_argument_list), maxperchunk)]
        res = None
        for chunk_args, chunk_kwargs in zip(chunks_args, chunks_kwargs):
            print("Processing new chunck size=", maxperchunk)
            pool = mp.Pool(10)
            #res = [pool.apply_async(function, args=args, kwds=kwds, error_callback=callback) for args, kwds in zip(chunk_args, chunk_kwargs)]
            res = [pool.apply_async(function, args=args, kwds=kwds, error_callback=callback) for args, kwds in zip(chunk_args, chunk_kwargs)]
            pool.close()
            pool.join()

        res_list = None
        try:
            res_list = [r.get() for r in res]
        except Exception as e:
            print("EXCEPTION")
            pass
        return res_list


        #return
def print_df_info(df):
    print("Dataframe info")

    print("Columns")
    for c in df.columns:
        print(f"\t{c}")

    print(f"Rows: {df.shape[0]}, Columns: {df.shape[1]}")



def save_fig(fig, path):
    """Wrapper to save figures

    Args:
        fig: figure (or anything that implements savefig)
        path: full path where to save
    """
    LOGGER.info("Save figure to %s", path)
    fig.savefig(path)


def make_2d_plot(df_input, columns, save_path, queries=None, **kwargs):
    """Make a 2d plot from 2 columns in a dataframe

    Args:
        df_pkl_path: path to pickled dataframe
        columns: tuple with 2 column names
        save_path: where to asve the plot
        queries: for each query 1 plot is added (optional)
        title: plot title (optional)
        **kwargs:
            title: str
            legend: tuple, legend title per query, if queries == None, str
    """
    
    LOGGER.info("Read dataframe from %s and save plot to %s", df_input, save_path)

    title = kwargs.pop("title", None)
    legend = kwargs.pop("legend", None)
    is_heatmap = kwargs.pop("is_heatmap", False)
    sub_query = kwargs.pop("sub_query", False)
    bins = kwargs.pop("bins", (200, 200))
    add_profile = kwargs.pop("add_profile", (False, False))

    df_all = None

    if isinstance(df_input, str):
        df_all = pickle.load(openfile(df_input, "rb"))
    else:
        df_all = df_input

    LOGGER.info("Found %i samples in the dataframe", df_all.shape[0])

    dfs = []
    if queries:
        if not legend:
            legend = [None] * len(queries)
        dfs = (df_all.query(query) for query in queries)
    else:
        dfs = (df_all,)
        legend = (legend,)
    fig, ax = plt.subplots(1, 1)

    has_legend = False
    for df, key in zip(dfs, legend):
        if is_heatmap:
            #data, xedges, yedges = np.histogram2d(df[columns[0]].values, df[columns[1]].values, bins=bins)
            data, xedges, yedges, image = ax.hist2d(df[columns[0]].values, df[columns[1]].values, bins=bins, cmap="coolwarm")
            plt.colorbar(image)

            y_bin_centers = [(yedges[i] + (yedges[i+1] - yedges[i]) / 2.) for i in range(0, len(yedges) - 1)]
            np_y_bin_centers = np.array(y_bin_centers)
            if add_profile[0]:
                # This means we do a profile plot
                profile_list_mean = []
                profile_list_std = []
                for d in data:
                    dnp = np.array(d)
                    if not np.any(dnp):
                        profile_list_mean.append(0)
                        profile_list_std.append(0)
                    else:
                        profile_list_mean.append(np.average(np_y_bin_centers, weights=dnp))
                        profile_list_std.append(sqrt(np.average((np_y_bin_centers - profile_list_mean[-1])**2, weights=dnp)))
                x_bin_centers = [(xedges[i] + (xedges[i+1] - xedges[i]) / 2.) for i in range(0, len(xedges) - 1)]
                ax.plot(x_bin_centers, profile_list_mean, color="black")
                ax.fill_between(x_bin_centers, [m - e for m, e in zip(profile_list_mean, profile_list_std)],
                                [m + e for m, e in zip(profile_list_mean, profile_list_std)], color="gray", alpha=0.8)

        else:
            ax.scatter(df[columns[0]].values, df[columns[1]].values, label=key if key else "", marker="x")
        has_legend = (has_legend or key)

    ax.set_xlabel(columns[0])
    ax.set_ylabel(columns[1])

    if has_legend and not is_heatmap:
        ax.legend(loc="best")

    if title:
        ax.set_title(title)

    fig.tight_layout()
    save_fig(fig, save_path)
    plt.close(fig)
    

def make_1d_plot(df_pkl_path, columns, save_path, queries=None, **kwargs):
    """Make a 2d plot from 2 columns in a dataframe

    Args:
        df_pkl_path: path to pickled dataframe
        columns: tuple with 2 column names
        save_path: where to asve the plot
        queries: for each query 1 plot is added (optional)
        title: plot title (optional)
        **kwargs:
            title: str
            legend: tuple, legend title per query, if queries == None, str
    """
    
    LOGGER.info("Read dataframe from %s and save plot to %s", df_pkl_path, save_path)

    titles = kwargs.pop("titles", None)
    legend = kwargs.pop("legend", None)


    if not titles:
        titles = [None] * len(columns)

    df_all = pickle.load(openfile(df_pkl_path, "rb"))
    print_df_info(df_all)


    LOGGER.info("Found %i samples in the dataframe", df_all.shape[0])

    dfs = []

    if queries:
        if not legend:
            legend = [None] * len(queries)
        dfs = [df_all.query(query) for query in queries]
    else:
        dfs = (df_all,)
        legend = (legend,)

    figsize = (30, 30)
    n_columns = 4
    n_rows = int(len(columns) / n_columns) + 1
    fig, axes_tmp = plt.subplots(n_rows, n_columns, figsize=figsize)


    # Flatten axes
    axes = []
    for ax_sub in axes_tmp:
        axes.extend(ax_sub)


    for ax, col, title in zip(axes, columns, titles):
        has_legend = False
        for df, key in zip(dfs, legend):
            ax.hist(df[col].values, label=key if key else "", bins=100, alpha=0.5)
            has_legend = (has_legend or key)

        ax.set_xlabel(col)
        ax.set_ylabel("counts")
        if has_legend:
            ax.legend(loc="best")

        if title:
            ax.set_title(title)

    save_fig(fig, save_path)
    plt.close(fig)



def skim_df(df_pkl_path, columns=None, query=None):
    """Skim dataframe

    Args:
        df_pkl_path: path to pickled dataframe
        columns: tuple with column names to extract
        save_path: where to asve the plot
        query: query the dataframe (optional)
    """
    
    #print("Read dataframe from %s", df_pkl_path)

    df = pickle.load(openfile(df_pkl_path, "rb"))

    #print(f"Found {df.shape[0]} samples in the dataframe")

    if query:
        df = df.query(query)

    if not columns:
        return df

    return df[columns]


#database = parse_yaml("/home/bvolkel/HF/MachineLearningHEP/machine_learning_hep/databases/database_ml_parameters_LcpK0spp_test.yml")
database = parse_yaml("/home/bvolkel/HF/MachineLearningHEP/machine_learning_hep/data/JetAnalysis/database_ml_parameters_LcpK0spp_20200301.yml")
database = database[list(database.keys())[0]]

mc_data = ("mc", "data")

ml_cut_presel = {}
ml_cut_sel = database["mlapplication"]["probcutoptimal"]
for case in mc_data:
    ml_cut_presel[case] = database["mlapplication"]["probcutpresel"][case]


pt_bins = list(zip(database["sel_skim_binmin"], database["sel_skim_binmax"]))
cut_MC_sig = database["ml"]["sel_sigml"]

training_variables = database["variables"]["var_training"]

file_name_base = "AnalysisResultsReco"
file_name_extension = ".pkl.lz4"
file_middle = database["var_binning"]


########################################
#
# V0M percentile
#
########################################

save_dir = "plots/v0m_perc"
if not osexists(save_dir):
    makedirs(save_dir)

query = "v0m_perc <= 100"
columns = ["perc_v0m", "v0m_eq_corr"]


for case in ("data",):
    for period, dir_applied in zip(database["multi"][case]["period"], database["multi"][case]["pkl_skimmed"]):
        args_list = []
        kwargs_list = []
        for pt_bin, training_vars in zip(pt_bins[4:5], training_variables[3:]):
            file_name = f"{file_name_base}_{file_middle}{pt_bin[0]}_{pt_bin[1]}{file_name_extension}"

            print(f"{dir_applied}/**/{file_name}")

            for f in  glob(f"{dir_applied}/**/{file_name}", recursive=True):
                #args_list.append((f, columns, query))
                args_list.append((f, columns))
                #kwargs_list.append({"query": query})


            dfs = parallelizer(skim_df, args_list, kwargs_list, 100)
            print(dfs)

            df_all = pd.concat(dfs)

            save_path = f"{file_name_base}{pt_bin[0]}_{pt_bin[1]}_{case}_{period}_v0m.png"
            save_path = osjoin(save_dir, save_path)
            make_2d_plot(df_all, columns, save_path, is_heatmap=True, add_profile=(True, False), bins=(100, 100))


exit(0)




########################################
#
# all training variables and a few more
#
########################################
"""
args_list = []
kwargs_list = []
# Need this as a numpy float32 since this is the type of nsigTOF and nsigTPC in the dataframe
# Otherwise the comparison fails
abs_max_tof = np.float32(10)
legend = (f"abs(nsigTOF) < {abs_max_tof} and {cut_MC_sig}",)
file_name_base = "AnalysisResultsReco"
file_name_extension = ".pkl.lz4"
file_middle = database["var_binning"]
case = "mc"


save_dir = "plots/MC_truth/TOF_smaller_10"
if not osexists(save_dir):
    makedirs(save_dir)

for case in ("mc",):
    for period, dir_applied in zip(database["multi"][case]["period"], database["multi"][case]["pkl_skimmed_merge_for_ml"]):
        for pt_bin, training_vars in zip(pt_bins, training_variables):
            file_name = f"{file_name_base}_{file_middle}{pt_bin[0]}_{pt_bin[1]}{file_name_extension}"
            file_name = osjoin(dir_applied, file_name)

            save_path = f"{file_name_base}{pt_bin[0]}_{pt_bin[1]}_{case}_{period}_TOF.png"
            save_path = osjoin(save_dir, save_path)

            cut = cut_MC_sig
            queries = (f"({cut}) & (abs(nsigTOF_Pr_0) < @abs_max_tof)",)

            plot_vars = training_vars + ["p_prong0", "p_prong1", "p_prong2", "pt_prong0", "pt_prong1", "pt_prong2", "phi_cand", "eta_cand", "pt_cand", "nsigTPC_Pr_0", "nsigTOF_Pr_0"]

            args_list.append((file_name, plot_vars, save_path, queries,))
            kwargs_list.append({"legend": legend, "titles": plot_vars})

parallelizer(make_1d_plot, args_list, kwargs_list, 100)

exit(0)
"""
####################################
#
# 
#
####################################


args_list = []
kwargs_list = []
mask_value = np.float32(-999)
abs_max_tof = np.float32(10)
var_names = ("pt_prong0", "nsigTOF_Pr_0")
case = "mc"
save_dir = "plots/MC_truth/pt_prong0_TOF"
if not osexists(save_dir):
    makedirs(save_dir)


for period, dir_applied in zip(database["multi"]["mc"]["period"], database["mlapplication"]["mc"]["pkl_skimmed_decmerged"]):
    for pt_bin, ml_cut_pre, ml_cut in zip(pt_bins, ml_cut_presel["mc"], ml_cut_sel):
        file_name = f"{file_name_base}{pt_bin[0]}_{pt_bin[1]}_{ml_cut_pre:.2f}{file_name_extension}"
        file_name = osjoin(dir_applied, file_name)

        save_path = f"{file_name_base}{pt_bin[0]}_{pt_bin[1]}_{ml_cut_pre:.2f}_{case}_{period}_truth_signal.png"
        #save_path = f"{file_name_base}{pt_bin[0]}_{pt_bin[1]}_{ml_cut_pre:.2f}_{case}_{period}_truth_tagged.png"
        save_path = osjoin(save_dir, save_path)

        #cut = f"y_test_probxgboost > {ml_cut}"
        cut = cut_MC_sig
        queries = (f"({cut}) & (abs(nsigTOF_Pr_0) < @abs_max_tof)",)

        n_bins = 100

        x_min = 0.
        x_max = 10.
        dist_x = (x_max - x_min) / n_bins
        y_min = -abs_max_tof
        y_max = abs_max_tof
        dist_y = (y_max - y_min) / n_bins

        bins = ([x_min + i * dist_x for i in range(n_bins)] + [x_max], [y_min + i * dist_y for i in range(n_bins)] + [y_max])

        args_list.append((file_name, var_names, save_path, queries,))
        kwargs_list.append({"title": f"{pt_bin[0]}_{pt_bin[1]}_{ml_cut_pre:.2f}_{case}_{period}, {cut}",
                            "is_heatmap": True,
                            "add_profile": (True, False),
                            "bins": (100, 100)})

parallelizer(make_2d_plot, args_list, kwargs_list, 100)
exit(0)

########################################
#
# all training variables and a few more
#
########################################

args_list = []
kwargs_list = []
# Need this as a numpy float32 since this is the type of nsigTOF and nsigTPC in the dataframe
# Otherwise the comparison fails
mask_value = np.float32(-999)
queries = ("nsigTOF_Pr_0 == @mask_value", "nsigTOF_Pr_0 != @mask_value", f"(nsigTOF_Pr_0 == @mask_value) & ({cut_MC_sig})", f"(nsigTOF_Pr_0 != @mask_value) & ({cut_MC_sig})",)
legend = ("nsigTOF == -999", "nsigTOF != -999", f"nsigTOF == -999 and {cut_MC_sig}", f"nsigTOF != -999 {cut_MC_sig}",)
file_name_base = "AnalysisResultsReco"
file_name_extension = ".pkl.lz4"
file_middle = database["var_binning"]
case = "mc"



save_dir = "plots/MC_truth"
if not osexists(save_dir):
    makedirs(save_dir)

for case in ("mc",):
    for period, dir_applied in zip(database["multi"][case]["period"], database["multi"][case]["pkl_skimmed_merge_for_ml"]):
        for pt_bin, training_vars in zip(pt_bins, training_variables):
            file_name = f"{file_name_base}_{file_middle}{pt_bin[0]}_{pt_bin[1]}{file_name_extension}"
            file_name = osjoin(dir_applied, file_name)

            save_path = f"{file_name_base}{pt_bin[0]}_{pt_bin[1]}_{case}_{period}_TOF.png"
            save_path = osjoin(save_dir, save_path)

            plot_vars = training_vars + ["p_prong0", "p_prong1", "p_prong2", "pt_prong0", "pt_prong1", "pt_prong2", "phi_cand", "eta_cand", "pt_cand"]

            args_list.append((file_name, plot_vars, save_path, queries,))
            kwargs_list.append({"legend": legend, "titles": plot_vars})

parallelizer(make_1d_plot, args_list, kwargs_list, 100)

exit(0)
