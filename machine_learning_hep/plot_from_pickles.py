import glob
import pickle
import numpy as np
import pandas as pd
import lz4.frame
from machine_learning_hep.utilities import openfile
from ROOT import TFile, TH1D, TCanvas, gROOT
from root_numpy import fill_hist
import multiprocessing as mp

var_name = "delta_mass_KK"

def fill(f, c, h):
    print(f"Process file {f}")
    print(f"Prob cut {c}")
    df = pickle.load(openfile(f, "rb")).query(f"y_test_probxgboost>{c}")
    #df = pickle.load(openfile(f, "rb")) #.query(f"y_test_probxgboost>{c}")
    df = df.query("ismcprompt==1")
    fill_hist(h, df[var_name])
    return h


def parallelizer(function, argument_list, maxperchunk):
        chunks = [argument_list[x:x+maxperchunk] \
                  for x in range(0, len(argument_list), maxperchunk)]
        res = None
        for chunk in chunks:
            print("Processing new chunck size=", maxperchunk)
            pool = mp.Pool(10)
            res = [pool.apply_async(function, args=chunk[i]) for i in range(len(chunk))]
            pool.close()
            pool.join()
        return [r.get() for r in res]


gROOT.SetBatch(True)

"""
files = ["/data/Derived/DskINT7HighMultwithJets/vAN-20191003_ROOT6-1_opt/pp_2016_data/258_20191004-0007/skpkldecmerged/*.lz4",
         "/data/Derived/DskINT7HighMultwithJets/vAN-20191003_ROOT6-1_opt/pp_2017_data/259_20191004-0008/skpkldecmerged/*.lz4",
         "/data/Derived/DskINT7HighMultwithJets/vAN-20191003_ROOT6-1_opt/pp_2018_data/260_20191004-0008/skpkldecmerged/*.lz4"]
"""
files = ["/data/Derived/DskINT7HighMultwithJets/vAN-20191003_ROOT6-1_opt/pp_2016_mc_prodDs/261_20191004-0007/skpkldecmerged/AnalysisResultsReco*.lz4",
         "/data/Derived/DskINT7HighMultwithJets/vAN-20191003_ROOT6-1_opt/pp_2017_mc_prodDs/261_20191004-0007/skpkldecmerged/AnalysisResultsReco*.lz4",
         "/data/Derived/DskINT7HighMultwithJets/vAN-20191003_ROOT6-1_opt/pp_2018_mc_prodDs/261_20191004-0007/skpkldecmerged/AnalysisResultsReco*.lz4"]
input_files = []

cut = [0.82,0.92,0.9,0.6,0.6]
all_cuts = []

for f in files:
    input_files = input_files + glob.glob(f)
    all_cuts = all_cuts + cut

print(input_files)
histos = []
for i in range(len(input_files)):
    histos.append(TH1D(f"{var_name}_{i}", var_name, 100, 0., 0.03))


args = []

for f, c, h in zip(input_files, all_cuts, histos):
    args.append((f, c, h,))

res = parallelizer(fill, args, 100)

histo = TH1D(var_name, var_name, 100, 0., 0.03)

for h in res:
    histo.Add(h)

histo.GetXaxis().SetTitle(var_name)
histo.GetYaxis().SetTitle("# entries")
c = TCanvas("name", "", 800, 800)
histo.Draw()
c.SaveAs("delta_mass_KK_MLopt.eps")
#c.SaveAs("delta_mass_KK_MLpresel.eps")
c.Close()

root_file = TFile.Open("delta_mass_KK_MLopt_MC.root", "recreate")
#root_file = TFile.Open("delta_mass_KK_MLpresel_MC.root", "recreate")
histo.Write()
root_file.Close()

