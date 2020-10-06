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

"""
Main script for doing the systematic studies. Standalone, so some parts similar to analyzer.py

At the moment includes: Cut variation and MC pT shape
The raw yield systematic is done within analyzer.py
"""
# pylint: disable=unused-wildcard-import, wildcard-import
# pylint: disable=no-name-in-module
# pylint: disable=import-error
from os.path import join, exists
from os import makedirs
import math
from array import *
import pickle
from copy import deepcopy

from ROOT import gROOT, gPad
from ROOT import TFile, TH1F, TCanvas, TLegend
from ROOT import kRed, kGreen, kBlack, kBlue, kOrange, kViolet, kAzure, kYellow
from ROOT import Double

from root_numpy import fill_hist

#from machine_learning_hep.utilities import selectdfrunlist
from machine_learning_hep.utilities import seldf_singlevar, openfile, make_file_path
from machine_learning_hep.utilities_plot import load_root_style_simple, load_root_style
from machine_learning_hep.utilities import mergerootfiles, get_timestamp_string
from machine_learning_hep.analysis.analyzer import Analyzer, AnalyzerAfterBurner

class SystematicsAfterBurner(AnalyzerAfterBurner):
    # pylint: disable=useless-super-delegation
    def __init__(self, database, case, typean):
        super().__init__(database, case, typean)

    def ml_cutvar_mass(self):

        tmp_merged = f"/data/tmp/hadd/{self.case}_{self.typean}/cutvar_mass/" \
                     f"{get_timestamp_string()}/"
        # This warning appears since it is set to None in the AnalyzerAfterBurner constructor
        # But it will have a meaningful value at this point
        # pylint: disable=not-an-iterable
        files_mass_cutvar = [syst.n_filemass_cutvar for syst in self.analyzers]

        filename_mass = self.datap["files_names"]["histofilename"].replace(".root", "_cutvar.root")
        merged_file = join(self.datap["analysis"][self.typean]["data"]["resultsallp"], "cutvar")
        merged_file = join(merged_file, filename_mass)

        mergerootfiles(files_mass_cutvar, merged_file, tmp_merged)

        # pylint: disable=fixme
        # Get the cut limits of the latest period and set them for the merged analyzer
        # FIXME This is not the entirely correct way as it does not correctly reflect
        #       the efficiency boundaries for the merged case
        # This warning appears since it is set to None in the AnalyzerAfterBurner constructor
        # But it will have a meaningful value at this point
        # pylint: disable=unsubscriptable-object
        self.analyzer_merged.min_cv_cut = self.analyzers[-1].min_cv_cut
        self.analyzer_merged.max_cv_cut = self.analyzers[-1].max_cv_cut

    def ml_cutvar_eff(self):

        tmp_merged = f"/data/tmp/hadd/{self.case}_{self.typean}/cutvar_eff/" \
                     f"{get_timestamp_string()}/"
        # This warning appears since it is set to None in the AnalyzerAfterBurner constructor
        # But it will have a meaningful value at this point
        # pylint: disable=not-an-iterable
        files_eff_cutvar = [syst.n_fileeff_cutvar for syst in self.analyzers]

        filename_eff = self.datap["files_names"]["efffilename"].replace(".root", "_cutvar.root")
        merged_file = join(self.datap["analysis"][self.typean]["data"]["resultsallp"], "cutvar")
        merged_file = join(merged_file, filename_eff)

        mergerootfiles(files_eff_cutvar, merged_file, tmp_merged)

    def mcptshape(self):

        tmp_merged = f"/data/tmp/hadd/{self.case}_{self.typean}/mcptshape_eff/" \
                     f"{get_timestamp_string()}/"
        # This warning appears since it is set to None in the AnalyzerAfterBurner constructor
        # But it will have a meaningful value at this point
        # pylint: disable=not-an-iterable
        files_eff_cutvar = [syst.n_fileeff_ptshape for syst in self.analyzers]

        filename_eff = self.datap["files_names"]["efffilename"].replace(".root", "_ptshape.root")
        merged_file = join(self.datap["analysis"][self.typean]["data"]["resultsallp"],
                           filename_eff)

        mergerootfiles(files_eff_cutvar, merged_file, tmp_merged)


# pylint: disable=too-many-lines
# pylint: disable=too-many-instance-attributes, too-many-statements, too-many-arguments
# pylint: disable=too-many-branches, too-many-nested-blocks
class Systematics(Analyzer):
    species = "systematics"

    cent_cv_cut = None
    min_cv_cut = None
    max_cv_cut = None

    def __init__(self, datap, case, typean, period, run_param, analyzer_merged, processer_mc, processer_data, analyzer):
        super().__init__(datap, case, typean, period)

        self.analyzer_merged = analyzer_merged
        self.nominal_processer_mc = processer_mc
        self.nominal_processer_data = processer_data
        self.nominal_analyzer = analyzer

        self.run_param = run_param
        self.p_period = datap["multi"]["data"]["period"][period] if period is not None \
                else "merged"
        self.d_pkl_decmerged_mc = datap["mlapplication"]["mc"]["pkl_skimmed_decmerged"][period] \
                if period is not None else ""
        self.d_pkl_decmerged_data = \
                datap["mlapplication"]["data"]["pkl_skimmed_decmerged"][period] \
                if period is not None else ""
        self.d_results = datap["analysis"][typean]["data"]["results"][period] \
                if period is not None else datap["analysis"][typean]["data"]["resultsallp"]

        self.v_var_binning = datap["var_binning"]
        #Binning used when skimming/training
        self.lpt_anbinmin = datap["sel_skim_binmin"]
        self.lpt_anbinmax = datap["sel_skim_binmax"]
        self.p_nptbins = len(datap["sel_skim_binmax"])
        #Analysis binning
        self.lpt_finbinmin = datap["analysis"][self.typean]["sel_an_binmin"]
        self.lpt_finbinmax = datap["analysis"][self.typean]["sel_an_binmax"]
        self.p_nptfinbins = len(self.lpt_finbinmin)
        self.bin_matching = datap["analysis"][self.typean]["binning_matching"]
        #Second analysis binning
        self.lvar2_binmin = datap["analysis"][self.typean]["sel_binmin2"]
        self.lvar2_binmax = datap["analysis"][self.typean]["sel_binmax2"]
        self.v_var2_binning = datap["analysis"][self.typean]["var_binning2"]
        self.v_var2_binning_gen = datap["analysis"][self.typean]["var_binning2_gen"]

        #ML model variables
        self.p_modelname = datap["mlapplication"]["modelname"]
        self.lpt_probcutfin = datap["mlapplication"]["probcutoptimal"]
        self.lpt_probcutpre_mc = datap["mlapplication"]["probcutpresel"]["mc"]
        self.lpt_probcutpre_data = datap["mlapplication"]["probcutpresel"]["data"]

        #Extra pre-selections
        self.triggerbit = datap["analysis"][self.typean]["triggerbit"]
        self.s_evtsel = datap["analysis"][self.typean]["evtsel"]
        self.s_presel_gen_eff = datap["analysis"][self.typean]["presel_gen_eff"]
        self.s_trigger_mc = datap["analysis"][self.typean]["triggersel"]["mc"]
        self.s_trigger_data = datap["analysis"][self.typean]["triggersel"]["data"]
        self.apply_weights = \
                datap["analysis"][self.typean]["triggersel"].get("usetriggcorrfunc", None) \
                is not None

        #Build names for input pickle files (data, mc_reco, mc_gen)
        self.n_reco = datap["files_names"]["namefile_reco"]
        self.n_gen = datap["files_names"]["namefile_gen"]
        self.lpt_recodec_mc = [self.n_reco.replace(".pkl", "%d_%d_%.2f.pkl" % \
                              (self.lpt_anbinmin[i], self.lpt_anbinmax[i], \
                               self.lpt_probcutpre_mc[i])) for i in range(self.p_nptbins)]
        self.lpt_recodec_data = [self.n_reco.replace(".pkl", "%d_%d_%.2f.pkl" % \
                                (self.lpt_anbinmin[i], self.lpt_anbinmax[i], \
                                 self.lpt_probcutpre_data[i])) for i in range(self.p_nptbins)]
        self.lpt_recodecmerged_mc = [join(self.d_pkl_decmerged_mc, self.lpt_recodec_mc[ipt])
                                     for ipt in range(self.p_nptbins)]
        self.lpt_recodecmerged_data = [join(self.d_pkl_decmerged_data, \
                                       self.lpt_recodec_data[ipt]) for ipt in range(self.p_nptbins)]

        self.lpt_gensk = [self.n_gen.replace(".pkl", "_%s%d_%d.pkl" % \
                          (self.v_var_binning, self.lpt_anbinmin[i], self.lpt_anbinmax[i])) \
                          for i in range(self.p_nptbins)]
        self.lpt_gendecmerged = [join(self.d_pkl_decmerged_mc, self.lpt_gensk[ipt]) \
                                 for ipt in range(self.p_nptbins)]

        #Build names for intermediate output ROOT files
        self.n_filemass = datap["files_names"]["histofilename"]
        self.n_fileeff = datap["files_names"]["efffilename"]
        self.n_filemass_cutvar = self.n_filemass.replace(".root", "_cutvar.root")
        self.n_fileeff_cutvar = self.n_fileeff.replace(".root", "_cutvar.root")
        self.n_fileeff_ptshape = self.n_fileeff.replace(".root", "_ptshape.root")
        #Build directories for intermediate output ROOT files
        self.d_results_cv = join(self.d_results, "cutvar")
        if not exists(self.d_results_cv):
            makedirs(self.d_results_cv)
        self.n_filemass_cutvar = join(self.d_results_cv, self.n_filemass_cutvar)
        self.n_fileeff_cutvar = join(self.d_results_cv, self.n_fileeff_cutvar)
        self.n_fileeff_ptshape = join(self.d_results, self.n_fileeff_ptshape)
        #Final file names for analyzer.py
        self.yields_filename_std = "yields"
        self.efficiency_filename_std = "efficiencies"
        self.cross_filename_std = "finalcross"
        #Final file names used for systematics
        self.yields_filename = "yields_cutvar"
        self.efficiency_filename = "efficiencies_cutvar"
        self.efficiency_filename_pt = "efficiencies_mcptshape"
        self.cross_filename = "finalcross_cutvar"
        self.ptspectra_filename = "ptspectra_for_weights"

        #Variables for cross section/corrected yield calculation (NB: not all corr are applied)
        self.f_evtnorm = join(self.d_results, "correctionsweights.root")
        self.p_indexhpt = datap["analysis"]["indexhptspectrum"]
        self.p_fd_method = datap["analysis"]["fd_method"]
        self.p_cctype = datap["analysis"]["cctype"]
        self.p_sigmav0 = datap["analysis"]["sigmav0"]
        self.p_bineff = datap["analysis"][self.typean]["usesinglebineff"]
        self.p_fprompt_from_mb = datap["analysis"][self.typean]["fprompt_from_mb"]
        self.p_triggereff = datap["analysis"][self.typean].get("triggereff", [1] * 10)
        self.p_triggereffunc = datap["analysis"][self.typean].get("triggereffunc", [0] * 10)
        self.p_inputfonllpred = datap["analysis"]["inputfonllpred"]

        #Variables for the systematic variations
        self.p_cutvar_minrange = datap["systematics"]["probvariation"]["cutvarminrange"]
        self.p_cutvar_maxrange = datap["systematics"]["probvariation"]["cutvarmaxrange"]
        self.p_ncutvar = datap["systematics"]["probvariation"]["ncutvar"]
        self.p_maxperccutvar = datap["systematics"]["probvariation"]["maxperccutvar"]
        self.p_fixedmean = datap["systematics"]["probvariation"]["fixedmean"]
        self.p_fixedsigma = datap["systematics"]["probvariation"]["fixedsigma"]
        self.p_weights = datap["systematics"]["mcptshape"]["weights"]
        self.p_weights_min_pt = datap["systematics"]["mcptshape"]["weights_min_pt"]
        self.p_weights_max_pt = datap["systematics"]["mcptshape"]["weights_max_pt"]
        self.p_weights_bins = datap["systematics"]["mcptshape"]["weights_bins"]
        # Require a minimum significance or a maximum chi2 for individual fits
        self.min_signif_fit = datap["systematics"]["probvariation"].get("min_signif_fit", -1.)
        self.max_red_chi2_fit = datap["systematics"]["probvariation"].get("max_red_chi2_fit", -1.)


        # Flag whether some internals methods have been executed
        self.done_mass = False
        self.done_eff = False
        self.done_fit = False

    def get_after_burner(self):
        return SystematicsAfterBurner(self.datap, self.case, self.typean, Systematic.)


    def define_cutvariation_limits(self):
        """
        Cut Variation: Defines the probability cuts based on a max percentage variation (set in DB)
        Produces N (set in DB) tighter and N looser probability cuts
        """

        # Only do this once for all
        if Systematics.cent_cv_cuts is not None:
            return

        self.logger.info("Defining systematic cut variations for period: %s", \
                         self.p_period)


        results_dirs_periods = [join(d, "tmp_ml_wp_limits") for d in datap["analyis"][self.typean]["mc"]["results"]]
        results_dir_all = join(datap["analyis"][self.typean]["mc"]["resultsallp"], "tmp_ml_wp_limits")

        # use multiprocesser here, prepare database
        datap = deepcopy(self.datap)

        datap["analysis"][self.typean]["mc"]["results"] = results_dirs_periods
        datap["analysis"][self.typean]["mc"]["resultsallp"] = results_dir_all

        for rdp in results_dirs_periods:
            if exists(rdp):
                continue
            makedirs(rdp)
        if not exists(results_dir_all):
            makedirs(results_dir_all)

        # MultiProcesser to cover all at once
        multi_processer_effs = MultiProcesser(self.case, self.nominal_processer_mc.__class__, datap,
                                              self.typean, self.nominal_processer_mc.run_param,
                                              "mc")

        # construct analyzer for all periods merged and use it for finding ML WP boundaries
        analyzer_effs = self.nominal_analyzer_merged.__class__(datap, self.case, self.typean, None)

        n_pt_bins = self.nominal_processer_mc.p_nptbins

        Systematics.cent_cv_cut = self.nominal_processer_mc.lpt_probcutfin

        ana_n_first_binning = analyzer_effs.p_nptbins
        ana_n_second_binning = analyzer_effs.p_nbin2

        bin_matching = self.nominal_analyzer_merged.bin_matching

        nominal_effs = [[None] * ana_n_first_binning for _ in ana_n_second_binning]
        for ibin1 in ana_n_first_binning:
            for ibin2 in ana_n_bins_binning:
                nominal_effs[ibin2][ibin1] = self.nominal_analyzer_merged.get_efficiency(ibin1, ibi2)

        # TODO This has to become a list of lists (mult and pT)
        # TODO extract central efficeincies
        Systematics.min_cv_cut = [[None] * ana_n_first_binning for _ in ana_n_second_binning]
        Systematics.max_cv_cut = [[None] * ana_n_first_binning for _ in ana_n_second_binning]

        ncutvar_temp = self.p_ncutvar * 2

        stepsmin = []
        stepsmax = []
        for ipt in range(n_pt_bins):

            stepsmin.append( \
              (self.lpt_probcutfin[ipt] - self.p_cutvar_minrange[ipt]) / ncutvar_temp)

            stepsmax.append( \
              (self.p_cutvar_maxrange[ipt] - self.lpt_probcutfin[ipt]) / ncutvar_temp)

        boundaries_found = [[0] * ana_n_first_binning for _ in ana_n_second_binning]
        for icv in range(ncutvar_temp):

            if sum([sum(bound) for bound in boundaries_found]) == ana_n_first_binning * ana_n_second_binning:
                break

            wps = [self.lpt_probcutfin[ipt] + stepsmin[ipt] for ipt in range(n_pt_bins)]
            wps_strings = ["y_test_prob%s>%s" % (self.p_modelname, wps[ipt]) \
                    for ipt in range(n_pt_bins)]
            # update processers and analyzer
            for proc in multi_processer_effs.process_listsample:
                proc.l_selml = wps_string
            analyzer_effs.lpt_probcutfin = wps

            multi_processer_effs.multi_efficiency()
            analyzer_effs.efficiency()

            # TODO I am here...
            for ibin1 in range(ana_n_first_binning):
                for ibin2 in range(ana_n_second_binning):
                    eff_new = analyzer_effs.get_efficiency(ibin1, ibin2)
                    if eff_new / nominal_effs[ibin2][ibin1] < 1 + self.p_maxperccutvar:
                        boundaries_found[ibin2][ibin1] = 1
                        continue
                    self.min_cv_cut[ibin2][ibin1] = wps[bin_matching[ibin1]]


        boundaries_found = [[0] * ana_n_first_binning for _ in ana_n_second_binning]
        for icv in range(ncutvar_temp):

            if sum([sum(bound) for bound in boundaries_found]) == ana_n_first_binning * ana_n_second_binning:
                break

            wps = [self.lpt_probcutfin[ipt] - stepsmax[ipt] for ipt in range(n_pt_bins)]
            wps_strings = ["y_test_prob%s>%s" % (self.p_modelname, wps[ipt]) \
                    for ipt in range(n_pt_bins)]
            # update processers and analyzer
            for proc in multi_processer_effs.process_listsample:
                proc.l_selml = wps_string
            analyzer_effs.lpt_probcutfin = wps

            multi_processer_effs.multi_efficiency()
            analyzer_effs.efficiency()

            # TODO I am here...
            for ibin1 in range(ana_n_first_binning):
                for ibin2 in range(ana_n_second_binning):
                    eff_new = analyzer_effs.get_efficiency(ibin1, ibin2)
                    if eff_new / nominal_effs[ibin2][ibin1] < 1 - self.p_maxperccutvar:
                        boundaries_found[ibin2][ibin1] = 1
                        continue
                    self.max_cv_cut[ibin2][ibin1] = wps[bin_matching[ibin1]]
        
        print("Limits for cut variation defined, based on eff %-var of: ", self.p_maxperccutvar)
        print("--Cut variation minimum: ", self.min_cv_cut)
        print("--Central probability cut: ", cent_cv_cut)
        print("--Cut variation maximum: ", self.max_cv_cut)


    def ml_cutvar_mass(self):
        """
        Cut Variation: Create ROOT file with mass histograms
        Histogram for each variation, for each pT bin, for each 2nd binning bin

        Similar as process_histomass_single(self, index) in processor.py
        """

        # Do this only for a single period
        if self.period is None:
            return

        # Define limits first
        self.define_cutvariation_limits()


        # list containing the ML cuts for each pT bin
        # [[cut1, cut2, cut3...], [cut1, cut2, cut3, ...], ...for all trials

        n_trials = 2 * self.p_ncutvar + 1
        ml_wps = [[] for _ in n_trials]

        print("Using run selection for mass histo for period", self.p_period)
        for ipt in range(self.p_nptfinbins):
            bin_id = self.bin_matching[ipt]

            stepsmin = (self.lpt_probcutfin[bin_id] - self.min_cv_cut[ipt]) / self.p_ncutvar
            stepsmax = (self.max_cv_cut[ipt] - self.lpt_probcutfin[bin_id]) / self.p_ncutvar
            
            for icv in range(self.p_ncutvar):
                lower_cut = Systematics.min_cv_cut[ipt] + icv * stepsmin
                upper_cut = self.lpt_probcutfin[bin_id] + (icv + 1) * stepsmax

                selml_cv_low = "y_test_prob%s>%s" % (self.p_modelname, lower_cut)
                selml_cv_up = "y_test_prob%s>%s" % (self.p_modelname, upper_cut)
                ml_wps[icv].append(selml_cv_low)
                ml_wps[self.p_ncutvar + icv + 1].append(selml_cv_up)

            ml_wps[self.p_ncutvar].append(self.lpt_probcutfin[bin_id])

        self.processers_syst = []

        for i, sel_pt in enumerate(ml_wps):
            # modify output path and create it
            new_output_mc = join(self.nominal_prcesser_mc.d_results, self.syst_out_dir, f"trial_{i}")
            new_output_data = join(self.nominal_prcesser_data.d_results, self.syst_out_dir, f"trial_{i}")
            if not exists(new_output_mc):
                makedirs(new_output_mc)
            if not exists(new_output_data):
                makedirs(new_output_data)
            # Make a new processer and append
            syst_processer_mc = self.processer_nominal_mc.__class__(self.nominal_processer_mc.case,
                                                                    self.nominal_processer_mc.datap,
                                                                    self.nominal_processer_mc.run_param,
                                                                    self.nominal_processer_mc.mcordata,
                                                                    self.nominal_processer_mc.p_maxfiles,
                                                                    self.nominal_processer_mc.d_root,
                                                                    self.nominal_processer_mc.d_pkl,
                                                                    self.nominal_processer_mc.d_pklsk,
                                                                    self.nominal_processer_mc.d_pkl_ml,
                                                                    self.nominal_processer_mc.p_period,
                                                                    self.nominal_processer_mc.i_period,
                                                                    self.nominal_processer_mc.p_chunksizeunp,
                                                                    self.nominal_processer_mc.p_chunksizesk,
                                                                    self.nominal_processer_mc.p_maxprocess,
                                                                    self.nominal_processer_mc.p_frac_merge,
                                                                    self.nominal_processer_mc.p_rd_merge,
                                                                    self.nominal_processer_mc.d_pkl_dec,
                                                                    self.nominal_processer_mc.d_pkl_decmerged,
                                                                    new_output_mc,
                                                                    self.nominal_processer_mc.typean,
                                                                    self.nominal_processer_mc.runlisttrigger,
                                                                    self.nominal_processer_mc.d_mcreweights)


            # Adjust ml cuts
            syst_processer_mc.lpt_probcutfin = sel_pt
            syst_processer_mc.l_selml = ["y_test_prob%s>%s" % (syst_processer_mc.p_modelname, syst_processer_mc.lpt_probcutfin[ipt]) \
                                           for ipt in range(syst_processer_mc.p_nptbins)]
            self.processers_mc_syst.append(syst_processer_mc)

            syst_processer_data = self.processer_nominal_data.__class__(self.nominal_processer_mc.case,
                                                                        self.nominal_processer_mc.datap,
                                                                        self.nominal_processer_mc.run_param,
                                                                        self.nominal_processer_mc.mcordata,
                                                                        self.nominal_processer_mc.p_maxfiles,
                                                                        self.nominal_processer_mc.d_root,
                                                                        self.nominal_processer_mc.d_pkl,
                                                                        self.nominal_processer_mc.d_pklsk,
                                                                        self.nominal_processer_mc.d_pkl_ml,
                                                                        self.nominal_processer_mc.p_period,
                                                                        self.nominal_processer_mc.i_period,
                                                                        self.nominal_processer_mc.p_chunksizeunp,
                                                                        self.nominal_processer_mc.p_chunksizesk,
                                                                        self.nominal_processer_mc.p_maxprocess,
                                                                        self.nominal_processer_mc.p_frac_merge,
                                                                        self.nominal_processer_mc.p_rd_merge,
                                                                        self.nominal_processer_mc.d_pkl_dec,
                                                                        self.nominal_processer_mc.d_pkl_decmerged,
                                                                        new_output_data,
                                                                        self.nominal_processer_mc.typean,
                                                                        self.nominal_processer_mc.runlisttrigger,
                                                                        self.nominal_processer_mc.d_mcreweights)


            # Adjust ml cuts
            syst_processer_data.lpt_probcutfin = sel_pt
            syst_processer_data.l_selml = ["y_test_prob%s>%s" % (syst_processer_data.p_modelname, syst_processer_data.lpt_probcutfin[ipt]) \
                                           for ipt in range(syst_processer_data.p_nptbins)]
            self.processers_data_syst.append(syst_processer_data)
            syst_processer_data.process_histomass()

        self.done_mass = True


    def ml_cutvar_eff(self):
        """
        Cut Variation: Create ROOT file with efficiencies
        Histogram for each variation, for each 2nd binning bin

        Similar as process_efficiency_single(self, index) in processor.py
        """

        # Do this only for a single period
        if self.period is None:
            return

        # Define limits first
        self.define_cutvariation_limits()

        for p in self.processers_mc_syst:
            p.process_efficiency()

        self.done_eff = True


    # pylint: disable=import-outside-toplevel
    def ml_cutvar_fit(self):
        """
        Cut Variation: Fit invariant mass histograms with AliHFInvMassFitter
        If requested, sigma+mean can be fixed to central fit

        Similar as fitter(self) in analyzer.py
        """

        if not self.nominal_analyzer:
            return

        # Define limits first
        if self.period is not None and not self.done_mass:
            self.logger.fatal("Cannot fit since mass histograms have not been produced yet.")

        # Put together all analyzers per trial
        self.analyzers_syst = []
        for i, sel_pt in enumerate(ml_wps):
            # modify output path and create it
            new_output_mc = join(self.nominal_analyzer.d_resultsallpmc, self.syst_out_dir, f"trial_{i}")
            new_output_data = join(self.nominal_analyzer.d_resultsallpdata, self.syst_out_dir, f"trial_{i}")

            datap == deepcopy(self.nominal_processer_data.datap)
            datap["mlapplication"]["probcutoptimal"] = sel_pt

            # update i_period result path
            typean = self.nominal_prcesser_data.typean
            i_period = self.nominal_processer_data.i_period
            if i_period is not None:
                datap["analysis"][typean]["mc"]["results"][i_period] = new_output_mc
                datap["analysis"][typean]["data"]["results"][i_period] = new_output_data
            else:
                datap["analysis"][typean]["mc"]["resultsallp"] = new_output_mc
                datap["analysis"][typean]["data"]["resultsallp"] = new_output_data


            # construct analyzer
            analyzer_new = self.nominal_analyzer.__class__(datap, self.nominal_analyzer.case,
                                                           typean, i_period)

            self.analyzers_syst.append(analyzer_new)
            analyzer_new.fit()
            analyzer.efficiency()
            analyzer.makenormyields()
            analyzer.plotternormyields()

        self.done_fit = True



    def ml_cutvar_cross(self):

        # Make normalized yields first
        self.cutvariation_makenormyields()
        # Then plot for each of these
        for name in ["histoSigmaCorr", "hDirectEffpt", "hFeedDownEffpt", "hRECpt"]:
            self.ml_cutvar_makeplots(name)


    def load_central_meansigma(self, imult):
        """
        Cut Variation: Get parameters (mean and sigma) from central fit
        """
        func_filename_std = make_file_path(self.d_results, self.yields_filename_std, \
                                           "root", None, [self.case, self.typean])

        massfile_std = TFile.Open(func_filename_std, "READ")
        means_histo = massfile_std.Get("hmeanss%d" % (imult))
        sigmas_histo = massfile_std.Get("hsigmas%d" % (imult))

        mean_for_data = []
        sigma_for_data = []
        for ipt in range(self.p_nptfinbins):

            mean_for_data.append(means_histo.GetBinContent(ipt + 1))
            sigma_for_data.append(sigmas_histo.GetBinContent(ipt + 1))

        massfile_std.Close()
        return mean_for_data, sigma_for_data


    def mcptshape_get_generated(self):
        """
        MC pT-shape: Get generated pT spectra from MC to define weights
        """
        fileout_name = make_file_path(self.d_results, self.ptspectra_filename, \
                                      "root", None, [self.typean, self.case])
        myfile = TFile(fileout_name, "RECREATE")

        for ibin2 in range(len(self.lvar2_binmin)):
            stringbin2 = "_%s_%.2f_%.2f" % (self.v_var2_binning_gen, \
                                        self.lvar2_binmin[ibin2], \
                                        self.lvar2_binmax[ibin2])

            h_gen_pr = TH1F("h_gen_pr" + stringbin2, "Prompt Generated in acceptance |y|<0.5", \
                            400, 0, 40)
            h_gen_fd = TH1F("h_gen_fd" + stringbin2, "FD Generated in acceptance |y|<0.5", \
                            400, 0, 40)

            for ipt in range(self.p_nptfinbins):
                bin_id = self.bin_matching[ipt]

                df_mc_gen = pickle.load(openfile(self.lpt_gendecmerged[bin_id], "rb"))
                df_mc_gen = df_mc_gen.query("abs(y_cand) < 0.5")

                df_mc_gen = seldf_singlevar(df_mc_gen, self.v_var_binning, \
                                     self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt])
                df_mc_gen = seldf_singlevar(df_mc_gen, self.v_var2_binning_gen, \
                                            self.lvar2_binmin[ibin2], self.lvar2_binmax[ibin2])

                df_gen_sel_pr = df_mc_gen[df_mc_gen.ismcprompt == 1]
                df_gen_sel_fd = df_mc_gen[df_mc_gen.ismcfd == 1]

                fill_hist(h_gen_pr, df_gen_sel_pr.pt_cand)
                fill_hist(h_gen_fd, df_gen_sel_fd.pt_cand)
            myfile.cd()
            h_gen_pr.Write()
            h_gen_fd.Write()
        myfile.Close()


    def mcptshape_build_efficiencies(self):
        """
        MC pT-shape: Create ROOT file with unweighted and weighted efficiencies
        Histogram for (un)weighted, for each 2nd binning bin

        Similar as process_efficiency_single(self, index) in processor.py
        """
        myfile = TFile.Open(self.n_fileeff_ptshape, "recreate")

        print("Using run selection for eff histo for period", self.p_period)
        for ibin2 in range(len(self.lvar2_binmin)):
            stringbin2 = "_%s_%.2f_%.2f" % (self.v_var2_binning_gen, \
                                        self.lvar2_binmin[ibin2], \
                                        self.lvar2_binmax[ibin2])

            n_bins = len(self.lpt_finbinmin)
            analysis_bin_lims_temp = self.lpt_finbinmin.copy()
            analysis_bin_lims_temp.append(self.lpt_finbinmax[n_bins-1])
            analysis_bin_lims = array('f', analysis_bin_lims_temp)

            h_gen_pr = TH1F("h_gen_pr" + stringbin2, \
                            "Prompt Generated in acceptance |y|<0.5", \
                            n_bins, analysis_bin_lims)
            h_presel_pr = TH1F("h_presel_pr" + stringbin2, \
                               "Prompt Reco in acc |#eta|<0.8 and sel", \
                               n_bins, analysis_bin_lims)
            h_sel_pr = TH1F("h_sel_pr" + stringbin2, \
                            "Prompt Reco and sel in acc |#eta|<0.8 and sel", \
                            n_bins, analysis_bin_lims)
            h_gen_fd = TH1F("h_gen_fd" + stringbin2, \
                            "FD Generated in acceptance |y|<0.5", \
                            n_bins, analysis_bin_lims)
            h_presel_fd = TH1F("h_presel_fd" + stringbin2, \
                               "FD Reco in acc |#eta|<0.8 and sel", \
                               n_bins, analysis_bin_lims)
            h_sel_fd = TH1F("h_sel_fd" + stringbin2, \
                            "FD Reco and sel in acc |#eta|<0.8 and sel", \
                            n_bins, analysis_bin_lims)

            bincounter = 0
            for ipt in range(self.p_nptfinbins):
                bin_id = self.bin_matching[ipt]
                selml = "y_test_prob%s>%s" % (self.p_modelname, self.lpt_probcutfin[bin_id])

                df_mc_reco = pickle.load(openfile(self.lpt_recodecmerged_mc[bin_id], "rb"))
                if self.s_evtsel is not None:
                    df_mc_reco = df_mc_reco.query(self.s_evtsel)
                if self.s_trigger_mc is not None:
                    df_mc_reco = df_mc_reco.query(self.s_trigger_mc)

                df_mc_gen = pickle.load(openfile(self.lpt_gendecmerged[bin_id], "rb"))
                df_mc_gen = df_mc_gen.query(self.s_presel_gen_eff)

                df_mc_reco = seldf_singlevar(df_mc_reco, self.v_var_binning, \
                                     self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt])
                df_mc_gen = seldf_singlevar(df_mc_gen, self.v_var_binning, \
                                     self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt])

                df_mc_reco = seldf_singlevar(df_mc_reco, self.v_var2_binning_gen, \
                                             self.lvar2_binmin[ibin2], self.lvar2_binmax[ibin2])
                df_mc_gen = seldf_singlevar(df_mc_gen, self.v_var2_binning_gen, \
                                            self.lvar2_binmin[ibin2], self.lvar2_binmax[ibin2])

                df_gen_sel_pr = df_mc_gen[df_mc_gen.ismcprompt == 1]
                df_reco_presel_pr = df_mc_reco[df_mc_reco.ismcprompt == 1]
                df_reco_sel_pr = None
                df_reco_sel_pr = df_reco_presel_pr.query(selml)

                df_gen_sel_fd = df_mc_gen[df_mc_gen.ismcfd == 1]
                df_reco_presel_fd = df_mc_reco[df_mc_reco.ismcfd == 1]
                df_reco_sel_fd = None
                df_reco_sel_fd = df_reco_presel_fd.query(selml)

                val = len(df_gen_sel_pr)
                err = math.sqrt(val)
                h_gen_pr.SetBinContent(bincounter + 1, val)
                h_gen_pr.SetBinError(bincounter + 1, err)
                val = len(df_reco_presel_pr)
                err = math.sqrt(val)
                h_presel_pr.SetBinContent(bincounter + 1, val)
                h_presel_pr.SetBinError(bincounter + 1, err)
                val = len(df_reco_sel_pr)
                err = math.sqrt(val)
                h_sel_pr.SetBinContent(bincounter + 1, val)
                h_sel_pr.SetBinError(bincounter + 1, err)

                val = len(df_gen_sel_fd)
                err = math.sqrt(val)
                h_gen_fd.SetBinContent(bincounter + 1, val)
                h_gen_fd.SetBinError(bincounter + 1, err)
                val = len(df_reco_presel_fd)
                err = math.sqrt(val)
                h_presel_fd.SetBinContent(bincounter + 1, val)
                h_presel_fd.SetBinError(bincounter + 1, err)
                val = len(df_reco_sel_fd)
                err = math.sqrt(val)
                h_sel_fd.SetBinContent(bincounter + 1, val)
                h_sel_fd.SetBinError(bincounter + 1, err)

                bincounter = bincounter + 1

            hw_gen_pr = TH1F("h_gen_pr_weight" + stringbin2, "Prompt Generated in acc |y|<0.5", \
                             n_bins, analysis_bin_lims)
            hw_presel_pr = TH1F("h_presel_pr_weight" + stringbin2, \
                                "Prompt Reco in acc |#eta|<0.8 and sel", \
                                 n_bins, analysis_bin_lims)
            hw_sel_pr = TH1F("h_sel_pr_weight" + stringbin2, \
                             "Prompt Reco and sel in acc |#eta|<0.8 and sel", \
                             n_bins, analysis_bin_lims)
            hw_gen_fd = TH1F("h_gen_fd_weight" + stringbin2, "FD Generated in acc |y|<0.5", \
                             n_bins, analysis_bin_lims)
            hw_presel_fd = TH1F("h_presel_fd_weight" + stringbin2, \
                                "FD Reco in acc |#eta|<0.8 and sel", \
                                n_bins, analysis_bin_lims)
            hw_sel_fd = TH1F("h_sel_fd_weight" + stringbin2, \
                             "FD Reco and sel in acc |#eta|<0.8 and sel", \
                             n_bins, analysis_bin_lims)

            bincounter = 0
            for ipt in range(self.p_nptfinbins):
                bin_id = self.bin_matching[ipt]
                selml = "y_test_prob%s>%s" % (self.p_modelname, self.lpt_probcutfin[bin_id])

                df_mc_reco = pickle.load(openfile(self.lpt_recodecmerged_mc[bin_id], "rb"))
                if self.s_evtsel is not None:
                    df_mc_reco = df_mc_reco.query(self.s_evtsel)
                if self.s_trigger_mc is not None:
                    df_mc_reco = df_mc_reco.query(self.s_trigger_mc)

                df_mc_gen = pickle.load(openfile(self.lpt_gendecmerged[bin_id], "rb"))
                df_mc_gen = df_mc_gen.query(self.s_presel_gen_eff)

                df_mc_reco = seldf_singlevar(df_mc_reco, self.v_var_binning, \
                                     self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt])
                df_mc_gen = seldf_singlevar(df_mc_gen, self.v_var_binning, \
                                     self.lpt_finbinmin[ipt], self.lpt_finbinmax[ipt])

                df_mc_reco = seldf_singlevar(df_mc_reco, self.v_var2_binning_gen, \
                                             self.lvar2_binmin[ibin2], self.lvar2_binmax[ibin2])
                df_mc_gen = seldf_singlevar(df_mc_gen, self.v_var2_binning_gen, \
                                            self.lvar2_binmin[ibin2], self.lvar2_binmax[ibin2])

                df_gen_sel_pr = df_mc_gen[df_mc_gen.ismcprompt == 1]
                df_reco_presel_pr = df_mc_reco[df_mc_reco.ismcprompt == 1]
                df_reco_sel_pr = None
                df_reco_sel_pr = df_reco_presel_pr.query(selml)

                df_gen_sel_fd = df_mc_gen[df_mc_gen.ismcfd == 1]
                df_reco_presel_fd = df_mc_reco[df_mc_reco.ismcfd == 1]
                df_reco_sel_fd = None
                df_reco_sel_fd = df_reco_presel_fd.query(selml)

                array_pt_gencand_gen = df_gen_sel_pr.pt_cand.values
                array_pt_recocand_reco_presel = df_reco_presel_pr.pt_cand.values
                array_pt_recocand_reco_sel = df_reco_sel_pr.pt_cand.values

                val, err = self.get_reweighted_count(array_pt_gencand_gen)
                hw_gen_pr.SetBinContent(bincounter + 1, val)
                hw_gen_pr.SetBinError(bincounter + 1, err)
                val, err = self.get_reweighted_count(array_pt_recocand_reco_presel)
                hw_presel_pr.SetBinContent(bincounter + 1, val)
                hw_presel_pr.SetBinError(bincounter + 1, err)
                val, err = self.get_reweighted_count(array_pt_recocand_reco_sel)
                hw_sel_pr.SetBinContent(bincounter + 1, val)
                hw_sel_pr.SetBinError(bincounter + 1, err)

                array_pt_gencand_genfd = df_gen_sel_fd.pt_cand.values
                array_pt_recocand_reco_preselfd = df_reco_presel_fd.pt_cand.values
                array_pt_recocand_reco_selfd = df_reco_sel_fd.pt_cand.values

                val, err = self.get_reweighted_count(array_pt_gencand_genfd)
                hw_gen_fd.SetBinContent(bincounter + 1, val)
                hw_gen_fd.SetBinError(bincounter + 1, err)
                val, err = self.get_reweighted_count(array_pt_recocand_reco_preselfd)
                hw_presel_fd.SetBinContent(bincounter + 1, val)
                hw_presel_fd.SetBinError(bincounter + 1, err)
                val, err = self.get_reweighted_count(array_pt_recocand_reco_selfd)
                hw_sel_fd.SetBinContent(bincounter + 1, val)
                hw_sel_fd.SetBinError(bincounter + 1, err)

                bincounter = bincounter + 1

            myfile.cd()
            h_gen_pr.Write()
            h_presel_pr.Write()
            h_sel_pr.Write()
            h_gen_fd.Write()
            h_presel_fd.Write()
            h_sel_fd.Write()
            hw_gen_pr.Write()
            hw_presel_pr.Write()
            hw_sel_pr.Write()
            hw_gen_fd.Write()
            hw_presel_fd.Write()
            hw_sel_fd.Write()
        myfile.Close()

    def mcptshape_efficiency(self):
        """
        MC pT-shape: Extract prompt and feeddown efficiencies
        Systematic = difference wrt 1 for ratio unweighted / weighted
        """
        load_root_style_simple()

        lfileeff = TFile.Open(self.n_fileeff_ptshape, "READ")

        fileout_name = make_file_path(self.d_results, self.efficiency_filename_pt, \
                                      "root", None, [self.typean, self.case])
        fileout = TFile(fileout_name, "RECREATE")

        for imult in range(len(self.lvar2_binmin)):

            stringbin2 = "_%s_%.2f_%.2f" % (self.v_var2_binning_gen, \
                                           self.lvar2_binmin[imult], \
                                           self.lvar2_binmax[imult])

            h_gen_pr = lfileeff.Get("h_gen_pr" + stringbin2)
            h_sel_pr = lfileeff.Get("h_sel_pr" + stringbin2)
            h_sel_pr.Divide(h_sel_pr, h_gen_pr, 1.0, 1.0, "B")

            h_gen_fd = lfileeff.Get("h_gen_fd" + stringbin2)
            h_sel_fd = lfileeff.Get("h_sel_fd" + stringbin2)
            h_sel_fd.Divide(h_sel_fd, h_gen_fd, 1.0, 1.0, "B")

            hw_gen_pr = lfileeff.Get("h_gen_pr_weight" + stringbin2)
            hw_sel_pr = lfileeff.Get("h_sel_pr_weight" + stringbin2)
            hw_sel_pr.Divide(hw_sel_pr, hw_gen_pr, 1.0, 1.0, "B")

            hw_gen_fd = lfileeff.Get("h_gen_fd_weight" + stringbin2)
            hw_sel_fd = lfileeff.Get("h_sel_fd_weight" + stringbin2)
            hw_sel_fd.Divide(hw_sel_fd, hw_gen_fd, 1.0, 1.0, "B")

            fileout.cd()
            h_sel_pr.SetName("eff_mult%d" % imult)
            h_sel_fd.SetName("eff_fd_mult%d" % imult)
            hw_sel_pr.SetName("eff_weight_mult%d" % imult)
            hw_sel_fd.SetName("eff_weight_fd_mult%d" % imult)
            h_sel_pr.Write()
            h_sel_fd.Write()
            hw_sel_pr.Write()
            hw_sel_fd.Write()
        fileout.Close()

    def mcptshape_makeplots(self):
        """
        MC pT shape: Make final plots.
        For the moment, value should be assigned by analyser
        """
        load_root_style()

        leg = TLegend(.15, .65, .85, .85)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.024)
        colours = [kBlack, kRed]
        markers = [20, 21]

        fileout_name = make_file_path(self.d_results, self.efficiency_filename_pt, \
                                      "root", None, [self.typean, self.case])

        f_fileout = TFile.Open(fileout_name)

        hweights = []
        hnoweights = []
        hfdweights = []
        hfdnoweights = []
        for imult in range(len(self.lvar2_binmin)):

            canv = TCanvas("systmcptshape_%d" % imult, '', 400, 400)
            plotname = "No weights / Weights"
            ptmax = self.lpt_finbinmax[-1] + 1
            canv.cd(1).DrawFrame(0, 0.85, ptmax, 1.15, \
                                 "%s %.2f < %s < %.2f;#it{p}_{T} (GeV/#it{c});Ratio %s" % \
                                 (self.typean, self.lvar2_binmin[imult], self.v_var2_binning, \
                                  self.lvar2_binmax[imult], plotname))

            hweights.append(f_fileout.Get("eff_weight_mult%d" % imult))
            hweights[imult].SetDirectory(0)
            hnoweights.append(f_fileout.Get("eff_mult%d" % imult))
            hnoweights[imult].SetDirectory(0)
            hfdweights.append(f_fileout.Get("eff_weight_fd_mult%d" % imult))
            hfdweights[imult].SetDirectory(0)
            hfdnoweights.append(f_fileout.Get("eff_fd_mult%d" % imult))
            hfdnoweights[imult].SetDirectory(0)

            hnoweights[imult].Divide(hnoweights[imult], hweights[imult], 1.0, 1.0, "B")
            hfdnoweights[imult].Divide(hfdnoweights[imult], hfdweights[imult], 1.0, 1.0, "B")

            hnoweights[imult].SetLineColor(colours[0])
            hnoweights[imult].SetMarkerColor(colours[0])
            hnoweights[imult].SetMarkerStyle(markers[0])
            hnoweights[imult].SetMarkerSize(0.8)
            hnoweights[imult].Draw("same")

            hfdnoweights[imult].SetLineColor(colours[1])
            hfdnoweights[imult].SetMarkerColor(colours[1])
            hfdnoweights[imult].SetMarkerStyle(markers[1])
            hfdnoweights[imult].SetMarkerSize(0.8)
            hfdnoweights[imult].Draw("same")

            if imult == 0:
                leg.AddEntry(hnoweights[imult], "Prompt", "LEP")
                leg.AddEntry(hfdnoweights[imult], "Feed-down", "LEP")

            leg.Draw()
            canv.SaveAs("%s/MCpTshape_Syst_mult%d.eps" % (self.d_results, imult))
        f_fileout.Close()


    def mcptshape(self):

        # Do only per period
        if self.period is not None:
            self.mcptshape_get_generated()
            self.mcptshape_build_efficiencies()

        # Do for all
        self.mcptshape_efficiency()
        self.mcptshape_makeplots()


    def get_reweighted_count(self, arraypt):
        """
        MC pT-shape: Reweight array of pTs from dataframe based on pT weights
        """
        weights = arraypt.copy()
        binwidth = (self.p_weights_max_pt - self.p_weights_min_pt)/self.p_weights_bins
        for j in range(weights.shape[0]):
            pt = arraypt[j]
            if pt - self.p_weights_min_pt < 0:
                self.logger.warning("pT_gen < minimum pT of weights!")
            ptbin_weights = int((pt - self.p_weights_min_pt)/binwidth)
            #improvement: make linear extrapolation with bins next to it
            weights[j] = self.p_weights[ptbin_weights]
        val = sum(weights)
        err = math.sqrt(val)
        return val, err

