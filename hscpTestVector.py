#!/usr/bin/env python3
import ROOT
from ROOT import TCanvas,TMultiGraph, TFile, TProfile, TNtuple, TH1F, TH2F,TH2D, TGraph, TGraphErrors, TEfficiency
from array import array
import json
import math
import re
import ctypes
import struct
import argparse
import sys
import ctypes


from DataFormats.FWLite import Events, Handle
from Analysis.HLTAnalyserPy.EvtData import EvtData, EvtHandles, add_product,get_objs

from Analysis.HLTAnalyserPy.CaloVectors import SetBin,FillHist,ActivityTab,SetActivityAround,PassThreshold,CreateCaloMap,Get12NeighboursLoop,CheckHighEMatrix,deltaR2,deltaR,testDeltaR2

from Analysis.HLTAnalyserPy.CaloVectors import Create3by3MatrixLoop,CreateCaloTabLoop,FindTrueSeed,FindAdjacentPair,IsInMatrix,FurthestInRatio,FillStepByStep,FillPdgIds,FindAllHSCP,FindChargedHSCP,CountIsoNum

from Analysis.HLTAnalyserPy.ArraysHSCP import ListCuts,idx_pdg_ch,idx_pdg_all,ListNames,hlt_trig_names

import Analysis.HLTAnalyserPy.CoreTools as CoreTools
import Analysis.HLTAnalyserPy.GenTools as GenTools
import Analysis.HLTAnalyserPy.HistTools as HistTools
import Analysis.HLTAnalyserPy.TrigTools as TrigTools



if __name__ == "__main__":

    CoreTools.load_fwlitelibs()

    parser = argparse.ArgumentParser(description='example e/gamma and/or muon HLT analyser')
    parser.add_argument('in_filenames',nargs="+",help='input filename')
    parser.add_argument('--prefix','-p',default='file:',help='file prefix')
    parser.add_argument('--l1menufile',required=True,help='l1menu file')
    parser.add_argument('--out','-o',default="output.root",help='output filename')
    parser.add_argument('eventmap',default=1,help='in which event you seek the map')
    parser.add_argument('minmipseed',default=1,help='how many MIP will the candidate leave in both CAL')
    parser.add_argument('nbmaps',default=20,help='how many eta-phi maps you want')
    parser.add_argument('minmipngh',default=1.5,help='minimal energy deposit inside neighbouring calo towers (in MIPs)')
    parser.add_argument('maxsumcut',default=10,help='maximal ecal+hcal deposit inside calo towers (in MIPs)')
    parser.add_argument('etacut',default=1,help='if we ask for a cut on calo towers eta')
    parser.add_argument('reversecheck',default=0,help='if we want a reverse check between all towers-hscp dr < 0.1')
    args = parser.parse_args()



   

    hfile = TFile( 'CaloStudiesHSCP_iso_vector_1_halfMIP.root', 'RECREATE', 'Output of hscp script to save HLT quantities' )

    SizeArrays = len(ListNames)

    EventToStudy = sys.argv[4]
    min_mip = float(sys.argv[5])
    min_mip_ngh = float(sys.argv[7])
    max_sum_mip = float(sys.argv[8])

    ask_cut = bool(int(sys.argv[9]))
    reverse_bool = bool(int(sys.argv[10]))

    print("We will study CaloTowers above ",min_mip, "MIP , and their neighbours above ", min_mip_ngh, " MIP , and sum (ECAL + HCAL) < ",max_sum_mip, " MIP")
    if ask_cut:
        print("We also ask all CaloTowers to have |eta| < 2.4")
    else:
        print("No requirements on Calo Towers")
    #============ CALO TOWER HISTOS ===========

    radius_idx_map=5


    nb_tower_per_event = TH1F('nb_tower_per_event', 'nb calo tower ',1000,0,1000)
    nb_tower_non_valid_hist = TH1F('nb calo tower non valid (for whole event)', 'nb calo towers missing',100000,0,100000)
    #ratio_mu_ecalhcal = TH1F ( 'ratio_ecal_hcal_MU', 'E_ecal / E_hcal ',100,0,1)
    recoCaloMet = TEfficiency("eff","reco::caloMET efficiency;reco::caloMET;#epsilon",100,0,1000)
    step_by_step_tower = TH1F('info nb tower step by step','tower number after each cut',6,0,6)


    nb_towers_8_above_trh = TH1F('number_of_neighbours_passing_cuts','tower nb passing cuts',9,0,9)
    nb_towers_8_below_trh = TH1F('number of neighbours (8 max) NOT passing cuts (same as above)','tower nb not passing cuts',9,0,9)
    nb_towers_8_all_trh = TH1F('number of neighbours per cell','nb of neibhours per cell',9,0,9)

    Sum_Energy_8_neighbours_threshold = TH1F('Sum_8_neighbors_tower','(Sum of X cells (below threshold) around cell of interest / (8- # cells)',400,0,20)

    Sum_Energy_GEN_HSCP_8_neighbours_threshold = TH1F('Sum_GEN_HSCP_8_neighbors_tower','(Sum of X cells (below threshold) around matched cell (HSCP GEN) / (8- # cells)',100,0,5)


    nb_towers_8_above_trh_HSCP_RECO_ISO = TH1F('HSCP_RECO_ISO_number_of_neighbours_passing_cuts','tower nb passing cuts',9,0,9)
    nb_towers_8_below_trh_HSCP_RECO_ISO = TH1F('HSCP_RECO_ISO_number_of_neighbours_NOT_passing_cuts','tower nb not passing cuts',9,0,9)
    nb_towers_8_all_trh_HSCP_RECO_ISO = TH1F('HSCP_RECO_ISO_number_of_neighbours_per_cell','# of neibhours per cell',9,0,9)

    Sum_Energy_8_neighbours_threshold_HSCP_RECO_ISO = TH1F('HSCP_RECO_ISO_Sum_8_neighbors_tower','(Sum of X cells (below threshold) around cell of interest iso(dr hscp calo tower min) / (8- # cells)',200,0,10)

    Sum_Energy_8_neighbours_threshold_HSCP_RECO_ISO_matched = TH1F('HSCP_RECO_ISO_Sum_8_neighbors_tower_matched','(Sum of X cells (below threshold) around cell of interest iso(dr hscp calo tower min) / (8- # cells)',200,0,10)

    Number_shift_seedOI = TH1F("Number_of_seed_shifts","# of shifts from initial seed to next one",10,0,10)

    ISO_HSCP_Number_shift_seedOI = TH1F("ISO_HSCP_Number_of_seed_shifts","# of shifts from initial seed to next one",10,0,10)

    HSCP_ISO_Number_non_physical_event = TH1F("Number_tower_more_3_ngh_hscp_iso","# of non-physical event",2,0,2)
    Number_non_physical_event = TH1F("Number_tower_more_3_ngh","# of non-physical event",2,0,2)
   
 

    if reverse_bool:
        E_Ecal_dr01 = TH1F("E_ECAL_dr01_genhscp_tower","Energy in Ecal for all calo towers - hscp dr < 0.1",80,0,4) 
        E_Hcal_dr01 = TH1F("E_HCAL_dr01_genhscp_tower","Energy in Hcal for all calo towers - hscp dr < 0.1",200,0,10) 
        E_Sum_dr01 = TH1F("E_SUM_dr01_genhscp_tower","Energy in Hcal+Ecal for all calo towers - hscp dr < 0.1",240,0,12) 
        E_Ratio_dr01 = TH1F("E_RATIO_dr01_genhscp_tower","Energy ratio ecal/hcal for all (calo towers - hscp) with dr < 0.1",100,0,1) 
        Nb_ngh_dr01 = TH1F("Nb_ngh_dr01_genhscp_tower","Number of neighbouring seeds (dr < 0.1) per event",10,0,10) 


    Number_neighbours_per_seed = TH1F('nb_neighbours_per_seed','# neighbours per seed',9,0,9)
    Nb_neighbour_within_2by2 = TH1F('nb_neighbours_within_2by2_all','# neighbours within 2x2 matrix',4,0,4)
    Number_HSCP_2by2_max_event = TH1F('nb_neighbours_within_2by2_max_event_matched_HSCP','nb HSCP',10,0,10)


    nb_matched_seed_2_charged = TH1F("nb_matched_Seed_when_2_charged","nb of matched seeds when having 2 charged HSCP", 4,0,4)
    nb_matched_seed_1_charged = TH1F("nb_matched_Seed_when_1_charged","nb of matched seeds when having 1 charged HSCP", 4,0,4)


    nb_matched_seed_2_charged_reco = TH1F("nb_matched_Seed_when_2_reco_charged","nb of matched seeds when having 2 reco charged HSCP", 4,0,4)
    nb_matched_seed_1_charged_reco = TH1F("nb_matched_Seed_when_1_reco_charged","nb of matched seeds when having 1 reco charged HSCP", 4,0,4)
    #----------------- RECO HSCP HISTOGRAMS ----------------------

    NB_RECO_HSCP_vs_matched_seed = TH2D("reco_hscp_vs_matched_seed","nb reco hscp vs nb matched seed per event", 3,0,3,5,0,5)
    NB_RECO_HSCP_vs_matched_seed.GetXaxis().SetTitle("# reco hscp")
    NB_RECO_HSCP_vs_matched_seed.GetYaxis().SetTitle("# matched seed")




    #----------------- GEN HSCP HISTOGRAMS ----------------------

    NB_GEN_HSCP_vs_matched_seed = TH2D("gen_hscp_vs_matched_seed","nb gen hscp vs nb matched seed per event", 3,0,3,5,0,5)
    NB_GEN_HSCP_vs_matched_seed.GetXaxis().SetTitle("# gen hscp")
    NB_GEN_HSCP_vs_matched_seed.GetYaxis().SetTitle("# matched seed")

    NB_GEN_HSCP_vs_NBTOWER = TH2D("gen_hscp_vs_tower","nb gen hscp vs nb tower per event", 3,0,3,15,0,15)
    NB_GEN_HSCP_vs_NBTOWER.GetXaxis().SetTitle("# gen hscp")
    NB_GEN_HSCP_vs_NBTOWER.GetYaxis().SetTitle("# tower per event")


    NB_GEN_HSCP_vs_NBTOWER_isovar03 = TH2D("gen_hscp_vs_tower_iso03","nb gen hscp vs nb tower iso v < 0.3", 3,0,3,15,0,15)
    NB_GEN_HSCP_vs_NBTOWER_isovar03.GetXaxis().SetTitle("# gen hscp")
    NB_GEN_HSCP_vs_NBTOWER_isovar03.GetYaxis().SetTitle("# tower with Iso v < 0.3")


    NB_GEN_HSCP_vs_NBTOWER_isovar05 = TH2D("gen_hscp_vs_tower_iso05","nb gen hscp vs nb tower iso v < 0.5", 3,0,3,15,0,15)
    NB_GEN_HSCP_vs_NBTOWER_isovar05.GetXaxis().SetTitle("# gen hscp")
    NB_GEN_HSCP_vs_NBTOWER_isovar05.GetYaxis().SetTitle("# tower with Iso v < 0.5")
 
 
    NB_GEN_HSCP_vs_NBTOWER_isovar07 = TH2D("gen_hscp_vs_tower_iso07","nb gen hscp vs nb tower iso v < 0.7", 3,0,3,15,0,15)
    NB_GEN_HSCP_vs_NBTOWER_isovar07.GetXaxis().SetTitle("# gen hscp")
    NB_GEN_HSCP_vs_NBTOWER_isovar07.GetYaxis().SetTitle("# tower with Iso v < 0.7")

    NB_GEN_HSCP_vs_NBTOWER_isovar09 = TH2D("gen_hscp_vs_tower_iso09","nb gen hscp vs nb tower iso v < 0.9", 3,0,3,15,0,15)
    NB_GEN_HSCP_vs_NBTOWER_isovar09.GetXaxis().SetTitle("# gen hscp")
    NB_GEN_HSCP_vs_NBTOWER_isovar09.GetYaxis().SetTitle("# tower with Iso v < 0.9")

    NB_RECO_HSCP_vs_NBTOWER = TH2D("reco_hscp_vs_tower","nb hscp vs nb tower per event", 15,0,15,25,0,25)
    NB_RECO_HSCP_vs_NBTOWER.GetXaxis().SetTitle("# hscp")
    NB_RECO_HSCP_vs_NBTOWER.GetYaxis().SetTitle("# tower per event")

    NB_CLEAN_RECO_HSCP_vs_NBTOWER = TH2D("clean_reco_hscp_vs_tower","nb clean reco hscp vs nb tower per event", 3,0,3,15,0,15)
    NB_CLEAN_RECO_HSCP_vs_NBTOWER.GetXaxis().SetTitle("#clean reco hscp")
    NB_CLEAN_RECO_HSCP_vs_NBTOWER.GetYaxis().SetTitle("# tower per event")


    NB_RECO_HSCP_vs_NBTOWER_isovar03 = TH2D("reco_hscp_vs_tower_iso_03","nb hscp vs nb tower per event with iso var < 0.3", 15,0,15,25,0,25)
    NB_RECO_HSCP_vs_NBTOWER_isovar03.GetXaxis().SetTitle("# hscp")
    NB_RECO_HSCP_vs_NBTOWER_isovar03.GetYaxis().SetTitle("# tower per event with iso var < 0.3")

    NB_RECO_cleanHSCP_vs_NBTOWER_isovar03 = TH2D("clean_reco_hscp_vs_tower_iso_03","nb clean hscp (pt > 50) vs nb tower per event with iso var < 0.3", 15,0,15,25,0,25)
    NB_RECO_cleanHSCP_vs_NBTOWER_isovar03.GetXaxis().SetTitle("# clean reco hscp")
    NB_RECO_cleanHSCP_vs_NBTOWER_isovar03.GetYaxis().SetTitle("# tower per event with iso var < 0.3")

    NB_RECO_HSCP_vs_NBTOWER_isovar05 = TH2D("reco_hscp_vs_tower_iso_05","nb hscp vs nb tower per event with iso var < 0.5", 15,0,15,25,0,25)
    NB_RECO_HSCP_vs_NBTOWER_isovar05.GetXaxis().SetTitle("# clean reco hscp")
    NB_RECO_HSCP_vs_NBTOWER_isovar05.GetYaxis().SetTitle("# tower per event with iso var < 0.5")

    NB_RECO_cleanHSCP_vs_NBTOWER_isovar05 = TH2D("clean_reco_hscp_vs_tower_iso_05","nb clean reco hscp (pt > 50) vs nb tower per event with iso var < 0.5", 15,0,15,25,0,25)
    NB_RECO_cleanHSCP_vs_NBTOWER_isovar05.GetXaxis().SetTitle("# clean reco hscp")
    NB_RECO_cleanHSCP_vs_NBTOWER_isovar05.GetYaxis().SetTitle("# tower per event with iso var < 0.5")

    NB_RECO_HSCP_vs_NBTOWER_isovar07 = TH2D("reco_hscp_vs_tower_iso_07","nb hscp vs nb tower per event with iso var < 0.7", 15,0,15,25,0,25)
    NB_RECO_HSCP_vs_NBTOWER_isovar07.GetXaxis().SetTitle("# hscp")
    NB_RECO_HSCP_vs_NBTOWER_isovar07.GetYaxis().SetTitle("# tower per event with iso var < 0.7")


    NB_RECO_cleanHSCP_vs_NBTOWER_isovar07 = TH2D("clean_reco_hscp_vs_tower_iso_07","nb clean hscp (pt > 50) vs nb tower per event with iso var < 0.7", 15,0,15,25,0,25)
    NB_RECO_cleanHSCP_vs_NBTOWER_isovar07.GetXaxis().SetTitle("# clean hscp")
    NB_RECO_cleanHSCP_vs_NBTOWER_isovar07.GetYaxis().SetTitle("# tower per event with iso var < 0.7")


    nb_matched_hscp_vs_nb_hscp = TH2D("matched_hscp_vs_hscp","nb matched hscp vs nb hscp", 15,0,15,15,0,15)
    nb_matched_hscp_vs_nb_hscp.GetXaxis().SetTitle("# hscp inside 2 by 2")
    nb_matched_hscp_vs_nb_hscp.GetYaxis().SetTitle("# hscp")

    NB_non_phys_vs_nb_tower_passing = TH2D("nb_non_phys_vs_tower","nb seed > 3 neighbours vs nb tower > threshold", 15,0,15,25,0,25)
    NB_non_phys_vs_nb_tower_passing.GetXaxis().SetTitle("# non physical event")
    NB_non_phys_vs_nb_tower_passing.GetYaxis().SetTitle("# tower > threshold")

    NB_match_tower_vs_nb_tower = TH2D("hscp_matched_iso_vs_tower","nb hscp matched and iso tower, vs nb tower per event", 10,0,10,10,0,10)



    NB_iso_below_07_per_event = TH1F("nb_below_iso_07_per_event","nb seed < iso 0.7 per event", 6,0,6)

    NB_iso_below_05_per_event = TH1F("nb_below_iso_05_per_event","nb seed < iso 0.5 per event", 6,0,6)
    NB_iso_below_03_per_event = TH1F("nb_below_iso_03_per_event","nb seed < iso 0.3 per event", 6,0,6)
    NB_iso_below_015_per_event = TH1F("nb_below_iso_015_per_event","nb seed < iso 0.15 per event", 6,0,6)


    NB_iso_below_vs_nb_seed = TH2D("nb_below_iso_vs_nb_seed","nb seed < iso 0.7 vs nb seed", 8,0,8,15,0,15)
    NB_iso_below_vs_nb_seed.GetXaxis().SetTitle("# seed IsoVar < 0.7")
    NB_iso_below_vs_nb_seed.GetYaxis().SetTitle("# seed")


    NB_iso_below_vs_nb_gen_hscp = TH2D("nb_below_iso_vs_nb_gen_hscp","nb seed < iso 0.7 vs nb hscp generator", 8,0,8,4,0,4)
    NB_iso_below_vs_nb_gen_hscp.GetXaxis().SetTitle("# seed IsoVar < 0.7")
    NB_iso_below_vs_nb_gen_hscp.GetYaxis().SetTitle("# gen hscp")


    Number_CH_GEN_HSCP_per_event = TH1F("nb_charged_generator_HSCP_per_event","nb generator HSCP per event",5,0,5)
    Number_ALL_GEN_HSCP_per_event = TH1F("nb_all_generator_HSCP_per_event","nb generator HSCP per event",5,0,5)

    Number_RECO_HSCP_per_event = TH1F("nb_HSCP_per_event","nb HSCP per event",15,0,15)
    Number_RECO_HSCP_cleaned_per_event = TH1F("nb_HSCP_per_event_post_cuts","nb HSCP per event after cuts on iso and pt",15,0,15)

    nb_towers_cut_per_Event = TH1F('number_of_towers_passing_cuts','tower nb passing thresholds',15,0,15)
    
    ratio_ecal_hcal_seeds = TH1F('ECAL_over_HCAL_final_seeds','ratio ecal/hcal for final seeds',40,0,0.4)

    Eff_Nseeds_over_nseedshscp = TH1F('Nseeds_over_Nseedhscp','efficiency # of seeds / # of seeds matching with an hscp',100,0,100)

    Xi2_choice_3 = TH1F("Xi2_3_in_line","smallest Xi2 when 3 neighbours are in line",100,0,0.5)
    Xi2_choice_2_and_1 = TH1F("Xi2_2_and_1_separate","smallest Xi2 when 2 together and 1 alone",100,0,0.5)

    Xi2_choice_2_and_1_hscp_iso = TH1F("HSCP_ISO_Xi2_2_and_1_separate","smallest Xi2 when 2 together and 1 alone with HSCP matching and iso",100,0,0.5)
    Xi2_choice_3_hscp_iso = TH1F ("HSCP_ISO_Xi2_3_in_line","smallest Xi2 when 3 neighbours are in line with hscp matching and iso",100,0,0.5)


    dr_min_matched_iso_hscp_tower = TH1F('min_dr_matched_gen_hscp_tower_iso','min dr between an hscp and a tower with iso var < 0.7',40,0,0.4)
    dr_min_iso_hscp_tower = TH1F('min_dr_hscp_tower_iso','min dr between an hscp and a tower with iso var < 0.7',1000,0,10)

    dr_min_matched_iso_reco_hscp_tower = TH1F('min_dr_matched_reco_hscp_tower_iso','min dr between a reco hscp and a tower with iso var < 0.7',40,0,0.4)
    
    dr_min_TEST_iso_hscp_tower = TH1F('min_dr_seed_gen_hscp_tower','min dr between all seeds and HSCPS',120,0,6)
    #============ END CALO TOWER HISTOS ===========
  
    #=============== START MET HISTOS ================
    eff_metfilter_vs_calo_met = TEfficiency("eff","efficiency of hltMET90 filter;hlt calo met;hltMET90;#epsilon",100,0,800)

    calomet_gen_hscp_raw = TH1F('caloMET_gen_hscp_raw','hlt Calo MET for events with at least 1 gen hscp - seed', 500,0,500)
    calomet_reco_hscp_raw = TH1F('caloMET_reco_hscp_raw','hlt Calo MET for events with at least 1 reco hscp - seed', 500,0,500)

    calomet_gen_hscp_matched = TH1F('caloMET_gen_hscp_matched','hlt Calo MET for matched gen hscp - seed', 500,0,500)
    calomet_reco_hscp_matched = TH1F('caloMET_reco_hscp_matched','hlt Calo MET for matched reco hscp - seed', 500,0,500)

    #=============== END MET HISTOS =============


    #============ GENERAL HISTOS ===========

    p_over_m_hscp = TH1F('P_over_M_hscp','ratio p/m = βγ for HSCPs',100,0,3)


    #============ END GENERAL HISTOS ===========


    number_maps = int(sys.argv[6])

    number_neighbours_max = 9

    Towers_neighboor_maps = [None] * number_maps
    for i in range(1,number_maps):
        trf = str(i)
        name = "Neighboor_activity_" +str(min_mip) + " mip central" + trf
        Towers_neighboor_maps[i] = TH2D(name,"eta-phi map tower activity",radius_idx_map,0,radius_idx_map,radius_idx_map,0,radius_idx_map)
        Towers_neighboor_maps[i].GetXaxis().SetTitle("Eta")
        Towers_neighboor_maps[i].GetYaxis().SetTitle("Phi")
    #================ END HISTOS ===============


    hlt_process = "HLT" #replace this with the HLT process name you want

    std_products = []
    #hlt products
    add_product(std_products,"trig_sum","trigger::TriggerEvent",f"hltTriggerSummaryAOD::{hlt_process}")
    add_product(std_products,"trig_res","edm::TriggerResults",f"TriggerResults::{hlt_process}")

    add_product(std_products,"l1algblk","BXVector<GlobalAlgBlk>","gtStage2Digis")
    add_product(std_products,"l1extblk","BXVector<GlobalExtBlk>","gtStage2Digis")

    #l1 objects
    add_product(std_products,"l1eg","BXVector<l1t::EGamma>","hltGtStage2Digis:EGamma")
    add_product(std_products,"l1etsum","BXVector<l1t::EtSum>","hltGtStage2Digis:EtSum")
    add_product(std_products,"l1jet","BXVector<l1t::Jet>","hltGtStage2Digis:Jet")
    add_product(std_products,"l1muon","BXVector<l1t::Muon>","gtStage2Digis:Muon")
    add_product(std_products,"l1egamma","BXVector<l1t::EGamma>","hltGtStage2Digis:EGamma")
    add_product(std_products,"l1tau","BXVector<l1t::Tau>","hltGtStage2Digis:Tau")
    # additional HLT products saved
    add_product(std_products,"calotower","edm::SortedCollection<CaloTower,edm::StrictWeakOrdering<CaloTower>>",f"hltTowerMakerForAll::{hlt_process}")
    add_product(std_products,"hscp","std::vector<susybsm::HSCParticle>","HSCParticleProducer")
    add_product(std_products,"hscpIso","edm::ValueMap<susybsm::HSCPIsolation>","HSCPIsolation:R03:HSCPAnalysis")
    add_product(std_products,"GenP","std::vector<reco::GenParticle>","genParticlesSkimmed")
    add_product(std_products,"tracks","std::vector<reco::Track>","hltMergedTracks")
    add_product(std_products,"dedx","edm::ValueMap<reco::DeDxData>","hltDeDxEstimatorProducer")

    add_product(std_products,"genPmc","vector<reco::GenParticle>","genParticles::SIM")
    add_product(std_products,"CaloMET","vector<reco::CaloMET>",f"hltMet::{hlt_process}")


    evtdata = EvtData(std_products,verbose=True)

    events = Events(CoreTools.get_filenames(args.in_filenames,args.prefix))


    l1menu_tree_file = ROOT.TFile.Open(args.l1menufile)
    l1menu_tree = ROOT.L1TUtmTriggerMenuRcd
    l1menu_tree.GetEntry(0)
    l1menu = l1menu_tree.L1TUtmTriggerMenu__

    l1name_to_indx = {x.second.getName() :x.second.getIndex() for x in l1menu.getAlgorithmMap()}

    #l1bitnr = l1name_to_indx["L1_SingleMu22"]

    ListBit = [0] * SizeArrays

    for x in range(0,SizeArrays):
        ListBit[x] = l1name_to_indx[ListNames[x]]

    nbevents = TH1F('nbevent', '# events',30000,20000,50000)

    print("This is a custom python script that allows the study of different aspects : L1 seed efficiencies, saving additional information about MET/MHT, energy clusters in ECAL/HCAL and more\n")
    print("Number of events to study : ",events.size())
    nbevents.Fill(events.size())

    nb_hscp_pre_iso, nb_hscp_after_iso = 0,0

    nb_cdt_2by2 , nb_cdt_out_2by2 , nb_cdt_2by2_hscp_iso , nb_cdt_out_2by2_hscp_iso = [0] * 4 , [0] * 4 , [0] * 4 , [0] * 4

    nb_seed_iso_matched , nb_seed_iso , nb_seed_iso_matched_reco = 0 , 0 , 0

    nb_evt_min_1_reco_hscp_presel = 0
    nb_evt_min_1_reco_hscp_raw = 0

    min_1_seed = 0
    
    totnb_tower_evt , totnb_tower_evt_hscp_iso , nb_ch_hscp , nb_ch_reco_hscp , nb_ch_clean_hscp , nb_ch_clean_reco_hscp , nb_matched_seed , nb_matched_reco_seed = 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0
    tot_nb_tower_evt_gen_min_1 , tot_nb_tower_evt_reco_min_1 = 0 , 0

    totnb_tower_evt_atleast_1_gen_hscp,totnb_tower_evt_atleast_1_reco_hscp = 0 , 0

    nb_belo_iso_tab_gen_tot , nb_belo_iso_tab_reco_tot = [0] * 5 , [0] * 5

    nb_tower_evt_cut , nb_belo_iso_tab_tot = [0] * 6 , [0] * 5
 
    denom_raw_evt=0
 
    or_met_ct_num_reco , met_filter_num , ct_filter_num , and_met_ct_num_reco , denom_presel_evt, denom_presel_evt_L1  = 0 , 0 , 0 , 0 , 0 , 0
    met_filter_num_L1 , ct_filter_num_l1 = 0,0
    met_filter_num_presel,met_filter_denom_presel = 0,0
    met_filter_denom_l1_raw,met_filter_num_l1_raw = 0,0
    or_met_ct_seed_num_reco , and_met_ct_seed_num_reco = 0,0 
    or_met_ct_seed_num_reco_nol1 , and_met_ct_seed_num_reco_nol1 = 0,0 
    or_met_ct_num_reco_nol1 , and_met_ct_num_reco_nol1 = 0,0
   
    or_first_l1_pfmet_num, or_first_l1_pfmet_denom = 0,0
    or_first_l1_pfmet_num_1_hscp_pe, or_first_l1_pfmet_denom_1_hscp_pe = 0,0
    or_first_l1_pfmet_num_presel, or_first_l1_pfmet_denom_presel = 0,0

 
    seed_filter_num = 0 
    seed_filter_num_l1 = 0
    for eventnr,event in enumerate(events):#events
        evtdata.get_handles(events)

        OR_L1_PFMET ,bool_min_1_seed ,min_1_seed_evt , seed_bool =  False, False, True, True
        OR_L1_PFMET_PRESEL = False
        nb_fake_evt , nb_tower_evt , nb_tower_evt_gen_min_1 , nb_tower_evt_reco_min_1 = 0 , 0 , 0 , 0
        nb_matched_pe, nb_matched_reco_pe = 0,0

        #we need to initialise the handles, must be called first for every event
        
        l1algblk = evtdata.get("l1algblk")

        l1dec_final = l1algblk.at(0,0).getAlgoDecisionFinal()
        l1dec_final_bxp1 = l1algblk.at(1,0).getAlgoDecisionFinal() if l1algblk.getLastBX()>=1 else None


        '''
        if eventnr == 10000:
            break
        '''


        '''
        l1muons_allbx = evtdata.get("l1muon")
        
        l1muons_bx0 = [l1muons_allbx.at(0,munr) for munr in range(0,l1muons_allbx.size(0))]
        l1muons_bx1 = [l1muons_allbx.at(1,munr) for munr in range(0,l1muons_allbx.size(1))] if l1muons_allbx.getLastBX()>=1 else []
        
        print(f"nr muons bx 1 {len(l1muons_bx1)}")
        for mu in l1muons_bx1:
            print(f" pt/eta/phi {mu.hwPt()} {mu.hwEta()} {mu.hwPhi()}")
        '''




        

        
        #print("Event # ", eventnr)
        if eventnr%1000 == 0 :
            print((eventnr/events.size())*100," %")


        hscparticle , hscp_iso , GenParticles , CaloTowers , GenParticlesMC , CaloMET = evtdata.get("hscp") , evtdata.get("hscpIso") , evtdata.get("GenP") , evtdata.get("calotower") , evtdata.get("genPmc"), evtdata.get("CaloMET")

        hltcalomet_objs = TrigTools.get_objs_passing_filter_aod(evtdata,"hltMET90")
        #print("Number of objects passing fitler hltMet90 :", len(hltcalomet_objs))
        #print("HLT Calo MET (from collection reco::calomet) = ", CaloMET[0].et())
        pass_hltmetfilter = False
        if len(hltcalomet_objs) != 0:
            pass_hltmetfilter = True

        eff_metfilter_vs_calo_met.Fill(pass_hltmetfilter,CaloMET[0].et())
        ''' 
        print("hlt calo met for this event :", CaloMET[0].et())
        for i in range(SizeArrays):
            print("L1 seed ", ListNames[i], " has bool = ", l1dec_final[ListBit[i]])
        '''
         
        #--------------- MET STUDY -------------------
        '''
        if CaloMET is not None:
            print("Size of CaloMET collection for this event :", len(CaloMET))
            for m in range(len(CaloMET)):
                print(CaloMET[m].et())

        '''

        #------------------ END MET STUDY -------------------


        all_pdg_id = []
        generator_hscp_ch,generator_hscp_clean_ch,generator_hscp_all,reco_hscp_clean_ch = 0,0,0,0

        FillPdgIds(GenParticlesMC,all_pdg_id)

        generator_hscp_ch , generator_hscp_all = FindChargedHSCP(all_pdg_id,idx_pdg_ch) , FindAllHSCP(all_pdg_id,idx_pdg_all)

        nb_ch_hscp += generator_hscp_ch

        Number_ALL_GEN_HSCP_per_event.Fill(generator_hscp_all)
        Number_CH_GEN_HSCP_per_event.Fill(generator_hscp_ch)
 
        indices = [i for i, x in enumerate(all_pdg_id) if x in idx_pdg_ch]       
        GEN_HSCPVector = []
        for p in range(len(indices)):
            p_over_m_hscp.Fill(GenParticlesMC[indices[p]].p()/1800)
            if GenParticlesMC[indices[p]].pt() > 55 and abs(GenParticlesMC[indices[p]].eta() ) < 2.1:                       
                GEN_HSCPVector.append((indices[p],GenParticlesMC[indices[p]].phi(),GenParticlesMC[indices[p]].eta()))

        nb_ch_clean_hscp += len(GEN_HSCPVector)
        generator_hscp_clean_ch = len(GEN_HSCPVector)

        tst_nb_dr01 = 0
        if reverse_bool:             
            for m in range(len(GEN_HSCPVector)):
                if CaloTowers is not None:
                    for i in range(CaloTowers.size()):
                        dr2 = deltaR(deltaR2(CaloTowers[i].eta(),CaloTowers[i].phi(),GEN_HSCPVector[m][2],GEN_HSCPVector[m][1]))
                        if dr2 < 0.1:
                            #print("Found 1 matching tower with hscp dr < 0.1 : em energy = ","%.2f" % CaloTowers[i].emEnergy(), " ,had energy =  ","%.2f" % CaloTowers[i].hadEnergy(), ", sum = ","%.2f" % (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy()))
                            tst_nb_dr01+=1
                            E_Ecal_dr01.Fill(CaloTowers[i].emEnergy())
                            E_Hcal_dr01.Fill(CaloTowers[i].hadEnergy())
                            E_Sum_dr01.Fill(CaloTowers[i].emEnergy() + CaloTowers[i].hadEnergy())
                            if CaloTowers[i].hadEnergy() != 0:
                                E_Ratio_dr01.Fill(CaloTowers[i].emEnergy()/CaloTowers[i].hadEnergy())
                            else:
                                E_Ratio_dr01.Fill(-5)
        
                Nb_ngh_dr01.Fill(tst_nb_dr01)



        all_pdg_id.clear()

        CaloVector = []
        idx_cvec = 0
        if CaloTowers is not None:
            nb_tower_per_event.Fill(CaloTowers.size())
            if ask_cut:
                for i in range(CaloTowers.size()):
                    if abs(CaloTowers[i].eta()) < 2.4:
                        CaloVector.append((idx_cvec,CaloTowers[i].iphi(),CaloTowers[i].ieta(),CaloTowers[i].phi(),CaloTowers[i].eta(),CaloTowers[i].emEnergy(),CaloTowers[i].hadEnergy()))
                        idx_cvec+=1
            else:
                for i in range(CaloTowers.size()):
                    CaloVector.append((idx_cvec,CaloTowers[i].iphi(),CaloTowers[i].ieta(),CaloTowers[i].phi(),CaloTowers[i].eta(),CaloTowers[i].emEnergy(),CaloTowers[i].hadEnergy()))
                    idx_cvec+=1

        else:
            nb_tower_non_valid_hist.Fill(1)

        CaloVectorAll = [] 
        CaloVectorAll = CaloVector[:]

        Number_RECO_HSCP_per_event.Fill(hscparticle.size())
        
        # RECO HSCP CANDIDATES CLEANING FOR EFFICIENCY N+1
        HSCPVector_raw = []
        HSCPVector = []
        for i in  range(hscparticle.size()):
            if hscparticle[i].trackRef().isNonnull():
                track_hscp = hscparticle[i].trackRef()
                iso_hscp = hscp_iso.get(track_hscp.key())
                HSCPVector_raw.append((i,track_hscp.phi(),track_hscp.eta()))
                if track_hscp.pt() > 55 and (iso_hscp.Get_ECAL_Energy() + iso_hscp.Get_HCAL_Energy())/track_hscp.p() < 0.3 and track_hscp.numberOfValidHits() > 8 and track_hscp.quality(2) and track_hscp.dxy() < 0.5 and track_hscp.dz() < 0.5 and abs(track_hscp.eta()) < 2.1 :
                    HSCPVector.append((i,track_hscp.phi(),track_hscp.eta()))
                

        Number_RECO_HSCP_cleaned_per_event.Fill(len(HSCPVector))

        reco_hscp_clean_ch = len(HSCPVector)
        nb_ch_clean_reco_hscp += len(HSCPVector)


        #END OF RECO HSCP CANDIDATES CLEANING
        or_first_l1_pfmet_denom+=1
        if (l1dec_final[ListBit[3]] == True or l1dec_final[ListBit[5]] == True or l1dec_final[ListBit[7]] == True or l1dec_final[ListBit[10]] == True or l1dec_final[ListBit[11]] == True):
            OR_L1_PFMET = True
            or_first_l1_pfmet_num +=1

        if len(HSCPVector_raw) !=0:
            or_first_l1_pfmet_num_1_hscp_pe +=1
            if (l1dec_final[ListBit[3]] == True or l1dec_final[ListBit[5]] == True or l1dec_final[ListBit[7]] == True or l1dec_final[ListBit[10]] == True or l1dec_final[ListBit[11]] == True):
                or_first_l1_pfmet_num_1_hscp_pe +=1

        if len(HSCPVector) !=0:
            or_first_l1_pfmet_denom_presel +=1
            if (l1dec_final[ListBit[3]] == True or l1dec_final[ListBit[5]] == True or l1dec_final[ListBit[7]] == True or l1dec_final[ListBit[10]] == True or l1dec_final[ListBit[11]] == True):
                or_first_l1_pfmet_num_presel +=1
                OR_L1_PFMET_PRESEL=True


        nb_hscp_evt=0


        if len(HSCPVector_raw) !=0:
           if CaloMET is not None:
               calomet_reco_hscp_raw.Fill(CaloMET[0].et()) 
        

        if len(GEN_HSCPVector) != 0:
            gen_vec = True
        else:
            gen_vec = False

        HSCPVector_eff = []
        HSCPVector_eff = HSCPVector[:]

        if len(HSCPVector) != 0:
            reco_vec = True
        else:
            reco_vec = False

        CT_matched = False
        if CaloVectorAll:
            nb_non_phys = 0
            nb_belo_iso_tab = [0] * 5
            nb_belo_iso_tab_gen , nb_belo_iso_tab_reco = [0] * 5 , [0] * 5

            for l in range(len(CaloVectorAll)-1,-1,-1):
                FillStepByStep(nb_tower_evt_cut,CaloVectorAll[l][5],CaloVectorAll[l][6])
                trf_eta,trf_phi,baryct_eta,baryct_phi, sum_ngh_above = 0,0,0,0,0

                cpt = True
                if CaloVectorAll[l][6] > 0:
                    p_all = PassThreshold("any",min_mip,CaloVectorAll[l][5],CaloVectorAll[l][6],max_sum_mip)
                    if p_all:
                        if bool_min_1_seed==False:
                            min_1_seed+=1
                            bool_min_1_seed = True
                        nb_tower_evt , totnb_tower_evt = nb_tower_evt+1, totnb_tower_evt+1
                        if gen_vec:
                            nb_tower_evt_gen_min_1 , tot_nb_tower_evt_gen_min_1  = nb_tower_evt_gen_min_1 + 1 , tot_nb_tower_evt_gen_min_1 + 1
                        if reco_vec:
                            tot_nb_tower_evt_reco_min_1 , nb_tower_evt_reco_min_1 = tot_nb_tower_evt_reco_min_1 + 1 , nb_tower_evt_reco_min_1 + 1

                        ratio_ecal_hcal_seeds.Fill(CaloVectorAll[l][5] / CaloVectorAll[l][6])

                        idx_phi,idx_eta,emEnergy,hadenergy,phi,eta,index = CaloVectorAll[l][1], CaloVectorAll[l][2],CaloVectorAll[l][5],CaloVectorAll[l][6],CaloVectorAll[l][3],CaloVectorAll[l][4],l
                        new_idx_phi,new_idx_eta,new_ratio_energy,new_phi,new_eta,new_index = FindTrueSeed(CaloVectorAll,idx_phi,idx_eta,(emEnergy/hadenergy),"both",phi,eta,index,max_sum_mip,min_mip)
                        all_nb_shift=0
                        while ((new_idx_phi != idx_phi) or (new_idx_eta != idx_eta)):
                            idx_phi,idx_eta,phi,eta,index = new_idx_phi,new_idx_eta,new_phi,new_eta,new_index
                            new_idx_phi,new_idx_eta,new_ratio_energy,new_phi,new_eta,new_index = FindTrueSeed(CaloVectorAll,idx_phi,idx_eta,new_ratio_energy,"both",new_phi,new_eta,new_index,max_sum_mip,min_mip)
                            all_nb_shift +=1

    
                        id_list,nb_all_ngh,nb_ngh_below,nb_ngh_above,sum_all_ngh_below = Create3by3MatrixLoop(CaloVectorAll,idx_phi,idx_eta,'both',min_mip_ngh,max_sum_mip)
                        sq_oi = 5
                        nb_towers_8_above_trh.Fill(nb_ngh_above)
                        Number_neighbours_per_seed.Fill((nb_ngh_above+nb_ngh_below))
                        Number_shift_seedOI.Fill(all_nb_shift)    
                        
                          
                        CaloVectorAll[new_index] = ((-999999,-999999,-999999,-999999,-999999,0,0))
                        '''
                        print("MAIN SEED has eta = ", eta, " and phi = ", phi,"\n")
                        print("This seed has  ", len(id_list), " neighbours ")
                        
                        for j in range (len(id_list)):
                            print("Candidate nb :",j, " has phi = ", id_list[j][3] , " and eta : ", id_list[j][4])
                        '''
 
                        if len(id_list) == 0:
                            nb_cdt_2by2[0] += 1
                            baryct_eta,baryct_phi = new_eta,new_phi
                            #print("0 neighbours, central seed has phi = ", new_phi, " and eta = ", new_eta)

                        elif len(id_list) == 1 or len(id_list) == 2 or len(id_list) == 3:
                            sq_oi = IsInMatrix(id_list,idx_phi,idx_eta)
                            if sq_oi == 1:
                                #print("We have ", len(id_list), " neighbours for that seed : ")
                                for u in range(len(id_list)):
                                   #print("Neighbour ",u+1, " -> eta = ", "%.2f" % id_list[u][4], " , phi = ", "%.2f" % id_list[u][3])
                                   trf_eta += id_list[u][4] * (id_list[u][5]+id_list[u][6])
                                   trf_phi += id_list[u][3] * (id_list[u][5]+id_list[u][6])
                                   sum_ngh_above += (id_list[u][5]+id_list[u][6]) 

                                baryct_eta,baryct_phi = (trf_eta/sum_ngh_above),(trf_phi/sum_ngh_above)
                                #print("Barycenter eta = ","%.2f" % baryct_eta, " and barycenter phi = ","%.2f" % baryct_phi)
                                
                                Nb_neighbour_within_2by2.Fill(len(id_list))
                                nb_cdt_2by2[len(id_list)] += 1
                                nb_towers_8_below_trh.Fill(nb_ngh_below)
                                nb_towers_8_all_trh.Fill(nb_ngh_below + nb_ngh_above)
                                for o in range(len(id_list)): 
                                    CaloVectorAll[id_list[o][0]] = ((-999999,-999999,-999999,-999999,-999999,0,0))
                                    #CaloVectorAll.remove((id_list[o]))
                         
                            else:
                                if len(id_list) == 1:
                                    sum_all_ngh_below += (id_list[0][5] + id_list[0][6])
                                    nb_ngh_below +=1

                                elif len(id_list) == 2:
                                    worst_cdt = FurthestInRatio(id_list,new_ratio_energy,1,idx_phi,idx_eta,Xi2_choice_2_and_1)
                                    best_cdt = FurthestInRatio(id_list,new_ratio_energy,2,idx_phi,idx_eta,Xi2_choice_2_and_1)
                                    nb_cdt_out_2by2[1] +=1
                                    for o in range(len(worst_cdt)):
                                        sum_all_ngh_below += (worst_cdt[o][5] + worst_cdt[o][6])
                                        nb_ngh_below +=1

                                    if len(best_cdt) != 0:
                                        CaloVectorAll[best_cdt[0][0]] = ((-999999,-999999,-999999,-999999,-999999,0,0))

                                elif len(id_list) == 3:
                                    worst_cdt_3 = FurthestInRatio(id_list,new_ratio_energy,1,idx_phi,idx_eta,Xi2_choice_3)
                                    best_cdt_3 = FurthestInRatio(id_list,new_ratio_energy,2,idx_phi,idx_eta,Xi2_choice_3)
                                    nb_cdt_out_2by2[2] +=1
                                    if len(worst_cdt_3) != 0:
                                        for i in range(len(worst_cdt_3)):
                                            sum_all_ngh_below += worst_cdt_3[i][5] + worst_cdt_3[i][6]
                                            nb_ngh_below +=1

                                    if len(best_cdt_3) != 0:
                                        for i in range(len(best_cdt_3)):
                                            #CaloVectorAll.remove((best_cdt_3[i]))
                                            CaloVectorAll[best_cdt_3[i][0]] = ((-999999,-999999,-999999,-999999,-999999,0,0))

                        else:
                            cpt = False
                            nb_non_phys +=1
                            Number_non_physical_event.Fill(1) 
   
                        if cpt:
                            if nb_ngh_above != 8:
                                Iso_var = sum_all_ngh_below/(8-nb_ngh_above)
                                #print("CaloTower ", l , " after ",all_nb_shift , " shifts")
                                #print("Computing iso var for seed",new_index ," id phi : ",idx_phi ," and id eta : ",idx_eta ," with", (8-nb_ngh_above), " neighbours outside thresholds, their sum = ","%.2f" % sum_all_ngh_below, " iso = ", Iso_var)
                                CountIsoNum(Iso_var,nb_belo_iso_tab)
                                CountIsoNum(Iso_var,nb_belo_iso_tab_tot)
                                if Iso_var < 1:
                                    nb_seed_iso += 1 
                                 
                                min_dr_nd = 999999
                                if gen_vec:
                                    totnb_tower_evt_atleast_1_gen_hscp += 1
                                    CountIsoNum(Iso_var,nb_belo_iso_tab_gen)
                                    CountIsoNum(Iso_var,nb_belo_iso_tab_gen_tot)
   
                                        
                                    for k in range(len(GEN_HSCPVector)):
                                        dr_nd = deltaR(deltaR2(new_eta,new_phi,GEN_HSCPVector[k][2],GEN_HSCPVector[k][1]))
                                        if dr_nd < min_dr_nd:
                                            min_dr_nd = dr_nd
                                            cdt_nd = GEN_HSCPVector[k]
    
                                    if min_dr_nd < 0.1:
                                        if CaloMET is not None:
                                            calomet_gen_hscp_matched.Fill(CaloMET[0].et())
                                        nb_matched_seed,nb_matched_pe,nb_seed_iso_matched = nb_matched_seed+1, nb_matched_pe+1, nb_seed_iso_matched+1
                                        dr_min_matched_iso_hscp_tower.Fill(min_dr_nd)
                                        GEN_HSCPVector.remove(cdt_nd)
                                        Sum_Energy_GEN_HSCP_8_neighbours_threshold.Fill(Iso_var)
    
                                    else:
                                        dr_min_iso_hscp_tower.Fill(min_dr_nd)
                                        Sum_Energy_8_neighbours_threshold.Fill(Iso_var)
                                    

                                min_dr_nd_reco = 999999
                                if reco_vec:
                                    totnb_tower_evt_atleast_1_reco_hscp += 1
                                    CountIsoNum(Iso_var,nb_belo_iso_tab_reco)
                                    CountIsoNum(Iso_var,nb_belo_iso_tab_reco_tot)
                                    for p in range(len(HSCPVector)):
                                        dr_nd_reco = deltaR(deltaR2(new_eta,new_phi,HSCPVector[p][2],HSCPVector[p][1]))
                                        if dr_nd_reco < min_dr_nd_reco:
                                            min_dr_nd_reco = dr_nd_reco
                                            cdt_nd_reco = HSCPVector[p]
 
                                    if min_dr_nd_reco < 0.1:
                                        if CaloMET is not None:
                                            calomet_reco_hscp_matched.Fill(CaloMET[0].et())

                                        CT_matched = True
                                        nb_matched_reco_seed,nb_matched_reco_pe,nb_seed_iso_matched_reco = nb_matched_reco_seed+1, nb_matched_reco_pe+1, nb_seed_iso_matched_reco+1
                                        dr_min_matched_iso_reco_hscp_tower.Fill(min_dr_nd_reco)
                                        HSCPVector.remove(cdt_nd_reco)
                                        Sum_Energy_8_neighbours_threshold_HSCP_RECO_ISO_matched.Fill(Iso_var)
    
                                    else:
                                        Sum_Energy_8_neighbours_threshold_HSCP_RECO_ISO.Fill(Iso_var)

       

                            else:
                                Sum_Energy_8_neighbours_threshold.Fill(-5)


 

            if generator_hscp_clean_ch == 2:
                nb_matched_seed_2_charged.Fill(nb_matched_pe)
            if generator_hscp_clean_ch == 1:
                nb_matched_seed_1_charged.Fill(nb_matched_pe)

            if reco_hscp_clean_ch == 2:
                nb_matched_seed_2_charged_reco.Fill(nb_matched_reco_pe)
            if reco_hscp_clean_ch == 1:
                nb_matched_seed_1_charged_reco.Fill(nb_matched_reco_pe)


            NB_iso_below_vs_nb_gen_hscp.Fill(nb_belo_iso_tab[3],generator_hscp_clean_ch)
            NB_iso_below_vs_nb_seed.Fill(nb_belo_iso_tab[3],nb_tower_evt)
            NB_iso_below_07_per_event.Fill(nb_belo_iso_tab[3])
            NB_iso_below_05_per_event.Fill(nb_belo_iso_tab[2])
            NB_iso_below_03_per_event.Fill(nb_belo_iso_tab[1])
            NB_iso_below_015_per_event.Fill(nb_belo_iso_tab[0])
    
            NB_non_phys_vs_nb_tower_passing.Fill(nb_non_phys,nb_tower_evt)
            CaloVectorAll.clear()
            #print("For Event nb ", eventnr, " there was ", nb_matched_pe , " matched seed for ", generator_hscp_ch, "charged hscp") 


            if reco_vec:
                NB_CLEAN_RECO_HSCP_vs_NBTOWER.Fill(reco_hscp_clean_ch,nb_tower_evt_reco_min_1)
                NB_RECO_HSCP_vs_matched_seed.Fill(reco_hscp_clean_ch,nb_matched_reco_pe)

                NB_RECO_cleanHSCP_vs_NBTOWER_isovar07.Fill(reco_hscp_clean_ch,nb_belo_iso_tab_reco[3])
                NB_RECO_cleanHSCP_vs_NBTOWER_isovar05.Fill(reco_hscp_clean_ch,nb_belo_iso_tab_reco[2])
                NB_RECO_cleanHSCP_vs_NBTOWER_isovar03.Fill(reco_hscp_clean_ch,nb_belo_iso_tab_reco[1])

                NB_RECO_HSCP_vs_NBTOWER_isovar03.Fill(hscparticle.size(),nb_belo_iso_tab_reco[1])
                NB_RECO_HSCP_vs_NBTOWER_isovar05.Fill(hscparticle.size(),nb_belo_iso_tab_reco[2])
                NB_RECO_HSCP_vs_NBTOWER_isovar07.Fill(hscparticle.size(),nb_belo_iso_tab_reco[3])
                NB_RECO_HSCP_vs_NBTOWER.Fill(hscparticle.size(),nb_tower_evt_reco_min_1)


            if gen_vec:
                NB_GEN_HSCP_vs_matched_seed.Fill(generator_hscp_clean_ch,nb_matched_pe)
                NB_GEN_HSCP_vs_NBTOWER.Fill(generator_hscp_clean_ch,nb_tower_evt_gen_min_1)

                NB_GEN_HSCP_vs_NBTOWER_isovar03.Fill(generator_hscp_clean_ch,nb_belo_iso_tab_gen[1])
                NB_GEN_HSCP_vs_NBTOWER_isovar05.Fill(generator_hscp_clean_ch,nb_belo_iso_tab_gen[2])
                NB_GEN_HSCP_vs_NBTOWER_isovar07.Fill(generator_hscp_clean_ch,nb_belo_iso_tab_gen[3])
                NB_GEN_HSCP_vs_NBTOWER_isovar09.Fill(generator_hscp_clean_ch,nb_belo_iso_tab_gen[4])

    
            nb_towers_cut_per_Event.Fill(nb_tower_evt)
        if len(HSCPVector_raw) != 0:
            nb_evt_min_1_reco_hscp_raw+=1
            denom_raw_evt+=1

            if CaloMET is not None:
                if OR_L1_PFMET:
                    met_filter_denom_l1_raw +=1
                    if len(hltcalomet_objs) != 0:
                    #if CaloMET[0].et() > 90:
                        met_filter_num_l1_raw += 1
                if len(hltcalomet_objs) != 0:         
                #if CaloMET[0].et() > 90:
                    met_filter_num+=1



        if len(HSCPVector_eff) !=0:
            nb_evt_min_1_reco_hscp_presel+=1 
            denom_presel_evt += 1  #Denominator is HSCP passing preselection AND first L1 seed
            #------------------------------------------------------------------------

            if OR_L1_PFMET:
                denom_presel_evt_L1 += 1
                if CT_matched:
                    ct_filter_num_l1+=1

            if bool_min_1_seed:
                seed_filter_num+=1
                if OR_L1_PFMET:
                    seed_filter_num_l1+=1

                if CaloMET is not None:
                    #if CaloMET[0].et() > 90:
                    if len(hltcalomet_objs) != 0:
                        and_met_ct_seed_num_reco_nol1+=1
                        if OR_L1_PFMET:
                            and_met_ct_seed_num_reco += 1

            if CT_matched:
                ct_filter_num+=1
                if CaloMET is not None:
                    #if CaloMET[0].et() > 90:
                    if len(hltcalomet_objs) != 0:
                        and_met_ct_num_reco_nol1+=1
                        if OR_L1_PFMET:
                            and_met_ct_num_reco+=1

            if CaloMET is None:
                if bool_min_1_seed:
                    or_met_ct_seed_num_reco_nol1+=1
                    if OR_L1_PFMET:
                        or_met_ct_seed_num_reco+=1
                       
                if CT_matched:
                    or_met_ct_num_reco_nol1+=1
                    if OR_L1_PFMET:
                        or_met_ct_num_reco+=1
            else:
                #if CaloMET[0].et() > 90: #Filter hltMet90 check
                if len(hltcalomet_objs) != 0:
                    met_filter_num_presel+=1
                    if OR_L1_PFMET:
                        met_filter_num_L1+=1

                if bool_min_1_seed or len(hltcalomet_objs) != 0: #CaloMET[0].et() > 90:
                    or_met_ct_seed_num_reco_nol1+=1
                    if OR_L1_PFMET:
                        or_met_ct_seed_num_reco+=1

                if CT_matched or len(hltcalomet_objs) != 0: #CaloMET[0].et() > 90:
                    or_met_ct_num_reco_nol1+=1
                    if OR_L1_PFMET:
                        or_met_ct_num_reco+=1


        GEN_HSCPVector.clear()




    # END OF LOOP ON EVENTS #

    for i in range(0,6):
        step_by_step_tower.SetBinContent(i+1,nb_tower_evt_cut[i])
        step_by_step_tower.GetXaxis().SetBinLabel(step_by_step_tower.GetXaxis().FindBin(i),str(ListCuts[i]))

    if nb_seed_iso: 
        Eff_Nseeds_over_nseedshscp.Fill((nb_seed_iso_matched/nb_seed_iso)*100)

    hfile.Write()
    print("We studied CaloTowers above ",min_mip, "MIP , and their neighbours above ", min_mip_ngh, " MIP , and sum (ECALL + HCAL) <",max_sum_mip, " MIP\n")


    err_eff_or_etmhf = math.sqrt(((or_first_l1_pfmet_num/or_first_l1_pfmet_denom) * ( 1 - (or_first_l1_pfmet_num/or_first_l1_pfmet_denom)))/or_first_l1_pfmet_denom)

    err_eff_or_etmhf_presel = math.sqrt(((or_first_l1_pfmet_num_presel/or_first_l1_pfmet_denom_presel) * ( 1 - (or_first_l1_pfmet_num_presel/or_first_l1_pfmet_denom_presel)))/or_first_l1_pfmet_denom_presel)
    print("ON ALL EVENTS WITH AT LEAST 1 RECO HSCP, NO PRESELECTION\n")

    print("Efficiency per event of L1_ETMHF100 OR L1_ETMHF110 OR L1_ETM150 OR L1_ETMHF120 OR L1_ETMHF150 ON ALL EVENTS= ", or_first_l1_pfmet_num, " / ", or_first_l1_pfmet_denom, " = ", (or_first_l1_pfmet_num/or_first_l1_pfmet_denom)*100, " % ± ", "%.2f" % (err_eff_or_etmhf*100) ,"% \n") 

    if denom_raw_evt !=0:
        eff_met = met_filter_num/denom_raw_evt
        err_eff_met = math.sqrt(((met_filter_num/denom_raw_evt) * ( 1 - (met_filter_num/denom_raw_evt)))/denom_raw_evt)
        print("hltMet90 filter efficiency = ",met_filter_num, " / ", denom_raw_evt, " = ","%.2f" % (eff_met*100) , " % ± ","%.2f" % (err_eff_met*100), " % \n")

    print("hltMet90 AFTER first L1 seed : ", met_filter_num_l1_raw, " / " , met_filter_denom_l1_raw, " = " , (met_filter_num_l1_raw/met_filter_denom_l1_raw)*100, " %")

    print("--------- EVENTS WITH 1 HSCP PASSING PRESELECTION ---------\n")
  
    print("Efficiency per event of L1_ETMHF100 OR L1_ETMHF110 OR L1_ETM150 OR L1_ETMHF120 OR L1_ETMHF150 PRESELECTION", or_first_l1_pfmet_num_presel, " / ", or_first_l1_pfmet_denom_presel, " = ", (or_first_l1_pfmet_num_presel/or_first_l1_pfmet_denom_presel)*100, " % ± ", "%.2f" % (err_eff_or_etmhf_presel*100) ,"% \n") 
    
    print("============ ON ALL HSCP PASSING PRESELECTION, NO L1 REQUIRED ============ \n")

    if denom_presel_evt != 0:
        eff_met_presel = met_filter_num_presel/denom_presel_evt
        err_eff_met_presel = math.sqrt(((met_filter_num_presel/denom_presel_evt) * ( 1 - (met_filter_num_presel/denom_presel_evt)))/denom_presel_evt)

        print("hltMet90 filter efficiency = ",met_filter_num_presel, " / ", denom_presel_evt, " = ","%.2f" % (eff_met_presel*100) , " % ± ","%.2f" % (err_eff_met_presel*100), " % \n")

        eff_ct_matched = ct_filter_num/denom_presel_evt
        err_eff_ct_matched = math.sqrt(((ct_filter_num/denom_presel_evt) * ( 1 - (ct_filter_num/denom_presel_evt)))/denom_presel_evt)
        print("CT matched filter efficiency = ",ct_filter_num, " / ", denom_presel_evt, " = ","%.2f" % (eff_ct_matched*100) , " % ± ","%.2f" % (err_eff_ct_matched*100) , " % \n")
        
        eff_seed_evt = seed_filter_num/denom_presel_evt
        err_eff_seed_evt = math.sqrt(((seed_filter_num/denom_presel_evt) * ( 1 - (seed_filter_num/denom_presel_evt)))/denom_presel_evt)
        print("CT 1 seed at least, no matching = ", seed_filter_num, " / ", denom_presel_evt, " = ","%.2f" % (eff_seed_evt*100), " % ± ","%.2f" % (err_eff_seed_evt*100)," % \n")

        print("------------ OR and AND efficiencies between hltMet90 and 1 seed at least per event ------------\n")

        
        print("hltMet90 || ct_01 matched = ",or_met_ct_seed_num_reco_nol1," / ",denom_presel_evt , " = ", (or_met_ct_seed_num_reco_nol1/denom_presel_evt)*100 , " % ± " , (math.sqrt((((or_met_ct_seed_num_reco_nol1/denom_presel_evt) * (1-(or_met_ct_seed_num_reco_nol1/denom_presel_evt)))/denom_presel_evt))*100) , " % \n")
    
        print("hltMet90 && ct_01 = ",and_met_ct_seed_num_reco_nol1," / ",denom_presel_evt , " = ", (and_met_ct_seed_num_reco_nol1/denom_presel_evt)*100 , " % ± " , math.sqrt((((and_met_ct_seed_num_reco_nol1/denom_presel_evt) * (1-(and_met_ct_seed_num_reco_nol1/denom_presel_evt)))/denom_presel_evt)))   
        print("\n")

        print("------------ OR and AND efficiencies between hltMet90 and matched CT ------------\n")
        print("\n")
    
        print("hltMet90 || ct_01 matched = ",or_met_ct_num_reco_nol1," / ",denom_presel_evt , " = ", (or_met_ct_num_reco_nol1/denom_presel_evt)*100 , " % ± " , math.sqrt((((or_met_ct_num_reco_nol1/denom_presel_evt) * (1-(or_met_ct_num_reco_nol1/denom_presel_evt)))/denom_presel_evt)))

        print("hltMet90 && ct_01 matched = ",and_met_ct_num_reco_nol1," / ",denom_presel_evt , " = ", (and_met_ct_num_reco_nol1/denom_presel_evt)*100 , " % ± " , math.sqrt((((and_met_ct_num_reco_nol1/denom_presel_evt) * (1-(and_met_ct_num_reco_nol1/denom_presel_evt)))/denom_presel_evt)))

        print("\n")
    else:
        print("Denominator = 0, efficiency = 0 %")

    print("============ ============ ============ ============ ============ ============\n")

    print("++++++++++++ ON ALL HSCP PASSING PRESELECTION + FIRST L1 FILTER REQUIRED ++++++++++++\n")

    if denom_presel_evt_L1 != 0:
        eff_met_l1 = met_filter_num_L1/denom_presel_evt_L1
        err_eff_met_l1 = math.sqrt(((met_filter_num_L1/denom_presel_evt_L1) * ( 1 - (met_filter_num_L1/denom_presel_evt_L1)))/denom_presel_evt_L1)
         
        print("hltMet90 filter efficiency = ",met_filter_num_L1, " / ", denom_presel_evt_L1, " = ","%.2f" % (eff_met_l1*100) , " % ± ","%.2f" % (err_eff_met_l1*100), " % \n")

        eff_ct_matched_l1 = ct_filter_num_l1/denom_presel_evt_L1
        err_eff_ct_matched_l1 = math.sqrt(((ct_filter_num_l1/denom_presel_evt_L1) * ( 1 - (ct_filter_num_l1/denom_presel_evt_L1)))/denom_presel_evt_L1)
        print("CT matched filter efficiency = ",ct_filter_num_l1, " / ", denom_presel_evt_L1, " = ","%.2f" % (eff_ct_matched_l1*100) , " % ± ","%.2f" % (err_eff_ct_matched_l1*100) , " \n")
        
        eff_seed_evt_l1 = seed_filter_num_l1/denom_presel_evt_L1
        err_eff_seed_evt_l1 = math.sqrt(((seed_filter_num_l1/denom_presel_evt_L1) * ( 1 - (seed_filter_num_l1/denom_presel_evt_L1)))/denom_presel_evt_L1)
        print("CT 1 seed at least, no matching = ", seed_filter_num_l1, " / ", denom_presel_evt_L1, " = ","%.2f" % (eff_seed_evt_l1*100), " % ± ","%.2f" % (err_eff_seed_evt_l1*100),"\n")

        print("------------ OR and AND efficiencies between hltMet90 and 1 seed at least per event ------------\n")

        print("hltMet90 || ct_01 matched = ",or_met_ct_seed_num_reco," / ",denom_presel_evt_L1 , " = ", (or_met_ct_seed_num_reco/denom_presel_evt_L1)*100 , " % ± " , math.sqrt((((or_met_ct_seed_num_reco/denom_presel_evt_L1) * (1-(or_met_ct_seed_num_reco/denom_presel_evt_L1)))/denom_presel_evt_L1)))
        print("\n")
    
        print("hltMet90 && ct_01 = ",and_met_ct_seed_num_reco," / ",denom_presel_evt_L1 , " = ", (and_met_ct_seed_num_reco/denom_presel_evt_L1)*100 , " % ± " , math.sqrt((((and_met_ct_seed_num_reco/denom_presel_evt_L1) * (1-(and_met_ct_seed_num_reco/denom_presel_evt_L1)))/denom_presel_evt_L1)))
        print("\n")

        print("------------ OR and AND efficiencies between hltMet90 and matched CT ------------\n")
    
        print("hltMet90 || ct_01 matched = ",or_met_ct_num_reco," / ",denom_presel_evt_L1 , " = ", (or_met_ct_num_reco/denom_presel_evt_L1)*100 , " % ± " , math.sqrt((((or_met_ct_num_reco/denom_presel_evt_L1) * (1-(or_met_ct_num_reco/denom_presel_evt_L1)))/denom_presel_evt_L1)))
        print("\n")
    
        print("hltMet90 && ct_01 = ",and_met_ct_num_reco," / ",denom_presel_evt_L1 , " = ", (and_met_ct_num_reco/denom_presel_evt_L1)*100 , " % ± " , math.sqrt((((and_met_ct_num_reco/denom_presel_evt_L1) * (1-(and_met_ct_num_reco/denom_presel_evt_L1)))/denom_presel_evt_L1)))
        print("\n")
    
    else:
        print("Denominator = 0, efficiency = 0 %")

    print("++++++++++++ ++++++++++++ ++++++++++++ ++++++++++++ ++++++++++++ ++++++++++++ \n")
    print("\n")

    print("Fraction nb event with 1 seed at least / nb event with 1 hscp passing presel at least = ",min_1_seed, " / " , nb_evt_min_1_reco_hscp_presel, " = ",(min_1_seed/nb_evt_min_1_reco_hscp_presel)*100 , " %" )
    print("\n")
   
    print("There were in total : ", totnb_tower_evt , " seeds passing thresholds \n")
    print("With cut IsoV < 0.1 :", nb_belo_iso_tab_tot[0])
    print("With cut IsoV < 0.3 :", nb_belo_iso_tab_tot[1])
    print("With cut IsoV < 0.5 :", nb_belo_iso_tab_tot[2])
    print("With cut IsoV < 0.7 :", nb_belo_iso_tab_tot[3])
    print("With cut IsoV < 0.9 :", nb_belo_iso_tab_tot[4])
    print("-----------------------------------------\n Events with at least 1 charged GEN HSCP \n")
    
    print("There was ", totnb_tower_evt_atleast_1_gen_hscp, " seeds with at least 1 charged gen HSCP")
    print("With cut IsoV < 0.1 :", nb_belo_iso_tab_gen_tot[0])
    print("With cut IsoV < 0.3 :", nb_belo_iso_tab_gen_tot[1])
    print("With cut IsoV < 0.5 :", nb_belo_iso_tab_gen_tot[2])
    print("With cut IsoV < 0.7 :", nb_belo_iso_tab_gen_tot[3])
    print("With cut IsoV < 0.9 :", nb_belo_iso_tab_gen_tot[4])
    print("\n") 
    print("-----------------------------------------\n Events with at least 1 charged RECO HSCP \n")
    print("There was ", totnb_tower_evt_atleast_1_reco_hscp, " seeds with at least 1 charged reco HSCP")
    print("With cut IsoV < 0.1 :", nb_belo_iso_tab_reco_tot[0])
    print("With cut IsoV < 0.3 :", nb_belo_iso_tab_reco_tot[1])
    print("With cut IsoV < 0.5 :", nb_belo_iso_tab_reco_tot[2])
    print("With cut IsoV < 0.7 :", nb_belo_iso_tab_reco_tot[3])
    print("With cut IsoV < 0.9 :", nb_belo_iso_tab_reco_tot[4])

    print("ON GENERATOR LEVEL : \n")
    print("----------------------------------------\n") 
    if totnb_tower_evt: 
        print("Efficiency (nb gen charged hscp cleaned/nb seed) = ", nb_ch_clean_hscp,"/",totnb_tower_evt, " = " ,(nb_ch_clean_hscp/totnb_tower_evt)*100, " %")
        print("Now efficiency (nb matched seed / nb seed) = ",nb_matched_seed, "/", totnb_tower_evt," = ", (nb_matched_seed/totnb_tower_evt)*100 , " %")
    if nb_ch_clean_hscp:
        print("Efficiency (nb matched seed/nb gen charged hscp cleaned) = ",nb_matched_seed, "/",nb_ch_clean_hscp, " = " ,(nb_matched_seed/nb_ch_clean_hscp)*100, " %")

    print("----------------------------------------\n") 
    print("ON RECO AFTER PRESEL LEVEL : \n")
    print("++++++++++++++++++++++++++++++++++++++++\n") 
 
    if totnb_tower_evt: 
        print("Efficiency (nb reco charged hscp cleaned/nb seed) = ", nb_ch_clean_reco_hscp,"/",totnb_tower_evt, " = " ,(nb_ch_clean_reco_hscp/totnb_tower_evt)*100, " %")
        print("Now efficiency (nb matched seed / nb seed) = ",nb_matched_reco_seed, "/", totnb_tower_evt," = ", (nb_matched_reco_seed/totnb_tower_evt)*100 , " %")
    if nb_ch_clean_reco_hscp:
        print("Efficiency (nb matched seed/nb reco charged hscp) = ",nb_matched_reco_seed, "/",nb_ch_clean_reco_hscp, " = " ,(nb_matched_reco_seed/nb_ch_clean_reco_hscp)*100, " %")

    print("++++++++++++++++++++++++++++++++++++++++\n") 
    print("ON ALL CALO TOWERS \n")
    print("--------------- DIRECT GOOD EVENTS ---------------")

    print(nb_cdt_2by2[0], " were solo seed, ", nb_cdt_2by2[1], " were seeds with 1 neighbours not in the diagonal", nb_cdt_2by2[2], " were seeds with 2 adjacent neighbours, and ",nb_cdt_2by2[3], " were seeds with 3 neighbours in an L shape (4 possibility since 4 corners in 3 by 3)")

    print("--------------- BAD EVENTS WHERE WE HAVE TO CHOOSE  ---------------")
    print(nb_cdt_out_2by2[1], " were seeds with 2 non-adjacent neighbours, pick the one with ratio ecal/hcal closer to initial seed, ", nb_cdt_out_2by2[2], " were seeds with 3 neighbours, where we have to chose (either 3 solo seeds, or 3 in a line or 2 adjacent and 1 not diagonal)\n")




    '''
    print("ON HSCP MATCHED + ISO CALO TOWERS \n")
    print("--------------- DIRECT GOOD EVENTS ---------------")

    print(nb_cdt_2by2_hscp_iso[0], " were solo seed, ", nb_cdt_2by2_hscp_iso[1], " were seeds with 1 neighbours not in the diagonal", nb_cdt_2by2_hscp_iso[2], " were seeds with 2 adjacent neighbours, and ",nb_cdt_2by2_hscp_iso[3], " were seeds with 3 neighbours in an L shape (4 possibility since 4 corners in 3 by 3)")

    print("--------------- BAD EVENTS WHERE WE HAVE TO CHOOSE  ---------------")
    print(nb_cdt_out_2by2_hscp_iso[1], " were seeds with 2 non-adjacent neighbours, pick the one with ratio ecal/hcal closer to initial seed, ", nb_cdt_out_2by2_hscp_iso[2], " were seeds with 3 neighbours, where we have to chose (either 3 solo seeds, or 3 in a line, or 2 adjacent and 1 not diagonal)")
    '''
    print("There was ",nb_hscp_pre_iso, " pair hscp-tower pre iso, and ", nb_hscp_after_iso, " pair hscp-tower aftert iso")
