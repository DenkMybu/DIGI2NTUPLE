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
import pandas as pd

from DataFormats.FWLite import Events, Handle
from Analysis.HLTAnalyserPy.EvtData import EvtData, EvtHandles, add_product,get_objs

from Analysis.HLTAnalyserPy.CaloVectors import SetBin,FillHist,ActivityTab,SetActivityAround,PassThreshold,CreateCaloMap,Get12NeighboursLoop,CheckHighEMatrix,deltaR2,deltaR

from Analysis.HLTAnalyserPy.CaloVectors import Create3by3MatrixLoop,CreateCaloTabLoop,FindTrueSeed,FindAdjacentPair,IsInMatrix,FurthestInRatio


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
    parser.add_argument('event to search for a muon-tower',default=1,help='in which event you seek the map')
    parser.add_argument('min nb of MIP in the ECAL',default=1,help='how many MIPS will the muon leave in the ECAL')
    parser.add_argument('number of maps you want',default=20,help='how many eta-phi maps you want')
    args = parser.parse_args()


    hfile = TFile( 'CaloStudiesHSCP_iso_vector.root', 'RECREATE', 'Output of hscp script to save HLT quantities' )
    EventToStudy = sys.argv[4]
    min_mip = int(sys.argv[5])

    #============ CALO TOWER HISTOS ===========

    radius_idx_map=5


    nb_tower_per_event = TH1F('nb_tower_per_event', '# calo tower ',1000,0,1000)
    nb_tower_non_valid_hist = TH1F('nb calo tower non valid (for whole event)', '# calo towers missing',100000,0,100000)
    #ratio_mu_ecalhcal = TH1F ( 'ratio_ecal_hcal_MU', 'E_ecal / E_hcal ',100,0,1)
    recoCaloMet = TEfficiency("eff","reco::caloMET efficiency;reco::caloMET;#epsilon",100,0,1000)
    nb_towers_8_above_trh = TH1F('number_of_neighbours_passing_cuts','tower nb passing cuts',9,0,9)
    nb_towers_8_below_trh = TH1F('number of neighbours (8 max) NOT passing cuts (same as above)','tower nb not passing cuts',9,0,9)
    nb_towers_8_all_trh = TH1F('number of neighbours per cell','# of neibhours per cell',9,0,9)
    Sum_Energy_8_neighbours_threshold = TH1F('Sum_8_neighbors_tower','(Sum of X cells (below threshold) around cell of interest / (8- # cells)',50,0,10)



    nb_towers_8_above_trh_HSCP = TH1F('HSCP_number_of_neighbours_passing_cuts','tower nb passing cuts',9,0,9)
    nb_towers_8_below_trh_HSCP = TH1F('HSCP number of neighbours (8 max) NOT passing cuts (same as above)','tower nb not passing cuts',9,0,9)
    nb_towers_8_all_trh_HSCP = TH1F('HSCP number of neighbours per cell','# of neibhours per cell',9,0,9)
    Sum_Energy_8_neighbours_threshold_HSCP = TH1F('HSCP_Sum_8_neighbors_tower','(Sum of X cells (below threshold) around cell of interest(dr hscp calo tower min) / (8- # cells)',50,0,10)

    nb_towers_8_above_trh_HSCP_ISO = TH1F('HSCP_ISO_number_of_neighbours_passing_cuts','tower nb passing cuts',9,0,9)
    nb_towers_8_below_trh_HSCP_ISO = TH1F('HSCP_ISO_number_of_neighbours_NOT_passing_cuts','tower nb not passing cuts',9,0,9)
    nb_towers_8_all_trh_HSCP_ISO = TH1F('HSCP_ISO_number_of_neighbours_per_cell','# of neibhours per cell',9,0,9)
    Sum_Energy_8_neighbours_threshold_HSCP_ISO = TH1F('HSCP_ISO_Sum_8_neighbors_tower','(Sum of X cells (below threshold) around cell of interest iso(dr hscp calo tower min) / (8- # cells)',50,0,10)


    Number_shift_seedOI = TH1F("Number_of_seed_shifts","# of shifts from initial seed to next one",10,0,10)
    HSCP_Number_shift_seedOI = TH1F("HSCP_Number_of_seed_shifts","# of shifts from initial seed to next one",10,0,10)
    ISO_HSCP_Number_shift_seedOI = TH1F("ISO_HSCP_Number_of_seed_shifts","# of shifts from initial seed to next one",10,0,10)

    HSCP_ISO_Number_non_physical_event = TH1F("Number_tower_more_3_ngh_hscp_iso","# of non-physical event",2,0,2)
    Number_non_physical_event = TH1F("Number_tower_more_3_ngh","# of non-physical event",2,0,2)
    

    Number_HSCP_2by2_max_event = TH1F('#_neighbours_within_2by2_max_event_matched_HSCP','# HSCP',10,0,10)
    Number_all_2by2_max_event = TH1F('#_neighbours_within_2by2_max_event_all','# HSCP',10,0,10)
 
    Number_fake_dr_sup_03 = TH1F('#_fake_tower','# fake tower (dr > 0.3 when passing all cuts)',10,0,10)

    NB_HSCP_vs_NBTOWER = TH2D("hscp_vs_tower","nb hscp vs nb tower per event", 15,0,15,25,0,25)
    NB_HSCP_vs_NBTOWER.GetXaxis().SetTitle("# hscp")
    NB_HSCP_vs_NBTOWER.GetZaxis().SetTitle("# tower per event")

    nb_matched_hscp_vs_nb_hscp = TH2D("matched_hscp_vs_hscp","nb matched hscp vs nb hscp", 15,0,15,15,0,15)
    nb_matched_hscp_vs_nb_hscp.GetXaxis().SetTitle("# hscp inside 2 by 2")
    nb_matched_hscp_vs_nb_hscp.GetYaxis().SetTitle("# hscp")

    NB_non_phys_vs_nb_tower_passing = TH2D("#_non_phys_vs_tower","nb seed > 3 neighbours vs nb tower > threshold", 15,0,15,25,0,25)
    NB_non_phys_vs_nb_tower_passing.GetXaxis().SetTitle("# non physical event")
    NB_non_phys_vs_nb_tower_passing.GetYaxis().SetTitle("# tower > threshold")

    NB_match_tower_vs_nb_tower = TH2D("hscp_matched_iso_vs_tower","nb hscp matched and iso tower, vs nb tower per event", 10,0,10,10,0,10)


    Number_HSCP_per_event = TH1F("#_HSCP_per_event","nb HSCP per event",15,0,15)
    Number_HSCP_CUTS_per_event = TH1F("#_HSCP_per_event_post_cuts","nb HSCP per event after cuts on iso and pt",15,0,15)
    nb_towers_cut_per_Event = TH1F('number_of_towers_passing_cuts','tower nb passing thresholds',50,0,50)
    
    ratio_ecal_hcal_seeds = TH1F('ECAL_over_HCAL_final_seeds','ratio ecal/hcal for final seeds',30,0,3)

    Eff_Nseeds_over_nseedshscp = TH1F('Nseeds_over_Nseedhscp','efficiency # of seeds / # of seeds matching with an hscp',100,0,100)

    Xi2_choice_3 = TH1F("Xi2_3_in_line","smallest Xi2 when 3 neighbours are in line",100,0,0.5)
    Xi2_choice_2_and_1 = TH1F("Xi2_2_and_1_separate","smallest Xi2 when 2 together and 1 alone",100,0,0.5)

    Xi2_choice_2_and_1_hscp_iso = TH1F("HSCP_ISO_Xi2_2_and_1_separate","smallest Xi2 when 2 together and 1 alone with HSCP matching and iso",100,0,0.5)
    Xi2_choice_3_hscp_iso = TH1F ("HSCP_ISO_Xi2_3_in_line","smallest Xi2 when 3 neighbours are in line with hscp matching and iso",100,0,0.5)
    #============ END CALO TOWER HISTOS ===========







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

    #l1 objects
    add_product(std_products,"l1eg","BXVector<l1t::EGamma>","hltGtStage2Digis:EGamma")
    add_product(std_products,"l1etsum","BXVector<l1t::EtSum>","hltGtStage2Digis:EtSum")
    add_product(std_products,"l1jet","BXVector<l1t::Jet>","hltGtStage2Digis:Jet")
    add_product(std_products,"l1muon","BXVector<l1t::Muon>","hltGtStage2Digis:Muon")
    add_product(std_products,"l1egamma","BXVector<l1t::EGamma>","hltGtStage2Digis:EGamma")
    add_product(std_products,"l1tau","BXVector<l1t::Tau>","hltGtStage2Digis:Tau")
    # additional HLT products saved
    add_product(std_products,"calotower","edm::SortedCollection<CaloTower,edm::StrictWeakOrdering<CaloTower>>",f"hltTowerMakerForAll::{hlt_process}")
    add_product(std_products,"hscp","std::vector<susybsm::HSCParticle>","HSCParticleProducer")
    add_product(std_products,"hscpIso","edm::ValueMap<susybsm::HSCPIsolation>","HSCPIsolation:R03:HSCPAnalysis")
    add_product(std_products,"GenP","std::vector<reco::GenParticle>","genParticlesSkimmed")
    add_product(std_products,"tracks","std::vector<reco::Track>","hltMergedTracks")
    add_product(std_products,"dedx","edm::ValueMap<reco::DeDxData>","hltDeDxEstimatorProducer")



    evtdata = EvtData(std_products,verbose=True)

    events = Events(CoreTools.get_filenames(args.in_filenames,args.prefix))
    nbevents = TH1F('nbevent', '# events',30000,20000,50000)
    print("This is a custom python script that allows the study of different aspects : L1 seed efficiencies, saving additional information about MET/MHT, energy clusters in ECAL/HCAL and more\n")
    print("Number of events to study : ",events.size())
    nbevents.Fill(events.size())

    print("Studying event ",EventToStudy)
    nb_hscp_pre_iso, nb_hscp_after_iso = 0,0

    nb_cdt_2by2,nb_cdt_out_2by2 = [0] * 4,[0] * 4
    nb_cdt_2by2_hscp_iso,nb_cdt_out_2by2_hscp_iso = [0] * 4,[0] * 4

    totnb_tower_evt,totnb_tower_evt_hscp_iso = 0,0
    for eventnr,event in enumerate(events):#events
        Pairs_idx_used_tower = []    
        nb_fake_evt,nb_tower_evt = 0,0
        nb_num_eff,nb_denom_eff = 0,0
        #we need to initialise the handles, must be called first for every event
        evtdata.get_handles(events)

        #print("Event # ", eventnr)
        if eventnr%1000 == 0 :
            print((eventnr/events.size())*100," %")

        #we can also do this for l1obj which the HLT filters write
        #these filters always begin with hltL1s
        l1met_objs = TrigTools.get_objs_passing_filter_aod(evtdata,"hltL1sAllETMHFSeeds")
        #tracks
        hlt_trk_p4s = TrigTools.get_objs_passing_filter_aod(evtdata,"hltTrk50Filter")
        #now the above is the basic HLT info you can expect to exist in AOD
        #ESTIMATING l1 SEED EFFICIENCIES
        l1singlemu18 = TrigTools.get_objs_passing_filter_aod(evtdata,"hltL1sSingleMu18")


        hscparticle = evtdata.get("hscp")
        hscp_iso = evtdata.get("hscpIso")
        GenParticles = evtdata.get("GenP")
        CaloTowers = evtdata.get("calotower")
        
        CaloVector = []
        CaloVectorAll = []    
        if CaloTowers is not None:
            nb_tower_per_event.Fill(CaloTowers.size())
            for i in range(CaloTowers.size()):
                CaloVector.append((i,CaloTowers[i].iphi(),CaloTowers[i].ieta(),CaloTowers[i].phi(),CaloTowers[i].eta(),CaloTowers[i].emEnergy(),CaloTowers[i].hadEnergy()))
        else:
            nb_tower_non_valid_hist.Fill(1)
    
        CaloVectorAll = CaloVector[:]


        Number_HSCP_per_event.Fill(hscparticle.size())

        HSCPVector = []
        for i in  range(hscparticle.size()):
            if hscparticle[i].trackRef().isNonnull():
                track_hscp = hscparticle[i].trackRef()
                iso_hscp = hscp_iso.get(track_hscp.key())
                if track_hscp.pt() > 50:
                    if iso_hscp.Get_TK_SumEt() < 50 and (iso_hscp.Get_ECAL_Energy() + iso_hscp.Get_HCAL_Energy())/track_hscp.p() < 0.3:
                        HSCPVector.append((i,track_hscp.phi(),track_hscp.eta()))
        #print("iso size : ",hscp_iso.size())

        Number_HSCP_CUTS_per_event.Fill(len(HSCPVector))
        #print("Size of HSCP collection : ", hscparticle.size())
        nb_hscp_evt=0
        if CaloVector:
            for i in range(hscparticle.size()):
                nb_tower_hscp_evt=0
                mindreg = 9999
                idx_dreg = 9999
                dreg = 9999
                if hscparticle[i].trackRef().isNonnull():
                    track = hscparticle[i].trackRef()
                    iso_got = hscp_iso.get(track.key())
                    #print("iso TK sumET : ", iso_got.Get_TK_SumEt())
                    if track.pt() > 50:
                        p_over_m_hscp.Fill(track.p()/1800)
                        for l in range(len(CaloVector)-1,0,-1):
                            if CaloVector[l][6] > 0:
                                if CaloVector[l][5] > 0.6 and CaloVector[l][6] > 2.5 and (CaloVector[l][5] / CaloVector[l][6]) < 0.3 and (CaloVector[l][5] + CaloVector[l][6]) < 15:
                                    dreg = deltaR(deltaR2(CaloVector[l][4],CaloVector[l][3],track.eta(),track.phi()))
                                    nb_tower_hscp_evt+=1
                                    if dreg < mindreg:
                                        mindreg = dreg
                                        idx_dreg = l
    
                        if mindreg < 0.3 and iso_got.Get_TK_SumEt() < 50 and (iso_got.Get_ECAL_Energy() + iso_got.Get_HCAL_Energy())/track.p() < 0.3:
                            nb_hscp_after_iso+=1
                            
                            idx_phi_hscp_iso,idx_eta_hscp_iso,emEnergy_hscp_iso,hadenergy_hscp_iso,phi_hscp_iso,eta_hscp_iso = CaloVector[idx_dreg][1], CaloVector[idx_dreg][2],CaloVector[idx_dreg][5],CaloVector[idx_dreg][6],CaloVector[idx_dreg][3],CaloVector[idx_dreg][4]
                            new_idx_phi_hscp_iso,new_idx_eta_hscp_iso,new_ratio_energy_hscp_iso,new_phi_hscp_iso,new_eta_hscp_iso = FindTrueSeed(CaloVector,idx_phi_hscp_iso,idx_eta_hscp_iso,(emEnergy_hscp_iso/hadenergy_hscp_iso),"both",phi_hscp_iso,eta_hscp_iso)
                            if new_idx_phi_hscp_iso != 999999:
                                ISO_HSCP_nb_shift=0
                                while ((new_idx_phi_hscp_iso != idx_phi_hscp_iso) or (new_idx_eta_hscp_iso != idx_eta_hscp_iso)):
                                    idx_phi_hscp_iso,idx_eta_hscp_iso,phi_hscp_iso,eta_hscp_iso = new_idx_phi_hscp_iso,new_idx_eta_hscp_iso,new_phi_hscp_iso,new_eta_hscp_iso
                                    new_idx_phi_hscp_iso,new_idx_eta_hscp_iso,new_ratio_energy_hscp_iso,new_phi_hscp_iso,new_eta_hscp_iso = FindTrueSeed(CaloVector,idx_phi_hscp_iso,idx_eta_hscp_iso,new_ratio_energy_hscp_iso,"both",new_phi_hscp_iso,new_eta_hscp_iso)

                                    ISO_HSCP_nb_shift +=1
    
    
                            if new_idx_phi_hscp_iso != 999999:
                                ISO_HSCP_Number_shift_seedOI.Fill(ISO_HSCP_nb_shift)        
                                id_list_hscp_iso,nb_all_ngh_hscp_iso,nb_ngh_below_hscp_iso,nb_ngh_above_hscp_iso,sum_all_ngh_below_hscp_iso = Create3by3MatrixLoop(CaloVector,idx_phi_hscp_iso,idx_eta_hscp_iso,'both',2)
                                sq_oi_hscp = 5

                                if len(id_list_hscp_iso) == 0:
                                    nb_cdt_2by2_hscp_iso[0] += 1
                                    nb_towers_8_below_trh_HSCP_ISO.Fill(nb_ngh_below_hscp_iso)
                                    nb_towers_8_all_trh_HSCP_ISO.Fill(nb_ngh_below_hscp_iso + nb_ngh_above_hscp_iso)
                                    if nb_ngh_below_hscp_iso != 8:
                                        Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(sum_all_ngh_below_hscp_iso/(8-nb_ngh_below_hscp_iso))
                                    else:
                                        Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(-1)

                                elif len(id_list_hscp_iso) == 1 or len(id_list_hscp_iso) == 2 or len(id_list_hscp_iso) == 3:
                                    nb_hscp_evt+=1
                                    Number_HSCP_2by2_max_event.Fill(len(id_list_hscp_iso)+1) 
                                    sq_oi_hscp = IsInMatrix(id_list_hscp_iso,idx_phi_hscp_iso,idx_eta_hscp_iso)
                                    if sq_oi_hscp == 1:
                                        nb_cdt_2by2_hscp_iso[len(id_list_hscp_iso)] +=1
                                        for o in range(len(id_list_hscp_iso)):
                                            CaloVector.remove((id_list_hscp_iso[o]))
    
                                        nb_towers_8_below_trh_HSCP_ISO.Fill(nb_ngh_below_hscp_iso)
                                        nb_towers_8_above_trh_HSCP_ISO.Fill(nb_ngh_above_hscp_iso)    
                                        nb_towers_8_all_trh_HSCP_ISO.Fill(nb_ngh_below_hscp_iso + nb_ngh_above_hscp_iso)
                                        if nb_ngh_below_hscp_iso != 8:
                                            Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(sum_all_ngh_below_hscp_iso/(8-nb_ngh_below_hscp_iso))
                                        else:
                                            Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(-1)
                                    else:
                                        if len(id_list_hscp_iso) == 1:
                                            sum_all_ngh_below_hscp_iso += (id_list_hscp_iso[0][5] + id_list_hscp_iso[0][6])
                                            nb_ngh_below_hscp_iso +=1
                                            if nb_ngh_below_hscp_iso != 8:
                                                Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(sum_all_ngh_below_hscp_iso/(8-nb_ngh_below_hscp_iso))
                                            else:
                                                Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(-1)

                                        elif len(id_list_hscp_iso) == 2:
                                            worst_cdt_hscp_iso = FurthestInRatio(id_list_hscp_iso,new_ratio_energy_hscp_iso,1,idx_phi_hscp_iso,idx_eta_hscp_iso,Xi2_choice_2_and_1_hscp_iso)
                                            best_cdt_hscp_iso = FurthestInRatio(id_list_hscp_iso,new_ratio_energy_hscp_iso,2,idx_phi_hscp_iso,idx_eta_hscp_iso,Xi2_choice_2_and_1_hscp_iso)
                                            nb_cdt_out_2by2_hscp_iso[1] +=1

                                            for o in range(len(worst_cdt_hscp_iso)):
                                                sum_all_ngh_below_hscp_iso += (worst_cdt_hscp_iso[o][5] + worst_cdt_hscp_iso[o][6])
                                                nb_ngh_below_hscp_iso +=1

                                            if nb_ngh_below_hscp_iso != 8:
                                                Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(sum_all_ngh_below_hscp_iso/(8-nb_ngh_below_hscp_iso))
                                            else:
                                                Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(-1)
                                            if len(best_cdt_hscp_iso) != 0:
                                                CaloVector.remove(best_cdt_hscp_iso[0])
    
                                        elif len(id_list_hscp_iso) == 3:
                                            worst_cdt_3_hscp_iso = FurthestInRatio(id_list_hscp_iso,new_ratio_energy_hscp_iso,1,idx_phi,idx_eta_hscp_iso,Xi2_choice_3_hscp_iso)
                                            best_cdt_3_hscp_iso = FurthestInRatio(id_list_hscp_iso,new_ratio_energy_hscp_iso,2,idx_phi_hscp_iso,idx_eta_hscp_iso,Xi2_choice_3_hscp_iso)
                                            nb_cdt_out_2by2_hscp_iso[2] +=1
                                            if len(worst_cdt_3_hscp_iso) != 0:
                                                for i in range(len(worst_cdt_3_hscp_iso)):
                                                    sum_all_ngh_below_hscp_iso += worst_cdt_3_hscp_iso[i][5] + worst_cdt_3_hscp_iso[i][6]
                                                    nb_ngh_below_hscp_iso +=1
                                                if nb_ngh_below_hscp_iso != 8:
                                                    Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(sum_all_ngh_below_hscp_iso/(8-nb_ngh_below_hscp_iso))
                                                else:
                                                    Sum_Energy_8_neighbours_threshold_HSCP_ISO.Fill(-1)                                                
                                            if len(best_cdt_3_hscp_iso) != 0:
                                                for i in range(len(best_cdt_3_hscp_iso)):
                                                    CaloVector.remove((best_cdt_3_hscp_iso[i]))
    
    
      
                                else:
                                    HSCP_ISO_Number_non_physical_event.Fill(1)
                                    
            nb_matched_hscp_vs_nb_hscp.Fill(nb_hscp_evt,hscparticle.size())
            CaloVector.clear()


        if CaloVectorAll:
            nb_non_phys = 0
            for l in range(len(CaloVectorAll)-1,0,-1):
                #print(len(CaloVectorAll))
                if CaloVectorAll[l][6] > 0:
                    if CaloVectorAll[l][5] > 0.6 and CaloVectorAll[l][6] > 2.5 and (CaloVectorAll[l][5] / CaloVectorAll[l][6]) < 0.3 and (CaloVectorAll[l][5] + CaloVectorAll[l][6]) < 15:
                        nb_tower_evt+=1
                        totnb_tower_evt+=1
                        ratio_ecal_hcal_seeds.Fill(CaloVectorAll[l][5] / CaloVectorAll[l][6])
                        nb_denom_eff +=1
                        idx_phi,idx_eta,emEnergy,hadenergy,phi,eta = CaloVectorAll[l][1], CaloVectorAll[l][2],CaloVectorAll[l][5],CaloVectorAll[l][6],CaloVectorAll[l][3],CaloVectorAll[l][4]

                        new_idx_phi,new_idx_eta,new_ratio_energy,new_phi,new_eta = FindTrueSeed(CaloVectorAll,idx_phi,idx_eta,(emEnergy/hadenergy),"both",phi,eta)

                        if new_idx_phi != 999999:
                            all_nb_shift=0
                            while ((new_idx_phi != idx_phi) or (new_idx_eta != idx_eta)):
                                idx_phi,idx_eta,phi,eta = new_idx_phi,new_idx_eta,new_phi,new_eta
                                new_idx_phi,new_idx_eta,new_ratio_energy,new_phi,new_eta = FindTrueSeed(CaloVectorAll,idx_phi,idx_eta,new_ratio_energy,"both",new_phi,new_eta)
                                all_nb_shift +=1
    
    
                        if new_idx_phi != 999999:
                            Number_shift_seedOI.Fill(all_nb_shift)
    
                            id_list,nb_all_ngh,nb_ngh_below,nb_ngh_above,sum_all_ngh_below = Create3by3MatrixLoop(CaloVectorAll,idx_phi,idx_eta,'both',2)
                            sq_oi = 5
                            nb_towers_8_above_trh.Fill(nb_ngh_above)
                            Number_all_2by2_max_event.Fill(len(id_list)+1)
                            min_dr_hscp = 999999
                            for k in range(len(HSCPVector)-1,0,-1):
                                dr_hscp = deltaR(deltaR2(new_eta,new_phi,HSCPVector[k][2],HSCPVector[k][1]))
                                if dr_hscp < min_dr_hscp:
                                    min_dr_hscp = dr_hscp
                                    id_hscp = HSCPVector[k][0]
                                    hscp_cdt = HSCPVector[k]

                            if min_dr_hscp < 0.2:
                                nb_num_eff +=1
                                HSCPVector.remove(hscp_cdt)

                            if len(id_list) == 0:
                                nb_cdt_2by2[0] += 1
                                nb_towers_8_below_trh.Fill(nb_ngh_below)
                                nb_towers_8_all_trh.Fill(nb_ngh_below + nb_ngh_above)
                                if nb_ngh_below != 8:
                                    Sum_Energy_8_neighbours_threshold.Fill(sum_all_ngh_below/(8-nb_ngh_below))
                                else:
                                    Sum_Energy_8_neighbours_threshold.Fill(-1)
                                

                            elif len(id_list) == 1 or len(id_list) == 2 or len(id_list) == 3:
                                sq_oi = IsInMatrix(id_list,idx_phi,idx_eta)
                                if sq_oi == 1:
                                        
                                    nb_cdt_2by2[len(id_list)] += 1
                                    nb_towers_8_below_trh.Fill(nb_ngh_below)
                                    nb_towers_8_all_trh.Fill(nb_ngh_below + nb_ngh_above)
                                    if nb_ngh_below != 8:
                                        Sum_Energy_8_neighbours_threshold.Fill(sum_all_ngh_below/(8-nb_ngh_below))
                                    else:
                                        Sum_Energy_8_neighbours_threshold.Fill(-1)

                                    for o in range(len(id_list)):
                                        CaloVectorAll.remove((id_list[o]))
                             
                                else:
                                    if len(id_list) == 1:
                                        sum_all_ngh_below += (id_list[0][5] + id_list[0][6])
                                        nb_ngh_below +=1
                                        if nb_ngh_below != 8:
                                            Sum_Energy_8_neighbours_threshold.Fill(sum_all_ngh_below/(8-nb_ngh_below))
                                        else:
                                            Sum_Energy_8_neighbours_threshold.Fill(-1)

                                    if len(id_list) == 2:
                                        worst_cdt = FurthestInRatio(id_list,new_ratio_energy,1,idx_phi,idx_eta,Xi2_choice_2_and_1)
                                        best_cdt = FurthestInRatio(id_list,new_ratio_energy,2,idx_phi,idx_eta,Xi2_choice_2_and_1)
                                        nb_cdt_out_2by2[1] +=1
                                        for o in range(len(worst_cdt)):
                                            sum_all_ngh_below += (worst_cdt[o][5] + worst_cdt[o][6])
                                            nb_ngh_below +=1

                                        if nb_ngh_below != 8:
                                            Sum_Energy_8_neighbours_threshold.Fill(sum_all_ngh_below/(8-nb_ngh_below))
                                        else:
                                            Sum_Energy_8_neighbours_threshold.Fill(-1)

                                        if len(best_cdt) != 0:
                                            CaloVectorAll.remove(best_cdt[0])

                                    elif len(id_list) == 3:
                                        worst_cdt_3 = FurthestInRatio(id_list,new_ratio_energy,1,idx_phi,idx_eta,Xi2_choice_3)
                                        best_cdt_3 = FurthestInRatio(id_list,new_ratio_energy,2,idx_phi,idx_eta,Xi2_choice_3)
                                        nb_cdt_out_2by2[2] +=1
                                        if len(worst_cdt_3) != 0:
                                            for i in range(len(worst_cdt_3)):
                                                sum_all_ngh_below += worst_cdt_3[i][5] + worst_cdt_3[i][6]
                                                nb_ngh_below +=1
                                            if nb_ngh_below != 8:
                                                Sum_Energy_8_neighbours_threshold.Fill(sum_all_ngh_below/(8-nb_ngh_below))
                                            else:
                                                Sum_Energy_8_neighbours_threshold.Fill(-1)                                                
                                        if len(best_cdt_3) != 0:
                                            for i in range(len(best_cdt_3)):
                                                CaloVectorAll.remove((best_cdt_3[i]))

                            else:
                                nb_non_phys +=1
                                Number_non_physical_event.Fill(1)

            if nb_denom_eff !=0:
                Eff_Nseeds_over_nseedshscp.Fill((nb_num_eff/nb_denom_eff)*100)

            NB_non_phys_vs_nb_tower_passing.Fill(nb_non_phys,nb_tower_evt)


            CaloVectorAll.clear()
            NB_HSCP_vs_NBTOWER.Fill(hscparticle.size(),nb_tower_evt)
            nb_towers_cut_per_Event.Fill(nb_tower_evt)
        
    hfile.Write()
    print("There were in total : ", totnb_tower_evt , " seeds passing thresholds")

    print("ON ALL CALO TOWERS \n")
    print("--------------- DIRECT GOOD EVENTS ---------------")

    print(nb_cdt_2by2[0], " were solo seed, ", nb_cdt_2by2[1], " were seeds with 1 neighbours not in the diagonal", nb_cdt_2by2[2], " were seeds with 2 adjacent neighbours, and ",nb_cdt_2by2[3], " were seeds with 3 neighbours in an L shape (4 possibility since 4 corners in 3 by 3)")

    print("--------------- BAD EVENTS WHERE WE HAVE TO CHOOSE  ---------------")
    print(nb_cdt_out_2by2[1], " were seeds with 2 non-adjacent neighbours, pick the one with ratio ecal/hcal closer to initial seed, ", nb_cdt_out_2by2[2], " were seeds with 3 neighbours, where we have to chose (either 3 solo seeds, or 3 in a line or 2 adjacent and 1 not diagonal)")

    print("ON HSCP MATCHED + ISO CALO TOWERS \n")
    print("--------------- DIRECT GOOD EVENTS ---------------")

    print(nb_cdt_2by2_hscp_iso[0], " were solo seed, ", nb_cdt_2by2_hscp_iso[1], " were seeds with 1 neighbours not in the diagonal", nb_cdt_2by2_hscp_iso[2], " were seeds with 2 adjacent neighbours, and ",nb_cdt_2by2_hscp_iso[3], " were seeds with 3 neighbours in an L shape (4 possibility since 4 corners in 3 by 3)")

    print("--------------- BAD EVENTS WHERE WE HAVE TO CHOOSE  ---------------")
    print(nb_cdt_out_2by2_hscp_iso[1], " were seeds with 2 non-adjacent neighbours, pick the one with ratio ecal/hcal closer to initial seed, ", nb_cdt_out_2by2_hscp_iso[2], " were seeds with 3 neighbours, where we have to chose (either 3 solo seeds, or 3 in a line, or 2 adjacent and 1 not diagonal)")

    print("There was ",nb_hscp_pre_iso, " pair hscp-tower pre iso, and ", nb_hscp_after_iso, " pair hscp-tower aftert iso")
