import math

def SetBin(hist,binx,biny,value):
    hist.SetBinContent(binx,biny,value)

def FillHist(hist,value):
    hist.Fill(value)


def ActivityTab(radius_idx_map,neighboor_cells):
    nb_wrong_neighboor = 0
    for i in range(radius_idx_map):
        for j in range(radius_idx_map):
                if neighboor_cells[i][j] == False or neighboor_cells[i][j] == 0:
                    nb_wrong_neighboor+=1

    return nb_wrong_neighboor

def SetActivityAround(radius_idx_map,neighboor_cells,hist):
    nb_wrong_neighboor = 0
    for i in range(radius_idx_map):
        for j in range(radius_idx_map):
            if neighboor_cells[i][j] == False:
                hist.SetBinContent(i+1,j+1,-1)
                nb_wrong_neighboor+=1
    return nb_wrong_neighboor

def PassThreshold(subsystem,threshold,emEnergy,hadEnergy,sum_cut):
    ecalmip,hcalmip = 0.27*threshold,1.25*threshold
    sum_mip = 1.52*sum_cut

    if (emEnergy > ecalmip) and (hadEnergy > hcalmip) and ((emEnergy/hadEnergy) < 0.4) and ((emEnergy/hadEnergy) > 0) and ((emEnergy+hadEnergy) < sum_mip):
        return True
    else:
        return False


def Create3by3MatrixLoop(CaloVector,id_central_phi,id_central_eta,calo,nbmipngh,sum_cut):
    nb_ngh_all,nb_ngh_below,nb_ngh_above,sum_emhad = 0,0,0,0
    idx_towers_above = [] #we go tower #, idx phi and idx eta
    if CaloVector:
        for i in range(len(CaloVector)):
            a = PassThreshold('ecal',nbmipngh,CaloVector[i][5],CaloVector[i][6],sum_cut)

            if calo == 'both':
                if CaloVector[i][1] == id_central_phi + 1:
                    if CaloVector[i][2] == id_central_eta - 1:
                        nb_ngh_all +=1
                        if a == False:
                            sum_emhad += (CaloVector[i][5]+CaloVector[i][6])
                            nb_ngh_below += 1
                        else:
                            nb_ngh_above +=1
                            idx_towers_above.append(CaloVector[i])


                    elif CaloVector[i][2] == id_central_eta:
                        nb_ngh_all +=1
                        if a == False:
                            sum_emhad += (CaloVector[i][5]+CaloVector[i][6])
                            nb_ngh_below += 1
                        else:
                            nb_ngh_above +=1
                            idx_towers_above.append(CaloVector[i])

        
                    elif CaloVector[i][2] == id_central_eta + 1:
                        nb_ngh_all +=1
                        if a == False:
                            sum_emhad += (CaloVector[i][5]+CaloVector[i][6])
                            nb_ngh_below += 1
                        else:
                            nb_ngh_above +=1
                            idx_towers_above.append(CaloVector[i])



                elif CaloVector[i][1] == id_central_phi:
                    if CaloVector[i][2] == id_central_eta - 1:
                        nb_ngh_all +=1
                        if a == False:
                            sum_emhad += (CaloVector[i][5]+CaloVector[i][6])
                            nb_ngh_below += 1
                        else:
                            idx_towers_above.append(CaloVector[i])
                            nb_ngh_above +=1

                    elif CaloVector[i][2] == id_central_eta + 1:
                        nb_ngh_all +=1
                        if a == False:
                            sum_emhad += (CaloVector[i][5]+CaloVector[i][6])
                            nb_ngh_below += 1
                        else:
                            idx_towers_above.append(CaloVector[i])
                            nb_ngh_above +=1


                elif CaloVector[i][1] == id_central_phi - 1:
                    if CaloVector[i][2] == id_central_eta - 1:
                        nb_ngh_all +=1
                        if a == False:
                            sum_emhad += (CaloVector[i][5]+CaloVector[i][6])
                            nb_ngh_below += 1
                        else:
                            idx_towers_above.append(CaloVector[i])
                            nb_ngh_above +=1


                    elif CaloVector[i][2] == id_central_eta:
                        nb_ngh_all +=1
                        if a == False:
                            sum_emhad += (CaloVector[i][5]+CaloVector[i][6])
                            nb_ngh_below += 1
                        else:
                            idx_towers_above.append(CaloVector[i])
                            nb_ngh_above +=1


        
                    elif CaloVector[i][2] == id_central_eta + 1:
                        nb_ngh_all +=1
                        if a == False:
                            sum_emhad += (CaloVector[i][5]+CaloVector[i][6])
                            nb_ngh_below += 1
                        else:
                            idx_towers_above.append(CaloVector[i])
                            nb_ngh_above +=1


    return (idx_towers_above,nb_ngh_all,nb_ngh_below,nb_ngh_above,sum_emhad)


def FurthestInRatio(id_list,ratio_en,mode,id_phi,id_eta,hist):
    list_duo = []

    if len(id_list) == 2:
        ratio_1 = id_list[0][5]/id_list[0][6]
        ratio_2 = id_list[1][5]/id_list[1][6]
        dist_seed_1 = abs((id_list[0][1] - id_phi) - (id_list[0][2] - id_eta))
        dist_seed_2 = abs((id_list[1][1] - id_phi) - (id_list[1][2] - id_eta))
        if dist_seed_1 != 1 and dist_seed_2 == 1:
            if mode == 1:
                list_duo.append(id_list[0])
            else:
                list_duo.append(id_list[1])
        elif dist_seed_2 != 1 and dist_seed_1 == 1:
            if mode == 1:
                list_duo.append(id_list[1])
            else:
                list_duo.append(id_list[0])

        elif dist_seed_1 == 1 and dist_seed_2 == 1:
            if abs(ratio_en - ratio_1) < abs(ratio_en - ratio_2):
                if mode == 1:
                    list_duo.append(id_list[1])
                else:
                    list_duo.append(id_list[0])

        elif dist_seed_1 != 1 and dist_seed_2 != 1:
            if mode == 1:
                list_duo.append(id_list[0])
                list_duo.append(id_list[1])



    elif len(id_list) == 3:
        ratio_1 = id_list[0][5]/id_list[0][6]
        ratio_2 = id_list[1][5]/id_list[1][6]
        ratio_3 = id_list[2][5]/id_list[2][6]
        dist_seed_1 = abs((id_list[0][1] - id_phi) - (id_list[0][2] - id_eta))
        dist_seed_2 = abs((id_list[1][1] - id_phi) - (id_list[1][2] - id_eta))
        dist_seed_3 = abs((id_list[2][1] - id_phi) - (id_list[2][2] - id_eta))
        dist12 = abs((id_list[0][1] - id_list[1][1]) - (id_list[0][2] - id_list[1][2]))
        dist13 = abs((id_list[0][1] - id_list[2][1]) - (id_list[0][2] - id_list[2][2]))
        dist23 = abs((id_list[1][1] - id_list[2][1]) - (id_list[1][2] - id_list[2][2]))

        if dist12 == 1 and dist23 == 1 and dist13 == 2 and (id_list[2][1] - id_list[0][1]) == 0:
            xi12 = ((ratio_1 - 0.2)**2 + (ratio_2 - 0.2)**2 )/ 2
            xi23 = ((ratio_2 - 0.2)**2 + (ratio_3 - 0.2)**2 )/ 2
            if xi12 < xi23:
                hist.Fill(xi12)
                if mode == 1:
                    list_duo.append(id_list[2])
                else:
                    list_duo.append(id_list[0])
                    list_duo.append(id_list[1])
            else:
                hist.Fill(xi23)
                if mode == 1:
                    list_duo.append(id_list[0])
                else:
                    list_duo.append(id_list[1])
                    list_duo.append(id_list[2])

        
        elif dist12 == 1 and dist13 == 1 and dist23 == 2 and (id_list[1][1] - id_list[2][1]) == 0:
            xi12 = ((ratio_1 - 0.2)**2 + (ratio_2 - 0.2)**2 )/ 2
            xi13 = ((ratio_1 - 0.2)**2 + (ratio_3 - 0.2)**2 )/ 2
            if xi12 < xi13:
                hist.Fill(xi12)
                if mode == 1:
                    list_duo.append(id_list[2])
                else:
                    list_duo.append(id_list[0])
                    list_duo.append(id_list[1])
            else:
                hist.Fill(xi13)
                if mode == 1:
                    list_duo.append(id_list[1])
                else:
                    list_duo.append(id_list[0])
                    list_duo.append(id_list[2])


        elif dist23 == 1 and dist13 == 1 and dist12 == 2 and (id_list[1][1] - id_list[0][1]) == 0:
            xi23 = ((ratio_2 - 0.2)**2 + (ratio_3 - 0.2)**2 )/ 2
            xi13 = ((ratio_1 - 0.2)**2 + (ratio_3 - 0.2)**2 )/ 2
            if xi23 < xi13:
                hist.Fill(xi23)
                if mode == 1:
                    list_duo.append(id_list[0])
                else:
                    list_duo.append(id_list[1])
                    list_duo.append(id_list[2])
            else:
                hist.Fill(xi13)
                if mode == 1:
                    list_duo.append(id_list[1])
                else:
                    list_duo.append(id_list[0])
                    list_duo.append(id_list[2])

        elif dist12 == 1 and dist13 != 1 and dist23 !=1:
            xi12 = ((ratio_1 - 0.2)**2 + (ratio_2 - 0.2)**2 )/ 2
            if dist_seed_3 == 1:
                if (ratio_3 - 0.2)**2 < xi12:
                    hist.Fill((ratio_3 - 0.2)**2)
                    if mode == 1:
                        list_duo.append(id_list[0])
                        list_duo.append(id_list[1])
                    else:
                        list_duo.append(id_list[2])
                else:
                    if mode == 1:
                        list_duo.append(id_list[2])
                    else:
                        list_duo.append(id_list[0])
                        list_duo.append(id_list[1])
                

            else:
                if mode == 1:
                    list_duo.append(id_list[2])
                else:
                    list_duo.append(id_list[0])
                    list_duo.append(id_list[1])

        elif dist13 == 1 and dist12 != 1 and dist23 !=1:
            xi13 = ((ratio_1 - 0.2)**2 + (ratio_3 - 0.2)**2 )/ 2
            if dist_seed_2 == 1:
                if (ratio_2 - 0.2)**2 < xi13:
                    hist.Fill((ratio_2 - 0.2)**2)
                    if mode == 1:
                        list_duo.append(id_list[0])
                        list_duo.append(id_list[2])
                    else:
                        list_duo.append(id_list[1])
                else:
                    if mode == 1:
                        list_duo.append(id_list[1])
                    else:
                        list_duo.append(id_list[0])
                        list_duo.append(id_list[2])
                

            else:
                if mode == 1:
                    list_duo.append(id_list[1])
                else:
                    list_duo.append(id_list[0])
                    list_duo.append(id_list[2])


        elif dist23 == 1 and dist13 != 1 and dist12 !=1:
            xi23 = ((ratio_2 - 0.2)**2 + (ratio_3 - 0.2)**2 )/ 2
            if dist_seed_1 == 1:
                if (ratio_1 - 0.2)**2 < xi23:
                    hist.Fill((ratio_1 - 0.2)**2)
                    if mode == 1:
                        list_duo.append(id_list[1])
                        list_duo.append(id_list[2])
                    else:
                        list_duo.append(id_list[0])
                else:
                    if mode == 1:
                        list_duo.append(id_list[0])
                    else:
                        list_duo.append(id_list[1])
                        list_duo.append(id_list[2])
            else:
                if mode == 1:
                    list_duo.append(id_list[0])
                else:
                    list_duo.append(id_list[1])
                    list_duo.append(id_list[2])



        elif dist12 != 1 and dist13 != 1 and dist23 !=1:
            if dist_seed_1 == 1 and dist_seed_2 != 1 and dist_seed_3 != 1:
                if mode == 1:
                    list_duo.append(id_list[1])
                    list_duo.append(id_list[2])
                else:
                    list_duo.append(id_list[0])
            elif dist_seed_2 == 1 and dist_seed_1 !=1 and dist_seed_3 != 1:
                if mode == 1:
                    list_duo.append(id_list[0])
                    list_duo.append(id_list[2])
                else:
                    list_duo.append(id_list[1])
         
            elif dist_seed_3 == 1 and dist_seed_1 !=1 and dist_seed_2 != 1:
                if mode == 1:
                    list_duo.append(id_list[0])
                    list_duo.append(id_list[1])
                else:
                    list_duo.append(id_list[2])
        
    return list_duo

def CountIsoNum(isov,isotab):
    if isov < 0.1:
        isotab[0] +=1
    if isov < 0.3:
        isotab[1] +=1
    if isov < 0.5:
        isotab[2] +=1
    if isov < 0.7:
        isotab[3] +=1
    if isov < 0.9:
        isotab[4] +=1


def FillStepByStep(tab,emEnergy,hadEnergy):
    tab[0] += 1
    if emEnergy > 0.27:
        tab[1] +=1
        if emEnergy > 0.54:
            tab[2] +=1
            if emEnergy > 0.81:
                tab[3] +=1
                if hadEnergy > 0:
                    tab[4] +=1
                    if hadEnergy > 1.25:
                        tab[5] +=1

def FillPdgIds(cpt_list,tab):
    for i in range(len(cpt_list)):
        tab.append(abs(cpt_list[i].pdgId()))

def FindChargedHSCP(all_ids,list_ch_id):
    nb_ch = 0
    for i in range(len(list_ch_id)):
        nb_ch += all_ids.count(list_ch_id[i])

    return nb_ch

def FindAllHSCP(all_ids,list_all_id):
    nb_all = 0
    for i in range(len(list_all_id)):
        nb_all += all_ids.count(list_all_id[i])

    return nb_all

def IsInMatrix(id_list,id_phi,id_eta):
    if len(id_list) == 1:
        dist_seed_1 = abs((id_list[0][1] - id_phi) - (id_list[0][2] - id_eta))
        if dist_seed_1 == 1:
            return 1
        else:
            return 5

    elif len(id_list) == 2:
        dist_seed_1 = abs((id_list[0][1] - id_phi) - (id_list[0][2] - id_eta))
        dist_seed_2 = abs((id_list[1][1] - id_phi) - (id_list[1][2] - id_eta))
        dist12 = abs((id_list[0][1] - id_list[1][1]) - (id_list[0][2] - id_list[1][2]))
        ratio_1 = id_list[0][5]/id_list[0][6]
        ratio_2 = id_list[1][5]/id_list[1][6]
        
        if dist12 == 1:
            return 1
        else:
            return 5

    elif len(id_list) == 3:
        dist_seed_1 = abs((id_list[0][1] - id_phi) - (id_list[0][2] - id_eta))
        dist_seed_2 = abs((id_list[1][1] - id_phi) - (id_list[1][2] - id_eta))
        dist_seed_3 = abs((id_list[2][1] - id_phi) - (id_list[2][2] - id_eta))
        dist12 = abs((id_list[0][1] - id_list[1][1]) - (id_list[0][2] - id_list[1][2]))
        dist13 = abs((id_list[0][1] - id_list[2][1]) - (id_list[0][2] - id_list[2][2]))
        dist23 = abs((id_list[1][1] - id_list[2][1]) - (id_list[1][2] - id_list[2][2]))
        if dist12 == 1 and dist23 == 1 and dist13 == 2 and abs(id_list[0][2] - id_list[2][2]) == 1 and abs(id_list[0][1] - id_list[2][1]) == 1:
            return 1

        elif dist12 == 1 and dist13==1 and dist23 == 2 and abs(id_list[1][2] - id_list[2][2]) == 1 and abs(id_list[1][1] - id_list[2][1]) == 1:
            return 1

        elif dist23 == 1 and dist13 == 1 and dist12 == 2 and abs(id_list[0][2] - id_list[1][2]) == 1 and abs(id_list[0][1] - id_list[1][1]) == 1:
            return 1
        else:
            return 5
    else:
        return 5
        
        


def FindAdjacentPair(phi1,eta1,phi2,eta2,phi_central,eta_central):
    if abs((phi1 - phi2) - (eta1 - eta2)) == 1:
        if phi1 > phi2:
            if (phi1 - phi_central) == 1 and (eta1 - eta_central) == 1:
                 return 3
            elif (phi1 - phi_central) == 0 and (eta1 - eta_central) == 1:
                return 4
            elif (phi1 - phi_central) == 1 and (eta1 - eta_central) == -1:
                return 2
            elif (phi1 - phi_central) == 0 and (eta1 - eta_central) == -1:
                return 1

        elif phi1 < phi2:
            if (phi2 - phi_central) == 1 and (eta2 - eta_central) == 1:
                return 3
            elif (phi2 - phi_central) == 0 and (eta2 - eta_central) == 1:
                return 4
            elif (phi2 - phi_central) == 1 and (eta2 - eta_central) == -1:
                return 2
            elif (phi2 - phi_central) == 0 and (eta2 - eta_central) == -1:
                return 1

        elif phi1 == phi2:                                            
            if eta1 > eta2:
                if (eta1 - eta_central ) == 0 and (phi1 - phi_central) == -1:
                    return 1
                elif (eta1 - eta_central ) == 0 and (phi1 - phi_central) == 1:
                    return 2
                elif (eta1 - eta_central ) == 1 and (phi1 - phi_central) == 1:
                    return 3
                elif (eta1 - eta_central ) == 1 and (phi1 - phi_central) == -1:
                    return 4

            elif eta1 < eta2:
                if (eta2 - eta_central ) == 0 and (phi2 - phi_central) == -1:
                    return 1
                elif (eta2 - eta_central ) == 0 and (phi2 - phi_central) == 1:
                    return 2
                elif (eta2 - eta_central ) == 1 and (phi2 - phi_central) == 1:
                    return 3
                elif (eta2 - eta_central ) == 1 and (phi2 - phi_central) == -1:
                    return 4

    elif abs(phi1-phi2) == 1 and abs(eta1-eta2)==1:

        if phi1 > phi2:
            if (phi1 - phi_central) == 0 and (eta1 - eta_central) == 1 and (phi2-phi_central) == -1 and (eta2-eta_central) == 0:
                return 4
            elif (phi1 - phi_central) == 0 and (eta1 - eta_central) == -1 and (phi2 - phi_central) == -1 and (eta2-eta_central) == 0:
                return 1
            elif (phi1 - phi_central) == 1 and (eta1 - eta_central) == 0 and (phi2 - phi_central) == 0 and (eta2-eta_central) == -1:
                return 2
            elif (phi1 - phi_central) == 1 and (eta1 - eta_central) == 0 and (phi2-phi_central) == 0 and (eta2-eta_central) == 1:
                return 3
        elif phi2 > phi1:
            if (phi1 - phi_central) == -1 and (eta1 - eta_central) == 0 and (phi2-phi_central) == 0 and (eta2-eta_central) == 1:
                return 4
            elif (phi1 - phi_central) == -1 and (eta1 - eta_central) == 0 and (phi2 - phi_central) == 0 and (eta2-eta_central) == -1:
                return 1
            elif (phi1 - phi_central) == 0 and (eta1 - eta_central) == 1 and (phi2 - phi_central) == 1 and (eta2-eta_central) == 0:
                return 3
            elif (phi1 - phi_central) == 0 and (eta1 - eta_central) == -1 and (phi2-phi_central) == 1 and (eta2-eta_central) == 0:
                return 2
        

    else:
        return 999999

def FindTrueSeed(CaloVector,id_central_phi,id_central_eta,ratio_energy,calo,phi_central,eta_central,indice,sum_cut,nbmipngh):
    #print("Find true seed called, id_phi : ", id_central_phi, " , id_eta :", id_central_eta, " , sum energy : ", sum_energy)
    new_highest = []
    new_highest.append(((abs(0.2 - ratio_energy)),id_central_phi,id_central_eta,ratio_energy,phi_central,eta_central,indice))
    if CaloVector:
        for i in range(len(CaloVector)):
            a = PassThreshold('ecal',nbmipngh,CaloVector[i][5],CaloVector[i][6],sum_cut)
            if calo == 'both':
                if CaloVector[i][1] == id_central_phi + 1:
                    if CaloVector[i][2] == id_central_eta - 1:
                        if a:
                            new_highest.append((abs(0.2 - (CaloVector[i][5]/CaloVector[i][6])),CaloVector[i][1],CaloVector[i][2],CaloVector[i][5]/CaloVector[i][6],CaloVector[i][3],CaloVector[i][4],i))


                    elif CaloVector[i][2] == id_central_eta:
                        if a:
                            new_highest.append((abs(0.2 - (CaloVector[i][5]/CaloVector[i][6])),CaloVector[i][1],CaloVector[i][2],CaloVector[i][5]/CaloVector[i][6],CaloVector[i][3],CaloVector[i][4],i))

                    elif CaloVector[i][2] == id_central_eta + 1:
                        if a:
                            new_highest.append((abs(0.2 - (CaloVector[i][5]/CaloVector[i][6])),CaloVector[i][1],CaloVector[i][2],CaloVector[i][5]/CaloVector[i][6],CaloVector[i][3],CaloVector[i][4],i))

                elif CaloVector[i][1] == id_central_phi:
                    if CaloVector[i][2] == id_central_eta - 1:
                        if a:
                            new_highest.append((abs(0.2 - (CaloVector[i][5]/CaloVector[i][6])),CaloVector[i][1],CaloVector[i][2],CaloVector[i][5]/CaloVector[i][6],CaloVector[i][3],CaloVector[i][4],i))

                    elif CaloVector[i][2] == id_central_eta + 1:
                        if a:
                            new_highest.append((abs(0.2 - (CaloVector[i][5]/CaloVector[i][6])),CaloVector[i][1],CaloVector[i][2],CaloVector[i][5]/CaloVector[i][6],CaloVector[i][3],CaloVector[i][4],i))

                elif CaloVector[i][1] == id_central_phi - 1:
                    if CaloVector[i][2] == id_central_eta - 1:
                        if a:
                            new_highest.append((abs(0.2 - (CaloVector[i][5]/CaloVector[i][6])),CaloVector[i][1],CaloVector[i][2],CaloVector[i][5]/CaloVector[i][6],CaloVector[i][3],CaloVector[i][4],i))

                    elif CaloVector[i][2] == id_central_eta:
                        if a:
                            new_highest.append((abs(0.2 - (CaloVector[i][5]/CaloVector[i][6])),CaloVector[i][1],CaloVector[i][2],CaloVector[i][5]/CaloVector[i][6],CaloVector[i][3],CaloVector[i][4],i))

                    elif CaloVector[i][2] == id_central_eta + 1:
                        if a:
                            new_highest.append((abs(0.2 - (CaloVector[i][5]/CaloVector[i][6])),CaloVector[i][1],CaloVector[i][2],CaloVector[i][5]/CaloVector[i][6],CaloVector[i][3],CaloVector[i][4],i))




    min_val = min(new_highest)
    min_idx = new_highest.index(min_val)

    if min_idx == 0:
        return (new_highest[min_idx][1],new_highest[0][2],new_highest[0][3],new_highest[min_idx][4],new_highest[min_idx][5],new_highest[min_idx][6])

    else: 
        new_id_phi,new_id_eta,new_ratio_energy,new_phi,new_eta,new_indice = min(new_highest)[1],min(new_highest)[2],min(new_highest)[3],min(new_highest)[4],min(new_highest)[5],min(new_highest)[6]
        #print("Returning new indexes -> phi :",new_id_phi, " , eta : ", new_id_eta, " energy : ", new_sum_energy)
        return (new_id_phi,new_id_eta,new_ratio_energy,new_phi,new_eta,new_indice)




def CreateCaloTabLoop(CaloTowers,id_central_phi,id_central_eta,neighboors):
    if CaloTowers is not None:
        for i in range(CaloTowers.size()):
            if CaloTowers[i].iphi() == id_central_phi+2:
                if CaloTowers[i].ieta() == id_central_eta-1:
                    neighboors[1][4] = True
                elif CaloTowers[i].ieta() == id_central_eta+1:
                    neighboors[3][4] = True
                elif CaloTowers[i].ieta() == id_central_eta:
                    neighboors[2][4] = True
                elif CaloTowers[i].ieta() == id_central_eta-2:
                    neighboors[0][4] = True
                elif CaloTowers[i].ieta() == id_central_eta+2:
                    neighboors[4][4] = True
        
            elif CaloTowers[i].iphi() == id_central_phi+1:
                if CaloTowers[i].ieta() == id_central_eta-1:
                    neighboors[1][3] = True
                elif CaloTowers[i].ieta() == id_central_eta+1:
                    neighboors[3][3] = True
                elif CaloTowers[i].ieta() == id_central_eta:
                    neighboors[3][2] = True
                elif CaloTowers[i].ieta() == id_central_eta-2:
                    neighboors[0][3] = True
                elif CaloTowers[i].ieta() == id_central_eta+2:
                    neighboors[4][3] = True
        
            elif CaloTowers[i].iphi() == id_central_phi:
                if CaloTowers[i].ieta() == id_central_eta-1:
                    neighboors[1][2] = True
                elif CaloTowers[i].ieta() == id_central_eta+1:
                    neighboors[3][2] = True
                elif CaloTowers[i].ieta() == id_central_eta-2:
                    neighboors[0][2] = True
                elif CaloTowers[i].ieta() == id_central_eta+2:
                    neighboors[4][2] = True
        
            elif CaloTowers[i].iphi() == id_central_phi-1:
                if CaloTowers[i].ieta() == id_central_eta-1:
                    neighboors[1][1] = True
                elif CaloTowers[i].ieta() == id_central_eta+1:
                    neighboors[3][1] = True
                elif CaloTowers[i].ieta() == id_central_eta:
                    neighboors[2][1] = True
                elif CaloTowers[i].ieta() == id_central_eta-2:
                    neighboors[0][1] = True
                elif CaloTowers[i].ieta() == id_central_eta+2:
                    neighboors[4][1] = True
        
            elif CaloTowers[i].iphi() == id_central_phi-2:
                if CaloTowers[i].ieta() == id_central_eta-1:
                    neighboors[1][0] = True
                elif CaloTowers[i].ieta() == id_central_eta+1:
                    neighboors[3][0] = True
                elif CaloTowers[i].ieta() == id_central_eta:
                    neighboors[2][0] = True
                elif CaloTowers[i].ieta() == id_central_eta-2:
                    neighboors[0][0] = True
                elif CaloTowers[i].ieta() == id_central_eta+2:
                    neighboors[4][0] = True




def CreateCaloMap(hist,id_central_phi,id_central_eta,id_phi,id_eta,emEnergy,hadEnergy,neighboors,calo):
    if id_phi == id_central_phi+2:
        if id_eta == id_central_eta-1:
            neighboors[1][4] = True
            if calo == 'ecal':
                SetBin(hist,2,5,emEnergy)
            elif calo == 'both':
                SetBin(hist,2,5,emEnergy+hadEnergy)

        elif id_eta == id_central_eta+1:
            neighboors[3][4] = True
            if calo == 'ecal':
                SetBin(hist,4,5,emEnergy)
            elif calo == 'both':
                SetBin(hist,4,5,emEnergy,+hadEnergy)

        elif id_eta == id_central_eta:
            neighboors[2][4] = True
            if calo == 'ecal':
                SetBin(hist,3,5,emEnergy)
            elif calo == 'both':
                SetBin(hist,3,5,emEnergy+hadEnergy)
                
        elif id_eta == id_central_eta-2:
            neighboors[0][4] = True
            if calo == 'ecal':
                SetBin(hist,1,5,emEnergy)
            elif calo == 'both':
                SetBin(hist,1,5,emEnergy,+hadEnergy)

        elif id_eta == id_central_eta+2:
            neighboors[4][4] = True
            if calo == 'ecal':
                SetBin(hist,5,5,emEnergy)
            elif calo == 'both':
                SetBin(hist,5,5,emEnergy,+hadEnergy)
        
    elif id_phi == id_central_phi+1:
        if id_eta == id_central_eta-1:
            neighboors[1][3] = True
            if calo == 'ecal':
                SetBin(hist,2,4,emEnergy)
            elif calo == 'both':
                SetBin(hist,2,4,emEnergy+hadEnergy)

        elif id_eta == id_central_eta+1:
            neighboors[3][3] = True
            if calo == 'ecal':
                SetBin(hist,4,4,emEnergy)
            elif calo == 'both':
                SetBin(hist,4,4,emEnergy,+hadEnergy)
   
        elif id_eta == id_central_eta:
            neighboors[3][2] = True
            if calo == 'ecal':
                SetBin(hist,4,3,emEnergy)
            elif calo == 'both':
                SetBin(hist,4,3,emEnergy+hadEnergy)

        elif id_eta == id_central_eta-2:
            neighboors[0][3] = True
            if calo == 'ecal':
                SetBin(hist,1,4,emEnergy)
            elif calo == 'both':
                SetBin(hist,1,4,emEnergy,+hadEnergy)

        elif id_eta == id_central_eta+2:
            neighboors[4][3] = True
            if calo == 'ecal':
                SetBin(hist,5,4,emEnergy)
            elif calo == 'both':
                SetBin(hist,5,4,emEnergy,+hadEnergy)


    elif id_phi == id_central_phi:
        if id_eta == id_central_eta-1:
            neighboors[1][2] = True
            if calo == 'ecal':
                SetBin(hist,2,3,emEnergy)
            elif calo == 'both':
                SetBin(hist,2,3,emEnergy+hadEnergy)


        elif id_eta == id_central_eta+1:
            neighboors[3][2] = True
            if calo == 'ecal':
                SetBin(hist,4,3,emEnergy)
            elif calo == 'both':
                SetBin(hist,4,3,emEnergy,+hadEnergy)

   
        elif id_eta == id_central_eta-2:
            neighboors[0][2] = True
            if calo == 'ecal':
                SetBin(hist,1,3,emEnergy)
            elif calo == 'both':
                SetBin(hist,1,3,emEnergy,+hadEnergy)

        elif id_eta == id_central_eta+2:
            neighboors[4][2] = True
            if calo == 'ecal':
                SetBin(hist,5,3,emEnergy)
            elif calo == 'both':
                SetBin(hist,5,3,emEnergy,+hadEnergy)


    elif id_phi == id_central_phi-1:
        if id_eta == id_central_eta-1:
            neighboors[1][1] = True
            if calo == 'ecal':
                SetBin(hist,2,2,emEnergy)
            elif calo == 'both':
                SetBin(hist,2,2,emEnergy+hadEnergy)


        elif id_eta == id_central_eta+1:
            neighboors[3][1] = True
            if calo == 'ecal':
                SetBin(hist,4,2,emEnergy)
            elif calo == 'both':
                SetBin(hist,4,2,emEnergy,+hadEnergy)

   
        elif id_eta == id_central_eta:
            neighboors[2][1] = True
            if calo == 'ecal':
                SetBin(hist,3,2,emEnergy)
            elif calo == 'both':
                SetBin(hist,3,2,emEnergy+hadEnergy)

        elif id_eta == id_central_eta-2:
            neighboors[0][1] = True
            if calo == 'ecal':
                SetBin(hist,1,2,emEnergy)
            elif calo == 'both':
                SetBin(hist,1,2,emEnergy,+hadEnergy)

        elif id_eta == id_central_eta+2:
            neighboors[4][1] = True
            if calo == 'ecal':
                SetBin(hist,5,2,emEnergy)
            elif calo == 'both':
                SetBin(hist,5,2,emEnergy,+hadEnergy)


    elif id_phi == id_central_phi-2:
        if id_eta == id_central_eta-1:
            neighboors[1][0] = True
            if calo == 'ecal':
                SetBin(hist,2,1,emEnergy)
            elif calo == 'both':
                SetBin(hist,2,1,emEnergy+hadEnergy)


        elif id_eta == id_central_eta+1:
            neighboors[3][0] = True
            if calo == 'ecal':
                SetBin(hist,4,1,emEnergy)
            elif calo == 'both':
                SetBin(hist,4,1,emEnergy,+hadEnergy)
   
        elif id_eta == id_central_eta:
            neighboors[2][0] = True
            if calo == 'ecal':
                SetBin(hist,3,1,emEnergy)
            elif calo == 'both':
                SetBin(hist,3,1,emEnergy+hadEnergy)

        elif id_eta == id_central_eta-2:
            neighboors[0][0] = True
            if calo == 'ecal':
                SetBin(hist,1,1,emEnergy)
            elif calo == 'both':
                SetBin(hist,1,1,emEnergy,+hadEnergy)

        elif id_eta == id_central_eta+2:
            neighboors[4][0] = True
            if calo == 'ecal':
                SetBin(hist,5,1,emEnergy)
            elif calo == 'both':
                SetBin(hist,5,1,emEnergy,+hadEnergy)

def returnIdMax(e1,e2,e3,e4):
    if max(e1,e2,e3,e4) == e4:
        return 4
    elif max(e1,e2,e3,e4) == e3:
        return 3
    elif max(e1,e2,e3,e4) == e2:
        return 2
    elif max(e1,e2,e3,e4) == e1:
        return 1
   

def Get12NeighboursLoop(CaloTowers,idmatrix,id_central_phi,id_central_eta):
    nb_ngh_all,sum_emhad = 0,0
    if CaloTowers is not None:
        for i in range(CaloTowers.size()):
            if idmatrix == 1:
                if CaloTowers[i].iphi() == id_central_phi + 1:
                    if CaloTowers[i].ieta() == id_central_eta - 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
                    elif CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

        
                elif CaloTowers[i].iphi() == id_central_phi:
                    if CaloTowers[i].ieta() == id_central_eta - 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1


        
                elif CaloTowers[i].iphi() == id_central_phi - 1:
                    if CaloTowers[i].ieta() == id_central_eta - 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1


        
                elif CaloTowers[i].iphi() == id_central_phi - 2:
                    if CaloTowers[i].ieta() == id_central_eta - 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
                    elif CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
                    elif CaloTowers[i].ieta() == id_central_eta:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
        
        
            elif idmatrix == 2:
                if CaloTowers[i].iphi() == id_central_phi + 2:
                    if CaloTowers[i].ieta() == id_central_eta - 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

        
                elif CaloTowers[i].iphi() == id_central_phi + 1:
                    if CaloTowers[i].ieta() == id_central_eta - 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
        
                elif CaloTowers[i].iphi() == id_central_phi:
                    if CaloTowers[i].ieta() == id_central_eta - 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1


        
                elif CaloTowers[i].iphi() == id_central_phi - 1:
                    if CaloTowers[i].ieta() == id_central_eta - 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

        
            elif idmatrix == 3:
                if CaloTowers[i].iphi() == id_central_phi + 2:
                    if CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
 
                    elif CaloTowers[i].ieta() == id_central_eta:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
 
                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

        
                elif CaloTowers[i].iphi() == id_central_phi + 1:
                    if CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

        
                elif CaloTowers[i].iphi() == id_central_phi:
                    if CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

        
                elif CaloTowers[i].iphi() == id_central_phi - 1:
                    if CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    if CaloTowers[i].ieta() == id_central_eta + 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

            elif idmatrix == 4:
                if CaloTowers[i].iphi() == id_central_phi + 1:
                    if CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

        
                elif CaloTowers[i].iphi() == id_central_phi:
                    if CaloTowers[i].ieta() == id_central_eta - 1:

                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1
                    elif CaloTowers[i].ieta() == id_central_eta + 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

        
                elif CaloTowers[i].iphi() == id_central_phi - 1:
                    if CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

        
                elif CaloTowers[i].iphi() == id_central_phi - 2:
                    if CaloTowers[i].ieta() == id_central_eta - 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    elif CaloTowers[i].ieta() == id_central_eta + 1:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

                    if CaloTowers[i].ieta() == id_central_eta + 2:
                        sum_emhad += (CaloTowers[i].emEnergy()+CaloTowers[i].hadEnergy())
                        nb_ngh_all += 1

    return (nb_ngh_all,sum_emhad)

#To change, put loop inside the func 
def CheckHighEMatrix(id_central_phi,id_central_eta,id_phi,id_eta,phi,eta,emEnergy,hadEnergy,corner):
    if corner == 1:
        low_eta,low_phi,high_eta,high_phi = 0,0,0,0
        if id_phi == id_central_phi-1:
            if id_eta == id_central_eta:
                low_phi = phi
                high_eta = eta
                return ((hadEnergy+emEnergy),low_phi,high_eta,2)
            elif id_eta == id_central_eta-1:
                return ((hadEnergy+emEnergy),-1,-1,-1)
            else:
                return (0,0,0,0)
        elif id_phi == id_central_phi:
            high_phi = phi
            if id_eta == id_central_eta - 1:
                low_eta = eta
                return ((hadEnergy+emEnergy),high_phi,low_eta,1)
            else:
                return (0,0,0,0)
        else:
            return (0,0,0,0)

    elif corner == 2:
        low_eta,low_phi,high_eta,high_phi = 0,0,0,0
        if id_phi == id_central_phi+1:
            if id_eta == id_central_eta:
                high_phi = phi
                high_eta = eta
                return ((hadEnergy+emEnergy),high_phi,high_eta,3)
            elif id_eta == id_central_eta-1:
                return ((hadEnergy+emEnergy),-1,-1,-1)
            else:
                return (0,0,0,0)
        elif id_phi == id_central_phi:
            if id_eta == id_central_eta - 1:
                low_phi = phi
                low_eta = eta
                return ((hadEnergy+emEnergy),low_phi,low_eta,4)
            else:
                return (0,0,0,0)
        else:
            return (0,0,0,0)

    elif corner == 3:
        low_eta,low_phi,high_eta,high_phi = 0,0,0,0
        if id_phi == id_central_phi+1:
            if id_eta == id_central_eta:
                high_phi = phi
                low_eta = eta
                return ((hadEnergy+emEnergy),high_phi,low_eta,1)
            elif id_eta == id_central_eta+1:
                return ((hadEnergy+emEnergy),-1,-1,-1)
            else:
                return (0,0,0,0)
        elif id_phi == id_central_phi:
            if id_eta == id_central_eta + 1:
                low_phi = phi
                high_eta = eta
                return ((hadEnergy+emEnergy),low_phi,high_eta,2)
            else:
                return (0,0,0,0)
        else:
            return (0,0,0,0)

    elif corner == 4:
        low_eta,low_phi,high_eta,high_phi = 0,0,0,0
        if id_phi == id_central_phi-1:
            if id_eta == id_central_eta:
                low_phi = phi
                low_eta = eta
                return ((hadEnergy+emEnergy),low_phi,low_eta,4)
            elif id_eta == id_central_eta + 1:
                return ((hadEnergy+emEnergy),-1,-1,-1)
            else:
                return (0,0,0,0)
        elif id_phi == id_central_phi:
            if id_eta == id_central_eta + 1:
                high_phi = phi
                high_eta = eta
                return ((hadEnergy+emEnergy),high_phi,high_eta,3)
            else:
                return (0,0,0,0)
        else:
            return (0,0,0,0)

def Preselection(genpart):
    yon = True
    if (genpart.eta() >= 2.1) or (genpart.eta() <= -2.1):
        yon = False
    if (genpart.pt() <= 55):
        yon = False
    if (yon == True):
        return 1
    else:
        return 0

def testDeltaR2(eta1,phi1,eta2,phi2):  #it is done like that in the CommonFunctions.h from SUSYBSMANALYSIS package
    deta = eta1 - eta2
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2 * math.pi
    while dphi <= -math.pi:
        dphi += 2 * math.pi
    return math.sqrt((deta*deta + dphi*dphi))



def deltaR2(eta1,phi1,eta2,phi2):
    diffphi = abs(phi1-phi2)
    if diffphi > math.pi:
        diffphi -= 2 * math.pi
    return (((eta1-eta2)*(eta1-eta2)) + (diffphi*diffphi))

def deltaR(delta):
    return(math.sqrt(delta))


