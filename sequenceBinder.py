#!/usr/bin/env python

import re
import argparse
import numpy as np

from datetime import datetime
startTime = datetime.now()


def atomicCharges(pdbList, pH):
    '''
    Get fraction of charge on ionisable atoms based on input pH,
    used in calcMultiBinding function.
    '''

    residueTypesDict = {
            'aliphatics': ['ALA', 'GLY', 'ILE', 'LEU', 'PRO', 'VAL'],
            'aromatics': ['PHE', 'TRP', 'TYR'],
            'acidics': [['ASP', 'CG', 4], ['GLU', 'CD', 4.4]],
            'basics': [['ARG', 'CZ', 12], ['LYS', 'NZ', 10.4], 
                    ['HIS', 'CE1', 6.3]],
            'positives': [['ARG', 'CZ', 12], ['LYS', 'NZ', 10.4]],
            'hydroxylics': ['SER', 'THR'],
            'sulfites': ['CYS', 'MET'],
            'amidics': ['ASN', 'GLN'],
            'terminals': [['NT', 7.5], ['CT', 3.8]]
    }


    sequence = []

    for atom in pdbList:
        for key, value in residueTypesDict.items():
            if key == 'acidics':
                for v in value:
                    pka = v[2]
                    if (atom[1], atom[0]) == (v[0], v[1]):
                        exp_here = pH - pka
                        ratio_neg = 10 ** exp_here
                        frac_neg = -(float(1) / float(1 + (float(1) / 
                                float(ratio_neg))))
                        sequence.append((atom[0], atom[1], atom[2], 
                                frac_neg, [atom[3], atom[4], atom[5]]))
                    else:
                        continue
            if key == 'basics':
                for v in value:
                    pka = v[2]
                    if (atom[1], atom[0]) == (v[0], v[1]):
                        exp_here = pka - pH
                        ratio_pos = 10 ** exp_here
                        frac_pos = float(1) / float(1 + (float(1) / 
                                float(ratio_pos)))
                        sequence.append((atom[0], atom[1], atom[2], frac_pos, 
                                [atom[3], atom[4], atom[5]]))
                    else:
                        continue
            if key == 'terminals':
                for v in value:
                    if (atom[0]) == (v[0]):
                        pka = v[1]
                        if v[0] == 'NT':
                            exp_here = pka - pH
                            ratio_pos = 10 ** exp_here
                            frac_pos = float(1) / float(1 + (float(1) / 
                                    float(ratio_pos)))
                            sequence.append((atom[0], atom[1], atom[2], 
                                    frac_pos, [atom[3], atom[4], atom[5]]))
                        if v[0] == 'CT':
                            exp_here = pH - pka
                            ratio_neg = 10 ** exp_here
                            frac_neg = -(float(1) / float(1 + (float(1) / 
                                    float(ratio_neg))))
                            sequence.append((atom[0], atom[1], atom[2], 
                                    frac_neg, [atom[3], atom[4], atom[5]]))
                    else:
                        continue


    return sequence





def calcMultiBinding(fileName, pdbList, num_pos_neg_points, *args, **kwargs):
    '''
    Bind counterions to charged multiple amino acids, where all charged 
    residues are clustered together, 
    instead of using patches from electrostatics. So effectively 
    we have 1 big pos patch and 1 big neg patch.
    '''



    pH = kwargs.get('pH', 7)
    conc_neg = kwargs.get('conc_neg', 'nul')
    conc_pos = kwargs.get('conc_pos', 'nul')
    neg_charge = kwargs.get('neg_charge', 1.0)
    pos_charge = kwargs.get('pos_charge', 1.0)


    if conc_pos == 'nul':
        conc_neg = float(conc_neg)
        conc_pos = float(neg_charge * conc_neg) / float(pos_charge)
    if conc_neg == 'nul':
        conc_pos = float(conc_pos)
        conc_neg = float(pos_charge * conc_pos) / float(neg_charge)

    ISCalc = 0.5 * (conc_pos * pos_charge * pos_charge  + 
            conc_neg * neg_charge * neg_charge)



    diel_rel = 78.4 #water dielectric for DH calcs
    diel_zero = 8.85419e-12
    kB = 1.38062e-23
    NA = 6.0219e23
    elec = 1.60219e-19
    T = 300


    es_constA = float(elec) / float(4 * np.pi * diel_zero * diel_rel * 1.0e-10)
    es_constB = elec * es_constA
    #a1 = NA * 1000 * ionstr
    #a2 = diel_zero * diel_rel * kB * T
    #a3 = elec
    #a4 = float(a3) * float((2 * a1) ** 0.5) / float(a2 ** 0.5)
    #debyeHuckel =  a4 * 1.0e-10
    RT_kJ = float(NA * kB * T) / float(1000)


    #for now work with no DH term

    ion_separation = 3.5
    Gbasic = float(es_constB) / float(ion_separation)
    Gbasic =  float(NA * Gbasic) / float(1000)
    Gbasic = -1 * Gbasic


    sequence = atomicCharges(pdbList, pH) ##get charge on protein residues


    pos_res_charges = 0
    neg_res_charges = 0
    tot_charge = 0


    for i in sequence:
        #print (i)
        residue_charge = i[3]
        if residue_charge > 0:
            pos_res_charges += residue_charge
            tot_charge += residue_charge
        elif residue_charge < 0:
            neg_res_charges += residue_charge
            tot_charge += residue_charge
        else:
            continue


    sequence = [pos_res_charges, neg_res_charges]


    start_charge = 0
    pos_start_charge = 0
    neg_start_charge = 0
    pos_bound_charge = 0
    neg_bound_charge = 0
    num_interactions = 0
    frac_pos_bound = 0
    frac_neg_bound = 0

    for i in sequence:
        residue_charge = i
        start_charge += residue_charge
        CI_bound_charge = 0

        ### bind positive counterions to negative group
        if residue_charge < 0:
            neg_start_charge += residue_charge
            reduced_valency = abs(i)
            Gint = 0
            while reduced_valency > 0:
                binding_charge = 0
                if reduced_valency < abs(pos_charge):
                    binding_charge = abs(reduced_valency)
                if reduced_valency > abs(pos_charge):
                    binding_charge = abs(pos_charge)

                Gint = binding_charge * Gbasic
                reduced_valency -= abs(pos_charge)

                Keq = np.exp(float(1.0 * Gint) / float(RT_kJ))
                one_over = 1.0 + float(Keq) / float(conc_pos)
                frac_bound = float(1) /  float(one_over)
                frac_pos_bound += frac_bound
                pos_bound_charge += frac_bound * pos_charge
                CI_bound_charge += frac_bound * pos_charge
                num_interactions += 1

        ### bind negative counterions to positive group
        if residue_charge > 0:
            pos_start_charge += residue_charge
            reduced_valency = abs(i)
            Gint = 0
            while reduced_valency > 0:
                binding_charge = 0
                if reduced_valency < abs(neg_charge):
                    binding_charge = abs(reduced_valency)
                if reduced_valency > abs(neg_charge):
                    binding_charge = abs(neg_charge)

                Gint = binding_charge * Gbasic
                reduced_valency -= abs(neg_charge)

                Keq = np.exp(float(1.0 * Gint) / float(RT_kJ))
                one_over = 1.0 + float(Keq) / float(conc_neg)
                frac_bound = float(1) /  float(one_over)
                frac_neg_bound += frac_bound
                neg_bound_charge += frac_bound * -abs(neg_charge)
                CI_bound_charge += frac_bound * -abs(neg_charge)
                num_interactions += 1



    tot_bound_charge_posneg = start_charge + neg_bound_charge + \
            pos_bound_charge



    ################### data for saving to file
    bindingInfo = [
            ['pH', pH],
            ['negativeIonConc', conc_neg],
            ['positiveIonConc', conc_pos],
            ['negativeIonCharge', neg_charge],
            ['positiveIonCharge', pos_charge],
            ['ionicStrength', ISCalc],
            ['proteinStartCharge', start_charge], 
            ['positiveStartCharge', pos_start_charge],
            ['negativeStartCharge', neg_start_charge],
            ['proteinEndCharge', tot_bound_charge_posneg],
            ['positiveEndCharge', pos_start_charge + neg_bound_charge], 
            ['negativeEndCharge', neg_start_charge + pos_bound_charge],
            ['numberOfBoundNegativeIons', frac_neg_bound],
            ['numberOfBoundPositiveIons', frac_pos_bound],
            ['positiveResiduesCount', len(num_pos_neg_points[0])],
            ['negativeResiduesCount', len(num_pos_neg_points[1])],
            ]



    data = open('%s_bindingInfo_pH%s_pos%s_neg%s_IS_%sM.csv' % 
            (fileName, pH, pos_charge, neg_charge, ISCalc), 'w')
    data.write("\n".join(['variable,value']) + "\n")

    for info in bindingInfo:
         data.write("\n".join(['%s,%s' % (info[0], info[1])]) + "\n")
    data.close()





def readPDB(fileName, *args, **kwargs):
    '''
    read pdb file and extract data. 
    '''

    pH = kwargs.get('pH', 7)

    pdbList = []
    resNumList = []
    residueNumDict = {}
    num_pos_neg_points = [[], []]


    with open(fileName) as data:
        count = 0
        patch_point_total = 0
        patchCoordsList = []
        patch_num = None

        for line in data:
            if line[0:6] == 'ATOM  ':
                atom_name = line[12:16].replace(' ','')
                resname = line[17:20].replace(' ','')
                chain = line[21:22].replace(' ','')
                resid = int(line[22:26].replace(' ',''))
                x = float(re.sub('[^-0-9.]', '', line[30:38]))
                y = float(re.sub('[^-0-9.]', '', line[38:46]))
                z = float(re.sub('[^-0-9.]', '', line[46:54]))

                try:
                    occupancy = float(line[54:60].replace(' ',''))
                    beta = float(line[60:66].replace(' ',''))
                except ValueError:
                    occupancy = None
                    beta = None

                pdbList.append([atom_name, resname, resid, x, y, z, 
                        occupancy, beta])


                if resname == 'ARG' and 'CZ' in atom_name and pH < 12:
                    num_pos_neg_points[0].append([resname, resid, 
                        atom_name])
                if resname == 'ASP' and 'CG' in atom_name and pH > 4:
                    num_pos_neg_points[1].append([resname, resid, 
                        atom_name])
                if resname == 'GLU' and 'CD' in atom_name and pH > 4.4:
                    num_pos_neg_points[1].append([resname, resid, 
                        atom_name])
                if resname == 'LYS' and 'NZ' in atom_name and pH < 10.4:
                    num_pos_neg_points[0].append([resname, resid, 
                        atom_name])

                if resname == 'HIS' and 'CE1' in atom_name and pH < 6:
                    num_pos_neg_points[0].append([resname, resid, 
                        atom_name])


                if resname not in residueNumDict and \
                        [resid, chain] not in resNumList:
                    residueNumDict[resname] = 1
                    resNumList.append([resid, chain])
                elif resname in residueNumDict and \
                        [resid, chain] not in resNumList:
                    residueNumDict[resname] += 1
                    resNumList.append([resid, chain])
                else:
                    continue


    data.close()



    if len(pdbList) != 0:
        return pdbList, residueNumDict, num_pos_neg_points



def outputData(fileName, *args, **kwargs):
    '''
    Use input data to run sequence binder and output protein
    residue info and binding info.

    '''

    pH = kwargs.get('pH', 7)
    conc_neg = kwargs.get('conc_neg', 'nul')
    conc_pos = kwargs.get('conc_pos', 'nul')
    neg_charge = kwargs.get('neg_charge', 1.0)
    pos_charge = kwargs.get('pos_charge', 1.0)
    res_counts = kwargs.get('res_counts', False)


    pdbList, residueNumDict, num_pos_neg_points = readPDB(fileName, pH=pH)

    fileName = fileName.split('.') 
    fileName = fileName[0]


    if res_counts == True: #save residue counts to file
        data = open('%s_residueCounts.csv' % (fileName), 'w')
        data.write("\n".join(['residue,count']) + "\n")
        for residue, count in residueNumDict.items():
            data.write("\n".join(['%s,%s' % (residue, count)]) + "\n")
        data.close()


    #calculate counterion binding
    calcMultiBinding(fileName, pdbList, num_pos_neg_points, pH=pH, 
            conc_neg=conc_neg, conc_pos=conc_pos, 
            neg_charge=neg_charge, pos_charge=pos_charge)





def main():


    try:
        usage = 'sequenceBinder.py [-h]'
        parser = argparse.ArgumentParser(description='Tool for binding '
                'counterions onto charged protein residues. '\
                '\nExample:'\
                '\npython sequenceBinder.py --fileName file.pdb --pH 7 '\
                '--conc_neg 0.1 --neg_charge 3 --pos_charge 1'\
                '\nOnly set concentration of positive or negative counterion.', 
                usage=usage, 
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        group = parser.add_argument_group('Options')
        group = parser.add_argument('-f', '--fileName', action='store', 
                default='file.pdb', help='input pdb list file name')
        group = parser.add_argument('-pH', '--pH', action='store', 
                type=float, default=7, help='input pH')
        group = parser.add_argument('-cn', '--conc_neg', action='store', 
                default='nul', 
                help='input concentration (Molar) of negative ion')
        group = parser.add_argument('-cp', '--conc_pos', 
                action='store', default='nul', 
                help='input concentration (Molar) of positive ion')
        group = parser.add_argument('-n', '--neg_charge', 
                action='store', type=float, default=1, 
                help='charge of negative ion')
        group = parser.add_argument('-p', '--pos_charge', 
                action='store', type=float, default=1, 
                help='charge of positive ion')
        group = parser.add_argument('-r', '--res_counts',
                action='store_true',
                help='include flag if residue counts are required')



        op = parser.parse_args()
    except argparse.ArgumentError:
        print ('ERROR - command line arguments are ill-defined, '\
                'please check the arguments')
        raise
        sys.exit(1)


    print (startTime)
       

    outputData(fileName=op.fileName, pH=op.pH, conc_neg=op.conc_neg, 
            conc_pos=op.conc_pos, neg_charge=op.neg_charge, 
            pos_charge=op.pos_charge, res_counts=op.res_counts)


    print(datetime.now() - startTime)


if __name__ == '__main__':
    main()
   

