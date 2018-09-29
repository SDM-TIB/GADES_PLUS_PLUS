#!/usr/bin/env python3

import requests, json, os, rdflib
from sem_sim_list import *
from statistics import mean
 
URL_SERVICE = "http://localhost:9000/similarity/initialize?model_1="
SERVICE_JACCARD = "http://localhost:9000/similarity/jaccard"
SERVICE_GADES = "http://localhost:9000/similarity/gades"

STAGE = "http://project-iasis.eu/vocab/tumorStageList"
DIAG_AGE = "http://project-iasis.eu/vocab/diagnosisAge"
BIOPSY = "http://project-iasis.eu/vocab/biopsyList"
ECOG = "http://project-iasis.eu/vocab/ecogList"
MONTHS_IN_TRE = "http://project-iasis.eu/vocab/monthsInTreatment"
ALPHA=0.6
MATRIX_FILE = "matrix.txt"
ENTITIES_FILE = "entities.txt"

SIM_OUTPUT = 'tmp/sim_output'
SIM_INPUT = 'tmp/sim_input'

class Entity:    
    ids = None
    monthsInTreatment = None

    ecog = []
    biopsy = []
    stage = []
    drugs = []
    comorbidities = []
    chemotherapyList = []
    tkiList = []
    immunotherapyList = []

    immunotherapy = False
    chemotherapy = False
    tki = False
    surgery = False
    antiangiogenic = False
    radiationtherapy = False
    
    systemicProgression = False
    localProgression = False
    brainMetastasis = False

    def __init__(self, id):
        self.ids = id
        
    def __str__(self):
        return "++++++++++++\nPatient: "+self.ids+"\n"+str(self.comorbidities)+"\t"+str(self.drugs)+"\t"+str(self.chemotherapy)+"\t"+str(self.chemotherapyList)+"\t"+str(self.tki)+"\t"+str(self.tkiList)+"\t"+str(self.immunotherapy)+"\t"+str(self.immunotherapyList)+"\t"+str(self.antiangiogenic)+"\t"+str(self.radiationtherapy)+"\t"+str(self.surgery)+"\t"+str(self.systemicProgression)+"\t"+str(self.localProgression)+"\t"+str(self.brainMetastasis)+"\t"+str(self.ecog)+"\t"+str(self.stage)+"\t"+str(self.biopsy)+"\t"+str(self.monthsInTreatment)

def load_nt_file(graph, filename):
    graph.parse(filename, format="nt")
    print("File with NT "+filename+" loaded!")

def load_gades_model(path_to_file_w_model):
    path_model = os.path.abspath(path_to_file_w_model)
    print("Loading the model, waiting for the GADES service")
    model = URL_SERVICE+path_model
    print(model)
    print("----------")
    r = requests.get(url = model)
    print("Respond of the GADES service: %s"%r)

"""    
def strip_entities(lentities):
    slentities = []
    for e in lentities:
        tok = e.split("/")
        slentities.append(tok[-1])
    return slentities
"""

def run_gades(file_pairs):
    with open(file_pairs) as fd:
        entities_to_compare = json.loads( fd.read() )
    """
    print("Computing the pairwise similarities, waiting for the GADES service")
    r = requests.post(SERVICE_JACCARD, json=entities_to_compare)
    print("Respond of the GADES service: %s"%r.text)
    djson =  r.json() 
    print("Results of JACCARD")
    for d in djson:
        print("URI 1", d['uri1'])
        print("URI 2", d['uri2'])
        print("Sim jaccard ", d['value'])
        print("--------------------")
    """
    print("Computing the pairwise similarities, waiting for the GADES service")
    r = requests.post(SERVICE_GADES, json=entities_to_compare)
    #print("Respond of the GADES service: %s"%r.text)
    print("Respond of the GADES service: %s "%r)
    djson =  r.json() 
    print("Results of GADES")
    """
    for d in djson:
        print("URI 1", d['uri1'])
        print("URI 2", d['uri2'])
        print("Sim gades ", d['value'])
        print("--------------------")
    """
    return djson

def load_molecules(filename):
    rdfGraph = rdflib.Graph()
    load_nt_file(rdfGraph, filename)
    return rdfGraph

def get_elements_to_compare(list_json):
    elements = set()
    for djson in list_json:
        for d in djson:
            elements.add(d['uri1'])
            elements.add(d['uri2'])
    lenities = {}
    for e in elements:
        oe = Entity(e) 
        lenities[e] = oe
    return lenities

def sim_dates(d1, d2):
    return 1- abs(d1-d2)/max(d1, d2)

def string_to_list(s):
    if s == "[]":
        return []
    else:
        ns = s[1:-1]
        return ns.split(";")

def tnorm(a, b):
    return a*b

def compute_average_set_sim(l1, l2, d_sim):
    total = []
    for s1 in l1:
        for s2 in l2:
            if s1 != s2:
                sim = d_sim[ (s1, s2) ]
                total.append(sim)
    return mean(total)

def get_entity_cancer_drugs(e):
    drugs = set()
    for c in e.chemotherapyList:
        drugs.add(c)
    for c in e.immunotherapyList:
        drugs.add(c)
    for c in e.tkiList:
        drugs.add(c)
    return drugs

def get_treatments_set(e):
    l = []
    if e.immunotherapy:
        l.append("immunotherapy")
    if e.chemotherapy:
        l.append("chemotherapy")
    if e.tki:
        l.append("tki")
    if e.surgery:
        l.append("surgery")
    if e.antiangiogenic:
        l.append("antiangiogenic")
    if e.radiationtherapy:
        l.append("radiationtherapy")
    return l

def get_progression_set(e):
    l = []
    if e.systemicProgression:
        l.append("systemicProgression")
    if e.localProgression:
        l.append("localProgression")
    if e.brainMetastasis:
        l.append("brainMetastasis")
    return l
    
def similarity(e1, e2, sim_comorbidities, sim_drugs, sim_drugs_c):
    if e1.ids == e2.ids:
        return 1.0

    if len(e1.ecog) != 0 and len(e2.ecog) != 0:
        if len(e1.ecog) != 1 and len(e2.ecog) != 1:
            ecog = sim_subseq(e1.ecog, e2.ecog)
        else:
            if len(e1.ecog) == 1 and len(e2.ecog) == 1:
                #print("Warning using soft tf-idf", (e1.ecog, e2.ecog))
                ecog = jaro_winkler(e1.ecog, e2.ecog)
            else:
                ecog = 0.0
    else:
        ecog = 0.0

    if len(e1.stage) != 0 and len(e2.stage) != 0:
        if len(e1.stage) != 1 and len(e2.stage) != 1:
            stage = sim_subseq(e1.stage, e2.stage)
        else:
            if len(e1.stage) == 1 and len(e2.stage) == 1:
                #print("Warning using jaro winkler", (e1.stage, e2.stage))
                stage = jaro_winkler(e1.stage, e2.stage)
            else:
                stage = 0.0
    else:
        stage = 0.0
        
    if len(e1.biopsy) != 0 and len(e2.biopsy) != 0:
        if len(e1.biopsy) != 1 and len(e2.biopsy) != 1:
            biopsy = sim_subseq(e1.biopsy, e2.biopsy)
        else:
            if len(e1.biopsy) == 1 and len(e2.biopsy) == 1:
                #print("Warning using jaro winkler", (e1.biopsy, e2.biopsy))
                biopsy = sim_jaccard(e1.biopsy, e2.biopsy)
            else:
                biopsy = 0.0
    else:
        biopsy = 0.0

    if e1.monthsInTreatment > 0.0 or e2.monthsInTreatment > 0.0:
        months = sim_dates(e1.monthsInTreatment, e2.monthsInTreatment)
    else:
        months = 0

    if (len(e1.comorbidities) > 0) and  (len(e2.comorbidities) > 0):
        sim_pair_comorbidities = compute_average_set_sim(e1.comorbidities, e2.comorbidities, sim_comorbidities)
    else:
        sim_pair_comorbidities = 0.0

    if (len(e1.drugs) > 0) and  (len(e2.drugs) > 0):
        sim_pair_drugs = compute_average_set_sim(e1.drugs, e2.drugs, sim_drugs)
    else:
        sim_pair_drugs = 0.0

    cancer_drugs_e1 = get_entity_cancer_drugs(e1)
    cancer_drugs_e2 = get_entity_cancer_drugs(e2)

    if (len(cancer_drugs_e1) > 0) and (len(cancer_drugs_e2) > 0):
        sim_pair_drugs_c = compute_average_set_sim(cancer_drugs_e1, cancer_drugs_e2, sim_drugs_c)
    else:
        sim_pair_drugs_c = 0.0

    treatments_e1 = get_treatments_set(e1)
    treatments_e2 = get_treatments_set(e2)
    treatments_pair_sim = sim_jaccard(treatments_e1, treatments_e2)

    progression_e1 = get_progression_set(e1)
    progression_e2 = get_progression_set(e2)
    progression_pair_sim = sim_jaccard(progression_e1, progression_e2)
    
    evolution = [ecog, stage, biopsy, sim_pair_comorbidities, sim_pair_drugs, sim_pair_drugs_c, treatments_pair_sim,  progression_pair_sim]
    
    return ((1.0-ALPHA)*mean(evolution) + months*ALPHA)/2.0
        
"""
def similarity(e1, e2):
    if e1.ids == e2.ids:
        return 1.0
    
    if len(e1.ecog) != 0 and len(e2.ecog) != 0:
        if len(e1.ecog) != 1 and len(e2.ecog) != 1:
            ecog = sim_subseq(e1.ecog, e2.ecog)
        else:
            if len(e1.ecog) == 1 and len(e2.ecog) == 1:
                #print("Warning using soft tf-idf", (e1.ecog, e2.ecog))
                ecog = jaro_winkler(e1.ecog, e2.ecog)
            else:
                ecog = 0.0
    else:
        ecog = 0.0
    if len(e1.stage) != 0 and len(e2.stage) != 0:
        if len(e1.stage) != 1 and len(e2.stage) != 1:
            stage = sim_subseq(e1.stage, e2.stage)
        else:
            if len(e1.stage) == 1 and len(e2.stage) == 1:
                #print("Warning using jaro winkler", (e1.stage, e2.stage))
                stage = jaro_winkler(e1.stage, e2.stage)
            else:
                stage = 0.0
    else:
        stage = 0.0
    if len(e1.biopsy) != 0 and len(e2.biopsy) != 0:
        biopsy = sim_jaccard(e1.biopsy, e2.biopsy)
    else:
        biopsy = 0.0
    if e1.monthsInTreatment > 0.0 or e2.monthsInTreatment > 0.0:
        months = sim_dates(e1.monthsInTreatment, e2.monthsInTreatment)
    else:
        months = 0
    #print(("ecog", ecog), ("stage", stage), ("biopsy", biopsy), ("dage", dage), ("months", months))
    #print(("ecog", ecog), ("stage", stage), ("months", months))
    #return tnorm( tnorm( tnorm( tnorm( tnorm(ecog, stage), biopsy ), dage), months ), gades )
    #return (ecog + stage + biopsy + dage + months + gades)/6
    #if e1.ids == 'http://project-iasis.eu/LCPatient/LCpatient_307cc02266869e778ed61b0341c15f51' and e2.ids == 'http://project-iasis.eu/LCPatient/LCpatient_a1d557bced94dad8691ff63b9f753626':
    
    #print(e1.ids, e2.ids)
    #print((e1.monthsInTreatment, e2.monthsInTreatment))
    #print((ecog, stage, months))
    #print((1.0-ALPHA)*((ecog + stage)/2), (months*ALPHA))
    #print(((1.0-ALPHA)*((ecog + stage)/2) + months*ALPHA)/2)
    #print("---")
    
    return ((1.0-ALPHA)*((ecog + stage)/2) + months*ALPHA)/2
    """
    
def print_entities(d_entities):
    for k in d_entities:
        #print("-----------------")
        entity = d_entities[k]
        #print(entity.ids, entity.ecog, entity.biopsy, entity.stage, entity.diagnosisAge, entity.monthsInTreatment)

def create_result_csv_file(results, outfilename):
    with open(outfilename, "w") as fd:
        for (e1, e2, sim) in results:
            fd.write(e1+'\t'+e2+'\t'+str(sim)+'\n')
"""
def create_semEP_node_input(elements, results):
    index_sim = {}
    for (e1, e2, sim) in results:
        key = e1+"---"+e2
        index_sim[key] = sim  
        #fd.write(e1+'\t'+e2+'\t'+str(sim)+'\n')
    with open("entities_gades.txt", "w") as fm:
        fm.write(str(len(elements))+"\n")
        for e in elements:
            fm.write(elements[e].ids+"\n")
    with open("sim_gades_matrix.txt", "w") as fm:
        s = ""
        fm.write(str(len(elements))+"\n")
        for e1 in elements:
            for e2 in elements:
                key = elements[e1].ids+"---"+elements[e2].ids
                sim = index_sim[key]
                s += str(sim)+" "
            s = s[:-1]+"\n"
        fm.write(s)


def get_short_id(e):
    tok = e.split("/")
    return tok[-1]
"""

def create_semEP_node_input(results):
    index_sim = {}
    elements_set = set() 
    for (e1, e2, sim) in results:
        key = e1+"---"+e2
        index_sim[key] = sim
        key = e2+"---"+e1
        index_sim[key] = sim
        elements_set.add(e1)
        elements_set.add(e2)
    #elements = strip_entities( list(elements_set) )
    elements = list(elements_set) 
    with open(ENTITIES_FILE, "w") as fm:
        fm.write(str(len(elements))+"\n")
        for e in elements:
            fm.write(e+"\n")
    with open(MATRIX_FILE, "w") as fm:
        s = ""
        fm.write(str(len(elements))+"\n")
        for e1 in elements:
            for e2 in elements:
                key = e1+"---"+e2
                sim = index_sim[key]
                s += str(sim)+" "
            s = s[:-1]+"\n"
        fm.write(s)

def run_gades_for_each_pairwise_file(path_to_file_w_model, list_file_pairs):
    list_json = []
    load_gades_model(path_to_file_w_model)
    for file_pairs in list_file_pairs:
        d_json = run_gades(file_pairs)
        list_json.append(d_json)
    return list_json
    
def run_gades_plus_plus(path_to_file_w_model, list_file_pairs, file_w_list):
    list_json = run_gades_for_each_pairwise_file(path_to_file_w_model, list_file_pairs)
    #print(list_json)
    elements = get_elements_to_compare(list_json)
    #print(elements)
    graph = load_molecules(file_w_list)
    get_entity_data(graph, elements)
    #print_entities(elements)
    results = compute_all_similarities(elements, list_json)
    create_result_csv_file(results)
    create_semEP_node_input(elements, results)

############################################################################



def compute_all_similarities(d_entities, sim_comorbidities, sim_drugs, sim_drugs_c):
    #patient_id_list = list(d_entities.keys())
    results = []
    cont = 0
    n = len(d_entities)
    for i in range(n):
        e1 = d_entities[i]
        for j in range(i, n):
            if cont > 0 and cont % 10 == 0:
                print("\nNumber of pairs so far: "+str(cont))
            e2 = d_entities[j]
            sim = similarity(e1, e2, sim_comorbidities, sim_drugs, sim_drugs_c)
            #assert( sim  == similarity(e2, e1) )
            #print("Entity 1")
            #print(e1.ids, e1.ecog, e1.biopsy, e1.stage, e1.diagnosisAge, e1.monthsInTreatment)
            #print(e1.ids, e1.ecog, e1.stage, e1.monthsInTreatment)
            #print("Entity 2")
            #print(e2.ids, e2.ecog, e2.biopsy, e2.stage, e2.diagnosisAge, e2.monthsInTreatment)
            #print(e2.ids, e2.ecog, e2.stage, e2.monthsInTreatment)
            #print("URI 1", uri1)
            #print("URI 2", uri2)
            #print("Sim Gades++ ", sim)
            #print("--------------------")
            results.append((e1.ids, e2.ids, sim))
            cont += 1 
    return results


def get_patient_data(graph):
    d_entities = {}
    for (s, p, o) in graph:
        subject = str(s)
        if subject not in d_entities:
            d_entities[subject] = Entity(subject)
        if str(p) == STAGE:
            entity = d_entities[subject]
            if str(o) != 'null':
                entity.stage = string_to_list(str(o))
            else:
                entity.stage = []
        elif str(p) == BIOPSY:
            entity = d_entities[subject]
            if str(o) != 'null':
                entity.biopsy = string_to_list(str(o))
            else:
                entity.biopsy = []
        elif str(p) == ECOG:
            entity = d_entities[subject]
            if str(o) != 'null':
                entity.ecog = string_to_list(str(o))
            else:
                entity.ecog = []
        elif str(p) == MONTHS_IN_TRE:
            entity = d_entities[subject]
            if str(o) != 'null':
                entity.monthsInTreatment = max(0, float(str(o)))
            else:
                entity.monthsInTreatment = 0.0
        else:
            #print(p)
            pass
    return d_entities

def call_ulms_similarity(f_input, f_output):
    current_path = os.path.dirname(os.path.realpath(__file__))
    print("Current Path: "+current_path)
    input_abspath = os.path.join(current_path, f_input)
    ouput_abspath = os.path.join(current_path, f_output)
    commd = "umls-similarity.pl --matrix --infile "+input_abspath+" --measure vector --getcache "+ouput_abspath
    print("Command to execute "+commd)
    return os.system(commd)

def compute_umls_similarities(l_entities, flag):
    f_input = SIM_INPUT+"_"+flag
    f_output = SIM_OUTPUT+"_"+flag
    with open(f_input, "w") as fd:
        for l in l_entities:
            fd.write(l+"\n")
    fd.close()
    print("Input file created")
    print("Computing the similarities between "+flag)
    r = call_ulms_similarity(f_input, f_output)
    print("UMLS sim result ",r)
    sim_dic = {}
    with open(f_output) as fd:
        for line in fd:
            tok = line.rstrip().split("<>")
            sim = abs(float(tok[0]))
            cui1 = tok[1]
            cui2 = tok[2]
            sim_dic[(cui1, cui2)] = sim  
    return sim_dic

def get_all_comorbidities(d_entities):
    comorbidities = set()
    for e in d_entities:
        for c in e.comorbidities:
            comorbidities.add(c)
    return comorbidities

def get_all_drugs(d_entities):
    drugs = set()
    for e in d_entities:
        for c in e.drugs:
            drugs.add(c)
    return drugs

def get_all_cancer_drugs(d_entities):
    drugs = set()
    for e in d_entities:
        for c in e.chemotherapyList:
            drugs.add(c)
        for c in e.immunotherapyList:
            drugs.add(c)
        for c in e.tkiList:
            drugs.add(c)
    return drugs

def load_patients_from_file(filename):
    with open(filename) as fd:
        header = fd.readline().rstrip().split("\t")
        print(header)
        print("Number of concepts: "+str(len(header)))
        d_entities = []
        for line in fd:
            tok = line.rstrip().split("\t")
            #assert( len(header) == len(tok) )
            if len(header) != len(tok):
                print(len(header), len(tok))
                sys.exit(1)
            entity = Entity(tok[0])
            entity.comorbidities = string_to_list(tok[1])
            entity.drugs = string_to_list(tok[2])
            entity.chemotherapy = True if tok[3] == "TRUE" else False
            entity.chemotherapyList = string_to_list(tok[4])
            entity.tki = True if tok[5] == "TRUE" else False
            entity.tkiList = string_to_list(tok[6])
            entity.immunotherapy = True if tok[7] == "TRUE" else False
            entity.immunotherapyList = string_to_list(tok[8])
            entity.antiangiogenic = True if tok[9] == "TRUE" else False
            entity.radiationtherapy = True if tok[10] == "TRUE" else False
            entity.surgery = True if tok[11] == "TRUE" else False
            entity.systemicProgression = True if tok[12] == "TRUE" else False
            entity.localProgression = True if tok[13] == "TRUE" else False
            entity.brainMetastasis = True if tok[14] == "TRUE" else False
            entity.ecog = string_to_list(tok[15])
            entity.stage = string_to_list(tok[16])
            entity.biopsy = string_to_list(tok[17])
            entity.monthsInTreatment = int(tok[18])
            print(entity)
            d_entities.append(entity)
        return d_entities
        
#def run_all_gades_plus_plus(nt_filename):
def run_all_gades_plus_plus(graph):
    #patients = get_data_patients(args[0])
    #graph = load_molecules(nt_filename)
    d_patients = get_patient_data(graph)
    #print("Number of Patients "+str(len(d_patients)))
    print("Computation of all pairs of similarity")
    similarities = compute_all_similarities(d_patients)
    #print(similarities)
    #create_result_csv_file(similarities, args[1])
    print("Building the semEP-Node input files")
    create_semEP_node_input(similarities)
    print("Done Gades Plus Plus")

def main(*args):
    """"
    #patients = get_data_patients(args[0])
    graph = load_molecules(args[0])
    d_patients = get_patient_data(graph)
    print("Number of Patients "+str(len(d_patients)))
    print("Computation of all pairs of similarity")
    similarities = compute_all_similarities(d_patients)
    #print(similarities)
    #create_result_csv_file(similarities, args[1])
    print("Building the semEP-Node input files")
    create_semEP_node_input(similarities)
    """
    d_entities = load_patients_from_file(args[0])
    comorbidities = get_all_comorbidities(d_entities)
    len("Comorbidities "+str(len(comorbidities)))
    print(comorbidities)
    sim_comorbidities = compute_umls_similarities(comorbidities, "Comorbidities")
    print(sim_comorbidities)
    
    drugs = get_all_drugs(d_entities)
    len("Drugs "+str(len(drugs)))
    print(drugs)
    sim_drugs = compute_umls_similarities(drugs, "Drugs")
    print(sim_drugs)

    drugs_c = get_all_cancer_drugs(d_entities)
    len("Cancer drugs "+str(len(drugs_c)))
    print(drugs_c)
    sim_drugs_c = compute_umls_similarities(drugs_c, "Cancer_Drugs")
    print(sim_drugs_c)
    similarities = compute_all_similarities(d_entities, sim_comorbidities, sim_drugs, sim_drugs_c)
    print(similarities)
    print("Building the semEP-Node input files")
    create_semEP_node_input(similarities)
    
    print("Done Gades Plus Plus")

if __name__ == "__main__":
    main(*sys.argv[1:])
