#!/usr/bin/env python3

import sys
from itertools import islice
from py_stringmatching import *

MINR = -2

class SimList:
    rank1 = None
    rank2 = None
    cf = None
    
    def __init__(self, l1, l2, seed=None):
        assert(len(l1) == len(l2))
        elems1 = l1[:]
        elems2 = l2[:]
        
        #self._get_rank_list(elems1, elems2)
        dr1 = self._get_position_indixes(elems1)
        dr2 = self._get_position_indixes(elems2)
        self._get_rank(dr1, dr2, elems1, elems2)
        repeated1 = self._get_repetitions(self.rank1)
        repeated2 = self._get_repetitions(self.rank2)
        new_rank1 = self._get_avg_ranks(repeated1)
        new_rank2 = self._get_avg_ranks(repeated2)
        self._set_ranks(self.rank1, new_rank1)
        self._set_ranks(self.rank2, new_rank2)
        self.cf = self._compute_factor(repeated1, repeated2)

    def _get_position_indixes(self, orig_seq):
        dicc_rank = {}
        elems = sorted(orig_seq, reverse=True)
        n = len(elems)
        for i in range(n):
            e = elems[i]
            if e not in dicc_rank:
                dicc_rank[e] = i + 1
        return dicc_rank

    def _get_rank(self, dic_rank1, dic_rank2, elems1, elems2):
        self.rank1 = []
        self.rank2 = []
        assert(len(elems1) == len(elems2))
        n = len(elems1)
        for i in range(n):
            pos1 = dic_rank1[elems1[i]]
            pos2 = dic_rank2[elems2[i]]
            self.rank1.append(pos1)
            self.rank2.append(pos2)
    
    def spearmanr_rho(self):
        assert(len(self.rank1) == len(self.rank2))
        n = len(self.rank1)
        total = 0
        for i in range(n):
            total += (self.rank1[i] - self.rank2[i])**2
        return 1 - 6*(total + (1/12)*self.cf)  / (n*(n**2 -1))
    """
    def _get_rank_list(self, elems1, elems2):
        assert(len(elems1) == len(elems2))
        ls1 = sorted(elems1, reverse=True)
        ls2 = sorted(elems2, reverse=True)
        n = len(elems2)
        self.rank1 = []
        self.rank2 = []
        for i in range(n):
            pos1 = ls1.index(elems1[i]) + 1
            pos2 = ls2.index(elems2[i]) + 1
            self.rank1.append(pos1)
            self.rank2.append(pos2)
    """
    def _get_repetitions(self, ranks):
        repetead = {}
        for r in ranks:
            if r not in repetead:
                repetead[r] = 1
            else:
                repetead[r] += 1
        return repetead

    def _get_avg_ranks(self, rep):
        new_pos = {}
        for key in rep:
            n = rep[key]
            if n > 1:
                acc = key
                for i in range(n-1):
                    acc += key + (i+1)
                new_pos[key]  = acc/n
        return new_pos

    def _set_ranks(self, rank, ties_rank):
        n = len(rank)
        for i in range(n):
            pos = rank[i]
            if pos in ties_rank:
                rank[i] = ties_rank[pos]

    def _compute_factor(self, repeated1, repeated2):
        cf = 0
        for key in repeated1:
            if repeated1[key] > 1:
                cf += repeated1[key] * (repeated1[key]**2 - 1)
        for key in repeated2:
            if repeated2[key] > 1:
                cf += repeated2[key] * (repeated2[key]**2 - 1)
        return cf
    
def window_list(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = list(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + [elem,]
        yield result

def spearmanr(s1, s2):
    cl = SimList(s1, s2)
    sim = cl.spearmanr_rho()
    del cl
    return sim 

def sim_subseq(s1, s2):
    if (len(s1) == len(s2)):
        bestSim = spearmanr(s1, s2)
    else:
        if (len(s1) > len(s2)):
            a1 = s1
            a2 = s2
        else:
            a1 = s2
            a2 = s1
        n = len(a2)
        bestSim = MINR
        for subseq in window_list(a1, n):
            #print(subseq, a2)
            sim = spearmanr(subseq, a2)
            if sim > bestSim:
                bestSim = sim
    #stfidf = soft_tf_idf(s1, s2)
    #return max(bestSim*stfidf, 0)
    return max(bestSim, 0)

def sim_subseq_t_norm_l(s1, s2):
    stfidf = soft_tf_idf(s1, s2)
    if (len(s1) == len(s2)):
        bestSim = spearmanr(s1, s2)
        result = lukasiewicz_t_norm(stfidf, bestSim) 
    else:
        result = stfidf
        if (len(s1) > len(s2)):
            a1 = s1
            a2 = s2
        else:
            a1 = s2
            a2 = s1
        n = len(a2)
        for subseq in window_list(a1, n):
            sim = spearmanr(subseq, a2)
            result = lukasiewicz_t_norm(result, sim)
    return result

def sim_subseq_t_norm_h(s1, s2):
    stfidf = soft_tf_idf(s1, s2)
    if (len(s1) == len(s2)):
        bestSim = spearmanr(s1, s2)
        result = hamacher_t_norm(stfidf, bestSim) 
    else:
        if (len(s1) > len(s2)):
            a1 = s1
            a2 = s2
        else:
            a1 = s2
            a2 = s1
        n = len(a2)
        result = stfidf
        for subseq in window_list(a1, n):
            sim = spearmanr(subseq, a2)
            print(result, sim)
            result = hamacher_t_norm(result, sim)
    return result

def lukasiewicz_t_norm(a, b):
    return max(0, a+b-1)

def hamacher_t_norm(a, b):
    if a == 0 and b == 0:
        r = 0
    else:
        r = (a*b)/(a+b-a*b)
    return r

def soft_tf_idf(s1, s2):
    soft_tfidf = SoftTfIdf([s1, s2])
    r = soft_tfidf.get_raw_score(s1, s2)
    del soft_tfidf
    return r

def jaro_winkler(s1, s2):
    str1 = ''.join(e+" " for e in s1)
    str2 = ''.join(e+" " for e in s2)
    jw = JaroWinkler()
    r = jw.get_raw_score(str1, str2)
    del jw
    return r 

def sim_jaccard(l1, l2):
    s = set(l1)
    t = set(l2)
    return len(s & t)/len(s | t)

def main(*args):
    """
    print("-----------------")
    l1 = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    l3 = ["9", "8", "7", "6", "5", "4", "3", "2", "1"]
    l2 = ["2", "1", "3", "4", "3", "3", "8", "2", "3"]
    for s in window(l1, 3):
        print(s, s[0])
    r = spearmanr([1,2,3,4,5],[5,6,7,8,7])
    print(r)
    print(r.correlation)
    l1 = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    l2 = ["2", "1", "3", "4", "3", "3", "8", "2", "3"]
    l1 = get_rank_list(l1)
    l2 = get_rank_list(l2)
    l1 = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    l2 = ["2", "1", "3", "4", "3", "3", "8", "2", "3"]
    r = spearmanr(l1,l2)
    print(r)
    print(r.correlation)
    r = spearmanr(l1,l2)
    print(r)
    print(r.correlation)
    l1 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    l2 = [2, 1, 3, 4, 3, 3, 8, 2, 3]
    r = spearmanr(l1,l2)
    print(r)
    print(r.correlation)
    l1 = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
    l2 = ["3", "4", "5"]
    l3 = ["9", "8", "7", "6", "5", "4", "3", "2", "1"]

    l1 = ["1", "2", "3", "4", "2", "2", "4", "5"]
    l2 = ["3", "4", "5"]
    l3 = ["2", "4", "2", "3", "3", "3", "5"]
    s  = sim_subseq(l1, l2)
    print(l1, l2,s)
    s  = sim_subseq(l1, l3)
    print(l1, l3, s)
    s  = sim_subseq(l2, l3)
    print(l2, l3, s)

    cl = SimList(l2,["2", "2", "2"])
    print("testing s_rho_nomrmalized")
    z = cl.spearmanr_w_r()    
    print("Spearmanr without repetions", z)
   
    p1 = [43,96,74,38,35,43,22,56,35,80]
    p2 = [30,94,84,13,30,18,30,41,48,95]

    cl = SimList(p1, p2)
    print("testing s_rho_nomrmalized")
    z = cl.spearmanr_rho()    
    
    print(p1, p2, sim_subseq(p1, p2))
    print([1,2,3], [3,2,1], sim_subseq([1,2,3], [3,2,1]))
    print([1,2,3], [1,2,3], sim_subseq([1,2,3], [1,2,3]))
    """

    """
    l1 = ["1", "2", "3", "4", "2", "2", "4", "5"]
    l2 = ["3", "4", "5"]
    l3 = ["2", "4", "2", "1", "3", "3", "5"]
    """
    # Case 1 
    #l1 = ["4", "4", "3", "3", "3", "2", "2", "1"]
    #l2 = ["3", "2", "2"]
    #l3 = ["4", "3", "2", "1"]
    #l4 = ["4", "4", "3", "3", "3" ]
    #lt = [l1, l2, l3, l4]

    """
    l1 = ["2", "3", "4"]
    l2 = ["2", "4"]
    lt = [l1, l2]
    """
    """
    l1 = ["3", "2", "4"]
    l2 = ["4", "3", "2"]
    lt = [l1, l2]
    """
    """
    l1 = ["4", "2", "2", "2", "3"]
    l2 = ["2", "4", "2", "3", "2"]
    lt = [l1, l2]
    """
    l1 = ["2", "3", "3", "3", "4", "5"]
    l2 = ["4", "3", "3", "3", "2"]
    lt = [l1, l2]

    n = len(lt)
    jw = "Jaro Winkler\n"
    st = "Soft tf idf\n"
    sl = "Sequence with t-norm lukasiewicz\n"
    sp = "Sequence with t-norm product\n"
    for i in range(n):
        l1 = lt[i]
        for j in range(n):
            l2 = lt[j]
            jw += str((l1, l2))+"\t"
            st += str((l1, l2))+"\t"
            sl += str((l1, l2))+"\t"
            sp += str((l1, l2))+"\t"
            print("--------------------------------")
            print((l1, l2))
            print("Jaro Winkler", jaro_winkler(l1, l2))
            jw += str(jaro_winkler(l1, l2)) + "\t"
            print("soft_tf_idf", soft_tf_idf(l1, l2))
            st += str(soft_tf_idf(l1, l2)) + "\t"
            print("sim_subseq lukasiewicz", sim_subseq_t_norm_l(l1, l2))
            sl += str(sim_subseq_t_norm_l(l1, l2)) + "\t"
            print("sim_subseq tnorm", sim_subseq(l1, l2))
            sp += str(sim_subseq(l1, l2)) + "\t"
        jw += "\n"
        st += "\n"
        sl += "\n"
        sp += "\n"

    jw +="\n"+st+"\n"+sl+"\n"+sp+"\n"
    with open("results_caso4.csv", "w") as fd:
        fd.write(jw)

    p1 = [43,96,74,38,35,43,22,56,35,80]
    p2 = [30,94,84,13,30,18,30,41,48,95]

    cl = SimList(p1, p2)
    print("testing s_rho_nomrmalized")
    z = cl.spearmanr_rho()    
    print(z)
        
if __name__ == "__main__":
    main(*sys.argv[1:])
