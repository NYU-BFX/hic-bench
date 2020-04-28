import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri
import numpy as np
from scipy.sparse.linalg import eigsh
import imagehash
from PIL import Image


def similarity_scc(map1, map2, TAD):
    r = robjects.r

    r('''
        library(hicrep)
        hicrep_similarity<-function(map1, map2, TAD){
          scc = c()
          for(i in 1:nrow(TAD)){
            start = TAD[i, 2]
            end = TAD[i, 3]
            hic1 = data.frame('chr', 1:(end - start+1), 2:(end - start+2), map1[start:end,start:end])
            hic2 = data.frame('chr', 1:(end - start+1), 2:(end - start+2), map2[start:end,start:end])
            processed <- prep(hic1, hic2, 1, 0, 5)
            scc.out = get.scc(processed, 1, 5)
            scc=c(scc,scc.out$scc)
          }
          TAD2 = data.frame(TAD,scc)
          return(TAD2)
        }
    ''')

    rpy2.robjects.numpy2ri.activate()
    pandas2ri.activate()

    nr, nc = TAD.shape
    TAD_r = r.matrix(TAD, nrow=nr, ncol=nc)
    r.assign("TAD1", TAD_r)

    hicrep_similarity = r('hicrep_similarity')
    scc = hicrep_similarity(map1, map2, TAD_r)
    scc = pandas2ri.ri2py(scc)
    scc = scc.values
    return scc


def similarity_Laplacian(map1, map2, TAD):
    def get_Laplacian(M):
        S = M.sum(1)
        i_nz = np.where(S > 0)[0]
        S = S[i_nz]
        M = (M[i_nz].T)[i_nz].T
        S = 1 / np.sqrt(S)
        M = S * M
        M = (S * M.T).T
        n = np.size(S)
        M = np.identity(n) - M
        M = (M + M.T) / 2
        return M

    def evec_distance(v1, v2):
        d1 = np.dot(v1 - v2, v1 - v2)
        d2 = np.dot(v1 + v2, v1 + v2)
        if d1 < d2:
            d = d1
        else:
            d = d2
        return np.sqrt(d)

    def get_ipr(evec):
        ipr = 1.0 / (evec * evec * evec * evec).sum()
        return ipr

    def get_reproducibility(M1, M2, num_evec=20):
        k1 = np.sign(M1).sum(1)
        d1 = np.diag(M1)
        kd1 = ~((k1 == 1) * (d1 > 0))
        k2 = np.sign(M2).sum(1)
        d2 = np.diag(M2)
        kd2 = ~((k2 == 1) * (d2 > 0))
        iz = np.nonzero((k1 + k2 > 0) * (kd1 > 0) * (kd2 > 0))[0]
        M1b = (M1[iz].T)[iz].T
        M2b = (M2[iz].T)[iz].T

        i_nz1 = np.where(M1b.sum(1) > 0)[0]
        i_nz2 = np.where(M2b.sum(1) > 0)[0]
        M1b_L = get_Laplacian(M1b)
        M2b_L = get_Laplacian(M2b)

        a1, b1 = eigsh(M1b_L, k=num_evec, which="SM")
        a2, b2 = eigsh(M2b_L, k=num_evec, which="SM")

        b1_extend = np.zeros((np.size(M1b, 0), num_evec))
        b2_extend = np.zeros((np.size(M2b, 0), num_evec))
        for i in range(num_evec):
            b1_extend[i_nz1, i] = b1[:, i]
            b2_extend[i_nz2, i] = b2[:, i]

        ipr_cut = 5
        ipr1 = np.zeros(num_evec)
        ipr2 = np.zeros(num_evec)
        for i in range(num_evec):
            ipr1[i] = get_ipr(b1_extend[:, i])
            ipr2[i] = get_ipr(b2_extend[:, i])

        b1_extend_eff = b1_extend[:, ipr1 > ipr_cut]
        b2_extend_eff = b2_extend[:, ipr2 > ipr_cut]
        num_evec_eff = min(np.size(b1_extend_eff, 1), np.size(b2_extend_eff, 1))

        evd = np.zeros(num_evec_eff)
        for i in range(num_evec_eff):
            evd[i] = evec_distance(b1_extend_eff[:, i], b2_extend_eff[:, i])

        Sd = evd.sum()
        l = np.sqrt(2)
        evs = abs(l - Sd / num_evec_eff) / l
        return evs

    def similarity(map1, map2, TAD_f):
        e = []

        for i in TAD_f:
            n1 = int(i[1])
            n2 = int(i[2])
            x1 = map1[n1:(n2 + 1), n1:(n2 + 1)]
            x2 = map2[n1:(n2 + 1), n1:(n2 + 1)]
            e.append(get_reproducibility(x1, x2))

        e = np.array(e)
        e = e.reshape(e.shape[0], 1)
        TAD_f = np.concatenate((TAD_f, e), 1)
        return TAD_f
    Laplacian = similarity(map1, map2, TAD)
    return Laplacian


def hash_similarity(map1, map2, TAD):
    e_ahash = []
    for i in TAD:
        n1 = int(i[1])
        n2 = int(i[2])
        x1 = map1[n1:(n2 + 1), n1:(n2 + 1)]
        x2 = map2[n1:(n2 + 1), n1:(n2 + 1)]
        x1 = x1.astype(np.uint8)
        x2 = x2.astype(np.uint8)
        x2 = x2 / np.mean(x2)
        img1 = Image.fromarray(x1)
        img2 = Image.fromarray(x2)
        hash1 = imagehash.average_hash(img1, 16)
        hash2 = imagehash.average_hash(img2, 16)
        d = 1 - (hash2 - hash1) / 256
        e_ahash.append(d)
    e_ahash = np.array(e_ahash)
    e_ahash = e_ahash.reshape(e_ahash.shape[0], 1)
    TAD = np.concatenate((TAD, e_ahash), 1)
    return TAD


