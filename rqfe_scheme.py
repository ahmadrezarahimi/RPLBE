from models import CRS, HelperSecretKey, MasterPublicKey, PublicKey
import utils
from petrelic.multiplicative.pairing import G1, G2, GT, Bn

"""
RQFE Scheme


TODO: Add description of the scheme.
TODO: Make it compatible with the databases
"""
def setup(L):
    L = 16
    n1 = 8
    n2 = 5
    F = []
    for i in range(L):
        F.append([[G1.order().random() for _ in range(n2)] for _ in range(n1)])
    crs = CRS(L, n1, n2, F)
    return crs



def keygen(crs,l):
    s = [G1.order().random() for _ in range(crs.n1)]
    w = G1.order().random()
    temp= []
    for k in range(crs.L):
        sF = utils.vec_mat_mul(s,crs.F[k])
        temp.append([crs.t[i]**sF[i] for i in range(crs.n2)])
    v = [G2.neutral_element() for k in range(crs.L)]
    for k in range(crs.L):
        for i in range(crs.n2):
            v[k] = v[k] * temp[k][i]

    dk = [v[k] * (crs.gamma[k]**w) for k in range(crs.L)]
    dk_pk_term = [dk[i] if i != l else G2.neutral_element() for i in range(crs.L)]
    pk = PublicKey([crs.g1 ** s[i] for i in range(crs.n1)],crs.g1 ** w, dk_pk_term)
    sk = dk[l]
    return pk,sk

def aggregate(crs, pks):
    w = G1.neutral_element()
    for pk in pks:
        w = w * pk.w

    s = [G1.neutral_element() for _ in range(crs.n1)]
    for pk in pks:
        for i in range(crs.n1):
            s[i] = s[i] * pk.s[i]

    h1 = [G2.neutral_element() for _ in range(crs.L)]
    for i in range(crs.L):
        for j in range(crs.L):
            if i != j:
                h1[i] = h1[i] * pks[j].dk[i]

    h2 = [crs.gamma[i] for i in range(crs.L)]

    hsk = [HelperSecretKey(h1[k], h2[k], crs.F[k]) for k in range(crs.L)]
    mpk = MasterPublicKey(s, w, crs.t)
    return mpk, hsk


def encrypt(crs, mpk, m):
    x, y = m

    alpha = G1.order().random()
    g1 = G1.generator()
    g2 = G2.generator()

    p = G1.order()
    a = G1.order().random()
    b = G1.order().random()
    c = G1.order().random()
    d = G1.order().random()

    det_M = a.mod_mul(d, p).mod_sub(b.mod_mul(c, p), p)  # det(M) = ad - bc
    assert det_M != Bn(0)
    det_M_inv = det_M.mod_inverse(p)
    M = [[a, b], [c, d]]
    M_inv_transpose = [[d.mod_mul(det_M_inv, p), (-c).mod(p).mod_mul(det_M_inv, p)],
                       [(-b).mod(p).mod_mul(det_M_inv, p), a.mod_mul(det_M_inv, p)]]

    c1 = g1 ** alpha
    c2 = mpk.w ** alpha  # This is g1^(alpha*w)
    c3 = [(((g1 ** (M_inv_transpose[0][0].mod_mul(x[i], p))) * (mpk.s[i] ** alpha) ** M_inv_transpose[0][1]),
           (g1 ** (M_inv_transpose[1][0].mod_mul(x[i], p)) * (mpk.s[i] ** alpha) ** M_inv_transpose[1][1])) for i in
          range(crs.n1)]

    c4 = [(((g2 ** (M[0][0].mod_mul(y[i], p))) * (mpk.t[i] ** (-1)) ** M[0][1]),
           (g2 ** (M[1][0].mod_mul(y[i], p)) * (mpk.t[i] ** (-1)) ** M[1][1])) for i in range(crs.n2)]

    return (c1, c2, c3, c4)


def decrypt(crs, sk, hsk, c):
    c1, c2, c3, c4 = c
    d0 = c1.pair(hsk.h1 * sk)
    d1 = c2.pair(hsk.h2)

    # d0 = c1.pair(sk)

    d2 = GT.neutral_element()

    for i in range(crs.n1):
        for j in range(crs.n2):
            d2 = d2 * (c3[i][0].pair(c4[j][0]) * c3[i][1].pair(c4[j][1])) ** (hsk.F[i][j])

    m = d0 * d2 * (d1 ** (-1))
    return m
