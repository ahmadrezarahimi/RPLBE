import math
from models import CRS,HelperSecretKeyRPLBE,MasterPublicKey,PublicKey
import utils
from petrelic.multiplicative.pairing import G1, G2,Bn


"""
RPLBE Scheme

This is the optimized RPLBE scheme represented in the paper.
"""

def setup(L):
    """
    Initalizes and returns a CRS object for the given parameter L representing the maximum number of slots.

    :param L: The value for L.
    :return: A CRS object with the specified parameters.
    """
    n1 = 2*int(math.sqrt(L))
    n2 = int(math.sqrt(L))+1
    crs = CRS(L, n1, n2, None)

    return crs
    # Create a new connection to the database
def aggregate(crs, pks):
    """
    Aggregate the given parameters to compute the master public key and helper secret keys.

    :param crs: The common reference string (CRS)
    :param pks: A list of public keys to be aggregated.
    :return: A tuple containing the master public key (mpk) and a list of helper secret keys (aux).
    """
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

    aux = [HelperSecretKeyRPLBE(h1[k], h2[k]) for k in range(crs.L)]
    mpk = MasterPublicKey(s, w, crs.t)
    return mpk, aux


def keygen(crs,l):
    """
    :param crs: The common reference string
    :param l: The index of the secret key to be generated
    :return: The public key and secret key
    """
    s_l = [G1.order().random() for _ in range(crs.n1)]
    w = G1.order().random()
    dk_l = [G2.neutral_element() for _ in range(crs.L)]

    sqrt_L = int(math.sqrt(crs.L))

    for k in range(crs.L):
        k1,k2 = utils.encode_number(k,crs.L)
        dk_l[k] = crs.t[0]**s_l[k1] * crs.t[k2+1]**s_l[k1+sqrt_L]*crs.gamma[k]**w
    dk_pk_term = [dk_l[i] if i != l else G2.neutral_element() for i in range(crs.L)]
    pk = PublicKey([crs.g1 ** s_l[i] for i in range(crs.n1)],crs.g1 ** w, dk_pk_term)
    sk = dk_l[l]
    return pk,sk


def encrypt(crs, mpk, m):
    """
    Encrypts a message using the given RPLBE parameters

    :param crs: The common reference string
    :type crs: CRS

    :param mpk: The master public key
    :type mpk: MasterPublicKey

    :param m: The message to be encrypted from {1,2}
    :type m: int

    :return: The encrypted ciphertext
    :rtype: tuple
    """
    x,y = utils.Z(m,crs.L)

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


def decrypt(crs, sk,ind, hsk, c):
    """
    :param crs: common reference string
    :param sk: secret key
    :param ind: index of the decryptor
    :param hsk: helper secret key
    :param c: encrypted message
    :return: decrypted message

    Decrypts the given encrypted message using the provided parameters and secret key.

    """
    c1, c2, c3, c4 = c
    d0 = c1.pair(hsk.h1 * sk)
    d1 = c2.pair(hsk.h2)

    i1,i2 = utils.encode_number(ind,crs.L)

    sqrt_L = int(math.sqrt(crs.L))

    d2 = (c3[i1][0].pair(c4[0][0]) * c3[i1][1].pair(c4[0][1])*
          c3[i1+sqrt_L][0].pair(c4[i2+1][0]) * c3[i1+sqrt_L][1].pair(c4[i2+1][1]))

    m = d0 * d2 * (d1 ** (-1))
    return m







