
from petrelic.multiplicative.pairing import G1, G2, GT,G1Element,G2Element

class CRS:
    """
    Initialize a CRS object.

    :param L: int, number of elements in gamma.
    :param n1: int, number of elements in g1.
    :param n2: int, number of elements in g2.
    :param F: Optional, default value is None. Any value to set F.
    :param set_gamma: Optional, default value is True. If True, gamma will be set.
    :param set_t: Optional, default value is True. If True, t will be set.
    :param set_F: Optional, default value is False. If True, F will be set.
    """
    def __init__(self,L,n1,n2,F = None,set_gamma = True,set_t = True,set_F = False):
        self.L = L
        self.n1 = n1
        self.n2 = n2
        self.g1 = G1.generator()
        self.g2 = G2.generator()

        if set_gamma:
            self.gamma = [self.g2**G2.order().random() for _ in range(L)]
        if set_t:
            self.t = [self.g2**G2.order().random() for _ in range(n2)]
        if set_F:
            self.F = F
    def set_gamma(self,gamma):
        self.gamma = gamma
    def set_t(self,t):
        self.t = t


class PublicKey:
    """
    Initializes a new instance of the PublicKey class.

    :param s: The 's' parameter for the public key.
    :type s: [G1Element]
    :param w: The 'w' parameter for the public key.
    :type w: Z_p element
    :param dk: The 'dk' parameter for the public key.
    :type dk: [G2Element]
    """
    def __init__(self,s,w,dk):
        self.s = s
        self.w = w
        self.dk = dk
        

class SecretKey:
    def __init__(self) -> None:
        pass

class MasterPublicKey:
    """
    MasterPublicKey represents a master public key.
    Initialize a MasterPublicKey instance with the given parameters.

    Parameters:
        - s [G1Element]: The s parameter.
        - w G1Element: The w parameter.
        - t [G2Element]: The t parameter.
    """
    def __init__(self,s,w,t):
        self.s = s
        self.w = w
        self.t = t
        

class HelperSecretKey:
    def __init__(self,h1,h2,F):
        self.h1 = h1
        self.h2 = h2
        self.F = F


class HelperSecretKeyRPLBE:
    """
    A class representing helper secret keys in the RPLBE system.

    :param h1: The first helper secret key parameter.
    :type h1: G2Element
    :param h2: The helper second secret key parameter.
    :type h2: G2Element
    """
    def __init__(self, h1, h2):
        self.h1 = h1
        self.h2 = h2


