import sqlite3
from functools import reduce
from petrelic.multiplicative.pairing import G1,G2,GT, Bn, G1Element, G2Element
from models import CRS, MasterPublicKey, HelperSecretKeyRPLBE
import msgpack
import matplotlib.pyplot as plt
import time
import pandas as pd

"""
Some utility functions that are used in the RPLBE scheme.
"""
def vec_mat_mul(a,M):

    """
    This function computes the vector matrix multiplication aM = b
    :param vector s:
    :param matrix f:
    :return:
    """
    n1 = len(M)
    n2 = len(M[0])

    b = [Bn(0) for _ in range(n2)]
    p = G1.order()

    for i in range(n1):
        for j in range(n2):
            b[j] = b[j].mod_add(a[i].mod_mul(M[i][j],p),p)
    return b

def inner_prod(a,b):
    """
    This function computes the inner product of two vectors
    :param vector a:
    :param vector b:
    :return:
    """
    p = G1.order()
    return reduce(lambda x,y: x.mod_add(y,p),[a[i].mod_mul(b[i],p) for i in range(len(a))])

def encode_number(k,L):
    """
    Encodes a number using given parameters.

    :param k: An integer representing the number to be encoded.
    :param L: An integer representing a perfect square.
    :return: A tuple (k1, k2) representing the encoded number.
    """
    # Assume L is a perfect square
    k_base_1 = k+1
    sqrt_L = int(L**0.5)
    k1 =  (k_base_1 // sqrt_L) if (k_base_1 % sqrt_L == 0)  else (k_base_1 // sqrt_L) + 1
    k2 =  k_base_1 % sqrt_L if (k_base_1 % sqrt_L != 0) else sqrt_L
    return (k1-1,k2-1)



def Z(m, L):
    """
    :param m: The value of m
    :param L: The value of L
    :return: A tuple containing the calculated values of x and y

    This method calculates the values of x and y based on the input parameters m and L.

    Example Usage:

    Z(2, 5) ouputs:
    ([Bn(0), Bn(2), Bn(2), Bn(2), Bn(2), Bn(2), Bn(2), Bn(2), Bn(2)], [Bn(1), Bn(1), Bn(1), Bn(1), Bn(1), Bn(1), Bn(1), Bn(1), Bn(1)])

    """
    sqrt_L = int(L**0.5)

    # Calculate v_tilde and v_hat
    v_tilde = [Bn(0)] + [Bn(m)] * (sqrt_L - 1)
    v_hat = [Bn(0)] * sqrt_L
    v_hat[0] = Bn(m)
    x = v_tilde + v_hat

    # Calculate y
    y = [Bn(1)] * (sqrt_L + 1)
    return x, y


"""
---------------------------------------------------------------------------------------------------------------------
Database Related Functions

These functions are used to save and load data from SQLite databases.
---------------------------------------------------------------------------------------------------------------------
"""
def save_crs(crs,filename='crs/crs.db'):
    """
    :param crs: The CRS object to be saved in the database.
    :param filename: The filename of the database to save the CRS object to. Default is 'crs/crs.db'.
    :return: None

    The save_crs method saves a CRS object to a SQLite database. The CRS object contains various attributes including L, n1, n2, gamma, g1, g2, and t. These attributes
    * are stored in different tables in the database.

    Example usage:
        crs = CRS(...) # create a CRS object
        save_crs(crs) # save the CRS object in the default database 'crs/crs.db'
        save_crs(crs, 'custom.db') # save the CRS object in a custom database 'custom.db'
    """
    conn = sqlite3.connect(filename)
    cur = conn.cursor()
    #save crs L,n1,n2
    cur.execute('CREATE TABLE IF NOT EXISTS crs_numbers (L INTEGER, n1 INTEGER, n2 INTEGER)')
    cur.execute('CREATE TABLE IF NOT EXISTS crs_g1 (id INTEGER PRIMARY KEY, value BLOB)')
    cur.execute('CREATE TABLE IF NOT EXISTS crs_g2 (id INTEGER PRIMARY KEY, value BLOB)')
    cur.execute('CREATE TABLE IF NOT EXISTS crs_gamma (id INTEGER PRIMARY KEY, value BLOB)')
    cur.execute('CREATE TABLE IF NOT EXISTS crs_t (id INTEGER PRIMARY KEY, value BLOB)')


    cur.execute("INSERT INTO crs_numbers (L, n1, n2) VALUES (?, ?, ?)", (crs.L, crs.n1, crs.n2))
    for i in range(crs.L):
        cur.execute("INSERT INTO crs_gamma (id, value) VALUES (?, ?)", (i, crs.gamma[i].to_binary()))
    for i in range(crs.n2):
        cur.execute("INSERT INTO crs_t (id, value) VALUES (?, ?)", (i, crs.t[i].to_binary()))

    cur.execute("INSERT INTO crs_g1 (id, value) VALUES (?, ?)", (0, crs.g1.to_binary()))
    cur.execute("INSERT INTO crs_g2 (id, value) VALUES (?, ?)", (0, crs.g2.to_binary()))


    conn.commit()

def load_crs(filename):
    """
    Load CRS from a SQLite database file.

    :param filename: The path of the SQLite database file.
    :return: The loaded CRS object.
    """
    conn = sqlite3.connect(filename)
    cur = conn.cursor()
    cur.execute("SELECT * FROM crs_numbers")
    row = cur.fetchone()
    L = row[0]
    n1 = row[1]
    n2 = row[2]
    cur.execute("SELECT * FROM crs_g1")
    row = cur.fetchone()
    g1 = G1Element.from_binary(row[1])
    cur.execute("SELECT * FROM crs_g2")
    row = cur.fetchone()
    g2 = G2Element.from_binary(row[1])
    cur.execute("SELECT * FROM crs_gamma")
    row = cur.fetchone()
    gamma = [G2Element.from_binary(row[1]) for _ in range(L)]
    cur.execute("SELECT * FROM crs_t")
    row = cur.fetchone()
    t = [G2Element.from_binary(row[1]) for _ in range(n2)]
    crs = CRS(L,n1,n2,None,False,False,False)
    crs.set_gamma(gamma)
    crs.set_t(t)
    crs.g1 = g1
    crs.g2 = g2
    return crs

def save_aux(aux,filename='aux.db'):
    """
    Save the auxiliary data to a SQLite database file.

    :param aux: The auxiliary data to be saved.
    :param filename: (Optional) The filename of the SQLite database file. Default is 'aux.db'.
    :return: None

    """
    conn = sqlite3.connect(filename)
    cur = conn.cursor()
    cur.execute('CREATE TABLE IF NOT EXISTS aux_h1 (id INTEGER PRIMARY KEY, value BLOB)')
    cur.execute('CREATE TABLE IF NOT EXISTS aux_h2 (id INTEGER PRIMARY KEY, value BLOB)')
    for i in range(len(aux)):
        cur.execute("INSERT INTO aux_h1 (id, value) VALUES (?, ?)", (i, aux[i].h1.to_binary()))
        cur.execute("INSERT INTO aux_h2 (id, value) VALUES (?, ?)", (i, aux[i].h2.to_binary()))
    conn.commit()
def load_aux(filename='aux.db'):
    """
    Load the auxiliary data from the given database file.

    :param filename: The filename of the SQLite database to load the auxiliary data from. Default value is 'aux.db'.
    :return: A tuple of two lists, where the first list contains the loaded data from the 'aux_h1' table, and the second list contains the loaded data from the 'aux_h2' table.

    """
    conn = sqlite3.connect(filename)
    cur = conn.cursor()
    cur.execute("SELECT * FROM aux_h1")
    row = cur.fetchone()
    h1 = [G2Element.from_binary(row[1]) for _ in range(len(row))]
    cur.execute("SELECT * FROM aux_h2")
    row = cur.fetchone()
    h2 = [G2Element.from_binary(row[1]) for _ in range(len(row))]
    return h1,h2

def load_hsk(index,filename='aux.db'):
    """
    Load a HelperSecretKeyRPLBE object from the database.

    :param index: The index of the object to load.
    :param filename: The filename of the database (default is 'aux.db').
    :return: A HelperSecretKeyRPLBE object loaded from the database.
    """
    conn = sqlite3.connect(filename)
    cur = conn.cursor()
    cur.execute("SELECT * FROM aux_h1 WHERE id = ?", (index,))
    row = cur.fetchone()
    h1 = G2Element.from_binary(row[1])
    cur.execute("SELECT * FROM aux_h2 WHERE id = ?", (index,))
    row = cur.fetchone()
    h2 = G2Element.from_binary(row[1])
    return HelperSecretKeyRPLBE(h1,h2)


def save_sks(sks,filename='sks.db'):
    """
    This method saves a list of sks objects to a SQLite database.

    sks = [sks_obj1, sks_obj2, sks_obj3]
    save_sks(sks, 'my_sks.db')
    """
    conn = sqlite3.connect(filename)
    cur = conn.cursor()
    cur.execute('CREATE TABLE IF NOT EXISTS sks (id INTEGER PRIMARY KEY, value BLOB)')
    for i in range(len(sks)):
        cur.execute("INSERT INTO sks (id, value) VALUES (?, ?)", (i, sks[i].to_binary()))
    conn.commit()

def load_sk(index,filename='sks.db'):
    """
    Load the secret key (SK) from the specified index in the sks.db file.

    :param index: The index of the secret key to load.
    :param filename: (Optional) The name of the sks.db file to use. Defaults to 'sks.db'.
    :return: The loaded secret key.
    """
    conn = sqlite3.connect(filename)
    cur = conn.cursor()
    cur.execute("SELECT * FROM sks WHERE id = ?", (index,))
    row = cur.fetchone()
    return G2Element.from_binary(row[1])
def save_mpk(mpk,filename='mpk.msgpack'):
    """
    Save the given mpk object as a binary MessagePack file.

    :param mpk: The mpk object to be saved.
    :type mpk: object

    :param filename: The name of the output file. Defaults to 'mpk.msgpack'.
    :type filename: str

    :return: None
    """
    data = {
        's': [s.to_binary() for s in mpk.s],
        'w': mpk.w.to_binary(),
        't': [t.to_binary() for t in mpk.t]
    }
    with open(filename, 'wb') as outfile:
        packed = msgpack.packb(data)
        outfile.write(packed)
def load_mpk(filename='mpk.msgpack'):
    """
    Load MasterPublicKey from a msgpack file.

    :param filename: The name of the msgpack file to load. Default is 'mpk.msgpack'.
    :return: A MasterPublicKey object.
    """
    with open(filename, 'rb') as data_file:
        data = msgpack.unpackb(data_file.read())
    s = [G1Element.from_binary(s) for s in data['s']]
    w = G1Element.from_binary(data['w'])
    t = [G2Element.from_binary(t) for t in data['t']]
    return MasterPublicKey(s,w,t)

def save_ct(c,filename='ct.msgpack'):
    """
    :param c: The ciphertext to be saved.
    :param filename: The name of the file to save the ciphertext to. Default is 'ct.msgpack'.
    :return: None

    Saves the given ciphertext to a file in msgpack format. The ciphertext is serialized in the following format:
    - c[0] is converted to binary and stored as 'c1'.
    - c[1] is converted to binary and stored as 'c2'.
    - Each element in c[2] is converted to binary and stored as a list in 'c3'.
    - Each element in c[3] is converted to binary and stored as a list in 'c4'.

    Example Usage:
    save_ct(ciphertext, 'encrypted.ct')
    """
    data = {
        'c1': c[0].to_binary(),
        'c2': c[1].to_binary(),
        'c3': [c[2][i][0].to_binary() for i in range(len(c[2]))],
        'c4': [c[3][i][0].to_binary() for i in range(len(c[3]))],
    }
    with open(filename, 'wb') as outfile:
        packed = msgpack.packb(data)
        outfile.write(packed)

"""
---------------------------------------------------------------------------------------------------------------------
Benchmarking functions
---------------------------------------------------------------------------------------------------------------------
"""

def pairings_times():
    """
    Computes the average time for exponentiation in G1, G2, GT, and pairing
    over 100 iterations.

    :return: None
    """
    # Initialize the elements
    g1_elem = G1.generator()
    g2_elem = G2.generator()
    gt_elem = GT.generator()

    # Variables to store the time
    exp_g1_time = 0
    exp_g2_time = 0
    exp_gt_time = 0
    pairing_time = 0

    # Repeat the operations 100 times
    for _ in range(100):
        # Exponentiation in G1
        start = time.perf_counter()
        _ = g1_elem ** 2
        end = time.perf_counter()
        exp_g1_time += (end - start) * 1000

        # Exponentiation in G2
        start = time.perf_counter()
        _ = g2_elem ** 2
        end = time.perf_counter()
        exp_g2_time += (end - start) * 1000

        # Exponentiation in GT
        start = time.perf_counter()
        _ = gt_elem ** 2
        end = time.perf_counter()
        exp_gt_time += (end - start) * 1000

        # Pairing
        start = time.perf_counter()
        _ = g1_elem.pair(g2_elem)
        end = time.perf_counter()
        pairing_time += (end - start) * 1000  # Convert to milliseconds

    # Compute the averages
    exp_g1_time_avg = exp_g1_time / 100
    exp_g2_time_avg = exp_g2_time / 100
    exp_gt_time_avg = exp_gt_time / 100
    pairing_time_avg = pairing_time / 100

    # Print the averages
    print("Average time for exponentiation in G1: ", round(exp_g1_time_avg, 4), "ms")
    print("Average time for exponentiation in G2: ", round(exp_g2_time_avg, 4), "ms")
    print("Average time for exponentiation in GT: ", round(exp_gt_time_avg, 4), "ms")
    print("Average time for pairing: ", round(pairing_time_avg, 4), "ms")



def convert_and_save():
    """
    Converts the time values in a CSV file to milliseconds and saves the updated data to a new CSV file.

    Usage:
        convert_and_save()
    """
    # Read the data from the csv file
    df = pd.read_csv('times.csv')

    # Convert time to milliseconds and round to 4 decimal places
    for col in df.columns:
        if col != 'L':
            df[col] = round(df[col] * 1000, 0)

    # Save the dataframe to new csv file
    df.to_csv('times_rounded.csv', index=False)


def plot_times():
  """
  Plot the time for each operation based on the data in 'times_rounded.csv'.

  :return: None
  """
  # Read the data from the new rounded csv file
  df = pd.read_csv('times_rounded.csv')

  # Set up the figure and axes
  fig, ax = plt.subplots()

  # Convert 'L' column to list for x axis
  x = df['L'].tolist()

  # Each line represents the avg time of each operation
  operations = ['Setup', 'KGen', 'Aggr', 'Enc', 'Dec']
  colors = ['b', 'g', 'r', 'c', 'm']

  for operation, color in zip(operations, colors):
      y = df[operation].tolist()
      ax.plot(x, y, color + 'o-', label=operation)  # Added 'o-' to the color code

  # Set the title, labels, and legend
  ax.set_xlabel('L')
  ax.set_ylabel('Time (ms, log scale)')
  ax.set_yscale('log')
  ax.set_xscale('log')
  ax.set_xticks([16, 64, 256, 1024])  # Set the scale values on X-axis
  ax.set_xticklabels(['16', '64', '256', '1024'])  # Set the labels for the scales on X-axis
  ax.legend(loc='upper left')

  # Save and show the plot
  plt.tight_layout()
  plt.savefig('plots/operations_time.png')
  plt.show()


"""
---------------------------------------------------------------------------------------------------------------------
Some functions that might be used later
---------------------------------------------------------------------------------------------------------------------
"""
# def test_relation(L, m):
#     i = 0
#     l=0
#     x, y = Z(m, L)
#     M_matrix = M(l, L)
#
#     # calculate xMy
#     temp = utils.vec_mat_mul(x, M_matrix)
#     xMy = utils.inner_prod(temp, y)
#
#     assert F(l, i, m) == xMy, f"Test failed for l={l}, i={i}, m={m}: F(l,i,m)={F(l, i, m)} but xMy={xMy}"


# def M(k, L):
#     # Size of the matrix
#     size = int(L**0.5)
#
#     # Getting k1 and k2
#     k1, k2 = encode_number(k, L)
#     # Initialize matrix with Bn(0)
#     matrix = [[Bn(0) for _ in range(size+1)] for _ in range(2*size)]
#
#     # Set specific position to Bn(1)
#     # Assume matrix indices start from 0
#     matrix[k1+size][k2+1] = Bn(1)
#     matrix[k1][0] = Bn(1)
#     return matrix