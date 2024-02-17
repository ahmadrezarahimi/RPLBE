import rplbe_scheme
import utils

import random
import pandas as pd

from petrelic.multiplicative.pairing import GT
import time


def create_dbs():
    """
    Runs the setup, keygen, and aggregate algorithms for RPLBE and saves the CRS, MPK, SKs, and AUXs to files.

    """
    columns = ['L', 'creation_time_crs', 'avg_keygen_time', 'aggregate_time', 'avg_encryption_time', 'avg_decryption_time']
    df = pd.DataFrame(columns=columns)  # Create new dataframe to store all data

    a = [2**4, 2**6,2**8, 2**10]

    for L in a:
        start_crs = time.perf_counter()
        crs = rplbe_scheme.setup(L)
        end_crs = time.perf_counter() - start_crs
        utils.save_crs(crs, "crs/crs_" + str(L) + ".db")

        pks, sks = [], []
        keygen_times = []
        for j in range(crs.L):
            start_kg = time.perf_counter()
            pk, sk = rplbe_scheme.keygen(crs, j)
            end_kg = time.perf_counter()
            pks.append(pk)
            sks.append(sk)
            keygen_times.append(end_kg - start_kg)

        avg_keygen_time = sum(keygen_times) / len(keygen_times)

        start_agg = time.perf_counter()
        mpk, hsk = rplbe_scheme.aggregate(crs, pks)
        end_agg = time.perf_counter() - start_agg

        df.loc[len(df)] = [L, end_crs, avg_keygen_time, end_agg, 'NA', 'NA']  # Append new row to dataframe

        utils.save_sks(sks, "sks/sks_" + str(L) + ".db")
        utils.save_mpk(mpk, "mpks/mpk_" + str(L) + ".msgpack")
        utils.save_aux(hsk, "auxs/aux_" + str(L) + ".db")

    df.to_csv('times.csv', index=False)


def benchmark_enc_dec():
    """
    Benchmark Encryption and Decryption

    This function benchmarks the encryption and decryption process using a predefined set of parameters.
    It performs 100 iterations of encryption and decryption and calculates the average
    * time for each operation. The results are then saved in a CSV file.
    """
    # L's values range
    a = [2 ** 4, 2 ** 6, 2 ** 8, 2 ** 10]
    for L in a:
        gt = GT.generator()
        crs = utils.load_crs(filename="crs/crs_" + str(L) + ".db")
        mpk = utils.load_mpk(filename="mpks/mpk_" + str(L) + ".msgpack")

        enc_times = []
        dec_times = []

        for _ in range(100):  # 100 iterations for encryption and decryption
            message = random.choice([1, 2])  # random message either 1 or 2
            decryption_index = random.randrange(0, L)  # random decryption index

            # Benchmark encryption
            start_enc = time.perf_counter()
            c = rplbe_scheme.encrypt(crs, mpk, message)
            end_enc = time.perf_counter()
            enc_times.append(end_enc - start_enc)  # store encryption time

            # Load keys for decryption
            hsk = utils.load_hsk(decryption_index, filename="auxs/aux_" + str(L) + ".db")
            sk = utils.load_sk(decryption_index, filename="sks/sks_" + str(L) + ".db")

            # Benchmark decryption
            start_dec = time.perf_counter()
            m = rplbe_scheme.decrypt(crs, sk, decryption_index, hsk, c)
            end_dec = time.perf_counter()
            dec_times.append(end_dec - start_dec)  # store decryption time

            assert m == gt ** message

        avg_encryption_time = sum(enc_times) / len(enc_times)  # Calculate average encryption time
        avg_decryption_time = sum(dec_times) / len(dec_times)  # Calculate average decryption time

        # Read existing csv file into pandas DataFrame
        df = pd.read_csv('times.csv')
        # Find row of current L and update encryption and decryption times
        df.loc[df['L'] == L, 'avg_encryption_time'] = avg_encryption_time
        df.loc[df['L'] == L, 'avg_decryption_time'] = avg_decryption_time
        # Save DataFrame back to csv
        df.to_csv('times.csv', index=False)




"""
Benchmarking
"""
# Uncomment if you want to create new CRS, MPK, SKs, and AUXs
# create_dbs()

# Uncomment if you want to benchmark encryption and decryption
# benchmark_enc_dec()
utils.plot_times()