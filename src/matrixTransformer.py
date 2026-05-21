#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Authors: Chris (chris@nohup.cc) & Young (young@nohup.cc)

import sys
import numpy as np
import math
import re
import fileinput

def average(matrixSum, seqLen):
    # average the summary of rows
    matrix_array = np.array(matrixSum)
    matrix_array = np.divide(matrix_array, seqLen)
    matrix_array_shp = np.shape(matrix_array)
    matrix_average = [np.reshape(matrix_array, (matrix_array_shp[0] * matrix_array_shp[1], ))]
    return matrix_average

def normalizePSSM(PSSM):
    PSSM = PSSM[:, 1:21]
    PSSM = PSSM.astype(float)
    seq_cn = np.shape(PSSM)[0]
    PSSM_norm = np.zeros((seq_cn, 20))  # Create a zero matrix for normalization
    mean_matrix = np.mean(PSSM, axis=1)
    std_matrix = np.std(PSSM, axis=1)

    for i in range(seq_cn):
        for j in range(20):
            if std_matrix[i] == 0.0:
                PSSM_norm[i][j] = PSSM[i][j] - mean_matrix[i]
            else:
                PSSM_norm[i][j] = (PSSM[i][j] - mean_matrix[i]) / std_matrix[i]
    return PSSM_norm

def window(PSSM, w_smth, w_slide):
    w_smth = int(w_smth)
    w_slide = int(w_slide)
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"

    PSSM = PSSM[:, 1:21].astype(float)
    seq_cn = np.shape(PSSM)[0]

    # Original PSSM
    PSSM_smth = np.zeros((seq_cn, 20))
    PSSM_orig = np.array(PSSM)

    # Section for PSSM_smth features
    PSSM_smth_full = pssm_smth(PSSM_orig, PSSM_smth, w_smth, seq_cn)
    PSSM_smth_final = np.zeros((w_slide, 20))

    for i in range(w_slide):
        PSSM_smth_final[i] = PSSM_smth_full[i]
    matrix_final = average(PSSM_smth_final, 1.0)
    return matrix_final

def pssm_smth(PSSM_orig, PSSM_smth, w_smth, l):
    for i in range(l):
        if i < (w_smth - 1) // 2:
            for j in range(i + (w_smth - 1) // 2 + 1):
                PSSM_smth[i] += PSSM_orig[j]
        elif i >= (l - (w_smth - 1) // 2):
            for j in range(i - (w_smth - 1) // 2, l):
                PSSM_smth[i] += PSSM_orig[j]
        else:
            for j in range(i - (w_smth - 1) // 2, i + (w_smth - 1) // 2 + 1):
                PSSM_smth[i] += PSSM_orig[j]
    return PSSM_smth

def handleRows(PSSM, SWITCH, COUNT):
    '''
    if SWITCH=0, we filter no element.
    if SWITCH=1, we filter all the negative elements.
    if SWITCH=2, we filter all the negative and positive elements greater than expected.
    '''
    '''
    if COUNT=20, we generate a 20-dimension vector.
    if COUNT=400, we generate a 400-dimension vector.
    '''
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"
    matrix_final = np.zeros((COUNT // 20, 20))
    seq_cn = 0

    PSSM_shape = np.shape(PSSM)
    for i in range(PSSM_shape[0]):
        seq_cn += 1
        str_vec = PSSM[i]
        str_vec_positive = np.array(list(map(int, str_vec[1:21])))  # Convert to list then array
        if SWITCH == 1:
            str_vec_positive[str_vec_positive < 0] = 0
        elif SWITCH == 2:
            str_vec_positive[str_vec_positive < 0] = 0
            str_vec_positive[str_vec_positive > 7] = 0

        if COUNT == 20:
            matrix_final[0] = list(map(sum, zip(str_vec_positive, matrix_final[0])))
        elif COUNT == 400:
            matrix_final[Amino_vec.index(str_vec[0])] = list(map(sum, zip(str_vec_positive, matrix_final[Amino_vec.index(str_vec[0])])))

    return matrix_final

def preHandleColumns(PSSM, STEP, PART, ID):
    if PART == 0:
        PSSM = PSSM[:, 1:21]
    elif PART == 1:
        PSSM = PSSM[:, 21:]
    PSSM = PSSM.astype(float)
    matrix_final = np.zeros((20, 20))
    seq_cn = np.shape(PSSM)[0]

    if ID == 0:
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - STEP):
                    matrix_final[i][j] += (PSSM[k][i] * PSSM[k + STEP][j])

    elif ID == 1:
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - STEP):
                    matrix_final[i][j] += ((PSSM[k][i] - PSSM[k + STEP][j]) * (PSSM[k][i] - PSSM[k + STEP][j]) / 4.0)

    return matrix_final

def handleTriColumns(PSSM):
    matrix_final = np.zeros((20, 20, 20))
    PSSM = PSSM[:, 21:].astype(float)
    seq_cn = np.shape(PSSM)[0]
    for m in range(20):
        for n in range(20):
            for r in range(20):
                for i in range(seq_cn - 2):
                    matrix_final[m][n][r] += (PSSM[i][m] * PSSM[i + 1][n] * PSSM[i + 2][r])
    matrix_final = np.divide(matrix_final, 1000000.0)
    matrix_final_shape = np.shape(matrix_final)
    matrix_result = [np.reshape(matrix_final, (matrix_final_shape[0] * matrix_final_shape[1] * matrix_final_shape[2], ))]
    return matrix_result

def handleMixed(PSSM, ALPHA):
    row1 = np.zeros(20)
    row2 = np.zeros(20)
    matrix_final = np.zeros((1, 40))

    PSSM_norm = normalizePSSM(PSSM)
    seq_cn = np.shape(PSSM)[0]
    for i in range(seq_cn):
        row1 += PSSM_norm[i]
    row1 = np.divide(row1, seq_cn)

    for j in range(20):
        for i in range(seq_cn - ALPHA):
            row2[j] += (PSSM_norm[i][j] - PSSM_norm[i + ALPHA][j]) ** 2
    row2 = np.divide(row2, seq_cn - ALPHA)

    row = np.hstack((row1, row2))
    matrix_final[0] = row
    return matrix_final

def handleMixed2(PSSM, ALPHA):
    row1 = np.zeros(40)
    row2 = np.zeros((20, 2 * ALPHA))
    matrix_final = np.zeros((1, 40 + 40 * ALPHA))

    PSSM_norm = normalizePSSM(PSSM)
    seq_cn = np.shape(PSSM)[0]
    for j in range(20):
        positive_count_1 = 0
        negative_count_1 = 0
        for i in range(seq_cn):
            if PSSM_norm[i][j] >= 0:
                positive_count_1 += 1
                row1[2 * j] += PSSM_norm[i][j]
            elif PSSM_norm[i][j] < 0:
                negative_count_1 += 1
                row1[2 * j + 1] += PSSM_norm[i][j]
        if positive_count_1 > 0:
            row1[2 * j] /= positive_count_1
        if negative_count_1 > 0:
            row1[2 * j + 1] /= negative_count_1

    for i in range(seq_cn - ALPHA):
        for j in range(20):
            row2[j][0] += (PSSM_norm[i][j] - PSSM_norm[i + ALPHA][j]) ** 2
            row2[j][1] += (PSSM_norm[i + ALPHA][j] - PSSM_norm[i][j]) ** 2

    for j in range(20):
        row2[j][0] /= (seq_cn - ALPHA)
        row2[j][1] /= (seq_cn - ALPHA)

    for j in range(20):
        row1 = np.hstack((row1, row2[j]))

    matrix_final[0] = row1
    return matrix_final

def pssm(file):
    matrix_final = []
    for line in fileinput.input(file):
        if not line.startswith("#"):
            line = line.strip().split("\t")
            matrix_final.append(line)
    matrix_final = np.array(matrix_final)
    return matrix_final

if __name__ == "__main__":
    pssm_file = sys.argv[1]
    PSSM = pssm(pssm_file)
    print(PSSM)

