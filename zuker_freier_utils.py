import numpy

Q_const = 0
M_const = 4.6
Qi_const = 0.4
m_const = 3

def FreierDangleLeft(sequence, x, a, b):
    five_prime = {"A":{"A":-0.3, "C":-0.5, "G":-0.2,"U":-0.3},
                  "C":{"A":-0.3, "C":-0.2, "G":-0.3,"U":-0.2},
                  "G":{"A":-0.4, "C":-0.2, "G":-0.0,"U":-0.2},
                  "U":{"A":-0.2, "C":-0.1, "G":-0.0,"U":-0.2}}
    return five_prime[sequence[x]][sequence[a]]

def FreierDangleRight(sequence, x, a, b):
    three_prime = {"A":{"A":-0.8, "C":-1.7, "G":-1.1,"U":-0.7},
                  "C":{"A":-0.5, "C":-0.8, "G":-0.4,"U":-0.1},
                  "G":{"A":-0.8, "C":-1.7, "G":-1.3,"U":-0.7},
                  "U":{"A":-0.6, "C":-1.2, "G":-0.6,"U":-0.1}}
    return three_prime[sequence[x]][sequence[a]]

def FreierStack(sequence, i, j):
    if sequence[i] == 'G' and sequence[j] == 'C':
        scores = {"A":{"A":-1.1, "C":-1.1, "G":-1.6,"U":-2.1},
                  "C":{"A":-1.3, "C":-0.6, "G":-3.4,"U":-0.5},
                  "G":{"A":-1.3, "C":-3.4, "G":-1.4,"U":-2.3},
                  "U":{"A":-2.3, "C":-0.5, "G":-1.4,"U":-0.7}}
        return scores[sequence[i + 1]][sequence[j - 1]]
    elif sequence[i] == 'C' and sequence[j] == 'G':
        scores = {"A":{"A":-1.9, "C":-1.0, "G":-1.9,"U":-1.7},
                  "C":{"A":-2.0, "C":-1.1, "G":-2.0,"U":-1.5},
                  "G":{"A":-1.9, "C":-2.9, "G":-1.9,"U":-1.9},
                  "U":{"A":-1.8, "C":-0.8, "G":-1.6,"U":-1.2}}
        return scores[sequence[i + 1]][sequence[j - 1]]
    elif sequence[i] == 'A' and sequence[j] == 'U':
        scores = {"A":{"A":-0.8, "C":-0.7, "G":-0.8,"U":-0.9},
                  "C":{"A":-1.0, "C":-0.7, "G":-1.7,"U":-0.8},
                  "G":{"A":-1.0, "C":-2.1, "G":-1.0,"U":-0.9},
                  "U":{"A":-0.9, "C":-0.7, "G":-0.9,"U":-0.8}}
        return scores[sequence[i + 1]][sequence[j - 1]]
    elif sequence[i] == 'U' and sequence[j] == 'A':
        scores = {"A":{"A":-1.0, "C":-0.7, "G":-1.1,"U":-0.9},
                  "C":{"A":-0.8, "C":-0.6, "G":-1.8,"U":-0.6},
                  "G":{"A":-1.1, "C":-1.8, "G":-1.2,"U":-1.0},
                  "U":{"A":-1.1, "C":-0.5, "G":-0.9,"U":-0.5}}
        return scores[sequence[i + 1]][sequence[j - 1]]
    return 0

def FreierIS1(i,j):
    global m_const
    length = j - i - 1
    if length < m_const:
        return float('inf')
    values = {3 : 7.4, 4 : 5.9, 5 : 4.4, 6 : 4.3, 7 : 4.1, 8 : 4.1, 9 : 4.2, 10 : 4.3, 12 : 4.9, 14 : 5.6, 16 : 6.1, 18 : 6.7, 20: 7.1, 25 : 8.1, 30 : 8.9}
    if length in values :
        return values[length]
    else :
        return round(4.0476 * numpy.log(length) - 4.9884, 1)

def FreierIS2(sequence, i, j, k, l) :
    is2= float('inf')
    if k - 1 == i and l + 1 == j:
        is2 = FreierStack(sequence, k, l)
    if k - 1 == i and l + 1 < j or k - 1 > i and l + 1 == j:
        is2 = FreierBulge(l - k - 1)
    else:
        is2 = FreierInternalLoop(l - k - 1)
    return is2

def FreierBulge(length):
    if length < 1:
        return float('inf')
    values = {1 : 3.3, 2 : 5.2, 3 : 6.0, 4 : 6.7, 5 : 7.4, 6 : 8.2, 7 : 9.1, 8 : 10.0, 9 : 3.1, 10 : 3.6, 12 : 4.4, 14 : 5.1, 16 : 5.6, 18 : 6.2, 20: 6.6, 25 : 7.6, 30 : 8.4}
    if length in values :
        return values[length]
    else :
        return round(4.3834 * numpy.log(length) + 0.8927, 1)

def FreierInternalLoop(length):
    if length < 2 :
        return float('inf')
    values = {2 : 0.8, 3 : 1.3, 4 : 1.7, 5 : 2.1, 6 : 2.5, 7 : 2.6, 8 : 2.8, 9 : 3.1, 10 : 3.6, 12 : 4.4, 14 : 5.1, 16 : 5.6, 18 : 6.2, 20: 6.6, 25 : 7.6, 30 : 8.4}
    if (length) in values :
        return values[length]
    else :
        return round(pow(length,2)*-0.0045 + 0.4185*(length)- 0.0003, 1)

def FreierPair(sequence, i, j):
    p = 0
    if sequence[i] == 'A' and sequence[j] == 'U' or sequence[i] == 'U' and sequence[j] == 'A' :
        p = -4.15
    elif sequence[i] == 'G' and sequence[j] == 'C' or sequence[i] == 'C' and sequence[j] == 'G' :
        p = -6.17
    else :
        p = 1
    return p
    