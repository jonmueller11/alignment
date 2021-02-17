from glob_alg import GlobalAlignment


sco = {'gap': 0, True: 1, False: 0}
alignment = GlobalAlignment('ACGTCATGCCCTGATGGT', 'AGTCCGGGCGCTTTAGCTAGC', sco)
print(alignment)