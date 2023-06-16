from sys import maxsize
from itertools import permutations

import numpy as np

# implementation of traveling Salesman Problem
# https://www.geeksforgeeks.org/traveling-salesman-problem-tsp-implementation/
# https://stackoverflow.com/questions/66084154/python-dynamic-programing-tsp-print-route
def travellingSalesmanProblem(graph, s):
  # keep all vertex other than the starting point
  vertex = [] 

  #traverse the diagram 
  for i in range(len(graph)): 
    if i != s: 
      vertex.append(i) 

  # keep minimum weight
  min_path = maxsize 
  next_permutation = permutations(vertex)
  best_path = []
  
  for i in next_permutation:
    print(i)
    # store current Path weight(cost) 
    current_pathweight = 0
    
    # compute current path weight 
    k = s 
    
    for j in i:
      current_pathweight += graph[k][j] 
      k = j 
    current_pathweight += graph[k][s] 
    
    # update minimum 
    if current_pathweight < min_path:
      min_path = current_pathweight
      best_path = [s]
      best_path.extend(list(i))
      best_path.append(s)
    
  return min_path, best_path

# https://www.geeksforgeeks.org/program-to-check-if-a-matrix-is-symmetric
# Returns true if mat is symmetric, else false
def isSymmetric(mat):
  N = len(mat)
  def transpose(mat, tr, N):
    for i in range(N):
        for j in range(N):
            tr[i][j] = mat[j][i]
     
  tr = [ [0 for j in range(len(mat[0])) ] for i in range(len(mat)) ]
  transpose(mat, tr, N)
  for i in range(N):
    for j in range(N):
      if (mat[i][j] != tr[i][j]):
        return False
  return True

# Grantham's Distance Matrix
# https://en.wikipedia.org/wiki/Amino_acid_replacement
original_grantham = [[110, 145, 74, 58, 99, 124, 56, 142, 155, 144, 112, 89, 68, 46, 121, 65, 80, 135, 177],
[102, 103, 71, 112, 96, 125, 97, 97, 77, 180, 29, 43, 86, 26, 96, 54, 91, 101],
[98, 92, 96, 32, 138, 5, 22, 36, 198, 99, 113, 153, 107, 172, 138, 15, 61],
[38, 27, 68, 42, 95, 114, 110, 169, 77, 76, 91, 103, 108, 93, 87, 147],
[58, 69, 59, 89, 103, 92, 149, 47, 42, 65, 78, 85, 65, 81, 128],
[64, 60, 94, 113, 112, 195, 86, 91, 111, 106, 126, 107, 84, 148],
[109, 29, 50, 55, 192, 84, 96, 133, 97, 152, 121, 21, 88],
[135, 153, 147, 159, 98, 87, 80, 127, 94, 98, 127, 184],
[21, 33, 198, 94, 109, 149, 102, 168, 134, 10, 61],
[22, 205, 100, 116, 158, 102, 177, 140, 28, 40],
[194, 83, 99, 143, 85, 160, 122, 36, 37],
[174, 154, 139, 202, 154, 170, 196, 215],
[24, 68, 32, 81, 40, 87, 115],
[46, 53, 61, 29, 101, 130],
[94, 23, 42, 142, 174],
[101, 56, 95, 110],
[45, 160, 181],
[126, 152],
[67]
]

#augment with zeros to maintain original shape
for idx, r in enumerate(original_grantham):
   zeros = [0] * idx
   original_grantham[idx] = zeros + r

grantham_adjacency_matrix = [ [0]*20 for i in range(20)] # 20 amino acids

# Transform the matrix to be an adjacency matrix with the below shape so we can treat it as a graph for the traveling salesman problem
    # Ser Arg Leu Pro Thr Ala Val Gly Ile Phe Tyr Cys His Gln Asn Lys Asp Glu Met Trp
# Ser 
# Arg
# Leu
# Pro
# Thr
# Ala
# Val
# Gly
# Ile
# Phe
# Tyr
# Cys
# His
# Gln
# Asn
# Lys
# Asp
# Glu
# Met
# Trp
for r_idx, r in enumerate(original_grantham):
  i = 0
  while i < r_idx:
    grantham_adjacency_matrix[r_idx][i] = original_grantham[i][r_idx - 1]
    i += 1
  c_idx = i
  while c_idx < len(r):
    grantham_adjacency_matrix[r_idx][c_idx + 1] = original_grantham[r_idx][c_idx]
    c_idx += 1
  grantham_adjacency_matrix[-1][r_idx] = original_grantham[r_idx][-1] # fill in the last Trp row that's not a row in the original matrix (only a col)

for r in grantham_adjacency_matrix:
  print(r)
print(isSymmetric(grantham_adjacency_matrix)) # the adj matrix should be symmetric

amidic_coords = np.ix_([13, 14],[13, 14])
aliphatic_coords = np.ix_([2, 3, 5, 6, 7, 8],[2, 3, 5, 6, 7, 8]) 
aromatic_coords = np.ix_([9, 10, 19],[9, 10, 19])
acidic_coords = np.ix_([16, 17],[16, 17]) 
basic_coords = np.ix_([1, 12, 15],[1, 12, 15]) 
hydroxylic_coords = np.ix_([0, 4], [0, 4])
sulphur_containing_coords = np.ix_([16, 17],[16, 17]) 
grantham_adjacency_matrix_np = np.asmatrix(grantham_adjacency_matrix)
print(grantham_adjacency_matrix_np[amidic_coords])
print(grantham_adjacency_matrix_np[aliphatic_coords])
print(grantham_adjacency_matrix_np[aromatic_coords])
print(grantham_adjacency_matrix_np[acidic_coords])
print(grantham_adjacency_matrix_np[basic_coords])
print(grantham_adjacency_matrix_np[hydroxylic_coords])
print(grantham_adjacency_matrix_np[sulphur_containing_coords])
print(travellingSalesmanProblem(grantham_adjacency_matrix_np[aliphatic_coords].tolist(), 0))

amino_acids_row_index = {
  "Ser": 0,
  "Arg": 1,
  "Leu": 2,
  "Pro": 3,
  "Thr": 4,
  "Ala": 5,
  "Val": 6,
  "Gly": 7,
  "Ile": 8,
  "Phe": 9,
  "Tyr": 10,
  "Cys": 11,
  "His": 12,
  "Gln": 13,
  "Asn": 14,
  "Lys": 15,
  "Asp": 16,
  "Glu": 17,
  "Met": 18,
  "Trp": 19
}



