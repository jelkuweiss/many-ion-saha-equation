import functions as saha

M = saha.augmented_matrix_builder(["H", "He"], [1662424176, 719087680], 1.549e+07)
#M = saha.augmented_matrix_builder(["H", "He","O"], [5000, 5000, 5000], 25000)
print(M)
s = saha.matrix_solver_for_ne(M, 3100599537)


#print(Saha.augmented_matrix_builder(["H", "He"], [5000, 5000], 25000))
#print(Saha.augmented_matrix_builder(["H", "C"], [5000, 5000], 25000))
#print(saha.augmented_matrix_builder(["H", "He", "C", "N", "O"], [5000, 5000, 5000, 200, 200], 25000))