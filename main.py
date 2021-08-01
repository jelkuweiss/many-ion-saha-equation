import functions as Saha

#M = Saha.augmented_matrix_builder(["H", "He"], [5000, 5000], 25000)
#print(M)
#s = Saha.matrix_solver_for_ne(M,10**40)


print(Saha.augmented_matrix_builder(["H", "He"], [5000, 5000], 25000))
print(Saha.augmented_matrix_builder(["H", "C"], [5000, 5000], 25000))
print(Saha.augmented_matrix_builder(["H", "He", "C", "N", "O"], [5000, 5000, 5000, 200, 200], 25000))