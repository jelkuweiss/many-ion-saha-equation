import functions as Saha


print(Saha.augmented_matrix_builder(["H", "He"], [5000, 5000], 25000))
print(Saha.augmented_matrix_builder(["H", "C"], [5000, 5000], 25000))


print(Saha.augmented_matrix_builder(["H", "He", "C", "N", "O"], [5000, 5000, 5000, 200, 200], 25000))