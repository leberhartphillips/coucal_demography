# matrix_structure() creates a matrix from a set of vital rates

matrix_structure <- expression(
  # top row of matrix
  0, NA, 0, NA,
  
  # second row of matrix
  (F_Nestling_survival * F_Groundling_survival * F_Fledgling_survival), F_Adult_survival, 0, 0,
  
  # third row of matrix
  0, NA, 0, NA,
  
  # fourth row of matrix
  0, 0, (M_Nestling_survival * M_Groundling_survival * M_Fledgling_survival), M_Adult_survival
  
)