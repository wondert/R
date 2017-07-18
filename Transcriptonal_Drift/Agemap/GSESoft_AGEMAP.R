# Download SOFT formatted family file
# In bash shell filter out lines using sed that start with ^ or #
# sed '/^#/ d' GSE7829_family.soft.txt > GSE7829_family.soft2.txt
# sed '/^\^/ d' GSE7829_family.soft2.txt > GSE7829_family.soft.cleaned.txt
filename_mat <- "GSE7829_family.soft.cleaned.txt"
matrix_table <- read.table(filename_mat, header = TRUE, row.names = 1, sep = "\t", blank.lines.skip = TRUE, comment.char = "!", fill = TRUE, skip = 0, nrows = -1)
N_rows <- length(matrix_table[, 1])
N_cols <- length(matrix_table[1, ])
cat("The data set", filename_mat, "has:", "\n", sep = " ")
cat("N_rows=", N_rows, "rows and \n", sep = " ")
cat("N_cols=", N_cols, "columns\n", sep = " ")
outputfilename <- paste(filename_mat, "-exprs.tab", sep = "")