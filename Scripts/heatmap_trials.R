set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("gene", seq_len(nr))
colnames(mat) = c(paste0("CC", 1:6), paste0("HKCC", 1:6), paste0("ctrl", 1:6), paste0("CL075", 1:6))

Heatmap(mat)

mat_log <- logCPM %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Symbol") %>% 
  left_join(HKCC_cc_logCPM,.) %>% 
  dplyr::select(-c(logFC, FDR)) %>% 
  tibble::column_to_rownames("Symbol")

mat_with_na = mat_log
na_index = sample(c(TRUE, FALSE), nrow(mat)*ncol(mat), replace = TRUE, prob = c(1, 9))
#mat_with_na[na_index] = NA
mat_with_na[c("ISG15", "DDX60", "SMAD1"),c("CC-1", "CC-2", "CC-3", "CC-4", "HKCC-1", "HKCC-2", "HKCC-3", "HKCC-4")] <- NA
Heatmap(mat_with_na, 
        cluster_rows = T,
        #cluster_columns = T,
        name = "mat", 
        na_col = "green",
        column_title = "a matrix with NA values")
