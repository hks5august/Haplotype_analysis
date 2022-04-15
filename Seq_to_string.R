#Convert sequnces into string 
x <- "ATCGCTGCAATTT"
n <- 1
substring(x, seq(1, nchar(x), n), seq(n, nchar(x) + n - 1, n)) 