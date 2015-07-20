source('~/Development/matchtm/match-tm.r')

# MATCH WEBSERVER V$HNF1_01 M00132
# TEST SEQUENCE: CTCCATGGGAGTTTCTGAAGAACCTTAGTTAATAATTTTCACAGCTGTGCAC
# MATCH: agTTAATaattttca core.match 1.00, matrix.match 0.948
li = 'PO      A      C      G      T  X
01      6      0     16      3      G
02      2      0     24      0      G
03      1      0      0     25      T
04      2      0      0     24      T
05     25      0      0      1      A
06     20      2      1      3      A
07      2      0      0     24      T
08      8      3      8      7      N
09     16      0      1      9      W
10      3      1      1     21      T
11      1      4      0     21      T
12     16      1      6      3      A
13      9     11      3      3      M
14      7     14      0      5      C
15     11      7      4      3      N'

li = strsplit(li, '\n')[[1]]
li = gsub('\\s', '\t', li)
mat = read.table(textConnection(li), header=T)
mat = mat[,-c(1, ncol(mat))]
mat = as.matrix(t(mat))

pwm = makePWM(mat)

s = matchScore(toupper('agTTAATaattttca'), pwm, both.strands = F)