# Set your working directory
setwd('~/Development/matchtm/')

# Source functions
source('match-tm.r')


###############################
##        DNA example        ##   
###############################
# Read DNA sequences
dna.seqs = readLines('sample_dna.txt')
# print(dna.seqs)

# Construct your PWM
dna.pwm = makePWM(dna.seqs, log.bg = F, priors = 'ecoli', seq.type = 'auto')

# Score some DNA sequences (can change method to "log")
dna.scores = matchScore(dna.seqs, dna.pwm)$mss
# print(dna.scores)

###############################
##        AA example         ##   
###############################

# Read AA sequences (AKT1 phosphorylation sites)
aa.seqs = readLines('sample_aa.txt')
# print(aa.seqs)

# Construct your PWM
aa.pwm = makePWM(aa.seqs, log.bg = F, priors = 'yeast', seq.type = 'auto')

# Score some DNA sequences (can change method to "log")
aa.scores = matchScore(aa.seqs, aa.pwm)$mss
# print(aa.scores)

###############################
##        PFM example        ##   
###############################

# Read position frequency matrix
pfm = as.matrix( read.table('sample_pfm.txt') )
# print(pfm)

# Construct your PWM
my.pwm = makePWM(pfm, log.bg = F, priors = 'ecoli', seq.type = 'auto')

# Score some DNA sequences (can change method to "log")
pfm.scores = matchScore('TGGCA', my.pwm)$mss
# print(pfm.scores)

