######################################################################
# Various numbers for the paper

######################################################################
# ABSTRACT
######################################################################

# number of mutants 
sum(data100[,2]+data100[,4]) # 986

# number of unique genes
sum(data100[,2]>0 | data100[,4]>0 | c(0,data100[-nrow(data100),4])>0) # 759



# est'd number of non-essential genes
nrow(data80)-mean(finalres$total) # 2375

# est'd percent non-essential genes
(1-mean(finalres$total)/nrow(data80))*100 # 56.5%

# confidence interval for number of non-essential genes
quantile(nrow(data80)-finalres$total, c(0.025, 0.975)) # 2091-2705

# est'd number of essential genes
mean(finalres$total) # 1829

# est'd percent essential genes
mean(finalres$total)/nrow(data80)*100 # 43.5%

# 95% CI for percent essential genes
quantile(finalres$total/nrow(data80)*100, c(0.025, 0.975)) # 35.7 - 50.3%

######################################################################
# RESULTS
######################################################################

# number of intragenic sites
sum(data100[,1]+data100[,3]) # 65,649

# genes with no insertion sites
noTA100 <- info[is.na(match(info[,1],rownames(data100))),]

# number of such
nrow(noTA100) # 16

# range of length of these
range(noTA100[,2]) # 96 - 300

# mean length
mean(noTA100[,2]) # 160

# median length
median(noTA100[,2]) # 132

# distribution of number of hits per gene
table(data100[,2]+data100[,4]+c(0,data100[-nrow(data100),4]))
# 0: 3475    # 4: 8
# 1: 590     # 5: 2
# 2: 129     # 6: 4
# 3: 25      # 7: 1 

# hit more than once = 169
sum(table(data100[,2]+data100[,4]+c(0,data100[-nrow(data100),4]))[-(1:2)])


# number of mutants with insertions in N-terminal 80% of gene
sum(data80[,2]+data80[,4]+c(0,data80[-nrow(data80),4])) # 758

# number of unique genes
sum(data80[,2]>0 | data80[,4]>0 | c(0,data80[-nrow(data80),4])>0) # 593

# number of genes with TA site within N-terminal 80%
nrow(data80) # 4204

# genes with no insertion sites
noTA80 <- info[is.na(match(info[,1],rownames(data80))),]

# number of such
nrow(noTA80) # 46

# range of length of these
range(noTA80[,2]) # 96 - 516

# mean length
mean(noTA80[,2]) # 182

# median length
median(noTA80[,2]) # 153

