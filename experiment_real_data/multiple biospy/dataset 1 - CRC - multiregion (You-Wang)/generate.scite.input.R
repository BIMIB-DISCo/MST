library(TRONCO)
options(scipen = 999)
script = ''
epos_level = 0.05
eneg_level = 0.05
steps = 900000

if (! dir.exists('scite_input')) {
    dir.create('scite_input')
}

if (! dir.exists('scite_input/datasets/')) {
    dir.create('scite_input/datasets/')
}

load('data.RData')
complete = data
geno = as.genotypes(complete)
geno = keysToNames(complete, geno)
geno = geno[-6,]
seed_level = round(runif(1) * 10000, 0)

scite.genotypes = t(geno)
nsample = ncol(scite.genotypes)
nmuts = nrow(scite.genotypes)
filename = paste0('datasets/complete.csv')
write.table(scite.genotypes, file = paste0('scite_input/', filename), quote = FALSE, sep = ' ', col.names = FALSE, row.names = FALSE)
script = paste0(script, '../scite -i ', filename,
    ' -n ', nmuts,
    ' -m ', nsample,
    ' -r 1 -l ', steps,
    ' -fd ', epos_level,
    ' -ad ', eneg_level,
    ' -seed ', 195,
    ' -max_treelist_size 1', ' \n')

load('data.merged.RData')
ind = data.merged
geno = as.genotypes(ind)
geno = keysToNames(ind, geno)
geno = geno[-6,]
seed_level = round(runif(1) * 10000, 0)

scite.genotypes = t(geno)
nsample = ncol(scite.genotypes)
nmuts = nrow(scite.genotypes)
filename = paste0('datasets/ind.csv')
write.table(scite.genotypes, file = paste0('scite_input/', filename), quote = FALSE, sep = ' ', col.names = FALSE, row.names = FALSE)
script = paste0(script, '../scite -i ', filename,
    ' -n ', nmuts,
    ' -m ', nsample,
    ' -r 1 -l ', steps,
    ' -fd ', epos_level,
    ' -ad ', eneg_level,
    ' -seed ', 1184,
    ' -max_treelist_size 1', ' \n')


script.filename = paste0('scite_input/scite.script.sh')
write(script, script.filename)


script = ''

load('data.merged.RData')
ind = data.merged
geno = as.genotypes(ind)
geno = keysToNames(ind, geno)
geno = geno[-6,]
steps = 100000
epos_level = 0.00000000001
eneg_level = 0.00000000001
set.seed(12345)

for (i in 1:100) {
    sampled.data =
        geno[sample(1:nrow(geno),
            size = nrow(geno),
            replace = TRUE), ]

    scite.genotypes = t(sampled.data)
    nsample = ncol(scite.genotypes)
    nmuts = nrow(scite.genotypes)
    filename = paste0('datasets/bootstrap_', i, '.csv')
    write.table(scite.genotypes, file = paste0('scite_input/', filename), quote = FALSE, sep = ' ', col.names = FALSE, row.names = FALSE)
    script = paste0(script, '../scite -i ', filename,
        ' -n ', nmuts,
        ' -m ', nsample,
        ' -r 1 -l ', steps,
        ' -fd ', epos_level,
        ' -ad ', eneg_level,
        ' -seed ', round(runif(1) * 10000, 0),
        ' -max_treelist_size 1', ' \n')
}

script.filename = paste0('scite_input/scite.script.bootstrap.sh')
write(script, script.filename)