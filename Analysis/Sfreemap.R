install.packages('devtools')
require(devtools)
install_github('dpasqualin/sfreemap')




install.packages("installr")

library(installr)

updateR()
packageStatus()
update.packages(checkBuilt=TRUE)
version
installed.packages()


require(sfreemap)
help(package = "sfreemap", help_type = "html")

mod.tip.dat <- array(dim = c(100,2))
mod.tip.dat[,1] <- tip.dat[,2]
mod.tip.dat[,2] <- tip.dat[,3]
colnames(mod.tip.dat) <- c('char1', 'char2')
rownames(mod.tip.dat) <- tip.dat[,1]
test <- sfreemap(tree = tree, tip_states = mod.tip.dat, model = 'ARD')
