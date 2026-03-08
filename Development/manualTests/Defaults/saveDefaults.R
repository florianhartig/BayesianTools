samp = getPossibleSamplerTypes()

defaults = list()

for(i in 1:(length(samp$BTname))){
  defaults[[i]] = applySettingsDefault(sampler = samp$BTname[i]) 
}

save(defaults, file = "./Development/manualTests/Defaults/defaults.Rdata")


library(BayesianTools)

load(file = "./Development/manualTests/Defaults/defaults.Rdata")

newDefaults = list()
samp = getPossibleSamplerTypes()

for(i in 1:(length(samp$BTname))){
  newDefaults[[i]] = applySettingsDefault(sampler = samp$BTname[i]) 
  
  newDefaults[[i]][[length(newDefaults[[i]])]] = NULL
  defaults[[i]][[length(defaults[[i]])]] = NULL
  
  print(samp$BTname[i])
  print("=======")
  print("not in new")
  print(defaults[[i]][which(! names(defaults[[i]]) %in% names(newDefaults[[i]]))])
  print("not in old")
  print(newDefaults[[i]][which(! names(newDefaults[[i]]) %in% names(defaults[[i]]))])
  print(" ")
}

