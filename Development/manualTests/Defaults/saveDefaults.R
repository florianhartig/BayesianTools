samp = getPossibleSamplerTypes()

defaults = list()

for(i in 1:(length(samp$BTname))){
  defaults[[i]] = applySettingsDefault(sampler = samp$BTname[i]) 
}

save(defaults, file = "./Development/manualTests/Defaults/defaults.Rdata")



