#======================================================================
#======================================================================
#==================  FILE TO RUN CELL CYCLE MODEL =====================
#======================================================================
#======================================================================


require(pracma)
require(reshape2)
require(ggplot2)
require(ggthemes)
require(ggbeeswarm)
require(ggpubr)

#======================================================================
#======================================================================
#=========================  UPLOAD THE DATA ===========================
#======================================================================
#======================================================================

#Text file should be specified according to the following columns:
#Sample Type  Day Gate_names.EdUposBrdUpos  Gate_names.EdUposBrdUneg  Gate_names.EdUnegBrdUpos  Gate_names.EdUnegBrdUneg
#where:
#Sample = name of the sample
#Type = sample type for when it is required to compare several groups
#Day = in case there is time series/kinetics
#Gate_names = name given to the pregatting processes before defining quadrants based on EdU and BrdU signal
#Values should be presented as percentages (between 0 and 100)

upload = function(textfile) {
  treated_data = read.table(textfile, header = TRUE, sep = "\t")
  return(treated_data)
}

convert = function(data) {
  data_frame = melt(data, id.vars = c("Sample", "Day", "Type"), variable.name = "Compartment")
  return(data_frame)
}

#FILENAME
filename = "complete_exp_data.txt"

#TIME ELAPSED BETWEEN THE TWO PULSES
t = 2

data = upload(filename)
data = convert(data)

#======================================================================
#======================================================================
#=========================== TRIM THE DATA ============================
#======================================================================
#======================================================================

compartments = c()
for (compartment in unique(data$Compartment[which(grepl("EdU", data$Compartment))])) {
  pregate_compartment = strsplit(compartment, ".EdU")[[1]][1]
  compartments = c(compartments, pregate_compartment)
}
rm(compartment)
compartments=unique(compartments)

data$value = data$value/100
for (sample in unique(data$Sample)) {
  for (compartment in compartments) {
    oldEdUnegBrdUneg = data$value[which(data$Sample == sample & 
                                          grepl(compartment, data$Compartment) &
                                          grepl("EdUnegBrdUneg",data$Compartment))]
    oldEdUnegBrdUpos = data$value[which(data$Sample == sample & 
                                          grepl(compartment, data$Compartment) &
                                          grepl("EdUnegBrdUpos",data$Compartment))]
    oldEdUposBrdUneg = data$value[which(data$Sample == sample & 
                                          grepl(compartment, data$Compartment) &
                                          grepl("EdUposBrdUneg",data$Compartment))]
    oldEdUposBrdUpos = data$value[which(data$Sample == sample & 
                                          grepl(compartment, data$Compartment) &
                                          grepl("EdUposBrdUpos",data$Compartment))]
    
    newEdUnegBrdUneg = oldEdUnegBrdUneg/(oldEdUnegBrdUneg+oldEdUnegBrdUpos)
    newEdUnegBrdUpos = oldEdUnegBrdUpos/(oldEdUnegBrdUneg+oldEdUnegBrdUpos)
    newEdUposBrdUneg = oldEdUposBrdUneg/(oldEdUposBrdUneg+oldEdUposBrdUpos)
    newEdUposBrdUpos = oldEdUposBrdUpos/(oldEdUposBrdUneg+oldEdUposBrdUpos)
    
    data$value[which(data$Sample == sample & 
                       grepl(compartment, data$Compartment) &
                       grepl("EdUnegBrdUneg",data$Compartment))] = newEdUnegBrdUneg
    
    data$value[which(data$Sample == sample & 
                       grepl(compartment, data$Compartment) &
                       grepl("EdUnegBrdUpos",data$Compartment))] = newEdUnegBrdUpos
    
    data$value[which(data$Sample == sample & 
                       grepl(compartment, data$Compartment) &
                       grepl("EdUposBrdUneg",data$Compartment))] = newEdUposBrdUneg
    
    data$value[which(data$Sample == sample & 
                       grepl(compartment, data$Compartment) &
                       grepl("EdUposBrdUpos",data$Compartment))] = newEdUposBrdUpos
    
  }
}
rm(sample)

rm(oldEdUnegBrdUneg)
rm(oldEdUnegBrdUpos)
rm(oldEdUposBrdUneg)
rm(oldEdUposBrdUpos)

rm(newEdUnegBrdUneg)
rm(newEdUnegBrdUpos)
rm(newEdUposBrdUneg)
rm(newEdUposBrdUpos)

data = droplevels(data)
#======================================================================
#======================================================================
#====================== MINIMIZATION FUNCTION =========================
#======================================================================
#======================================================================

minimize = function(x) {
  
  u1 = x[1]
  u2 = x[2]
  
  #imposing restrictions in the values of the parameters
  if (u1 < 0 |  u2 < 0 | u1 >1 |  u2 >1) {
    return(c(runif(n=1, min=10000, max=1000000),
             runif(n=1, min=10000, max=1000000)))
  }
  
  EdUposBrdUpos = u2/(u1+u2) + u1/(u1+u2)*exp(-(u1+u2)*t)
  EdUposBrdUneg = 1 - EdUposBrdUpos
  EdUnegBrdUpos = u2/(u1+u2) - u2/(u1+u2)*exp(-(u1+u2)*t)
  
  residualsEdUposBrdUpos = EdUposBrdUpos - expEdUposBrdUpos
  residualsEdUposBrdUneg = EdUposBrdUneg - expEdUposBrdUneg
  residualsEdUnegBrdUpos = EdUnegBrdUpos - expEdUnegBrdUpos
  
  residuals = c(residualsEdUposBrdUneg,residualsEdUposBrdUpos,residualsEdUnegBrdUpos)
  return(residuals)
}

#======================================================================
#======================================================================
#========================== MODEL FITTING =============================
#======================================================================
#======================================================================
for (sample in unique(data$Sample)) {
  for (compartment in compartments) {
    sub_data = data[which(data$Sample == sample & grepl(compartment, data$Compartment)),]
    day = unique(sub_data$Day)
    type = paste(unique(sub_data$Type))
    
    expEdUposBrdUpos = sub_data$value[which(grepl("EdUposBrdUpos", sub_data$Compartment))]
    expEdUposBrdUneg = sub_data$value[which(grepl("EdUposBrdUneg", sub_data$Compartment))]
    expEdUnegBrdUpos = sub_data$value[which(grepl("EdUnegBrdUpos", sub_data$Compartment))]
    
    estimates = lsqnonlin(x0 = c(0.5,0.5), 
                          fun = minimize)
    
    u1 = estimates$x[1]
    u2 = estimates$x[2]
    residuals = estimates$ssq
    newrow = c(sample, day, type, compartment, u1, u2, residuals)
    model = rbind(model, newrow)
  }
}


rm(expEdUposBrdUpos)
rm(expEdUposBrdUneg)
rm(expEdUnegBrdUpos)
rm(newrow)
rm(u1)
rm(u2)
rm(residuals)
rm(day)
rm(type)
rm(estimates)
rm(sample)
rm(compartment)
rm(sub_data)

model = data.frame(model)
colnames(model) = c("Sample", "Day", "Type", "Compartment",  "u1", "u2", "residuals")
model$u1 = as.numeric(paste(model$u1))
model$u2 = as.numeric(paste(model$u2))
model$Sphase = 1/model$u1
model$G2MG1phase = 1/model$u2
model$CC = model$Sphase+model$G2MG1phase
rownames(model) = NULL

write.table(model, file = paste(filename, "_model.txt", sep = ""), sep = "\t")

#======================================================================
#======================================================================
#=============================== END ==================================
#======================================================================
#======================================================================