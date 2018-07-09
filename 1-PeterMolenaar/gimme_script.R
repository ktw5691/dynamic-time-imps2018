#install gimme
install.packages("gimme", dependencies = TRUE)

#load gimme
library(gimme)

#tell gimme where your input and output folders are located
###CHANGE WORKING DIRECTORY PATH
setwd("C:\\Users\\pxm21\\Desktop\\GIMME Symposium")

#run gimme without subgrouping
###CHANGE DATA AND OUT DIRECTORY
gimmeSEM(data = "play_data",
out = "gimmer_output2",
sep = "",
header = FALSE,
ar = FALSE,
plot = TRUE,
subgroup = FALSE,
paths = NULL,
groupcutoff = .75,
subcutoff = .5) 

#run gimme with subgrouping
###CHANGE DATA AND OUT DIRECTORY
gimmeSEM(data = "play_data",
out = "subgroup_output",
sep = "",
header = FALSE,
ar = TRUE,
plot = TRUE,
subgroup = TRUE,
paths = NULL,
groupcutoff = .75,
subcutoff = .5)
