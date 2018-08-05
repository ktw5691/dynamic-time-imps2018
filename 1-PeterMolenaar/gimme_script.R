library(pacman)
p_load(gimme)

# Run gimme without subgrouping
### NOTE: CHANGE DATA AND OUT DIRECTORY IF NEEDED
gimmeSEM(
  data = "1-PeterMolenaar/play_data",
  out = "gimmer_output2",
  sep = "",
  header = FALSE,
  #ar = FALSE,
  ar = TRUE, # Recommended by Peter Molenaar to set ar = TRUE
  plot = TRUE,
  subgroup = FALSE,
  paths = NULL,
  groupcutoff = .75,
  subcutoff = .5)

# Run gimme with subgrouping
### NOTE: CHANGE DATA AND OUT DIRECTORY IF NEEDED
gimmeSEM(
  data = "play_data",out = "subgroup_output",
  sep = "",
  header = FALSE,
  ar = TRUE,
  plot = TRUE,
  subgroup = TRUE,
  paths = NULL,
  groupcutoff = .75,
  subcutoff = .5)
