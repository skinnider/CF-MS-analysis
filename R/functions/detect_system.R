# detect system
system = 'sockeye'
node = Sys.info()[["nodename"]]
if (grepl("cedar", node)) {
  system = 'cedar'
} else if (grepl("frontend|worker", node)) {
  system = 'elasti'
}
# set up base directory
if (system == 'cedar') {
  base_dir = "/home/skinnim/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis"
  args$allocation = 'rrg-ljfoster-ab'
} else if (system == 'sockeye') {
  base_dir = "/scratch/st-ljfoster-1/CF-MS-analysis"
  args$allocation = 'st-ljfoster-1'
} else if (system == 'elasti') {
  base_dir = '/home/ubuntu/projects/rrg-ljfoster-ab/skinnim/CF-MS-analysis'
  args$allocation = 'root'
}
