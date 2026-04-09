# Load TopoTools and read the topology
package require topotools
topo readlammpsdata loaded_npt_final.lmp full

# Load the trajectory into the same molecule
mol addfile traj_loaded.lammpstrj type lammpstrj waitfor all

