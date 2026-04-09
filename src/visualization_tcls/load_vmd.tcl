package require topotools
package require pbctools

mol new loaded_npt_final.lmp type lammpsdata waitfor all
mol addfile traj_loaded.lammpstrj type lammpstrj waitfor all

set m [molinfo top]

pbc unwrap -all
pbc wrap -all -compound fragment -center com -centersel "type 1 2"
mol reanalyze top

mol reanalyze $m

mol delrep 0 $m

mol representation Lines 2.0
mol selection "type 1 2"
mol color Name
mol addrep $m

mol representation CPK 0.8 0.3 12.0 12.0
mol selection "type 3 4 5 6"
mol color Type
mol addrep $m
