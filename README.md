# Electron-Defect Interaction Calculator

## Installation

 EDIC compiles on top of Quantum Espresso. To install, copy folder to the main folder of QE 6.8, then run make.
 ./configure --with-hdf5  --with-hdf5-include=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include

# from old files
calc_m-ml-nc.f90                      : rj code from EDSC
calc_m-mnl-soc.f90                    : rj code from EDSC
calc_m.f90                            : rj code from EDSC
calcmdefect-all.f90                   : subroutine from calcmdefect.f90
calcmdefect-read_perturb_file.f90     : subroutine from calcmdefect.f90
calcmdefect-get_vloc_colin.f90        : subroutine from calcmdefect.f90
calcmdefect-get_gind_rhoandpsi_gw.f90 : subroutine from calcmdefect.f90 with small change for header etc
calcmdefect-getvloc.f90               : subroutine from calcmdefect.f90
calcmdefect-gw_eps_init.f90           : subroutine from calcmdefect.f90
calcmdefect-h5gw_read.f90             : subroutine from calcmdefect.f90
calcmdefect-mcharge-lfa.f90           : subroutine from calcmdefect.f90
calcmdefect-mcharge.f90               : subroutine from calcmdefect.f90
calcmdefect-ml-nc.f90                 : subroutine from calcmdefect.f90
calcmdefect-ml-rs.f90                 : subroutine from calcmdefect.f90
calcmdefect-mnl-ks.f90                : subroutine from calcmdefect.f90
calcmdefect-mnl-nc.f90                : subroutine from calcmdefect.f90
calcmdefect-mnl-soc.f90               : subroutine from calcmdefect.f90
calcmdefect.f90                       : original 

# new files
edic.f90                              : new program file
edic_mod.f90                          : all new mod file
get_gind_rhoandpsi_gw.f90             : adapt from calcmdefect-get_gind_rhoandpsi_gw.f90 
getvloc.f90                           : adapt from calcmdefect-getvloc.f90
gw_eps_init.f90                       : adapt from calcmdefect-gw_eps_init.f90   
h5gw_read.f90                         : adapt from calcmdefect-h5gw_read.f90     
init_us_2_sc.f90                      : init supercell 
initialization.f90                    : integrated into edic.f90
makefile                              : 
mcharge-lfa.f90                       : adapt from calcmdefect-mcharge-lfa.f90    
mcharge.f90                           : adapt from calcmdefect-mcharge.f90       
ml-nc.f90                             : adapt from calcmdefect-ml-nc.f90         
ml-rs.f90                             : adapt from calcmdefect-ml-rs.f90         
mnl-ks.f90                            : adapt from calcmdefect-mnl-ks.f90        
mnl-nc.f90                            : adapt from calcmdefect-mnl-nc.f90        
mnl-soc.f90                           : adapt from calcmdefect-mnl-soc.f90       
read_perturb_file.f90                 : combine EDSC read_perturb_file.f90 get_vloc_colin.f90
sh










