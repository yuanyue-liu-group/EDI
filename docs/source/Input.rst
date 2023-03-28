Input
=====

.. _installation:

An example input file is shown below:

.. code-block:: console

    &calcmcontrol
    prefix='mos2',
    outdir='dout/'
    lvacalign=.true.
    vac_idx=0
    lcorealign=.false.
    core_v_d=0.0
    core_v_p=0.0 
    wt_filename='wt.dat'
    klist_filename='scfklist.dat'
    ev_filename='v.dat'
    noncolin =.true.
    lspinorb =.true.
    calcmlocal = .true.
    calcmnonlocal = .true.
    V_d_filename='./V_d.dat'
    Bxc_1_d_filename='./Bxc_1_d.dat'
    Bxc_2_d_filename='./Bxc_2_d.dat'
    Bxc_3_d_filename='./Bxc_3_d.dat'
    V_p_filename='./V_p.dat'
    Bxc_1_p_filename='./Bxc_1_p.dat'
    Bxc_2_p_filename='./Bxc_2_p.dat'
    Bxc_3_p_filename='./Bxc_3_p.dat'
    calcmcharge=.true.
    mcharge_dolfa=.true.
    qeh_eps_filename='./eps.dat'
    eps_type='gw'
    eps_type='qeh'
    eps_type='tf'
    doqeh=.true.
    dogwfull=.true.
    dogwdiag=.true.
    k0screen_read=0.27
    gw_epsmat_filename='./epsmat.h5'
    gw_eps0mat_filename='./eps0mat.h5'
    /


Variables
------------

The detailed meaning of all variables are listed below in the table:

====================      ======================================
variable                     meaning                            
====================      ======================================
prefix                       qe prefix
outdir                       qe outdir
lvacalign                    vacuum aligment
vac_idx                     vacuum alignment location
lcorealign                   core alignment
core_v_d                    core alignment value for defect
core_v_p                    core alignment value for pristine
wt_filename                  weight file
klist_filename               k point list file
ev_filename                  
noncolin                     non-colinear calculation
lspinorb                     spin-orbit calculation
calcmlocal                   calculate local part M
calcmnonlocal                   calculate non-local part M
V_d_filename                 defect system local potential 
Bxc_1_d_filename             defect system magnetic field
Bxc_2_d_filename             defect system magnetic field
Bxc_3_d_filename             defect system magnetic field
V_p_filename                 pristine system local potential
Bxc_1_p_filename             pristine system magnetic field
Bxc_2_p_filename             pristine system magnetic field
Bxc_3_p_filename             pristine system magnetic field
calcmcharge                 charged defect
mcharge_dolfa               LFA approximation
qeh_eps_filename            dielectric function file from QEH
doqeh                       use QEH dielectric function 
dogwfull                    use BGW dielectric function
k0screen_read               Lindhard model carrier screening
gw_epsmat_filename          BGW dielectric function file
gw_eps0mat_filename          BGW dielectric function file
====================      ======================================


