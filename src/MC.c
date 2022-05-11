/*MC.cfg*/
ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄ type of Monte Carlo simulation ÄÄÄÄÄÄÄÄÄÄÄÄ¿
³ (1) - Diffusion M.C.       ³ (0) - Variational M.C.      ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÁÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³  MC= 0                                                   ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ algorithm ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³ (0) - quadratic            ³ (0) - Move all              ³
³ (1) - linear (E at drift)  ³ (1) - Move one by one       ³
³ (2) - linear (E at drift/2)³                             ³
³ (3) - pseudopotential      ³                             ³
³ (4) - Metropolis (at drift)³                             ³
³ (5) - Metropolis (drift/2) ³                             ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÁÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³  SmartMC= 1                                              ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ walkers ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³  Npop= 1 00     - anticipated number of walkers          ³
³               1 - means absence of branching             ³
³  NwalkersMax= 5 0   - maximal number of walkers          ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄ store coord for par optimization ÄÄÄÄÄÄÄÄÄÄÄ´
³  optimization= 0            Niter_store= 1000            ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ parameters ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³  var_par_array= 0 let the branching choose best var.par. ³
ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ

ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ general parameters ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
³ N= 120       - number of particles                       ³
³ Ndens= 120       - used to calculate box size            ³
³ Nspin= 4      - number of spins                          ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ Homogeneous system ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³ n= 0.2        - density                                  ³
³ bilayer_width= 0.6                                       ³
ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ

ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ boundary conditions ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
³ (0) - no boundary conditions (3D trap)                   ³
³ (1) - periodic boundary conditions (infinite system)     ³
³ (2) - periodic boundary condition in z direction         ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³ boundary condition parameter is defined in main.h        ³
ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ

ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ trial wave functions ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
³ Trial wavefunction type is defined in main.h             ³
³ grid_trial= 1000 0       - number of points in trial w.f.³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³   Rpar= 0.0       maximal distance for same layer        ³
³   Rpar12= 2.4      maximal distance for differnt layers  ³
³              if equal to zero, is set to L/2             ³
³   Epotcutoff= 1000 0                                     ³
³   R0par= 5.0  Effective repulsive interaction spin1=spin2³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³ One-body term: (1) oscillator (2) lattice (3) Crystal    ³
³ alpha_x_up= 0              alpha_x_dn= 0  alpha_Rx= 0.25 ³
³ alpha_y_up= 0              alpha_y_dn= 0  alpha_Ry= 0.25 ³
³ alpha_z_up= 0.             alpha_z_dn= 0. alpha_Rz= 0.0  ³
ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ

ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ information while run ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
³ video= 1  text (0) / video (1) mode                      ³
³  in text mode the ammount of information                 ³
³  (0) minimal (1) normal (2) high (3) highest             ³
³ verbosity= 0                                             ³
ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ

ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ Initial configuration ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
³ (1) generate new coordinates or (0) load from file       ³
³ generate_new_coordinates= 1                              ³
³ in the latter case specify the filename                  ³
³  file_particles= in3Dprev.in                             ³
ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ

ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ timesteps ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
º          DMC               º             VMC             º
ÇÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ×ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¶
º dt= 1e-3                   º    acceptance_rate= 50      º
º                            º or set it to 0 and specify  º
º                            º     dt_all= 0.01            º
º                            º     dt_one= 0.01            º
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ grid sizes ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³  gridOBDM= 50                                            ³
³  gridTBDM= 30                                            ³
³  gridTBDM_MATRIX= 10                                     ³
³  gridOBDM_MATRIX= 20                                     ³
³  gridPD= 100                                             ³
³  gridRD= 100                                             ³
³  gridSk= 56                                              ³
³  gridNk= 20                                              ³
³  gridSD= 100                                             ³
³  gridg3= 50                                              ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ files ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³  file_particles= in3Dprev.in  - coordinates              ³
³  file_wf=      inwf.in      - trial wavefunction         ³
³  file_energ= oute.dat        - energy                   ³
³  file_OBDM= outdr.dat         - OBDM nondiagonal elements³
³  file_OBDM_MATRIX= outobdm.dat - one body density matrix ³
³  file_PD=  outpd.dat          - pair distribution        ³
³  file_RD=  outrd.dat          - radial distribution      ³
³  file_PDz= outpdz.dat         - pair distribution (z)    ³
³  file_RDz= outrdz.dat         - radial distribution (z)  ³
³  file_R2=  outr2.dat          - sqrt<r^2>                ³
³  file_z2=  outz2.dat          - sqrt<z^2>                ³
³  file_PD_pure= outpdp.dat     - pair distribution pure   ³
³  file_PD_accum= outpda.dat     - pair distribution pure  ³
³  file_R2_pure= outr2pur.dat   - sqrt<r^2> pure           ³
³  file_z2_pure= outz2pur.dat   - sqrt<z^2> pure           ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄ quantities to be measured ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³ measure_energy= 1                                        ³
³ measure_RadDistr= 0                                      ³
³ measure_PairDistr= 1                                     ³
³ measure_g3= 0                                            ³
³ measure_Sk= 0                                            ³
³ measure_Nk= 0                                            ³
³ measure_Lind= 0  (Lindemann ratio)                       ³
³ measure_PairDistr_pure= 0                                ³
³ measure_PairDistr_accum= 0                               ³
³ measure_OBDM= 0                                          ³
³ measure_TBDM= 0                                          ³
³ measure_TBDM_MATRIX= 0                                   ³
³ measure_OBDM_MATRIX= 0                                   ³
³ measure_R2= 0                                            ³
³ measure_Sk_pure= 0                                       ³
³ measure_OP= 0  (order parameter)                         ³
³ measure_SD= 0                                            ³
ÃÄÄÄÄÄÄÄ number of McMillan points in OBDM measurement ÄÄÄÄ´
³ McMillan_points= 10                                      ³
ÃÄÄÄÄÄÄÄ width of the strip used in RD measurement ÄÄÄÄÄÄÄÄ´
³ RDmax= 10     length of the array                        ³
³ RDzmax= 10                                               ³
³ RDwidth=  10                                             ³
³ RDzwidth= 10                                             ³
³ PDmax= 0   (0 means L/2)                                 ³
³ OBDMmax= 0   (0 means L/2)                               ³
³ TBDMmax= 0   (0 means L/2)                               ³
³ TBDM_MATRIXmax= 0   (0 means L/2)                        ³
³ OBDM_MATRIXmax= 5   (0 means L/2)                        ³
³ SDspacing= 100                                           ³
ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ Measurements ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´
³ gridPD_pure_block= 10 00                                 ³
³ blck_heating= 0 - number of blocks spend for heating    ³
³ blck= 2    number of blocks spend for doing measurements³
³ Niter= 5 000         - number of iterations                  ³
³ Nmeasure= 1 0     do measurement each Nmeasurent iterations³
³ file_append= 0                                           ³
ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
