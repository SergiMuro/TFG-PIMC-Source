/*MC.cfg*/
��������������� type of Monte Carlo simulation �����������Ŀ
� (1) - Diffusion M.C.       � (0) - Variational M.C.      �
����������������������������������������������������������Ĵ
�  MC= 0                                                   �
������������������������ algorithm �����������������������Ĵ
� (0) - quadratic            � (0) - Move all              �
� (1) - linear (E at drift)  � (1) - Move one by one       �
� (2) - linear (E at drift/2)�                             �
� (3) - pseudopotential      �                             �
� (4) - Metropolis (at drift)�                             �
� (5) - Metropolis (drift/2) �                             �
����������������������������������������������������������Ĵ
�  SmartMC= 1                                              �
�������������������������� walkers �����������������������Ĵ
�  Npop= 1 00     - anticipated number of walkers          �
�               1 - means absence of branching             �
�  NwalkersMax= 5 0   - maximal number of walkers          �
�������������� store coord for par optimization ����������Ĵ
�  optimization= 0            Niter_store= 1000            �
����������������������� parameters �����������������������Ĵ
�  var_par_array= 0 let the branching choose best var.par. �
������������������������������������������������������������

������������������� general parameters �������������������Ŀ
� N= 120       - number of particles                       �
� Ndens= 120       - used to calculate box size            �
� Nspin= 4      - number of spins                          �
����������������� Homogeneous system ���������������������Ĵ
� n= 0.2        - density                                  �
� bilayer_width= 0.6                                       �
������������������������������������������������������������

������������������� boundary conditions ������������������Ŀ
� (0) - no boundary conditions (3D trap)                   �
� (1) - periodic boundary conditions (infinite system)     �
� (2) - periodic boundary condition in z direction         �
����������������������������������������������������������Ĵ
� boundary condition parameter is defined in main.h        �
������������������������������������������������������������

������������������� trial wave functions �����������������Ŀ
� Trial wavefunction type is defined in main.h             �
� grid_trial= 1000 0       - number of points in trial w.f.�
����������������������������������������������������������Ĵ
�   Rpar= 0.0       maximal distance for same layer        �
�   Rpar12= 2.4      maximal distance for differnt layers  �
�              if equal to zero, is set to L/2             �
�   Epotcutoff= 1000 0                                     �
�   R0par= 5.0  Effective repulsive interaction spin1=spin2�
����������������������������������������������������������Ĵ
� One-body term: (1) oscillator (2) lattice (3) Crystal    �
� alpha_x_up= 0              alpha_x_dn= 0  alpha_Rx= 0.25 �
� alpha_y_up= 0              alpha_y_dn= 0  alpha_Ry= 0.25 �
� alpha_z_up= 0.             alpha_z_dn= 0. alpha_Rz= 0.0  �
������������������������������������������������������������

���������������� information while run �������������������Ŀ
� video= 1  text (0) / video (1) mode                      �
�  in text mode the ammount of information                 �
�  (0) minimal (1) normal (2) high (3) highest             �
� verbosity= 0                                             �
������������������������������������������������������������

����������������� Initial configuration ������������������Ŀ
� (1) generate new coordinates or (0) load from file       �
� generate_new_coordinates= 1                              �
� in the latter case specify the filename                  �
�  file_particles= in3Dprev.in                             �
������������������������������������������������������������

����������������������� timesteps ������������������������Ŀ
�          DMC               �             VMC             �
����������������������������������������������������������Ķ
� dt= 1e-3                   �    acceptance_rate= 50      �
�                            � or set it to 0 and specify  �
�                            �     dt_all= 0.01            �
�                            �     dt_one= 0.01            �
���������������������� grid sizes ������������������������Ĵ
�  gridOBDM= 50                                            �
�  gridTBDM= 30                                            �
�  gridTBDM_MATRIX= 10                                     �
�  gridOBDM_MATRIX= 20                                     �
�  gridPD= 100                                             �
�  gridRD= 100                                             �
�  gridSk= 56                                              �
�  gridNk= 20                                              �
�  gridSD= 100                                             �
�  gridg3= 50                                              �
����������������������� files ����������������������������Ĵ
�  file_particles= in3Dprev.in  - coordinates              �
�  file_wf=      inwf.in      - trial wavefunction         �
�  file_energ= oute.dat        - energy                   �
�  file_OBDM= outdr.dat         - OBDM nondiagonal elements�
�  file_OBDM_MATRIX= outobdm.dat - one body density matrix �
�  file_PD=  outpd.dat          - pair distribution        �
�  file_RD=  outrd.dat          - radial distribution      �
�  file_PDz= outpdz.dat         - pair distribution (z)    �
�  file_RDz= outrdz.dat         - radial distribution (z)  �
�  file_R2=  outr2.dat          - sqrt<r^2>                �
�  file_z2=  outz2.dat          - sqrt<z^2>                �
�  file_PD_pure= outpdp.dat     - pair distribution pure   �
�  file_PD_accum= outpda.dat     - pair distribution pure  �
�  file_R2_pure= outr2pur.dat   - sqrt<r^2> pure           �
�  file_z2_pure= outz2pur.dat   - sqrt<z^2> pure           �
��������������� quantities to be measured ����������������Ĵ
� measure_energy= 1                                        �
� measure_RadDistr= 0                                      �
� measure_PairDistr= 1                                     �
� measure_g3= 0                                            �
� measure_Sk= 0                                            �
� measure_Nk= 0                                            �
� measure_Lind= 0  (Lindemann ratio)                       �
� measure_PairDistr_pure= 0                                �
� measure_PairDistr_accum= 0                               �
� measure_OBDM= 0                                          �
� measure_TBDM= 0                                          �
� measure_TBDM_MATRIX= 0                                   �
� measure_OBDM_MATRIX= 0                                   �
� measure_R2= 0                                            �
� measure_Sk_pure= 0                                       �
� measure_OP= 0  (order parameter)                         �
� measure_SD= 0                                            �
�������� number of McMillan points in OBDM measurement ���Ĵ
� McMillan_points= 10                                      �
�������� width of the strip used in RD measurement �������Ĵ
� RDmax= 10     length of the array                        �
� RDzmax= 10                                               �
� RDwidth=  10                                             �
� RDzwidth= 10                                             �
� PDmax= 0   (0 means L/2)                                 �
� OBDMmax= 0   (0 means L/2)                               �
� TBDMmax= 0   (0 means L/2)                               �
� TBDM_MATRIXmax= 0   (0 means L/2)                        �
� OBDM_MATRIXmax= 5   (0 means L/2)                        �
� SDspacing= 100                                           �
��������������������� Measurements �����������������������Ĵ
� gridPD_pure_block= 10 00                                 �
� blck_heating= 0 - number of blocks spend for heating    �
� blck= 2    number of blocks spend for doing measurements�
� Niter= 5 000         - number of iterations                  �
� Nmeasure= 1 0     do measurement each Nmeasurent iterations�
� file_append= 0                                           �
������������������������������������������������������������
