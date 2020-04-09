! 09/12/2011 Alex
! commvar.f90 -- module containing definitions of global variables. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

module commvar
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    

!!!  Global Constants:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Selected Boundary conditions:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), parameter :: periodic_bc = 1    ! constant to denote the periodic boundary conditions 
integer (I4B), parameter :: exact_bc = 2       ! constant to denote the exact boundary conditions when the 
                                               ! distribution function is known exactly at the boundary
integer (I4B), parameter :: diffusive_reflection_bc = 3  ! constant to denote the diffusive reflective boundary conditions 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Global Variables:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Variables to describe zero's level of refinement (So far only refinement in the variable u will be 
!! used. All constants in x are the zero level constants. In U we will have levels of refinement. Parameters of 
!! the refinement will be build in, calculated automatically of solicited somehow separately. However, this file will 
!! define constants to create the mesh of level zero refinement. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

integer (I4B), dimension (:), allocatable :: k_b_list, k_c_list, su_list, sv_list, sw_list, Mu_list, &
 Mv_list, Mw_list, N_list, d_c_list, d_b_list ! arrays of different values of k_c,k_b,s,M,N   

integer (I4B) :: su,sv,sw,Mu,Mv,Mw ! M is the number of velocity cells of level 0, s is the number of gauss points per level zero cell.

integer (I4B) :: k_c,k_b,N ! N is the number of physical cells, k_c is order of the highest Legendre polynomials 
						! used in the basis inside the cell  k_d is the number of Gauss-Lobatto nodes used on  
						! the t=const interface of the cell. 
integer (I4B) :: d_c,d_b ! d_c is the order of the Legendre basis in time variable used inside the cell (for the local solutions/test functions).
						! d_b is the number of Gauss-Lobatto nodes on the x=const cell interface
integer (I4B) :: max_deg_xt ! this is the max degree of the two-dimensional basis in x and t. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (DP), dimension (:), allocatable ::  x_gauss_n, x_gauss_w ! gauss nodes and weights for calculation of spectral coefficients or integrals in u 
real (DP), dimension (:), allocatable ::  u_gauss_n, u_gauss_w ! Gauss nodes and weights for calculation of spectral coefficients or integrals in x
real (DP), dimension (:), allocatable ::  t_gauss_n, t_gauss_w ! Gauss nodes and weights for calculation of spectral coefficients or integrals in t

!
real (DP), dimension (:), allocatable ::  x_lobat_n, x_lobat_w ! Gauss-Lobatto nodes and weights for calculation of spectral coefficients or integrals in u 
real (DP), dimension (:), allocatable ::  t_lobat_n, t_lobat_w ! Gauss-Lobatto nodes and weights for calculation of spectral coefficients or integrals in t
!! NOTE: orders of gauss lobatto quadratures in t and in x are controlled by d_b and k_b. 

real (DP) :: u_L, u_R, v_L, v_R, w_L, w_R ! left and right endpoints of the interval in U,V,and W 
real (DP), dimension (:), allocatable ::   grids_u, grids_v, grids_w 
integer (I4B), dimension (:), allocatable :: grids_cap_u,grids_cap_v,grids_cap_w
     ! grids_cap_u, grids_cap_v and _w contain the total number of 1D meshpoints (including both enpoints, e.g., one cell has two points 
     ! two cells have three and so on. recall that the 3D grid is a cartesian product of the 1d meshes.      
     ! grids_u, grids_v, and grids_w contain 1D meshes for each of the grids. Meshes are following each other, 
     ! to find the appropriate records use the grid capacity array "grids_cap_u/_v/_w
real (DP), dimension (:), allocatable :: cells_lu, cells_lv, cells_lw, cells_ru, cells_rv, cells_rw
integer (I4B), dimension (:), allocatable :: cells_pgrid, cells_cgrid, cells_refu, cells_refv, &
                                      cells_refw, cells_gow, cells_gou, cells_gov
integer (I4B), dimension (:), allocatable :: cells_ugi, cells_vgi, cells_wgi ! relative addresses on the grid... 
     ! cells are the main object of the discretization, they form an ierarchical grid. 
     ! cells_pgrid is the number of the parent grid. 
     ! cells_cgrid is -1 if the cell is not refined or is equal to the number of the child grid. 
     ! cells_lu/_lv/_lw and cells_ru/_rv/_rw are the corners of the cells
     ! cells_refu/_refv/_refw are the refinement coefficients 
     ! cells_gou/_gow /_gov Cell's Gauss order. --- so far will be the same in both u, v, and w 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
real (DP), dimension (:), allocatable ::  nodes_u, nodes_v, nodes_w, nodes_gwts
integer (I4B), dimension (:), allocatable :: nodes_pcell, nodes_ui, nodes_vi, nodes_wi
     ! velocity nodes of the nodal-DG discretization. nodes are listed for each cell in the order cells are listed
     ! nodes_pcell is the number of the parent cell
     ! nodes_u/_v/_w are the coordinates of the node. 
     ! nodes_gwts are the values of products of weights of gaussian quadrature formulas to be used in 3D integration
     ! nodes_ui/_vi/_wi are the values of the local indices that determine addresses within the cell of the 1D nodal 
     ! points in each variable associated with the node
 integer (I4B), dimension (:), allocatable :: nodes_Ashift,nodes_dui,nodes_dvi,nodes_dwi,nodes_phican    
     ! this array contains the shift in the A index to be used for periodic grids/selkf similar basis functions
     ! nodes_dui/_dvi/_dwi contains the shift in cell idex for each node. = Shift in the index to go from the cell that has the node to the canonical cell. 
     ! nodes_phican contains the number of the canonical node on the canonical cell that corresponds to the given node.
!!!!!!!!!!! Parameters of time and space discretization 
  			 	
real (DP) :: x_L, x_R ! left and right endpoints of the interval in X
real (DP), dimension (:), allocatable ::  xmesh_l, xmesh_r ! left and right endpoints of the cells in x (usually, the length of xmesh = N)
real (DP), dimension (:), allocatable ::  xgnodes ! Gauss nodal points for the cells in x. listed for all cells from 
real (DP), dimension (:), allocatable ::  xlnodes ! Gauss-Lobatto nodal points for the interfaces t=const in 
 ! cells in x. 
real (DP) :: t_L, t_R ! left and right endpoints of the interval in T
real (DP) :: e_time ! current time 
real (DP) :: dt ! the timestep
real (DP), dimension (:), allocatable ::  tlnodes  ! Gauss-Lobatto nodal points for the interfaces x=const in


logical :: mesh_x_uniform, mesh_u_uniform, mesh_v_uniform, mesh_w_uniform       ! mesh is uniform in x and u (true or false) (applied only to meshes of level 0)

!!!!!!
! Choice of BCs and the sourse function
integer (I4B) :: selected_rightbcond, selected_leftbcond ! variable to denote the choice of BCs
integer (I4B) :: selected_exact_sol ! variable to denote the choice of the exact solution (if available) 
!!!!!!
! The names of computed solution file and the directory to store it
character (len=50) :: current_solution_dir  ! name of the directory (relative to the local tree?) where the solution is stored
character (len=20) :: current_solution_base_name ! base name for the solution file -- more information about order, mesh, and time is added to it automatically
character (len=20) :: current_Aoperator_base_name ! base name for the files containing Aarrays -- more information about order, mesh, and time is added to it automatically
!!!!
integer (I4B) :: num_save_solution, num_eval_error ! how many records of solution from intial time to final time
                                                   ! how many times to caclulate the error
!!!!
! varaibles related to generation of non-unifom meshes: 

integer (I4B) :: x_nonuniform_mesh_type, u_nonuniform_mesh_type, v_nonuniform_mesh_type, w_nonuniform_mesh_type ! parameter to describe non-uniform mesh in x and u 

! supported types of nonuminform mesh: 
! 1 --- the mesh is build based on the gauss nodes used for the integration of moments, moments_x_gauss_nodes and 
!       moments_u_gauss_nodes. intervals [x_left,x_right], [u_left,u_right] is divided in subintervals 
!       as to have gauss nodes at centerpoints. Some extra points need to be introduced to make this possible.
! 2 --- this mesh is to be used with the diffusive boundary conditions. The velocity of the wall (currently only u_{w}=0)
!       will be included in the mesh. Also, the cell near the walls will be 1/8 - 1/4 - 1/2        
! 3 --  for variable u -- this mesh will have small cells surrounding u=0 as prescribed by parameters
!       sml -- small cell levels amd smr -- small cell refinement factor 
! 3 --- This is a non-niform mesh in "x" with cells near wall be 1/4-1/2. Currenlty is not supported for meshes in "u"

! parameters related to the integration of moments in u and integration in x
!
integer (I4B) :: moments_x_gauss_order, moments_u_gauss_order ! order of gausss quadrateures in x and u to calculated (total) moments
integer (I4B) :: moments_refine_x, moments_refine_u ! coefficients of mesh refinement in x and u for the evaluation of moments
real (DP), dimension (:), allocatable :: moments_x_gauss_nodes,   moments_u_gauss_nodes     ! arrays to keep gauss nodes for evaluation of moments 
real (DP), dimension (:), allocatable :: moments_x_gauss_weights, moments_u_gauss_weights ! arrays to keep gauss weights for evaluation of moments  

real (DP), dimension (:,:), allocatable :: g_nds_all, g_wts_all ! arrays to store vaules of used gauss nodes and weights. 
                                                                 ! g_nds_all(i,j) i is the number of the node and j is the order of the gauss formula 
integer (I4B) :: g_max_order=10 ! maximum gauss order that will be used ... 

! variables describing the gas:
real (DP) :: gasR,gasTref,gasalpha,gasmuref ! is the normal gas constant, gas reference temperature, gas reference viscosity and alpha constant
                                           ! --- need to be read from the parameter.dat
                                           
        !! CHANGE TO READ FROM PARAMETERS.DAT
                                           
real (DP) :: mol_mass = 6.63d-26 ! variable to keep molecular mass 

        !! CHANGE TO READ FROM PARAMETERS.DAT
        
! variables for temperatures on the walls -- will use in Diffusive BCs.
real (DP) :: T_w_left, T_w_right ! temperature of the left and right wall. Only used with the diffusive BCs

! to check whether restart is used.
character (len=20) :: restart_time_txt  ! time of last saved solution/restart time in text format
logical :: need_to_restart                ! tells if the simulations need to be restarted 
 
!!!!!!!!!!!!! operator A and the related staff !!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: A_capphi ! A_capphi(i) gives the number of non-zero entries in A(\xi,\xi_{1},\varphi^{j}_{p})
                                       ! for each basis function with number i
integer (I4B), dimension (:), allocatable :: A_xi,A_xi1 !for each non-zero entry of A this one keeps the index of velocities that produced that non-zero entries. 
                                       ! We assume that all velocities are indexed with one index
integer (I4B), dimension (:), allocatable :: A_phi ! this is the index of the basis function fro which the non-zero entry was computed.                                        
real (DP), dimension (:), allocatable :: A ! this is the collision information operator A(\xi,\xi_{1},\phi^{j}_{p})
                                !  A(i) i -- is the index of nonzero entires, we need to use other A-arrrays to restore 
                                ! what velocities and wnat basis function this index correcpods  
                                ! for example, A_sind_xi(i) gives the first velocity, A_sind_xi1 gives the second velcity
                                ! and A_phi gives the index of the used basis function. 
real (DP) :: Trad ! cutoff radius for A-array (it is actually diameter of the collisio sphere -- distance between xi and xi1)
real (DP) :: ErrEps, ErrChi, min_sens ! errors in integrals in chi and epsilon and cutoff level for A
! Notice that min_sense will have to uses and they normally require values different in order of magnitude. 
! The historically first use of min_sense is in the subroutine for the evaluation of $A$. These values should be decresing 
! with resolution and essentially be equal to ErrChi.
! The second use is for after the truncation of pre-cmoputed A. 
 
integer (I4B), dimension (:), allocatable :: I1_list ! list of basis functions for which A-array should be computed
                                     
!!!!!!!!!!!!!!! Parameters of Hard shepre model
real (DP) :: mol_diam !         molecular diameter                            

!!!!!!!!!!!!!!! Parameters of Dimensionless Reduction
! C_inf = termal velocity in m/s 
! L_inf = characteristic length in m
! N_inf = the total number of molecules in the volume
! T_inf =  the normalization for time is selected from the condition T_inf*C_inf = L_inf calculated automatically
!!!!!!!!!!!!!!!
real (DP) :: T_inf, C_inf, L_inf, N_inf

!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), allocatable :: f ! the distribution function -- nodal values. 
real (DP), dimension (:), allocatable :: f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5 ! the distribution function -- nodal values. 

!!!!!!!!!!!!!!!!!!!
integer (I2B) :: Num_OMP_threads ! number of threads for OpenMP. 

!!!!!!!!!!!!!!!!!!!
integer (I4B) :: rkmts ! The variable to keep the order of the Runge-Kutta method and the Adams-Bashworth MTS schemes

!!!!!!!!!!!!!!!!!!!
integer (I4B) :: run_mode ! This is a variable that keeps track in which mode is the simulations: 
                          ! run_mode=0 is the simulation is strongly non-equilibrium. Set decomp_lev < |f_M-f|
                          ! run_mode=1 is when the simulation are close to equilibrium, that is that the local 
                          !            velocity distribution function is close to a maxwellian, however, it is 
                          !            not quite converged to the Maxwellian yet.  linear_lev < |f_M - f| < decomp_lev.
                          !            Max_decomp decides when the solution need to be decomposed into a Maxwell and $$ 
                          !
                          ! run_mode=2 is the linearized regime. In this regime the function is a perturbation of 
                          !            a Maxwellian and the norm of the pertrurbation is small so that quadratic term in the 
                          !            projection of the collision integral can be neglected. |f_{M}-f|<linear_lev, 
                          !            min_linear is the parameter that decides whether the addition is small enough. 
                          ! run_mode=-777 this is a signal to stop the exectution of an MPI exectution
                          !              
                          
real (DP) :: linear_lev, decomp_lev ! these are parameters detemining at which level linearization or perturbative decomposition starts. 
								    ! make sure that linear_lev is in the order of magnitude smaller than decomp_lev                         
!!!!!!!!!!!!!!!!!!!
real (DP) :: LocDens, LocUbar, LocVbar, LocWbar, LocTempr ! these are the variables to hold the local values of the macroparameters
!!!!!!!!!!!!!!!!!!!
! These are variables related to the linearization of the collision operator
real (DP), dimension (:), allocatable :: fm ! this is a storage to keep the local maxwellian
real (DP), dimension (:,:), allocatable :: fmA, BfmA! this is to keep the linearized operator, BfmA -- is the consolidated linearized operator
integer :: fmA_flag=0,BfmA_flag=0! these are flags if the arrays are allocated -- the size function sometimes messes up..

!!!!!!!!!!!!!!!!!!!                          

!!!!!!!!!! These variables are invoked to save the relative L1 error of the deviation of the distribution function from the local Maxwellian
real (DP), dimension (:), allocatable :: L1_a,L1_t ! L1_a keeps the error, L1_t keeps the time moments that correspond to the erro values
integer (I4B) :: L1_count ! This is the counter to keep track of the records
real (DP) :: L1_record_period,L1_next_time_record,L1_err_eval_period,L1_next_time_eval ! some scrap timing variables.

!!!!!!!!!!!!!!!!!!! variables for chunking the A-arrays
integer (I4B) :: numchnks ! the number of chunks where A arrays are saved. 


!!!!!!!!!!!!!!!!!!! MPI Variables for the evaluation of collision operator !!!!!!!!!
integer :: MPI_LINEAR_COLL_WORLD
integer :: MPI_LINEAR_COLL_WORLD_GROUP ! Communications and commpnicator group handles for creating a communicator
                                        ! That will contain processes participating in consolidated evaluation of linearied operator.  
integer, dimension(:), allocatable :: MPI_Acopy_UNIVS, MPI_Acopy_UNIVS_GROUP ! This array will hold pointers to universes that will receive a single copy of A_array.
integer, dimension(:), allocatable :: MPI_LinFmA_UNIVS, MPI_LinFmA_UNIVS_GROUP ! This array will hold pointers to universes that will receive a single copy of A_array.

                                          ! the first element will keep the number of receives and the other elements will keep the number of the universes
integer (I4B), dimension (:), allocatable :: procs_nodes_wld ! this array will store nodes that are assigned for a processor 
                                                               ! on which the program is running. the index runs trhought the nodes assigned for the 
                                                               ! evaluation. It is expected that the node knows the necessary peices of arrys A

integer (I4B), dimension (:), allocatable :: lin_proc_nodes_wld ! this array will collect the information about the local workload of the each processor in the 
                                                              ! group of processors responsible for the evaluation of the linearized collision operator
                                                              ! the index runs over assigned nodes. 
integer (I4B) :: My_Acopy_Univ_indx ! these keep the index of the Acopy universe that contain this proceesor   
integer  :: Lin_Self_Send_FmA_Flag ! this flag =0 if the master node of a Acopy Universe also appears as a processor recieving 
									! components of the linearized operator from this universe.   

integer (I4B), dimension (:,:), allocatable :: nodes_proc_groups ! this array will contain the beaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b
                     
integer (I4B), dimension(:), allocatable :: LinFmA_Univs_callorder ! entries of this array determine the order in wich LinFmA universes enter the 
                     ! collective communication. The entires are numbers from 1 to num_Acopies that are permuted in some way. 
                                                               
integer (I4B), dimension (:), allocatable ::  linprocndswld_Acopy_univs, linprocndswld_BfmA_addr, linprocndswld_Acopy_addr
        ! These arrays are located on the group of processors computing linearized solutions
        ! they store local information about what components from what universice the processor is 
        ! expecting as well as the information about the addresss of the components in the 
        ! broadcasted arrays 
        ! receive_proc_lin and send_process_lin
        !! in these arrays the information is in the following format:
        !! (# number of the LinFmA_recv_univ universe where the components of the linearized operator will come from, 
        !!  # total number of nodes for which the componets of the linearizeed operator need to be recieved from this universe
        !!  #node1, #node2, and so on, #nodeN, 
        !! then continue onto the next universe 
        !! the stucture of the array linprocndswld_Acopy_addr is the same as of linprocndswld_Acopy_univs
        !! exception that instead of the nodes, we have their indices in the corresponding Acopy universe's 
        !! Acopy_wrkld array. This array will be useful for the preparation of the linearized operator
        !! linprocndswld_BfmA_addr has the same structure except the records are indices of nodes from array LinFmA_recv_univ
        !! as they show up in the array BfmA
        !!
        !!
integer (I4B), dimension (:), allocatable :: Acopy_wrkld 
        !! This array is very similar to the array procs_nodes_wld, with the exception that is applies to the ]
        !! Acopy_universe rather than to individual processor. 
        !! Specifically, each Acopy universe will have some nodes assigned to it. All processors in this universe 
        !! will have a copy of this array containing the nodes that are assigned to their Acopy universe
        !! no nodes belong to two Acopy universes at the same time and no two processors belong to the same Acopy 
        !! universes, to there will be no ambiguity
integer (I4B), dimension (:), allocatable :: procs_nodes_wld_Acopy_addr
        !! each slave processor will contribute components of the linearized operator 
        !! in each Acopy universe an array will be created (FmA) containing the components of the linearized operator 
        !! the rows of this array will correspond to the nodes listed in the Acopy_wrkld array
        !! now Acopy_wrkld_addr will contain indices of nodes from the procs_nodes_wld on this processor in the 
        !! its A copy Universe's array Acopy_wrkld. The length of the  procs_nodes_wld_Acopy_addr is tha same as      
        !! array procs_nodes_wld. 
       
integer (I4B) :: num_lin_proc ! this variable will keep the numer of processors dedicated to the evaluation of the linearized operator. This number must be less than or equal to the 
                              ! total number of processors in the MPI universe. 
integer (I4B) :: num_Acopies ! this variable will store the desired number of times A is copied in memory between the processor. 
                             ! this is a good approach is A does not take a lot of storage, so say it fits in one node memory, say 16 procs.  
                             ! however there are many nodes and essentially we chunk th enodes between groups of processors. each gorup will have a copy of entire A                              
integer (I4B) :: num_linstep_wait,lin_no_update ! this is the numeber of steps to wait until start computing linear solutio using consolidated linearied collision operator                             
              ! lin_no_update is the counter to see how many times we computed linearied collision operator without updating it.
integer :: Acopy_irank ! MPI variable for universe Acopy. Gives the rank of this processor in the Acopy universe it belongs. Does not make sense on master node i.e. not valid with irank=0!!!                


!!!!!!!!!!!!!!!!!! There are some diagnostic variables !!!!!!!!!!!!!!
real (DP), dimension(:), allocatable :: MomDerivRelaxT ! array to store values of the derivatives of the moments and the relaxation speeds
integer :: MomDerivRelaxT_flag = 0 ! flag to denote if MomDerivRelaxT has been created 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              
end module commvar