!
!  miscset.f90
! 
!  Alex 09/13/2011
!
!  Miscellaneous setup for the UWVFbgk code. 
!  constains generation of meshes and more. 
!!!!!!!!!!!!!!!!!!!!!!1

module miscset

use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetDGblzmGnodes
! 
! This subroutine sets up the notes used in various integrations and 
! also used as ordimnates/nodal values 
!
! This subroutine uses parameters set by SetUWbgkParams
!
! This subroutine sets arrays used in other programs -- must go before 
! SetDGVblzmmesh
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetDGblzmGnodes
use commvar, only: x_gauss_n, x_gauss_w, t_gauss_n, t_gauss_w, &
             x_lobat_n, x_lobat_w, t_lobat_n, t_lobat_w, k_c, k_b, d_c, d_b, &
             g_nds_all, g_wts_all, g_max_order, moments_u_gauss_nodes, moments_u_gauss_weights, &
             moments_x_gauss_nodes,moments_x_gauss_weights,moments_u_gauss_order,moments_x_gauss_order
use gaussian_mod
use gaussquad, only: gauss_rule

integer (I4B) :: loc_alloc_stat, info  ! local variable to keep the allocation status
real (DP), dimension (1:k_b) :: xb ! scratch array to use in the subroutine generating gauss-lobatto nodes. 
real (DP), dimension (1:d_b) :: tb ! scratch array to use in the subroutine to generate gauss-Lobatto nodes. 
real (DP), dimension (:), allocatable :: ugn, ugw ! scratch arrays to create gauss nodes... 
integer (I4B) :: i ! local counter

!!! Prepare the Gaussian nodes and weights in x variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
  allocate (x_gauss_n(1:max(k_c,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (x_gauss_n)"
     end if 
     !
  allocate (x_gauss_w(1:max(k_c,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (x_gauss_w)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  if (k_c == 1) then 
   x_gauss_n(:)= (/ 0.0_DP /)
   x_gauss_w(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),x_gauss_n,x_gauss_w)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
  end if 
!!!
!!! Prepare the Gaussian nodes and weights in u variables !!!!!!!!!!!!!!!!!!!!!!
!   first we need to allocate the arrays to store gauss nodes and weights: 
!
  allocate (g_nds_all(1:g_max_order, 1:g_max_order), stat=loc_alloc_stat) 
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (g_nds_all)"
     end if 
     !
  allocate (g_wts_all(1:g_max_order, 1:g_max_order), stat=loc_alloc_stat) 
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (g_wts_all)"
     end if 
     !
!!! reset the arrays to zeros:
g_nds_all=0
g_wts_all=0
!!! 
g_nds_all(1,1)=0.0_DP
g_wts_all(1,1)=2.0_DP
do i=2,g_max_order 
 ! will need two scrap arrays... 
 allocate (ugn(1:i), stat=loc_alloc_stat) 
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (ugn), i=", i
     end if 
 allocate (ugw(1:i), stat=loc_alloc_stat) 
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (ugw), i=", i
     end if 
call GauLeg(Real(-1,DP),Real(1,DP),ugn,ugw)
g_nds_all(1:i,i)=ugn
g_wts_all(1:i,i)=ugw
deallocate(ugn,ugw)
end do 
!!! Prepare the Gaussian nodes and weights in t variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
  allocate (t_gauss_n(1:max(d_c,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (t_gauss_n)"
     end if 
     !
  allocate (t_gauss_w(1:max(d_c,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (t_gauss_w)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  if (d_c == 1) then 
   t_gauss_n(:)= (/ 0.0_DP /)
   t_gauss_w(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),t_gauss_n,t_gauss_w)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
  end if 
!!!!!!!
! 
! Next we will set the Gauss-Lobatto nodes used in the nodal formulation. 
!!! Allocate arrays for gauss-lobatto nodes and weights
  allocate (t_lobat_n(1:max(d_b,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (t_lobat_n)"
     end if 
     !
  allocate (t_lobat_w(1:max(d_b,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (t_lobat_w)"
     end if 
   call gauss_rule(0, d_b, t_lobat_n, t_lobat_w, tb, 0.0_DP, 0.0_DP, 'B', info)  
!!! Allocate arrays for gauss-lobatto nodes and weights
  allocate (x_lobat_n(1:max(k_b,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (x_lobat_n)"
     end if 
     !
  allocate (x_lobat_w(1:max(k_b,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (x_lobat_w)"
     end if 

   call gauss_rule(0, k_b, x_lobat_n, x_lobat_w, xb, 0.0_DP, 0.0_DP, 'B', info)
   ! We set the variables of moments_x_gauss_nodes, moments_u_gauss_nodes, moments_x_gauss_weights, moments_u_gauss_weights
   allocate (moments_x_gauss_nodes(1:moments_x_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (moments_x_gauss_nodes)"
     end if 
     !
  allocate (moments_x_gauss_weights(1:moments_x_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (moments_x_gauss_nodes)"
     end if 
 if (moments_x_gauss_order == 1) then 
   moments_x_gauss_nodes(:)= (/ 0.0_DP /)
   moments_x_gauss_weights(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),moments_x_gauss_nodes,moments_x_gauss_weights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
  end if 
 allocate (moments_u_gauss_nodes(1:moments_u_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (moments_u_gauss_nodes)"
     end if 
     !
  allocate (moments_u_gauss_weights(1:moments_u_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (moments_u_gauss_nodes)"
     end if 
 if (moments_u_gauss_order == 1) then 
   moments_u_gauss_nodes(:)= (/ 0.0_DP /)
   moments_u_gauss_weights(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),moments_u_gauss_nodes,moments_u_gauss_weights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
  end if 

end subroutine SetDGblzmGnodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetDGVblzmmesh 
!
! This subroutine sets up the meshes. Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters 
! N,M, mesh_x_uniform, mesh_u_uniform, moments_x_gauss_order, 
! moments_u_gauss_order are selected                               !!!!! 
! CALL this SUBROUTINE before other subroutines use varaibles N,M.
!
! DO not call before arrays moments_x_gauss_nodes,moments_x_gauss_weights_x,moments_x_gauss_order
!          moments_u_gauss_order,moments_u_gauss_nodes,moments_u_gauss_weights are selected! 
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetDGVblzmmesh 
use commvar, only:  xmesh_l, xmesh_r, N, Mu, Mv, Mw, x_L, x_R, u_L, u_R, &
                   v_L, v_R, w_L, w_R, grids_cap_u, grids_cap_v, grids_cap_w, grids_u, grids_v, grids_w, & 
                   mesh_x_uniform, mesh_u_uniform, x_nonuniform_mesh_type, u_nonuniform_mesh_type,&
                   moments_x_gauss_nodes, moments_u_gauss_nodes,mesh_v_uniform, mesh_w_uniform, & 
                   v_nonuniform_mesh_type,w_nonuniform_mesh_type

intrinsic Min, Max, MinVal
                   
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) ::  i,j       ! local counter 
integer (I4B) :: sml,smr   ! parameters for umesh with small cells around zero sml -- levels of refinment and smr -- ratio of refinement
integer (I4B) :: mxgo,mugo      ! local variables to keep the order of the gauss integration for the moments
real (DP) :: dx, du ! local temp variables to store 1/3 smallest distance between nodes.
 ! first the mesh in x
 if (mesh_x_uniform) then 
  allocate (xmesh_l(1:N), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
  stop
  end if
     !
  xmesh_l = (/ (x_L + (x_R - x_L )/Real(N,DP)*i, i=0,N-1) /)
  !
  allocate (xmesh_r(1:N), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
  stop
  end if
     !
  xmesh_r = (/ (x_L + (x_R - x_L )/Real(N,DP)*i, i=1,N) /)
  !
 else
  select case (x_nonuniform_mesh_type)
  case (1) ! use "the gauss nodes for evaluating moments in x" to set up the "mesh in x"
           ! interval [x_left,x_right] is divided in subintervals as to have gauss nodes at centerpoints 
           ! (and some extra points need to be introduced to make this possible). 
  mxgo = size(moments_x_gauss_nodes) ! temp remember order of the guass method for integration of moments .... 
  dx=min(minval(moments_x_gauss_nodes(2:mxgo) - moments_x_gauss_nodes(1:mxgo-1)),&
                 moments_x_gauss_nodes(1)+Real(1,DP), Real(1,DP)-moments_x_gauss_nodes(mxgo))/Real(3,DP)
  !! allocate memory for the nodes in x
     allocate (xmesh_l(1:mxgo*2+1), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
     stop
     end if
     !
!! allocate memory for the nodes in x
     allocate (xmesh_r(1:mxgo*2+1), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
     stop
     end if
     !
  !! Now we set up the mesh in x:
  xmesh_l(1)=x_L
  do i=1,mxgo
  xmesh_l(2*i) = (x_R+x_L)/Real(2,DP) + (x_R-x_L)/Real(2,DP)*(moments_x_gauss_nodes(i) - dx)
  xmesh_l(2*i+1) = (x_R+x_L)/Real(2,DP) + (x_R-x_L)/Real(2,DP)*(moments_x_gauss_nodes(i) + dx)
  end do 
  xmesh_r(1:2*mxgo)=xmesh_l(2:2*mxgo+1)
  xmesh_r(mxgo*2+1) = x_R
  N=mxgo*2+1 ! the new value of (N) in case somebody uses it directly... 
  !! mesh in x is ready
  case (2) ! This mesh is to be used with diffusive boundary conditions
           ! and will make 1/8-1/4-1/2 refinement near the walls
           ! Only Static wall is implemented at this time! 
           ! Requires N>=6 !!! 
  !! allocate memory for the nodes in x
  if (N<6) then 
    print *, "SetDGVblzmmesh: For the type 2 nonuniform mesh in x the number of mesh points has to be >=6 (xmesh_l)"
    stop  
  end if  
    allocate (xmesh_l(1:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
    stop
    end if
    !
    allocate (xmesh_r(1:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
    stop
    end if
    !
    dx=(x_R-x_L)/(Real(N-4,DP)-Real(1,DP)/Real(4,DP))
  !! Now we will set up the mesh in variable x:
    xmesh_l(1)=x_L
    xmesh_l(2)=x_L+dx/Real(8,DP)
    xmesh_l(3)=xmesh_l(2)+dx/Real(4,DP)
    xmesh_l(4)=xmesh_l(3)+dx/Real(2,DP)
     do i=5,N-3
      xmesh_l(i)=xmesh_l(i-1)+dx
     end do 
    xmesh_r(N)=x_R
    xmesh_l(N)=xmesh_r(N)-dx/Real(8,DP)
    xmesh_l(N-1)=xmesh_l(N)-dx/Real(4,DP)
    xmesh_l(N-2)=xmesh_l(N-1)-dx/Real(2,DP)
    xmesh_r(1:N-1)=xmesh_l(2:N)         
  !! end of the non-uniform mesh type 2. 
  case (3) ! This mesh is to be used with diffusive boundary conditions
           ! and will make 1/4-1/2 refinement near the walls
           ! Only Static wall is implemented at this time! 
           ! Requires N>=4 !!! 
  !! allocate memory for the nodes in x
  if (N<4) then 
    print *, "SetDGVblzmmesh: For the type 3 nonuniform mesh in x the number of mesh points has to be >=4 (xmesh_l)"
    stop  
  end if  
    allocate (xmesh_l(1:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
    stop
    end if
    !
    allocate (xmesh_r(1:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
    stop
    end if
    !
    dx=(x_R-x_L)/(Real(N-2,DP)-Real(1,DP)/Real(2,DP))
  !! Now we will set up the mesh in variable x:
    xmesh_l(1)=x_L
    xmesh_l(2)=x_L+dx/Real(4,DP)
    xmesh_l(3)=xmesh_l(2)+dx/Real(2,DP)
     do i=4,N-2
      xmesh_l(i)=xmesh_l(i-1)+dx
     end do 
    xmesh_r(N)=x_R
    xmesh_l(N)=xmesh_r(N)-dx/Real(4,DP)
    xmesh_l(N-1)=xmesh_l(N)-dx/Real(2,DP)
    xmesh_r(1:N-1)=xmesh_l(2:N)
  !!! end of the non-uniform mesh type 3. 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  case default 
    allocate (xmesh_l(1:3), stat=loc_alloc_stat)
     !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
    stop
    end if
     !
    allocate (xmesh_r(1:3), stat=loc_alloc_stat)
     !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
    stop
    end if
     ! 
    xmesh_l = (/ (x_L + (x_R-x_L)/Real(3,DP)*i, i=0,2) /)
    xmesh_r(1:2)=xmesh_l(2:3)
    xmesh_r(3) = x_R 
    N=3 ! the new value of (N) in case somebody uses it directly... 
    print *, "SetDGVblzmmesh: unsupported choice of non-uniform mesh. x_nonuniform_mesh_type= ", x_nonuniform_mesh_type
end select
 end if  
 ! now we will setup the level zero meshes in u ( grids_u/_v/_w and grids_cap_u/_v/_w )  
 ! first the mesh in u
  allocate (grids_u(1:Mu+1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_u)"
  stop
  end if
     !
  allocate (grids_cap_u(1:1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_cap_u)"
  stop
  end if
     !   
  call OldCreateUmesh(grids_u,grids_cap_u,mesh_u_uniform,u_nonuniform_mesh_type,u_L,u_R)
  ! then mesh in v 
  allocate (grids_v(1:Mv+1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_v)"
  stop
  end if
     !
  allocate (grids_cap_v(1:1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_cap_v)"
  stop
  end if
     !   
  call OldCreateUmesh(grids_v,grids_cap_v,mesh_v_uniform, v_nonuniform_mesh_type,v_L,v_R)
  ! now the mesh in w
  allocate (grids_w(1:Mw+1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_w)"
  stop
  end if
     !
  allocate (grids_cap_w(1:1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_cap_w)"
  stop
  end if
     !   
  call OldCreateUmesh(grids_w,grids_cap_w,mesh_w_uniform,w_nonuniform_mesh_type,w_L,w_R)
     !
 end subroutine SetDGVblzmmesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  Subroutine OldCreateUmesh 
!  This subroutine is used to generate 1D meshes in the velocity variables 
!  

subroutine OldCreateUmesh(ugrid,ugrid_cap,mesh_u_uniform,u_nonuniform_mesh_type,u_L,u_R) 
!
real (DP), dimension (:), intent (out)  :: ugrid
integer (I4B), dimension (:), intent (out)  :: ugrid_cap
logical, intent(in) :: mesh_u_uniform
integer (I4B), intent (in) :: u_nonuniform_mesh_type
real (DP), intent (out) :: u_L, u_R
! local variables: 
integer (I4B) :: M ! number of cells in the grid
integer (I4B) :: i,j ! local counters
!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: smr, sml ! refinement parameters for mesh of type 3
real (DP) :: du 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
M=size(ugrid,1)-1
ugrid_cap(1) = M+1
! 
if (mesh_u_uniform) then 
  ugrid = (/ (u_L + (u_R-u_L)/Real(M,DP)*i, i=0,M) /) 
  else
  select case (u_nonuniform_mesh_type)
  case (1) ! Unsupported  
     print *, "OldCreateUmesh: Unsupported type of mesh in velocity variable (u_nonuniform_mesh_type=1)"
     stop
  case (2) ! This mesh is to be used with diffusive boundary conditions
           ! and will include the velocity of the wall, u_{w}=0 as a mesh point
           ! Only Static wall is implemented at this time! 
           ! 
 ! first, let us check that u_wall =0  is in the velocity interval                     
  if ((u_L > 0.0_DP) .or. (u_R < 0.0_DP)) then 
    u_L = min (0.0_DP,u_L)
    u_R = max (0.0_DP,u_R)
  end if                    
  ! now we start to fill in the mesh:
  i=-M
  j=1
  ugrid(1)=u_L
  do while ( u_R - (u_R-u_L)/Real(M,DP)*i > (u_R-u_L)/Real(M,DP)/2 - 0.00000000001_DP )    
  if (i*(u_R-u_L)/Real(M,DP)-u_L > (u_R-u_L)/Real(M,DP)/2) then 
  j=j+1
  ugrid(j)=i*(u_R-u_L)/Real(M,DP)
  end if 
  i=i+1
  end do 
  ugrid(M+1)=u_R
  ! now we have filled in the mesh:
 case (3) ! This mesh is to be used with diffusive boundary conditions
           ! and will include the velocity of the wall, u_{w}=0 as a mesh point
           ! also it will surround $u=0$ with small cells. The sizes of the sell are obtained 
           ! by refining the size of biggest cells $smr^{sml}$ times. Here $smr$ -- is the 
           ! refinement factor and $sml$ is the total level of reinements. 
           ! Only Static wall is implemented at this time! 
           ! 
 ! First we set values for $smr$ and $sml$ 
 smr=10
 sml=2
 if (2*sml > M-2) then
  print *, "SetDGVblzmmesh: selected number of refinements 'sml' of mesh is too large for the given M"
  stop
 end if
   ! now we need to evaluate the size of the largest cell. 
 du=(u_R-u_L)/((M-2*sml) + 2.0_DP/Real(smr,DP)* &
             (1.0_DP/(Real(smr,DP)**sml)-1.0_DP)/(1.0_DP/Real(smr,DP)-1.0_DP))
 ! now we start to fill in the mesh:
  j=1
  ugrid(1) = u_L
  do while (ugrid(j) < - du*(1.0_DP/(Real(smr,DP)**sml) - 1.0_DP)/(1.0_DP/Real(smr,DP)-1.0_DP)/Real(smr,DP) - 0.00000000001_DP )    
  j=j+1
  ugrid(j)=ugrid(j-1)+du
  end do
  ! Now we fill in the small cells
  do i=1,sml
  j=j+1
  ugrid(j) = ugrid(j-1) + du/(Real(smr,DP)**i)
  end do 
  ugrid(j)=0.0_DP ! Just to have it exact --- otherwise it is not zero, but 10^-13
  do i=sml,1,-1
  j=j+1
  ugrid(j) = ugrid(j-1) + du/(Real(smr,DP)**i)
  end do 
  ! now the rest of the cells.
  do while (ugrid(j) < u_R - du - 0.00000000001_DP)
  j = j+1
  ugrid(j) = ugrid(j-1)+du
  end do
  ugrid(M+1) = u_R
 case default 
     print *, "SetDGVblzmmesh: Unknown type of mesh"
     stop
  end select
 end if  
end subroutine OldCreateUmesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set3DCellsR_DGVblzm
!
! This subroutine creates a 3D velocity mesh. 
!
! The subroutine looks up the grids arrays and creates 3D cells from the grids arrays. 
! If more than one grid is defined, grids are treated as non-related. Their ierarchy is ignored -- actually it is 
! not known from the grids array rigth away and needs to be discovered by a painful process. So we just do not do it. 
! instead, we would like to use this subroutine to create zero level mesh and we will hope that we only have one grid at this point.
! 
! Later, another subroutine may be called to create embedded grids. 
!
! Depends on the main program. Before calling Make sure su,sv,sw are defined!
! Also make sure that 1D grids in velocity are defined. 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Set3DCellsR_DGVblzm

use commvar, only: grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,su,sv,sw
!
integer (I4B) :: iu,iv,iw, ig ! some local counters

integer (I4B) :: mm,mx,igu,igv,igw ! local scrap counters,  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creating cells from whaterwer is already in the grids -- refinement will be applied later 
!

! we assume that at least one grid of level zero is defined...
if (size(grids_u,1) < 1) then 
  print *, "Set3DCellsR_DGVblzm: Error, there seem to be no 1D grids_u/_v/_w defined."
  stop 
end if 
mm=(grids_cap_u(1)-1)*(grids_cap_v(1)-1)*(grids_cap_w(1)-1)
if (mm < 1) then 
  print *, "Set3DCellsR_DGVblzm: Error, at least one of the 1D grids_u/_v/_w defined."
  stop 
end if 
! create the cells arrays... 
call AllocCellsArrsDGV(mm)
! now we will need to populate the cells: 
call FillCellsArrsDGV(1,mm, grids_u(1:grids_cap_u(1)),grids_v(1:grids_cap_v(1)),grids_w(1:grids_cap_w(1)),1,su,sv,sw)
! set some counters in case we have more stuff to process...
mx = mm  ! this will accumulate the total number of the recorded cells ...
igu=grids_cap_u(1)
igv=grids_cap_v(1)
igw=grids_cap_w(1)
! now if there is more grids... we needs to continue to create cells... 
do ig=2,size(grids_cap_u,1)   ! all grids_u/_v/_w arrays are of the same length 
 ! for each grid create cells...  
 mm=(grids_cap_u(ig)-1)*(grids_cap_v(ig)-1)*(grids_cap_w(ig)-1) !! this is how many new cells will be on this grid
 ! extend the cells arrays to fit more cells ...
 call ExtendCellsArrsDGV(mx+mm)
 ! now we will need to populate the cells: 
 call FillCellsArrsDGV(mx+1,mx+mm,grids_u(igu+1:igu+grids_cap_u(ig)),grids_v(igv+1:igv+grids_cap_v(ig)), & 
             grids_w(igw+1:igw+grids_cap_w(ig)),ig,su,sv,sw)
 ! set some counters in case we have more stuff to process...
 mx = mx+mm  ! this will accumulate the total number of the recorded cells ...
 igu=igu+grids_cap_u(ig)
 igv=igv+grids_cap_v(ig)
 igw=igw+grids_cap_w(ig)
 !
end do 

end subroutine Set3DCellsR_DGVblzm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! AllocCellsArrsDGV
!! This is a "macro" subrouting to save some space in another subroutine. 
!! given the sizes, it allocates the arrays and returns emply arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine AllocCellsArrsDGV(M)
use commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

integer (I4B), intent (in) :: M ! The size of the arrays...  
!
integer (I4B) :: loc_alloc_stat
!!!!!!!!!!!!!!
!
allocate (cells_pgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_pgrid)"
  stop
  end if
allocate (cells_cgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_cgrid)"
  stop
  end if
allocate (cells_lu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_lu)"
  stop
  end if
allocate (cells_lv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_lv)"
  stop
  end if
allocate (cells_lw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_lw)"
  stop
  end if
allocate (cells_ru(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_ru)"
  stop
  end if
allocate (cells_rv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_rv)"
  stop
  end if
allocate (cells_rw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_rw)"
  stop
  end if
allocate (cells_refu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_refu)"
  stop
  end if
allocate (cells_refv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_refv)"
  stop
  end if
allocate (cells_refw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_refw)"
  stop
  end if
allocate (cells_gou(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_gou)"
  stop
  end if
allocate (cells_gov(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_gov)"
  stop
  end if
allocate (cells_gow(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_gow)"
  stop
  end if
end subroutine AllocCellsArrsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ExtendCellsArrsDGV
!! This is a "macro" subroutine to save some space in another subroutine. 
!! given the sizes, this subroutine extends the cells arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtendCellsArrsDGV(M)
use commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

integer (I4B), intent (in) :: M ! The new size of the cells arrays...  
!
integer (I4B) :: loc_alloc_stat
integer (I4B) :: Md ! a crap variable to keed the old size 
!!!!!!!!!!!!!!
real (DP), dimension(:), allocatable ::cpg,ccg,clu,clv,clw,cru,crv,crw,&
                                       crefu,crefv,crefw,cgou,cgov,cgow  

Md=size(cells_pgrid, 1)
if (Md>M) then 
 print *, "ExtendCellsArrsDGV: new size for arrays cells is smaller then the old one... "
 stop
end if
!!!!
!!!! 
allocate (cpg(1:Md),ccg(1:Md),clu(1:Md),clv(1:Md),clw(1:Md),cru(1:Md),crv(1:Md),crw(1:Md), stat=loc_alloc_stat)
!
if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variables cpg,ccg,clu,clv,clw,cru,crv,crw"
  stop
  end if
allocate (crefu(1:Md),crefv(1:Md),crefw(1:Md),cgou(1:Md),cgov(1:Md),cgow(1:Md),stat=loc_alloc_stat)
!
if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variables cells_refw, cells_gow, cells_gou, cells_gov"
  stop
  end if
! save the cells arrays in the temp arrays... 
cpg=cells_pgrid
ccg=cells_cgrid
clu=cells_lu 
clv=cells_lv 
clw=cells_lw
cru=cells_ru
crv=cells_rv
crw=cells_rw
crefu=cells_refu
crefv=cells_refv
crefw=cells_refw 
cgow=cells_gow
cgou=cells_gou 
cgov=cells_gov
! deallocate the old cells arrays... 
deallocate (cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov)
! now we allocate with the new size M   
allocate (cells_pgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_pgrid)"
  stop
  end if
allocate (cells_cgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_cgrid)"
  stop
  end if
allocate (cells_lu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_lu)"
  stop
  end if
allocate (cells_lv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_lv)"
  stop
  end if
allocate (cells_lw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_lw)"
  stop
  end if
allocate (cells_ru(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_ru)"
  stop
  end if
allocate (cells_rv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_rv)"
  stop
  end if
allocate (cells_rw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_rw)"
  stop
  end if
allocate (cells_refu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_refu)"
  stop
  end if
allocate (cells_refv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV Allocation error for variable (cells_refv)"
  stop
  end if
allocate (cells_refw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_refw)"
  stop
  end if
allocate (cells_gou(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_gou)"
  stop
  end if
allocate (cells_gov(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_gov)"
  stop
  end if
allocate (cells_gow(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_gow)"
  stop
  end if
!! and the inverse assignment:
cells_pgrid(1:Md)=cpg
cells_cgrid(1:Md)=ccg
cells_lu(1:Md)=clu
cells_lv(1:Md)=clv
cells_lw(1:Md)=clw
cells_ru(1:Md)=cru
cells_rv(1:Md)=crv
cells_rw(1:Md)=crw
cells_refu(1:Md)=crefu
cells_refv(1:Md)=crefv
cells_refw(1:Md)=crefw
cells_gow(1:Md)=cgow
cells_gou(1:Md)=cgou
cells_gov(1:Md)=cgov
!! the rest will be filled outside. 
deallocate(cpg,ccg,clu,clv,clw,cru,crv,crw,crefu,crefv,crefw,cgou,cgov,cgow)
  
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine deallocates the cells arrays 
!

subroutine  DeAllocCellsArrsDGV
use commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

deallocate(cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov)

end subroutine DeAllocCellsArrsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FillCellsArrsDGV
! 
! This subroutine populates the certain portion of the cells_arrays.
! 
! Expects to have cells_arrays to exists and to be of a proper size
! 
! cells_cgrid = number of the grid. If the cell is not refined then the value is -1, if the cell is refined than 
!               cells_cgrid gives the number of the grid. 

subroutine FillCellsArrsDGV(ibeg,iend,umesh,vmesh,wmesh,pgrid,gou,gov,gow)

use commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

integer (I4B), intent(in) :: ibeg,iend ! beginning and end of the range to fill the cells array
real (DP), dimension (0:), intent (in) :: umesh,vmesh,wmesh ! 1D meshes to be used
integer (I4B), intent (in) :: pgrid ! the number of the parent grid --- grid where these cells belong. 
integer (I4B), intent (in) :: gou,gov,gow ! the max order of gauss lagrange basis to be assigned to the cells -- all cells in this range will be assigned the same order in all 
! 
integer (I4B) :: iu,iv,iw,mm ! local counters
! first some elementary checks
if ((size(cells_pgrid,1)< ibeg) .or. (size(cells_pgrid,1)< iend)) then 
 print *, "FillCellsArrsDGV: supplied range does not work for the cells arrays. "
 stop
end if
if ((ibeg > iend) .or. ((iend-ibeg + 1) .neqv. (size(vmesh,1)-1)*(size(wmesh,1)-1)*(size(umesh,1)-1)))   then 
 print *, "FillCellsArrsDGV: supplied range does not work: ibeg>iend or iend-ibeg + 1 .neq. size(vmesh)*size(wmesh)*size(umesh). "
 stop
end if  
! ok, let us populate the cells
! some properties will be the same for the entire range: 
cells_pgrid(ibeg:iend) = pgrid
cells_cgrid(ibeg:iend) = -1
cells_refu(ibeg:iend) = 1
cells_refv(ibeg:iend) = 1
cells_refw(ibeg:iend) = 1
cells_gou(ibeg:iend) = gou
cells_gov(ibeg:iend) = gov
cells_gow(ibeg:iend) = gow
! some properties will be different: 
mm=ibeg
do iu=0,size(umesh,1)-2
do iv=0,size(vmesh,1)-2
do iw=0,size(wmesh,1)-2
   cells_lu(mm) = umesh(iu)
   cells_lv(mm) = vmesh(iv)
   cells_lw(mm) = wmesh(iw)
   cells_ru(mm) = umesh(iu+1)
   cells_rv(mm) = vmesh(iv+1)
   cells_rw(mm) = wmesh(iw+1)
   mm=mm+1
end do 
end do 
end do 
end subroutine FillCellsArrsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! CellsRefineDGV
!! 
!! This subrotine refines cells that are on the supplied list
!!
!! Depends on the main program. The cells arrays must be set up before calling this subroutine
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CellsRefineDGV(reflist,refu,refv,refw)

use commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w

integer (I4B), dimension (:), intent (in) :: reflist  !list of cells to refine at this step
integer (I4B), intent(in) :: refu,refv,refw ! refinement factors in dimension u, v, and w
!!!!!!!!
!! IMPORTANT: 
!!
integer (I4B) :: go_inc=0 !! THIS IS A PARAMETER: the refined cells will have orders of the parent + go_inc
!!
!!!!!!!! 
integer (I4B) :: i,ic,mm ! local counters 
integer (I4B) :: isc ! the number of the cell that is being refined
integer (I4B) :: cgrid ! a scrap number to keep the number of the new grid
integer (I4B) :: gou,gov,gow ! the scrap variables to keep the gauss order of the basis
integer (I4B) :: loc_alloc_stat  ! local variable to keep the allocation status

real (DP), dimension (:), allocatable :: dumu,dumv,dumw 
! a quick check for errors: 
if ((refu*refv*refw < 1 ) .or. (size(reflist,1)<1)) then 
  print *,"CellsRefineDGV: error in incoming parameters, refu*refv*refw < 1 or size(reflist,1)<1"
end if 
! also, we need to have a shorthand for how many cells we will add
mm=(refu)*(refv)*(refw)
!
do i=1,size(reflist,1)
  isc = reflist(i)
  ! We need to construnct the new refined grids
  allocate (dumu(1:refu+1),dumv(1:refv+1),dumw(1:refw+1), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "CellsRefineDGV: Allocation error for variable (dumu),(dumv),or (dumw)"
  stop
  end if
     !
  dumu = (/ (cells_lu(isc) + (cells_ru(isc) - cells_lu(isc))/Real(refu,DP)*i, i=0,refu) /)
  dumv = (/ (cells_lv(isc) + (cells_rv(isc) - cells_lv(isc))/Real(refv,DP)*i, i=0,refv) /)
  dumw = (/ (cells_lw(isc) + (cells_rw(isc) - cells_lw(isc))/Real(refw,DP)*i, i=0,refw) /)
  ! Now we need to extend the grids arrays to include the new grids... 
  call ExtendGridsDGV(dumu,dumv,dumw,cgrid)
  ! Now we mark the cell with number (isc) to have a child grid...
  cells_cgrid(isc) = cgrid
  cells_refu(isc)= refu
  cells_refv(isc)= refv
  cells_refw(isc)= refw
  ! We also want to know what is the order of the Lagrange basis on that cell: 
  gou=cells_gou(isc)
  gov=cells_gov(isc)
  gow=cells_gow(isc)
  ! we need to know how long was cells before we added new: 
  ic=size(cells_lu,1)
  ! Now we need to extend the cells arrays to include the new cells... 
  call ExtendCellsArrsDGV(ic+mm)
  ! now we need to populate the new cells
  call FillCellsArrsDGV(ic+1,ic+mm,dumu,dumv,dumw,cgrid,gou+go_inc,gov+go_inc,gow+go_inc)
  ! ready to do the next cell 
end do 


end subroutine CellsRefineDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ExtendGridsDGV
!
! This subroutine extends the grids_cap_u/_v/_w and the grids_u/_v/_w arrays
! Depends on the main program. 
!
! dumu,dumv,dumw -- these are the new grids that need to be added to the old grids... 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtendGridsDGV(dumu,dumv,dumw,cgrid)

use commvar, only: grids_cap_u, grids_cap_v, grids_cap_w, grids_u, grids_v, grids_w

real (DP), dimension (:), intent (in) :: dumu,dumv,dumw ! the new velocity values that will be added to grids
integer (I4B), intent (out) :: cgrid ! the number of the child grid that will be created 
!
integer (I4B) :: iu,iv,iw ! some scrap variables... 
integer (I4B) :: loc_alloc_stat  ! local variable to keep the allocation status

real (DP), dimension (:), allocatable :: tmp_u,tmp_v,tmp_w ! temporary arrays to keep the old grids

! a quick check for errors
if ((size(dumu,1)< 3 ) .or. (size(dumv,1)< 3) .or. (size(dumw,1)< 3 )) then 
  print *,"ExtendGridsDGV: error in incoming parameters, at least one of the arrays does not require a refinement < 1"
end if 
! create the dump arrays to save current grids_cap_u/v/w:
iu=size(grids_cap_u,1)
allocate (tmp_u(1:iu),tmp_v(1:iu),tmp_w(1:iu), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendGridsDGV: Allocation error 1 for variable (tmp_u),(tmp_v),or (tmp_w)"
  stop
  end if
! save the old grids_cap arrays
tmp_u=grids_cap_u
tmp_v=grids_cap_v
tmp_w=grids_cap_w
! deallocate the grids_cap 
deallocate (grids_cap_u,grids_cap_v,grids_cap_w)
! allocate them 1 cell longer... 
allocate (grids_cap_u(1:iu+1),grids_cap_v(1:iu+1),grids_cap_w(1:iu+1), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendGridsDGV: Allocation error for variable (grids_cap_u/_v/_w)"
  stop
  end if
! restore the old grids_cap arrays
grids_cap_u(1:iu)=tmp_u
grids_cap_v(1:iu)=tmp_v
grids_cap_w(1:iu)=tmp_w
! record the new cells
grids_cap_u(iu+1)=size(dumu,1)
grids_cap_v(iu+1)=size(dumv,1)
grids_cap_w(iu+1)=size(dumw,1)
! done with grids_cap arrays
cgrid = iu+1 ! this will be passed out as the number of the new created grid
deallocate(tmp_u,tmp_v,tmp_w)
! now we extend the grids_u/_v/_w arrays... 
iu=size(grids_u,1)
iv=size(grids_v,1)
iw=size(grids_w,1)
allocate (tmp_u(1:iu),tmp_v(1:iv),tmp_w(1:iw), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendGridsDGV: Allocation error 2 for variable  (tmp_u),(tmp_v),or (tmp_w)"
  stop
  end if
! save the old grids_u/_v/_w arrays
tmp_u=grids_u
tmp_v=grids_v
tmp_w=grids_w
! deallocate the grids_cap 
deallocate(grids_u,grids_v,grids_w)
! allocate them again but longer... 
allocate (grids_u(1:iu+size(dumu,1)),grids_v(1:iv+size(dumv,1)),grids_w(1:iw+size(dumw,1)), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendGridsDGV: Allocation error for variable (grids_cap_u/_v/_w)"
  stop
  end if
! restore the old grids_cap arrays
grids_u(1:iu)=tmp_u
grids_v(1:iv)=tmp_v
grids_w(1:iw)=tmp_w
! record the new cells
grids_u(iu+1:iu+size(dumu,1)) = dumu
grids_v(iv+1:iv+size(dumv,1)) = dumv
grids_w(iw+1:iw+size(dumw,1)) = dumw
! done with grids_cap arrays
deallocate(tmp_u,tmp_v,tmp_w)
end subroutine ExtendGridsDGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine SetNodesDGV 
!
! This subroutine sets all nodes in velocity that will ever be used. 
!
! nodes_cells contains the # of the cell where this node belongs
! nodex_x/_y/_z contains the coordinates of the nodes.
! nodes_weights contains the product of the weights of gauss quadrature for each node so that the 
!               functions can be evaluated on the nodes and multiplied by the weigth and summed. 
!
! grids_cap_u -- stores the grid capasity 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetNodesDGV

use commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w, nodes_gwts,&
				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
                   g_wts_all, nodes_ui, nodes_vi, nodes_wi
                   
intrinsic SUM 

integer (I4B) :: ic,im,jm,iu,iv,iw ! some scrap variables... 
real (DP) :: ru,qu,rv,qv,rw,qw ! some scrap variables ...

integer (I4B) :: loc_alloc_stat  ! local variable to keep the allocation status

im=0
do ic=1,size(cells_pgrid,1)
 if (cells_cgrid(ic) .eq. -1) then 
  im=im+cells_gou(ic)*cells_gov(ic)*cells_gow(ic)
 end if
end do 
!
! We now allocating the nodes arrays
allocate (nodes_pcell(1:im),nodes_ui(1:im),nodes_vi(1:im),nodes_wi(1:im),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetNodesArDGV: Allocation error for variable (nodes_pcell,nodes_ui/_vi/_wi)"
  stop
  end if

allocate (nodes_u(1:im),nodes_v(1:im),nodes_w(1:im),nodes_gwts(1:im),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetNodesArDGV: Allocation error for variables (nodes_u/_v/_w/_gwts)"
  stop
  end if
! Now we will populate them. 
jm=1
do ic=1,size(cells_pgrid,1)
if (cells_cgrid(ic) .eq. -1) then 
 do iu=1,cells_gou(ic)
   ru=g_nds_all(iu,cells_gou(ic))*(cells_ru(ic)-cells_lu(ic))/2.0_DP + (cells_ru(ic)+cells_lu(ic))/2.0_DP
   qu=g_wts_all(iu,cells_gou(ic))
 do iv=1,cells_gov(ic)
   rv=g_nds_all(iv,cells_gov(ic))*(cells_rv(ic)-cells_lv(ic))/2.0_DP + (cells_rv(ic)+cells_lv(ic))/2.0_DP
   qv=g_wts_all(iv,cells_gov(ic))
 do iw=1,cells_gow(ic)
   rw=g_nds_all(iw,cells_gow(ic))*(cells_rw(ic)-cells_lw(ic))/2.0_DP + (cells_rw(ic)+cells_lw(ic))/2.0_DP
   qw=g_wts_all(iw,cells_gow(ic))
! ready to fill the nodes! 
nodes_pcell(jm) = ic   
nodes_u(jm)= ru 
nodes_v(jm)= rv
nodes_w(jm)= rw
nodes_gwts(jm)=qu*qv*qw*(cells_ru(ic)-cells_lu(ic))*(cells_rv(ic)-cells_lv(ic))*(cells_rw(ic)-cells_lw(ic))/8.0_DP
nodes_ui(jm)=iu
nodes_vi(jm)=iv
nodes_wi(jm)=iw
jm=jm+1
 end do 
 end do 
 end do 
end if  
end do 
!
if (jm-1 .ne. im) then 
  print *, "SetNodesArDGV: Error, jm-1 .neq. im. Nodes arrays may be populated wrong."
  stop
  end if

end subroutine SetNodesDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetCellsUGIAshiftArrays
!
!
! This subroutine is setting up the arrays cells_gui, _gvi, _gwi that 
! store the relative integer adress of the cell on the grid. 
!
! Only wil work if there is only one grid in the mesh.
!
! Also this subroutine sets up an arrays that gives the adress shift for different nodes to 
! read staff from the A-Array.
!
! Note that all these arrays can be created in other places of the algorithm, but we did nto think about them 
! earlier. So we put them here.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetCellsUGINdsSftArrs(ccell)

use commvar, only: nodes_pcell, nodes_ui, nodes_vi, nodes_wi, cells_ugi, cells_vgi, cells_wgi, &
             cells_pgrid, grids_cap_u, grids_cap_v, grids_cap_w, nodes_Ashift, nodes_phican,&
             cells_gou,cells_gov,cells_gow,nodes_dui,nodes_dvi,nodes_dwi,A_capphi
!!!!!!!!!!
integer (I4B), intent(in) :: ccell ! the number of the canonical cell -- a cell in the middle of the domain for which A is computed 
!!!!!!!!!!
integer (I4B) :: i,Ashift,phishift,my_cshift ! scrap variable...
integer (I4B) :: iphi ! the number of the basis function. Recall that each node generates a basis function.  
integer (I4B) :: phicell ! the number of the cell where phi belongs
integer (I4B) :: gou,gov,gow ! scrap variables to keep the number of nodes in a cell and the shift in nodes to the canonical cell
integer (I4B) :: pgcu,pgcv,pgcw ! scrap variables to keep the number of cells on the grid in each dimension
integer (I4B) :: ccell_ugi,ccell_vgi,ccell_wgi ! scrap variables to keep the integer address of the canonical cell on the grid
integer :: loc_alloc_stat ! scrap variable to keep the allocation status
!!!!!!!!!!!!!!!!!!!!
if (ccell < 0) then  ! a quick check for garbage in data: 
 print *,"EvalCollisionPeriodicA_DGV: Error. The provided value of (ccell) is negative or invalid. Stop."
end if 
!!!!
if (size(grids_cap_u,1) /= 1) then 
     print *, "SetCellsUGIAshiftArrays: Error. The number of grids must be 1. Stop"
     stop
end if  
!!! First, we need to calculate some useful constants:
gou=cells_gou(ccell)
gov=cells_gov(ccell)
gow=cells_gow(ccell)
my_cshift=(ccell-1)*gou*gov*gow ! this will be the number of the node right before the first node in canonical cell
! first we find the number of cells in each dimension
pgcu=grids_cap_u(cells_pgrid(ccell))-1
pgcv=grids_cap_v(cells_pgrid(ccell))-1
pgcw=grids_cap_w(cells_pgrid(ccell))-1
! now we will calculate the relative integer address of all cells on the grid.
! first, we create a space where to store this information: 
i=size(cells_pgrid,1)
allocate (cells_ugi(1:i),cells_vgi(1:i),cells_wgi(1:i),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetCellsUGIAshiftArrays: Allocation error for variable (cells_ugi/_vgi/_wgi)"
  stop
  end if
! Nex we will calculate the adresses
! IMPORTANT! Notice that this algorithm uses the convention that when enumerating cells the outside loop is in x, the 
! middle loop is in y and the innermost is in z.  

do phicell = 1,pgcu*pgcv*pgcw ! loop in all cells
!
ccell_ugi=1
ccell_vgi=1
ccell_wgi=1
i=1
do while (i<phicell)  
 i=i+1
 if (ccell_wgi == pgcw) then 
   ccell_wgi = 1 
   if (ccell_vgi == pgcv) then 
     ccell_vgi = 1
     if (ccell_ugi == pgcu) then  
      print *,"EvalCollisionPeriodicA_DGV: Error the index of the canonical cell is not fouind. Value of (gui) is too big."
      stop
     else 
      ccell_ugi=ccell_ugi+1
     end if
   else 
     ccell_vgi=ccell_vgi+1
   end if     
 else
   ccell_wgi = ccell_wgi+1
 end if           
end do 
! finshed finding the cells addresses, now recording:
cells_ugi(phicell) = ccell_ugi
cells_vgi(phicell) = ccell_vgi
cells_wgi(phicell) = ccell_wgi
end do 
! finished finding the local adress of the cells

! Next find the nodes_Ashift and nodes_dui/_dvi/_dwi arrays. 
! First we allocated them... 
i=size(nodes_vi,1)
allocate (nodes_Ashift(1:i),nodes_dui(1:i),nodes_dvi(1:i),nodes_dwi(1:i),nodes_phican(1:i),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetCellsUGIAshiftArrays: Allocation error for variable (nodes_Ashift,nodes_dui/dvi/dwi)"
  stop
  end if
! next we populate the cells_ugi/_vgi/wgi -arrays... 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iphi=1,size(nodes_vi,1)
    phicell = nodes_pcell(iphi)
    if (cells_pgrid(phicell) /= cells_pgrid(ccell)) then 
     print *, "SetCellsUGIAshiftArrays: Error. The parent cell of phi belongs to a different grid than the canonical cell. Stop"
     stop
    end if 
    ! next we record the shifts in integer index in from the cell containting phi to the canonical
    nodes_dui(iphi)=cells_ugi(phicell) - cells_ugi(ccell)
    nodes_dvi(iphi)=cells_vgi(phicell) - cells_vgi(ccell)
    nodes_dwi(iphi)=cells_wgi(phicell) - cells_wgi(ccell)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now we need to calculate the nodes_Ashift number and the nodes_phican number 
    ! First we calculate the number of the node on the canonical cell that corresponds to the node with number iphi 
    phishift = my_cshift + (nodes_ui(iphi)-1)*gov*gow + (nodes_vi(iphi)-1)*gow + nodes_wi(iphi) ! this is the number of the node
    nodes_phican(iphi) = phishift ! we will record this number. 
    ! Now we will calculate the Ashift number
    Ashift = sum(A_capphi(1:phishift-1))
    nodes_Ashift(iphi)=Ashift
    !!!!! all done
end do 
end subroutine SetCellsUGINdsSftArrs


end module miscset
