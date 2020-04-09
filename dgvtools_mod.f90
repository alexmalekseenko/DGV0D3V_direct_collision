!
! dgv_mod.f90
! 
! This module contains subroutines for working with Discontinuous Galerkin velocity discretization. 
! It contains definitions and routines for working with the basis functions 
! 
! Subroutines of this module complement the subroutines of setting up grids and cells and nodes in miscset
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dgvtools_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none


 interface lagrbasfun
     module procedure lagrbasfun, lagrbasfun_xvec
   end interface

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! lagrbasfun
!
! This subroutine evaluates the lagrange basis function given the x and the 
! array of nodes kappa
! 
! \Prod_{j\neq i} \frac{kappa(j)-x}{kappa(j)-kappa(i)}
!

function lagrbasfun(i,x,kappa) result (y)

integer (I4B) :: i ! the number of the node where the lagrange basis function is one. 
real (DP) :: x ! the point where the function needs to be evaluated
real (DP), dimension(:) :: kappa ! the nodes of the lagrange basis functions
real (DP) :: y ! the value of the function 
!
integer (I4B) :: j ! some local counter


!!!!!!!!!!!!!1
! a quick consistancy check
if ((i<1) .or. (i>size(kappa,1))) then 
 print *," lagrbasfun: error. (i<1) .or. (i>size(kappa,1)) no Lagrange basis function with this number."
 stop
endif  
!
y=1.0_DP
do j=1,i-1
y = y*(kappa(j)-x)/(kappa(j)-kappa(i))
enddo 
do j=i+1,size(kappa,1)
y = y*(kappa(j)-x)/(kappa(j)-kappa(i))
enddo 

end function lagrbasfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! lagrbasfun
!
! This subroutine evaluates the lagrange basis function given the x and the 
! array of nodes kappa
! 
! \Prod_{j\neq i} \frac{kappa(j)-x}{kappa(j)-kappa(i)}
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function lagrbasfun_xvec(i,x,kappa) result (y)

integer (I4B) :: i ! the number of the node where the lagrange basis function is one. 
real (DP), dimension(:) :: x ! the point where the function needs to be evaluated
real (DP), dimension(:) :: kappa ! the nodes of the lagrange basis functions
real (DP), dimension(1:size(x,1)) :: y ! the value of the function 
!
integer (I4B) :: j ! some local counter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a quick consistancy check
if ((i<1) .or. (i>size(kappa,1))) then 
 print *," lagrbasfun: error. (i<1) .or. (i>size(kappa,1)) no Lagrange basis function with this number."
 stop
endif  
!
y=1.0_DP
do j=1,i-1
y = y*(kappa(j)-x)/(kappa(j)-kappa(i))
enddo 
do j=i+1,size(kappa,1)
y = y*(kappa(j)-x)/(kappa(j)-kappa(i))
enddo 

end function lagrbasfun_xvec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! EvalLagrBasisFunByNdsDGblzm
! 
! This subroutine evaluates the basis function identified by a node 
! for a given a value of velocity. To do this, first the cell is found 
! to which this velocity belongs. Then the cell is found to which the basis function belongs and 
! the numbers of the 1D-basis functions are identified. Then if the basis function is defined on this 
! cell, then the basis function is evaluated using the function lagrbasfun. If the velocity and the 
! basis function belong to different cells then the value of the basis function is zero
!
! This function directly uses arrays from commvar.f90. Specifically, grids, cells and nodes.
!
! u,v,w = components of the velocity where the basis function needs to be evaluated. 
! 
! i = the number of the node associated with this basis function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!

function EvalLagrBasisFunByNdsDGblzm(u,v,w,i) result (y)

use commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w,&
				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
                   nodes_ui, nodes_vi, nodes_wi

real (DP) :: u,v,w ! the components of the given velocity
integer (I4B) :: i  ! the number of the node in the nodes array. this the number of the node associated with 
                    ! the given basis function. Thus basis functions may be numbered using nodes... 
real (DP) ::  y ! the value of the basis function on the given velocity                    

!!!!!!!!!!!!!!                      
integer (I4B) :: j ! a counter -- usually denotes the number of the grid we are working on. 
integer (I4B) :: gi_v,gi_u,gi_w ! counters to help count nodes in the grids array
integer (I4B) :: celli, cellp ! these will store the number of the cell where the velocity belongs and the number of the 
            ! parent cell where the basis function belongs. 
integer (I4B) :: ui,vi,wi ! are the local numbers on the particular grids where u belongs. Because we have ierarchicval grids, 
             ! there may be different grids that have ui,vi,wi defined. If one of these numbers is zero -- then the 
             ! velocity u,v,w is not on the grid. Unless a mistake was done in setting the grids either all three of them 
             ! are zero -- that means that the velocity is not on this (level zero grid) or all three are non-zero  -- means that
             ! the velocity has been found. Mized results may indicate a breach in data integrity. 
             !!!!!
integer (I4B) :: jj ! a counter
real (DP) :: unor,vnor,wnor ! scrap variables to store the normalized velus of the components u,v,w             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! first, given a velocity, we need to identify the cell where it came from and its 1D numbers 
! we start looking in course cells first. If a cell is found that has that velocity we check if the cell is refined. 
! if it is, then we look into the corresponding grid and seek the smaller cell that has it and so on. 
!!!!!!!!
!
j=1
!! set up the shift in grids_u/_v/_w corresponding to the first grid (zero shift):
gi_u=0
gi_v=0
gi_w=0
do while (j<=size(grids_cap_u,1))
  call finduiviwi(grids_u(gi_u+1:gi_u+grids_cap_u(j)),grids_v(gi_v+1:gi_v+grids_cap_v(j)),&
                       grids_w(gi_w+1:gi_w+grids_cap_w(j)),u,v,w,ui,vi,wi)
  if (ui*vi*wi /= 0) then 
     ! now we need to find the cell that correspond to this grid and these ui,vi,wi
     celli=1
     do jj=1,j-1
      celli = celli + (grids_cap_u(jj)-1)*(grids_cap_v(jj)-1)*(grids_cap_w(jj)-1)
     end do
     celli=celli + (ui-2)*(grids_cap_v(j)-1)*(grids_cap_w(j)-1)+(vi-2)*(grids_cap_w(j)-1)+(wi-2)
     !  check if this cell is refined
     if ((cells_refu(celli)>1) .or. (cells_refv(celli)>1) .or. (cells_refw(celli)>1)) then  
        j=cells_cgrid(celli)
        ! set up the shift in grids_u/_v/_w corresponding to the j'th grid:
        gi_u=0
        gi_v=0
        gi_w=0
        do jj=1,j-1
         gi_u = gi_u + grids_cap_u(jj)
         gi_v = gi_v + grids_cap_v(jj)
         gi_w = gi_w + grids_cap_w(jj)
        end do
     else 
        exit
     endif
  else 
     celli=0
     exit
  endif                      
enddo
! now "celli" either is the number of the cell or 0 (0 means that the velocity in not on any cell)
! next step is to find the numbers that will help us compute the value of the basis function
! first we need to know to what cell in u,v,w this basis function belongs. This is simple becasue this 
! information is stored in the nodes arrays -- nodes_pcell
cellp=nodes_pcell(i)
! A quick check, if the velocity is not on the cell than the value of the basis function is zero and we are all done
if ((cellp /= celli) .or. (celli == 0)) then 
y=0.0_DP
                    else 
! if the velocity is on the cell, then we need to compute the value of $y$:
y=1.0_DP    
! next we need to know the three local indices that tell what velocity nodal values correspond to this 
! basis function. this is also simple since this information is also stored in the Nodes Arrays.
unor = ( u - (cells_ru(cellp) + cells_lu(cellp))/2.0_DP )/(cells_ru(cellp) - cells_lu(cellp))*2.0_DP 
y=y*lagrbasfun(nodes_ui(i),unor,g_nds_all(:cells_gou(cellp),cells_gou(cellp)))
!
vnor = ( v - (cells_rv(cellp) + cells_lv(cellp))/2.0_DP )/(cells_rv(cellp) - cells_lv(cellp))*2.0_DP 
y=y*lagrbasfun(nodes_vi(i),vnor,g_nds_all(:cells_gov(cellp),cells_gov(cellp)))
!
wnor = ( w - (cells_rw(cellp) + cells_lw(cellp))/2.0_DP )/(cells_rw(cellp) - cells_lw(cellp))*2.0_DP 
y=y*lagrbasfun(nodes_wi(i),wnor,g_nds_all(:cells_gow(cellp),cells_gow(cellp)))
!   
!!!!!!!!!!!!
end if
 
end function EvalLagrBasisFunByNdsDGblzm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! EvalLagrBasisFunByNdsDGblzmEZ
! 
! This subroutine evaluates the basis function identified by a node 
! for a given a value of velocity. To do this, first the cell is found 
! to which this velocity belongs. We check if the velocity falls within the 
! cell where the basis function is defined. The 
! numbers of the 1D-basis functions are identified. Then the basis function is 
! evaluated using the function lagrbasfun. If the velocity and the 
! basis function belong to different cells then the value of the basis function is zero
!
! This function directly uses arrays from commvar.f90. Specifically, grids, cells and nodes.
!
! u,v,w = components of the velocity where the basis function needs to be evaluated. 
! 
! i = the number of the node associated with this basis function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!

function EvalLagrBasisFunByNdsDGblzmEZ(u,v,w,i) result (y)

use commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w,&
				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
                   nodes_ui, nodes_vi, nodes_wi

real (DP) :: u,v,w ! the components of the given velocity
integer (I4B) :: i  ! the number of the node in the nodes array. this the number of the node associated with 
                    ! the given basis function. Thus basis functions may be numbered using nodes... 
real (DP) ::  y ! the value of the basis function on the given velocity                    

!!!!!!!!!!!!!!                      
integer (I4B) :: cellp ! these will store the number of the cell where the velocity belongs and the number of the 
            ! parent cell where the basis function belongs. 
real (DP) :: lu,ru,lv,rv,lw,rw,unor,vnor,wnor ! scrap variables to store the normalized velus of the components u,v,w             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! first, given the number of the node/basis function, we know the cell where the function belongs:
!!!(rem to restore)!!!cellp= nodes_pcell(i) ! this is the number of the nodes parent cell. 
! Now we check is the velocity is within the cell bounds: 
!!!(rem to restore)!!!lu=cells_lu(cellp); ru=cells_ru(cellp) 
!!!(rem to restore)!!!lv=cells_lv(cellp); rv=cells_rv(cellp) 
!!!(rem to restore)!!!lw=cells_lw(cellp); rw=cells_rw(cellp) 
!!!(rem to restore)!!!if ((u<lu) .or. (ru<u) .or. (v<lv) .or. (rv<v) .or. (w<lw) .or. (rw<w)) then 
!!!(rem to restore)!!!          y=0.0_DP
!!!(rem to restore)!!!else 
! if we are here then the velocity belongs to the cell where the basis function is defined:
! next step is to find the numbers that will help us compute the value of the basis function
! first we need to know to what cell in u,v,w this basis function belongs. This is simple becasue this 
! information is stored in the nodes arrays -- nodes_pcell (see above)

!!! Uncomment this line if want the zeros's moment 
!y=1.0_DP 
!!! Uncomment one of the three lines for the first moment:
!y = u
!y = v
!y = w 
!!! Uncomment the next line to get the second moment
y=u*u + v*v + w*w  
! next we need to know the three local indices that tell what velocity nodal values correspond to this 
! basis function. this is also simple since this information is also stored in the Nodes Arrays.
!!!(rem to restore)!!!unor = ( u - (ru+lu)/2.0_DP )/(ru - lu)*2.0_DP 
!!!(rem to restore)!!!y=y*lagrbasfun(nodes_ui(i),unor,g_nds_all(:cells_gou(cellp),cells_gou(cellp)))
!
!!!(rem to restore)!!!vnor = ( v - (rv+lv)/2.0_DP )/(rv - lv)*2.0_DP 
!!!(rem to restore)!!!y=y*lagrbasfun(nodes_vi(i),vnor,g_nds_all(:cells_gov(cellp),cells_gov(cellp)))
!
!!!(rem to restore)!!!wnor = ( w - (rw+lw)/2.0_DP )/(rw - lw)*2.0_DP 
!!!(rem to restore)!!!y=y*lagrbasfun(nodes_wi(i),wnor,g_nds_all(:cells_gow(cellp),cells_gow(cellp)))
!   
!!!!!!!!!!!!
!!!(rem to restore)!!!end if

end function EvalLagrBasisFunByNdsDGblzmEZ


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! finduiviwi(ugrid,vgrid,wgrid,u,v,w,ui,vi,wi)
! 
! This subroutine lookes trhough three one-dimensional meshes and checks whether u is in the mesh ugrid, 
! v is in the mesh vgrid, w is in the mesh wgrid. If it is, it returns the numbe greater than or equal to 1
! if a zero is returned, the number is not on the mesh
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine finduiviwi(ugrid,vgrid,wgrid,u,v,w,ui,vi,wi)

real (DP), dimension (:), intent (in) :: ugrid,vgrid,wgrid ! the grids in the three variables..
real (DP), intent (in) :: u,v,w ! the components of the velocity in question
integer (I4B), intent (out) :: ui,vi,wi  !  the coordinates of the right endpoint in each dimension. if =0  then the point is not there.. 
!!
real (DP) :: epsill=10d-11 ! a small parameter to help curb the effects of the round off error..  
!!

if ((u < ugrid(1)-epsill) .or. (u > ugrid(size(ugrid,1)) + epsill)) then 
  ui=0
else
  ui=1
  do while ((u >=ugrid(ui)-epsill) .and. (ui<size(ugrid,1)))
  ui=ui+1
  enddo
endif
!
if ((v < vgrid(1)-epsill) .or. (v > vgrid(size(vgrid,1)) + epsill)) then 
  vi=0
else
  vi=1
  do while ((v >= vgrid(vi)-epsill) .and. (vi<size(ugrid,1)))
  vi=vi+1
  enddo
endif
!
if ((w < wgrid(1)-epsill) .or. (w > wgrid(size(wgrid,1)) + epsill)) then 
  wi=0
else
  wi=1
  do while ((w >= wgrid(wi)-epsill) .and. (wi<size(wgrid,1)))
  wi=wi+1
  enddo
endif
end subroutine finduiviwi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! SetA_DGV
!
! This subroutine sets up the collision information operator A(\xi,\xi_1,\varphi^{j}_{p})
! This is a symmetric operator in $\xi$ and $\xi_{1}$, therefore we only are interested in the 
! records for unique unordered pairs (\xi,\xi_{1}). Because all velocities are indexed by one index, say i, 
! we only interested in producing records for pairs $(\xi_{i},\xi_{j})$, $i>j$. Notice that the case \xi=\xi_{i} 
! gives value 0.
!   
! This subroutine depends on the main program
! 
! You should consult the notes n101011.tex for detail on the formulas. 
!
! Evaluation of each entry of A involves a two dimensional integral. If the entry is smaller than some given number, 
! the entry is neglected/nullified by force
! 
! The structure of the arrays related to the operator A: 
! A_capphi(i) gives the number of non-zero entries in A(\xi,\xi_{1},\varphi^{j}_{p}) for each basis function with number i
! for each non-zero entry of A this one keeps the index of velocities that produced that non-zero entries. 
! A(i) i -- is the index of nonzero entires, we need to use other A-arrrays to restore 
! what velocities and wnat basis function this index correcpods  
! for example, A_sind_xi(i) gives the first velocity, A_sind_xi1 gives the second velcity
! and A_phi gives the index of the used basis function. 
!                                    
! A(\xi,\xi_{1};\varphi^{j}_{p})=\frac{d^2 |g|}{8}  \int_{0}^{\pi} \int_{0}^{2\pi}
!(\varphi^{j}_{p}(\xi')+\varphi^{j}_{p}(\xi'_{1}))
! d\varepsilon\, \sin \chi d \chi - \frac{d^2 |g| \pi }{2} (\varphi^{j}_{p}(\xi)+\varphi^{j}_{p}(\xi_{1})) 
!
! Takes 
! Trad  === This variable determines the cut off radius for $A$. It should be based on the estimated size of 
!                               ! non-trivial support for the distribution function
! ErrChi == the parameter defining how accurate the evaluation of the integral in Chi  should be. 
! ErrEps == the parameter defining how accurate the evaluation of the integral in Epsilon should be. 
! min_sens == minimum accepted value for A, values below will be neglected
! I1_list == List of basis functions for which A-array is evaluated.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetA_DGV
use commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w, nodes_ui, nodes_vi, nodes_wi, &
                   nodes_gwts,A_capphi,A_xi,A_xi1,A_phi,A,&
                   cells_lu, cells_lv, cells_lw, cells_ru, cells_rv, cells_rw,&
                   Trad, ErrChi, ErrEps, min_sens, I1_list,Num_OMP_threads

use gaussian_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), parameter :: pi25DT = 3.141592653589793238462643d0
integer (I4B), parameter :: NdsChi = 10000 ! The max number of cells for integration in chi and epsilon
real (DP) :: ires_e=100.0_DP ! this parameter will determine the resolution for integration in \varepsilon
real (DP) :: ires_t=100.0_DP ! this parameter will determine the resolution in t.
                                   ! the number of points is dictated by the radius of the collision sphere, which is |g|/2
                                   ! the bigger is the radius, the more points is needed to integrate over the sphere.
                                   ! these parameters will represent a lengh of the arc that is desired for application of 
                                   ! a gauss quadrature rule. Then the angular intervals will be broken in portions no bigger than 
                                   ! ires_t(_e)/|g| to ensure sufficient angular resolution
!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (2*NdsChi) :: NodesChi1, NodesChi2 
real (DP), dimension (3*NdsChi) :: FuncChi1, FuncChi2 ! Arrays to store cells and functin values on cells Two copies are needed. Integration in Chi
integer (I4B), dimension (NdsChi) :: CellsChiRef1,CellsChiRef2 ! Arrays to keep refinement flags
!!
integer (I4B) :: A_ct, phi_ct, loc_alloc_stat ! scrap variables. A_ct is a counter of how many records have been created for A
							!phi_ct counts how many records have been created for this particular basis function.  
integer (I4B) :: An ! this varable will keep the current size of the array A.
integer (I4B) :: ni ! this one will keep the size of the nodes array
integer (I4B) :: i1,i2,i3,d,e ! local counters  
integer (I4B) :: pci ! local integer
real (DP), dimension (8) :: uu,vv,ww ! coordinates of the 8 vertices of the support of the basis function
real (DP) :: dphi, dsph ! diameters for the circumscribed sphere for basis function and the collision shpere
real (DP) :: xiu,xiv,xiw,xi1u,xi1v,xi1w,xiupr,xivpr,&
             xiwpr,xi1upr,xi1vpr,xi1wpr,ugu,ugv,ugw ! scrap variables to keep the values of the pre and post collision velocities 
real (DP) :: ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3 ! more screap variables
real (DP) :: ku,kv,kw ! coordinates of the post-coll g'
real (DP) :: dist2, g1,g2 ! useful variables 
integer (I4B)  :: Nchi1 ! numbers of subdivisions in epsilon and chi
real (DP) :: Atemp,quad,quad1,quad2,es,my_int1,varphi_xi,varphi_xi1 ! scrap variables.
logical :: chiadaptflag
!
integer :: iiii ! a test variable to play with OpenMP runtime functions calls
!!!!!!!!!!!!!!!!!!!!!!!!!! Interface for OpenMP runtime libraries !!!!!!!!!!!!!!!!!!!
!interface 
! function omp_get_thread_num() result (y)
!   integer :: y 
! end function omp_get_thread_num
! function omp_get_num_threads() result (y)
!  integer :: y 
! end function omp_get_num_threads 
! function omp_get_num_procs() result (y)
!  integer :: y 
! end function omp_get_num_procs
! function omp_get_stack_size () result (y)
!  use nrtype
!  integer (I2B) :: y
! end function omp_get_stack_size
!end interface  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initial Allocation of the A-arrays 
ni=size(nodes_pcell,1)
allocate (A_capphi(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetA_DGV: Allocation error for variable  A_capphi"
  stop
  end if
An=10*ni
allocate (A_xi(1:An), A_xi1(1:An), A(1:An), A_phi(1:An), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetA_DGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
A_xi=0; A_xi1=0; A=0; A_phi=0; A_capphi=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
A_ct=0 ! in the beginning there is zero records
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!TEMP: do it just for one velocity nodes:  do i1=1,ni  ! loop in basis functions 
i1=I1_list(1)
!!
phi_ct=0 ! in the beginning there is zero records
pci=nodes_pcell(i1)
!
uu(1)=cells_lu(pci); vv(1)=cells_lv(pci); ww(1)=cells_lw(pci) ! the first one holds the lower left left corner
uu(2)=cells_ru(pci); vv(2)=cells_rv(pci); ww(2)=cells_rw(pci) ! the second holds the upper right right corner
uu(3)=cells_lu(pci); vv(3)=cells_lv(pci); ww(3)=cells_rw(pci)
uu(4)=cells_lu(pci); vv(4)=cells_rv(pci); ww(4)=cells_lw(pci)
uu(5)=cells_lu(pci); vv(5)=cells_rv(pci); ww(5)=cells_rw(pci)
uu(6)=cells_ru(pci); vv(6)=cells_lv(pci); ww(6)=cells_lw(pci)
uu(7)=cells_ru(pci); vv(7)=cells_lv(pci); ww(7)=cells_rw(pci)
uu(8)=cells_ru(pci); vv(8)=cells_rv(pci); ww(8)=cells_lw(pci)
dphi = sqrt((uu(2)-uu(1))**2+(vv(2)-vv(1))**2+(ww(2)-ww(1))**2)
!! we set the resolution parameters based on the size of the velocity cell.
ires_t = max(uu(2)-uu(1),vv(2)-vv(1),ww(2)-ww(1))/32
ires_e = ires_t
!!!
! Add openMP parallel directives here....  
!!!
! iiii=omp_get_num_procs
! OpenMP set the number of threads: 

!call omp_set_num_threads(Num_OMP_threads)
!$OMP PARALLEL DO PRIVATE(xiu,xiv,xiw,xi1u,xi1v,xi1w,dsph,dist2,ugu,ugv,ugw, & 
!$OMP    ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,Nchi1,NodesChi1,FuncChi1,d, & 
!$OMP    CellsChiRef1,Atemp,ChiAdaptFlag,my_int1,NodesChi2,FuncChi2,CellsChiRef2,e, &
!$OMP    quad,quad1,quad2,es,i2,i3,iiii,varphi_xi,varphi_xi1) NUM_THREADS(Num_OMP_threads) &
!$OMP    SCHEDULE(DYNAMIC, 5)  
do i2=1,ni  ! loop in velocity 1
!!!!
   xiu=nodes_u(i2); xiv=nodes_v(i2); xiw=nodes_w(i2)
   varphi_xi = EvalLagrBasisFunByNdsDGblzmEZ(xiu,xiv,xiw,i1) ! this is the value of the basis function on the first velocity -- will be passed to integrator
do i3=i2+1,ni ! loop in velocity 2
   xi1u=nodes_u(i3); xi1v=nodes_v(i3); xi1w=nodes_w(i3)
   varphi_xi1 = EvalLagrBasisFunByNdsDGblzmEZ(xi1u,xi1v,xi1w,i1) ! this is the value of the basis function on the second velocity -- -- will be passed to integrator
!!!!!!!!!!!! evaluation of A !!!!!!!!!!!!!!!
 ! evaluate the diamieter of the collision sphere and the double distance from the center of the basis function support to the center of the collision sphere
   dsph = sqrt((xiu-xi1u)**2+(xiv-xi1v)**2+(xiw-xi1w)**2)
   dist2 = sqrt((xiu+xi1u-uu(1)-uu(2))**2+(xiv+xi1v-vv(1)-vv(2))**2+(xiw+xi1w-ww(1)-ww(2))**2)
 ! Now we calculate some useful quantities:
   ugu = (xiu-xi1u)/dsph  ! need to introduce a unit vector in the direction of g
   ugv = (xiv-xi1v)/dsph  
   ugw = (xiw-xi1w)/dsph
   ! The components of the vector (ugu,ugv,ugw) should be orthogonal to vecotor g. 
   ! there is only one case when these componets are degenerate. is it when g is parallel to (1,1,1)
   ! in this case, we need to cook a non-trivial vecotor. The next if statement takes case of that: 
   if (abs(ugw-ugv)+abs(ugu-ugw)+abs(ugv-ugu) < 1.0d-6) then 
   ! first, we create two orthogonal non-trivial vectors
   ugux2 = ugv 
   ugvx2 = -ugu 
   ugwx2 = 0
   ugux3 = ugu*ugw
   ugvx3 = ugv*ugw
   ugwx3 = -(ugu)**2-(ugv)**2
   ! 
   g1=sqrt(ugux2**2+ugvx2**2)
   g2=sqrt(ugux3**2+ugvx3**2+ugwx3**2)
   !
   else 
   !first we create the components of two non-trivial orthogonal vectors: 
   ugux2 = (ugw-ugv)
   ugvx2 = (ugu-ugw)
   ugwx2 = (ugv-ugu)
   ugux3 = ((ugv)**2-(ugv)*(ugu)-(ugw)*(ugu)+(ugw)**2)
   ugvx3 = ((ugw)**2-(ugw)*(ugv)-(ugu)*(ugv)+(ugu)**2)
   ugwx3 = ((ugu)**2-(ugu)*(ugw)-(ugv)*(ugw)+(ugv)**2)
   !  
   g1 = sqrt(ugux2**2 + ugvx2**2 + ugwx2**2)
   g2 = sqrt(ugux3**2+ugvx3**2+ugwx3**2)
   end if
!  First we estimate if A is zero by cheking the overlap of spheres
   if ((dist2 > dsph+dphi) .or. (dist2+dphi < dsph) .or. (dsph > Trad)) then 
        cycle ! collision shpere does not hit the support of the basis function continue with the next velocity
              ! or the collision shpere is so large that it is not possible for distribution function to be nonzero for both velocities 
   end if 
 ! if we got here then the collision sphere has some possible overlap with the support of basis function
 ! To evaluate the integral, we will implement adaptive quadrature in both directions using Simpson's rule. 
 !
 ! begin integration in \chi
 ! first, we set up the initial mesh
   Nchi1=FLOOR(pi25DT*dsph/ires_t/2.0_DP)+1
   if (Nchi1>NdsChi) then 
     Nchi1=NdsChi
     print *,"SetA_DGV: Warning! Number of nodes (NdsChi) for integration in epsilon gives insufficient resolution"  
   end if 
   ! we set the initial mesh
   do d=1,Nchi1
    NodesChi1(2*(d-1)+1)=(d-1)*pi25DT/Real(Nchi1,DP)
    NodesChi1(2*d)=d*pi25DT/Real(Nchi1,DP)
   end do
   ! now we evaluate the integrand on the initial mesh
   FuncChi1=0
   do d=1,Nchi1
   ! each d gives one cell
   ! First, we evaluate the integrand on the left node of the cell
    if ((d>1) .and. (NodesChi1(2*d-1) == NodesChi1(2*d-2))) then 
     FuncChi1(3*d-2)= FuncChi1(3*d-3) ! if the node is repeating, the value has been computed already
    else 
     FuncChi1(3*d-2) = A_IntEpsilon(NodesChi1(2*d-1),ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
    end if
   ! This takes care of the left node in the cell... 
   !  
   !Now we evaluate the solution on the right node of the cell and at the center. 
   FuncChi1(3*d) = A_IntEpsilon(NodesChi1(2*d),ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
   FuncChi1(3*d-1) = A_IntEpsilon((NodesChi1(2*d-1) + NodesChi1(2*d))/2.0_DP,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,&
              ugu,ugv,ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
   ! repeat for all cells... 
   end do      
   ! finally, we mark all cells for the refinement in the beginning... 
   CellsChiRef1=0 
   CellsChiRef1(1:Nchi1)=1
   ! 
   !Next we adaptively refine and evaluate the integral:    
   Atemp=0
   ChiAdaptFlag =.true.
   my_int1 = 100*ErrChi ! initialize with something big for the intial iteration. 
   do while (ChiAdaptFlag)
   ChiAdaptFlag =.false.
   ! Now we will go over the array of the integrand values and will cross 
   ! out cells where the integrand is zero, and those that are not marked for the refinment 
   ! on those cells we summ the integrand and keep the sum.
   ! Cells marked for the refinement we divide into halves 
   NodesChi2=0
   FuncChi2=0
   CellsChiRef2=0
   e=0
   if (my_int1>ErrChi) then  ! this checks that the integral over the mesh NodesChi1 will be small anyway
   my_int1=0 ! this will keep the estimate of the integral on the mesh NodesChi1 so that we stop resolving when this number is small
   do d=1,Nchi1
   if ((abs(FuncChi1(3*d-2))+abs(FuncChi1(3*d-1))+abs(FuncChi1(3*d)) > 1.0d-15) .and. (CellsChiRef1(d)==1)) then 
    if (e > NdsChi-2) then 
     print *,"SetA_DGV: Number of nodes in chi is too big (e>NdsChi-2)"
     stop
    end if
    ! save the nodes and the integrand values
    NodesChi2(2*e+1)=NodesChi1(2*d-1)
    FuncChi2(3*e+1)=FuncChi1(3*d-2)
    NodesChi2(2*e+4)=NodesChi1(2*d)
    FuncChi2(3*e+6)=FuncChi1(3*d)
    ! Evaluate the midpoint node: 
    NodesChi2(2*e+2)=(NodesChi1(2*d-1)+NodesChi1(2*d))/2
    NodesChi2(2*e+3)=NodesChi2(2*e+2)
    ! Save the midpoint value of the integrand:
    FuncChi2(3*e+3)=FuncChi1(3*d-1)
    FuncChi2(3*e+4)=FuncChi2(3*e+3)
    ! Evaluate the integrand on the new nodes
    FuncChi2(3*e+2) = A_IntEpsilon((NodesChi2(2*e+2)+NodesChi2(2*e+1))/2.0_DP,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,&
               ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
    FuncChi2(3*e+5) = A_IntEpsilon((NodesChi2(2*e+4)+NodesChi2(2*e+3))/2.0_DP,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,&
               ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
    ! Now it is time to compare the quadratures and to mark cells for the refinement 
    quad  = (FuncChi1(3*d-2)+4*FuncChi1(3*d-1)+FuncChi1(3*d))*(NodesChi1(2*d)-NodesChi1(2*d-1))/6
    quad1 = (FuncChi2(3*e+1)+4*FuncChi2(3*e+2)+FuncChi2(3*e+3))*(NodesChi2(2*e+2)-NodesChi2(2*e+1))/6
    quad2 = (FuncChi2(3*e+4)+4*FuncChi2(3*e+5)+FuncChi2(3*e+6))*(NodesChi2(2*e+4)-NodesChi2(2*e+3))/6
    es = abs(quad1+quad2-quad)/15/(NodesChi1(2*d)-NodesChi1(2*d-1)) ! error indicator
    my_int1=my_int1+abs(quad1)+abs(quad2) ! this quantity estimates the integral over the entire mesh NodesChi2 that will be NodesChi1 soon.. 
    if (es > ErrChi) then 
     ChiAdaptFlag = .true.
     CellsChiRef2(e+1) = 1
     CellsChiRef2(e+2) = 1
    else
     CellsChiRef2(e+1) = 0
     CellsChiRef2(e+2) = 0
     Atemp = Atemp + quad1 + quad2 
    end if
    e=e+2
   ! end refinement 
   end if 
   end do 
   !
   end if ! end of the check that the integral over the mesh NodesChi1 will be small anyway
   ! Finally, we need to replace the first mesh with the refined mesh 
   FuncChi1=0
   FuncChi1=FuncChi2
   NodesChi1=0
   NodesChi1=NodesChi2
   CellsChiRef1=0
   CellsChiRef1=CellsChiRef2
   Nchi1=e
   end do
   !!!!!!!!!!!!!!!!!! Looks like we are done with the integration in chi !!!!! 
   ! now it is time to check if the value of A is non-zero. If it is bigger than some specified level of 
   ! minimum sensitivity, then is it added to storage, Otherwise we continue to the next velocity #2.
   if ( ABS(Atemp) > min_sens ) then     
!$omp critical    
     A_ct = A_ct+1               ! we count the record for A
     phi_ct=phi_ct + 1           ! we count the record for this basis function. 
     if (An < A_ct) then
        call ExtendAarraysDGV(An,ni)
     end if  
      ! now we need to add the volume elements for xi and xi1 
      ! ATTENTION: molecular diameter mol_diam is removed in the dimensionless formulation)
      ! ATTENTION: for evaluation of A operator using the dimensional code use mol_diam = 1.0:
      ! ATTENTION: OLD (dimenional) code Atemp = Atemp*mod_diam^2*(nodes_gwts(i3))**(nodes_gwts(i2))
      ! Dimensionless code. Molecular diamter is accounted for in the spatial operator. 
      ! xi1:
      Atemp = Atemp*(nodes_gwts(i3)) ! this takes care of the volume elements and the weight for xi1 
      ! xi			
      Atemp = Atemp*(nodes_gwts(i2)) ! this takes care of the volume elements and the weight for xi 
      !		
      dsph=dsph/8.0_DP  ! this takes care of the fraction 1/8 still need to add |g|
      A(A_ct) = Atemp*dsph ! this takes care of |g|/8
      A_xi(A_ct) = i2
      A_xi1(A_ct) = i3 
      A_phi(A_ct) = i1
!!
!iiii = omp_get_thread_num()
!print *, "I2=", i2, "Thread", iiii, "A_ct=", A_ct, i3        
!$omp end critical           
   end if                    
!!!!!!!!!!!! end evaluation of A !!!!!!!!!!!
end do  ! END LOOP in I3
! add this two lines to track progress and Print a hello message to Check if parallel stuff worked... !!!!
!iiii = omp_get_thread_num()
!print *, "I2=", i2, "Thread", iiii, "A_ct=%i8", A_ct    
! 
end do ! END LOOP IN I2
A_capphi(i1)=phi_ct
print *, "Set_A i1=", i1, "A_ct=", A_ct
!! end do !! TEMPORARY do it for just one node...  
!!!!!!!!!!!!
call ShrinkAarraysDGV(A_ct)
!!!!!!!!!!!!
end subroutine SetA_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ExtendAarraysDGV(An,ni)
!
! This subroutine extends arrays A, A_phi,A_xi,A_xi1 for additional ni records

subroutine ExtendAarraysDGV(An,ni)

use commvar, only: A,A_phi,A_xi,A_xi1

integer (I4B), intent (out) :: An ! An is the length of the arrays that will be updated.  
integer (I4B), intent (in) :: ni ! ni is the number of records to add

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B) , dimension (:), allocatable :: A_iscr  ! integer scrap array
integer (I4B) :: nn ! scrap variable
!
integer (I4B) :: loc_alloc_stat
!
nn=size(A,1)
An=nn+ni
! extending A ... 
allocate (A_rscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_rscr, nn=", nn
  stop
  end if
A_rscr=A
deallocate(A)
allocate (A(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A, nn=", nn+ni
  stop
  end if
A(1:nn)=A_rscr
deallocate(A_rscr)
! end extending A

! extending A_xi
allocate(A_iscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_iscr, nn=", nn
  stop
  end if
A_iscr=A_xi
deallocate(A_xi)
allocate (A_xi(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_xi, nn=", nn+ni
  stop
  end if
A_xi(1:nn)=A_iscr
! extending A_xi1
A_iscr=A_xi1
deallocate(A_xi1)
allocate (A_xi1(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_xi1, nn=", nn+ni
  stop
  end if
A_xi1(1:nn)=A_iscr
!extending A_phi
A_iscr=A_phi
deallocate(A_phi)
allocate (A_phi(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_phi, nn=", nn+ni
  stop
  end if
A_phi(1:nn)=A_iscr
! end extending A_phi
deallocate(A_iscr)

end subroutine ExtendAarraysDGV

subroutine ShrinkAarraysDGV(ni)

use commvar, only: A,A_phi,A_xi,A_xi1

integer (I4B), intent (in) :: ni ! ni is the new size of the arrays to be shrunk

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B), dimension (:), allocatable :: A_iscr  ! integer scrap array
!
integer (I4B) :: loc_alloc_stat 

! shrinking A ... 
allocate (A_rscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_rscr, nn=", ni
  stop
  end if
A_rscr=A(1:ni)
deallocate(A)
allocate (A(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A, nn=", ni
  stop
  end if
A=A_rscr
deallocate(A_rscr)
! end shrinking A

! shrinking A_xi
allocate (A_iscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_iscr, nn=", ni
  stop
  end if
A_iscr=A_xi(1:ni)
deallocate(A_xi)
allocate (A_xi(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_xi, nn=", ni
  stop
  end if
A_xi=A_iscr
! extending A_xi1
A_iscr=A_xi1(1:ni)
deallocate(A_xi1)
allocate (A_xi1(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_xi1, nn=", ni
  stop
  end if
A_xi1=A_iscr
!extending A_phi
A_iscr=A_phi(1:ni)
deallocate(A_phi)
allocate (A_phi(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_phi, nn=", ni
  stop
  end if
A_phi=A_iscr
! end extending A_phi
deallocate(A_iscr)
end subroutine ShrinkAarraysDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! IphiDGV
!
! This subroutine evaluates the moment of the collision integral for basis functions. 
! The basis functions are defined for each velocity node. The Collision Information Operator
! is calculated for all (or few basis functions). 
! The subroutine will take calculate the moment for each basis function recorded in A
!
! Depends on the main program. A-arrays must be defined as well as the disftribution function f
! 
!!!!!!!!!!!!!!!!!!!!!!!!

function IphiDGV(f) result (Iphi)

use commvar, only: A,A_capphi,A_xi,A_xi1,A_phi,nodes_gwts,mol_diam,L_inf,N_inf

real (DP), dimension(:) :: f ! the main variable -- distribution function -- one component for each velocity node
real (DP), dimension (size(A_capphi,1)) :: Iphi ! the result of integration
real (DP), dimension(:), allocatable :: A_SP, F_SP
!!! 
integer (I4B) :: i,j, A_ct,my_count ! scrap variables,
!!!
allocate (A_SP(size(A,1)), F_SP(size(f,1)))
A_SP = A
f_SP = f
!!!
my_count=0
Iphi=0
A_ct=0
do i=1,size(A_capphi,1) ! loop in basis functions --- one node=one basis function
   Iphi(i)=0
   do j=1+A_ct,A_ct+A_capphi(i)
   Iphi(i)=Iphi(i)+A_SP(j)*f_SP(A_xi(j))*f_SP(A_xi1(j))
   if (ABS(A_SP(j)*f_SP(A_xi(j))*f_SP(A_xi1(j)))<1.0d-8) then 
   my_count=my_count+1
   end if 
   end do 
   Iphi(i)=2*Iphi(i)/nodes_gwts(i)*((mol_diam/L_inf)**2*N_inf)
   A_ct=A_ct+A_capphi(i)
end do 
!
deallocate (A_SP, f_SP)
end function IphiDGV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! IphiDGV_decomp
!
! This is a diagnistics subroutine 
!
! This subroutine evaluates the moment of the collision integral for basis functions. 
! using the decomposition method. It is to test the advantage of decomposition -- if any 
! the solution is split into a maxwellian adn the rest. The collision operator is evcaluated.
! The basis functions are defined for each velocity node. The Collision Information Operator
! is calculated for all (or few basis functions). 
! The subroutine will take calculate the moment for each basis function recorded in A
!
! Depends on the main program. A-arrays must be defined as well as the disftribution function f
! 
!!!!!!!!!!!!!!!!!!!!!!!!

function IphiDGV_decomp(f) result (Iphi)

use commvar, only: A,A_capphi,A_xi,A_xi1,A_phi,nodes_gwts,mol_diam,L_inf,N_inf,nodes_u, nodes_v, nodes_w

use distributions_mod



real (DP), dimension(:) :: f ! the main variable -- distribution function -- one component for each velocity node
real (DP), dimension (size(A_capphi,1)) :: Iphi ! the result of integration

!!! 
integer (I4B) :: i,j, A_ct ! scrap variables,
real (DP), dimension(size(f,1)) :: fmaxwellNew,Df
real (DP) :: LocDens,LocUbar,LocVbar,LocWbar,LocTempr 
!!!!!!!
call MassCheckRec (f,LocDens,LocUbar,LocVbar,LocWbar,LocTempr) ! First, we evaluate the macroparamters of the solution.
fMaxwellNew = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_u, nodes_v, nodes_w) ! now we populate the maxwellian with the same macroparamters.
Df=f-fMaxwellNew ! evaluate the perturbation from the maxwellian. 

Iphi=0
A_ct=0
do i=1,size(A_capphi,1) ! loop in basis functions --- one node=one basis function
   Iphi(i)=0
   do j=1+A_ct,A_ct+A_capphi(i)
   Iphi(i)=Iphi(i) + A(j)*fMaxwellNew(A_xi(j))*Df(A_xi1(j))+A(j)*Df(A_xi(j))*fMaxwellNew(A_xi1(j)) + A(j)*Df(A_xi(j))*Df(A_xi1(j))
   end do 
   Iphi(i)=2*Iphi(i)/nodes_gwts(i)*((mol_diam/L_inf)**2*N_inf)
   A_ct=A_ct+A_capphi(i)
end do 
!
end function IphiDGV_decomp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!
!EvalPhiPostCollIntegrand
!
! This function helps to calculate the post collision velocities and calls evaluation of the basis function
! on the post collision velocitis. The function is mainly to help coding -- to reduce the number of repeating lines. 
!
!

function EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
             ugux3,ugvx3,ugwx3,epsil,sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1) result (y)
!
!
real (DP), intent (in) :: sinchi, coschi, g1,g2                           ! some useful numbers 
real (DP), intent (in) :: ugu,ugv,ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3 ! useful coefficients
real (DP), intent (in) :: epsil ! the angle for which the integral must be evaluated.
real (DP), intent (in) :: dsph,xiu,xiv,xiw,xi1u,xi1v,xi1w ! the pre-collision velocities 
integer (I4B), intent (in) :: i1 ! the number of the basis fucntion to evaluate
real (DP), intent (in) :: varphi_xi,varphi_xi1 ! the values of the basis function on xi and xi1
real (DP) :: y ! the result

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP) :: sinchicoseps, sinchisineps ! useful quantities
real (DP) :: xiupr,xivpr,xiwpr,xi1upr,xi1vpr,xi1wpr ! components of post-collision velocity
real (DP) :: z,ku,kv,kw ! compoents of the unit vector in the direction of post collision relative speed
!!!!!!
     sinchicoseps = cos(epsil)*sinchi/g1
     sinchisineps = sin(epsil)*sinchi/g2
   ! Now it is time to evaluate post-collision velocities... 
   ! First we evaluate the direction of $g'$ : 
     ku = (ugu)*coschi - ugux2*sinchicoseps - ugux3*sinchisineps
     kv = (ugv)*coschi - ugvx2*sinchicoseps - ugvx3*sinchisineps
     kw = (ugw)*coschi - ugwx2*sinchicoseps - ugwx3*sinchisineps
   ! Now we calculate the post collision velocities: 
   !
     xiupr=((xiu + xi1u) + dsph*ku)/2.0_DP
     xi1upr=((xiu + xi1u) - dsph*ku)/2.0_DP
     xivpr=((xiv + xi1v) + dsph*kv)/2.0_DP
     xi1vpr=((xiv + xi1v) - dsph*kv)/2.0_DP
     xiwpr=((xiw + xi1w) + dsph*kw)/2.0_DP
     xi1wpr=((xiw + xi1w) - dsph*kw)/2.0_DP

   ! now we evaluate basis the basis functions on these velocities 
     z = EvalLagrBasisFunByNdsDGblzmEZ(xiupr,xivpr,xiwpr,i1) + &
       EvalLagrBasisFunByNdsDGblzmEZ(xi1upr,xi1vpr,xi1wpr,i1) - varphi_xi - varphi_xi1 
     y=z  
end function EvalPhiPostCollIntegrand     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A_IntEpsilon
! This function is to shorten the code. It is pretty much a piece of code that is cut out
! and paste in the form of a function to make the rest of the code read easier. 
! This portion implements the integration in \varepsilon in the evaluation of operator A
!
! This code evaluates \sin\chi \int_{0}^{2\pi} (\varphi(\xi')+\varphi(\xi'_{1}) d\varepsilon 
!
!
! The result is the integral for this particular angle \chi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function A_IntEpsilon(chi,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1) result (y)
             
real (DP), intent (in) :: chi,g1,g2,ires_e                    ! some useful numbers 
real (DP), intent (in) :: ugu,ugv,ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3 ! useful coefficients
real (DP), intent (in) :: dsph,xiu,xiv,xiw,xi1u,xi1v,xi1w ! the pre-collision velocities 
integer (I4B) :: i1 ! the number of the basis function which is to evaluate. 
real (DP), intent (in) :: ErrEps ! the max set error of integral evauation. 
real (DP), intent (in) :: varphi_xi,varphi_xi1 ! the values of the basis function on xi and xi1
!
real (DP) :: y ! the result
             

!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP) :: sinchi,coschi ! useful scrap constants to remember..
real (DP), parameter :: pi25DT = 3.141592653589793238462643d0 
integer (I4B), parameter :: NdsEps=10000 ! The max number of cells for integration in chi and epsilon

real (DP), dimension (2*NdsEps) :: NodesEps1, NodesEps2
real (DP), dimension (3*NdsEps) :: FuncEps1, FuncEps2 ! Arrays to store cells and functin values on cells Two copies are needed. Integration in Eps
integer (I4B), dimension (NdsEps) :: CellsEpsRef1,CellsEpsRef2 ! Arrays to keep refinement flags
real (DP) :: Atemp_temp,quad,quad1,quad2,es,my_int1 ! scrap variables.
!!!!!
integer (I4B) :: Neps1,g ! variables to keep the number of integration cells
integer (I4B) :: f ! local counter
logical :: EpsAdaptFlag ! a logical variable to tell if any of the refinements need to be done. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   sinchi=sin(chi) ! will be useful to remember this one... 
   coschi=cos(chi) ! will be useful to remember this one too... 
! Begin Integration in epsilon ...  
   Neps1=FLOOR(pi25DT*dsph*sinchi/ires_e)+1
   if (Neps1>NdsEps) then 
     Neps1=NdsEps
     print *,"SetA_DGV: Warning! Number of nodes (NdsEps) for integration in epsilon gives insufficient resolution"  
   end if 
   ! we set the initial mesh
   do f=1,Neps1
    NodesEps1(2*(f-1)+1)=(f-1)*2*pi25DT/Real(Neps1,DP)
    NodesEps1(2*f)=f*2*pi25DT/Real(Neps1,DP)
   end do
   ! now we evaluate the integrand on the initial mesh
   FuncEps1=0
   do f=1,Neps1
   ! each f gives one cell
   ! First, we evaluate the integrand on the left node of the cell
    if ((f>1) .and. (NodesEps1(2*f-1) == NodesEps1(2*f-2))) then 
     FuncEps1(3*f-2)= FuncEps1(3*f-3) ! if the node is repeating, the value has been computed already
    else 
     FuncEps1(3*f-2) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
                ugux3,ugvx3,ugwx3,NodesEps1(2*f-1),sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1) 
    end if
   ! This takes care of the left node in the cell... 
   !  
   !Now we evaluate the solution on the right node of the cell and at the center. 
   FuncEps1(3*f) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,NodesEps1(2*f),sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1) 
   FuncEps1(3*f-1) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,(NodesEps1(2*f)+NodesEps1(2*f-1))/2,sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1)
   ! repeat for all cells... 
   end do      
   ! finally, we mark all cells for the refinement in the beginning... 
   CellsEpsRef1=0 
   CellsEpsRef1(1:Neps1)=1
   ! 
   !Next we adaptively refine and evaluate the integral:    
   my_int1=10*ErrEps !! set it to a large number fo the rist run
   EpsAdaptFlag =.true.
   Atemp_temp=0
   do while (EpsAdaptFlag)
   EpsAdaptFlag =.false.
   ! Now we will go over the array of the integrand values and will cross 
   ! out cells where the integrand is zero, and those that are not marked for the refinment 
   ! on those cells we summ the integrand and keep the sum.
   ! Cells marked for the refinement we divide into halves 
   NodesEps2=0
   FuncEps2=0
   CellsEpsRef2=0
   g=0
   if (my_int1 > ErrEps) then   ! check is the integral over NodesEps1 is expected to be small. If it is small then the refinement is not performed
   my_int1=0 ! this will keep the estimate of the integral on the mesh NodesEps1 so that we stop resolving when this number is small
   do f=1,Neps1
   if ((abs(FuncEps1(3*f-2))+abs(FuncEps1(3*f-1))+abs(FuncEps1(3*f)) > 1.0d-15) .and. (CellsEpsRef1(f)==1)) then 
    if (g > NdsEps-2) then 
     print *,"SetA_DGV: Number of nodes in epsilon is too big (g>NdsEps-2)"
     stop
    end if
    ! save the nodes and the integrand values
    NodesEps2(2*g+1)=NodesEps1(2*f-1)
    FuncEps2(3*g+1)=FuncEps1(3*f-2)
    NodesEps2(2*g+4)=NodesEps1(2*f)
    FuncEps2(3*g+6)=FuncEps1(3*f)
    ! Evaluate the midpoint node: 
    NodesEps2(2*g+2)=(NodesEps1(2*f-1)+NodesEps1(2*f))/2
    NodesEps2(2*g+3)=NodesEps2(2*g+2)
    ! Save the midpoint value of the integrand:
    FuncEps2(3*g+3)=FuncEps1(3*f-1)
    FuncEps2(3*g+4)=FuncEps2(3*g+3)
    ! Evaluate the integrand on the new nodes
    FuncEps2(3*g+2) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
    ugux3,ugvx3,ugwx3,(NodesEps2(2*g+2)+NodesEps2(2*g+1))/2.0_DP,sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1) 
    FuncEps2(3*g+5) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
    ugux3,ugvx3,ugwx3,(NodesEps2(2*g+4)+NodesEps2(2*g+3))/2.0_DP,sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1) 
    ! Now it is time to compare the quadratures and to mark cells for the refinement 
    quad  = (FuncEps1(3*f-2)+4*FuncEps1(3*f-1)+FuncEps1(3*f))*(NodesEps1(2*f)-NodesEps1(2*f-1))/6
    quad1 = (FuncEps2(3*g+1)+4*FuncEps2(3*g+2)+FuncEps2(3*g+3))*(NodesEps2(2*g+2)-NodesEps2(2*g+1))/6
    quad2 = (FuncEps2(3*g+4)+4*FuncEps2(3*g+5)+FuncEps2(3*g+6))*(NodesEps2(2*g+4)-NodesEps2(2*g+3))/6
    es = abs(quad1+quad2-quad)/15/(NodesEps1(2*f)-NodesEps1(2*f-1)) ! local error indicator for simpson's rule.
    my_int1 = my_int1 + abs(quad1) + abs(quad2) ! This will calculate the estimate from above to the total integral on the current mesh 
    if (es > ErrEps) then 
     EpsAdaptFlag = .true.
     CellsEpsRef2(g+1) = 1
     CellsEpsRef2(g+2) = 1
    else
     CellsEpsRef2(g+1) = 0
     CellsEpsRef2(g+2) = 0
     Atemp_temp = Atemp_temp + quad1 + quad2 
    end if                   
    g=g+2
   ! end refinement 
   end if 
   end do 
   !
   end if ! end check if the integral over NodesEps1 is too small to worry about it. 
   ! Finally, we need to replace the first mesh with the refined mesh 
   FuncEps1=0
   FuncEps1=FuncEps2
   NodesEps1=0
   NodesEps1=NodesEps2
   CellsEpsRef1=0
   CellsEpsRef1=CellsEpsRef2
   Neps1=g
   end do
   !!!!!!!!!!!!!!!!!! Looks like we are done with the integration in epsilon !!!!! 
   y=Atemp_temp*sinchi
end function A_IntEpsilon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MassCheck
! 
! this fucntions will try to calculate mass of the disctribution function. This is a debug subroutine
! 
!!!!!!!!!!!!!!!!!
function MassCheck (f) result (y)
use commvar, only: nodes_gwts,nodes_u,nodes_v,nodes_w,gasR

real (DP), dimension (:), intent (in) :: f !the vector of nodal values of the distribution function
real (DP) :: y ! the result (the mass) 
!!!!!!!!!!!!!!!!!!
real (DP) :: n,ubar,vbar,wbar,temp ! number density, av_v
 
integer (I4B) :: i ! scrap index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check mass
n=sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check momentum 
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
! check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!
y=ubar ! return this... 

end function MassCheck  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FindCellContainsPoint_DGV
!
! This subroutine find the number of the active cell in the cell arrays that containes
! the point (u,v,w)
!
! the function returns zero if the velocity point is not on any cells. 
!
!
! the subroutine accesses arrays
! nodes_pcell, nodes_u, nodes_v, nodes_w,&
!				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
!                  cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
!                  cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
!                  grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
!                  nodes_ui, nodes_vi, nodes_wi
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function FindCellContainsPoint_DGV(u,v,w) result (y)

use commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w,&
				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
                   nodes_ui, nodes_vi, nodes_wi

real (DP) :: u,v,w ! the components of the given velocity point
integer (I4B) ::  y ! the value of the index in the cell arrays that correcponds to the active cell containing point (u,v,w)

!!!!!!!!!!!!!!                      
integer (I4B) :: j ! a counter -- usually denotes the number of the grid we are working on. 
integer (I4B) :: gi_v,gi_u,gi_w ! counters to help count nodes in the grids array
integer (I4B) :: celli, cellp ! these will store the number of the cell where the velocity belongs and the number of the 
            ! parent cell where the basis function belongs. 
integer (I4B) :: ui,vi,wi ! are the local numbers on the particular grids where u belongs. Because we have ierarchicval grids, 
             ! there may be different grids that have ui,vi,wi defined. If one of these numbers is zero -- then the 
             ! velocity u,v,w is not on the grid. Unless a mistake was done in setting the grids either all three of them 
             ! are zero -- that means that the velocity is not on this (level zero grid) or all three are non-zero  -- means that
             ! the velocity has been found. Mized results may indicate a breach in data integrity. 
             !!!!!
integer (I4B) :: jj ! a counter
real (DP) :: unor,vnor,wnor ! scrap variables to store the normalized velus of the components u,v,w             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! first, given a velocity, we need to identify the cell where it came from and its 1D numbers 
! we start looking in course cells first. If a cell is found that has that velocity we check if the cell is refined. 
! if it is, then we look into the corresponding grid and seek the smaller cell that has it and so on. 
!!!!!!!!
!
j=1
!! set up the shift in grids_u/_v/_w corresponding to the first grid (zero shift):
gi_u=0
gi_v=0
gi_w=0
do while (j<=size(grids_cap_u,1))
  call finduiviwi(grids_u(gi_u+1:gi_u+grids_cap_u(j)),grids_v(gi_v+1:gi_v+grids_cap_v(j)),&
                       grids_w(gi_w+1:gi_w+grids_cap_w(j)),u,v,w,ui,vi,wi)
  if (ui*vi*wi /= 0) then 
     ! now we need to find the cell that correspond to this grid and these ui,vi,wi
     celli=1
     do jj=1,j-1
      celli = celli + (grids_cap_u(jj)-1)*(grids_cap_v(jj)-1)*(grids_cap_w(jj)-1)
     end do
     celli=celli + (ui-2)*(grids_cap_v(j)-1)*(grids_cap_w(j)-1)+(vi-2)*(grids_cap_w(j)-1)+(wi-2)
     !  check if this cell is refined
     if ((cells_refu(celli)>1) .or. (cells_refv(celli)>1) .or. (cells_refw(celli)>1)) then  
        j=cells_cgrid(celli)
        ! set up the shift in grids_u/_v/_w corresponding to the j'th grid:
        gi_u=0
        gi_v=0
        gi_w=0
        do jj=1,j-1
         gi_u = gi_u + grids_cap_u(jj)
         gi_v = gi_v + grids_cap_v(jj)
         gi_w = gi_w + grids_cap_w(jj)
        end do
     else 
        exit
     endif
  else 
     celli=0
     exit
  endif                      
enddo
! now "celli" either is the number of the cell or 0 (0 means that the velocity in not on any cell)
! we return this number 
y=celli
! 
end function FindCellContainsPoint_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FindI1sByCellNum_DGV
!
! This functions looks at the nodes arrays and finds the numbers of the nodes that
! belong to cell with the number i
! 
! The numbers will be recorded in the array of results with the first number indicating the number of useful records/
! 
! For example: if the result array has a zero in its first place, then the cell has no nodes -- which is not expected, really... 
! 
! if the array has number 5 in its firs place, than elements 2--6 contain the I1s --- these numbers will be used to build the A-array.
! 
!  if there is no cell with the number i -- the program retrurns zero records and prints a worning
!  if the result array is too short to fit all the numbers i1, the program prints the error message and stops.
!
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FindI1sByCellNum_DGV(i,y)

use commvar, only: nodes_pcell,cells_lu
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), intent (in) :: i ! the number of the cell where I1s need to be looked at. 
integer (I4B), dimension (:), intent (out) :: y ! the numbers of the nodes (I1 -- in our sleng...) that belong to the cell with number i 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (i4B) :: j,k,pcell_i,sizey ! local counters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sizey=size(y,1)
y=0 ! nullify the thing... Just In Case 
k=0 ! in the beginning there is no records///
if (i>size(cells_lu,1)) then  ! check if i is too big
 y=0 !
 print *,"FindI1sByCellNum_DGV: the cell with the provided number i do not exist -- i is too big. Returned zero records"
else ! if i is not too big, 
 do j=1,size(nodes_pcell,1)
  pcell_i=nodes_pcell(j) 
  if (pcell_i==i) then ! check if the node belongs to cell i 
   k=k+1                ! if it does, check if can still records in y
   if (k+1<=sizey) then  
    y(1+k)=j     !record
   else  ! othersize print the error message and stop
    print *,"FindI1sByCellNum_DGV: The size fo the result array is too small. Can not put all I1s. Stop"
    stop
   end if
  end if  
 end do 
y(1)=k ! the first records is reserved for the number of found I1s
end if 
end subroutine FindI1sByCellNum_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TruncAarrsTrhld_DGV (trhld)
!
! This subroutine trancates the A-arrays. All entries of A that are below the provided 
! treshhold value (trhld) are being cut, arrays Axi, Axi1, Aphi, A_capphi are updated correspondingly... 
!
! ATTENTION: array nodes_Ashift need to be re-created after A has been truncated.
!
!
!!!!!!!!!!!!!!!!!!!
subroutine TruncAarrsTrhls_DGV (trhld)

use commvar, only: A, A_xi, A_xi1, A_phi, A_capphi 

!!!!
real (DP), intent (in) :: trhld ! the treshhold at which to cut A:
!!!
integer (I4B) :: nold,nnew ! are the old and new sizes of the A-arrays
integer (I4B) :: i,j,phicap ! scrap indices
!!!
real (DP), dimension (:), allocatable :: Ascr ! real scrap array
integer (I4B), dimension (:), allocatable :: A_xiscr, A_xi1scr, A_phiscr  ! integer scrap arrays
!
integer (I4B) :: loc_alloc_stat 

nold = size(A,1)
! shrinking A ... 
allocate (Ascr(1:nold), A_xiscr(1:nold), A_xi1scr(1:nold), A_phiscr(1:nold), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "TruncAarrsTrhls_DGV: Allocation error for variables  Ascr, A_xiscr, _xi1scr, _phiscr" 
  stop
  end if
! Save the A arrays .... 
Ascr=A; A_xiscr=A_xi; A_xi1scr=A_xi1; A_phiscr=A_phi
! Now clear the A-arrays;  
A=0.0_DP; A_xi=0; A_xi1=0; A_phi=0
! Now start to fill them in, but omitting all records that are below the treshhold value
nold=0 ! this will count the original records of A  === similar to A_ct index
nnew=0 ! this will count the records after the truncation. 
do i = 1,size(A_capphi,1) ! loop in basis functions --- one node=one basis function
   phicap = 0
   do j = 1+nold,nold+A_capphi(i)
    if (ABS(Ascr(j))>trhld) then  ! if records in A array are above treshhold, they are saved, otherwise ignored
     phicap = phicap+1
     nnew = nnew+1 
     A(nnew)=Ascr(j)
     A_xi(nnew)=A_xiscr(j)
     A_xi1(nnew)=A_xi1scr(j)
     A_phi(nnew)=A_phiscr(j)
    end if 
   end do 
   nold = nold+A_capphi(i) ! shift the start index to the place where record sfo rthe next basis function start.
   A_capphi(i) = phicap ! Update the A_capphi --- it now stores the number of records after the truncation
end do 
deallocate(Ascr,A_xiscr,A_xi1scr,A_phiscr) !

call ShrinkAarraysDGV(nnew)

end subroutine TruncAarrsTrhls_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TruncAarrsRadius_DGV (rad)
!
! This subroutine trancates the A-arrays based on the distance between xi and xi1. If the distance is larger than 
! the specified radus (rad), the record is discarded from A and arrays Axi, Axi1, Aphi, A_capphi are updated correspondingly... 
!
! ATTENTION: array nodes_Ashift needs to be re-created after A has been truncated.
!
!
!!!!!!!!!!!!!!!!!!!
subroutine TruncAarrsRadius_DGV (rad)

use commvar, only: A, A_xi, A_xi1, A_phi, A_capphi, nodes_u, nodes_v, nodes_w 

!!!!
real (DP), intent (in) :: rad ! the treshhold at which to cut A:
!!!
integer (I4B) :: nold,nnew ! are the old and new sizes of the A-arrays
integer (I4B) :: i,j,phicap,ixi,ixi1 ! scrap indices
!!!
real (DP), dimension (:), allocatable :: Ascr ! real scrap array
integer (I4B), dimension (:), allocatable :: A_xiscr, A_xi1scr, A_phiscr  ! integer scrap arrays
real (DP) :: dsq, radsq ! scrap variables to keep distance squared and radius squared
!
integer (I4B) :: loc_alloc_stat 

nold = size(A,1)
! shrinking A ... 
allocate (Ascr(1:nold), A_xiscr(1:nold), A_xi1scr(1:nold), A_phiscr(1:nold), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "TruncAarrsTrhls_DGV: Allocation error for variables  Ascr, A_xiscr, _xi1scr, _phiscr" 
  stop
  end if
! Save the A arrays .... 
Ascr=A; A_xiscr=A_xi; A_xi1scr=A_xi1; A_phiscr=A_phi
! Now clear the A-arrays;  
A=0.0_DP; A_xi=0; A_xi1=0; A_phi=0
! Now start to fill them in, but omitting all records that are below the treshhold value
nold=0 ! this will count the original records of A  === similar to A_ct index
nnew=0 ! this will count the records after the truncation. 
radsq=rad**2
do i = 1,size(A_capphi,1) ! loop in basis functions --- one node=one basis function
   phicap = 0
   do j = 1+nold,nold+A_capphi(i)
    ixi=A_xiscr(j)
    ixi1=A_xi1scr(j)
    dsq = (nodes_u(ixi)-nodes_u(ixi1))**2 + (nodes_v(ixi)-nodes_v(ixi1))**2 + (nodes_w(ixi)-nodes_w(ixi1))**2
    if (dsq <= radsq) then  ! if vectors xi and xi1 are at or closer than distance rad, the record in A is saved, otherwise ignored
     phicap = phicap+1
     nnew = nnew+1 
     A(nnew)=Ascr(j)
     A_xi(nnew)=A_xiscr(j)
     A_xi1(nnew)=A_xi1scr(j)
     A_phi(nnew)=A_phiscr(j)
    end if 
   end do 
   nold = nold+A_capphi(i) ! shift the start index to the place where record sfo rthe next basis function start.
   A_capphi(i) = phicap ! Update the A_capphi --- it now stores the number of records after the truncation
end do 
deallocate(Ascr,A_xiscr,A_xi1scr,A_phiscr) !

call ShrinkAarraysDGV(nnew)

end subroutine TruncAarrsRadius_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MassCheckRec
! 
! this fucntions will try to calculate mass of the disctribution function. This is a debug subroutine
! 
!!!!!!!!!!!!!!!!!
subroutine MassCheckRec (f,n,ubar,vbar,wbar,temp)
use commvar, only: nodes_gwts,nodes_u,nodes_v,nodes_w,gasR

real (DP), dimension (:), intent (in) :: f !the vector of nodal values of the distribution function
real (DP), intent (out) :: n,ubar,vbar,wbar,temp ! number density, av_v
!!!!!!!!!!!!!!!!!!
 
integer (I4B) :: i ! scrap index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check mass
n=sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check momentum 
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
! check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!
end subroutine MassCheckRec  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MassCheckRecPlus
! 
! this fucntions will try to calculate mass of the disctribution function. This is a debug subroutine
! This one also calculase momentum, temperature and the directional temperature 
!!!!!!!!!!!!!!!!!
subroutine MassCheckRecPlus (f,n,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w)
use commvar, only: nodes_gwts,nodes_u,nodes_v,nodes_w,gasR

real (DP), dimension (:), intent (in) :: f !the vector of nodal values of the distribution function
real (DP), intent (out) :: n,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w ! number density, av_v
!!!!!!!!!!!!!!!!!!
 
integer (I4B) :: i ! scrap index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check mass
n=sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check momentum 
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
! check temperature 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
! check directional temperatures 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp_u = sum(f*nodes_gwts*(nodes_u-ubar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tx
temp_v = sum(f*nodes_gwts*(nodes_v-vbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Ty
temp_w = sum(f*nodes_gwts*(nodes_w-wbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tz
!!!!!!!!!!
end subroutine MassCheckRecPlus  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MassCheckRecHighCMoments
! 
! this fucntions will try to calculate mass of the disctribution function. This is a debug subroutine
! This one also calculates momentum, temperature and the directional temperature as well as additional moments
! and will update the moments history
!!!!!!!!!!!!!!!!!
subroutine MassCheckRecHighCMoments (f,n,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w,mom3_u,mom3_v,mom3_w,mom4_u,mom4_v,&
   mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w)
use commvar, only: nodes_gwts,nodes_u,nodes_v,nodes_w,gasR

real (DP), dimension (:), intent (in) :: f !the vector of nodal values of the distribution function
real (DP), intent (out) :: n,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w,mom3_u,mom3_v,mom3_w
real (DP), intent (out) :: mom4_u,mom4_v,mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w ! number density, av_v
!!!!!!!!!!!!!!!!!!

integer (I4B) :: i ! scrap index

!!!!!!!!!
! check mass
n=sum(f*nodes_gwts)
!!!!!!!!!
! check momentum 
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
! check temperature 
!!!!!!!!!
temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
! check directional temperatures 
!!!!!!!!!
temp_u = sum(f*nodes_gwts*(nodes_u-ubar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tx
temp_v = sum(f*nodes_gwts*(nodes_v-vbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Ty
temp_w = sum(f*nodes_gwts*(nodes_w-wbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tz
!!!!!!!!!
! Moments for Q = u^3
mom3_u = sum(f*nodes_gwts*(nodes_u-ubar)**3)/n
mom3_v = sum(f*nodes_gwts*(nodes_v-vbar)**3)/n
mom3_w = sum(f*nodes_gwts*(nodes_w-wbar)**3)/n
!!!!!!!!!
! Moments for Q = u^4
mom4_u = sum(f*nodes_gwts*(nodes_u-ubar)**4)/n
mom4_v = sum(f*nodes_gwts*(nodes_v-vbar)**4)/n
mom4_w = sum(f*nodes_gwts*(nodes_w-wbar)**4)/n
!!!!!!!!!
! Moments for Q = u^5
mom5_u = sum(f*nodes_gwts*(nodes_u-ubar)**5)/n
mom5_v = sum(f*nodes_gwts*(nodes_v-vbar)**5)/n
mom5_w = sum(f*nodes_gwts*(nodes_w-wbar)**5)/n
!!!!!!!!!
! Moments for Q = u^6
mom6_u = sum(f*nodes_gwts*(nodes_u-ubar)**6)/n
mom6_v = sum(f*nodes_gwts*(nodes_v-vbar)**6)/n
mom6_w = sum(f*nodes_gwts*(nodes_w-wbar)**6)/n
!!!!!!!!!
end subroutine MassCheckRecHighCMoments

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RecordMomDerivRelaxTimes_DGV
!! 
!!
!! This is a diagnostics subroutine that takes records derivatives and relaxation speeds of
!! a selected group of moments. 
!! 
!! The derivatives of moments are evaluated for the problem of spatially homogeneous relaxation. 
!! The derivatives are thus computed by taking the moment of the collision integral. Also, derivatives of the local maxwellian are computed
!! for true problem of spatial relaxatiopn this derivative is zero, however, due to the numerical erros, the 
!! maxwellian will change. These errors will tell us how bad is the violation of the conservation laws 
!! 
!! The relaxation speeds will be computed from the following definition: 
!! \nu_{\phi} = -\partial_{t} \ln | f_{\varphi}(t) - f_{\varphi}^M (t)|
!! 
!! 
!! The records will be stored in an array - MomDerivRelaxT with the first two 
!! numbers being (1) the number of saved records -- how many times the derivatives and speeds were evaluated. 
!! (2) the number of columns in the record. The first column will the the time 
!! then there will go the columns for the derivatives and then for speeds. The number of columns is then 1+2*m where m 
!! is the number of moments. 
!! 
!! MomDerivRelaxT -- array to keep the data. Size of this array is (1+2*m)*num_eval_error 
!!				
!! MomDerivRelaxT_flag -- this variable will keep the status of the allocation of 
!!
!! On the first call, the MomDerivRelaxT array will be allocated. Then adata will recorder in this array, adding one record 
!! every time the subroutine is called. The same subrouine will deallocated the array after the exectution is finished.  
!! In the MPI implementatio the subroutine will run on the master node. 
!! 
!! variables: 
!! 
!! from commvar: MomDerivRelaxT, MomDerivRelaxT_flag
!! from commvar: frhs1 -- last evaluated collision operator
!! from commvar: f1 -- value of f that corresponds to that frhs1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RecordMomDerivRelaxTimes_DGV(curr_time)

use commvar, only: MomDerivRelaxT,MomDerivRelaxT_flag,f1,frhs1,num_eval_error,t_R,nodes_u,nodes_v,nodes_w,&
                   rkmts,nodes_gwts 
use distributions_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), intent (in) :: curr_time ! the currect time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: i,j,k ! scrap indices 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: loc_alloc_stat, info  ! local variable to keep the allocation status
real (DP), parameter :: epps=0.01_DP ! threshhold parameters for evaluation of the relaxation times
real (DP) :: tresh ! this is a scrap variable to use in evaluation of relaxation times
real (DP) :: n,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w,mom3_u,mom3_v,mom3_w ! scrap variables to store the moments
real (DP) :: mom4_u,mom4_v,mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w ! number density, av_v
real (DP) :: mn,mubar,mvbar,mwbar,mtemp,mtemp_u,mtemp_v,mtemp_w,mmom3_u,mmom3_v,mmom3_w ! scrap variables to store the difference in the moments between the solution and the local maxwellian
real (DP) :: mmom4_u,mmom4_v,mmom4_w,mmom5_u,mmom5_v,mmom5_w,mmom6_u,mmom6_v,mmom6_w ! number density, av_v
real (DP) :: dn,dubar,dvbar,dwbar,dtemp,dtemp_u,dtemp_v,dtemp_w,dmom3_u,dmom3_v,dmom3_w ! scrap vaiables to store derivatives of the moments...
real (DP) :: dmom4_u,dmom4_v,dmom4_w,dmom5_u,dmom5_v,dmom5_w,dmom6_u,dmom6_v,dmom6_w !  
real (DP) :: taun,tauubar,tauvbar,tauwbar,tautemp,tautemp_u,tautemp_v,tautemp_w,taumom3_u,taumom3_v,taumom3_w ! scrap vaiables to store relaxation times for the moments...
real (DP) :: taumom4_u,taumom4_v,taumom4_w,taumom5_u,taumom5_v,taumom5_w,taumom6_u,taumom6_v,taumom6_w !  
real (DP), dimension (size(f1,1)) :: fm,derfm ! this is the storage for the local Maxwellian and the derivative of the local maxwellian 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! First the program will check if the array to store the moment derivatives has been created yet. 
! If it was not, it will create the array.
if (MomDerivRelaxT_flag /= 1) then 
  k = 1 + 2*20 !1+3+1+5x3 = 20 -- number of saved moments
             ! k is the total number of columents in the array.
  allocate (MomDerivRelaxT(1:2+k*(num_eval_error+3+rkmts+1)), stat=loc_alloc_stat) ! we added rkmts so that we can save on each of the step of the prepare MTS
  if (loc_alloc_stat >0) then 
   print *, "RecordMomDerivRelaxTimes_DGV: Allocation error for variable (MomDerivRelaxT)"
   stop
  end if              
  MomDerivRelaxT(1)=0 ! this cell stores the number of records
  MomDerivRelaxT(2) = k  ! this cell stored the number of columns
  MomDerivRelaxT_flag = 1  ! set the flag to 1 to indicate that the array was created 
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we compute the moments, we will need them  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!
! check mass
n=sum(f1*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!
! check momentum 
!!!!!!!!!!!!!!!!!!!!!!!!!
ubar=sum(f1*nodes_gwts*nodes_u)/n
vbar=sum(f1*nodes_gwts*nodes_v)/n
wbar=sum(f1*nodes_gwts*nodes_w)/n
! check temperature 
!!!!!!!!!!!!!!!!!!!!!!!!!
temp = sum(f1*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!!!!!!!!!
! check directional temperatures 
!!!!!!!!!
temp_u = sum(f1*nodes_gwts*(nodes_u-ubar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tx
temp_v = sum(f1*nodes_gwts*(nodes_v-vbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Ty
temp_w = sum(f1*nodes_gwts*(nodes_w-wbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tz
!!!!!!!!!
! Moments for Q = u^3
mom3_u = sum(f1*nodes_gwts*(nodes_u-ubar)**3)/n
mom3_v = sum(f1*nodes_gwts*(nodes_v-vbar)**3)/n
mom3_w = sum(f1*nodes_gwts*(nodes_w-wbar)**3)/n
!!!!!!!!!
! Moments for Q = u^4
mom4_u = sum(f1*nodes_gwts*(nodes_u-ubar)**4)/n
mom4_v = sum(f1*nodes_gwts*(nodes_v-vbar)**4)/n
mom4_w = sum(f1*nodes_gwts*(nodes_w-wbar)**4)/n
!!!!!!!!!
! Moments for Q = u^5
mom5_u = sum(f1*nodes_gwts*(nodes_u-ubar)**5)/n
mom5_v = sum(f1*nodes_gwts*(nodes_v-vbar)**5)/n
mom5_w = sum(f1*nodes_gwts*(nodes_w-wbar)**5)/n
!!!!!!!!!
! Moments for Q = u^6
mom6_u = sum(f1*nodes_gwts*(nodes_u-ubar)**6)/n
mom6_v = sum(f1*nodes_gwts*(nodes_v-vbar)**6)/n
mom6_w = sum(f1*nodes_gwts*(nodes_w-wbar)**6)/n
!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will calculate the derivatives of the moments and will save them into a whole bunch of scrap variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Notice that the derivatives of the first five invariants must be zero in spatially homogeneous problem. The values of the derivatives 
! tell us how much the conservation laws are violated. 
!!!!!!!!!!!!!!!!!!!!!!!!!
! derivative of mass
dn = sum(frhs1*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!
! derivative of momentum 
!!!!!!!!!!!!!!!!!!!!!!!!!
dubar = (sum(frhs1*nodes_gwts*nodes_u) - dn*ubar)/n
dvbar = (sum(frhs1*nodes_gwts*nodes_v) - dn*vbar)/n
dwbar = (sum(frhs1*nodes_gwts*nodes_w) - dn*wbar)/n
! derivative of temperature 
!!!!!!!!!!!!!!!!!!!!!!!!!
dtemp = ( sum(frhs1*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2) &
              -f1*2.0_DP*((nodes_u-ubar)*dubar+(nodes_v-vbar)*dvbar+(nodes_w-wbar)*dwbar))/3.0_DP*2.0_DP &
                -dn*temp )/n ! derivative of dimensionless temperature
!!!!!!!!!
! BEGIN COMMENTED CODE. ALEX 06022014
! IN THE COMMENTED LINES THE DERIVATIVES of the moments are evaluated 
! Assuming that there are non-trivial derivatives of the first five moments. However, 
! In the exact solution to spatially homogeneous relaxqtion the derivatives of the first five moments are zero
! Nozero derivatives are due exclusively to the numerical errors. We suspect that these numerical errors contaminate 
! the derivatives of moments whose values are near zero??? Odd moments that tend to zero or may be some other moments. 
! Thefore we will comment these lines and replace them with formuls where the derivatives of the first five moments are
! neglected.  
!
! derivatives of directional temperatures 
!!!!!!!!!
!dtemp_u =( sum(frhs1*nodes_gwts*(nodes_u-ubar)**2 - f1*2.0_DP*(nodes_u-ubar)*dubar)/3.0_DP*2.0_DP &
!              -dn*temp_u )/n  ! derivative of dimensionless directional temperature Tx
!dtemp_v =( sum(frhs1*nodes_gwts*(nodes_v-vbar)**2 - f1*2.0_DP*(nodes_v-vbar)*dvbar)/3.0_DP*2.0_DP &
!              -dn*temp_v )/n  ! derivative of dimensionless directional temperature Ty
!dtemp_w =( sum(frhs1*nodes_gwts*(nodes_w-wbar)**2 - f1*2.0_DP*(nodes_w-wbar)*dwbar)/3.0_DP*2.0_DP &
!              -dn*temp_w )/n  ! derivative of dimensionless directional temperature Tz
!!!!!!!!!!
!! Derivatives of Moments for Q = u^3
!dmom3_u = ( sum(frhs1*nodes_gwts*(nodes_u-ubar)**3 - f1*3.0_DP*((nodes_u-ubar)**2)*dubar) - dn*mom3_u )/n
!dmom3_v = ( sum(frhs1*nodes_gwts*(nodes_v-vbar)**3 - f1*3.0_DP*((nodes_v-vbar)**2)*dvbar) - dn*mom3_v )/n
!dmom3_w = ( sum(frhs1*nodes_gwts*(nodes_w-wbar)**3 - f1*3.0_DP*((nodes_w-wbar)**2)*dwbar) - dn*mom3_w )/n
!!!!!!!!!!
!! Derivatives of Moments for Q = u^4
!dmom4_u = ( sum(frhs1*nodes_gwts*(nodes_u-ubar)**4 - f1*4.0_DP*((nodes_u-ubar)**3)*dubar) - dn*mom4_u )/n
!dmom4_v = ( sum(frhs1*nodes_gwts*(nodes_v-vbar)**4 - f1*4.0_DP*((nodes_v-vbar)**3)*dvbar) - dn*mom4_v )/n
!dmom4_w = ( sum(frhs1*nodes_gwts*(nodes_w-wbar)**4 - f1*4.0_DP*((nodes_w-wbar)**3)*dwbar) - dn*mom4_w )/n
!!!!!!!!!
! Moments for Q = u^5
!dmom5_u = ( sum(frhs1*nodes_gwts*(nodes_u-ubar)**5 - f1*5.0_DP*((nodes_u-ubar)**4)*dubar) - dn*mom5_u )/n
!dmom5_v = ( sum(frhs1*nodes_gwts*(nodes_v-vbar)**5 - f1*5.0_DP*((nodes_v-vbar)**4)*dvbar) - dn*mom5_v )/n
!dmom5_w = ( sum(frhs1*nodes_gwts*(nodes_w-wbar)**5 - f1*5.0_DP*((nodes_w-wbar)**4)*dwbar) - dn*mom5_w )/n
!!!!!!!!!
! Moments for Q = u^6
!dmom6_u = ( sum(frhs1*nodes_gwts*(nodes_u-ubar)**6 - f1*6.0_DP*((nodes_u-ubar)**5)*dubar) - dn*mom6_u )/n
!dmom6_v = ( sum(frhs1*nodes_gwts*(nodes_v-vbar)**6 - f1*6.0_DP*((nodes_v-vbar)**5)*dvbar) - dn*mom6_v )/n
!dmom6_w = ( sum(frhs1*nodes_gwts*(nodes_w-wbar)**6 - f1*6.0_DP*((nodes_w-wbar)**5)*dwbar) - dn*mom6_w )/n
!!!!!!!!!
! END COMMENTED LINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! derivatives of directional temperatures 
!!!!!!!!!
dtemp_u =( sum(frhs1*nodes_gwts*(nodes_u-ubar)**2))/3.0_DP*2.0_DP/n  ! derivative of dimensionless directional temperature Tx. changes in the first four moments neglected  
dtemp_v =( sum(frhs1*nodes_gwts*(nodes_v-vbar)**2))/3.0_DP*2.0_DP/n  ! derivative of dimensionless directional temperature Ty. changes in the first four moments neglected
dtemp_w =( sum(frhs1*nodes_gwts*(nodes_w-wbar)**2))/3.0_DP*2.0_DP/n  ! derivative of dimensionless directional temperature Tz. changes in the first four moments neglected
!!!!!!!!!
! Derivatives of Moments for Q = u^3
dmom3_u = ( sum(frhs1*nodes_gwts*(nodes_u-ubar)**3))/n ! changes in the first four moments neglected
dmom3_v = ( sum(frhs1*nodes_gwts*(nodes_v-vbar)**3))/n ! changes in the first four moments neglected
dmom3_w = ( sum(frhs1*nodes_gwts*(nodes_w-wbar)**3))/n ! changes in the first four moments neglected
!!!!!!!!!
! Derivatives of Moments for Q = u^4
dmom4_u = ( sum(frhs1*nodes_gwts*(nodes_u-ubar)**4))/n ! changes in the first four moments neglected
dmom4_v = ( sum(frhs1*nodes_gwts*(nodes_v-vbar)**4))/n ! changes in the first four moments neglected
dmom4_w = ( sum(frhs1*nodes_gwts*(nodes_w-wbar)**4))/n ! changes in the first four moments neglected
!!!!!!!!!
! Moments for Q = u^5
dmom5_u = ( sum(frhs1*nodes_gwts*(nodes_u-ubar)**5))/n ! changes in the first four moments neglected
dmom5_v = ( sum(frhs1*nodes_gwts*(nodes_v-vbar)**5))/n ! changes in the first four moments neglected
dmom5_w = ( sum(frhs1*nodes_gwts*(nodes_w-wbar)**5))/n ! changes in the first four moments neglected
!!!!!!!!!
! Moments for Q = u^6
dmom6_u = ( sum(frhs1*nodes_gwts*(nodes_u-ubar)**6))/n ! changes in the first four moments neglected
dmom6_v = ( sum(frhs1*nodes_gwts*(nodes_v-vbar)**6))/n ! changes in the first four moments neglected
dmom6_w = ( sum(frhs1*nodes_gwts*(nodes_w-wbar)**6))/n ! changes in the first four moments neglected
!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will calculate the local maxwellian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fm = maxwelveldist(temp,ubar,vbar,wbar,n,nodes_u,nodes_v,nodes_w) ! now we have the maxwellian with the same macroparamters.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we compute the difference in the moments between the function and the local maxwellian, we will need them  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODIFIED ALEX: 06/02/2014
!! explicitly will compute the moments of the steady state it and subtract them. Previously, momements of f1-fm were computed
!!!!!!!!!!!!!!!!!!!!!!!!!
! difference in mass
mn=sum((f1-fm)*nodes_gwts) ! this one is a measure of how accurate is the evaluation of macroparameters.
!!!!!!!!!!!!!!!!!!!!!!!!!
! check momentum 
!!!!!!!!!!!!!!!!!!!!!!!!!
mubar=sum((f1-fm)*nodes_gwts*nodes_u)/n  !! this one is a measure of how accurate is the evaluation of macroparameters.
mvbar=sum((f1-fm)*nodes_gwts*nodes_v)/n  ! this one is a measure of how accurate is the evaluation of macroparameters.
mwbar=sum((f1-fm)*nodes_gwts*nodes_w)/n  ! this one is a measure of how accurate is the evaluation of macroparameters.
! check temperature 
!!!!!!!!!!!!!!!!!!!!!!!!!
mtemp = sum((f1-fm)*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature ! this one is a measure of how accurate is the evaluation of macroparameters.
!!!!!!!!!
! check directional temperatures 
!!!!!!!!!
mtemp_u = sum(f1*nodes_gwts*(nodes_u-ubar)**2)/n/3.0_DP*2.0_DP - temp/3.0_DP ! dimensionless directional temperature Tx
mtemp_v = sum(f1*nodes_gwts*(nodes_v-vbar)**2)/n/3.0_DP*2.0_DP - temp/3.0_DP ! dimensionless directional temperature Ty
mtemp_w = sum(f1*nodes_gwts*(nodes_w-wbar)**2)/n/3.0_DP*2.0_DP - temp/3.0_DP! dimensionless directional temperature Tz
!!!!!!!!!
! Moments for Q = u^3 ! try without fm -- the corresponding moment is zero anyway... 
mmom3_u = sum(f1*nodes_gwts*(nodes_u-ubar)**3)/n
mmom3_v = sum(f1*nodes_gwts*(nodes_v-vbar)**3)/n
mmom3_w = sum(f1*nodes_gwts*(nodes_w-wbar)**3)/n
!!!!!!!!!
! Moments for Q = u^4
mmom4_u = sum(f1*nodes_gwts*(nodes_u-ubar)**4)/n - 3.0_DP*temp**2/4.0_DP
mmom4_v = sum(f1*nodes_gwts*(nodes_v-vbar)**4)/n - 3.0_DP*temp**2/4.0_DP
mmom4_w = sum(f1*nodes_gwts*(nodes_w-wbar)**4)/n - 3.0_DP*temp**2/4.0_DP 
!!!!!!!!!
! Moments for Q = u^5 ! try without fm -- the corresponding moment is zero anyway
mmom5_u = sum(f1*nodes_gwts*(nodes_u-ubar)**5)/n 
mmom5_v = sum(f1*nodes_gwts*(nodes_v-vbar)**5)/n
mmom5_w = sum(f1*nodes_gwts*(nodes_w-wbar)**5)/n
!!!!!!!!!
! Moments for Q = u^6
mmom6_u = sum(f1*nodes_gwts*(nodes_u-ubar)**6)/n - 15.0_DP*temp**3/8.0_DP
mmom6_v = sum(f1*nodes_gwts*(nodes_v-vbar)**6)/n - 15.0_DP*temp**3/8.0_DP
mmom6_w = sum(f1*nodes_gwts*(nodes_w-wbar)**6)/n - 15.0_DP*temp**3/8.0_DP
!!!!!!!!!

!! END modified Alex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we compute the relaxation times/speeds for the moments
! the formula for the relaxation times is
! 
! \nu_{\phi} = -\partial_{t} \ln |f_{\phi} - f_{\phi}^{M}| , \partial_{t} f_{\phi}^{M}=0
!
! thus
! 
! \nu_{\phi} = -\partial_{t} f_{\phi}/(f_{\phi} - f_{\phi}^{M}) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
tresh=abs(dn*epps) ! this is treshhold for sensitivity of the evaluation of relaxation times. 
! if the difference in $f_{\varphi}-f^{M}_{\varphi}$ falls below this treshhold, then 
! we values on the relaxation time probably are not reliable. in this case the relaxation time are equated to zero... 
!!!
if (abs(mn) > tresh) then 
 taun = dn/mn
else 
 taun = 0.0_DP
end if 
!!! 
if (abs(mubar) > tresh) then 
 tauubar = -dubar/mubar 
else 
 tauubar = 0.0_DP
end if 
if (abs(mvbar) > tresh) then 
 tauvbar = -dvbar/mvbar 
else 
 tauvbar = 0.0_DP
end if 
if (abs(mwbar) > tresh) then 
 tauwbar = -dwbar/mwbar 
else 
 tauwbar = 0.0_DP
end if 
!!!
if (abs(mtemp) > tresh) then 
 tautemp = -dtemp/mtemp 
else 
 tautemp = 0.0_DP
end if 
!!!
if (abs(mtemp_u) > tresh) then 
 tautemp_u = -dtemp_u/mtemp_u 
else 
 tautemp_u = 0.0_DP
end if 
if (abs(mtemp_v) > tresh) then 
 tautemp_v = -dtemp_v/mtemp_v 
else 
 tautemp_v = 0.0_DP
end if 
if (abs(mtemp_w) > tresh) then 
 tautemp_w = -dtemp_w/mtemp_w 
else 
 tautemp_w = 0.0_DP
end if 
!!!
if (abs(mmom3_u) > tresh) then 
 taumom3_u = -dmom3_u/mmom3_u 
else 
 taumom3_u = 0.0_DP
end if 
if (abs(mmom3_v) > tresh) then 
 taumom3_v = -dmom3_v/mmom3_v 
else 
 taumom3_v = 0.0_DP
end if 
if (abs(mmom3_w) > tresh) then 
 taumom3_w = -dmom3_w/mmom3_w 
else 
 taumom3_w = 0.0_DP
end if 
!!!
if (abs(mmom4_u) > tresh) then 
 taumom4_u = -dmom4_u/mmom4_u 
else 
 taumom4_u = 0.0_DP
end if 
if (abs(mmom4_v) > tresh) then 
 taumom4_v = -dmom4_v/mmom4_v 
else 
 taumom4_v = 0.0_DP
end if 
if (abs(mmom4_w) > tresh) then 
 taumom4_w = -dmom4_w/mmom4_w 
else 
 taumom4_w = 0.0_DP
end if 
!!! 
if (abs(mmom5_u) > tresh) then 
 taumom5_u = -dmom5_u/mmom5_u 
else 
 taumom5_u = 0.0_DP
end if 
if (abs(mmom5_v) > tresh) then 
 taumom5_v = -dmom5_v/mmom5_v 
else 
 taumom5_v = 0.0_DP
end if 
if (abs(mmom5_w) > tresh) then 
 taumom5_w = -dmom5_w/mmom5_w 
else 
 taumom5_w = 0.0_DP
end if 
!!!
if (abs(mmom6_u) > tresh) then 
 taumom6_u = -dmom6_u/mmom6_u 
else 
 taumom6_u = 0.0_DP
end if 
if (abs(mmom6_v) > tresh) then 
 taumom6_v = -dmom6_v/mmom6_v 
else 
 taumom6_v = 0.0_DP
end if 
if (abs(mmom6_w) > tresh) then 
 taumom6_w = -dmom6_w/mmom6_w 
else 
 taumom6_w = 0.0_DP
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Next we will need to put the 
!! calculated quatities in the storage array... 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i = int(MomDerivRelaxT(1)) ! Number of records
k = int(MomDerivRelaxT(2)) ! number of columns in the record
j = 2+i*k ! now the index j points to the last filled cell. The first new record can be placed in the following cell.
! A quick check if there is enought room for storing the computed values. 
if (j+k > size(MomDerivRelaxT,1)) then 
 print *, "RecordMomDerivRelaxTimes_DGV: Not enough memory in (MomDerivRelaxT) to store moments derivs. and relax. times."
else 
 MomDerivRelaxT(1) = Real(i)+1.0    
 !!! Time
 MomDerivRelaxT(j+1)=curr_time
 !!!! First will go the derivatives of macroparameters. 
 MomDerivRelaxT(j+2)=dn
 !!!
 MomDerivRelaxT(j+3)=dubar
 MomDerivRelaxT(j+4)=dvbar
 MomDerivRelaxT(j+5)=dwbar
 !!!
 MomDerivRelaxT(j+6)=dtemp
 !!!
 MomDerivRelaxT(j+7)=dtemp_u
 MomDerivRelaxT(j+8)=dtemp_v
 MomDerivRelaxT(j+9)=dtemp_w
 !!!
 MomDerivRelaxT(j+10)=dmom3_u
 MomDerivRelaxT(j+11)=dmom3_v
 MomDerivRelaxT(j+12)=dmom3_w
 !!!
 MomDerivRelaxT(j+13)=dmom4_u
 MomDerivRelaxT(j+14)=dmom4_v
 MomDerivRelaxT(j+15)=dmom4_w
 !!!
 MomDerivRelaxT(j+16)=dmom5_u
 MomDerivRelaxT(j+17)=dmom5_v
 MomDerivRelaxT(j+18)=dmom5_w
 !!!
 MomDerivRelaxT(j+19)=dmom6_u
 MomDerivRelaxT(j+20)=dmom6_v
 MomDerivRelaxT(j+21)=dmom6_w
 !!!
 !! Next go the frequences for the same moments
 !!!
 MomDerivRelaxT(j+22)=taun
 !!!
 MomDerivRelaxT(j+23)=tauubar
 MomDerivRelaxT(j+24)=tauvbar
 MomDerivRelaxT(j+25)=tauwbar
 !!!
 MomDerivRelaxT(j+26)=tautemp
 !!!
 MomDerivRelaxT(j+27)=tautemp_u
 MomDerivRelaxT(j+28)=tautemp_v
 MomDerivRelaxT(j+29)=tautemp_w
 !!!
 MomDerivRelaxT(j+30)=taumom3_u
 MomDerivRelaxT(j+31)=taumom3_v
 MomDerivRelaxT(j+32)=taumom3_w
 !!!
 MomDerivRelaxT(j+33)=taumom4_u
 MomDerivRelaxT(j+34)=taumom4_v
 MomDerivRelaxT(j+35)=taumom4_w
 !!!
 MomDerivRelaxT(j+36)=taumom5_u
 MomDerivRelaxT(j+37)=taumom5_v
 MomDerivRelaxT(j+38)=taumom5_w
 !!!
 MomDerivRelaxT(j+39)=taumom6_u
 MomDerivRelaxT(j+40)=taumom6_v
 MomDerivRelaxT(j+41)=taumom6_w
 !!!!!!!!!!!!!!!!!!!!
end if  

!!! if the time(current time) is bigger than the final time t_R, then the program will deallocate the array -- this is a cleanup call
if (curr_time > t_R) then
 deallocate(MomDerivRelaxT)
 MomDerivRelaxT_flag = 0 
end if    
!!!!!!!!!!!!
end subroutine RecordMomDerivRelaxTimes_DGV



   
end module dgvtools_mod
