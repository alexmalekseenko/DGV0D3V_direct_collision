!
! collision_mod.f90
!
! in this module, the subroutines involved in the evaluation of the collision integral are listed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module collision_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none



contains 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicA_DGV (f,fc) 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated on one cell. The values on other cells are
! obtained by a transposition from operator A on the canonical cell --- which should be selected at the 
! center fot he mesh (or close to the center)
! 
! STILL NEED TO ADD THE COLLISION FREQUENCY...
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicA_DGV(f,fc)

use commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: k,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
real (DP), dimension(:), allocatable :: Aphi ! scrap array to keep the portion of A
real (DP), dimension (:), allocatable :: fxi,fxi1 !scrap array to keep the f(xi) and f(xi1) 
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift = nodes_Ashift(iphi)
   phicap = A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   allocate(fxi(1:phicap),fxi1(1:phicap),Aphi(1:phicap),stat=loc_alloc_stat)
    !
   if (loc_alloc_stat >0) then 
   print *, "EvalCollisionPeriodicA_DGV: Allocation error for variables (fxi,fxi1,Aphi)"
   stop
   end if
   !!!!!!!!!!!!!!!!!
   fxi=0;fxi1=0 ! nullify the arrays
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   k=0 ! counter -- will be used to compress the zero records in arrays... 
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             k=k+1
             fxi(k)=f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) )
           !  zz=((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
           !   (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1)
             fxi1(k)=f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) )
             Aphi(k)=A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do   
   !! the arrays are set. Proceed to evaluate the matrix product   
   fc(iphi)=2*sum(Aphi(1:k)*fxi(1:k)*fxi1(1:k))/nodes_gwts(iphi)     
   deallocate (fxi,fxi1,Aphi)   
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicA_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicA_DGV (f,fc) 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated on one cell. The values on other cells are
! obtained by a transposition from operator A on the canonical cell --- which should be selected at the 
! center fot he mesh (or close to the center)
! 
!
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicAPlus_DGV(f,fc)

use commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             fc(iphi)=fc(iphi)+2*f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) )*&
                               f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) )*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=fc(iphi)/nodes_gwts(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicAPlus_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicMixedTermsA_DGV (f,fm,fc) 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated on one cell. The values on other cells are
! obtained by a transposition from operator A on the canonical cell --- which should be selected at the 
! center fot he mesh (or close to the center)
! 
! This SUBROUTINE IS TO BE USED IN THE DECOMPOSITION MODE when storing derivative is not 
! feasible. It evaluates the cross term \int\int f_{i} fm_{j} A^{ij}_{k} 
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicMixedTermsA_DGV(f,fm,fc)

use commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution-maxwellian at the current time step. 
real (DP), dimension (:), intent (in) :: fm ! the components of the maxwellian at the current time step.
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) 
             fc(iphi)=fc(iphi)+( fm(xi_j)*f(xi1_j) + f(xi_j)*fm(xi1_j) )*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=2*fc(iphi)/nodes_gwts(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicMixedTermsA_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionLinear(f,fc) 
!
! This subroutine evaluates the linearized collision kernel. To call this subroutine, 
! a vector of the linearized kernel must be prepared. the vector is 2\int f_{m}(v)A(v,v_1,phi) dv
! it is calcuated in the subroutine PrepareFMA_DGV
! 
! 
! Here the arrays fmA has the linearizised operator in it
! accessed directly from the commvar
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine EvalCollisionLinear(f,fc)        ! This evaluates the linear part

use commvar, only: fmA

intrinsic MATMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fc = 2*MATMUL(f,fmA) ! this is a contraction with the gradient , i.e., the derivative of the collisio operator at some maxwellian

end subroutine EvalCollisionLinear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine PrepareFMA_DGV (fm)
!
!
! This subroutine prepares the derivative of the collision integral, 2\int f_{M}(v) A(v,v_1) dv
! 
! fm is the current local maxwellian
! 
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated on one cell. The values on other cells are
! obtained by a transposition from operator A on the canonical cell --- which should be selected at the 
! center fot he mesh (or close to the center)
!
! Before calling this function, make sure that  Nodes_Ashift, Nodes_ccan and other supplementary arrays are set before calling the 
! the evaluation of the linearized collision operator.

subroutine PrepareFMA_DGV
use commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,fm,fmA,nodes_gwts

!!!!!!!!!!! ALL IMPORTANT VARIABLES ARE ACCESSED DIRECTLY !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the linearized collision operator: fmA=\int_{R^3} f_{M}(v)A(v,v_1,\phi)dv.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fmA=0 ! nullify the result before computing anything... 
!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) 
             fmA(xi_j,iphi)  = fmA(xi_j,iphi)  + fm(xi1_j)*A(Ashift+j)
             fmA(xi1_j,iphi) = fmA(xi1_j,iphi) + fm(xi_j)*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do   
   !! Need to divide by the volume element for the node with the number iphi 
   !! this division appears in the formula for the collision operator
   fmA(:,iphi)=fmA(:,iphi)/nodes_gwts(iphi)
   !! the arrays are set. Proceed to evaluate the linearized solution  
end do ! End of the main loop in nodal points
!
end subroutine PrepareFMA_DGV 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicA_MPI_DGV (f,fc,irank) 
!
! This is an MPI analog of the subroutine. On the master processor, only the collective operation is called. On the slave processor, 
! the evaluation of the portion of the collision operator is performed. 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated on one cell. The values on other cells are
! obtained by a transposition from operator A on the canonical cell --- which should be selected at the 
! center fot he mesh (or close to the center)
! 
!
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicAPlus_MPI_DGV(f,fc)

use commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts,procs_nodes_wld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
integer :: irank ! the number of the process on which the software is running
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (DP), dimension (:), allocatable :: fc_buff ! 
integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: i,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status
! 
integer :: ierr ! variables for MPI Calls

!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
do i=1,size(procs_nodes_wld,1)  ! Loop in the nodal points
   iphi=procs_nodes_wld(i)
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             fc(iphi)=fc(iphi)+2*f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) )*&
                               f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) )*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=fc(iphi)/nodes_gwts(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicAPlus_MPI_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine PrepareFMA_MPI_DGV (fm)
!
! this is the MPI extension of the above subroutine
!
! This subroutine prepares the derivative of the collision integral, 2\int f_{M}(v) A(v,v_1) dv
! 
! fm is the current local maxwellian
! 
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated on one cell. The values on other cells are
! obtained by a transposition from operator A on the canonical cell --- which should be selected at the 
! center fot he mesh (or close to the center)
!
! Before calling this function, make sure that  Nodes_Ashift, Nodes_ccan and other supplementary arrays are set before calling the 
! the evaluation of the linearized collision operator.
!
! It is assumed that the fmA array of the right size exists...
!
subroutine PrepareFMA_MPI_DGV(irank)
use commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,fm,fmA,nodes_gwts,procs_nodes_wld,&
                   fmA_flag

!!!!!!!!!!! ALL IMPORTANT VARIABLES ARE ACCESSED DIRECTLY !!!!!!!!!!!!!

integer, intent (in):: irank ! the number of the processor in the MPI universe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j,i! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the linearized collision operator: fmA=\int_{R^3} f_{M}(v)A(v,v_1,\phi)dv.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check if fmA have a wrong size then deallocate it and allocate with the right size
if (fmA_flag==0) then 
 allocate (fmA(1:size(nodes_ui,1),1:size(procs_nodes_wld,1)), stat=loc_alloc_stat)
 !
 if (loc_alloc_stat >0) then 
  print *, "PrepareFMA_MPI_DGV: Allocation error for variables (fmA) on proicess", irank, ". Stop."
  stop
 end if
 fmA_flag = 1  ! this falg is on if the array FmA was allocated
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fmA=0 ! nullify the result before computing anything... 
!
!
do i=1,size(procs_nodes_wld,1)  ! Loop in the nodal points
   iphi=procs_nodes_wld(i)
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) 
             fmA(xi_j,i)  = fmA(xi_j,i)  + fm(xi1_j)*A(Ashift+j)
             fmA(xi1_j,i) = fmA(xi1_j,i) + fm(xi_j)*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do   
   !! Need to divide by the volume element for the node with the number iphi 
   !! this division appears in the formula for the collision operator
   fmA(:,i)=fmA(:,i)/nodes_gwts(iphi)
   !! the arrays are set. Proceed to evaluate the linearized solution  
end do ! End of the main loop in nodal points
!
end subroutine PrepareFMA_MPI_DGV 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine PrepConsolLinColl_highperf_DGV(irank,Acopy_irank)
!
! this is the MPI subroutine to prepare the linearized operator on the group of 
! processors that will deal with the evaluation of the collision operator
!
! This subroutine is similar to PrepareFMA_MPI_DGV. The derivative is compiuted for the nodes given in the 
! array procs_nodes_wld. 
! 
! Once the components are computed -- they are sent out to the master processor of the corresponding Acopy universe
! 
! This subroutine prepares the derivative of the collision integral, 2\int f_{M}(v) A(v,v_1) dv
! 
! fm is the current local maxwellian
! 
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated on one cell. The values on other cells are
! obtained by a transposition from operator A on the canonical cell --- which should be selected at the 
! center fot he mesh (or close to the center)
!
! Before calling this function, make sure that  Nodes_Ashift, Nodes_ccan and other supplementary arrays are set before calling the 
! the evaluation of the linearized collision operator.
!
! It is assumed that the fmA array of the right size exists...
!
!
! !! The new algorithm for evaluatin the collision operator is as follows:
 !!
 !! A. Slave Processors  compute the components of the linearized operator.
 !! B. Slave processors place these componets in the copy of an array 
 !!    FmA whose rows go over all nodes in the array Acopy_Wrkld that are the same for 
 !!    all processors in the same Acopy universe
 !!    to help them stack the components the array procs_nodes_wkld_Acopy_Addr array us used. 
 !!    procs_nodes_wld_Acopy_addr has the same amount of elements as procs_nodes_wld
 !!    The array procs_nodes_wld_Acopy_addr contains indices of the nodes from 
 !!    procs_nodes_wld in Acopy_Wrkld
 !!    The lead processor of Acopy univ. gathers the FmA to it
 !! 
 !! C. The lead proc of Acopy universe is also the lead processor of the matching LinFmA universe 
 !!    This processor broadcasts FmA to processors involved in the valuation of the collision oeprator. 
 !!    For each Acopy universe there is a LinFmA_Univ universe that combines all processors that need 
 !!    to receive data from processors invloved in Acopy_univ. 
 !! D. The processors pick the data from the boradsact and build the local copy of the linearized operator.
 !!    To do that three arrays are used:
 !!    lin_proc_nodes_wld --- contains the componets of the linearized operator assigned to each processor involved 
 !!    in the evaluation of the linearized operator. 
 !!    linprocndswld_Acopy_univs, --- this array contains informaion about arriving components of the linearized operator:
 !!     # of the universe A copy that will send the components
 !!     # total number of the nodes (each node is one row of the lineaired collision operator
 !!     # Nodes .... 
 !!    then repeat
 !!    linprocndswld_Acopy_addr has the same stucture except instead of nodes, there will be their indices in the corresponding 
 !!    Acopy_Wrkld array
 !!
 !! E. The processors involced in the evaluation of the linearized colliion operator create a copy of BFmA where they store the 
 !!    componets of the consolidated operator.  
 !!    the processors evaluate the parts of the consolidated collision operator and stack the results in the temporary solution 
 !!    array. Finally master processrs gathers all componets of the temp. solution array by a collective communication call. 
 !! 
 !!!!!!!!!!!!!!!!!!!!!
!
!!! the order in which the LinFmA universes exchange the components is determined in the array LinFmA_Univs_callorder.
!! this array contains the sequence in which the LinFmA universes are called
!! It contains numbers from 1 to num_Acopies that are permuted in some way
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

subroutine PrepConsolLinColl_highperf_DGV(irank)         
use commvar, only: fmA,fmA_flag,BfmA,BfmA_flag,&
                   procs_nodes_wld,lin_proc_nodes_wld,Acopy_wrkld,& 
                   procs_nodes_wld_Acopy_addr,LinFmA_Univs_callorder,&
                   num_lin_proc,nodes_u,num_Acopies,nodes_proc_groups,&
                   MPI_Acopy_UNIVS,My_Acopy_Univ_indx,&
                   !MPI_Acopy_UNIVS_GROUP, MPI_LinFmA_UNIVS_GROUP,& ! theese do not seem to be needed either...
                   MPI_LinFmA_UNIVS,linprocndswld_Acopy_univs, linprocndswld_Acopy_addr,&
                   linprocndswld_BfmA_addr,Acopy_irank

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!! 

!!!!!!!!!!! ALL IMPORTANT VARIABLES ARE ACCESSED DIRECTLY !!!!!!!!!!!!!

integer, intent (in):: irank ! the number of the processor in the MPI universe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(DP), dimension(:), allocatable :: buffer
!!!1
integer (I4B) :: i,j,ii,jj,nn,jjj,buff_len,BfmA_addr,Acopy_addr !indices
integer (I4B) :: proc, nnds, send_node, active_Lin_FmA,send_Lin_FmA! useful variables 
integer :: loc_alloc_stat ! variable to keep the allocation status

logical :: is_in_Acopy_univ, is_Acopy_master ! two flags to signal that the Lin_FmA/Acopy universe is in the list of sendres and that the current  processor is the master of the next Lin_FmA

real (DP), dimension (:,:), allocatable :: Acopy_FmA,LinUniv_FmA! this is to keep the linearized operator, BfmA -- is the consolidated linearized operator
integer :: Acopy_FmA_flag=0, LinUniv_FmA_flag=0 ! these are flags if the arrays are allocated -- the size function sometimes messes up..
!! MPI VARIABLES:
integer, dimension (MPI_STATUS_SIZE) :: istatus
integer :: ierr,sendproc ! status variable 
integer :: ixx, root ! the value of the calling processor
integer :: Acopy_UNIV_GROUP,Acopy_UNIV,LinFmA_UNIV_GROUP,LinFmA_UNIV ! scrap pointers to Acopy and LinFmA universes

!! END MPI VARIABLES

!!! need this for consistency check:
integer :: My_FmA_rank  ! rank of the processor in an FmA_universe
!!! end of the consistency check variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Effectively, we should not evaluate any new derivatives. We just need to send out the ones that are already computed. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (irank < num_lin_proc) then ! the consolidated lineariez operator is only created on these processes   
  ! we check if the consolidated linearized operator has not been created yet, or may be it was deallocated 
  ! if it is allocated, it must have the correct size. 
  if (BfmA_flag == 0) then
   BfmA_flag=1 ! set the flag to the allocation status 1 - allocated
   nn = size(nodes_u,1) ! the total number of nodes in velocity
   allocate (BfmA(1:nn,1:size(lin_proc_nodes_wld,1)), stat=loc_alloc_stat)
   if (loc_alloc_stat>0) then 
    print *, "PrepConsolLinColl_highperf_DGV: Allocation error for variables (BfmA) on process", irank, ". Stop."
    stop
   end if
  end if
  ! 
  BfmA=0 ! reset the consolidated lienarized array.
  !   
end if  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! During Step 1. all slave processes send the componets of the linearized operator to the master processor of their Acopy universe
!! 
!!!!!!
if (irank > 0) then ! Only do it on the slave processors
 ! first, we create the Acopy_FmA array which consolidates the components of the linearized operator 
 ! inside the Acopy Universe. 
 nn = size(nodes_u,1) ! the total number of nodes in velocity
 ! actually, this array should only be allocated once in the call of the subroutine
 ! so the next check is not needed.
 if (Acopy_FmA_flag == 1) then 
  deallocate(Acopy_FmA)
  Acopy_FmA_flag = 0
 end if 
 ! create the Acopy_FmA arrat
 allocate (Acopy_FmA(1:nn,1:size(Acopy_Wrkld,1)), stat=loc_alloc_stat)
 if (loc_alloc_stat>0) then 
  print *, "PrepConsolLinColl_highperf_DGV: Allocation error for variables (Acopy_FmA) on process", irank, ". Stop."
  stop
 end if
 Acopy_FmA_flag=1 ! set the flag to the allocation status 1 - allocated
 Acopy_FmA=0 ! nullify newly created array. It is important that we have it initially zero.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! Now let us prepare the Acopy_FmA for collective operation.
 !! We record components of the linearized operator into the local copy of Acopy_FmA array
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do i = 1, size(procs_nodes_wld,1)
   j = procs_nodes_wld_Acopy_addr(i) ! determine the index of this components in the array Acopy_FmA (rows correspond to order of nodes in Acopy_Wrkld)
   do ii=1,nn
    Acopy_FmA(ii,j)= fmA(ii,i) ! this will take components of the linearized operator that are stored on this processor and move them to the consolidated operator in this Acopy universe 
   end do 
 end do 
 !
 ! we release the memory of the fmA array. There should not be such an array on irank=0, but there will be 
 ! on other processors...  copy the contents of the scrap array into the fmA array. fmA need to be re-allocated to
 deallocate (fmA)
 fmA_flag = 0 ! 0 means that the array fmA is not allocated
 !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
 !! The array Acopy_FmA is now ready for consolidation on the master processor of this Acopy Universe.
 ! shound not need this: Acopy_UNIV_GROUP = MPI_Acopy_UNIVS_GROUP(My_Acopy_Univ_indx) !set up pointers to communicator group 
 Acopy_UNIV = MPI_Acopy_UNIVS(My_Acopy_Univ_indx)             ! and the communicator
 !!!
 ii = size(Acopy_Wrkld,1)
 buff_len = nn*ii 
 ! set up buffer on the processor with Acopy_irenk =0
 ! on other processors the buffer is ignored in the collective operation:
 if (Acopy_irank==0) then 
  allocate (buffer(1:buff_len), stat=loc_alloc_stat)  !
  if (loc_alloc_stat >0) then 
   print *, "PrepConsolLinColl_highperf_DGV: Allocation error for variable (buffer). Process", irank, ". Stop"
   stop
  end if
 end if
 !!  make a collective operation call in the local Acopy universe. 
 call mpi_reduce (Acopy_FmA,buffer,buff_len,MPI_DOUBLE_PRECISION,MPI_SUM,0,Acopy_UNIV,ierr) !
 if (ierr /= 0 ) then 
  print *,"PrepConsolLinColl_highperf_DGV: slave processor", irank, "MPI_REDUCE to proc 0 returned error", ierr
  stop
 end if
 !!  now the master processor on the Acopy_Univ has the consolidated Acopy_FmA in the buffer. we need to copy it in the Acopy_FmA
 if (Acopy_irank==0) then 
  do j=1,ii
   do i=1,nn
    Acopy_FmA(i,j) = buffer((j-1)*nn+i) ! should not there be a simpler way to do it? 
   end do 
  end do 
  deallocate (buffer) ! release the memory used by the buffer array
 else 
  deallocate(Acopy_FmA)
  Acopy_FmA_flag = 0
 end if      
 !!! Now Acopy_FmA is ready on the master nodes of all Acopy universes. 
end if
!! We are ready to start send the components of the linearized operator out.... 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! The componets of the linearized operator will be broadcasted to LinFmA universes.
!!! 
!!! The universes will take turns broadcasting the data
!!!
!!! The order in which broadcasts will be called is determined by the array LinFmA_Univs_callorder
!!!
!!! Each process will go down the array LinFmA_Univs_callorder and match its entries (nubers of the FmA universes) 
!!! against the records in the local array linprocndswld_Acopy_univs
!!! if a universe is found in linprocndswld_Acopy_univs with the current call number, than the process enters into collective communication 
!!! for this universe. 
!!! if the record for this universe is not found, then we still check if the process My_Acopy_indx is the same as the universe in question.
!!! If it is, than we check if Acopy_irank of this processor is zero. This would mean that this process is a master node of 
!!! the Acopy universe that needs to be communicating next. In this case, the process performs the collective communcation: the 
!!! broadcast of Acopy_FmA to its LinFmA universe.
!!! 
!!! Now it may be possible to both find the universe number in linprocndswld_Acopy_univs and to have 
!!! My_Acopy index to be equal to the number and Acopy_irank=0. In this case the process needs to send the data to itself and 
!!! the flag Lin_Self_Send_FmA_Flag should be set to 1 by the preparatory subroutine. 
!!! 
!!! The difference between just sending and just recieving (via boradcast) is that the array Acopy_FmA is already setup in the processor 
!!! with Acopy_irank=0. While on the receiveing processors it needs to be created. 
!!!
!!! It is is self-send then there is no need to create the array Acopy_FmA. However, after the broadcast, we need to get out the components of the Acopy_FmA.
!!! 
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

do i=1,num_Acopies ! run over all universes 
 active_Lin_FmA = LinFmA_Univs_callorder(i) ! determine which Lin_FmA universe is speaking next 
 ! resent the flags
 is_in_Acopy_univ = .false. 
 is_Acopy_master = .false.
 ! now we test the linprocndswld_Acopy_univs if the active_LinFmA is present there. 
 ! we only may have meaningful arrays on the processors participating in the evaluation of the 
 ! consolidated linear operator
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (irank < num_lin_proc) then ! so if the flag is_in_Acopy_univ is true, then irank<num_lin_proc
  j=1
  do while (j<size(linprocndswld_Acopy_univs,1))
   send_Lin_FmA = linprocndswld_Acopy_univs(j) ! pick a Lin_FmA universe that sends staff
   nnds = linprocndswld_Acopy_univs(j+1)
   if (send_Lin_FmA == active_Lin_FmA) then 
    is_in_Acopy_univ = .true.  ! the process receives data from the active Acopy_universe.
    exit
   end if  
   j=j+nnds+2
  end do 
 end if 
 ! now let us check if this process is the master node of the active Lin_FmA universe
 if ((active_Lin_FmA == My_Acopy_Univ_indx) .and. (Acopy_irank == 0))  then  !My_Acopy_Univ_indx=0 on the master prociessor, Acopy_irank=-1 on the master processor
  is_Acopy_master = .true. ! this processor is the master of the active Lin_fmA_Universe 
 end if 
 ! now we need to do the send-recieve call. 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! according the to values of flags is_in_Acopy_univ and is_Acopy_master, different things are performed.
 ! if is_in_Acopy_univ = true and is_Acopy_master = false, then we need to create LinUniv_FmA array of proper size
 ! then call the communication subroutine. After the broadcast, we need to pick out a few components out of this array
 !
 ! if is_in_Acopy_univ = true and is_Acopy_master = true, we have Acopy_FmA array ready on this processor, but we still need 
 ! to broadcast it. After (or even before) the broadcast, we need to pick out a few components out of this array
 !
 ! if is_in_Acopy_univ = false and is_Acopy_master = true, then we need to broadcast the Acopy_FmA array and this will be it. 
 !  
 ! if is_in_Acopy_univ = false and is_Acopy_master = false, we mode on to the different sending universe,
 ! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if ((is_in_Acopy_univ == .true.) .and. (is_Acopy_master == .false.)) then 
   ! We will now prepare the temporary array for the communication of the componets of linearized operator.
   ! it should not be necessary, but in case the array already exists, we need to destroy it 
   if (LinUniv_FmA_flag == 1) then 
    deallocate(LinUniv_FmA)
    LinUniv_FmA_flag = 0
   end if 
   ! we need to allocate the array LinUniv_FmA_flag to receive components from the active Lin_FmA iniverse
   nn = size(nodes_u,1)
   nnds = nodes_proc_groups(2,active_Lin_FmA) - nodes_proc_groups(1,active_Lin_FmA) + 1
   allocate (LinUniv_FmA(1:nn,1:nnds), stat=loc_alloc_stat)
   if (loc_alloc_stat>0) then 
    print *, "PrepConsolLinColl_highperf_DGV: Allocation error 1 for (Acopy_FmA) on process", irank, ". Stop."
    stop
   end if
   LinUniv_FmA_flag = 1 ! set the flag to the allocation status 1 - allocated
   LinUniv_FmA = 0 ! nullify newly created array. It is important that we have it initially zero. 
 end if 
 ! now we are ready to send of receive or both 
 ! this is the situation of self-send.... 
 if ((is_Acopy_master == .true.) .or. (is_in_Acopy_univ == .true.)) then ! in all three cases, there will be a boradcast: need to prepare:    
   ! should not need this ...: LinFmA_Univ_GROUP = MPI_LinFmA_UNIVS_GROUP(active_Lin_FmA) !set up pointers to the next communicator group 
   LinFmA_Univ = MPI_LinFmA_UNIVS(active_Lin_FmA)             ! and the communicator 
   nn = size(nodes_u,1)
   nnds = nodes_proc_groups(2,active_Lin_FmA) - nodes_proc_groups(1,active_Lin_FmA) + 1
   nnds=nnds*nn ! this is the size of the boradcast....
 end if 
 ! Next we see if we need to do a broadcast to send Acopy_FmA
 if (is_Acopy_master == .true.) then ! this is the situation of a send including a self-send. 
  call mpi_bcast (Acopy_FmA,nnds,MPI_DOUBLE_PRECISION,0,LinFmA_Univ,ierr)  ! boradcast Acopy_FmA from the maste processor of an active Lin_FmA Universe
  if (ierr /= 0 ) then 
   print *,"PrepConsolLinColl_highperf_DGV: bcast of Acopy_FmA from ", irank, "Acopy", My_Acopy_Univ_indx," returned error. stop", ierr
   stop
  end if
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 
  !! A consistency check -- can comment later -- the master processor of an Acopy universe also much be a master processor onf 
  My_FmA_rank=-1
  call mpi_comm_rank(LinFmA_Univ,My_FmA_rank,ierr)
  if ((Acopy_irank == 0) .and. (My_FmA_rank /= 0)) then  ! if we are here, were are on the master mode of a sending FmA_universe
    print *, "proc", irank, "master of FmA Univ. is not the master of Acopy univ. My_FmA_rank", My_FmA_rank, "Acopy_irank",Acopy_irank
    stop
  end if   
  !!! end consistency check 
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! now if there is no self-send, then we can go ahead and destroy the array Acopy_FmA
  if (is_in_Acopy_univ == .true.) then 
   ! this is self-send, so we just need to extract data from Acopy_FmA 
   ! by this time the consolidated linearized operator must have been created and 
   ! it is ready to be populated. 
   ! essentially, all it takes is to get the data in Acopy_FmA and ADD it on the top of whatever is in BfmA.
   j=1
   do while (j<size(linprocndswld_Acopy_univs,1))
    send_Lin_FmA = linprocndswld_Acopy_univs(j) ! pick a Lin_FmA universe that sends staff
    jjj = linprocndswld_Acopy_univs(j+1)
    if (send_Lin_FmA == active_Lin_FmA) then 
     do jj=1,jjj
      Acopy_addr = linprocndswld_Acopy_addr(j+1+jj)
      BfmA_Addr = linprocndswld_BfmA_addr(j+1+jj)
      do ii = 1,nn
       BfmA(ii,BfmA_Addr) = BfmA(ii,BfmA_Addr) + Acopy_FmA(ii,Acopy_addr)
      end do 
     end do
     exit !We do not expect more than one message at a time
    end if  
    j=j+jjj+2
   end do  
  end if
  ! we will now deallocate the Acopy_FmA to free the memory
  deallocate (Acopy_FmA)
  Acopy_FmA_flag = 0
  !  
 end if 
 ! Next we see if we need to receive a copy of Acopy_FmA 
 if ((is_in_Acopy_univ == .true.) .and. (is_Acopy_master == .false.)) then  ! need to recive data 
  ! but there was no self -send 
  ! in this case, we need to call collective communication to recieve the data: 
  call mpi_bcast (LinUniv_FmA,nnds,MPI_DOUBLE_PRECISION,0,LinFmA_Univ,ierr)  ! recive Acopy_FmA from the master processor 
  if (ierr /= 0 ) then 
   print *,"PrepConsolLinColl_highperf_DGV: bcast of Acopy_FmA from ", irank, "Acopy", My_Acopy_Univ_indx," returned error. stop", ierr
   stop
  end if
  ! once the operator is received, we need to get componets out of it... 
  ! by this time the consolidated linearized operator must have been created and 
  ! it is ready to be populated. 
  ! essentially, all it takes is to get the data in Acopy_FmA and ADD it on the top of whatever is in BfmA.
  j=1
  do while (j<size(linprocndswld_Acopy_univs,1))
   send_Lin_FmA = linprocndswld_Acopy_univs(j) ! pick a Lin_FmA universe that sends staff
   jjj = linprocndswld_Acopy_univs(j+1)
   if (send_Lin_FmA == active_Lin_FmA) then 
    do jj=1,jjj
     Acopy_addr = linprocndswld_Acopy_addr(j+1+jj)
     BfmA_Addr = linprocndswld_BfmA_addr(j+1+jj)
     do ii = 1,nn
      BfmA(ii,BfmA_Addr) = BfmA(ii,BfmA_Addr) + LinUniv_FmA(ii,Acopy_addr)
     end do 
    end do
    exit !We do not expect more than one message at a time
   end if  
   j=j+jjj+2
  end do  
  ! we will now deallocate the LinUniv_FmA to fee the memory
  deallocate (LinUniv_FmA)
  LinUniv_FmA_flag = 0
 ! done with the non-self-send receive
 end if 
! done with this LinFmA universe -- moving onto next 
end do   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine PrepConsolLinColl_highperf_DGV 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionLinear_MPI_DGV
!
! This subroutine evaluates the linearized collision kernel. To call this subroutine, 
! a vector of the linearized kernel must be prepared. the vector is 2\int f_{m}(v)A(v,v_1,phi) dv
! it is calcuated in the subroutine PrepareFMA_DGV
! 
! 
! Here the arrays fmA has the linearizised operator in it
! accessed directly from the commvar
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine EvalCollisionLinear_MPI_DGV(f,fc,workload_nodes)        ! This evaluates the linear part

use commvar, only: fmA

intrinsic MATMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
integer (I4B), dimension (:), intent (in) :: workload_nodes ! these are the nodes to which fmA corresponds: in MPI 
                                 ! implementation each processor has only a portion of the linearized collision 
                                 ! operator fmA. These nodes are listed in workload_nodes array. 
                                 ! the seocnd index in fmA runs over these nodes.
                                 !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (1:size(workload_nodes,1)) :: fc_temp ! scrap array
integer (I4B) :: i ! index 

fc = 0 ! we nullify the result just in case 
fc_temp = 2*MATMUL(f,fmA) ! this is a contraction with the gradient , i.e., the derivative of the collisio operator at some maxwellian
!now we need to place the results into their places: 
do i=1,size(workload_nodes,1)
 fc(workload_nodes(i)) = fc_temp(i)
end do
! the linearized portion of collision operator is ready.
end subroutine EvalCollisionLinear_MPI_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalConsolCollisionLinear_MPI_DGV
!
! This subroutine evaluates the linearized collision kernel USING THE CONOSLIDATED COLLISION OPERATOR. 
! Other than the different collision kernel this subroutine is identical the the one above 
! 
! To call this subroutine, 
! a vector of the linearized kernel must be prepared. the vector is 2\int f_{m}(v)A(v,v_1,phi) dv
! it is calcuated in the subroutine PrepConsolLinColl_DGV 
! 
! 
! Here the arrays BfmA has the linearizised operator in it
! accessed directly from the commvar
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine EvalConsolCollisionLinear_MPI_DGV(f,fc,workload_nodes)        ! This evaluates the linear part

use commvar, only: BfmA

intrinsic MATMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
integer (I4B), dimension (:), intent (in) :: workload_nodes ! these are the nodes to which fmA corresponds: in MPI 
                                 ! implementation each processor has only a portion of the linearized collision 
                                 ! operator fmA. These nodes are listed in workload_nodes array. 
                                 ! the seocnd index in fmA runs over these nodes.
                                 !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (1:size(workload_nodes,1)) :: fc_temp ! scrap array
integer (I4B) :: i ! index 

fc = 0 ! we nullify the result just in case 
fc_temp = 2*MATMUL(f,BfmA) ! this is a contraction with the gradient , i.e., the derivative of the collisio operator at some maxwellian
!now we need to place the results into their places: 
do i=1,size(workload_nodes,1)
 fc(workload_nodes(i)) = fc_temp(i)
end do
! the linearized portion of collision operator is ready.
end subroutine EvalConsolCollisionLinear_MPI_DGV

end module collision_mod
