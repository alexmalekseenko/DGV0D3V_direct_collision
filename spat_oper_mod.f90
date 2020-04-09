!
!  spat_oper_mod.f90
!
! This module contains routines that are involved in the evaluation of the spatial operator  
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module spat_oper_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none



contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SpatialOperSH_DGV (f) 
!
! This operator evaluates the spatial operator of the Boltzmann equation 
! in the spatially homogeneous case
!
! COMMENT: STILL NEED TO ADD COLLISION FREQUENCY! 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SpatialOperSH_DGV (f,fcol,dt)

use commvar, only: run_mode,mol_diam,L_inf,N_inf

use collision_mod

real (DP), dimension (:), intent (in) :: f ! the solution at the current  
real (DP), dimension (:), intent (out) :: fcol ! value of the right side  
real (DP), intent (in) :: dt ! The time step 
!!!!!!!!!!!
real (DP), dimension (size(fcol,1)) :: fcol_scr,Df ! scatch variable to keep the right side and the perturbation part
real (DP) :: coef_temp ! Scrap variable 

! WARNING: Make sure that Nodes_Ashift, Nodes_ccan and other supplementary arrays are set before calling the 
! the evaluation of the collision operator.

call CheckSolutionMode_DGV(f,Df) ! This will check the local distribution function against a maxwellian with the same 5 macroparameters.

select case (run_mode) ! run_mode is set outside, when solution is evaluated for closeness to a maxwellian...  

case (0) ! run_mode=0 means we are very far from a Maxwellian. In this case we just call the collision operator
  !!! Two choises for the call of collision operator: Uncomment only one of them!  These procedures do the same, but slightly different in implementation
  !!! call EvalCollisionPeriodicA_DGV(f,fcol) ! This one uses intermediate arrays and is slower
  call EvalCollisionPeriodicAPlus_DGV(f,fcol) ! This one is a little faster that the one above... 
  !!!!
case (1) ! this is the non-linear perturbation mode: we need both the linear part and the non-linear 
  call EvalCollisionLinear(Df,fcol_scr)        ! This evaluates the linear part
  call EvalCollisionPeriodicAPlus_DGV(Df,fcol) ! this evaluates the non-linear part 
  fcol = fcol+fcol_scr
  ! 
case(2) ! this is the linear mode, we pretty much neglect the quadratic part..
  call EvalCollisionLinear(Df,fcol)        ! This evaluates the linear part
case default 
    print *, "SpatialOperSH_DGV: cannot process value of run_mode", run_mode
    stop
end select
!! Now we are addiing the dimensionless coefficient mol_diam^2*N_inf/(L_inf)^2
coef_temp = (mol_diam/L_inf)**2*N_inf*dt
fcol = fcol*coef_temp
!! 
end subroutine SpatialOperSH_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SpatialOperSH_MPI_DGV (f) 
!
! This is an MPI enhansed version of the previous subroutine
!
! This operator evaluates the spatial operator of the Boltzmann equation 
! in the spatially homogeneous case
!
! COMMENT: STILL NEED TO ADD COLLISION FREQUENCY! 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SpatialOperSH_MPI_DGV (f,fcol,dt)

use commvar, only: run_mode, mol_diam, L_inf, N_inf,lin_no_update,num_linstep_wait,&
                   lin_proc_nodes_wld,MPI_LINEAR_COLL_WORLD,nodes_u,nodes_v,nodes_w,&
                   BfmA,BfmA_flag   

use collision_mod
use distributions_mod
use dgvtools_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!! 

real (DP), dimension (:), intent (in) :: f ! the solution at the current  
real (DP), dimension (:), intent (out) :: fcol ! value of the right side  
real (DP), intent (in) :: dt ! The time step  
!!!!!!!!!!!
real (DP) :: coef_temp ! Scrap variable 
real (DP), dimension (1:size(fcol,1)) :: fcol_scr ! scrap variable tenporary value of the right side  
real (DP), dimension (1:size(f,1)) :: Df ! scrap variable to keep the difference f-f_M (perturbation from maxwellian)
real (DP), dimension (1:size(f,1)) :: fm ! scrap variable to keep the f_M (the local maxwellian)
real (DP) :: ndens,ubar,vbar,wbar,tempr ! scrap variables to keep macroparameters

!! MPI VARIABLES:
integer :: ierr,irank ! status variable 
integer (I4B) :: nmsg ! long integer to store the length of the message

! WARNING: Make sure that Nodes_Ashift, Nodes_ccan and other supplementary arrays are set before calling the 
! the evaluation of the collision operator.
! has an MPI broadcast in it  
call CheckSolutionMode_MPI_DGV(f,Df) ! This will check the local distribution function against a maxwellian with the same 5 macroparameters.
! end have MPI broadcast 
! Next , we broadcast the current solution to all nodes
nmsg=size(f,1)
call mpi_bcast (f,nmsg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if (ierr /= 0 ) then 
   print *,"SpatialOperSH_MPI_DGV: MPI boradcast of f from proc 0 returned error", ierr
   stop
  end if   
! after this broadcast the slave modes will evaluate the solution. We just need to assemble it here:
fcol = 0.0_DP     ! nullify the results first, just in case being paranoid
fcol_scr = 0.0_DP ! nullify the results first, just in case being paranoid
!!! in all cases below the solution needs to be reduced to the master node via mpi_reduce

if ((run_mode == 0) .or. (run_mode == 1) .or. ((run_mode == 2) .and. (lin_no_update <= num_linstep_wait ))) then 
   !!!!!
   ! REMARK: lin_no_update on the master node is modified in CheckSolutionMode_MPI_DGV
   !!!!!
   ! A little block to release the memory used by the consolidated linerar operator
   if (BfmA_flag==1) then 
    deallocate(BfmA)  ! release the memory used by the consolidated linear operator. 
    BfmA_flag=0 ! set the flag into the "not allocated= 0" 
   end if 
   ! Now we consolidate the solution on the master node:
   nmsg = size(fcol,1) ! long integer to store the length of the message
   call mpi_reduce (fcol_scr,fcol,nmsg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr) !
   if (ierr /= 0 ) then 
     print *,"SpatialOperSH_MPI_DGV:: master node MPI_REDUCE returned error", ierr
     stop
   end if
end if
!! Next we describe what happens when the evaluation of the linearized collision operator is 
!! consolidated on a group of processors. These processors are clustered in a special world: 
!!  MPI_LINEAR_COLL_WORLD -- these are the processors that consolidate the evaluation of the linearized 
!! collision operator. Processors with numbers 0 <= num_lin_proc - 1 go into this world
!! Essentially, we just need to call the linearied evaluation on each of the processors and then assemble 
!! the result on the master node by doing a reduce operation. 
!! Note that the master node is in the MPI_LINEAR_COLL_WORLD world so, we will do the next do next lines 
!! whenever the lin_no_update => num_linstep_wait
if ((run_mode == 2) .and. (lin_no_update == num_linstep_wait+1)) then 
 call PrepConsolLinColl_highperf_DGV(0) !! call the subroutine that sets up the consolidated linearized operator... 
end if 
!! next we call the evaluation of the collision part
if ((run_mode == 2) .and. (lin_no_update > num_linstep_wait)) then
 ! first we evaluate the part that belongs to the master node: 
 call MassCheckRec (f,ndens,ubar,vbar,wbar,tempr) ! First, we evaluate the macroparamters of the solution.
 fm = maxwelveldist(tempr,ubar,vbar,wbar,ndens,nodes_u, nodes_v, nodes_w) ! now we populate the maxwellian with the same macroparamters.
 !! Now we call the evaluation of the portion of the collision operator
 Df=f-fm ! update the perturbatin of the Maxwellian 
 call EvalConsolCollisionLinear_MPI_DGV(Df,fcol_scr,lin_proc_nodes_wld)  ! This evaluates the linear part 
 !! next we compile our part with other parts wia MPI-REDUCE
 nmsg=size(fcol,1) ! long integer to store the length of the message
 call mpi_reduce (fcol_scr,fcol,nmsg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_LINEAR_COLL_WORLD,ierr) !
    if (ierr /= 0 ) then 
      call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on... 
      print *,"SpatialOperSH_MPI_DGV:: slave processor", irank, "MPI_REDUCE to proc 0 returned error", ierr
      stop
    end if
end if 
!! Done with the consolidation of the linearied operator
  
!! Now we are addiing the dimensionless coefficient mol_diam^2*N_inf/(L_inf)^2
coef_temp = (mol_diam/L_inf)**2*N_inf*dt
fcol = fcol*coef_temp
!! 
end subroutine SpatialOperSH_MPI_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine CheckSolutionMode_DGV(f,Df) 
! 
! This subroutine evaluates macroparameters of the solution. Then it subtracts 
! from the solution a Maxwellian with the computed macroparameters. This difference is placed into Df
! then a norm of Df is assessed. If |Df| is large, (rum_mode=0) then we are in a strongly non-equilibrium 
! regime and noting is done. 
! if the |Df| is moderate, than we switch to mode 1 (run_mode=1) then linear and quadratic contributions are computed separately
! if |Df| is small, we switch to run_mode=2 and only the linear part of collision opertator is computed.   
! 
! This subroutine also checks the solution and calles for the update of the linearized operator fmA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckSolutionMode_DGV(f,Df)

use commvar, only: run_mode, LocDens, LocUbar, LocVbar, LocWbar, LocTempr, fm, nodes_u, nodes_v, nodes_w, &
                   decomp_lev, linear_lev,nodes_gwts

use distributions_mod
use dgvtools_mod
use collision_mod
use readwrite

Intrinsic SUM, ABS

real (DP), dimension (:), intent(in) :: f ! the solution on the current state
real (DP), dimension (:), intent(out)::Df ! the difference between the solutio and the local Maxwellian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real (DP), dimension (size(fm,1)) :: fMaxwellNew  ! a scratch storage for the maxwellian 
real (DP) :: L1_err, L1_err_fm ! a sctratch variables for the errors

 
call MassCheckRec (f,LocDens,LocUbar,LocVbar,LocWbar,LocTempr) ! First, we evaluate the macroparamters of the solution.
fMaxwellNew = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_u, nodes_v, nodes_w) ! now we populate the maxwellian with the same macroparamters.
L1_err = SUM(ABS(f-fMaxwellNew)*nodes_gwts)/LocDens! evaluate the relative L1-norm of the difference (L1 norm of the maxwellian should be density)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! will try to track the deviation from the Maxwellian 
! For that purpose we will create a global allocatable array that will store the times and the values of L1_err 
! It will periodically damp it to disk... all on its own... 
!!!!!!!!!!!!
call KeepTrackL1_err_DGV(L1_err)
!!!!!!!!!!!!
if (L1_err < decomp_lev) then 
   select case (run_mode) 
   case(0) ! if we were in a strongly non-linear regime on the previous time step, then do this: 
      fm=fMaxwellNew
      call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_{m}(v)A(v,v_1,phi) dv
      if (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2                      ! linear regime
        print *, "Switching to linearized regime from strongly non-linear"
      else
        run_mode=1                      ! decomposition regime
        print *, "Switching to decomposition regime from strongly non-linear"
      end if   
   case(1) ! if already were in the decomposition regime
      L1_err_fm=SUM(ABS(fm-fMaxwellNew)*nodes_gwts)/LocDens   ! next we check is the linearized operator needs to be updated
      if (L1_err_fm > linear_lev/2.0_DP ) then ! we use the linear_level as a measure of two perturbations being close
        fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated: 
        call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_{m}(v)A(v,v_1,phi) dv
      end if
      if (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode = 2                    ! linear regime
        print *, "Switching to linearized regime from decomposition regime"
      end if   
   case (2) ! if were in the linear regime.. 
      L1_err_fm=SUM(ABS(fm-fMaxwellNew)*nodes_gwts)/LocDens   ! next we check is the linearized operator needs to be updated
      if (L1_err_fm > linear_lev/2.0_DP) then ! we use the linear_level as a measure of two perturbations being close
        fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated: 
        call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_{m}(v)A(v,v_1,phi) dv
      end if
      if (L1_err > linear_lev) then     ! Detemine the new regime:
        run_mode = 1                    ! linear regime
        print *, "Switching to decomposition regime from linear regime"
      end if
   case default
      print *, "CheckSolutionMode_DGV: Error, can no tprocess the value of run_mode=", run_mode    
      stop
   end select       
   Df=f-fm                           ! evaluate the perturbatio nfrom maxwellian
else    
  if (run_mode==1) then 
  print *, "Switching from decomposition mode to strongly non-linear mode."    
  end if 
  if (run_mode==2) then 
  print *, "Switching from linear mode to strongly non-linear mode."    
  end if 
  run_mode=0
end if 
end subroutine CheckSolutionMode_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine CheckSolutionMode_MPI_DGV(f,Df)
!
! This is an MPI version of the above subroutine. Essentially, only three subroutine calls had changed. All the rest is the same.  
! 
! This subroutine evaluates macroparameters of the solution. Then it subtracts 
! from the solution a Maxwellian with the computed macroparameters. This difference is placed into Df
! then a norm of Df is assessed. If |Df| is large, (rum_mode=0) then we are in a strongly non-equilibrium 
! regime and noting is done. 
! if the |Df| is moderate, than we switch to mode 1 (run_mode=1) then linear and quadratic contributions are computed separately
! if |Df| is small, we switch to run_mode=2 and only the linear part of collision opertator is computed.   
! 
! This subroutine also checks the solution and calles for the update of the linearized operator fmA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckSolutionMode_MPI_DGV(f,Df)

use commvar, only: run_mode, LocDens, LocUbar, LocVbar, LocWbar, LocTempr, fm, nodes_u, nodes_v, nodes_w, &
                   decomp_lev, linear_lev,nodes_gwts,lin_no_update,num_linstep_wait

use distributions_mod
use dgvtools_mod
use collision_mod
use readwrite

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!! 

Intrinsic SUM, ABS

real (DP), dimension (:), intent(in) :: f ! the solution on the current state
real (DP), dimension (:), intent(out)::Df ! the difference between the solutio and the local Maxwellian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real (DP), dimension (size(fm,1)) :: fMaxwellNew  ! a scratch storage for the maxwellian 
real (DP) :: L1_err, L1_err_fm ! a sctratch variables for the errors
 
!! MPI VARIABLES:
integer :: ierr,irank ! status variable 

 
call MassCheckRec (f,LocDens,LocUbar,LocVbar,LocWbar,LocTempr) ! First, we evaluate the macroparamters of the solution.
fMaxwellNew = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_u, nodes_v, nodes_w) ! now we populate the maxwellian with the same macroparamters.
L1_err = SUM(ABS(f-fMaxwellNew)*nodes_gwts)/LocDens! evaluate the relative L1-norm of the difference (L1 norm of the maxwellian should be density)
!!!! Begin debug 
  print *,"L1_err", L1_err
!!!!! end debug
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! will try to track the deviation from the Maxwellian 
! For that purpose we will create a global allocatable array that will store the times and the values of L1_err 
! It will periodically damp it to disk... all on its own... 
!!!!!!!!!!!! Decomposition err treshold must be smaller than the linearization error threshold. Here is a quick check to ensure it
if (decomp_lev<=linear_lev) then 
 print *,"CheckSolutionMode_MPI_DGV: Error: Decomposition treshold is lower than linearization treshhold. Stop."
 stop
end if 
!!!!!!!!!!!! 
call KeepTrackL1_err_DGV(L1_err)
!!!!!!!!!!!!
if (L1_err < decomp_lev) then 
   select case (run_mode) 
   case(0) ! if we were in a strongly non-linear regime on the previous time step, then do this: 
      if (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2                      ! linear regime
        print *, "Switching to linearized regime from strongly non-linear"
        ! Then we send the transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
        call mpi_bcast (21,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 21 means linearized regime and new linearized operator is needed
        if (ierr /= 0 ) then 
         print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
         stop
        end if  
      else
        run_mode=1                      ! decomposition regime
        print *, "Switching to decomposition regime from strongly non-linear"
        ! Then we send the transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
        call mpi_bcast (11,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 11 means decomposition regime and new linearized operator is needed
        if (ierr /= 0 ) then 
         print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
         stop
        end if
      end if  
      fm=fMaxwellNew
      ! we do not broadcast fm therefore fm needs to be re-calulated on the slave processor
   case(1) ! if already were in the decomposition regime
      L1_err_fm=SUM(ABS(fm-fMaxwellNew)*nodes_gwts)/LocDens   ! next we check is the linearized operator needs to be updated
      if (L1_err_fm > linear_lev/2.0_DP ) then ! we use the linear_level as a measure of two perturbations being close
        fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated: 
        if (L1_err < linear_lev) then     ! Detemine the new regime:
         run_mode = 2                    ! linear regime
         print *, "Switching to the linearized regime from the decomposition regime."
         ! We send a transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
         call mpi_bcast (21,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 21 means linearized regime and new linearized operator is needed
         if (ierr /= 0 ) then 
          print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
          stop
         end if
        else 
         ! We send a transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
         call mpi_bcast (11,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 11 means decomposition regime and new linearized operator is needed
         if (ierr /= 0 ) then 
          print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
          stop
         end if
        end if   
      else 
        if (L1_err < linear_lev) then     ! Detemine the new regime:
         run_mode = 2                    ! linear regime
         print *, "Switching to the linearized regime from the decomposition regime."
         ! We send a transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
         call mpi_bcast (2,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 2 means linearized regime and new linearized operator is NOT needed
         if (ierr /= 0 ) then 
          print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
          stop
         end if
         !! Now if CheckSolutionMode_MPI_DGV broadcaseted run_mode = 2 --- that is the linearized collision operator needs to be 
         !! evaluated, and the linearized operator need not to be updated, the software starts to count how many consequtive steps 
         !! the run_mode=2 is boradcasted. Normally, the solution is assembled in the spatial operator subroutine from the slave 
         !! processes. However, if the solution is reaching Maxwellian and it is getting closer to a particular Maxwellian during 
         !! many time steps, then solution stays in the linearized regime and also no recalculating of the linearized operator is 
         !! necessary. Because contracting the linearized collision operator with the solution is a lot faster 
         !! than the evaluating the linearized operator and evaluating the full collision operator, it is desired to consolidate 
         !! the evaluatin of the linearized operator on fewer porocessors. 
         !! for that purpose, pieces of the linearized collision integral will be send to a dedicated group of processors as identified by local 
         !! arrays receive_lin_procs,send_lin_procs. The consolidation operation will begin after run_mode=2 has been transmitted 
         !! num_linstep_wait times. -- This will correspond to num_linstep_wait consequive evaluations of linearized collision term without 
         !! updating the linearized operator. Everytime solution switches to any other regime, the counting is interrupted. 
         lin_no_update=lin_no_update+1 ! this controls the lin_no_update on master node!
        else 
         ! We send a transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
         call mpi_bcast (1,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 1 means decomposition regime and new linearized operator is NOT needed
         if (ierr /= 0 ) then 
          print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
          stop
         end if
        end if   
      end if
   case (2) ! if were in the linear regime.. 
      L1_err_fm=SUM(ABS(fm-fMaxwellNew)*nodes_gwts)/LocDens   ! next we check is the linearized operator needs to be updated
      if (L1_err_fm > linear_lev/2.0_DP ) then ! we use the linear_level as a measure of two perturbations being close
        ! nullify the linearization with no update counter
        lin_no_update=0
        !
        fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated: 
        if (L1_err >= linear_lev) then    ! Detemine the new regime:
         run_mode = 1                    ! decomposition regime
         print *, "Switching to the decomposion regime from the linearized regime."
         ! We send a transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
         call mpi_bcast (11,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 11 means decomposition regime and new linearized operator is needed
         if (ierr /= 0 ) then 
          print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
          stop
         end if
        else 
         ! We send a transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
         print *, "mode 21: need to update the linear operator"
         call mpi_bcast (21,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 21 means linearized regime and new linearized operator is needed
         if (ierr /= 0 ) then 
          print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
          stop
         end if
        end if   
      else 
        if (L1_err >= linear_lev) then     ! Detemine the new regime:
         run_mode = 1                    ! decomposition regime
         print *, "Switching to the decomposition from the linearized regime."
         ! now we need to check is the consolidated operator was used in the linearized regime. If it was, then 
         ! all local fmA arrays were deallocated. In this case, we need to allocate them again. Therefore we 
         ! transmit run_mode=11 in this case. If the consolidation was not used, then fmA exist and can be used
         if (lin_no_update > num_linstep_wait) then ! this is if the consolidation has been invoked... 
          ! We send a transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
          call mpi_bcast (11,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 1 means decomposition regime and new linearized operator is NOT needed
          if (ierr /= 0 ) then 
           print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
           stop
          end if
         else  ! this is if the colsolidation has not been envoked
          ! We send a transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
          call mpi_bcast (1,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 1 means decomposition regime and new linearized operator is NOT needed
          if (ierr /= 0 ) then 
           print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
           stop
          end if
         end if 
         !! nullify the linearization with no update counter
         lin_no_update=0
         !!
        else 
         ! We send a transmission about the new run_mode and that we need to caclulate the linearized operator to slave processes
         call mpi_bcast (2,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! 2 means linearized regime and new linearized operator is NOT needed
         if (ierr /= 0 ) then 
          print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
          stop
         end if
         !! Now if CheckSolutionMode_MPI_DGV broadcaseted run_mode = 2 --- that is the linearized collision operator needs to be 
         !! evaluated, and the linearized operator need not to be updated, the software starts to count how many consequtive steps 
         !! the run_mode=2 is boradcasted. Normally, the solution is assembled in the spatial operator subroutine from the slave 
         !! processes. However, if the solution is reaching Maxwellian and it is getting closer to a particular Maxwellian during 
         !! many time steps, then solution stays in the linearized regime and also no recalculating of the linearized operator is 
         !! necessary. Because contracting the linearized collision operator with the solution is a lot faster 
         !! than the evaluating the linearized operator and evaluating the full collision operator, it is desired to consolidate 
         !! the evaluatin of the linearized operator on fewer porocessors. 
         !! for that purpose, pieces of the linearized collision integral will be send to a dedicated group of processors as identified by local 
         !! arrays receive_lin_procs,send_lin_procs. The consolidation operation will begin after run_mode=2 has been transmitted 
         !! num_linstep_wait times. -- This will correspond to num_linstep_wait consequive evaluations of linearized collision term without 
         !! updating the linearized operator. Everytime solution switches to any other regime, the counting is interrupted. 
         lin_no_update=lin_no_update+1
        end if   
      end if
      
   case default
      print *, "CheckSolutionMode_MPI_DGV: Error, can not process the value of run_mode=", run_mode    
      stop
   end select       
   Df=f-fm                           ! evaluate the perturbatio nfrom maxwellian
else    
  if (run_mode==1) then 
   print *, "Switching from decomposition mode to strongly non-linear mode."    
  end if 
  if (run_mode==2) then 
   print *, "Switching from linear mode to strongly non-linear mode."    
   !!! nullify the linearization with no update counter
   lin_no_update=0
   !!!
  end if 
  run_mode=0
  ! Then we send the transmission about the new run_mode to slave processes
  call mpi_bcast (run_mode,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr) 
  if (ierr /= 0 ) then 
   print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
   stop
  end if 
end if 
end subroutine CheckSolutionMode_MPI_DGV


end module spat_oper_mod 