!
!  Collection of subroutines for 3D nodal-DG velocity discretization of the Boltzmann collision integral 
!  Will be integrated with CLAWPACK and with other codes 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program main
      
      use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
      use readwrite ! contains useful subroutines for reading/writing. in particular the reading of the problem data file
      use miscset ! contains various setup subroutines, defines useful arrays, variables etc.
      use commvar, only: k_b,k_c,d_c,d_b,N,Mv,Mu,Mw,su,sv,sw,nodes_u,nodes_v,nodes_w,f,A_capphi,&
                  Mu_list,Mv_list,Mw_list,su_list,sv_list,sw_list,e_time,dt,t_L,t_R,num_save_solution,num_eval_error,rkmts,&
                  min_sens,Trad,fm,fmA,BfmA,run_mode,procs_nodes_wld,num_lin_proc,lin_proc_nodes_wld,lin_no_update,&
                  mpi_linear_coll_world,num_linstep_wait,fmA_flag,BfmA_flag,restart_time_txt,need_to_restart,&
                  g_max_order
      
      use dgvtools_mod
      use collision_mod
      use sf02, only:f_1D3D
      use time_integr_mod
      use distributions_mod
      use miscmpiset
      
      
      
      implicit none
      
      !!! MPI !!! 
      include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
      !!! END MPI !!!  

      
     ! Variables
 
 real (DP) :: test1 ! scrap variables
 real (DP), dimension(:), allocatable :: fcol,fcol_scr, Df ! scrap variable to keep results of the evaluation of collision operator 
 integer (I4B) :: gg,loc_alloc_stat ! to keep them allocation status
 !!!!!!!!!!!
 real :: pr_time_1, pr_time_2 ! variable to calculate the processor time 
 !!!!!!!!!!! SAVING SOLUTION AND SUCH !!!!!!!!!!!!!!!
 integer (I4B) :: rec, time_step, err_count
 real (DP) :: next_time_file_record, file_record_period, next_time_error_eval, error_eval_period
 character (len=15) :: suff  ! the string to hold the suffix
 real (DP)  :: ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w ! number density, av_v
 real (DP)  :: mom3_u,mom3_v,mom3_w,mom4_u,mom4_v,mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w ! higher moments 
 real (DP), dimension (1:1000) :: time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a ! scrap arrays to keep the macroapameters in.
 real (DP), dimension (1:1000) :: mom3_u_a,mom3_v_a,mom3_w_a,mom4_u_a,mom4_v_a,mom4_w_a,mom5_u_a,mom5_v_a,mom5_w_a! scrap arrays to keep the macroapameters in.
 real (DP), dimension (1:1000) :: mom6_u_a,mom6_v_a,mom6_w_a ! scrap arrays to keep the macroapameters in.
 
 !!!!!!!!!!! MPI variables go below... 
 integer :: ierr,irank ! variables for MPI Calls
 integer (I4B) :: nmsg ! long integer to store the length of the message
 !!!!!!!!!!!
      ! Body
 !!! MPI !!!
 !!! INITIALIZE THE MPI environment
 call mpi_init(ierr)
 !!! Now the MPI routines and variables will make sense ..
 !!! MPI FORK
 call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on... 
 !!! Now we will check what processor we are on. The master processor  
 if (irank == 0) then 
     call SetUWbgkParams("parameters.dat",0) ! if the second variable is not zero than the subroutine returns printouts
 else 
     call SetUWbgkParams("parameters.dat",1) ! if the second variable =0 then the printouts are not returend 
 end if 
 !!! END MPI FORK
      ! Read the problem parameters from the data file and set the problem variables. 
      !call SetUWbgkParams("parameters.dat")
      ! set the order of used gauss-lobatto nodes in x and t
      k_b=7 
      d_b=8
      d_c=2
      k_c=2

      call SetDGblzmGnodes ! sets some basic arrays of gauss and gauss-lobatto nodes. 
      ! set the dimensions of the level zero meshes... 
      N=3
      Mv=Mv_list(1)
      Mu=Mu_list(1)
      Mw=Mw_list(1)
      su=su_list(1)
      sv=sv_list(1)
      sw=sw_list(1)
      !! Consistency check: the the max allowed order of gauss quadrature is g_max_order
      !! if the user supplied values are bigger than that, the program will issue a warning and will stop.
      if (max(su,sv,sw)>g_max_order) then 
       print *, "The supplied order of Gauss order exceeds allowed g_max_order (see in commonvar.f90). Stop"
       stop
      end if 
      !! End consistency check:
      !
      call SetDGVblzmmesh  ! sets one-dimensional meshes in X and one-dimensional grids in u,v, and w
      call Set3DCellsR_DGVblzm
      !call CellsRefineDGV((/ 14 /),2,2,2)
      !call CellsRefineDGV((/ 27 + 14 /),3,3,3)
      call SetNodesDGV
      if (irank == 0) then 
       call WriteCellsDGV
       call WriteGridsDGV
       call WriteNodesDGV
      end if 
      !!! call WriteI1DGV(0.0_DP,0.0_DP,0.0_DP)
      !!! stop
      !!! MPI FORK !!!
      call PrepMPIcollsnDGV_highperf(irank)  ! distributed work and sets up some arrays that will be used in parallel evaluation of the collision integral
      !!! END MPI FORK !!!
     
     !!! Now all the preparations have been completed and the program is ready to compute the collision integral. 
     !!! What will happen now is that the master node will broadcast the solution to the slave processors. 
     !!! The slave processes will take this solution and compute a portion of the collision integral. 
     !!! Then the master node will perform a collective operation to combine all peices into one.
     !!! THen the master node will perform the time step and update the solution. Aftert that the process is repeated 
     !!! until the desired time is reaches. 
     
     !!! To realize this algorithm we will make another fork. The master process will have the usual time integratio nroutines and all that. 
     !!! the solution boradcast will be the fisrt operation then, the master node will call the collective operation in the evaluation 
     !!! of the spatial part. 
          
     !!! The slave process will run in the loop untill a signal is received that the time integration is complete. 
     !!! the loop will consist of three operations: recived the updated f, compute the portion of the collision operator and sending 
     !!! out the collision operator
     
     !!! Allocate arrays for the solution, for the local maxwllian and for the linearizaed operator. 
     allocate (f(1:size(nodes_u,1)),fcol(1:size(nodes_u,1)),fm(1:size(nodes_u,1)),Df(1:size(nodes_u,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DGVblzm: Allocation error for variables (f), (fcol), (fm)"
     end if 
     !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !!!!!!!!!!!!!!!!!!!!!!!    
     !! INITIALIZE THE RUN_MODE to ZERO -- in the beginning the regime is strongly non-linear
     run_mode=0
     err_count=0
     lin_no_update=0 ! reset the counter for linear steps on the same derivative ...
     fmA_flag = 0; BfmA_flag = 0 ! reset the allocatable arrays status flag to "not-allocated"=0
     !!! MPI FORK !!!!
 if (irank == 0) then 
     !!! Do this on the master processor
     !!!!!!!!!!!!!!!!!!!!! TIME EVOLUTION !!!!
     rkmts = 5
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! begin restart fork
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (need_to_restart) then
       call RestartSolutionDGV(restart_time_txt)
       Read (restart_time_txt, fmt = * ,  iostat = loc_alloc_stat) e_time ! set the current time to be equal to the time of restore  
       err_count=0
     else 
      rkmts = 5
      e_time = t_L ! initial time
      f = f_1D3D(0.0_DP, nodes_u, nodes_v, nodes_w, e_time) ! set up initial data
      !!!!!!!!!!!
      ! Write the initial data on disk.
      !!!!!!!!!!!
      suff = "i"
      call WriteSol_init_DGV (suff)  
      ! compute the macroparameters of the initial data and save it on disk
      err_count=0
      !!!call MassCheckRecPlus (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w)
      call MassCheckRecHighCMoments (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w,mom3_u,mom3_v,mom3_w,mom4_u,mom4_v,&
                                       mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w)
      err_count=err_count+1
      if (err_count<1000) then
        time_a(err_count) = e_time 
        ndens_a(err_count) = ndens
        ubar_a(err_count) = ubar
        vbar_a(err_count) = vbar
        wbar_a(err_count) = wbar
        tempr_a(err_count) = tempr
        tempr_u_a(err_count) = tempr_u
        tempr_v_a(err_count) = tempr_v
        tempr_w_a(err_count) = tempr_w
        mom3_u_a(err_count) = mom3_u
        mom3_v_a(err_count) = mom3_v
        mom3_w_a(err_count) = mom3_w
        mom4_u_a(err_count) = mom4_u
        mom4_v_a(err_count) = mom4_v
        mom4_w_a(err_count) = mom4_w
        mom5_u_a(err_count) = mom5_u
        mom5_v_a(err_count) = mom5_v
        mom5_w_a(err_count) = mom5_w
        mom6_u_a(err_count) = mom6_u
        mom6_v_a(err_count) = mom6_v
        mom6_w_a(err_count) = mom6_w
      else
        print *, "the number of error evaluations exceeded the storage, err_count >= 1000"
      end if 
      !!!!!save macroparameters on hard drive
      !call WriteErrPlus_DGV (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
      !                        vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
      !                        tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count))
      call WriteErrHighCMoments_DGV (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
                              vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
                              tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count),&
                              mom3_u_a(1:err_count),mom3_v_a(1:err_count),mom3_w_a(1:err_count),mom4_u_a(1:err_count),&
                              mom4_v_a(1:err_count),mom4_w_a(1:err_count),mom5_u_a(1:err_count),mom5_v_a(1:err_count),&
                              mom5_w_a(1:err_count),mom6_u_a(1:err_count),mom6_v_a(1:err_count),mom6_w_a(1:err_count))
      !!!!! The initial data  and the initial macroparameters were saved on the hard drive    !!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      print *, "Got here...", pr_time_2 - pr_time_1 ! Report the account of time
      
      
      call cpu_time (pr_time_1) ! Start the account of time
      ! prepare MTS steps.
      time_step=0
      !!! AN MPI FORK IS HIDDEN HERE
      call PrepareMTS_SH_DGV (f, dt)
      !!! END OF THE MPI FORK   
      call cpu_time (pr_time_2) ! End the account of time
      print *, "Processor time lapsed in seconds for Prepare MTS:", pr_time_2 - pr_time_1 ! Report the account of time
      !!!!!!
      !! TESTING MPI SPEEDUP
      !! comment the next 8 lines if not testing speedup.
      !!! broadcast an  end of work code to the processors... 
      !!call mpi_bcast (-777,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! -777  means end of work
      !!if (ierr /= 0 ) then 
      !! print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
      !!  stop
      !!end if 
      !! call mpi_finalize(ierr)
      !! stop
      !!
      !! END Testing MPI SPEEDUP
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! now all the arrays are ready for the time evolution.
      ! save the solution at the initial stage:
      !write (suff, "(I4)") rec
      write (suff, "(F14.10)") e_time 
      call WriteSol_DGV (suff)  
      !!!!!evaluate macroparameters after the MTS prepare stage
      !call MassCheckRecPlus (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w)
      call MassCheckRecHighCMoments (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w,mom3_u,mom3_v,mom3_w,mom4_u,mom4_v,&
                                       mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w)
      err_count=err_count+1 
      if (err_count<1000) then ! record the macroparameters
        time_a(err_count) = e_time 
        ndens_a(err_count) = ndens
        ubar_a(err_count) = ubar
        vbar_a(err_count) = vbar
        wbar_a(err_count) = wbar
        tempr_a(err_count) = tempr
        tempr_u_a(err_count) = tempr_u
        tempr_v_a(err_count) = tempr_v
        tempr_w_a(err_count) = tempr_w
        mom3_u_a(err_count) = mom3_u
        mom3_v_a(err_count) = mom3_v
        mom3_w_a(err_count) = mom3_w
        mom4_u_a(err_count) = mom4_u
        mom4_v_a(err_count) = mom4_v
        mom4_w_a(err_count) = mom4_w
        mom5_u_a(err_count) = mom5_u
        mom5_v_a(err_count) = mom5_v
        mom5_w_a(err_count) = mom5_w
        mom6_u_a(err_count) = mom6_u
        mom6_v_a(err_count) = mom6_v
        mom6_w_a(err_count) = mom6_w
      else
        print *, "the number of error evaluations exceeded the storage, err_count >= 1000"
      end if
      !!!!!save macroparameters on hard drive
      !call WriteErrPlus_DGV (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
      !                        vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
      !                        tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count))
      call WriteErrHighCMoments_DGV (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
                              vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
                              tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count),&
                              mom3_u_a(1:err_count),mom3_v_a(1:err_count),mom3_w_a(1:err_count),mom4_u_a(1:err_count),&
                              mom4_v_a(1:err_count),mom4_w_a(1:err_count),mom5_u_a(1:err_count),mom5_v_a(1:err_count),&
                              mom5_w_a(1:err_count),mom6_u_a(1:err_count),mom6_v_a(1:err_count),mom6_w_a(1:err_count))
      !!!!! The MTS prepare stage data and the macroparameters were saved on the hard drive   !!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     end if
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
     ! END OF RESTART FORK   
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
     ! measure time:  
     call cpu_time (pr_time_1) ! Start the account of time
     ! set up error evaluation and solution recording 
     file_record_period = (t_R-t_L)/Real(num_save_solution,DP)
     next_time_file_record = e_time + file_record_period
     error_eval_period = (t_R-t_L)/Real(num_eval_error,DP)
     next_time_error_eval = e_time + error_eval_period
     ! THE TIME LOOP
     do while (e_time <= t_R+0.00000000001_dp)
      time_step=time_step+1
      !!! AN MPI FORK IS HIDDEN HERE
      call TimeIntegratorMTS_SH_MPI_DGV(f,dt)
      !!! END OF AN MPI FORK
      ! check if it is time to evaluate error or any other quantity (like mass, bulk vel and tempr)
      if (e_time >= next_time_error_eval) then
       next_time_error_eval = e_time + error_eval_period
       ! evaluate mass momenum and tempretaure of the solution 
       ! end eevaluate mass momenum and tempretaure of the solution
       !call MassCheckRecPlus (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w)
       call MassCheckRecHighCMoments (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w,mom3_u,mom3_v,mom3_w,mom4_u,mom4_v,&
                                       mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w)
       err_count=err_count+1 
       if (err_count<1000) then
        time_a(err_count) = e_time 
        ndens_a(err_count) = ndens
        ubar_a(err_count) = ubar
        vbar_a(err_count) = vbar
        wbar_a(err_count) = wbar
        tempr_a(err_count) = tempr
        tempr_u_a(err_count) = tempr_u
        tempr_v_a(err_count) = tempr_v
        tempr_w_a(err_count) = tempr_w
        mom3_u_a(err_count) = mom3_u
        mom3_v_a(err_count) = mom3_v
        mom3_w_a(err_count) = mom3_w
        mom4_u_a(err_count) = mom4_u
        mom4_v_a(err_count) = mom4_v
        mom4_w_a(err_count) = mom4_w
        mom5_u_a(err_count) = mom5_u
        mom5_v_a(err_count) = mom5_v
        mom5_w_a(err_count) = mom5_w
        mom6_u_a(err_count) = mom6_u
        mom6_v_a(err_count) = mom6_v
        mom6_w_a(err_count) = mom6_w
       else
        print *, "the number of error evaluations exceeded the storage, err_count >= 1000"
       end if  
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! DIAGNOSTIC call to record a derivatives of macroparameters and the relaxation times of macroparameters. 
       ! The function will use the updated frhs1 directly from commvar and the old f can be found in f1 directly from commvar; 
       ! The time corresponds to the time step of f1. 
       !!!!
       call RecordMomDerivRelaxTimes_DGV(e_time-dt)
       !!!!  
       ! END DIAGNOSTIC
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      end if 
      ! Check if it is time to save the solution 
      if (e_time >= next_time_file_record) then 
        next_time_file_record = next_time_file_record + file_record_period
      ! write the solution of disk
       !rec=rec+1
       !write (suff, "(I4)") rec
       write (suff, "(F14.10)") e_time 
       call WriteSol_DGV (suff)
       !call WriteErrPlus_DGV (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
       !                       vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
       !                       tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count))
       call WriteErrHighCMoments_DGV (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
                              vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
                              tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count),&
                              mom3_u_a(1:err_count),mom3_v_a(1:err_count),mom3_w_a(1:err_count),mom4_u_a(1:err_count),&
                              mom4_v_a(1:err_count),mom4_w_a(1:err_count),mom5_u_a(1:err_count),mom5_v_a(1:err_count),&
                              mom5_w_a(1:err_count),mom6_u_a(1:err_count),mom6_v_a(1:err_count),mom6_w_a(1:err_count))
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! DIAGNOSTIC call to write on the hard drive the derivatives of macroparameters and the relaxation times of macroparameters. 
       !!!!
       call WriteMomDerivRelaxTimes_DGV
       !!!!  
       ! END DIAGNOSTIC
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                       
      end if
      print *, "time step", time_step     
     end do
     call cpu_time (pr_time_2) ! End the account of time
     print *, "Processor time lapsed in seconds for main run:", pr_time_2 - pr_time_1 ! Report the account of time
     !!! broadcast an  end of work code to the processors... 
     call mpi_bcast (-777,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)  ! -777  means end of work
     if (ierr /= 0 ) then 
       print *,"CheckSolutionMode_MPI_DGV: MPI boradcast of run_mode from proc 0 returned error", ierr
       stop
     end if 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else  ! MPI FORK
     ! do this on the slave processor... 
     ! need this scrape variables:
     allocate (fcol_scr(1:size(nodes_u,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DGVblzm: Allocation error for variables (fcol_scr)"
     end if 
     !
     ! Next the slave mode goes in the wait mode for the directives from the master node and exists there untill recieves the exist code -777
     do while (run_mode /= -777) 
       ! First we wait for the transmission of the run_mode
        call mpi_bcast (run_mode,1,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
        if (ierr /= 0 ) then 
          print *,"DGVblzmPL: slave processor", irank, "MPI boradcast of run_mode from proc 0 returned error", ierr
        stop
        end if
        ! if the run_mode /= -777 we will need the solution 
        if (run_mode /= -777) then    
         ! we wait for the transmission of the solution 
         nmsg=size(f,1)
         call mpi_bcast (f,nmsg,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         if (ierr /= 0 ) then 
          print *,"DGVblzmPL: slave processor", irank, "MPI boradcast of f from proc 0 returned error", ierr
         stop
         end if   
        else 
         print *, "run_mode = -777. stop next."
        end if 
        ! nullify the scrap variables:
        fcol=0
        fcol_scr=0
        !
        !
        select case (run_mode)
         case (0) ! mode zero means we are in non-linear regime and collision integral needs to be computed. 
          ! this little block releases memory used by the linearized operators
          if (BfmA_flag==1) then 
            deallocate(BfmA)  ! release the memory used by the consolidated linear operator. 
            BfmA_flag=0 ! set the flag into the "not allocated= 0" 
          end if 
          if (fmA_flag==1) then 
            deallocate(fmA)  ! release the memory used by the consolidated linear operator. 
            fmA_flag=0 ! set the flag into the "not allocated= 0" 
          end if 
          ! end of the memory release block...     
          lin_no_update = 0 ! reset the counter for linerizaed steps with no update of the derivative
          ! this is a fully nonlinear case, we need to evaluate the portion of the collision integral
          call EvalCollisionPeriodicAPlus_MPI_DGV(f,fcol) !MPI FORK HERE
          ! see for the consolidated statement on mpi_reduce. Done the case of full nonlinear collision operator  
         case (1) !
          ! this is the case when decomposition is used but no re-caclulation 
          ! of the linearized operator is necessary 
          ! we do not need to recalculate the Maxwellian (fm) -- it must be good to go somewhere in the memory:
          lin_no_update = 0 ! reset the counter for linerizaed steps with no update of the derivative 
          Df=f-fm ! update the perturbatin of the Maxwellian 
          call EvalCollisionLinear_MPI_DGV(Df,fcol_scr,procs_nodes_wld)  ! This evaluates the linear part
          call EvalCollisionPeriodicAPlus_MPI_DGV(Df,fcol) ! this evaluates the non-linear part 
          fcol = fcol+fcol_scr ! this adds the linear and nonlinear parts
          ! see below for adding a dimensionless coefficient and for sending the solution out
         case (11)
          ! this is the case when decomposition is used but re-caclulation 
          ! of the linearized operator is necessary
          ! this little block releases memory used by the linearized operators
          if (BfmA_flag==1) then 
            deallocate(BfmA)  ! release the memory used by the consolidated linear operator. 
            BfmA_flag=0 ! set the flag into the "not allocated= 0" 
          end if 
          ! end of the memory release block...     
          lin_no_update = 0 ! reset the counter for linerizaed steps with no update of the derivative
          call MassCheckRec (f,ndens,ubar,vbar,wbar,tempr) ! First, we evaluate the macroparamters of the solution.
          fm = maxwelveldist(tempr,ubar,vbar,wbar,ndens,nodes_u, nodes_v, nodes_w) ! now we populate the maxwellian with the same macroparamters.
		  call PrepareFMA_MPI_DGV(irank) ! this will prepare the linearized operator  fmA = 2\int f_{m}(v)A(v,v_1,phi) dv
          !! Now we call the evaluation of the portion of the collision operator
          Df=f-fm ! update the perturbatin of the Maxwellian 
          !!! Uncomment next line if need to remove using derivative of collision operator because of the storage limitaions.
          !!! call EvalCollisionPeriodicMixedTermsA_DGV(Df,fm,fcol_scr)
          call EvalCollisionLinear_MPI_DGV(Df,fcol_scr,procs_nodes_wld)  ! This evaluates the linear part using the derivative
          !!!! 
          call EvalCollisionPeriodicAPlus_MPI_DGV(Df,fcol) ! this evaluates the non-linear part 
          fcol = fcol+fcol_scr ! this adds the linear and nonlinear parts
          ! see below for adding a dimensionless coefficient and for sending the solution out
          !!!!!!!! DEBUG !!!!!!!!!!!!!!!!!!!!!     
         case (2)
          ! this is the case when the solution is linearized, however, old linear operator can be used.
          lin_no_update = lin_no_update + 1      
          ! if the number of linear steps is within the maximum allowed, evaluate the linearized solution and send it out
          if (lin_no_update <= num_linstep_wait) then 
            ! evaluate the linearized solution
            Df=f-fm ! update the perturbatin of the Maxwellian 
            call EvalCollisionLinear_MPI_DGV(Df,fcol,procs_nodes_wld)  ! This evaluates the linear part
          end if 
          if (lin_no_update == num_linstep_wait+1) then 
            ! else need to send out the pieces for the derivative for consolidation and if this is one of the 
            ! processors that performs the consolidated evaluation of linear collision operator, then this will prepare fmA as well
            ! all is accomplished by the call of the subroutine:
            call PrepConsolLinColl_highperf_DGV(irank)  !! call the subroutine that sets up the consolidated linearized operator... 
          end if  
          if ((irank < num_lin_proc) .and. (lin_no_update > num_linstep_wait)) then 
            ! evaluate the linearized solution
            Df=f-fm ! update the perturbatin of the Maxwellian 
            call EvalConsolCollisionLinear_MPI_DGV(Df,fcol,lin_proc_nodes_wld)  ! This evaluates the linear part
          end if 
         case (21) 
           ! this is the case when the solution is linearized and the new linearized operator needs to be computed 
           ! this little block releases memory used by the linearized operators
           if (BfmA_flag==1) then 
            deallocate(BfmA)  ! release the memory used by the consolidated linear operator. 
            BfmA_flag=0 ! set the flag into the "not allocated= 0" 
           end if 
           ! end of the memory release block...     
           lin_no_update = 0 ! reset the counter for linerizaed steps with no update of the derivative
           call MassCheckRec (f,ndens,ubar,vbar,wbar,tempr) ! First, we evaluate the macroparamters of the solution.
           fm = maxwelveldist(tempr,ubar,vbar,wbar,ndens,nodes_u, nodes_v, nodes_w) ! now we populate the maxwellian with the same macroparamters.
		   call PrepareFMA_MPI_DGV(irank) ! this will prepare the linearized operator  fmA = 2\int f_{m}(v)A(v,v_1,phi) dv
           !! Now we call the evaluation of the portion of the collision operator
           Df=f-fm ! update the perturbatin of the Maxwellian 
           call EvalCollisionLinear_MPI_DGV(Df,fcol,procs_nodes_wld)  ! This evaluates the linear part
         case (-777)
           print *,"DGVblzmPL: slave processor", irank, "end of work signal recieved. Exit now" 
         case default
           print *,"DGVblzmPL: slave processor", irank, "Exit Error: unknown case of run_mode received, run_mode=", run_mode 
           stop
        end select
        !!! These lines are consolidated from the sace select to minimize code repetition:
        !!! in all cases below the solution needs to be sent to the master node via mpi_reduce
        if ((run_mode == 0) .or. (run_mode == 1)  .or. (run_mode == 11)  .or. (run_mode == 21) &
                   .or. ((run_mode == 2) .and. (lin_no_update <= num_linstep_wait))) then 
          ! Now we send back the portion of the collision integral
          nmsg=size(f,1)
          call mpi_reduce (fcol,f,nmsg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr) ! f for the receive buffer should be ignored
          if (ierr /= 0 ) then 
           print *,"DGVblzmPL: slave processor", irank, "MPI_REDUCE to proc 0 returned error", ierr
           stop
          end if
        end if
        if (((run_mode == 2) .and. (lin_no_update > num_linstep_wait)) .and. (irank < num_lin_proc)) then 
          ! Now we send back the portion of the collision integral
           nmsg=size(f,1)
          call mpi_reduce (fcol,f,nmsg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_LINEAR_COLL_WORLD,ierr) ! f for the receive buffer should be ignored
          if (ierr /= 0 ) then 
           print *,"DGVblzmPL: slave processor", irank, "MPI_REDUCE to proc 0 returned error", ierr
           stop
          end if
        end if 
        !! We sent the solution out
     end do    
      print *, "slave process", irank, "exiting the main loop. Bye."
     
 end if  
 !!! END MPI FORK
 !!!!!!!!!!!!!!!!!!!!!

call mpi_finalize(ierr)


end program main
