!
!  readwrite.f90
!
! Alex 09/12/11
!
! This modle contains useful procedures for reading and writing from the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!1

module readwrite
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

contains

subroutine SetUWbgkParams(pfname,slt)
use commvar
intrinsic Real,Index, Len, Scan, Trim

character (len=*), intent (in) :: pfname ! Name of the file where the parameters are stored 
integer, intent (in) :: slt ! parameter detemining if the prinout is generated. If = 0 -- no printout generated. 
character (len=132) :: line              ! string to keep one record 
character (len=10) :: fmtchar            ! string to keep format 
character (len=50) :: line_head          ! string to keep the line header 
integer (I4B) :: code_file, code_line          ! get the return code from READ
integer (I4B) :: m_count                            ! dump counter 
integer (I4B) :: pos                 ! to store position within the string 
integer (I4B), dimension (20) :: i_bulk  ! to store temporarily integers that has been read from the parameter file 
real (DP), dimension (20) :: r_bulk      ! to store temporarily reals that has been  read from the parameter file 
integer (I4B) :: loc_alloc_stat ! some dump varaible
!
open (15, file=pfname, position ="REWIND", action="READ") ! open the file for reading
code_file=0
do while (code_file == 0)
 read (15, "(A)", iostat=code_file)  line                               ! read one line
 pos = Scan(line, "=")
  if ((pos /= 0) .and. (line(1:1)/="!")) then
  write (fmtchar, "(I3)") pos-1 
  read (line, "(A"//trim(fmtchar)//")", iostat = code_line ) line_head  
  line_head = trim(line_head)
  line(:) = line(pos+1:)   ! remove the heading from the line 
  !
  pos = Scan(line, "!")
  if (pos > 1) then 
  line(:) = trim(line(:pos-1))   ! remove any comments from the line
  end if 
  ! 
  select case (line_head)
   case ("degree of local Legendre basis in x")           ! ready to set up local order in variable x
    !!! First we read the paramters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the parameters we have just read !!! 
    if (m_count > 0) then 
     allocate (k_c_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (k_c_list)"
     stop
     end if
     !
     k_c_list=i_bulk(1:m_count)
    end if  
    if (slt==0) then 
    print *, "k_c_list=", k_c_list 
    end if
   case ("number of G-L boundary nodes in x")           ! ready to set up boundary order in variable x
    !!! First we read the paramters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count,2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the parameters we have just read !!! 
    if (m_count > 0) then 
     allocate (k_b_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (k_b_list)"
     stop
     end if
     !
     k_b_list=i_bulk(1:m_count)
    end if
    if (slt==0) then   
    print *, "k_b_list=", k_b_list 
    end if 
   case ("degree of local Legendre basis in t")           ! ready to set up local order in variable t
    !!! First we read the paramters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the parameters we have just read !!! 
    if (m_count > 0) then 
     allocate (d_c_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (d_c_list)"
     stop
     end if
     !
     d_c_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then  
    print *, "d_c_list=", d_c_list 
    end if 
   case ("number of G-L boundary nodes in t")           ! ready to set up boundary order in variable t
    !!! First we read the paramters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count,2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the parameters we have just read !!! 
    if (m_count > 0) then 
     allocate (d_b_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (d_b_list)"
     stop
     end if
     !
     d_b_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then  
    print *, "d_b_list=", d_b_list  
    end if
   case ("left endpoint in x")         ! ready to set up left endpoint in x
    !!! First we read the parameters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     x_L=r_bulk(1) ! all other inputs values are ignored 
    else 
     x_L=0.0_DP ! set up to the default value 
    end if
    if (slt==0) then   
    print *, "x_L=", x_L
    end if  
   case ("right endpoint in x")         ! ready to set up right endpoint in x
    !!! First we read the parameters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     x_R=r_bulk(1) ! all other input values are ignored 
    else 
     x_R=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
    print *, "x_R=", x_R
    end if 
   case ("left endpoint in u")         ! ready to set up left endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     u_L=r_bulk(1) ! all other inputs values are ignored 
    else 
     u_L=0.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
    print *, "u_L=", u_L
    end if 
   case ("right endpoint in u")         ! ready to set up right endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     u_R=r_bulk(1) ! all other input values are ignored 
    else 
     u_R=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
    print *, "u_R=", u_R
    end if 
   case ("left endpoint in v")         ! ready to set up left endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     v_L=r_bulk(1) ! all other inputs values are ignored 
    else 
     v_L=0.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
    print *, "v_L=", v_L
    end if 
   case ("right endpoint in v")         ! ready to set up right endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     v_R=r_bulk(1) ! all other input values are ignored 
    else 
     v_R=1.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
    print *, "v_R=", v_R
    end if 
   case ("left endpoint in w")         ! ready to set up left endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     w_L=r_bulk(1) ! all other inputs values are ignored 
    else 
     w_L=0.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
    print *, "w_L=", w_L
    end if 
   case ("right endpoint in w")         ! ready to set up right endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     w_R=r_bulk(1) ! all other input values are ignored 
    else 
     w_R=1.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
    print *, "w_R=", w_R
    end if 
   case ("uniform mesh in x")         ! ready to set up mesh in x is uniform parameter 
    !!! We read the parameter from the input line 
   if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     mesh_x_uniform = .TRUE.
    else
     mesh_x_uniform = .FALSE.
    end if 
    if (slt==0) then 
    print *, "mesh_x_uniform=", mesh_x_uniform
    end if 
   case ("uniform mesh in u")         ! ready to set up mesh in u is uniform parameter 
    !!! We read the parameter from the input line 
    if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     mesh_u_uniform = .TRUE.
    else
     mesh_u_uniform = .FALSE.
    end if
    if (slt==0) then 
    print *, "mesh_u_uniform=", mesh_u_uniform 
    end if 
   case ("uniform mesh in v")         ! ready to set up mesh in u is uniform parameter 
    !!! We read the parameter from the input line 
    if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     mesh_v_uniform = .TRUE.
    else
     mesh_v_uniform = .FALSE.
    end if
    if (slt==0) then 
    print *, "mesh_v_uniform=", mesh_v_uniform 
    end if 
   case ("uniform mesh in w")         ! ready to set up mesh in u is uniform parameter 
    !!! We read the parameter from the input line 
    if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     mesh_w_uniform = .TRUE.
    else
     mesh_w_uniform = .FALSE.
    end if
    if (slt==0) then 
    print *, "mesh_w_uniform=", mesh_w_uniform 
    end if 
   case ("number of cells in x")           ! ready to set up the number of cells in x
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (N_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (N_list)"
     stop
     end if
     !
     N_list=i_bulk(1:m_count)
    end if  
    if (slt==0) then 
    print *, "N_list=", N_list
    end if 
   case ("number of cells in u")           ! ready to set up the number of cells in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (Mu_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (Mu_list)"
     stop
     end if
     !
     Mu_list=i_bulk(1:m_count)
    end if
    if (slt==0) then   
    print *, "Mu_list=", Mu_list
    end if 
   case ("degree of local Lagrange basis in u")           ! ready to set up order in variable u 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (su_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (su_list)"
     stop
     end if
     !
     su_list=i_bulk(1:m_count)
    end if
    if (slt==0) then   
    print *, "su_list=", su_list 
    end if 
   case ("number of cells in v")           ! ready to set up the number of cells in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (Mv_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (Mv_list)"
     stop
     end if
     !
     Mv_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then  
    print *, "Mv_list=", Mv_list
    end if 
   case ("degree of local Lagrange basis in v")           ! ready to set up order in variable u 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (sv_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (sv_list)"
     stop
     end if
     !
     sv_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then  
    print *, "sv_list=", sv_list 
    end if 
   case ("number of cells in w")           ! ready to set up the number of cells in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (Mw_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (Mw_list)"
     stop
     end if
     !
     Mw_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then  
    print *, "Mw_list=", Mw_list
    end if 
   case ("degree of local Lagrange basis in w")           ! ready to set up order in variable u 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (sw_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (sw_list)"
     stop
     end if
     !
     sw_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then  
    print *, "sw_list=", sw_list 
    end if 
   case ("conditions on the left boundary")         ! ready to set up selected type of boundary conditions
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    selected_leftbcond=i_bulk(1) ! all other inputs values are ignored 
    else 
    selected_leftbcond=1 ! set up to the default value 
    end if 
    if (slt==0) then   
    print *, "selected_leftbcond=", selected_leftbcond  
    end if 
   case ("conditions on the right boundary")         ! ready to set up selected type of boundary conditions
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    selected_rightbcond=i_bulk(1) ! all other inputs values are ignored 
    else 
    selected_rightbcond=1 ! set up to the default value 
    end if 
    if (slt==0) then  
    print *, "selected_rightbcond=", selected_rightbcond  
    end if 
   case ("type of exact solution")         ! ready to set up selected type of boundary data
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    selected_exact_sol=i_bulk(1) ! all other inputs values are ignored 
    else 
    selected_exact_sol=1 ! set up to the default value 
    end if  
    if (slt==0) then  
    print *, "selected_exact_sol=", selected_exact_sol  
    end if
   case ("current solution save directory")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_solution_dir = "solution/"
    else
     current_solution_dir = Trim(Adjustl(line))
    end if 
    if (slt==0) then
    print *, "current_solution_dir=", current_solution_dir 
    end if
   case ("current solution base name")         ! ready to set base name (first part fo the name) to store various soltuion files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_solution_base_name = "bgk1d"
    else
     current_solution_base_name = Trim(Adjustl(line))
    end if 
    if (slt==0) then
    print *, "current_solution_base_name=", current_solution_base_name 
    end if 
   case ("current A operator base name")         ! ready to set base name (the first part of the name) of the file that has the A operator in it
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_Aoperator_base_name = "noname"
    else
     current_Aoperator_base_name = Trim(Adjustl(line))
    end if
    if (slt==0) then 
    print *, "current_Aoperator_base_name=", current_Aoperator_base_name  
    end if
   case ("initial time")       ! ready to set up the initial time
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     t_L=r_bulk(1) ! all other input values are ignored 
    else 
     t_L=0.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "t_L=", t_L
    end if 
   case ("final time")         ! ready to set up the end time
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     t_R=r_bulk(1) ! all other input values are ignored 
    else 
     t_R=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "t_R=", t_R
    end if 
   case ("instances to evaluate error")  ! ready to set up the number of error evaluations
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1000) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_eval_error=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_eval_error=1000 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "num_eval_error=", num_eval_error 
    end if 
   case ("instances to save solution")  ! ready to set up the number solution recordings
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 10) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_save_solution=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_save_solution=10 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "num_save_solution=", num_save_solution 
    end if 
   case ("type of nonuniform mesh in x")  ! ready to set up the type of nonuniform mesh in x
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    x_nonuniform_mesh_type=i_bulk(1) ! all other inputs values are ignored 
    else 
    x_nonuniform_mesh_type=1 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "x_nonuniform_mesh_type=", x_nonuniform_mesh_type 
    end if 
   case ("type of nonuniform mesh in u")  ! ready to set up the type of nonuniform mesh in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    u_nonuniform_mesh_type=i_bulk(1) ! all other inputs values are ignored 
    else 
    u_nonuniform_mesh_type=1 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "u_nonuniform_mesh_type=", u_nonuniform_mesh_type
    end if 
   case ("type of nonuniform mesh in v")  ! ready to set up the type of nonuniform mesh in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    v_nonuniform_mesh_type=i_bulk(1) ! all other inputs values are ignored 
    else 
    v_nonuniform_mesh_type=1 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "v_nonuniform_mesh_type=", v_nonuniform_mesh_type
    end if 
   case ("type of nonuniform mesh in w")  ! ready to set up the type of nonuniform mesh in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    w_nonuniform_mesh_type=i_bulk(1) ! all other inputs values are ignored 
    else 
    w_nonuniform_mesh_type=1 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "w_nonuniform_mesh_type=", w_nonuniform_mesh_type
    end if 
   case ("gauss order for moments in x")  ! ready to set up the order of gauss method in x variable for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_x_gauss_order=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_x_gauss_order=1 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "moments_x_gauss_order=", moments_x_gauss_order 
    end if 
   case ("mesh refinement in x for moments")  ! ready to set up the coefficient of mesh refinement in x for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_refine_x=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_refine_x=1 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "moments_refine_x=", moments_refine_x
    end if 
   case ("gauss order for moments in u")  ! ready to set up the order of gauss method in u variable for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_u_gauss_order=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_u_gauss_order=1 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "moments_u_gauss_order=", moments_u_gauss_order 
    end if 
   case ("mesh refinement in u for moments")  ! ready to set up the coefficient of mesh refinement in u for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_refine_u=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_refine_u=1 ! set up to the default value 
    end if  
    if (slt==0) then
    print *, "moments_refine_u=", moments_refine_u
    end if 
    ! we are ready to read the ordinary gas constant, reference teperature and viscosity and alpha constant 
    ! for the experiment --- will need in the collision term
  case ("the ordinary gas constant")         ! ready to set up the ordinary gas constant
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 1.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     gasR=r_bulk(1) ! all other inputs values are ignored 
    else 
     gasR=1.0_DP ! set up to the default value (garbage! gasR=1 has no meaning!!)
    end if 
    if (slt==0) then 
    print *, "gasR=", gasR
    end if 
    ! we are ready to read the gas reference temperature for the experiment --- will need in the collision term
  case ("gas reference temperature")         ! ready to set up gas reference temperature
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 1.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     gasTref=r_bulk(1) ! all other inputs values are ignored 
    else 
     gasTref=1.0_DP ! set up to the default value (garbage! gasTref=1 has no meaning!!)
    end if
    if (slt==0) then  
    print *, "gasTref=", gasTref 
    end if 
    ! we are ready to read the gas reference viscosity for the experiment --- will need in the collision term
  case ("gas reference viscosity")         ! ready to set up gas reference viscosity
        !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 1.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     gasmuref=r_bulk(1) ! all other inputs values are ignored 
    else 
     gasmuref=1.0_DP ! set up to the default value (garbage! gasmuref=1 has no meaning!!)
    end if 
    if (slt==0) then 
    print *, "gasmuref=", gasmuref 
    end if 
    ! we are ready to read the gas "alpha" in the gas state law for the experiment --- will need in the collision term
  case ("gas alpha constant")         ! ready to set up gas alpha constant
        !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 1.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     gasalpha=r_bulk(1) ! all other inputs values are ignored 
    else 
     gasalpha=1.0_DP ! set up to the default value (garbage! gasalpha=1 has no meaning!!)
    end if 
    if (slt==0) then 
    print *, "gasalpha=", gasalpha 
    end if 
    ! we are ready to read the temerature of the left wall. This is to be used with diffusive BCs 
  case ("temperature of the left wall")         ! ready to set up the temperature of the left wall
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 300.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     T_w_left=r_bulk(1) ! all other inputs values are ignored 
    else 
     T_w_left=300.0_DP ! set up to the default value (garbage! T_w_left=300K has no particular meaning!!)
    end if  
    if (slt==0) then
    print *, "T_w_left=", T_w_left
    end if 
    ! we are ready to read the temerature of the right wall. This is to be used with diffusive BCs 
  case ("temperature of the right wall")         ! ready to set up the temperature of the right wall
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 300.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     T_w_right=r_bulk(1) ! all other inputs values are ignored 
    else 
     T_w_right=300.0_DP ! set up to the default value (garbage! T_w_right=300K has no particular meaning!!)
    end if  
    if (slt==0) then
    print *, "T_w_right=", T_w_right
    end if 
case ("solution restart")         ! ready to set up mesh in x is uniform parameter 
    !!! We read the parameter from the input line 
    if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     need_to_restart = .TRUE.
    else
     need_to_restart = .FALSE.
    end if 
    if (slt==0) then
    print *, "need_to_restart=", need_to_restart
    end if 
case ("restart time")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     restart_time_txt = "0.00000000"
    else
     restart_time_txt = Trim(Adjustl(line))
    end if
    if (slt==0) then 
    print *, "restart_time_txt=", restart_time_txt 
    end if 
case ("cutoff radius of A-array")       ! ready to set up the cutoff radius for A-array
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     Trad=r_bulk(1) ! all other input values are ignored 
    else 
     Trad=100.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "Trad=", Trad 
    end if 
case ("error in integral in Chi")       ! ready to set up the error of integral evaluation in Chi in A-array
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     ErrChi=r_bulk(1) ! all other input values are ignored 
    else 
     ErrChi=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "ErrChi=", ErrChi 
    end if 
case ("error in integral in Epsilon")       ! ready to set up the error in integral in Epsilon for evaluation of A-array
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     ErrEps=r_bulk(1) ! all other input values are ignored 
    else 
     ErrEps=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "ErrEps=", ErrEps 
    end if 
case ("cutoff values for A-array")       ! ready to set up the cutoff values of A-array
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     min_sens=r_bulk(1) ! all other input values are ignored 
    else 
     min_sens=1.0_DP ! set up to the default value 
    end if  
    if (slt==0) then
    print *, "min_sens=", min_sens
    end if 
  case ("list of basis functions for A-array")           ! ready to set up the list of basis functions for evaluation of A-array 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (I1_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (I1_list)"
     stop
     end if
     !
     I1_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then 
    print *, "I1_list=", I1_list 
    end if 
 case ("number of OMP threads")  ! ready to set up the order of gauss method in x variable for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    Num_OMP_threads=i_bulk(1) ! all other inputs values are ignored 
    else 
    Num_OMP_threads=1 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "Num_OMP_threads=", Num_OMP_threads    
    end if 
!!!!! end added Alex 01/23/12 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
case ("error maxwell linearization")       ! ready to set up the cut off for linearization
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     linear_lev=r_bulk(1) ! all other input values are ignored 
    else 
     linear_lev=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "linear_lev=", linear_lev 
    end if 
case ("error maxwell decomposition")       ! ready to set up the cutoff valiue for non-linear decomposition
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     decomp_lev=r_bulk(1) ! all other input values are ignored 
    else 
     decomp_lev=100.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "decomp_lev=", decomp_lev 
    end if
case ("ref termal velocity")       ! ready to set up the reference thermal velocity unsed in dimensionless reduction
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     C_inf=r_bulk(1) ! all other input values are ignored 
    else 
     C_inf=1.0_DP ! set up to the default value 
    end if
    if (slt==0) then  
    print *, "C_inf=", C_inf     
    end if 
case ("ref characteristic length")       ! ready to set up the reference characteristic length used in dimensionless reduction
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     L_inf=r_bulk(1) ! all other input values are ignored 
    else 
     L_inf=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "L_inf=", L_inf     
    end if 
    if ((C_inf > 0) .and. (L_inf > 0)) then
     T_inf=L_inf/C_inf 
     if (slt==0) then 
     print *, "T_inf=", T_inf
     end if 
    else 
     T_inf=1.0_DP
     if (slt==0) then 
     print *, "T_inf=", T_inf
     end if 
    end if
case ("ref number of molecules")       ! ready to set up the reference total number of molecules in the volume used in dimensionless reduction
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     N_inf=r_bulk(1) ! all other input values are ignored 
    else 
     N_inf=1.0_DP ! set up to the default value 
    end if
    if (slt==0) then
    print *, "N_inf=", N_inf
    end if      
case ("ref molecular diameter")       ! ready to set up the molecular diameter for hard scpere model used in dimensionless reduction
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     mol_diam=r_bulk(1) ! all other input values are ignored 
    else 
     mol_diam=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "mol_diam=", mol_diam     
    end if 
case ("time step")         ! ready to set up the time step dt 
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     dt=r_bulk(1) ! all other input values are ignored 
    else 
     dt=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "dt=", dt
    end if     
   case ("number of chunks for Aarrays")  ! ready to set up the number of chunks in A arrays 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    numchnks=i_bulk(1) ! all other inputs values are ignored 
    else 
    numchnks=0 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "numchnks=", numchnks 
    end if 
   case ("number of procs for linear problem")  ! ready to set up the number of chunks in A arrays 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_lin_proc=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_lin_proc=1 ! set up to the default value 
    end if 
    if (num_lin_proc < 1) then
     num_lin_proc = 1 
     if (slt==0) then 
      print *,  "num_lin_proc=", num_lin_proc, "Specified value is too small. Reset to 1."
     end if
    end if 
    if (slt==0) then 
    print *, "num_lin_proc=", num_lin_proc 
    end if 
   case ("no update linearization steps wait")  ! ready to set up the number of chunks in A arrays 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_linstep_wait=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_linstep_wait=1 ! set up to the default value 
    end if 
    if (num_linstep_wait < 1) then
     num_linstep_wait =1 
     if (slt==0) then 
      print *,  "num_linstep_wait=", num_linstep_wait, "Specified value is too small. Reset to 1."
     end if
    end if 
    if (slt==0) then 
    print *, "num_linstep_wait=", num_linstep_wait
    end if 
   case ("number of mpi copies of A")  ! ready to set up the number of chunks in A arrays 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_Acopies=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_Acopies=0 ! set up to the default value 
    end if 
    if (slt==0) then 
    print *, "num_Acopies=", num_Acopies
    end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   case default
    if (slt==0) then 
    print *, "Can not process:" // line
    end if 
  end select
 else
 if (slt==0) then
 print *, line
 end if 
 end if   
 end do 
close(15)
end subroutine SetUWbgkParams 

subroutine ReadIntegersFromLine (line_head,line,i_bulk, m, defval)
character (len=*), intent (in) :: line_head          ! line with parameter name 
character (len=*), intent (in) :: line               ! line with numbers to be processed    
integer (I4B), intent (in) :: defval 
integer (I4B), dimension (:), intent (out) :: i_bulk ! storage for the varaibles 
integer (I4B), intent (out) :: m                           ! counter how many records has been read   
integer (I4B) :: code_line ! to use with Read statement 
!
character (len=10) :: fmtchar ! to keep format
!
code_line = 0
m=-1
do while (code_line == 0 )
   write (fmtchar,"(I3)") m
   m=m+1
   Read ( line, fmt = * ,  iostat = code_line) i_bulk(1:m+1)
   if ((code_line > 0) .and. (m==0)) then 
     print *, "Set1DHOVParameters: Error read line", code_line, & 
            "Parameter "//line_head//" may have an empty or improper numeric record. I use the default value", defval 
     i_bulk(m+1)=defval                       ! the default value for the parameter
     m=1
   end if 
end do
end subroutine ReadIntegersFromLine

subroutine ReadRealsFromLine (line_head,line,r_bulk, m, defval)
character (len=*), intent (in) :: line_head          ! line with parameter name 
character (len=*), intent (in) :: line               ! line with numbers to be processed    
real (DP), intent (in) :: defval 
real (DP), dimension (:), intent (out) :: r_bulk ! storage for the varaibles 
integer (I4B), intent (out) :: m                           ! counter how many records has been read   
integer (I4B) :: code_line ! to use with Read statement 
!
character (len=10) :: fmtchar ! to keep format
!
code_line = 0
m=-1
do while (code_line == 0 )
   write (fmtchar,"(I3)") m
   m=m+1
   Read ( line, fmt = * ,  iostat = code_line) r_bulk(1:m+1)
   if ((code_line > 0) .and. (m==0)) then 
     print *, "Set1DHOVParameters: Error read line", code_line, & 
            "Parameter "//line_head//" may have an empty or improper numeric record. I use the default value", defval 
     r_bulk(m+1)=defval                       ! the default value for the parameter
     m=1
   end if 
end do
end subroutine ReadRealsFromLine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteGridsDGV
! 
! This subroutine writes the grids_arrays on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WriteGridsDGV

use commvar, only: grids_u, grids_v, grids_w, grids_cap_u, grids_cap_v, grids_cap_w
!
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm ! is the scrap variable to keep te size of the cells arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_grids.dat"  ! this will keep the cells array of the problem
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
mm=size(grids_cap_u,1)                   
! we record the size of the arrays
write (15) mm
! now goes the error itself
write(15) grids_cap_u, grids_cap_v, grids_cap_w
!
mm=size(grids_u,1)                   
! we record the size of the arrays
write (15) mm
! now goes the error itself
write(15) grids_u
!
mm=size(grids_v,1)                   
!
write (15) mm
! now goes the error itself
write(15) grids_v
!
mm=size(grids_w,1)                   
!
write (15) mm
! now goes the error itself
write(15) grids_w
!
close (15)
!
end subroutine WriteGridsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteCellsDGV
! 
! This subroutine writes the cells_arrays on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WriteCellsDGV

use commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov
!
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm ! is the scrap variable to keep te size of the cells arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_cells.dat"  ! this will keep the cells array of the problem
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
mm=size(cells_pgrid,1)                   
! we record the size of the arrays
write (15) mm
! now goes the error itself
write(15) cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

!
close (15)
!
end subroutine WriteCellsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteNodesDGV
! 
! This subroutine writes the nodes arrays on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WriteNodesDGV

use commvar, only: nodes_u, nodes_v, nodes_w, nodes_gwts, & 
                   nodes_pcell, nodes_ui, nodes_vi, nodes_wi
!                   
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm ! is the scrap variable to keep te size of the cells arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_nodes.dat"  ! this will keep the cells array of the problem
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
mm=size(nodes_pcell,1)                   
! we record the size of the arrays
write (15) mm
! now goes the error itself
write(15)  nodes_u,nodes_v,nodes_w,nodes_gwts,&
           nodes_pcell,nodes_ui,nodes_vi, nodes_wi
!
close (15)
!
end subroutine WriteNodesDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadAarraysDGV
! 
! This subroutine reads A-arrays from the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ReadAarraysDGV

use commvar, only: A,A_capphi,A_xi,A_xi1,A_phi
!                   
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
real :: zz ! some scrap variable

! A quck check if the Arrays are already allocated
if (size(A,1)>0) then 
   print *,"ReadAarraysDGV: A arrays are already allocated. Exit."
   stop
end if    

! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperDGV(file_name)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)

! We now need to prepare the storage for the data.
allocate (A_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_capphi"
  stop
  end if
allocate (A_xi(1:nn), A_xi1(1:nn), A(1:nn), A_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
! now we open the file, again, to populate the A-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	A_capphi,A,A_xi,A_xi1,A_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
close (15)
!
zz=maxval(A)
zz=minval(A)
end subroutine ReadAarraysDGV
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadAarraysChnksDGV
! 
! This subroutine reads A-arrays from the disk when A-arrays are stored on the hard drive in chunks, 
! it assembled them 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ReadAarraysChnksDGV(nchnks)

use commvar, only: A,A_capphi,A_xi,A_xi1,A_phi
use dgvtools_mod
!                   
intrinsic Trim
!
integer (I4B), intent (in) :: nchnks ! the number of chunks of data  
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn,An,i ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
!
real (DP), dimension (:), allocatable :: Achnks
integer (I4B), dimension (:), allocatable :: Achnks_capphi, Achnks_xi, Achnks_xi1, Achnks_phi

! A quck check if the Arrays are already allocated
if (size(A,1)>0) then 
   print *,"ReadAarraysDGV: A arrays are already allocated. Exit."
   stop
end if    
!! Now the chunks of the A-arrays will be read and pieced together in the sequence imbedded in the filennames
!! IT IS IMPORTNAT THAT THE INFORMATION IN STORED IN THE CHUNKS IN THE RIGHT ORDER. 
!! In particular, A_capphi array will list the number of records for nodes. It will be expected that records in A 
!! are ordered in the increasing I1, or A_phi. 
!! 
!! The first chunk will go directly to A-arrays.... 
!! 
! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperChnksDGV(file_name,0)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)

! We now need to prepare the storage for the data.
allocate (A_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_capphi"
  stop
  end if
allocate (A_xi(1:nn), A_xi1(1:nn), A(1:nn), A_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
! now we open the file, again, to populate the A-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	A_capphi,A,A_xi,A_xi1,A_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
close (15)
!
! We have read the first chunk directly into A-arrays. 
!
do i=2,nchnks ! This array will loop until all chunks are read
 ! first, we prepare the file name to read the A-arrays
 call MakeBaseNameAoperChnksDGV(file_name,i-1)
 file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
 !
 ! now we open the file for reading and read some stuff from it: 
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 close (15)

 ! We now need to prepare the storage for the data.
 allocate (Achnks_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_capphi"
  stop
  end if
 allocate (Achnks_xi(1:nn), Achnks_xi1(1:nn), Achnks(1:nn), Achnks_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
 ! now we open the file, again, to populate the A-arrays
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 read (15, iostat=code_line)	Achnks_capphi,Achnks,Achnks_xi,Achnks_xi1,Achnks_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
 close (15)
 !
 ! now we need to record the new chunk into the A-arrays. For that we will need to extend the A-arrays (see below)
 ! first he A_capphi array -- because we expect that I1 do not overlap in chunks and if so, we expect the records be kept 
 ! in the increasing order of I1's --- I1s are stored in A_phi. Because of all that, we can just add Achnks_capphi to A_capphi:
 A_capphi = A_capphi + Achnks_capphi
 ! Now we need to extend the A arrays to add the new records:
 An=size(A,1) ! Let us find out how lond the arrays are now
 call ExtendAarraysDGV(An,nn)  ! we make room for the new records, An is the new size 
 ! Now we add the new reords
 A(An-nn+1:An)=Achnks
 A_xi(An-nn+1:An)=Achnks_xi
 A_xi1(An-nn+1:An)=Achnks_xi1
 A_phi(An-nn+1:An)=Achnks_phi
 ! now we need to destroy the chunks arrays:
 deallocate (Achnks,Achnks_xi,Achnks_xi1,Achnks_phi,Achnks_capphi)
 ! ready to reead the next chunk
end do 

end subroutine ReadAarraysChnksDGV
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadAarraysChnksOfstNrecDGV
! 
! This subroutine reads a given number of records of A-arrays from the disk when A-arrays are stored on 
! the hard drive in chunks. This subroutine will open chunk with the number fstchnk, it will skip the first 
! (ofst) records, and reads exactly numrec records into A arrays. If the chunk is finished but not eanough records have been read, 
! the software moves to the next chunk untill the desired number of records achieved or untill all chunks were read. 
!
! fstchnk -- the number of the first chunk 
! ofst  -- the number of records to skip in the first chunk (presumably these records are read at another process.
! numrec -- the total number of records to read into A-arrays, 
! nchnks -- the total number of chunks of A.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ReadAarraysChnksOfstNrecDGV(fstchnk,ofst,numrec,nchnks)

use commvar, only: A,A_capphi,A_xi,A_xi1,A_phi
use dgvtools_mod
!                   
intrinsic Trim
!
integer (I4B), intent (in) :: fstchnk ! the number of first chunks of data  
integer (I4B), intent (in) :: ofst ! the number of records to skip in the first chunk   
integer (I4B), intent (in) :: numrec ! the number of records that needs to be read  
integer (I4B), intent (in) :: nchnks ! the total number of chunks in A  
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn,i,j,Act, Act_need ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
!
real (DP), dimension (:), allocatable :: Achnks
integer (I4B), dimension (:), allocatable :: Achnks_capphi, Achnks_xi, Achnks_xi1, Achnks_phi

! A quck check if the Arrays are already allocated
if (size(A,1)>0) then 
   print *,"ReadAarraysChnksOfstNrecDGV: A arrays are already allocated. Exit."
   stop
end if    
!! Now the chunks of the A-arrays will be read and pieced together in the sequence imbedded in the filennames
!! IT IS IMPORTNAT THAT THE INFORMATION IN STORED IN THE CHUNKS IN THE RIGHT ORDER. 
!! In particular, A_capphi array will list the number of records for nodes. It will be expected that records in A 
!! are ordered in the increasing I1, or A_phi. 
!! 
!! The first chunk will go directly to A-arrays.... 
!! 
! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperChnksDGV(file_name,fstchnk-1)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)
! check if the read was successful:
if (code_line/=0) then 
 print *,  "ReadAarraysChnksOfstNrecDGV: read error for chunk of A-array. Possibly wrong name or damaged file. Stop"
 stop
end if

! We now need to prepare the storage for the data.
allocate (A_capphi(1:mm), Achnks_capphi(1:mm), stat = loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysChnksOfstNrecDGV: Allocation error for variable  A_capphi"
  stop
  end if
allocate (A_xi(1:numrec), A_xi1(1:numrec), A(1:numrec), A_phi(1:numrec), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysChnksOfstNrecDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
 allocate (Achnks_xi(1:nn), Achnks_xi1(1:nn), Achnks(1:nn), Achnks_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysChnksOfstNrecDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
! now we open the file, again, to populate the A-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	Achnks_capphi,Achnks,Achnks_xi,Achnks_xi1,Achnks_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
close (15)
!
! We have read the first chunk directly into Achnks-arrays.
! Now we skip the number of records ofst and put the rest in the permanent A-Arrays.
! first, determine the number of records to read: 
A_capphi = 0
Act = 0
Act_need = min(nn-ofst,numrec-Act)  ! this variable keeps track of how many records needs to be added to the array from this chunk of A 
A(1:Act_need)=Achnks(ofst+1:ofst+Act_need)
A_xi(1:Act_need)=Achnks_xi(ofst+1:ofst+Act_need)
A_xi1(1:Act_need)=Achnks_xi1(ofst+1:ofst+Act_need)
A_phi(1:Act_need)=Achnks_phi(ofst+1:ofst+Act_need)
Act=Act_need ! now we remember how many records we added to the A arrays... 
!
!  Now we need to write an A_capphi that corresponds to what is in A
i=1
do while ((ofst > sum(Achnks_capphi(1:i))) .and. (i<= mm)) ! we need to skip ofst records in this array.
 i=i+1
end do                            ! when this difference is negative, we went over the offset number of recods
if (i>mm) then                    ! a quick check if we had enough records 
 print *, "ReadAarraysChnksOfstNrecDGV: may be rerror in the work array -- no records in the chunk after offset taken"
 stop
end if
                                 
if (sum(Achnks_capphi(1:i)) - ofst < Act_need) then            ! now we need to see if the records for the basis function i is longer than the total number of 
 A_capphi(i) = sum(Achnks_capphi(1:i)) - ofst            !  records that will be read from this chunk. The appropriate number of records is acconted for.
 Act_need = Act_need - sum(Achnks_capphi(1:i)) + ofst
else 
 A_capphi(i) = Act_need  
 Act_need = 0
end if 
 
do while ((Act_need > 0) .and. (i<mm))
 i=i+1
 if (Achnks_capphi(i) >= Act_need) then 
   A_capphi(i) = Act_need
   Act_need = 0
 else 
   A_capphi(i) = Achnks_capphi(i) 
   Act_need = Act_need - Achnks_capphi(i)
 end if   
end do 
if ((i>=mm).and.(Act_need >0)) then                    ! a quick check if we had enough records 
 print *, "ReadAarraysChnksOfstNrecDGV: may be error in the work array -- no records in the chunk after offset taken"
 stop
end if     
deallocate (Achnks,Achnks_xi,Achnks_xi1,Achnks_phi,Achnks_capphi)                            
! We finished dealing with the first chunk.

! Now if we still have insufficent records in A, we will try to read the next chunk and us records from there..
j = fstchnk+1 ! point the index to the next chunk (remember that chunks are numbered from zero) 

do while ((Act < numrec) .and. (j <= nchnks)) ! This array will loop until enough records is read or until all chunks are read 
 ! first, we prepare the file name to read the A-arrays
 call MakeBaseNameAoperChnksDGV(file_name,j-1) ! the names of chunks are numbered from zero, therefore it is j-1
 file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
 !
 ! now we open the file for reading and read some stuff from it: 
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 close (15)
 
 ! We now need to prepare the storage for the data.
 allocate (Achnks_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_capphi"
  stop
  end if
 allocate (Achnks_xi(1:nn), Achnks_xi1(1:nn), Achnks(1:nn), Achnks_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
 ! now we open the file, again, to populate the A-arrays
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 read (15, iostat=code_line)	Achnks_capphi,Achnks,Achnks_xi,Achnks_xi1,Achnks_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
 close (15)
 !
 ! now we need to record the new chunk into the A-arrays. 
 ! We need to put only enough records to meet the total number of records requirement
 Act_need = min(nn, numrec - Act)  ! this calculates the number of records to be saved fron this chunk.  There is no offset aftert the first chunk
 A(Act+1:Act+Act_need) = Achnks(1:Act_need)
 A_xi(Act+1:Act+Act_need) = Achnks_xi(1:Act_need)
 A_xi1(Act+1:Act+Act_need) = Achnks_xi1(1:Act_need)
 A_phi(Act+1:Act+Act_need) = Achnks_phi(1:Act_need)
 Act=Act+Act_need ! now we remember how many records have been added to the A arrays so far
 ! No wwe need to update the A_caphi array. 
 i=0              ! reset the counter 
 do while ((Act_need > 0) .and. (i<mm)) ! will loop untill all records were taken account of or while ran out of records
  i=i+1
  if (Achnks_capphi(i) >= Act_need) then 
   A_capphi(i) = A_capphi(i) + Act_need 
   Act_need = 0
  else 
   A_capphi(i) = A_capphi(i) + Achnks_capphi(i) 
   Act_need = Act_need - Achnks_capphi(i)
  end if   
 end do
 if ((i>=mm).and.(Act_need >0)) then                    ! a quick check if we had enough records 
  print *, "ReadAarraysChnksOfstNrecDGV: may be rerror in the work array -- no records in the chunk after offset taken"
  stop
 end if     
 ! now we need to destroy the chunks arrays:
 deallocate (Achnks,Achnks_xi,Achnks_xi1,Achnks_phi,Achnks_capphi)
 ! ready to reead the next chunk
 j=j+1 ! Advance the number of chunk by one
end do 
! Next we need to set up the supplementary "shift" arrays that will be used for integration of the right side...... 
end subroutine ReadAarraysChnksOfstNrecDGV
!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteAarraysDGV
! 
! This subroutine writes the nodes arrays on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WriteAarraysDGV

use commvar, only: A,A_capphi,A_xi,A_xi1,A_phi
!                   
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn ! is the scrap variable to keep te size of the cells arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
mm=size(A_capphi,1)
nn=size(A,1)                   
! we record the size of the arrays
write (15) mm,nn
! now goes the error itself
write(15)  A_capphi,A,A_xi,A_xi1,A_phi
!                 
close (15)
!
end subroutine WriteAarraysDGV 
 
 
                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameDGV 
!
! This function creates a recognizable name for the file.
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters k_c, su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, x_nonuniform_mesh_type, u_nonuniform_mesh_type, mesh_x_uniform, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameDGV (file_name)

use commvar, only: k_c, su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_solution_base_name, &
          Mu,Mv,Mw,N,x_nonuniform_mesh_type,u_nonuniform_mesh_type,mesh_x_uniform,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform
!
intrinsic Trim, Min, Adjustl
!
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
file_name = trim(current_solution_dir)//trim(current_solution_base_name)
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k_c
file_name = trim(file_name)//trim(Adjustl(parchar))//"kc"
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"XU"
else
write (parchar, "(I1)") x_nonuniform_mesh_type
file_name = trim(file_name)//"X"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mu
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mu"
if (mesh_u_uniform) then 
 file_name = trim(file_name)//"UU"
else
 write (parchar, "(I1)") u_nonuniform_mesh_type
 file_name = trim(file_name)//"U"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mv
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mv"
if (mesh_v_uniform) then 
 file_name = trim(file_name)//"VU"
else
 write (parchar, "(I1)") v_nonuniform_mesh_type
 file_name = trim(file_name)//"V"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mw
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mw"
if (mesh_w_uniform) then 
 file_name = trim(file_name)//"WU"
else
 write (parchar, "(I1)") w_nonuniform_mesh_type
 file_name = trim(file_name)//"W"//trim(Adjustl(parchar))
end if 
!
end subroutine MakeBaseNameDGV
!!!!! End Modified Alex 10/18/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameChnksDGV 
!
! This function creates a recognizable name for the file. It should be used when data is writted in several chunks
! The file name will be formed by adding the chunk number 00X to the base name.
!
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters k_c, su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, x_nonuniform_mesh_type, u_nonuniform_mesh_type, mesh_x_uniform, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameChnksDGV (file_name,i)

use commvar, only: k_c, su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_solution_base_name, &
          Mu,Mv,Mw,N,x_nonuniform_mesh_type,u_nonuniform_mesh_type,mesh_x_uniform,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform
!
intrinsic Trim, Min, Adjustl
!
integer (I4B) :: i ! the index of the chunk (starts from zero)
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
write (parchar, "(I3.3)") i
file_name = trim(current_solution_dir)//trim(current_solution_base_name)//"ch"//trim(Adjustl(parchar))//"_"
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k_c
file_name = trim(file_name)//trim(Adjustl(parchar))//"kc"
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"XU"
else
write (parchar, "(I1)") x_nonuniform_mesh_type
file_name = trim(file_name)//"X"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mu
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mu"
if (mesh_u_uniform) then 
 file_name = trim(file_name)//"UU"
else
 write (parchar, "(I1)") u_nonuniform_mesh_type
 file_name = trim(file_name)//"U"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mv
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mv"
if (mesh_v_uniform) then 
 file_name = trim(file_name)//"VU"
else
 write (parchar, "(I1)") v_nonuniform_mesh_type
 file_name = trim(file_name)//"V"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mw
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mw"
if (mesh_w_uniform) then 
 file_name = trim(file_name)//"WU"
else
 write (parchar, "(I1)") w_nonuniform_mesh_type
 file_name = trim(file_name)//"W"//trim(Adjustl(parchar))
end if 
!
end subroutine MakeBaseNameChnksDGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteDGV2dArray(A,suffix)
!
! This subroutine dupms on the hard drive the 2D array. 
! A is the array of doubles 
! suffix is a the string containing a combination of letters that 
! will be attached to the name fo the file... 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteDGV2dArray(A,suffix)

!use commvar, only: nodes_u, nodes_v, nodes_w, nodes_gwts, & 
!                   nodes_pcell, nodes_ui, nodes_vi, nodes_wi
!                   
intrinsic Trim
!
real (DP), dimension (:,:), intent (in) :: A ! the incoming array of doubles
character (len=10) :: suffix ! the string to hold the suffix
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: m1,m2 ! scrap variables to keep te size of the 2d arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_"//trim(Adjustl(suffix))//".dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
m1=size(A,1)
m2=size(A,2)                   
! we record the sizes of the array
write (15) m1,m2
! now goes the array itself
write(15)  A
!
close (15)
!
end subroutine WriteDGV2dArray

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteI1DGV
!
! This subroutine will write a sequence of node indices to be later used in the evaluation of A operator. 
! 
! I1 is the array of integers that needs to be written on the disk. Notice that the first element in the array 
!  I1(1) contains the total number of records. So then only 2..I1(1)+1 cells in the array are importnat   
!
! u,v,w are the components of the velocity point. Am active cell in velocity space will be found that containes 
! this velocity point and indices of all nodes form this cell will be listed.  
!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteI1DGV(u,v,w)
!
use dgvtools_mod
!@
intrinsic Trim
!!!!!
real (DP), intent (in) :: u,v,w !  the components of the velocity
!!!!!!!!!!!!!!!!!!!!!!
!
integer (I4B) :: celli ! the number of the cell that contains point (u,v,w)
character (len=132) :: file_name ! the variable to store the file name 
integer, parameter :: nnn=5 !! Number of entries in the row 
integer (I4B), dimension (300) :: I1 ! an array to store the found nodes' indices.
integer (I4B) :: j,k,n ! some counters  
character (len=20) :: parchar    ! string to keep char value of the paameter 
character (len=10) :: FMT1 ! the format string 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we find the number of the active cell that contains point (u,v,w)
celli = FindCellContainsPoint_DGV(u,v,w)
! we find the indices of all nodes that belong to cell with the number icell
call FindI1sByCellNum_DGV(celli,I1)
! Now it is time to record this infomation 

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_I1"
write (parchar, "(F4.1)") u
file_name = trim(file_name)//trim(Adjustl(parchar))//"u"
write (parchar, "(F4.1)") v
file_name = trim(file_name)//trim(Adjustl(parchar))//"v"
write (parchar, "(F4.1)") w
file_name = trim(file_name)//trim(Adjustl(parchar))//"w.txt"  ! this file will keep the array 
!
! now we open the file in formatted write regime for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the numbers of nodes belonging to the cell containing velocity point:"
write(15, *) "u=", u, "v=", v, "w=", w
write(15, *) "I1= "
!prepare the format string to print nnn entries in a row
write (parchar, "(I)") nnn
FMT1="("//trim(Adjustl(parchar))//"(I5, A))" 
! next we print array using nnn entries in a row
n=I1(1)
k=0
do while (k + nnn < n)
  write (15, FMT1) (I1(1+j), "," , j=k+1,k+nnn) 
  k=k+nnn
end do 
write (15,"(3(I5, "","" ))") I1(2+k:1+n)
close (15)
end subroutine WriteI1DGV 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteSol_DGV(suff)
!
! This subroutine will damp on the disk the current state of the solution along with some 
! parameters.
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteSol_DGV (suff) 

use commvar, only:  f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5,&
                    e_time,dt,rkmts

intrinsic Trim
!
character (len=15) :: suff ! the string to hold the suffix
!
character (len=132) :: file_name ! the variable to store the file name 

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(suff))//"_sol.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")

write (15) e_time,dt ! record the current time and the used dt
write (15) rkmts ! record the order of the time integrator.
write (15) size(f,1) ! record the length of the arrays.
! Next we save the solution .... 
SELECT CASE(rkmts)
	CASE(1)
write (15) f
    CASE(2)
write (15) f,f1,frhs1
    CASE(3)
write (15) f,f1,f2,frhs1,frhs2
    CASE(4)
write (15) f,f1,f2,f3,frhs1,frhs2,frhs3
    CASE(5)
write (15) f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4
    CASE default 
PRINT *, "Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT    
! end save 
close (15)
end subroutine WriteSol_DGV  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteSol_init_DGV(suff)
!
! THIS IS A COPY OF THE ABOVE SUBROUTINE THAT ONLY SAVES f 
!
! SHOUD ONLY BE USED TO SAVE THE INITIAL DATA>>>
!
! This subroutine will damp on the disk the current state of the solution along with some 
! parameters.
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteSol_init_DGV (suff) 

use commvar, only:  f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5,&
                    e_time,dt,rkmts

intrinsic Trim
!
character (len=15) :: suff ! the string to hold the suffix
!
character (len=132) :: file_name ! the variable to store the file name 

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(suff))//"_sol.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")

write (15) e_time,dt ! record the current time and the used dt
write (15) rkmts ! record the order of the time integrator.
write (15) size(f,1) ! record the length of the arrays.
! Next we save the solution .... 
SELECT CASE(rkmts)
	CASE(1)
write (15) f
    CASE(2)
write (15) f
    CASE(3)
write (15) f
    CASE(4)
write (15) f
    CASE(5)
write (15) f
    CASE default 
PRINT *, "Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT    
! end save 
close (15)
end subroutine WriteSol_init_DGV  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadSol_DGV(suff)
!
! This subroutine will read the solution from the disk.
! the solution 
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadSol_DGV (suff) 

use commvar, only:  f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5,&
                    e_time,dt,rkmts

intrinsic Trim
!
character (len=15) :: suff ! the string to hold the suffix
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: m ! the length of the arrays
integer (I4B) :: loc_alloc_stat ! some dump varaible


! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(suff))//"_sol.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")

read (15) e_time,dt ! record the current time and the used dt
read (15) rkmts ! record the order of the time integrator.
read (15) m ! record the length of the arrays.
! A quick check that the storages are the correct size: 
if (m /= size(f,1)) then 
 print *, "ReadSol_DGV: trying to read a storage of wrong size, m/= size(f,1)"
 close (15)
 stop
end if  
! Next we save the solution .... 
SELECT CASE(rkmts)
	CASE(1)
read (15) f ! the space for this one needs to be allocated already ... 
    CASE(2)
! We now need to prepare the storage for the data.
allocate (f1(1:m), frhs1(1:m), frhs2(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,frhs1,frhs2"
  close(15)
  stop
  end if
read (15) f,f1,frhs1
    CASE(3)
! We now need to prepare the storage for the data.
allocate (f1(1:m),f2(1:m),frhs1(1:m),frhs2(1:m),frhs3(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,f2,frhs1,frhs2,frhs3"
  close(15)
  stop
  end if
read (15) f,f1,f2,frhs1,frhs2
    CASE(4)
  allocate (f1(1:m),f2(1:m),f3(1:m),frhs1(1:m),frhs2(1:m),frhs3(1:m),frhs4(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,f2,f3,frhs1,frhs2,frhs3,frhs4"
  close(15)
  stop
  end if
read (15) f,f1,f2,f3,frhs1,frhs2,frhs3
    CASE(5)
  allocate (f1(1:m),f2(1:m),f3(1:m),f4(1:m),frhs1(1:m),frhs2(1:m),frhs3(1:m),frhs4(1:m),frhs5(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5"
  close(15)
  stop
  end if
read (15) f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4
    CASE default 
PRINT *, "Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT    
! end save 
close (15)
end subroutine ReadSol_DGV  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErr_DGV (ndens_a,ubar_a,vbar_a,wbar_a,tempr_a)
!
! This subroutine will damp on the disk some error arrays
! the solution 
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteErr_DGV (time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a)
!
intrinsic Trim
!!!!!
real (DP), dimension (:), intent (in) :: time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a ! these arrays contain time when the errors are recorder and the errors themselves
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: n ! some counters  
character (len=10) :: FMT1 ! the format string 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_macroerr.txt"  ! this file will keep the array 
!
! now we open the file in formatted write regime for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the values of density, average velocity and temperature, and the time stamp:"
write(15, *) "Time		Ndens		Ubar		Vbar		Wbar"
!prepare the format string to print nnn entries in a row
FMT1="(F18.15 "","" F18.5 "","" F18.5 "","" F18.5 "","" F18.5 "","" F18.5)" 
! next we print array using nnn entries in a row
do n=1,size(ndens_a,1)
  write (15,"(F18.7, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15)" ) time_a(n), "," , ndens_a(n),",", &
              ubar_a(n),",", vbar_a(n),",", wbar_a(n),",", tempr_a(n)
end do 
close (15)
end subroutine WriteErr_DGV 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErrPlus_DGV (ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a)
!
! This subroutine will damp on the disk some error arrays
! the solution 
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteErrPlus_DGV (time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a)
!
intrinsic Trim
!!!!!
real (DP), dimension (:), intent (in) :: time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a ! these arrays contain time when the errors are recorder and the errors themselves
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: n ! some counters  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_macroerr.txt"  ! this file will keep the array 
!
! now we open the file in formatted write regime for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the values of density, average velocity and temperature, and the time stamp:"
write(15, *) "Time		Ndens		Ubar		Vbar		Wbar    T     Tx     Ty    Tz"
! next we print array using nnn entries in a row
do n=1,size(ndens_a,1)
  write (15,"(F18.14, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15,  A2, F18.15)" ) &
        time_a(n), "," , ndens_a(n),",", ubar_a(n),",", vbar_a(n),",", wbar_a(n),",", tempr_a(n), &
        ",", tempr_u_a(n),",", tempr_v_a(n),",", tempr_w_a(n)
end do 
close (15)
end subroutine WriteErrPlus_DGV 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErrHighCMoments_DGV (ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a)
!
! This subroutine will damp on the disk some error arrays
! the solution 
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteErrHighCMoments_DGV (time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a,&
                              mom3_u_a,mom3_v_a,mom3_w_a,mom4_u_a,mom4_v_a,mom4_w_a,mom5_u_a,mom5_v_a,&
                              mom5_w_a,mom6_u_a,mom6_v_a,mom6_w_a)
!
intrinsic Trim
!!!!!
real (DP), dimension (:), intent (in) :: time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a ! these arrays contain time when the errors are recorder and the errors themselves
real (DP), dimension (:), intent (in) :: mom3_u_a,mom3_v_a,mom3_w_a,mom4_u_a,mom4_v_a,mom4_w_a,mom5_u_a,mom5_v_a
real (DP), dimension (:), intent (in) :: mom5_w_a,mom6_u_a,mom6_v_a,mom6_w_a
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: n ! some counters  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_macroerr.txt"  ! this file will keep the array 
!
! now we open the file in formatted write regime for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the values of density, average velocity and temperature, and the time stamp:"
write(15, *) "Time , Ndens , Ubar , Vbar , Wbar , T , Tx , Ty , Tz , (u-ubar)^3 ,(v-vbar)^3 , (w-wbar)^3 , ",& 
     "(u-ubar)^4 , (v-vbar)^4 , (w-wbar)^4, (u-ubar)^5 , (v-vbar)^5 , (w-wbar)^5 , (u-ubar)^6 , (v-vbar)^6 , (w-wbar)^6" 
! next we print array using nnn entries in a row
do n=1,size(ndens_a,1)
  write (15,"(F18.14, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15,  A2, F18.15, "//&
            "A2, F18.14, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2,"//&
            "F18.14, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15 )") &
        time_a(n), "," , ndens_a(n),",", ubar_a(n),",", vbar_a(n),",", wbar_a(n),",", tempr_a(n), &
        ",", tempr_u_a(n),",",tempr_v_a(n),",", tempr_w_a(n),",",mom3_u_a(n),",",mom3_v_a(n), &
        ",", mom3_w_a(n),",",mom4_u_a(n),",", mom4_v_a(n),",",mom4_w_a(n),",", mom5_u_a(n),&
        ",", mom5_v_a(n),",",mom5_w_a(n),",", mom6_u_a(n),",",mom6_v_a(n),",", mom6_w_a(n)
end do 
close (15)
end subroutine WriteErrHighCMoments_DGV 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KeepTrackL1_err_DGV(L1_err)
! 
! This subroutine saves the values of the relative L1 norm of the deviation of the current solution from the local maxwellian 
! 
! Takes the parameters from commvar.mod
!!!!!!!!!!!!!!!!!!

subroutine KeepTrackL1_err_DGV(L1_err)

use commvar, only: L1_a, L1_t, e_time, L1_count, t_L, t_R,num_save_solution,&
                    L1_record_period,L1_next_time_record,L1_err_eval_period,L1_next_time_eval
intrinsic Trim                    

real (DP), intent (in) :: L1_err ! relative L1 error of the deviation of distribution function from the local maxwellian computed elsewhere
integer (I4B), parameter :: c = 1000 ! the total number of records in the error file
integer :: loc_alloc_stat
integer (I4B) :: n ! scrap counter
character (len=132) :: file_name ! the variable to store the file name 

!!!!!!!!!!!!!!!!!!!!!!!!
! check is the error arrays are allocated. if not -- allocate them and reset the counter
if (size(L1_a,1)<2) then 
 allocate (L1_a(1:c),L1_t(1:c),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "KeepTrackL1_err_DGV: Allocation error for variable (L1_a(1:c),L1_t(1:c))"
  stop
  end if
L1_count = 1
L1_a(1) = L1_err
L1_t(1) = e_time 
!!!!!!!!
! we also set up the time intervals for evaluation and saves
L1_record_period = (t_R-t_L)/Real(num_save_solution,DP)
L1_next_time_record = e_time + L1_record_period
L1_err_eval_period = (t_R-t_L)/Real(c,DP)
L1_next_time_eval = e_time + L1_err_eval_period
end if 
! once the arrays are allocated, we periodically save the solution 
if (e_time >= L1_next_time_eval) then ! check if it is time to record
if (L1_count < c) then
       L1_count=L1_count+1 
       L1_next_time_eval = e_time + L1_err_eval_period
       L1_a(L1_count) = L1_err
       L1_t(L1_count) = e_time 
end if       
end if 
if (e_time >= L1_next_time_record) then ! check if it is time to save on disk.
 L1_next_time_record = e_time + L1_record_period 
 ! Writing the error on disk: 
 ! first, we prepare the file name to store the solution
 call MakeBaseNameDGV(file_name)
 file_name = trim(Adjustl(file_name))//"_L1err.txt"  ! this file will keep the array 
 !
 ! now we open the file in formatted write regime for record and save some stuff in it: 
 open (15, file=file_name, position ="REWIND", action="WRITE", &
                    form="FORMATTED")  
 write(15, *) "This file contains the values of relative L1 error of deviation from local maxwellian and the time stamp:"
 write(15, *) "Time		L1_err"
 ! next we print array using nnn entries in a row
 do n=1,L1_count
   write (15,"(F18.14, A2, F18.15)" ) L1_t(n), "," , L1_a(n)  
 end do 
 close (15)
end if 
end subroutine KeepTrackL1_err_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameAoperDGV 
!
! This function creates a recognizable name for the file to store A-operator in one file.
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters k_c, su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, x_nonuniform_mesh_type, u_nonuniform_mesh_type, mesh_x_uniform, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameAoperDGV (file_name)

use commvar, only: k_c, su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_Aoperator_base_name, &
          Mu,Mv,Mw,N,x_nonuniform_mesh_type,u_nonuniform_mesh_type,mesh_x_uniform,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform
!
intrinsic Trim, Min, Adjustl
!
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
file_name = trim(current_solution_dir)//trim(current_Aoperator_base_name)
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k_c
file_name = trim(file_name)//trim(Adjustl(parchar))//"kc"
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"XU"
else
write (parchar, "(I1)") x_nonuniform_mesh_type
file_name = trim(file_name)//"X"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mu
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mu"
if (mesh_u_uniform) then 
 file_name = trim(file_name)//"UU"
else
 write (parchar, "(I1)") u_nonuniform_mesh_type
 file_name = trim(file_name)//"U"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mv
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mv"
if (mesh_v_uniform) then 
 file_name = trim(file_name)//"VU"
else
 write (parchar, "(I1)") v_nonuniform_mesh_type
 file_name = trim(file_name)//"V"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mw
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mw"
if (mesh_w_uniform) then 
 file_name = trim(file_name)//"WU"
else
 write (parchar, "(I1)") w_nonuniform_mesh_type
 file_name = trim(file_name)//"W"//trim(Adjustl(parchar))
end if 
!
end subroutine MakeBaseNameAoperDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameAoperChnksDGV 
!
! This function creates a recognizable name for the files to store chunks of operator A. It should be used when data is writted in several chunks
! The file name will be formed by adding the chunk number 00X to the base name.
!
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters k_c, su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, x_nonuniform_mesh_type, u_nonuniform_mesh_type, mesh_x_uniform, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameAoperChnksDGV (file_name,i)

use commvar, only: k_c, su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_Aoperator_base_name, &
          Mu,Mv,Mw,N,x_nonuniform_mesh_type,u_nonuniform_mesh_type,mesh_x_uniform,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform
!
intrinsic Trim, Min, Adjustl
!
integer (I4B) :: i ! the index of the chunk (starts from zero)
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
write (parchar, "(I3.3)") i
file_name = trim(current_solution_dir)//trim(current_Aoperator_base_name)//"ch"//trim(Adjustl(parchar))//"_"
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k_c
file_name = trim(file_name)//trim(Adjustl(parchar))//"kc"
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"XU"
else
write (parchar, "(I1)") x_nonuniform_mesh_type
file_name = trim(file_name)//"X"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mu
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mu"
if (mesh_u_uniform) then 
 file_name = trim(file_name)//"UU"
else
 write (parchar, "(I1)") u_nonuniform_mesh_type
 file_name = trim(file_name)//"U"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mv
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mv"
if (mesh_v_uniform) then 
 file_name = trim(file_name)//"VU"
else
 write (parchar, "(I1)") v_nonuniform_mesh_type
 file_name = trim(file_name)//"V"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mw
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mw"
if (mesh_w_uniform) then 
 file_name = trim(file_name)//"WU"
else
 write (parchar, "(I1)") w_nonuniform_mesh_type
 file_name = trim(file_name)//"W"//trim(Adjustl(parchar))
end if 
!
end subroutine MakeBaseNameAoperChnksDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ScanAarrsChnks4Acapphi(A_capphi_aggreg,Act_aggreg)
! 
! This subroutine reads A-arrays from the disk when A-arrays are stored on the hard drive in chunks, 
! It determines the number of records of A arrays stored in each chunk. This information is recorded in 
! the array chnks_Act
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ScanAarrsChnks4Acapphi(chunks_Act)

use dgvtools_mod
!                   
intrinsic Trim
!
integer (I4B), dimension (:), intent (out) :: chunks_Act
!

character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn,i ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) :: nchnks
!
integer (I4B), dimension (:), allocatable :: Achnks_capphi

nchnks=size(chunks_Act,1) ! determine the number of chunks

!! Now the chunks of the A-arrays will be read and their A_capphi pieced together
!! IT IS IMPORTNAT THAT THE INFORMATION IN STORED IN THE CHUNKS IN THE RIGHT ORDER. 
!! 
! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperChnksDGV(file_name,0)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)

! We now need to prepare the storage for the data.
allocate (Achnks_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ScanAarrsChnks4Acapphi: Allocation error for variables Achnks_capphi"
  stop
  end if

! now we open the file, again, to populate the A_capphi arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	Achnks_capphi ! read the array		A_capphi (other arrays should be still there : A,A_xi,A_xi1,A_phi)
close (15)
!
! We have scanned the first chunk of A-arrays and read the A_capphi for the first chunk. Now we need to move it to A_capphi_aggreg
chunks_Act(1)= sum(Achnks_capphi)
! 
! Now we need to proceed with other chuncks. 

do i=2,nchnks ! This array will loop until all chunks are read
 ! first, we prepare the file name to read the A-arrays
 call MakeBaseNameAoperChnksDGV(file_name,i-1)
 file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
 !
 ! now we open the file for reading and read some stuff from it: 
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 read (15, iostat=code_line)	Achnks_capphi ! read the arrays A_capphi					
 close (15)
 !
 ! now we need to add the number of records to A_capphi . 
 chunks_Act(i) = sum(Achnks_capphi)
 ! ready to reead the next chunk
end do 
deallocate (Achnks_capphi)
!
end subroutine ScanAarrsChnks4Acapphi
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RestartSolutionDGV
!
! This procedure reads a solution from the disk. 
!
! This subroutine allocates the storage for the solution  
!
! the directory path is set through parameter current_solution_dir
! the name is composed from the curr_sol_base_name and a values of different parameters 
!
! This procedure is dependent on the main program! Do not call before all parameters are 
! set especially f -- this has to have been allocated correctly, its dimensions will be used to 
! set dimensions of other arrays. If the dimensions are wrong, the arrays will not read correctly  
! 
!
! most of the variables are taken directly from 
! commvar.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RestartSolutionDGV(restart_time_txt)
use commvar, only: f, f1, f2, f3, f4, frhs1, frhs2, frhs3, frhs4, &
							frhs5, rkmts,dt,e_time
!
intrinsic Trim, Adjustl
!
character (len=20) :: restart_time_txt ! the time at which the solution need to be restored
character (len=132) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the parameter 
integer (I4B) :: dump1,dump2     ! string to read some garbage bits from the file
integer :: m ! local counters
integer :: code_line !! local variable to keep the status of reading from file
integer :: loc_alloc_stat ! to keep allocation status
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! First, we need to allocate the storages
SELECT CASE(rkmts)
	CASE(1)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:1: Allocation error for variable (frhs1)"
		stop
	  END IF
	CASE(2)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),frhs2(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:2: Allocation error for variable (frhs1,frhs2,f1)"
		stop
	  END IF
	CASE(3)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:3: Allocation error for variable (frhs1,f1)"
		stop
	  END IF
	ALLOCATE(frhs2(1:size(f,1)),frhs3(1:size(f,1)),f2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:4: Allocation error for variable (frhs2,frhs3,f2)"
		stop
	  END IF
	CASE(4) 
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:5: Allocation error for variable (frhs1,f1)"
		stop
	  END IF
	ALLOCATE(frhs2(1:size(f,1)),f2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:6: Allocation error for variable (frhs2,f2)"
		stop
	  END IF
	  ALLOCATE(frhs3(1:size(f,1)),frhs4(1:size(f,1)),f3(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:7: Allocation error for variable (frhs3,frhs4,f3)"
		stop
	  END IF
    CASE(5)
      !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:8: Allocation error for variable (frhs1,f1)"
		stop
	  END IF
	  ALLOCATE(frhs2(1:size(f,1)),f2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:9: Allocation error for variable (frhs2,f2)"
		stop
	  END IF
	ALLOCATE(frhs3(1:size(f,1)),f3(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:10: Allocation error for variable (frhs3,f3)"
		stop
	  END IF
	  ALLOCATE(frhs4(1:size(f,1)),frhs5(1:size(f,1)),f4(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionDGV:11: Allocation error for variable (frhs4,frhs5,f4)"
		stop
	  END IF
	CASE default
			PRINT *, "RestartSolutionDGV:12: The value of (rkmts) must be from 1 to 5. No such RK or MTS methods implemented"
			stop
END SELECT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Second, we prepare the file name to read the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(restart_time_txt))//"_sol.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")

read (15) e_time,dt ! record the current time and the used dt
read (15) rkmts ! record the order of the time integrator.
read (15) m ! record the length of the arrays.
! A quick check that the storages are the correct size: 
if (m /= size(f,1)) then 
 print *, "RestartSolutionDGV:13: trying to read a storage of wrong size, m/= size(f,1)"
 close (15)
 stop
end if  
! Next we save the solution .... 
SELECT CASE(rkmts)
	CASE(1)
read (15) f ! the space for this one needs to be allocated already ... 
    CASE(2)
read (15) f,f1,frhs1
    CASE(3)
read (15) f,f1,f2,frhs1,frhs2
    CASE(4)
read (15) f,f1,f2,f3,frhs1,frhs2,frhs3
    CASE(5)
read (15) f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4
    CASE default 
PRINT *, "RestartSolutionDGV:14:Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT    
! end save 
close (15)
!
end subroutine RestartSolutionDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteMomDerivRelaxTimes_DGV 
!
! This subroutine will dump on the disk values fo the derivatives for the moments and values of the relaxation times for the same moments. 
! there will be two files created. The file with the suffix _dermom will have derivatives of the moments
! and the file with the suffix _reltimes will contain the relaxation times 
!  
! the file that contains the data is accessed directly from the commvar
!
! This subroutine expects that 20 moments are recorded in all. There will be a quick chsch for that.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteMomDerivRelaxTimes_DGV 

use commvar, only: MomDerivRelaxT,MomDerivRelaxT_flag 
!!!!
intrinsic Trim
!!!!!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: i,n,k ! some counters  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first we check is the array exists. Otherwise we send our an error message and stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (MomDerivRelaxT_flag /= 1) then 
 print *, "WriteMomDerivRelaxTimes_DGV: the flag for array MomDerivRelaxT indicates array does not exist. stop"
 stop
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! second check: we will try to read 40 + 1 variables from the array. We need to check that the information 
! is set up correspondingly. The second cell in the array MomDerivRelaxT gives the number of 
k = int(MomDerivRelaxT(2))
if (k /= 41) then 
 print *, "WriteMomDerivRelaxTimes_DGV: number of variables in each record is inconsistent with the expected value. stop"
 stop
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will write two files. The first file will contain the values of the moment's derivatives:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first the file with the derivatives of the moments.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the derivatives of the moments.
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_dermom.txt"  ! this file will keep the array 
!
! now we open the file in formatted write regime for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This file contains values of time derivatives of density, average velocity and temperature, and the time stamp:"
write(15, *) "Time , Ndens , Ubar , Vbar , Wbar , T , Tx , Ty , Tz , (u-ubar)^3 ,(v-vbar)^3 , (w-wbar)^3 , ",& 
     "(u-ubar)^4 , (v-vbar)^4 , (w-wbar)^4, (u-ubar)^5 , (v-vbar)^5 , (w-wbar)^5 , (u-ubar)^6 , (v-vbar)^6 , (w-wbar)^6" 
! next we print array using nnn entries in a row
do i=0,int(MomDerivRelaxT(1))-1
  n=2+i*k
  write (15,"(E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, "//&
            "A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2,"//&
            "E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6 )") &
        MomDerivRelaxT(n+1), "," , MomDerivRelaxT(n+2), "," , MomDerivRelaxT(n+3), "," , MomDerivRelaxT(n+4), "," , &
        MomDerivRelaxT(n+5), "," , MomDerivRelaxT(n+6), "," , MomDerivRelaxT(n+7), "," , MomDerivRelaxT(n+8), "," , &
        MomDerivRelaxT(n+9),  "," , MomDerivRelaxT(n+10), "," , MomDerivRelaxT(n+11), "," , MomDerivRelaxT(n+12), "," , &
        MomDerivRelaxT(n+13), "," , MomDerivRelaxT(n+14), "," , MomDerivRelaxT(n+15), "," , MomDerivRelaxT(n+16), "," , &
        MomDerivRelaxT(n+17), "," , MomDerivRelaxT(n+18), "," , MomDerivRelaxT(n+19), "," , MomDerivRelaxT(n+20), "," , &
        MomDerivRelaxT(n+21)
end do 
close (15)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! second the file with the relaxation times 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the derivatives of the moments.
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_relaxtimes.txt"  ! this file will keep the array 
!
! now we open the file in formatted write regime for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This file contains values of relaxation times for density, average velocity and temperature, and the time stamp:"
write(15, *) "Time , Ndens , Ubar , Vbar , Wbar , T , Tx , Ty , Tz , (u-ubar)^3 ,(v-vbar)^3 , (w-wbar)^3 , ",& 
     "(u-ubar)^4 , (v-vbar)^4 , (w-wbar)^4, (u-ubar)^5 , (v-vbar)^5 , (w-wbar)^5 , (u-ubar)^6 , (v-vbar)^6 , (w-wbar)^6" 
! next we print array using nnn entries in a row
do i=0,int(MomDerivRelaxT(1))-1
  n=2+i*k
  write (15,"(E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6,  A2, E18.6, "//&
            "A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2,"//&
            "E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6, A2, E18.6 )") &
        MomDerivRelaxT(n+1), "," , MomDerivRelaxT(n+22), "," , MomDerivRelaxT(n+23), "," , MomDerivRelaxT(n+24), "," , &
        MomDerivRelaxT(n+25), "," , MomDerivRelaxT(n+26), "," , MomDerivRelaxT(n+27), "," , MomDerivRelaxT(n+28), "," , &
        MomDerivRelaxT(n+29), "," , MomDerivRelaxT(n+30), "," , MomDerivRelaxT(n+31), "," , MomDerivRelaxT(n+32), "," , &
        MomDerivRelaxT(n+33), "," , MomDerivRelaxT(n+34), "," , MomDerivRelaxT(n+35), "," , MomDerivRelaxT(n+36), "," , &
        MomDerivRelaxT(n+37), "," , MomDerivRelaxT(n+38), "," , MomDerivRelaxT(n+39), "," , MomDerivRelaxT(n+40), "," , &
        MomDerivRelaxT(n+41)
end do 
close (15)
end subroutine WriteMomDerivRelaxTimes_DGV 


end module readwrite
