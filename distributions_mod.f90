!! THis module makes a bimodal distribution using two maxwellians 
!! one maxwellian on the left and another on the right and some sort of linear? switching in between... 
!! IN 3D in velocity space and 1D in physical space...

!DIMENSIONLESS SEE NOTES FOR THE REDUCTION FORMULAS
module distributions_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 

implicit none
   interface maxwelveldist
     module procedure maxwelveldist, maxwelveldist_T_vector, maxwelveldist_u_vectors, &
                                     maxwelveldist_T_vector_u_vectors
   end interface
   
!!!!!!!!!!!!! parameters and global variables 
real (DP), parameter, private :: pi25DT = 3.141592653589793238462643d0
!!!!!!!!!!


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! maxwelveldist (T,u_0,n,u) result (y)
! 
! This function evauates the 1D maxwellian equilibrium distribution with given temperature and average velocity
! Temperature and average velocity can be arrays (corresponding to different points in x variable
!
! This function evaluates 
! f2_{M}(t,x,u)=(2\pi RT(t,x))^{-1/2} \exp(-\frac{(u-\bar{u})^2}{2RT}) 
!
! This is a reloadable function
!!!!!!!!!!!!!!!!!!!!!!!
!
! this is the copy when T is vector and u is scalar
function maxwelveldist_T_vector_u_vectors (T,u_0,v_0,w_0,n,u,v,w) result (y)
   real (DP), dimension (:), intent (in) :: T    ! temperature parameter (may depend on x)
   real (DP), dimension (:), intent (in) :: n    ! density parameter (may depend on x)
   real (DP), dimension (:), intent (in) :: u_0,v_0,w_0  ! average velocity (may depend on x)
   !!! T,n,u_0 must be the same size !!! 
   real (DP), dimension (:), intent (in) :: u,v,w    ! value of the velocity variable 
   real (DP), dimension (size(T),size(u))  :: y ! value of the dencity for this values of u and T   
!!!
   real (DP), dimension (size(T))  :: beta ! local variable  to keep temporary results    
   integer (I4B) :: i ! local counter 
!!!     
beta=sqrt(pi25DT*T)*(pi25DT*T)
do i=1,size(u)
y(:,i) = n*exp(-((u(i)-u_0)*(u(i)-u_0)+(v(i)-v_0)*(v(i)-v_0)+&
                  (w(i)-w_0)*(w(i)-w_0))/max(T,0.0000001_DP))/beta
end do 
end function maxwelveldist_T_vector_u_vectors

! this is the copy when T is vector and u is scalar
function maxwelveldist_T_vector (T,u_0,v_0,w_0,n,u,v,w) result (y)
   real (DP), dimension (:), intent (in) :: T    ! temperature parameter (may depend on x)
   real (DP), dimension (:), intent (in) :: n    ! density parameter (may depend on x)
   real (DP), dimension (:), intent (in) :: u_0, v_0,w_0  ! average velocity (may depend on x)
   !!! T,n,u_0 must be the same size !!! 
   real (DP), intent (in) :: u,v,w    ! value of the velocity variable 
   real (DP), dimension (size(T))  :: y ! value of the dencity for this values of u and T   
!!!
   real (DP), dimension (size(T))  :: beta ! local variable  to keep temporary results    
beta=sqrt(pi25DT*T)*(pi25DT*T)
y = n*exp(-((u-u_0)*(u-u_0)+(v-v_0)*(v-v_0)+&
                  (w-w_0)*(w-w_0))/max(T,0.0000001_DP))/beta
end function maxwelveldist_T_vector

! this is the copy when T is scalar and u is vector
function maxwelveldist_u_vectors (T,u_0,v_0,w_0,n,u,v,w) result (y)
   real (DP), intent (in) :: T    ! temperature parameter (scalar)
   real (DP), intent (in) :: u_0,v_0,w_0  ! average velocity (scalar)
   real (DP), intent (in) :: n    ! density parameter (scalar)
   real (DP), dimension (:), intent (in) :: u,v,w    ! values of the velocity variable 
   real (DP), dimension (size(u))  :: y ! value of the dencity for this values of u and T   
!!!
   real (DP)  :: beta ! local variable  to keep temporary results    
beta=sqrt(pi25DT*T)*(pi25DT*T)
y = n*exp(-((u-u_0)*(u-u_0)+(v-v_0)*(v-v_0)+&
                  (w-w_0)*(w-w_0))/max(T,0.0000001_DP))/beta
end function maxwelveldist_u_vectors

! this is the copy when both T and u are scalars
function maxwelveldist (T,u_0,v_0,w_0,n,u,v,w) result (y)
   real (DP), intent (in) :: T    ! temperature parameter (may depend on x)
   real (DP), intent (in) :: u_0,v_0,w_0  ! average velocity (scalar)
   real (DP), intent (in) :: n    ! density parameter (scalar)
   real (DP), intent (in) :: u,v,w    ! value of the velocity variable 
   real (DP)  :: y ! value of the dencity for this values of u and T   
!!!
   real (DP)  :: beta ! local variable  to keep temporary results    
beta=sqrt(pi25DT*T)*(pi25DT*T)
y = n*exp(-((u-u_0)*(u-u_0)+(v-v_0)*(v-v_0)+&
                  (w-w_0)*(w-w_0))/max(T,0.0000001_DP))/beta
end function maxwelveldist 

end module distributions_mod