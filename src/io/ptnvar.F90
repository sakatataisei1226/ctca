!
  subroutine ptnvar(pmax, pomega, pmax_tmp, t, t_max)
!
!   ____________________________________________________________
!
!               S U B R O U T I N E   P T N V A L
!   ____________________________________________________________
!
!-------------------- 
  implicit none
  real(kind=8),intent(in) :: pmax, t, t_max, pomega
  real(kind=8),intent(out) :: pmax_tmp 


!------------------- step function
      if(t_max.eq.0.0) then
        pmax_tmp = pmax 
        return
      end if 


!------------------- linear increase of potential(to ftime)
      if(pomega.eq.0.0) then
        if(t.lt.t_max) then
          pmax_tmp = pmax/t_max * t
        else
          pmax_tmp = pmax
        end if
      end if
      print*, 'pmax_tmp  =  ', pmax_tmp


  return
  end
