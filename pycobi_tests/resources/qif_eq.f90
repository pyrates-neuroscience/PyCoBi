module qif_equations

double precision :: PI = 4.0*atan(1.0)

contains


subroutine qif(t,y,dy,tau,Delta,I_ext,eta,weight)

implicit none

double precision, intent(in) :: t
double precision, intent(in) :: y(2)
double precision :: r
double precision :: v
double precision :: r_in
double precision, intent(inout) :: dy(2)
double precision, intent(in) :: tau
double precision, intent(in) :: Delta
double precision, intent(in) :: I_ext
double precision, intent(in) :: eta
double precision, intent(in) :: weight


r = y(1)
v = y(2)
r_in = r*weight

dy(1) = (Delta/(pi*tau) + 2.0*r*v)/tau
dy(2) = (I_ext + eta - pi**2*r**2*tau**2 + r_in*tau + v**2)/tau

end subroutine


end module


subroutine func(ndim,y,icp,args,ijac,dy,dfdu,dfdp)

use qif_equations
implicit none
integer, intent(in) :: ndim, icp(*), ijac
double precision, intent(in) :: y(ndim), args(*)
double precision, intent(out) :: dy(ndim)
double precision, intent(inout) :: dfdu(ndim,ndim), dfdp(ndim,*)

call qif(args(14), y, dy, args(1), args(2), args(3), args(4), args(5))

end subroutine func


subroutine stpnt(ndim, y, args, t)

implicit None
integer, intent(in) :: ndim
double precision, intent(inout) :: y(ndim), args(*)
double precision, intent(in) :: t

args(1) = 1.0  ! tau
args(2) = 1.0  ! Delta
args(3) = 0.0  ! I_ext
args(4) = -5.0  ! eta
args(5) = 15.0  ! weight
y(1) = 0.01  ! r
y(2) = -2.0  ! v

end subroutine stpnt



subroutine bcnd
end subroutine bcnd


subroutine icnd
end subroutine icnd


subroutine fopt
end subroutine fopt


subroutine pvls
end subroutine pvls
