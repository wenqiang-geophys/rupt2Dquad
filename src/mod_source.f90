module mod_source

implicit none

contains

function svf_ricker(t,fc,tdelay)

real*8 :: svf_ricker
real*8 :: f0,r,rr,t,fc,tdelay
real*8,parameter :: pi = 3.141593

f0 = sqrt(pi)/2.
r = pi * fc * (t-tdelay)
rr = r**2

svf_ricker = r*(3.-2.*rr)*exp(-rr)*f0*pi*fc

end function


end module
