program modelreproc
implicit none
! find the kind of a high precision variable by finding the kind of 1.0d0
integer, parameter :: dp=kind(1.0d0)
! Define parameters
real (kind=dp), parameter :: h = 6.626068D0
! Planck’s constant (factors of 10 removed)
real (kind=dp), parameter :: k = 1.3806503D0
! Boltzmann constant (factors of 10 removed)
real (kind=dp), parameter :: c = 299792458.0D0
! Speed of light
real (kind=dp), parameter :: pi = 3.14159265D0
! Pi
real (kind=dp), parameter :: massp = 1.67262158D0
! Proton Mass (factors of 10 removed)
145

real (kind=dp), parameter :: Grav = 6.67300D0
! Gravitational constant (factors of 10 removed)
real (kind=dp), parameter :: StefanBoltz = 5.670400D0
! Stefan Boltzmann constant (factors of 10 removed)
real (kind=dp), parameter :: Thomson = 6.6524D0
! Thomson cross-section (factors of 10 removed)
real (kind=dp), parameter :: MSolar = 1.98892D0
! Solar Mass (factors of 10 removed)
real (kind=dp), parameter :: Parsec = 3.08568025D0
! Parsec (factors of 10 removed)
! Define variables
real (kind=dp), dimension(2) :: rin, eff
! Inner disc radius, efficiency
real (kind=dp), dimension(5) :: hx
! height of X-ray source (rin)
real (kind=dp), dimension(25) :: Mdot
! Accretion rate
real (kind=dp), dimension(30) :: Lx
! X-ray luminosity (W)
real (kind=dp), dimension(21) :: Massbh
! Blackhole mass
real (kind=dp), dimension(6001) :: radius, Temp, Flux, lagr
! Radius, temperature, flux, lagradius i.e. inc=0
real (kind=dp) :: w2cf, m2cf, w1cf, ucf, bcf, vcf
! Cumulative fluxes(w2,m2,w1,u,b,v)
real (kind=dp) :: w2flux, m2flux, w1flux, uflux, bflux, vflux
! Fluxes(w2,m2,w1,u,b,v)
real (kind=dp), dimension(10000) :: wavelength, prefactor
 
 
! Wavelength(lambda), For blackbody calculations
real (kind=dp), dimension(10000) :: Fluxlambda, Flambda
! Flux(lambda)
real (kind=dp), dimension(10000) :: TotFluxlambda
! TotalFlux(lambda)
real (kind=dp), dimension(10000,81) :: TFlambda
! TotalFlux(lambda,lag)
real (kind=dp) :: expnt, expntdenom
! For blackbody calculations
real (kind=dp), dimension(641) :: wavelengtha, wavelengthb
real (kind=dp), dimension(7,641) :: effarea
! Wavelength (1,n) and effarea for w2, m2, w1, u, b and v (2,7)
real (kind=dp) :: constantsa, constantsb, constantsc, distfactor
! Constants in hc/lambda k t, reproc constants,
disc flux constants, correction for distance
real (kind=dp) :: factora, factorb, factorc
! Factors in front of disc, reproc terms
real (kind=dp), dimension(81) :: lag
! Lag in seconds
real (kind=dp) :: laga, lagb
! Lag for use in integration
real (kind=dp), dimension(6) :: inclination, cosinc, sininc
! Inclination, cos(inclination), sin(inclination)
real (kind=dp) :: dhr, dlr
real (kind=dp) :: areaterms, deltar, distconstants
real (kind=dp) :: Router
! Outer radius ie when temp < 1750K

 
real (kind=dp), dimension(6) :: rad, temper
! Radius, temperature within lag calc
! Filenames
character(len=100) :: fnames, fnamet
! Spectrum and temperature profile
character(len=100) :: fnamer, fname
! Lag and cumulative flux
integer :: fileno, fileno2
! Filenos for calling files
integer :: mbh, md, hxray, l, ri, r, i, j, lambda, lagi, inc, rstep
! Counting integers
! Read in Effective Areas
open(21, file=’uw2effarea.dat’)
open(22, file=’um2effarea.dat’)
open(23, file=’uw1effarea.dat’)
open(24, file=’uuueffarea.dat’)
open(25, file=’ubbeffarea.dat’)
open(26, file=’uvveffarea.dat’)
open(27, file=’modelreproclag.log’)
write(*,*) ’Testing...’
wavelengtha=0.0D0
wavelengthb=0.0D0
effarea=0.0D0
do j=1,6
do i=1,641
read(j+20,*,end=100) wavelengtha(i), wavelengthb(i), effarea(j+1,i)
end do
100 close(j+20)
write(*,*) j
end do
! Put effective area in metres^2 and wavelengths in 10^-7 m
do i=1,641
 
 
effarea(1,i)=(wavelengtha(i)+wavelengthb(i))/2000.0D0
effarea(2,i)=effarea(2,i)/10000.0D0
effarea(3,i)=effarea(3,i)/10000.0D0
effarea(4,i)=effarea(4,i)/10000.0D0
effarea(5,i)=effarea(5,i)/10000.0D0
effarea(6,i)=effarea(6,i)/10000.0D0
effarea(7,i)=effarea(7,i)/10000.0D0
end do
write(*,*) effarea(1,300), effarea(2,300), effarea(3,300),
effarea(4,300), effarea(5,300), effarea(6,300), effarea(7,300)
! Calculate wavelengths in units of 10^-7 m
do lambda=1,10000
wavelength(lambda)=lambda*1.0D0 ! Angstroms
wavelength(lambda)=wavelength(lambda)/1000.0D0 ! 10^-7 m
prefactor(lambda) = (2.0D0*c) / (wavelength(lambda)**4.0D0)
! Prefactor for Planck formula
end do
write(*,*) wavelength(10), prefactor(10)
! Define input values
! Black hole masses
do mbh=1,11
Massbh(mbh)=mbh-1.0D0 ! Range 0-10
Massbh(mbh)=Massbh(mbh)/2.0D0 ! 0-5
Massbh(mbh)=10000.0D0*10.0D0**Massbh(mbh) !10^4 -10^9
write(*,*) Massbh(mbh)
end do
! Accretion rate/Eddington
do md=1,13
Mdot(md)=md-3.0D0
Mdot(md)=-0.5D0*Mdot(md) ! 1 - -5
Mdot(md)=10.0D0**Mdot(md) ! 10^1 - 10^-5
write(*,*) Mdot(md)
end do
! X-ray source height in units of inner radius
do hxray=1,5
hx(hxray)=hxray-1.0D0 ! 0-4
hx(hxray)=hx(hxray)/2.0D0 ! 0-2

 
hx(hxray)=10.0D0**hx(hxray) ! 1-10^2
write(*,*) hx(hxray)
end do
! X-ray source luminosity (J)
Lx(1)=0.0D0
do l=2,16
Lx(l)=l-2.0D0 ! 0-28
Lx(l)=Lx(l)*0.5D0 ! 0-7
Lx(l)=Lx(l)+32.0D0 ! 10^32 - 10^39
! 0.7 takes into account albedo of 0.3
Lx(l)=0.7D0*10.0D0**Lx(l)
write(*,*) Lx(l)
end do
! Inner disc radius and efficiency for Kerr(1) and Schwarzschild(2)
rin(1)=1.235D0 ! Kerr inner disc radius
rin(2)=6.0D0 ! Schwarzschild inner disc radius
eff(1)=0.32D0 ! Accretion efficiency
eff(2)=0.06D0 ! Accretion efficiency
! Calculate lags for evaluation
do lagi=1,81
lag(lagi)=lagi-21.0D0 ! -20 - 60
lag(lagi)=lag(lagi)*0.1D0 ! -2 - 6
lag(lagi)=10.0D0**lag(lagi)
! 0.01 - 1,000,000 log scale
end do
! Inclination angles (rads)
do inc=1,6
inclination(inc)=inc-1.0D0 ! 0-5
inclination(inc)=inclination(inc)*pi/12.0D0
! 0 - 75 degrees in radians
cosinc(inc)=cos(inclination(inc))
sininc(inc)=sin(inclination(inc))
end do
! hc/(lambda k T) constant parts
constantsa=2.0D0*Msolar*Grav
! 2 G Msolar / 10^19
 
 
constantsa=constantsa**2.0D0
! (2 G Msolar)^2 / 10^38
constantsa=constantsa*StefanBoltz
! StefanBoltz (2 G Msolar)^2 / 10^30
constantsa=constantsa**0.25D0
! StefanBoltz^0.25 sqrt(2 G Msolar) / 10^7.5
constantsa=h*constantsa/k
! (h/k) StefanBoltz^0.25 sqrt(2 G Msolar) * 10^3.5
constantsa=constantsa*10.0D0**-3.5D0
! (h/k) StefanBoltz^0.25 sqrt(2 G Msolar)
! Constant parts of the reprocessed flux
constantsb=1.0D0/pi
! Constant parts of the intrinsic disc flux
constantsc=6.0D0*massp*c*Grav*Msolar/Thomson
! (6 mp c Grav Msolar / Thomson) * 10^-21
constantsc=constantsc*10.0D0**21.0D0
! (6 mp c Grav Msolar / Thomson)
! Distance correction factor constants
distconstants=Grav*Msolar/c**2.0D0
! Rg for Sun / 10^19
distconstants=distconstants/(Parsec*10.0D0**3.0D0)
! Rg/Mpc
distconstants=distconstants**2.0D0
! (Rg/Mpc)^2
! This 4pi comes from the surface of a sphere,
doesn’t cancel with a pi from the disc because
only using the scale size here
distconstants=distconstants/(4.0D0*pi)
! Rg(solar)^2/(4 pi Mpc^2)
write(*,*) constantsa, constantsb, constantsc, distconstants
do r=1,6001
radius(r)=r*1.0D0 ! 1-6001
radius(r)=radius(r)-1.0D0
radius(r)=radius(r)*0.001D0 ! 0.001-6
radius(r)=10.0D0**radius(r) ! 10^0.001-10^6
end do

 
do ri=1,1
do mbh=1,11
! Overall factors in hc/kT
and Mass/inner radius distance correction factor
factora=Massbh(mbh)*Rin(ri)
! M Rin
distfactor=factora**2.0D0
! (M Rin)^2 for distance factor ie Rg(sun) --> Rin(mbh)
factora=factora**0.5D0
! sqrt(M Rin) Mass and inner radius correction factor
do md=1,13
! Factor in intrinsic flux
factorc=mdot(md)*Massbh(mbh)/(eff(ri)*Rin(ri))
! accretion rate/efficiency correction factor
do hxray=1,5
do l=13,16
inc=1
j=350+mbh*50
! Define filenames
if (md<3.5) then
i=50*(3-md)
if (ri==1) then
write(fnamet,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Kerr4/Temps’,j,’+’,i,’H’,hxray,’L’,l,’.dat’
write(fnames,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Kerr4/Spec’,j,’+’,i,’H’,hxray,’L’,l,’.dat’
write(fname,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Kerr4/Total’,j,’+’,i,’H’,hxray,’L’,l,’.dat’
write(fnamer,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Kerr4/Lag’,j,’+’,i,’H’,hxray,’L’,l,’.dat’
end if
if (ri==2) then
write(fnamet,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Schwarzschild4/Temps’,j,’+’,i,’H’,hxray,’L’,l,’.dat’
write(fnames,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Schwarzschild4/Spec’,j,’+’,i,’H’,hxray,’L’,l,’.dat’
write(fname,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Schwarzschild4/Total’,j,’+’,i,’H’,hxray,’L’,l,’.dat’
write(fnamer,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Schwarzschild4/Lag’,j,’+’,i,’H’,hxray,’L’,l,’.dat’
 
 
end if
else
i=50*(md-3)
if (ri==1) then
write(fnamet,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Kerr4/Temps’,j,’-’,i,’H’,hxray,’L’,l,’.dat’
write(fnames,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Kerr4/Spec’,j,’-’,i,’H’,hxray,’L’,l,’.dat’
write(fname,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Kerr4/Total’,j,’-’,i,’H’,hxray,’L’,l,’.dat’
write(fnamer,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Kerr4/Lag’,j,’-’,i,’H’,hxray,’L’,l,’.dat’
end if
if (ri==2) then
write(fnamet,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Schwarzschild4/Temps’,j,’-’,i,’H’,hxray,’L’,l,’.dat’
write(fnames,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Schwarzschild4/Spec’,j,’-’,i,’H’,hxray,’L’,l,’.dat’
write(fname,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Schwarzschild4/Total’,j,’-’,i,’H’,hxray,’L’,l,’.dat’
write(fnamer,’(a,i0,a,i3.3,a,i0,a,i0,a)’)
’Schwarzschild4/Lag’,j,’-’,i,’H’,hxray,’L’,l,’.dat’
end if
end if
! Factor in reprocessed flux
factorb=hx(hxray)*Lx(l)
! Set variables to zero
write(*,*) fnamet
write(*,*) fnames
write(*,*) fname
open(30, file=fnamet)
open(31, file=fnames)
open(32, file=fname)
write(*,*)fnamer
open(33, file=fnamer)
! Reset TotFlux to zero
do lambda=1,10000
TotFluxlambda(lambda)=0.0D0
end do
do r=1,6001
Temp(r)=((1-(1/radius(r))**0.5D0)*factorc*constantsc)

 
/(radius(r))**3.0D0
! Part of temperature due to accretion in disc
Temp(r)=Temp(r)+(factorb*constantsb)
/((hx(hxray)**2.0D0+radius(r)**2.0D0)**1.5D0)
! + temperature due to reprocessing
Temp(r)=Temp(r)**0.25D0
! ^0.25
Temp(r)=Temp(r)/(factora*constantsa)
! kT/hc
Temp(r)=c*Temp(r)*h/(k*(10.0D0**11.0D0))
! T
write(30,*) radius(r)*Rin(ri), Temp(r)
if (Temp(r) < 1750.0D0) then
if (Temp(r-1) > 1750.0D0) then
Router = Radius(r)
! Outer extent of disc due to dust sublimation
write(*,*) Router
end if
else
if (r > 1.5) then
do lambda=1,10000
Fluxlambda(lambda)=exp(h*c/((k*(Temp(r)+Temp(r-1))/2.0D0)
*wavelength(lambda)*10.0D0**4.0D0))
! 10^4 from 10^11 from h/k and 10!-7 wavelength
Fluxlambda(lambda)=pi*(radius(r)**2.0D0-radius(r-1)**2.0D0)
/(Fluxlambda(lambda)-1.0D0)
! Area within annulus
Fluxlambda(lambda)=Fluxlambda(lambda)
*prefactor(lambda)*10.0D0**28.0D0
! Prefactor and wavelength correction (10^7)^4
Fluxlambda(lambda)=Fluxlambda(lambda)*distfactor*distconstants
! Correction for distance
TotFluxlambda(lambda)=TotFluxlambda(lambda)+Fluxlambda(lambda)
! Cumulative
end do
w2cf=0.0D0
m2cf=0.0D0
w1cf=0.0D0
ucf=0.0D0
bcf=0.0D0
vcf=0.0D0
 
 
w2flux=0.0D0
m2flux=0.0D0
w1flux=0.0D0
uflux=0.0D0
bflux=0.0D0
vflux=0.0D0
do i=1,641
do lambda=1585+10*i,1595+10*i
if (lambda == 1585+10*i .or. lambda == 1595+10*i) then
w2cf=w2cf+TotFluxlambda(lambda)*effarea(2,i)*10.0**-10.0D0/2.0D0
! x effarea and correct for counts
per angstrom instead of metres
m2cf=m2cf+TotFluxlambda(lambda)*effarea(3,i)*10.0**-10.0D0/2.0D0
! x effarea and correct for counts
per angstrom instead of metres
w1cf=w1cf+TotFluxlambda(lambda)*effarea(4,i)*10.0**-10.0D0/2.0D0
! x effarea and correct for counts
per angstrom instead of metres
ucf=ucf+TotFluxlambda(lambda)*effarea(5,i)*10.0**-10.0D0/2.0D0
! x effarea and correct for counts
per angstrom instead of metres
bcf=bcf+TotFluxlambda(lambda)*effarea(6,i)*10.0**-10.0D0/2.0D0
! x effarea and correct for counts
per angstrom instead of metres
vcf=vcf+TotFluxlambda(lambda)*effarea(7,i)*10.0**-10.0D0/2.0D0
! x effarea and correct for counts
per angstrom instead of metres
w2flux=w2flux+Fluxlambda(lambda)*effarea(2,i)*10.0**-10.0D0/2.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
m2flux=m2flux+Fluxlambda(lambda)*effarea(3,i)*10.0**-10.0D0/2.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
w1flux=w1flux+Fluxlambda(lambda)*effarea(4,i)*10.0**-10.0D0/2.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
uflux=uflux+Fluxlambda(lambda)*effarea(5,i)*10.0**-10.0D0/2.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
bflux=bflux+Fluxlambda(lambda)*effarea(6,i)*10.0**-10.0D0/2.0D0
! flux x effarea and correct for counts

 
per angstrom instead of metres
vflux=vflux+Fluxlambda(lambda)*effarea(7,i)*10.0**-10.0D0/2.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
else
w2cf=w2cf+TotFluxlambda(lambda)*effarea(2,i)*10.0**-10.0D0
! x effarea and correct for counts
per angstrom instead of metres
m2cf=m2cf+TotFluxlambda(lambda)*effarea(3,i)*10.0**-10.0D0
! x effarea and correct for counts
per angstrom instead of metres
w1cf=w1cf+TotFluxlambda(lambda)*effarea(4,i)*10.0**-10.0D0
! x effarea and correct for counts
per angstrom instead of metres
ucf=ucf+TotFluxlambda(lambda)*effarea(5,i)*10.0**-10.0D0
! x effarea and correct for counts
per angstrom instead of metres
bcf=bcf+TotFluxlambda(lambda)*effarea(6,i)*10.0**-10.0D0
! x effarea and correct for counts
per angstrom instead of metres
vcf=vcf+TotFluxlambda(lambda)*effarea(7,i)*10.0**-10.0D0
! x effarea and correct for counts
per angstrom instead of metres
w2flux=w2flux+Fluxlambda(lambda)*effarea(2,i)*10.0**-10.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
m2flux=m2flux+Fluxlambda(lambda)*effarea(3,i)*10.0**-10.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
w1flux=w1flux+Fluxlambda(lambda)*effarea(4,i)*10.0**-10.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
uflux=uflux+Fluxlambda(lambda)*effarea(5,i)*10.0**-10.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
bflux=bflux+Fluxlambda(lambda)*effarea(6,i)*10.0**-10.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
vflux=vflux+Fluxlambda(lambda)*effarea(7,i)*10.0**-10.0D0
! flux x effarea and correct for counts
per angstrom instead of metres
 
 
end if
end do
end do
write(32,*) radius(r)*Rin(ri), w2cf, m2cf, w1cf, ucf, bcf, vcf
end if
 r(r)=sqrt(hx(hxray)**2.0D0+radius(r)**2.0D0)+hx(hxray)
lagr(r)=Rin(ri)*lagr(r)*Grav*Msolar*Massbh(mbh)*10.0D0**19.0D0/(c**3.0D0)
if (r > 1.5) then
w2flux=w2flux/(lagr(r)-lagr(r-1))
m2flux=m2flux/(lagr(r)-lagr(r-1))
w1flux=w1flux/(lagr(r)-lagr(r-1))
uflux=uflux/(lagr(r)-lagr(r-1))
bflux=bflux/(lagr(r)-lagr(r-1))
vflux=vflux/(lagr(r)-lagr(r-1))
write(33,*) inc, lagr(r), w2flux, m2flux, w1flux, uflux, bflux, vflux
end if
end if
end do
do lambda=1,10000
write(31,*) lambda,TotFluxlambda(lambda)
end do
close(30)
close(31)
close(32)
do inc=2,6
do lambda=1,10000
do lagi=1,81
TFlambda(lambda,lagi)=0.0D0
end do
end do
write(*,*) inc
do r=2,6001
if (radius(r)<router) then
do rstep=1,6
rad(rstep)=rstep-1.0D0
rad(rstep)=rad(rstep)/5.0D0
temper(rstep)=rad(rstep)*(temp(r)-temp(r-1))
temper(rstep)=temper(rstep)+temp(r-1)
rad(rstep)=rad(rstep)*(radius(r)-radius(r-1))
rad(rstep)=rad(rstep)+radius(r-1)
deltar=(radius(r)-radius(r-1))/5.0D0

 
do lambda=1,10000
Fluxlambda(lambda)=0.0D0
end do
do lambda=
5,8005
Fluxlambda(lambda)
=exp(h*c/(k*temper(rstep)*wavelength(lambda)*10.0D0**4.0D0))-1.0D0
Fluxlambda(lambda)=prefactor(lambda)*10.0D0**28.0D0
/Fluxlambda(lambda)
if (rstep==1) then
Fluxlambda(lambda)=Fluxlambda(lambda)/2.0D0
end if
if (rstep==6) then
Fluxlambda(lambda)=Fluxlambda(lambda)/2.0D0
end if
end do
do lagi=1,81
if (lagi==1) then
laga=0.0D0
else
laga=lag(lagi-1)
end if
lagb=lag(lagi)
laga=laga/rin(ri) ! Lag lower/inner radius
lagb=lagb/rin(ri) ! Lag upper/inner radius
laga=laga/Grav ! Lagl/(Rin G) / 10^11
lagb=lagb/Grav ! Lagu/(Rin G) / 10^11
laga=laga/Msolar ! Lagl/(Msun Rin G) * 10^19
lagb=lagb/Msolar ! Lagu/(Msun Rin G) * 10^19
laga=laga/Massbh(mbh) ! Lagl/(Mbh Rin G) * 10^19
lagb=lagb/Massbh(mbh) ! Lagu/(Mbh Rin G) * 10^19
! Lag in terms of distance to inner disc radius
laga=(laga*c**3.0D0)/10.0D0**19.0D0 ! Lagl*c/(Rin Rg(mbh))
lagb=(lagb*c**3.0D0)/10.0D0**19.0D0 ! Lagl*c/(Rin Rg(mbh))
dhr=sqrt(hx(hxray)**2.0D0+rad(rstep)**2.0D0)+hx(hxray)*cosinc(inc)
dlr=dhr-rad(rstep)*sininc(inc)
dhr=dhr+rad(rstep)*sininc(inc)
if (lagb>dlr) then
if (laga<dhr) then
if (lagb >= dhr) then
if (laga <= dlr) then
areaterms=2.0D0*pi*rad(rstep)
do lambda=1595,8005
Flambda(lambda)=Fluxlambda(lambda)*areaterms*deltar
TFlambda(lambda,lagi)=TFlambda(lambda,lagi)+Flambda(lambda)
end do
else
areaterms=(pi/2.0D0)-atan(sqrt((laga-dlr)/(dhr-laga)))
areaterms=areaterms*4.0D0*rad(rstep)
do lambda=1595,8005
Flambda(lambda)=Fluxlambda(lambda)*areaterms*deltar
TFlambda(lambda,lagi)=TFlambda(lambda,lagi)+Flambda(lambda)
end do
end if
else
if (laga <= dlr) then
areaterms=atan((sqrt((lagb-dlr)/(dhr-lagb))))
areaterms=areaterms*4.0D0*rad(rstep)
do lambda=1595,8005
Flambda(lambda)=Fluxlambda(lambda)*areaterms*deltar
TFlambda(lambda,lagi)=TFlambda(lambda,lagi)+Flambda(lambda)
end do
else
areaterms=atan((sqrt((lagb-dlr)/(dhr-lagb))))
-atan(sqrt((laga-dlr)/(dhr-laga)))
areaterms=areaterms*4.0D0*rad(rstep)
do lambda=1595,8005
Flambda(lambda)=Fluxlambda(lambda)*areaterms*deltar
TFlambda(lambda,lagi)=TFlambda(lambda,lagi)+Flambda(lambda)
end do
end if
end if
end if
end if
end do
end do
end if
end do
do lagi=1,81
w2flux=0.0D0
m2flux=0.0D0
w1flux=0.0D0

 
uflux=0.0D0
bflux=0.0D0
vflux=0.0D0
do i=1,641
do lambda=1585+10*i,1595+10*i
if (lambda == 1585+10*i .or. lambda == 1595+10*i) then
w2flux=w2flux+TFlambda(lambda,lagi)
*effarea(2,i)*10.0**-10.0D0/2.0D0
m2flux=m2flux+TFlambda(lambda,lagi)
*effarea(3,i)*10.0**-10.0D0/2.0D0
w1flux=w1flux+TFlambda(lambda,lagi)
*effarea(4,i)*10.0**-10.0D0/2.0D0
uflux=uflux+TFlambda(lambda,lagi)
*effarea(5,i)*10.0**-10.0D0/2.0D0
bflux=bflux+TFlambda(lambda,lagi)
*effarea(6,i)*10.0**-10.0D0/2.0D0
vflux=vflux+TFlambda(lambda,lagi)
*effarea(7,i)*10.0**-10.0D0/2.0D0
else
w2flux=w2flux+TFlambda(lambda,lagi)
*effarea(2,i)*10.0**-10.0D0
m2flux=m2flux+TFlambda(lambda,lagi)
*effarea(3,i)*10.0**-10.0D0
w1flux=w1flux+TFlambda(lambda,lagi)
*effarea(4,i)*10.0**-10.0D0
uflux=uflux+TFlambda(lambda,lagi)
*effarea(5,i)*10.0**-10.0D0
bflux=bflux+TFlambda(lambda,lagi)
*effarea(6,i)*10.0**-10.0D0
vflux=vflux+TFlambda(lambda,lagi)
*effarea(7,i)*10.0**-10.0D0
end if
end do
end do
if (lagi==1) then
w2flux=cosinc(inc)*w2flux*distfactor*distconstants/lag(lagi)
m2flux=cosinc(inc)*m2flux*distfactor*distconstants/lag(lagi)
w1flux=cosinc(inc)*w1flux*distfactor*distconstants/lag(lagi)
uflux=cosinc(inc)*uflux*distfactor*distconstants/lag(lagi)
bflux=cosinc(inc)*bflux*distfactor*distconstants/lag(lagi)
vflux=cosinc(inc)*vflux*distfactor*distconstants/lag(lagi) 
else
w2flux=cosinc(inc)*w2flux*distfactor*distconstants
/(lag(lagi)-lag(lagi-1))
m2flux=cosinc(inc)*m2flux*distfactor*distconstants
/(lag(lagi)-lag(lagi-1))
w1flux=cosinc(inc)*w1flux*distfactor*distconstants
/(lag(lagi)-lag(lagi-1))
uflux=cosinc(inc)*uflux*distfactor*distconstants
/(lag(lagi)-lag(lagi-1))
bflux=cosinc(inc)*bflux*distfactor*distconstants
/(lag(lagi)-lag(lagi-1))
vflux=cosinc(inc)*vflux*distfactor*distconstants
/(lag(lagi)-lag(lagi-1))
end if
write(33,*) inc, lag(lagi), w2flux,
m2flux, w1flux, uflux, bflux, vflux
end do
end do
close(33)
end do
end do
end do
end do
end do
end program modelreproc