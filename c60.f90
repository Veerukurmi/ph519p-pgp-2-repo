! Aim:To find out the phase shift for c60 fullerene without polarization
 !All the calculation is done in atomic units (a.u.)
! 1 Hartree= 27.2114 eV
! 2m/hbra^2=alpha=2.0 hbar=1, m=1, c=1 (in atomic units)

program phase_shift

implicit none

integer, parameter :: dp=selected_real_kind(14,200)
real(dp) , allocatable:: u(:), r(:),  f(:), Uo(:)
! u for radial wavefunction, r for radial distance, f has usual meaning in Numerov's method
integer:: i, mesh
real(dp), parameter:: alpha=2.0_dp, pi=3.141592653589793238_dp, Ht=27.2114_dp
!rmin and  is minimum and maximum radial distances respectively
real(dp) ::rmin, rmax, dr, energy, delx2
!lamx is the maximum value of angular momentum  quantum number 
integer :: lmax,l
!phaseshift is the phase shift
real(dp):: phaseshift,kk,tandelta
! lambda4 is 1/4 of the wavelength
!k is the wave vector
real(dp):: k,  lambda4, r2,  r1
! r1, r2 are the two random points beyons the potential
real(dp):: A,r0,del
! r0 and del are the inner radius and thickness of the well respectively
real(dp):: nl,jl
!nl and jl are the Bessel functions
integer:: i1,i2
!i1 and i2 are the corresponding mesh number to r1 and r2
real(dp):: rmax1,rmax2, xmin, xmax
integer :: j, n
!eup is maximum energy, elw is lower energy, and de is energy step size
real(dp):: eup, elw, de
 del=1.9_dp	
 r0=5.8_dp
 xmin=1.0_dp     
 xmax=r0+del+16.0_dp
 !mesh is the number of grid points
 mesh=3000
 !to read the value of l
 print*,"Enter the value of l :" 
 read(*,*) l
 eup=0.8_dp                   
 elw=0.05_dp                  
 de=0.001_dp  				 
 ! n is  the number of different energy values
!allocation of arrays
 allocate(r(0:mesh),u(0:mesh),f(0:mesh),Uo(0:mesh))
 
 n=nint((eup-elw)/de)+1    
! do loop for energy is changing in every interation and corresponding cross section value we are getting 	
do j =1, n
    
	energy=elw+(j-1)*de
    k=sqrt(energy*alpha)
	lambda4 =pi/k/2.0_dp    
	rmax1=xmax+lambda4       
	dr=(rmax1-xmin)/mesh
	r2 = rmax1
	i2 = mesh
	r1 = r2 - dr*nint(lambda4/dr)     !integer times dr
	i1 = i2 -nint(lambda4/dr)

	
	do i = 0,mesh
	    r(i) =xmin + float(i) * dr
              
	end do
		    
	!Numerov's method		

	delx2=dr*dr/12.0_dp

	!do loop for potential array for corresponding r values

	do i=0,mesh
			
			              
		if (r(i)<r0) then
			     Uo(i)=0.0_dp
	
		elseif (r(i)>r0+del) then 
			     Uo(i)=0.0_dp

			              		  
		else 
	            Uo(i)=-0.302_dp

   
			              		 
		endif
				f(i)=1.0_dp-((delx2)*(alpha))*((l*(l+1))/(alpha*r(i)*r(i))+Uo(i)-energy)
			
             
       end do
                     
            
            ! putting two initial values of wave function using asymptotic solution
             	
			u(0) =(r(0))**(l+1)     
			u(1) =(r(1))**(l+1)    
      
            !iterations for wavefunction by using Numerov's method
			do i =1,mesh-1
				u(i+1)=((12.0_dp-10.0_dp*f(i))*u(i)-f(i-1)*u(i-1))/f(i+1)
			end do
	       !conditions for  matching and and finding the value of kk=r1u(kr2)/r2u(kr1)and then phaseshift(phase shift)

			kk = (r1*u(i2))/(r2*u(i1))

			tandelta=(kk*jl(l,k,r1) - jl(l,k,r2) ) / (kk*nl(l, k,r1) - nl(l,k,r2) )
			
			!phase shifts for different values of l
			
			phaseshift = atan(tandelta)  
		
    !to get the data of phase shift at different energy values
              open(18, file='phaseshift.dat')
	             write(18,*) energy*Ht, phaseshift       

			
 
 end do
 
end program phase_shift


! to calculate the value of  Bessel functions for asymptotic solutions of the wavefunction
!not good  for large value of l

function jl (l,k,r)

  implicit none
  
  integer, parameter :: dp = selected_real_kind(14,200)
  integer :: l,lm
  real(dp) :: k,r, jl,jm1, jp1
  
  jm1= (cos (k*r)) / (k*r)
  jl = (sin (k*r)) / (k*r)
  
  do lm =0, l-1
    jp1= (2*lm+1)/(k*r)*jl-jm1
    jm1= jl
    jl = jp1
 end do

end function jl

function nl (l,k,r )

  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer :: l,lm
  real(dp) :: k,r, nl,nm1, np1

  nm1= (sin (k*r))/(k*r)
  nl =-(cos (k*r))/(k*r)
  
  do lm =0, l-1
    np1=(2*lm+1)/(k*r)*nl-nm1
    nm1=nl
     nl=np1
  end do
  
end function nl


