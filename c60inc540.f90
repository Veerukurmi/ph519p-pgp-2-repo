! Aim:To find out the sigma section  for sthe potential(C60@C540)
! for calculated parameters C60@C540 , Uo1=-0.302 Uo2=-0.404
 !All the calculationn is done in atomic units (a.u.)
! 1 Hartree= 27.2114 eV
! 2m/hbra^2=alpha=2.0
! hbar=1, m=1, c=1 (in atomic units)
!If we need the time delay from this code, then we can just convert the phase shift, energy values, in array form
program phase_shift

implicit none

integer, parameter :: dp=selected_real_kind(14,200)
real(dp) , allocatable:: u(:), r(:),  f(:), Uo(:)
integer:: i,j, mesh,ii
real(dp), parameter:: alpha=2.0_dp, pi=3.141592653589793238_dp, Ht=27.2114_dp
!rmin and rmax  are the  minimum and the maximum radial distances respectively
real(dp) ::rmin, rmax, dr, energy, delx2
!lamx is the maximum value of angular momentum quantum number            
integer :: l
!deltal is the phase shift
real(dp):: phaseshift,kk,tandelta,sigma
! quarterlambda is 1/4 of the wavelength
!k is the wave vector
real(dp):: k,  lambda4, r2,  r1
! r1, r2 are the two random points beyond the potential
real(dp):: r01,r02,del1,del2
!r01, r02, del1, del2 are the inner radius of c60 and c240 and thicknesses of the wells
real(dp):: nl,jl
! nl and jl are Bessel functions
integer:: i1,i2
!i1 and i2 are the corresponding mesh number to r1 and r2
real(dp):: rmax1,rmax2, xmin, xmax
integer ::  n
!eup is maximum energy, elw is lower energy, and de is energy step size
real(dp):: eup, elw, de
  character (len=80) :: fileout,fileout1
 del1=1.9_dp	
 del2=1.9
 r01=5.8_dp
 r02=18.85_dp
 xmin=1.0_dp     
 xmax=r02+del2+16.0_dp
 !here we are taking a rough value 16 because we want to take any two points that are beyond the potential range. So, we can taking any number like 2,3,5,etc.
 dr=0.01_dp
 !no. of grid points
 mesh=3000
 !to read the value of angular momentum number
 print*,"Enter the value of l :" 
 read(*,*) l
 eup=0.8_dp                   
 elw=0.05_dp                  
 de=0.001_dp  				   
 n=nint((eup-elw)/de)+1        
 ! this n is  the number of different energy values
 
 
 
!allocation of arrays	
 allocate(r(0:mesh),u(0:mesh),f(0:mesh),Uo(0:mesh))
 
 
 !for output files name 
print*, "enter the name of file for phaseshift:"
read (*,*) fileout
!print*, "enter the name of file for cross section:"
!read (*,*) fileout1
   
 ! do loop fo energies	
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
   
            
	! Numerov's method
            
	delx2=dr*dr/12.0_dp
            
	do i=0, mesh
			
			              
			        if (r(i)>r01.and.r(i)<(r01+del1)) then
			              		   Uo(i)=-0.302_dp
			              		   	  	
			             			elseif(r(i)>r02.and.r(i)<(r02+del2))then
			              		   		Uo(i)=-0.404_dp
				       	else 
	              		 	Uo(i)=0.0_dp

   
			              		 
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

			!conditions for  matching and and finding the value of kk=r1u(kr2)/r2u(kr1)and then deltal(phase shift)

			kk = (r1*u(i2))/(r2*u(i1))			
			tandelta=(kk*jl(l,k,r1) - jl(l,k,r2) ) / (kk*nl(l, k,r1) - nl(l,k,r2) )
			
			!to calculate the phase shift for different l values 
			
			phaseshift = atan(tandelta)  
		   
		    !to get the data of phase shift at different energy values
 	    	  open (18,file=fileout, status='unknown', form='formatted')
	             write(18,*) energy*Ht, phaseshift      
	             
	             !for cross section
	            ! sigma = ((4*pi)/(k**2))*(2*l+1)*sin(phaseshift)**2   
	               
	             !for partial cross section
	              ! open (19,file=fileout1, status='unknown', form='formatted')
	             !write(19,*) energy*Ht, sigma
   
          
	
	
 enddo
 
!to extract the potential

!open(12, file="potential.dat")		
		!do j=0,mesh
        
		!	   write (12,*) r(j), Uo(j)
		!		end do	 

 
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


