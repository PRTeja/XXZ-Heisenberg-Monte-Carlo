PROGRAM hei
implicit none
integer*8::lattice, num, nnum, ii,jj,kk, g1
integer*8, parameter::ndiv = 30, MCcycles =75*10**2, Eqcycles = 75*10**4, ndiv2 = 60, nbins = 200, errorbins = 300
real*8,ALLOCATABLE, DIMENSION(:,:,:):: spin
real*8,parameter::J1 = 1.D0, zeronum = 0.d0, pi = acos(-1.D0)!, J2  = 1.D0!, delta = 1.D0
real*8:: delta, J2, AFMbinder4, amag, Pmax, in_x, in_y, out_z, inplane, magerr, time0, time1
integer*8::N,i,j,k, i0, counter, i22, nall, N2, N3,bins,i0t, int1 = MCcycles/errorbins, i33
real*8::t0,t1, T, energy, beta, chiral, ch0, energy0, magnav, stripeav, Enav, En2av, chix, chiy, chiz,  suscchir, Cvav
real*8:: chi1, chi2, chi3, rands, stripe, binder, stagmag, stagav, stagchi, den, magz,ds, chirav, AFMav
character(len=120)::str, str2, title, title2, title10, title11, title12,  title13, title14, title19, title15, title16
integer*8,dimension(nbins)::  P_En
real*8,dimension(MCcycles)::magn_data, stripe_data, chiral_data, En_data, stagmag_data,En2_data,AFMdata,magxdata,magydata,magzdata
real*8,DIMENSION(3)::magn,magnetisation, field, stmag
real*8, dimension(MCcycles/errorbins):: magn_temp, stripe_temp, chiral_temp, En_temp, stagmag_temp
real*8:: mag_err, stripe_err, chiral_err, En_err, stagmag_err, total_err1,total_err2,total_err3, chi_in, chi_out, Cv_err
real*8, dimension(errorbins):: mag_avs, stripe_avs, chiral_avs, En_avs, stagmag_avs, Cv_avs, chi1_avs, chi2_avs, chi3_avs, En2_avs
real*8, dimension(errorbins):: inchi_avs, outchi_avs, SF_avs, AFMbind1, AFMbind2, Enbind, AFMbind0
real*8, dimension(MCcycles):: Gtdata, magqtdata, smagn2data
real*8:: Gt, magqt, Gtav, magqtav, smagn2av,xit,xis,xitcon,xiscon,Gtcon,AFMbinderr1,AFMbinderr2,Enbinderr, AFMbinderr0

!variables
real*8:: SF0, SF1, SFav, SFerr, kval(3) 
real*8:: magn0, magn1, magn2, magn4, stripe0, stripe1, stripe2, stripe4, chir0, chir1, chir2, chir4, En1, En2, En4, Gt1, magqt1
real*8:: Cv, AFMBinder1, AFMBinder2, AFMBinder3, EnBinder, E1sum, E2sum, E1av, E2av,stiff,chi1_err,chi2_err,chi3_err,inchi_err
real*8:: AFMchi1, AFMchi2, AFMchi3, magnx1, magnx2, magny1, magny2, magnz1, magnz2,in_chi,out_chi,inchiav,outchiav,outchi_err
nall =120
k=10
nnum = MCcycles/errorbins
!$OMP PARALLEL DO default(firstprivate)
do i22=31,31

write(str,*)i22
lattice=120
!read(*,*)lattice

call cpu_time(time0)
!kval = 2*pi/lattice
kval(:)  = pi
den = float(lattice*lattice)							 ! L^d
ALLOCATE(spin(lattice,lattice,3))					 

J2 = 0.7d0
delta = 0.d0
do i33 = 1,1

t0 = 0.420d0
t1 = 0.480d0 

title = 'details'//str!2//'_'//str
title2 =''
do i=1,len(title)   
   if (title(i:i).ne.' ')title2=trim(title2)//trim(title(i:i))   
end do
title10 =  (title2)
	Open(2, file =title10)
	write(2,*) 'T_range ', t0, t1, ndiv
	write(2,*) 'delta_range ', delta!1.2d0, 0.d0, ndiv2
	write(2,*) 'J2_range ', J2
	write(2,*) 'Lattice  ', lattice
	write(2,*) 'MC_steps ', MCcycles
	close(2)


title = 'diagram'//str!2//'_'//str
title2 =''
do i=1,len(title)   
   if (title(i:i).ne.' ')title2=trim(title2)//trim(title(i:i))   
end do
title11 =  (title2)
Open(1, file =title11)
write(1,*) 'T ','J1 ', 'J2 ', 'delta ', 'stag_mag ','stripe ','chiral ', 'Energy ','Cv ', 'SF '
close(1)

title = 'susceptibilities'//str!2//'_'//str
title2 =''
do i=1,len(title)   
   if (title(i:i).ne.' ')title2=trim(title2)//trim(title(i:i))   
end do
title12 =  (title2)
open(14, file = title12)
write(14,*)'T ','J1 ','J2 ','delta ','stag_mag_chi ','stripe_chi ','chiral_chi ', 'mod_in ','mod_out'
close(14)

title = 'Binder_parameters'//str!2//'_'//str
title2 =''
do i=1,len(title)   
   if (title(i:i).ne.' ')title2=trim(title2)//trim(title(i:i))   
end do
title13 =  (title2)
Open(16, file =title13)
write(16,*)'T ', 'J1 ', 'J2 ', 'delta ','stagmag_binder ','stripe_binder ','ising_binder ','En_binder '
close(16)

title = 'errors'//str
title2 =''

do i=1,len(title)   
   if (title(i:i).ne.' ')title2=trim(title2)//trim(title(i:i))   
end do
title14 =  (title2)
Open(18, file =title14)
write(18,*)'T ','J1 ','J2 ','delta ','magn_error ','stripe_error ','chiral_error ','EnErr ', 'Cv_err '
	close(18)

title = 'chi_errors'//str
title2 =''
do i=1,len(title)   
   if (title(i:i).ne.' ')title2=trim(title2)//trim(title(i:i))   
end do
title19 =  (title2)
Open(19, file =title19)
write(19,*)'T ','J1 ','J2 ','delta ','magn_error ','stripe_error ','chiral_error ','EnErr ', 'Cv_err '
	close(19)

title = 'Binderr'//str
title2 =''
do i=1,len(title)   
   if (title(i:i).ne.' ')title2=trim(title2)//trim(title(i:i))   
end do
title15 =  (title2)
Open(20, file =title15)
write(20,*)'T ','J1 ','J2 ','delta ', 'neel_binder ', 'stripe_binder ', 'ising_binder ', 'energy_binder '
	close(20)

title = 'spin'//str
title2 =''
do i=1,len(title)   
   if (title(i:i).ne.' ')title2=trim(title2)//trim(title(i:i))   
end do
title16 =  (title2)


close(1)
close(4)
close(14)
close(16)
close(19)
close(18)

CALL RANDOMNUMB(nall, rands,k)
CALL INITIALIZE(spin, lattice)
call callenergy(lattice, spin, J1, J2, delta, energy0, E1sum, E2sum)	

PRINT*, "ORDERED ENERGY = ", energy0/den

call readvariables(lattice, spin,stripe, stagmag, amag, magz, magn, chiral)	
!print*, lattice, energy0, stripe, stagmag, amag, chiral


DO N = 1,1!ndiv2+1


P_En(:)= 0.d0

!T = 0.1 ! To find out the ground state config
!do ii = 1,200000
!				call heatbath(lattice, spin, counter, energy0,T,J1,J2,delta, pi)
			!	CALL METROPOLIS(lattice, spin, counter, energy0,T,J1,J2,delta, P_En)
!				call overrelax(lattice, spin, counter, energy,T,J1,J2,delta, P_En, 8)
!end do

dO N2 = 1,1

	write(str2,*)N2
	T = t1 - (i22-1)*(t1-t0)/ndiv
	beta = 1.D0/T
		
	do i0 = 1, 100000
		 call heatbath(lattice, spin, counter, energy0,T,J1,J2,delta, pi)
	end do

	DO i0 = 1, Eqcycles
		!CALL METROPOLIS(lattice, spin, counter, energy0,T,J1,J2,delta, P_En)
		call heatbath(lattice, spin, counter, energy0,T,J1,J2,delta, pi)
		call overrelax(lattice, spin, counter, energy0,T,J1,J2,delta, P_En, 15)
!		call overrelax(lattice, spin, counter, energy,T,J1,J2,delta, P_En)
!		call overrelax(lattice, spin, counter, energy,T,J1,J2,delta, P_En)
!		call overrelax(lattice, spin, counter, energy,T,J1,J2,delta, P_En)
!		call overrelax(lattice, spin, counter, energy,T,J1,J2,delta, P_En)
		if (mod(i0,Eqcycles/10) .eq. 2) then
			
			call overrelax(lattice, spin, counter, energy,T,J1,J2,delta, P_En, 40)

		endif
	END DO	

!	open(2, file = title16)
	do ii = 1, lattice
		do jj = 1,lattice
			do kk = 1,lattice
!				write(2,*) ii, jj, spin(ii,jj,:)
			end do
		end do
	end do
!	close(2)

	
DO N3 = 1, errorbins	
	
	P_En(:) = 0.d0
	call callenergy(lattice, spin, J1, J2, delta, energy0, E1sum, E2sum)

!	print*, 'energy1 = ', energy0/den, stripe/den, stagmag/den, lattice, den, J1, J2, delta

	!CALL RANDOMNUMB(nall, rands,k)
	
	counter=0
	chiral=0.d0
	
	magn0   = 0.d0
	magn1   = 0.d0 
	magn2   = 0.d0 
	magn4   = 0.d0 
	stripe0 = 0.d0 
	stripe1 = 0.d0 
	stripe2 = 0.d0 
	stripe4 = 0.d0 
	chir0   = 0.d0 
	chir1   = 0.d0 
	chir2   = 0.d0 
		chir4   = 0.d0 
	En1     = 0.d0 
	En2     = 0.d0 
	En4     = 0.d0
	Gt1     = 0.d0
	magqt1  = 0.d0
	magnx1  = 0.d0
	magnx2  = 0.d0
	magny1  = 0.d0
	magny2  = 0.d0
	magnz1  = 0.d0
	magnz2  = 0.d0
	SF1   = 0.d0
	
	

	
	DO i0 = 1, MCcycles
		!CALL METROPOLIS(lattice, spin, counter, energy0,T,J1,J2,delta, P_En)
		call heatbath(lattice, spin, counter, energy0,T,J1,J2,delta, pi)

		call overrelax(lattice, spin, counter, energy0,T,J1,J2,delta, P_En, 10)
		
!		call overrelax(lattice, spin, counter, energy,T,J1,J2,delta, P_En)

		call readvariables(lattice, spin, stripe,stagmag,amag,magz, magn, chiral)

!		call Sfval(spin, lattice,kval, SF0)
		!call corr_fun(lattice, spin, Gt, magqt)
		
		
!		SF1   = SF1   + abs(SF0)
		magn1 = magn1 + abs(stagmag)
		magn2 = magn2 + stagmag*stagmag
		magn4 = magn4 + stagmag*stagmag*stagmag*stagmag
		
		magnx1 = magnx1 + abs(magn(1))
		magnx2 = magnx2 + magn(1)**2
		magny1 = magny1 + magn(2)
		magny2 = magny2 + magn(2)**2
		magnz1 = magnz1 + magn(3)
		magnz2 = magnz2 + magn(3)**2
		
		stripe1 = stripe1 + abs(stripe)
		stripe2 = stripe2 + stripe*stripe
		stripe4 = stripe4 + stripe*stripe*stripe*stripe
		
		chir1 = chir1 + abs(chiral)
		chir2 = chir2 + chiral**2
		chir4 = chir4 + chiral**4
		
		En1   = En1 + abs(energy0)
		En2   = En2 + energy0*energy0
		En4   = En4 + energy0*energy0*energy0*energy0		

		E1av = E1av + E1sum*E1sum
		E2av = E2av + E2sum				
		
!		Gt1   = Gt1 + Gt
!		magqt1= magqt1 + magqt
			
	END DO	
	
	call cpu_time(time1)

	
	magn1 =abs( magn1/float(MCcycles))
	mag_avs(N3) = magn1
	
	SF_avs(N3)  = SF1/float(MCcycles)

	magn2 = magn2/float(MCcycles)
	magn4 = magn4/float(MCcycles)
	stripe1 = abs(stripe1/float(MCcycles))
	stripe_avs(N3) =  stripe1
	
	stripe2 = stripe2/float(MCcycles)
	stripe4 = stripe4/float(MCcycles)
	chir1 = abs(chir1/float(MCcycles))
	chiral_avs(N3) =  chir1
	
	chir2 = chir2/float(MCcycles)
	chir4 = chir4/float(MCcycles)
	En1   = abs(En1/float(MCcycles))
	En_avs(N3) =  En1
	En2   = En2/float(MCcycles)
	En4   = En4/float(MCcycles)
!	Gt1   = Gt1/float(MCcycles)
	magqt1= magqt1/float(MCcycles)
	
	E1av = E1av/float(MCcycles)

	E2av = E2av/float(MCcycles)
	stiff = E2av/den - (E1av/T)/den
 	
	magnx1 = magnx1/float(MCcycles)
	magnx2 = magnx2/float(MCcycles)
	magny1 = magny1/float(MCcycles)
	magny2 = magny2/float(MCcycles)
	magnz1 = magnz1/float(MCcycles)
	magnz2 = magnz2/float(MCcycles)
	
	
!	Gtcon  = Gt1 - magqt1**2
!	xit = sqrt(magn2 - Gt1)/(Gt1*kval*kval)
!	xitcon = sqrt(abs(magn2 - magn1**2 - Gtcon)/(Gtcon*kval*kval))
	
	Cv = (En2 - En1**2)*beta*beta
	Cv_avs(N3) = Cv
	

	AFMBind0(N3) = 1 - magn4/(3.D0*magn2*magn2)
	AFMBind1(N3) = 1 - stripe4/(3.D0*stripe2*stripe2)
	AFMBind2(N3) = 1 - chir4/(3.D0*chir2*chir2)
	EnBind(N3)   = 1 - En4/(3.D0*En2*En2)
	
	AFMchi1  =  (magn2 - magn1*magn1)*beta
	chi1_avs(N3) =  AFMchi1
	AFMchi2  =  (stripe2 - stripe1*stripe1)*beta
	chi2_avs(N3) =  AFMchi2
	AFMchi3  =  (chir2 - chir1*chir1)*beta	
	chi3_avs(N3) =  AFMChi3
	
	out_chi = (magnz2 - magnz1**2)
	in_chi  = (magnx2 - magnx1**2 + magny2 - magny1**2)
	
	inchi_avs(N3) = in_chi
	outchi_avs(N3) = out_chi

	print*, N3, errorbins, T, lattice, En1, En2, Cv, EnBind(N3), AFMBind0(N3), AFMBind1(N3), stripe1, chir1, &
	AFMchi2, AFMchi3, in_chi, out_chi, time1-time0
end DO	
	title = 'spin'//str2
	title2 =''
	do i=1,len(title)   
	   if (title(i:i).ne.' ')title2=trim(title2)//trim(title(i:i))   
	end do
	open(50, file = title2)
	do ii = 1, lattice
		do jj = 1,lattice
			do kk = 1,lattice
				write(50,*) ii, jj, kk, spin(ii,jj,:)
			end do
		end do
	end do
	close(50)

	call averages(mag_avs, errorbins, magnav)	
	call averages(stripe_avs, errorbins, stripeav)
	call averages(chiral_avs, errorbins, chirav)
	call averages(En_avs, errorbins, Enav)
	call averages(Cv_avs, errorbins, Cvav)
	call averages(chi1_avs, errorbins, AFMchi1)
	call averages(chi2_avs, errorbins, AFMchi2)
	call averages(chi3_avs, errorbins, AFMchi3)
	call averages(inchi_avs, errorbins, inchiav)
	call averages(outchi_avs, errorbins, outchiav)
	call averages(SF_avs, errorbins, SFav)

	call averages(AFMBind0, errorbins, AFMBinder1)
	call averages(AFMBind1, errorbins, AFMBinder2)
	call averages(AFMBind2, errorbins, AFMBinder3)	
	call averages(EnBind, errorbins, EnBinder)	
	
	!mag_avs, stripe_avs, chiral_avs, En_avs, stagmag_avs, Cv_avs, chi1_avs, chi2_avs, chi3_avs
	!mag_err, stripe_err, chiral_err, En_err,
	call Error_value(mag_avs, errorbins, mag_err)
	call Error_value(stripe_avs, errorbins, stripe_err)
	call Error_value(chiral_avs, errorbins, chiral_err)
	call Error_value(En_avs, errorbins, En_err)
	call Error_value(Cv_avs, errorbins, Cv_err)
	call Error_value(chi1_avs, errorbins, chi1_err)
	call Error_value(chi2_avs, errorbins, chi2_err)
	call Error_value(chi3_avs, errorbins, chi3_err)
	call Error_value(inchi_avs, errorbins, inchi_err)
	call Error_value(outchi_avs, errorbins, outchi_err)
	call Error_value(SF_avs, errorbins, SFerr)
	
	call Error_value(AFMBind0, errorbins, AFMbinderr0)
	call Error_value(AFMBind1, errorbins, AFMbinderr1)
	call Error_value(AFMBind2, errorbins, AFMbinderr2)
	call Error_value(EnBind, errorbins, Enbinderr)
	
	call cpu_time(time1)

	print*, lattice, T, J1, J2, delta, magn1/den, stripe1/den, chir1/den, En1/den, Cv/den, xit/lattice, time1 - time0

	open(1,  file = title11, position ='append')
	open(14, file = title12, position ='append')
	open(16, file = title13, position ='append')
	open(18, file = title14, position ='append')
	open(19, file = title19, position ='append')
	open(25, file = title15, position ='append')
	write(1,*)T, J1,J2, delta, magnav/den,  stripeav/den, chirav/den, Enav/den, Cvav/den, SFav/(den*den), counter/(MCcycles*den)
	write(14,*)T, J1, J2, delta, AFMchi1/den, AFMchi2/den, AFMchi3/den, inchiav/den, outchiav/den
	write(16,*)T, J1, J2, delta, AFMBinder1, AFMBinder2, AFMBinder3, EnBinder
	write(18,*)T, J1, J2, delta, mag_err/den, stripe_err/den, chiral_err/den, En_err/den, Cv_err/den, SFerr/(den*den),zeronum
	write(19,*)T, J1, J2, delta, chi1_err/den, chi2_err/den, chi3_err/den, inchi_err/den, outchi_err/den
	write(25,*)T, j1, J2, delta, AFMbinderr0, AFMbinderr1, AFMbinderr2, Enbinderr
	close(1)
	close(14)
	close(16)
	close(18)
	close(19)
	close(25)
end do

END DO	

end do
	close(4)
	close(1)
	close(14)
	close(16)
deallocate(spin)
end do
!$OMP END PARALLEL DO

END PROGRAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INITIALIZE(spin, lattice)
integer*8::lattice, i,j, k, i0, nall
real*8,dimension(lattice,lattice,3)::spin
real*8::s, x,y,z, norm, rands
!real*8, dimension(100)::rands
real*8,PARAMETER::pi = acos(-1.d0)
nall =100
i0=1
CALL RANDOMNUMB(nall, rands,i0)

DO i = 1,lattice
	DO j = 1,lattice
!	DO k = 1,lattice
	!CALL RANDOM_NUMBER(theta)	 
	!CALL RANDOM_NUMBER(phi)	 		
	norm = 1.
	do while (norm >0.5)
		CALL RANDOM_NUMBER(x); x= x-0.5
		CALL RANDOM_NUMBER(y); y= y-0.5
		CALL RANDOM_NUMBER(z); z= z-0.5	
		norm = sqrt(x*x + y*y + z*z) 	
	end do
	x= x/norm; y=y/norm; z=z/norm	!spin(i,j,1) = cos(2*pi*theta)*sin(2*pi*phi)
	spin(i,j,1) = x!0
	spin(i,j,2) = y!0
	spin(i,j,3) = z!(-1)**(i+j)
!	end do
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine METROPOLIS(lattice, spin, counter, energy,T,J1,J2,delta, P_En)
implicit none
integer*8::lattice
real*8,DIMENSION(lattice,lattice,3)::spin
integer*8::i,j,i0,k,k2, x,y,z, counter, Ecount
real*8:: J1,J2, delta, energy, T, beta, chiral, ch0, ch1
real*8:: theta, phi, px, py, pz, test, dE, sx,sy,sz, norm, ds
real*8,DIMENSION(3):: sold, snew, magn, field , delspin, stagmag

integer*8, dimension(200)::P_En
!real*8,PARAMETER::pi = acos(-1.d0)
!nall =100
i0=1

!CALL RANDOMNUMB(nall, rands,i0)
DO i=1,lattice**2

	CALL RANDOM_NUMBER(px)
	CALL RANDOM_NUMBER(py)
!	CALL RANDOM_NUMBER(pz)
	x = 1+INT(px*lattice)
	y = 1+INT(py*lattice)
!	z = 1+ int(pz*lattice)
	norm = 1.
	do while (norm >0.5)
		CALL RANDOM_NUMBER(sx); sx = sx-0.5!sx= spin(x,y,1) + ds*(sx-0.5)*2
		CALL RANDOM_NUMBER(sy); sy = sy-0.5!sy= spin(x,y,2) + ds*(sy-0.5)*2!sy-0.5
		CALL RANDOM_NUMBER(sz); sz = sz-0.5!sz= spin(x,y,3) + ds*(sz-0.5)*2!sz-0.5	
		norm = sqrt(sx*sx + sy*sy + sz*sz) 			
	end do
!	sx = spin(x,y,1) + 0.1*sx
!	sy = spin(x,y,2) + 0.1*sy
!	sz = spin(x,y,3) + 0.1*sz
!	norm = sqrt(sx*sx + sy*sy + sz*sz) 			
	snew(1) = sx/norm
	snew(2) = sy/norm
	snew(3) = sz/norm

	sold(:) = spin(x,y,:)
	delspin(:) =  snew(:) - spin(x,y,:)

	CALL	callfield(spin, lattice, x,y,field,J1,J2,delta)
	dE =  DOT_PRODUCT(field,delspin)

	IF (dE<0) THEN
		DO i0 =1,3
			spin(x,y,i0) = snew(i0)
		END DO
		energy = energy + dE
		counter = counter + 1		
	ELSE
		CALL RANDOM_NUMBER(test)	
		IF (test<EXP(-dE/T)) THEN
			DO i0 =1,3
				spin(x,y,i0) = snew(i0)			
			END DO
			energy = energy + dE
			counter = counter + 1		
		END IF
	END IF
	
!	Ecount  = int(abs(energy*100/float(lattice*lattice)))
!	if (Ecount <2000) then
!		P_En(Ecount) =  P_En(Ecount) + 1
!	else
!		print*, " - - - Energy OVERFLOW - - -", energy
!	end if

END DO

END SUBROUTINE

Subroutine heatbath(lattice, spin, counter, energy,T,J1,J2,delta, pi)

integer*8::lattice
real*8,DIMENSION(lattice,lattice,3)::spin
integer*8::i,j,i0,k,k2, x,y,z, counter, Ecount
real*8:: J1,J2, delta, energy, T, beta, chiral, ch0, ch1, h, cosx1, sinx1, x2, cosx2, sinx2,cosy1, siny1, cosy2, siny2
real*8:: theta, phi, px, py, pz, test, dE, sx,sy,sz, norm, ds, pi
real*8,DIMENSION(3):: sold, snew, magn, field , delspin, stagmag

integer*8, dimension(200)::P_En
!real*8,PARAMETER::pi = acos(-1.d0)
!nall =100
i0=1

!CALL RANDOMNUMB(nall, rands,i0)
DO i=1,lattice**2

	CALL RANDOM_NUMBER(px)
	CALL RANDOM_NUMBER(py)
!	CALL RANDOM_NUMBER(pz)
	x = 1+INT(px*lattice)
	y = 1+INT(py*lattice)
!	z = 1+ int(pz*lattice)
	norm = 1.
	sold(:) = spin(x,y,:)
	
	CALL	callfield(spin, lattice, x,y,field,J1,J2,delta)
	!field = field/T
	h =  sqrt(dot_product(field, field))
	!py=  exp(h)
	cosx1 = 2
	
	do while(abs(cosx1) > 1.d0)
		call random_number(px)
		!cosx1=1.+1./h*log(px);
		!cosx1 = (1/h)*(log(py*(1-px) + px/py))
		cosx1=-1. - (log(1.+px*(exp(-2.*h/T)-1.)) )*T/h;
!		cosx1= -1 + (log(1 + px*( exp(-2*h) - 1.)))/h
	end do
	!print*, cosx1, px, py, h

	sinx1=sqrt(1.-cosx1*cosx1);
	call random_number(px)
	x2=pi*(px*2.d0 - 1.d0);
	cosx2=cos(x2);
	if (x2>0) then
		sinx2 = sqrt(1.-cosx2*cosx2)
	else
		sinx2 = -sqrt(1.-cosx2*cosx2)
	endif

	cosy1=field(3)/h;
	siny1=sqrt(1.-cosy1*cosy1);
	cosy2=field(1)/(h*siny1);
	siny2=field(2)/(h*siny1);

	snew(3)=cosx1*cosy1 - cosx2*sinx1*siny1;
	snew(1)=cosy2*(cosx2*cosy1*sinx1 + cosx1*siny1) - sinx1*sinx2*siny2;
	snew(2)=cosy2*sinx1*sinx2 + (cosx2*cosy1*sinx1 + cosx1*siny1)*siny2;
	
	spin(x,y,:) = snew(:)
	energy = energy + dot_product(field, snew - sold)
	!print*, dot_product(snew, snew), cosx1, sinx1, x2, cosx2, sinx2,cosy1, siny1, cosy2, siny2
END DO

END SUBROUTINE


subroutine overrelax(lattice, spin, counter, energy,T,J1,J2,delta, P_En, g1)
implicit none
integer :: g1
integer*8::lattice
real*8,DIMENSION(lattice,lattice,3)::spin
integer*8::i,j,i0,k,k2, x,y,z, counter, Ecount, P_En(200)
real*8:: J1,J2, delta, energy, T, beta, chiral, ch0, ch1
real*8:: theta, phi, px, py, pz, test, dE, sx,sy,sz, norm, ds
real*8,DIMENSION(3):: sold, snew, magn, field , delspin, stagmag

DO i=1,g1*lattice**2
	CALL RANDOM_NUMBER(px)
	CALL RANDOM_NUMBER(py)
!	CALL RANDOM_NUMBER(pz)
	x = 1+INT(px*lattice)
	y = 1+INT(py*lattice)
!	z = 1+ int(pz*lattice)
	sold(:) =  spin(x,y,:)
	CALL	callfield(spin, lattice, x,y,field,J1,J2,delta)
	
	spin(x,y,:) =  2*field(:)*(dot_product(field,sold)/dot_product(field, field)) - sold(:)

end do


end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE callfield(spin, lattice, x,y,field,J1,J2,delta)
implicit none
integer*8::lattice
Real*8,DIMENSION(lattice,lattice,3)::spin
Real*8,DIMENSION(3)::field
integer*8::i,n,k, x,y, xleft,xright,yup,ydown, L
Real*8::J1,J2, delta
L =lattice
xleft = x-1
xright =x+1
yup = y-1
ydown = y+1
IF (x .EQ. 1) THEN
	xleft = lattice
ELSE IF (x .EQ. lattice) THEN       
	xright = 1
END IF

IF (y .EQ. 1) THEN
	yup = lattice
ELSE IF (y .EQ. lattice) THEN       
	ydown = 1
END IF

	field(:) = J1*(spin(x, yup,:) + spin(x, ydown,:) + spin(xright, y,:) + spin(xleft, y,:))

	field(:) = field(:) + J2*(spin(xright, yup,:) + spin(xright, ydown,:) + spin(xleft, ydown,:) + spin(xleft, yup,:))
	
	field(3) = delta*field(3)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine callenergy(lattice, spin, J1, J2, delta, energy, E1sum, E2sum)
integer*8	:: lattice, x,y,z
real*8		:: spin(lattice, lattice, 3), energy, field(3), sp(3), J1, J2, delta, E1sum, E2sum

energy = 0.d0
E1sum  = 0.d0
E2sum  = 0.d0


do x= 1, lattice
do y= 1, lattice
!do z= 1, lattice
	sp(:)  = spin(x,y,:)
	call  callfield(spin, lattice, x,y,field,J1,J2,delta)!surroudings(spin, lattice, field,x,y,z)
	energy = energy + 0.5d0*dot_product(field, sp)
	!print*, 'field', field
!end do
end do
end do
!print*, energy

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE readvariables(lattice, spin, stripe, stagmag, amag, magz, magn, chiral)
implicit none
integer*8::lattice, i,j,k, x,y,z
real*8,DIMENSION(lattice,lattice,3)::spin
real*8,DIMENSION(3)::magn, stagmag1, stripe1, stripe2, stagmag3, zmag, stripe3
real*8::chiral, chiral0, s, stripe, stagmag, magz, amag, AFmagn0, AFmag, dot

magn(:) =0.d0	
stripe1(:) =0.d0
stripe2(:) =0.d0
stripe3(:) = 0.d0
stagmag3(:) =0.d0

!magz=0.d0
chiral =0.D0
AFmag =0.d0
dot = 0.d0
DO x=1,lattice
	DO y=1,lattice
!	do z=1,lattice
!		magz = magz + spin(x,y,3)*((-1)**(x+y))
		call chiralxy(lattice, spin, x,y, chiral0)
		!call NeelMagn(spin, lattice, x,y,AFmagn0)
		!call cosphi(spin, lattice, x,y, dot)
		chiral = chiral + chiral0
!		AFmag =  AFmag + dot*0.5d0!AFmagn0
		DO i=1,3
			magn(i) =  magn(i) + ((-1)**(x+y))*spin(x,y,i)
			!stagmag3(i) = stagmag3(i) + ((-1)**(x+y))*spin(x,y,i)!(1)**(x+y)*spin(x,y,i)
			stripe1(i) = stripe1(i) +  ((-1)**(x))*spin(x,y,i)
			stripe2(i) = stripe2(i) +  ((-1)**(y))*spin(x,y,i)
!			stripe3(i) = stripe3(i) +  ((-1)**z)*spin(x,y,z,i)
		END DO
!	end do
	END DO
END DO

magz =  sqrt(magn(3)*magn(3))!sqrt(dot_product(magn,magn))
stripe = Sqrt(	Dot_Product(stripe1,stripe1) + Dot_Product(stripe2,stripe2)	)
stagmag = sqrt(dot_product(magn,magn))!sqrt(dot_product(stagmag3,stagmag3))
magn(:) = sqrt(stripe1(:)**2 + stripe2(:)**2)

!amag = (AFmag)

END SUBROUTINE


subroutine Sfval(spin, lattice,kval, Sval)
		integer*8:: x,y,z, x1, y1, styp1, styp2, lattice
		real*8   :: Sval, kval(3), dot, r1(3), r2(3), r12(3), dot2, sp1(3), sp2(3), realx, realy, realz, imgx, imgy, imgz
		real*8, dimension(lattice, lattice, lattice, 3)   :: spin!, R_val
		real*8   :: rand, cos_i, sin_i
				
		realx = 0.d0
		realy = 0.d0
		realz = 0.d0
		imgx  = 0.d0
		imgy  = 0.d0
		imgz  = 0.d0
		Sval = 0.d0
		do x = 1, lattice
		r1(1)= x
		do y = 1, lattice	
		r1(2)= y
		do z = 1, lattice
!			print*, x,y,z
				r1(3)    = z
!			do styp1 = 1, sublattice
				!r1(:) = R_val(x,y,z,:)
				sp1(:) = spin(x,y,z,:)
						dot    = dot_product(r1, kval)
						cos_i  = cos(dot)
						sin_i  = sin(dot)
						realx  = realx + sp1(1)*cos_i
						realy  = realy + sp1(2)*cos_i
						realz  = realz + sp1(3)*cos_i
						imgx   = imgx  + sp1(1)*sin_i
						imgy   = imgy  + sp1(2)*sin_i
						imgz   = imgz  + sp1(3)*sin_i
!				print*, r1, R_val(x,y,z,1,1)
!			end do
		
		end do	
		end do
		end do
		Sval  = (realx**2 + realy**2 +realz**2) + (imgx**2 + imgy**2 + imgz**2)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE chiralxy(lattice, spin, x,y, chiral)
implicit none
integer*8::lattice, i0
real*8,DIMENSION(lattice,lattice,3)::spin
integer*8:: x,y, xleft,xright,yup,ydown, L
real*8, DIMENSION(3):: s1, s2
real*8::dot, chiral

xleft = x-1
xright =x+1
yup = y-1
ydown = y+1
IF (x .EQ. 1) THEN
	xleft = lattice
ELSE IF (x .EQ. lattice) THEN       
	xright = 1
END IF

IF (y .EQ. 1) THEN
	yup = lattice
ELSE IF (y .EQ. lattice) THEN       
	ydown = 1
END IF

DO i0 = 1,3
	s1(i0) = spin(x,y,i0) - spin(xright,ydown,i0)
	s2(i0) = spin(xright,y,i0) - spin(x,ydown,i0)
END DO
dot =  DOT_PRODUCT(s1,s2)
chiral = dot/4.0d0

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cosphi(spin, lattice, x,y, dot)
implicit none
integer*8::lattice, i,j,k, x,y, xleft,xright,yup,ydown, L
real*8,DIMENSION(lattice,lattice,3)::spin
real*8:: dot
xleft = x-1
xright =x+1
yup = y-1
ydown = y+1
IF (x .EQ. 1) THEN
	xleft = lattice
ELSE IF (x .EQ. lattice) THEN       
	xright = 1
END IF

IF (y .EQ. 1) THEN
	yup = lattice
ELSE IF (y .EQ. lattice) THEN       
	ydown = 1
END IF

dot =       (( spin(x,y,1)*spin(xright,y,1)+ spin(x,y,2)*spin(xright,y,2)+ spin(x,y,3)*spin(xright,y,3) ))
dot = dot + (( spin(x,y,1)*spin(x,ydown,1) + spin(x,y,2)*spin(x,ydown,2) + spin(x,y,3)*spin(x,ydown,3)  ))
dot = dot + (( spin(x,y,1)*spin(xleft,y,1) + spin(x,y,2)*spin(xleft,y,2) + spin(x,y,3)*spin(xleft,y,3)  ))
dot = dot + (( spin(x,y,1)*spin(x,yup,1)   + spin(x,y,2)*spin(x,yup,2)   + spin(x,y,3)*spin(x,yup,3)    ))

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Savespin(spin,lattice,k)
implicit none
integer*8::lattice
integer*8::i,j,k,n 
real*8, Dimension(lattice,lattice,3) :: spin
Do i = 1,lattice
	Do j=1,lattice
		write(k,*) spin(i,j,1), spin(i,j,2), spin(i,j,3)
	END Do
END Do

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine averages(magn_data, MCcycles, magnav)
implicit none
integer*8:: MCcycles, i
real*8, dimension(MCcycles)::magn_data
real*8::magnav
magnav=0.
do i =1,MCcycles
	magnav = magnav + abs(magn_data(i))
end do
magnav = magnav/float(MCcycles)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Error_value(magn_data, MCcycles, error)
implicit none
integer*8:: MCcycles, i
real*8, dimension(MCcycles)::magn_data
real*8::chi, magsum, magsqsum, T, mnn, error
magsum=0.; magsqsum=0.
do i=1,MCcycles
	mnn = abs(magn_data(i))
	magsum = magsum + mnn
	magsqsum = magsqsum + (mnn)**2
end do
error = sqrt((magsqsum/float(MCcycles) - (magsum/float(MCcycles))**2)/float(MCcycles-1))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RANDOMNUMB(nall, rands,i0)
implicit none
integer:: n1, i
real*8,dimension(nall) :: rands0, arr
integer,dimension(nall)::iseed
integer*8 :: clock, clock2, nall, i0
real*8::i2, rands

n1=nall
call random_seed(size = n1)
call system_clock(COUNT=clock)

iseed = clock + 10000*i0 + 100 * [(i, i = 0,nall-1)]

call random_seed(PUT = iseed)
n1=nall
call random_number(rands0) 
rands = rands0(1)
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
