!*************************************************************************
!This is the main program of the whole procedures
!*************************************************************************
program cindo 
  !
  IMPLICIT NONE
  !
  integer::i,j,natom,natpar,ok,iham,nelec,ionicity,nbast,na,nb,la,lb,m,ibast,&
       jbast,mm,imethod,nalpha,nbeta,nocc,ncount,dipol,NITER,nzmax,maxit,&
       idamp,iall,vecu,ioptics,imullik,iplot,ndir,nplane,nplot,nplot_u,&
       nplot_d,ilogplot
  integer,allocatable::nc(:),lc(:),iatom(:),nbasis(:),mc(:),lcmax(:),&
       nc_atm(:),idir(:,:)
  parameter(natpar=110)
  real(kind=8)::coeff1(0:2),coeff2(0:1),coeff(0:3),coef(0:4,0:4),enuc,&
       xdamp,conv,buf,rmax,rmin,dr
  real(kind=8)::lwidth,omega1,omega2,domega,scale
  real(kind=8),allocatable::zn(:),rnuc(:,:),zeta(:),rovp(:),ycoef(:,:),&
       zcoef(:,:),sovp(:),smat(:),g(:,:),zet_atm(:),hcore(:),ezn(:),huck(:),&
       dip(:,:),g1(:),f2(:),eneg(:,:),bondpar(:),dipr(:,:),zval(:),eng_at(:)  
  character(len=80)::line
  real(kind=8)::xmin,xmax,dx,ymin,ymax,dy
  integer,allocatable::iorbplt(:),iorbplt_u(:),iorbplt_d(:),iaxis(:)
  !
  character(len=2)::atomname(110)=(/'H ', 'He', 'Li', 'Be', 'B ', 'C ', &
       'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', &
       'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', & 
       'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', & 
       'Rb', 'Sr', 'Y ',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', & 
       'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', &
       'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', & 
       'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', & 
       'Ir', 'Pt', 'Au', 'Hg', 'Ti', 'Pb', 'Bi', 'Po', 'At', 'Rn', & 
       'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', & 
       'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', & 
       'Bh', 'Hs', 'Mt', 'Ds'/) 
  character(len=80)::hamilt
  !
  !***************************************************************************
  !
  ! The total  number of atoms in the molecule
  !
  ! 
! read(*,*)line
  PRINT*, ''
  PRINT*, 'Welcom to the CINDO program!'
  PRINT*, ''
  PRINT*, 'Non-profit use and redistribution without modification are approved.'
  PRINT*, 'Please read https://github.org/brhr-iwao/cindo_windows/readme.md'
  PRINT*, 'for the detailed licence.'
  PRINT*, ''
  PRINT*, 'If you publish a paper which result from using the program,'
  PRINT*, 'cite "S. Sahu and A. Shukla, Fortran 90 implementation of '
  PRINT*, 'the Hartree-Fock approach within the CNDO/2 and INDO models,'
  PRINT*,  'Computer Physics Communications, 180, 724, 2009"'
  PRINT*, '(doi.org/10.1016/j.cpc.2008.11.004)'
  PRINT*, ''
   PRINT*, 'Enter the number of atoms in the molecule/cluster as integer' 
  read*,natom
  allocate(zn(natom),ezn(natom),stat=ok)
  if(ok/=0)then
     print*,'Stopping: allocation failed zn & ezn'
     stop
  end if
  !WRITE(*,*) 'Number of atoms =', natom
  !PRINT '("Number of atoms = ", I2)', natom
  !
  ! Read the nuclear charges

  !
  ! read(*,*)line

  PRINT*, ''
  PRINT '("Enter the nuclear charge of each atom as space-separated", i2, " floating point numbers")', natom
  PRINT*, 'E.g., "1.0 1.0" for H2, "5.0 5.0" for B2, and so on.'
  read*,(zn(i),i=1,natom)
  ! do i=1, natom
    ! PRINT*, 'Enter the nuclear charge of atom' , i
    ! PRINT '( "Enter the nuclear charge of atom", i2) ', i
    ! PRINT*, '(1 for H, 2 for He and so on)'
  !  read*, zn(i)
  ! ENDDO
  !
  ! Nuclear coordinates
  !
  allocate(rnuc(3,natom),stat=ok)
  if(ok/=0)then
     print*,'allocation failed for rnuc'
     stop
  end if
  !
  nzmax=18
  !
  allocate(f2(nzmax),eneg(nzmax,3),g1(nzmax),bondpar(nzmax),eng_at(nzmax),&
       stat=ok)
  if(ok/=0)then
     print*,'allocation failed for empirical parameters'
     stop
  end if
  !
  PRINT*, ''
  PRINT*,  'Enter the nuclear coordinates (Angstrom) one by one.'
  ! read(*,*)line
  do i=1,natom
     PRINT '("Enter the coordinate of atom", i2, " in Angstrom")', i
     PRINT '("as space-separated three floating point numbers. e.g. 0.0 0.0 1.0 ")'
     read*,(rnuc(j,i),j=1,3)
  !PRINT*, 'Enter the x coordinate (Angstrom) of atom', i
  !PRINT '("Enter the x coordinate (Angstrom) of atom", i2)', i
  !read*, rnuc(1,i)
  !PRINT*, 'Enter the y coordinate (Angstrom) of atom', i
 ! PRINT '("Enter the y coordinate (Angstrom) of atom", i2)', i
  !read*, rnuc(2,i)
  !PRINT*, 'Enter the z coordinate (Angstrom) of atom', i
  !PRINT '("Enter the z coordinate (Angstrom) of atom", i2)', i
  !read*, rnuc(3,i)
     !
     ! Convert atomic coordinates from Angstroms to atomic units
     !   
     do j=1,3
        rnuc(j,i)=rnuc(j,i)/0.5291772083d0
!        rnuc(j,i)=rnuc(j,i)/0.529167d0
     end do
     !
  end do

  !
  do i=1,natom
     if(zn(i)==1)then
        ezn(i)=1
     else
        if(zn(i)>=3.and. zn(i)<=10)then
           ezn(i)=zn(i)-2
        else
           if(zn(i)>=11.and. zn(i)<=17)then
              ezn(i)=zn(i)-10
           else
              if(zn(i)>17)then
                 print*,"Stopping,the program stops with this atomic #"
                 stop
              end if
           end if
        end if
     end if
  end do

  !print out the atom names
  !
  write(*,*)'All the calculated values are in atomic units'
  write(*,*)'Conversion factors:'
  write(*,*)'0.5291772083 Ang = 1 A.U.(Bohr)'
  write(*,*)'27.21 eV = 1 A.U(Hartree).'

  !
  print*,' '
  write(*,1)
1 format(T5,"Atom  #",T26,"Atom Name",T44,"Atomic Coordinates (a.u.)")
  write(*,2)
2 format(T47,"X       Y        Z")
  do i=1,natom
     write(*,3)i, atomname(nint(zn(i))), (rnuc(j,i),j=1,3)
  end do
3 format(3X,I4.1,21X,A2,11X,3F8.3)
  !
  ! Semi-empirical methods to be employed
  !
  print*,' '
  !
  ! Read the type of Hamiltonian to use
  ! 
  ! read(*,*)line
  PRINT*,''
  PRINT*, 'Enter the Hamiltonian to be used. i.e. CNDO or INDO'
  read(*,*)hamilt
  hamilt=adjustl(hamilt) 
  if(hamilt(1:4)=='CNDO')then
     print*,"CNDO/2 model will be used"
     iham=1
  elseif(hamilt(1:4)=='INDO')then
     print*,"INDO model will be used"
     iham=2
  else
     print*,'Hamiltonian you gave was: ',hamilt
     print*,"This is not a valid Hamiltonian"
  end if
  !
  ! Read the method to employ
  ! 
  !read(*,*)line
  PRINT*,''
  PRINT*, 'Choose the method RHF or UHF'
  read(*,*)hamilt
  hamilt=adjustl(hamilt) 
  if(hamilt(1:3)=='RHF')then
     print*,"RHF method will be used"
     imethod=1
  elseif(hamilt(1:3)=='UHF')then
     print*,"UHF method will be used"
     imethod=2
 !    read(*,*)line
     PRINT*, 'Your choice is UHF.'
     PRINT*, 'Enter the numbers of alpha spin and beta spin electrons respectively'
     PRINT*, 'as two integers separated by space.  e.g. 4 2'
     read(*,*)nalpha,nbeta
     print*,"nalpha,nbeta",nalpha,nbeta
  else
     print*,'Method you gave was: ',hamilt
     print*,"This is not a valid method"
  end if
  !

  !
  !Restriction on atomic nmber to be used
  !
  do i=1,natom
     if(zn(i)>18)then
        print*,"Stopping,progrogram stops with this atomic #"
        stop
     endif
     if(zn(i)>10)then
        if(iham==2)then
           print*,"Stopping,with this atomic # INDO cannot be employed"
           stop
        endif
     end if
  end do
  ! 
  !
  ! Read the ionicity of the system
  !
  print*,' '
  !read*,line
  PRINT*, 'Enter the ionicity of the system. e.g. 0 for a neutral, -1 for a single anion.'
  read*,ionicity
  nelec=0
  do i=1,natom 
     nelec=nelec+nint(ezn(i))
  end do
  nelec=nelec-ionicity
  write(*,4)nelec
4 format(1x,'Number of valence electrons in the system =',1x,I4.1)
  print*,' '
  !
  ! Compute the number of occupied orbitals
  ! nocc: is the number of occupied orbitals for the RHF method
  ! nalpha/nbeta are the number of up/down spin orbitals for the UHF method
  !
  if(imethod==1)then
     if(mod(nelec,2)>0)then
        PRINT*, ''
        ! print*,'The system has odd number of electrons=',nelec
	PRINT '("The system has odd number of electrons=", i2)', nelec
        print*,'So your choice of RHF method is not possible'
        stop
     end if
     nocc=nelec/2
  else
     ncount=nalpha+nbeta
     if(nelec/=ncount)then
        !print*,'nalpha,nbeta',nalpha,nbeta
	PRINT*, ''
        !print*,'For the UHF method, you gave nalpha, nbeta=',nalpha,nbeta
	print '("For the UHF method, you gave nalpha, nbeta=", i2, t1, i2)',nalpha,nbeta
        print*,'Error: nalpha,nbeta donot add up to nelectron=',nelec
        stop
     end if
  end if
  !
  !Number of basis functions to be used in the system
  !
  allocate(nbasis(natom),lcmax(natom),stat=ok)
  if(ok/=0)then
     print*,"allocation fails for nbasis,lcmax"
     stop
  end if
  do i=1,natom
     if(zn(i)<=2)then
        nbasis(i)=1
        lcmax(i)=0
     else
        if(zn(i)>=3.and. zn(i)<=10)then
           nbasis(i)=4
           lcmax(i)=1
        else
           if(zn(i)>=11.and. zn(i)<=17)then
              nbasis(i)=9
              lcmax(i)=2
           else
              if(zn(i)>17)then
                 print*,"Stopping,the program stops with this atomic #"
                 stop
              end if
           end if
        end if
     end if
  end do
  !

  nbast=0
  do i=1,natom
     nbast=nbast+nbasis(i)
  end do
  write(*,5)nbast
5 format(1x,'Total no. of basis functions in the system=',1x,I4.1)


  !

  allocate(zeta(nbast),lc(nbast),nc(nbast),iatom(nbast),mc(nbast),&
       rovp((nbast*(nbast+1))/2),sovp((nbast*(nbast+1))/2),&
       smat((nbast*(nbast+1))/2),g(natom,natom),zet_atm(natom),&
       nc_atm(natom),hcore((nbast*(nbast+1))/2),&
       huck((nbast*(nbast+1))/2),dip(3,(nbast*(nbast+1))/2),&
       dipr(3,(nbast*(nbast+1))/2),stat=ok)
  if(ok/=0)then
     print*,"allocation fails for zeta,lc, nc,iatom,rovp,sovp"
     stop
  end if

  !
  !read(*,*)line
  PRINT*, 'Enter the convergence threshold and max number of iteration' 
  PRINT*, 'as space-separated floating point number and integer. e.g. 1.d-7 100'
  read(*,*)conv,maxit
  !
  idamp=0 
  !read(*,*)line
  PRINT*, 'To facillitate the convergence, Would you like to do Fock matrix mixing?'
  PRINT*, 'If you wish to do Fock matrix mixing, enter DAMP. Otherwise enter NODAMP.'
  read(*,*)hamilt
  hamilt=adjustl(hamilt)
  ! 
! 
  if(hamilt(1:4)=='DAMP')then     
     idamp=1
     PRINT*, 'Your choice is DAMP.'
     PRINT*, 'Enter the mixing parameter as a floating point number between 0.0 and 1.0'
     PRINT*, '(The ith iteration result R(i) = param*R(i)+(1.0-param)*R_(i-1))'
     read*,xdamp
     print*, 'Fock matrix mixing (damping) will be done'
     write(*, 3)xdamp
  end if
 !
6 format(1x,'Fock matrix mixing coefficient=',2x,f7.4)
  !
  imullik=0
  !read(*,*)line
  PRINT*, 'Would you like to calculate Mulliken charge?'
  PRINT*, 'If you wish to calculate it, enter MULLIK. Otherwise enter NOMULLIK.'
  read(*,*)hamilt
  hamilt=adjustl(hamilt)
  if(hamilt(1:6)=='MULLIK')then
     imullik=1
  end if
  !
  iplot=0
  !read(*,*)line
  PRINT*, 'Would you like to plot the orbitals/densities as a function of spatial coordinates?'
  PRINT*, 'Enter either of ORBPLT, DENPLT, BOTHPLT or NOPLT'
  read(*,*)hamilt
  print*,' '
  hamilt=adjustl(hamilt) 
  if(hamilt(1:6)=='ORBPLT')then
     print*,"Orbitals will be plotted"
     iplot=1
  elseif(hamilt(1:6)=='DENPLT')then
     print*,"Charge density will be plotted"
     iplot=2
  elseif(hamilt(1:7)=='BOTHPLT')then
     print*,"Both orbitals and charge density plot will be plotted"
     iplot=3
  end if
  !
  if(iplot>0)then
     !
     ilogplot=0
     !
     PRINT*, 'Your choice is DENPLT.'
     PRINT*, 'You can print out the log10 of the density which can be useful in making contour plot.'
     PRINT*, 'Enter 1D or 1DLOG or 2D or 2DLOG'
     read(*,*)hamilt
     hamilt=adjustl(hamilt)
     if(hamilt(1:2)=='2D')iplot=iplot+3
     if(hamilt(3:5)=='LOG')ilogplot=1
     if(ilogplot>0)iplot=7
     if((iplot==1.or.iplot==4).and.ilogplot==1)then
        print*," Semi_emp:Logarithim plot is only allowed for density"
        print*,"So stopping due to input error"
        stop
     end if
     !
     if(iplot<=3)then
        !
        ! One dimensional plots with given direction
        !
	PRINT*, 'Enter the starting, ending position and step size for the specified directions, respectively,'
	PRINT*, 'as space-separated floating point numbers, e.g. -5.d0 5.d0 0.5d0'
        read(*,*)rmin,rmax,dr
	PRINT*, 'If you choosed the 1D plot, enter 1.'
	PRINT*, 'Otherwise enter 2 (i.e. you choosed 2D plot.)'
        read(*,*)ndir
        allocate(idir(3,ndir),stat=ok)
        if(ok/=0)then
           print*,'Program semi_emp: ndir allocation fails'
           stop
        end if
        print*,'Plots will be along the direction(s):'
        !
	PRINT*, 'If you choosed the 1D plot, specify the direction vector as space-separeted three integers.'
	PRINT*, 'e.g., 1 0 0'
	PRINT*, 'Otherwise, that is, you choose the 2D plot, specify the two direction vectors as space-separeted six integers '
	PRINT*, 'e.g., 1 0 0 0 1 0'
        do i=1,ndir
           read(*,*)(idir(j,i),j=1,3)
           write(*,7),(idir(j,i),j=1,3)
        end do
     else
        !
        !2D Plotting (xy-plane, for a given value of z)
        !
        !
	PRINT*, 'Enter starting, ending position and step size for x-direction, e.g. -5.d0 5.d0 0.5d0'
        read(*,*)xmin,xmax,dx
	PRINT*, 'Enter starting, ending position and step size for y-direction, e.g. -5.d0 5.d0 0.5d0'
        read(*,*)ymin,ymax,dy
	PRINT*, 'Enter number of plane(s) mutually pararell on which the plots are needed. e.g. 3'
	PRINT*, '(These planes are characterized by the above two coordinates.)'
        read(*,*)nplane
        allocate(iaxis(nplane),zval(nplane),stat=ok)
        if(ok/=0)then
           print*,'Allocation failure for iaxis and zval'
           stop
        end if
        ! 
	PRINT*,'Enter the index/indices interger for the cartesian axis which is perpendicular to the each planes. e.g. 1 2 3'
        read*,(iaxis(i),i=1,nplane)
	PRINT*, 'Enter the location of the each planes in terms of the value of the cartesian coordinate.'
        read*,(zval(i),i=1,nplane)
        !
     end if
     ! a RHF case
     if(imethod==1)then
        PRINT*, 'How many orbitals will be plotted?'
	PRINT*, 'If you wish to plot all the occupied orbitals, add negative sign to the orbital number.'
	PRINT*, 'Otherwise you wish to plot a part of orbitals, leave the orbital number as positive value.'
        read*,nplot
        if(nplot<=0)then
           nplot=nelec/2
           allocate(iorbplt(nplot),stat=ok)
           if(ok /= 0)then
              WRITE(*,*)'allocation fails for iorbplt'
              stop
           end if
           !
           do i=1,nplot
              iorbplt(i)=i
           end do
        else
           allocate(iorbplt(nplot),stat=ok)
           if(ok /= 0)then
              WRITE(*,*)'Allocation fails for iorbplt'
              stop
           end if
	   PRINT*, 'Enter indices of orbitals which you wish to plot as space-separated numbers. e.g 1 2 3'
           read*,(iorbplt(i),i=1,nplot)
        end if
        !
     ! an UHF case
     elseif(imethod==2)then
        PRINT*, 'This is a UHF calculation.'
	PRINT*, 'Enter the numbers of upper-spin and lower-spin orbitals to be plotted respectively'
	PRINT*, 'as space-separetd two integers.'
        read*,nplot_u,nplot_d
        if(nplot_u<=0)then
           nplot_u=nalpha
           allocate(iorbplt_u(nplot_u),stat=ok)
           if(ok /= 0)then
              WRITE(*,*)'Allocation fails for iorbplt_u'
              stop
           end if
           do i=1,nplot_u
              iorbplt_u(i)=i
           end do
        else
           allocate(iorbplt_u(nplot_u),stat=ok)
           if(ok /= 0)then
              WRITE(*,*)'Allocation fails for iorbplt_u'
              stop
           end if
           read*,(iorbplt_u(i),i=1,nplot_u)
        end if
        if(nplot_d<=0)then
           nplot_d=nbeta
           allocate(iorbplt_d(nplot_d),stat=ok)
           if(ok /= 0)then
              WRITE(*,*)'Allocation fails for iorbplt_d'
              stop
           end if
           do i=1,nplot_d
              iorbplt_d(i)=i
           end do
        else
           allocate(iorbplt_d(nplot_d),stat=ok)
           if(ok /= 0)then
              WRITE(*,*)'Allocation fails for iorbplt_d'
              stop
           end if
           read*,(iorbplt_d(i),i=1,nplot_d)
        end if
     end if
  end if
!
 7 format(2x,3i3.1)
  !
  dipol=0
  ! read(*,*)line
  PRINT*, 'If you wish to calculate the dipole moment, enter DIPOLE. Otherwise enter NODIPPOLE.'
  read(*,*)hamilt
  hamilt=adjustl(hamilt)
  if(hamilt(1:6)=='DIPOLE')then     
     dipol=1
  else
     print*,'Method you gave was: ',hamilt
     print*,"This is not a valid method"
  end if
  !


  !
  ioptics=0
  ! a RHF case
  if(imethod==1)then
    ! read(*,*)line
    PRINT*, 'If you wish to calculate the linear optical absorption (which is applicable only for RHF calculation),'
    PRINT*, 'Enter OPTICS. Otherwise enter NOOPTICS.'
    PRINT*,'If you calculate the linear optical absorption, the calculation result is written in a filed name spectrum.dat'
    read(*,*)hamilt
    hamilt=adjustl(hamilt)
    if(hamilt(1:6)=='OPTICS')then
       ioptics=1
    end if
  end if  
!
  if(ioptics==1) THEN
  PRINT*, 'Enter the plotting frequency step, minmum and maximum limit of frequencies to be plotted'
  PRINT*, 'and line width of  the excited state as space-separeted five floating-point numbers.'
  PRINT*, 'e.g. 0.25 0.0 10.d0 0.01 1.d0'
  read*,lwidth,omega1,omega2,domega,scale  
 ENDIF  

  !

  !
  !
  !
  call basegen(zeta,zn,nc,lc,mc,iatom,natom,nbast,nbasis,zet_atm,nc_atm,&
       bondpar,eng_at,eneg,g1,f2,nzmax,iham)

  !
  ! get the factorials
  !
  call factcal 

  call assoc_legndre

  ! 
  !
  !
  call redovint(rovp,sovp,smat,nc,lc,lcmax,mc,rnuc,zeta,iatom,nbast,&
       nbasis,natom)
  !
  call coul_int(natom,g,rnuc,zet_atm,nc_atm) 
  !            
  call core_int(nbast,natom,zn,ezn,iatom,iham,lc,hcore,huck,smat,g, bondpar,&
  eneg,g1,f2,nzmax,rnuc,enuc)
  !
  iall=0
  !
  if(ioptics>0)iall=1
  if(iplot>0)iall=1
  !
  if(imethod==1)then
     call scf_rhf(natom,nbast,nbasis,iatom,zn,huck,hcore,g,&
          xdamp,conv,maxit,nelec,nocc,iham,iall,idamp,g1,f2,nzmax,enuc,&
          imullik,ezn,eng_at)
     !
     !
  elseif(imethod==2)then
     !print*,'iall',iall
     call scf_uhf(natom,nbast,nbasis,iatom,zn,huck,hcore,g,&
          xdamp,conv,maxit,nelec,nalpha,nbeta,iham,iall,idamp,&
          g1,f2,nzmax,enuc,imullik,ezn,eng_at)
  else
     print*,"No other methods have been employed"
     stop
  end if
  !
  call dipint(nbast,iatom,nc,lc,mc,zeta,dip)
  !

  call  property(ilogplot,nplane,imethod,dipol,iall,iplot,natom,nbasis,&
       ezn,rnuc,nbast,iatom,lc,mc,nc,zeta,ioptics,lwidth,omega1,&
       omega2,scale,nocc,nalpha,nbeta,domega,dip,ndir,idir,zval,&
       rmax,rmin,dr,xmax,xmin,dx,ymax,ymin,&
       dy,iorbplt,iorbplt_u,iorbplt_d,nplot,nplot_u,nplot_d,iaxis)
  ! 
  !
  !             
end program cindo 
 
