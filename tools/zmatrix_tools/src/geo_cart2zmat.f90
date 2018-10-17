program geo_cart2zmat
  
  implicit none
  
  integer,parameter  ::  dp = selected_real_kind(2*precision(1.0)) 
  integer,parameter  ::  maxatoms=100
  integer,parameter  ::  unit1=13
  integer,parameter  ::  unit2=14

  real(kind=dp),dimension(3,maxatoms)  ::  czcart=0.0_dp
  real(kind=dp),dimension(3,maxatoms)  ::  czint=0.0_dp
  integer,dimension(3,maxatoms)        ::  izcmat=0
  integer  ::  i,k
  integer  ::  read_info,tmp_int
  integer  ::  natoms1,natoms2
  character(len=200)  ::  file1,file2
 
  call getarg(n=1,buffer=file1)
  call getarg(n=2,buffer=file2)

  open(unit=unit1,file=file1,status='old',action='read')
  open(unit=unit2,file=file2,status='old',action='read')

  do i=1,maxatoms
     read(unit1,*,iostat=read_info) tmp_int, (czcart(k,i),k=1,3)
     if(read_info.ne.0) then
        natoms1=i-1
        exit
     endif
  enddo

  
  read(unit2,*,iostat=read_info) 
  read(unit2,*,iostat=read_info) tmp_int, izcmat(1,2)
  read(unit2,*,iostat=read_info) tmp_int, (izcmat(k,3),k=1,2)
  do i=4,maxatoms
     read(unit2,*,iostat=read_info) tmp_int, (izcmat(k,i),k=1,3)
     if(read_info.ne.0) then
        natoms2=i-1
        exit
     endif
  enddo

  if(natoms1.ne.natoms2) then
     write(6,*) 'ERROR IN INPUT DATA FILES' ;  stop     
  endif

  call cart2zmat(natoms1,czcart,izcmat,czint)

  write(6,'(''Atom'',t9,''bond length'',t21,''bond angle'',t34,''dihed. angle'',t50,''Connection matrix'')')
  write(6,'(1x,i2,a57)') 1, ' '
  write(6,'(1x,i2,4x,(f11.7,2x),26x,4x,(i2,2x),a6)') 2, czint(1,2), izcmat(1,2), ' '
  write(6,'(1x,i2,4x,(f11.7,2x),(f11.6,2x),13x,4x,2(i2,2x),a2)') 3,  (czint(k,3),k=1,2), (izcmat(k,3),k=1,2), ' '
  do i=4,natoms1
     write(6,'(1x,i2,4x,(f11.7,2x),2(f11.6,2x),4x,3(i2,2x))') i, (czint(k,i),k=1,3), (izcmat(k,i),k=1,3)
  enddo

  !write(6,'(1x,i2)') 1
  !write(6,'(1x,i2,3x,(i2,3x),10x,2x,(f11.7,3x))') 2, izcmat(1,2), czint(1,2)
  !write(6,'(1x,i2,3x,2(i2,3x),5x,2x,(f11.7,3x),(f12.6,3x))') 3, (izcmat(k,3),k=1,2), (czint(k,3),k=1,2)
  !do i=4,natoms1
  !   write(6,'(1x,i2,3x,3(i2,3x),2x,(f11.7,3x),2(f12.6,3x))') i, (izcmat(k,i),k=1,3), (czint(k,i),k=1,3)
  !enddo



end program geo_cart2zmat
