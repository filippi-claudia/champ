program geo_zmat2cart
  
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

  read(unit1,*,iostat=read_info) 
  read(unit1,*,iostat=read_info) tmp_int, czint(1,2)
  read(unit1,*,iostat=read_info) tmp_int, (czint(k,3),k=1,2)
  do i=4,maxatoms
     read(unit1,*,iostat=read_info) tmp_int, (czint(k,i),k=1,3)
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
     write(6,*) 'ERROR IN DATA FILES' ; stop
  endif

  call zmat2cart(natoms1,izcmat,czint,czcart)

  write(6,'(''Atom'',t12,''x'',t25,''y'',t38,''z'')')
  do i=1,natoms1
     write(6,'(1x,i2,4x,3(f11.7,2x))') i, (czcart(k,i),k=1,3)
  enddo



end program geo_zmat2cart
