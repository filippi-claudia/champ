program geo_zmat2cart_rc
  
  implicit none
  
  integer,parameter  ::  dp = selected_real_kind(2*precision(1.0)) 
  integer,parameter  ::  maxatoms=40
  integer,parameter  ::  unit1=13
  integer,parameter  ::  unit2=14
  integer,parameter  ::  unit3=15

  real(kind=dp),dimension(3,maxatoms)  ::  czcart=0.0_dp
  real(kind=dp),dimension(3,maxatoms)  ::  czint=0.0_dp
  real(kind=dp),dimension(3,3)         ::  czcart_rc=0.0_dp
  integer,dimension(3,maxatoms)        ::  izcmat=0
  character(len=2),dimension(maxatoms) ::  atomlabel
  integer  ::  i,k
  integer  ::  read_info,tmp_int
  integer  ::  natoms1,natoms2
  character(len=50)   ::  tmp_char
  character(len=200)  ::  file1,file2,file3
  character(len=57)   ::  comment1,comment2,comment3
 
  call getarg(n=1,buffer=file1)
  call getarg(n=2,buffer=file2)
  call getarg(n=3,buffer=file3)



  open(unit=unit1,file=file1,status='old',action='read')
  open(unit=unit2,file=file2,status='old',action='read')
  open(unit=unit3,file=file3,status='old',action='read')

  read(unit1,*,iostat=read_info) natoms1
  read(unit1,*,iostat=read_info) comment1
  read(unit1,*,iostat=read_info) atomlabel(1)
  read(unit1,*,iostat=read_info) atomlabel(2), czint(1,2)
  read(unit1,*,iostat=read_info) atomlabel(3), (czint(k,3),k=1,2)
  do i=4,natoms1
     read(unit1,*,iostat=read_info) atomlabel(i), (czint(k,i),k=1,3)
     if(read_info.ne.0 .and. i.ne.natoms1) then
        write(6,*) 'ERROR IN XYZ INPUT FILE' ; stop
     endif
  enddo

  read(unit2,*) natoms2
  if (natoms1.ne.natoms2) then
     write(6,*) 'ERROR IN INPUT FILES' ; stop
  endif
  read(unit2,*) comment2
  read(unit2,*) tmp_char
  read(unit2,*) tmp_char, izcmat(1,2)
  read(unit2,*) tmp_char, (izcmat(k,3),k=1,2)
  do i=4,natoms1
     read(unit2,*,iostat=read_info) tmp_char, (izcmat(k,i),k=1,3)
     if(read_info.ne.0 .and. i.ne.natoms1) then
        write(6,*) 'ERROR IN  CONNECTION MATRIX INPUT FILE' ; stop
     endif
  enddo

  read(unit3,*,iostat=read_info) tmp_int
  read(unit3,*,iostat=read_info) tmp_char
  do i=1,3
     read(unit3,*,iostat=read_info) tmp_char, (czcart_rc(k,i),k=1,3)
  enddo


  call zmat2cart_rc(natoms1,izcmat,czint,czcart,czcart_rc)

  if (natoms1.ge.10) then
     write(6,'(1x,i2,a56)') natoms1, ' '
  else
     write(6,'(1x,i1,a57)') natoms1, ' '
  endif
  write(6,'(2x,a57)') comment1
  do i=1,natoms1
     write(6,'(1x,a2,4x,3(f11.7,2x))') atomlabel(i), (czcart(k,i),k=1,3)
  enddo


end program geo_zmat2cart_rc
