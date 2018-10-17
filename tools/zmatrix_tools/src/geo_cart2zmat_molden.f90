program geo_cart2zmat
  
  implicit none
  
  integer,parameter  ::  dp = selected_real_kind(2*precision(1.0)) 
  integer,parameter  ::  maxatoms=40
  integer,parameter  ::  unit1=13
  integer,parameter  ::  unit2=14

  real(kind=dp),dimension(3,maxatoms)   ::  czcart=0.0_dp
  real(kind=dp),dimension(3,maxatoms)   ::  czint=0.0_dp
  integer,dimension(3,maxatoms)         ::  izcmat=0
  character(len=2),dimension(maxatoms)  ::  atomlabel
  integer  ::  i,k
  integer  ::  read_info,tmp_int
  integer  ::  natoms1,natoms2
  character(len=50)  ::  tmp_char
  character(len=200)  ::  file1,file2
  character(len=57)  ::  comment1,comment2
 

  call getarg(n=1,buffer=file1)
  call getarg(n=2,buffer=file2)


  open(unit=unit1,file=file1,status='old',action='read')
  open(unit=unit2,file=file2,status='old',action='read')

  read(unit1,*) natoms1
  read(unit1,*) comment1
  do i=1,natoms1
     read(unit1,*,iostat=read_info) atomlabel(i), (czcart(k,i),k=1,3)
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


  call cart2zmat(natoms1,czcart,izcmat,czint)

  if (natoms1.ge.10) then
     write(6,'(1x,i2,a56)') natoms1, ' '
  else
     write(6,'(1x,i1,a57)') natoms1, ' '
  endif
  write(6,'(2x,a57)') comment1
  write(6,'(1x,a2,a56)') atomlabel(1), ' '
  write(6,'(1x,a2,4x,(f11.7,2x),26x,4x,(i2,2x),a6)') atomlabel(2), czint(1,2), izcmat(1,2), ' '
  write(6,'(1x,a2,4x,(f11.7,2x),(f11.6,2x),13x,4x,2(i2,2x),a2)') atomlabel(3),  (czint(k,3),k=1,2), (izcmat(k,3),k=1,2), ' '
  do i=4,natoms1
     write(6,'(1x,a2,4x,(f11.7,2x),2(f11.6,2x),4x,3(i2,2x))') atomlabel(i), (czint(k,i),k=1,3), (izcmat(k,i),k=1,3)
  enddo

  !write(6,'(1x,i2)') 1
  !write(6,'(1x,i2,3x,(i2,3x),10x,2x,(f11.7,3x))') 2, izcmat(1,2), czint(1,2)
  !write(6,'(1x,i2,3x,2(i2,3x),5x,2x,(f11.7,3x),(f12.6,3x))') 3, (izcmat(k,3),k=1,2), (czint(k,3),k=1,2)
  !do i=4,natoms1
  !   write(6,'(1x,i2,3x,3(i2,3x),2x,(f11.7,3x),2(f12.6,3x))') i, (izcmat(k,i),k=1,3), (czint(k,i),k=1,3)
  !enddo



end program geo_cart2zmat
