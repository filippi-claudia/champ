subroutine get_masses(atomtype, mass)

implicit none
  
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = selected_real_kind(15, 307)
  character,intent(in) :: atomtype
  real(kind=dp), dimension(15):: masses
  real(kind=dp) :: mass, elec_mass

 elec_mass = 1822.88 
 masses(1) = 1.0079   !Hydorgen
 masses(2) = 12.0107  !Carbon
 masses(3) = 14.0067  !Nitrogen
 masses(4) = 15.9994  !Oxygen
 masses(5) = 32.065   !Sulfur

! write(6,*) "ATOMTYPE = ", atomtype
 mass = 0
 if (atomtype== "H") then
   mass = masses(1)
 elseif(atomtype== "C") then
   mass= masses(2) 
 elseif(atomtype=="N") then
   mass= masses(3)
 elseif(atomtype == "O") then
   mass=masses(4) 
 elseif(atomtype == "S") then
   mass= masses(5) 
 else
  call fatal_error("ATOM TYPE NOT IN THE LIST")
 endif

 mass = mass*elec_mass

 return
end

