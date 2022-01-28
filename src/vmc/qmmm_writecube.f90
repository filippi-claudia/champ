module qmmm_writecube_mod
contains
	subroutine qmmm_writecube(filename,title,n_atoms,x0,n_x,n_y,n_z, &
                delta,x_atom,id_atom,chrg_atom,pot)
        
        implicit none
        character (len=80) :: filename,title(2)
        integer :: n_x,n_y,n_z
        integer :: n_atoms
        integer :: i,ix,iy
        integer :: id_atom(n_atoms)
        double precision :: x0(3),delta(3),tmp(3)
        double precision :: x_atom(n_atoms,3)
        double precision :: chrg_atom(n_atoms),pot(n_x,n_y,n_z)
      
        tmp(:)=0.d0 
        open (12,file=filename,form='formatted',status='unknown')
        write (12,'(a80)') title(1)
        write (12,'(a80)') title(2)
        write (12,'(i5,3f12.6)') n_atoms,x0
        write (12,'(i5,3f12.6)') n_x,delta(1),tmp(2),tmp(3)
        write (12,'(i5,3f12.6)') n_y,tmp(1),delta(2),tmp(3)
        write (12,'(i5,3f12.6)') n_z,tmp(1),tmp(2),delta(3)
        do i=1,n_atoms
          write(12,'(i5,4f12.6)') id_atom(i),chrg_atom(i),x_atom(i,:)
        enddo
        do ix=1,n_x
          do iy=1,n_y
            write(12,'(6e13.5)') pot(ix,iy,:)
          enddo
        enddo

       write (*,*)
       write(*,*) "Writing cube file ..."
       write (*,'(a15,3i7)') 'Dimensions  :' ,n_x,n_y,n_z
       write (*,'(a10,3f12.6)') 'Mesh  : ',delta

        return
        end
end module
