      module nodes_distance_mod
      contains
      subroutine nodes_distance(v,distance_node,iflag)
c Written by Claudia Filippi

      use system, only: nelec
      use velocity_jastrow, only: vj, vjn
      use precision_kinds, only: dp
      implicit none

      integer :: i, iflag, k
      real(dp) :: distance_node
      real(dp), dimension(3, *) :: v
      real(dp), dimension(3, nelec) :: vdonly










      if(iflag.eq.0) then
        do k=1,3
          do i=1,nelec
            vdonly(k,i)=v(k,i)-vjn(k,i)
          enddo
        enddo
       else
        do k=1,3
          do i=1,nelec
            vdonly(k,i)=v(k,i)-vj(k,i)
          enddo
        enddo
      endif

      distance_node=0.d0
      do k=1,3
        do i=1,nelec
          distance_node=distance_node+vdonly(k,i)*vdonly(k,i)
        enddo
      enddo
      distance_node=1.d0/dsqrt(distance_node)

      return
      end

      function rnorm_nodes_num(distance_node,epsilon)

      use precision_kinds, only: dp
      implicit none


      real(dp) :: distance_node, dum, epsilon, rnorm_nodes_num

      if(distance_node.ge.epsilon) then
        rnorm_nodes_num=distance_node
       else
        dum=distance_node/epsilon
        rnorm_nodes_num=epsilon*dum**dum
      endif

      return
      end
      end module
