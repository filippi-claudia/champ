      subroutine nodes_distance(v,distance_node,iflag)
c Written by Claudia Filippi

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use const, only: nelec
      use velocity_jastrow, only: vj, vjn
      implicit real*8(a-h,o-z)






      parameter(one=1.d0)



      dimension v(3,*),vdonly(3,MELEC)

      if(iflag.eq.0) then
        do 10 k=1,3
          do 10 i=1,nelec
  10        vdonly(k,i)=v(k,i)-vjn(k,i)
       else
        do 20 k=1,3
          do 20 i=1,nelec
  20        vdonly(k,i)=v(k,i)-vj(k,i)
      endif

      distance_node=0.d0
      do 100 k=1,3
        do 100 i=1,nelec
  100     distance_node=distance_node+vdonly(k,i)*vdonly(k,i)
      distance_node=1.d0/dsqrt(distance_node)

      return
      end

      function rnorm_nodes_num(distance_node,epsilon)

      implicit real*8(a-h,o-z)

      if(distance_node.ge.epsilon) then
        rnorm_nodes_num=distance_node
       else
        dum=distance_node/epsilon
        rnorm_nodes_num=epsilon*dum**dum
      endif

      return
      end
