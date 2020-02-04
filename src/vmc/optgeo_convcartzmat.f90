  
  ! Omar Valsson (o.valsson@tnw.utwente.nl)
  !
  ! Set of subroutines for converting between 
  ! cartesian coordinanes and internial coordinaets
  ! (Z matrix).


module misc_bond_func
  implicit none
  integer, parameter, private  ::  dp = selected_real_kind(2*precision(1.0))
  real(kind=dp), parameter, private  ::  pi = 3.14159265358979323846_dp
  real(kind=dp), parameter, public  ::  rad2deg = 180.0_dp/pi
  real(kind=dp), parameter, public  ::  deg2rad = pi/180.0_dp

  private
  public distance
  public bond_angle
  public torsion_angle
  public vector_length
  public cross_product

  contains
  
  function distance(cA,cB)
    implicit none
    real(kind=dp)  ::  distance
    real(kind=dp), dimension(1:3), intent(in)  ::  cA,cB
    
    distance = sqrt(dot_product(cA-cB,cA-cB))
  end function distance

  function bond_angle(cA,cB,cC)
    implicit none
    real(kind=dp)  ::  bond_angle
    real(kind=dp), dimension(1:3), intent(in)  ::  cA,cB,cC
    real(kind=dp)  ::  cos_angle
    real(kind=dp), dimension(1:3)  ::  vecBA,vecBC

    vecBA=cA-cB; vecBC=cC-cB
    cos_angle = dot_product(vecBA,vecBC)/(distance(cA,cB)*distance(cB,cC))
    !if((cos_angle + 1.0_dp) .lt. 0.0001) then
    !   bond_angle = 180.0_dp
    !else
    !TODO JF try in radians
    !bond_angle = acos(cos_angle)*rad2deg
    bond_angle = acos(cos_angle)
    !endif
  end function bond_angle

  function torsion_angle(cA,cB,cC,cD)
  ! Taken from the torsion_angle subroutine from bo.c in
  ! Babel 1.6 
    implicit none
    real(kind=dp)  ::  torsion_angle
    real(kind=dp), dimension(1:3), intent(in)  ::  cA,cB,cC,cD
    real(kind=dp), dimension(1:3)  ::  vecBA,vecCB,vecDC,norm1,norm2
    real(kind=dp)  ::  cos_angle

    vecBA=cA-cB; vecCB=cB-cC; vecDC=cC-cD
    norm1 = cross_product(vecCB,vecBA)
    norm2 = cross_product(vecDC,vecCB)   

    cos_angle = dot_product(norm1,norm2)/ &
        & (vector_length(norm1)*vector_length(norm2))
    
    if (cos_angle > 1.0_dp) then
       cos_angle = 1.0_dp
    endif
    if (cos_angle < -1.0_dp) then
       cos_angle = -1.0_dp
    endif
   
    !torsion_angle = acos(cos_angle)*rad2deg
    torsion_angle = acos(cos_angle)
    !
    if (dot_product(vecBA,norm2) < 0.0_dp) then
    !   torsion_angle = 360.0_dp - torsion_angle
      torsion_angle = 2d0 * pi - torsion_angle
    endif

    !if (torsion_angle > 180.0_dp) then
    !   torsion_angle = torsion_angle - 360.0_dp
    if (torsion_angle > pi) then
      torsion_angle = torsion_angle - 2d0 * pi
    endif

  end function torsion_angle

  function cross_product(vA,vB)
    implicit none
    real(kind=dp), dimension(1:3)  ::  cross_product
    real(kind=dp), dimension(1:3), intent(in)  ::  vA,vB

    cross_product(1) = vA(2)*vB(3) - vA(3)*vB(2)
    cross_product(2) = vA(3)*vB(1) - vA(1)*vB(3)
    cross_product(3) = vA(1)*vB(2) - vA(2)*vB(1)
  end function cross_product

  function vector_length(vA)
    implicit none
    real(kind=dp)  ::  vector_length
    real(kind=dp), dimension(1:3), intent(in)  ::  vA

    vector_length = sqrt(dot_product(vA,vA))
  end function vector_length
end module misc_bond_func


subroutine cart2zmat(natoms,coord_cart,conn_zmatrix,coord_zmatrix)
  use misc_bond_func 
  implicit none
  integer, parameter  ::  dp = selected_real_kind(2*precision(1.0))
  integer, intent(in)  ::  natoms
  real(kind=dp), dimension(1:3,1:natoms), intent(in)  ::  coord_cart
  integer, dimension(1:3,1:natoms), intent(in)  ::  conn_zmatrix
  real(kind=dp), dimension(1:3,1:natoms), intent(out)  ::  coord_zmatrix
  integer  ::  nA,nB,nC,nD
  real(kind=dp), dimension(1:3)  ::  cA,cB,cC,cD
  integer  ::  i

  coord_zmatrix = 0.0_dp
  ! Atom 1
  nA = 1
  coord_zmatrix(1,nA) = 0.0_dp
  coord_zmatrix(2,nA) = 0.0_dp
  coord_zmatrix(3,nA) = 0.0_dp
  if (natoms .eq. 1) return
  ! Atom 2
  nA = 2
  nB = conn_zmatrix(1,nA)
  if(nB.ge.nA) then
     write(6,*) 'zmat2cart: error in conn_zmatrix.' ; return
  endif
  if(nB .eq. 0) return
  cA = coord_cart(:,nA)
  cB = coord_cart(:,nB)
  coord_zmatrix(1,nA) = distance(cA,cB)
  coord_zmatrix(2,nA) = 0.0_dp
  coord_zmatrix(3,nA) = 0.0_dp
  if (natoms .eq. 2) return
  ! Atom 3
  nA = 3
  nB = conn_zmatrix(1,nA)
  nC = conn_zmatrix(2,nA)
  if(nB.ge.nA .or. nC.ge.nA) then
     write(6,*) 'zmat2cart: error in conn_zmatrix.' ; return
  endif
  if(nB .eq. 0 .or. nC .eq. 0) return
  cA = coord_cart(:,nA)
  cB = coord_cart(:,nB)
  cC = coord_cart(:,nC)
  coord_zmatrix(1,nA) = distance(cA,cB)
  coord_zmatrix(2,nA) = bond_angle(cA,cB,cC)
  coord_zmatrix(3,nA) = 0.0_dp
  if (natoms .eq. 3) return
  ! Atoms  = > 4
  do i = 4,natoms
     nA = i
     nB = conn_zmatrix(1,nA)
     nC = conn_zmatrix(2,nA)
     nD = conn_zmatrix(3,nA)
     if(nB.ge.nA .or. nC.ge.nA .or. nD.ge.nA) then
        write(6,*) 'zmat2cart: error in conn_zmatrix.' ; return
     endif
     if (nB.eq.0 .or. nC.eq.0 .or. nD.eq.0) return
     cA = coord_cart(:,nA)
     cB = coord_cart(:,nB)
     cC = coord_cart(:,nC)
     cD = coord_cart(:,nD)
     coord_zmatrix(1,nA) = distance(cA,cB)
     coord_zmatrix(2,nA) = bond_angle(cA,cB,cC)
     coord_zmatrix(3,nA) = torsion_angle(cA,cB,cC,cD)
  enddo    
end subroutine cart2zmat


subroutine zmat2cart(natoms,conn_zmatrix,coord_zmatrix,coord_cart)
  ! this is taken from int_to_cart subroutine in intcart.c 
  ! in Babel 1.6
  use misc_bond_func 
  implicit none
  integer, parameter  ::  dp = selected_real_kind(2*precision(1.0))
  integer, intent(in)  ::  natoms
  integer, dimension(1:3,1:natoms), intent(in)  ::  conn_zmatrix
  real(kind=dp), dimension(1:3,1:natoms), intent(in)  ::  coord_zmatrix
  real(kind=dp), dimension(1:3,1:natoms), intent(out)  ::  coord_cart
  integer  ::  nA,nB,nC,nD
  real(kind=dp), dimension(1:3)  ::  cA,cB,cC,cD
  real(kind=dp)  ::  dist,angle,dihed
  real(kind=dp)  ::  cosa,sina,cosd,sind
  real(kind=dp)  ::  xa,ya,za,xb,yb,zb,xd,yd,zd
  real(kind=dp)  ::  xpa,ypa,zqa,xpd,ypd,zpd,xqd,yqd,zqd
  real(kind=dp)  ::  cosph,sinph,costh,sinth,coskh,sinkh
  real(kind=dp)  ::  rbc,xyb,yza,tmp1_dp
  logical  ::  rotate_yaxis
  integer  ::  i

  coord_cart = 0.0_dp
  ! Atom 1
  nA = 1
  coord_cart(1,nA) = 0.0_dp
  coord_cart(2,nA) = 0.0_dp
  coord_cart(3,nA) = 0.0_dp
  if (natoms .eq. 1) return
  ! Atom 2
  nA = 2
  nB = conn_zmatrix(1,nA)
  if(nB.ge.nA) then
     write(6,*) 'zmat2cart: error in conn_zmatrix.' ; return
  endif
  if(nB .eq. 0) return
  coord_cart(1,nA) = coord_zmatrix(1,nA)
  coord_cart(2,nA) = 0.0_dp
  coord_cart(3,nA) = 0.0_dp
  if (natoms .eq. 2) return
  ! Atom 3
  nA = 3
  nB = conn_zmatrix(1,nA)
  nC = conn_zmatrix(2,nA)
  if(nB.ge.nA .or. nC.ge.nA) then
     write(6,*) 'zmat2cart: error in conn_zmatrix.' ; return
  endif
  if(nB .eq. 0 .or. nC .eq. 0) return
  dist = coord_zmatrix(1,nA)
  angle = coord_zmatrix(2,nA)*deg2rad
  cosa = cos(angle)
  sina = sin(angle)
  coord_cart(1,nA) = coord_cart(1,nB) - (-1)**nB*dist*cosa
  coord_cart(2,nA) = dist*sina
  coord_cart(3,nA) = 0.0_dp
  if (natoms .eq. 3) return
  ! Atoms  = > 4
  do i = 4,natoms
     nA = i
     nB = conn_zmatrix(1,nA)
     nC = conn_zmatrix(2,nA)
     nD = conn_zmatrix(3,nA)
     if(nB.ge.nA .or. nC.ge.nA .or. nD.ge.nA) then
        write(6,*) 'zmat2cart: error in conn_zmatrix.' ; return
     endif
     if (nB.eq.0 .or. nC.eq.0 .or. nD.eq.0) return
     dist = coord_zmatrix(1,nA)
     angle = coord_zmatrix(2,nA)*deg2rad
     dihed = coord_zmatrix(3,nA)*deg2rad
     cB = coord_cart(:,nB)
     cC = coord_cart(:,nC)
     cD = coord_cart(:,nD)
     
     xb = cC(1) - cB(1)
     yb = cC(2) - cB(2)
     zb = cC(3) - cB(3)
     
     rbc = xb**2 + yb**2 + zb**2
     if(rbc.lt.0.0001_dp) then
        write(6,*) 'zmat2cart: error rbc.lt.0001.'; return
     endif
     rbc = 1.0_dp/sqrt(rbc)

     cosa = cos(angle)
     sina = sin(angle)

     if(abs(cosa).ge.0.999999) then
        ! Colinear
        coord_cart(1,nA) = coord_cart(1,nB)+dist*rbc*cosa*xb
        coord_cart(2,nA) = coord_cart(2,nB)+dist*rbc*cosa*yb
        coord_cart(3,nA) = coord_cart(3,nB)+dist*rbc*cosa*zb
     else
        xa = cD(1) - cB(1)
        ya = cD(2) - cB(2)
        za = cD(3) - cB(3)

        sind = -sin(dihed)
        cosd = cos(dihed)
        xd = dist*cosa
        yd = dist*sina*cosd
        zd = dist*sina*sind
        
        xyb = sqrt(xb**2 + yb**2)
        if(xyb .lt. 0.1) then
           ! Rotate about the y-axis
           tmp1_dp=za; za=-xa; xa=tmp1_dp
           tmp1_dp=zb; zb=-xb; xb=tmp1_dp
           xyb = sqrt(xb**2 + yb**2)
           rotate_yaxis = .true.
        else
           rotate_yaxis = .false.
        endif
        
        costh = xb/xyb
        sinth = yb/xyb
        xpa = costh*xa + sinth*ya
        ypa = costh*ya - sinth*xa
        
        sinph = zb*rbc    
        cosph = sqrt(1.0_dp - sinph**2)
        zqa = cosph*za - sinph*xpa
       
        yza = sqrt(ypa**2 + zqa**2)
        if(yza .gt. 1.0e-10_dp) then
           coskh = ypa/yza
           sinkh = zqa/yza
           
           ypd = coskh*yd - sinkh*zd
           zpd = coskh*zd + sinkh*yd
        else
           ! coskh = 1.0
           ! sinkh =  0.0
           ypd = yd
           zpd = zd
        endif
       
        xpd = cosph*xd  - sinph*zpd
        zqd = cosph*zpd + sinph*xd
        xqd = costh*xpd - sinth*ypd
        yqd = costh*ypd + sinth*xpd
 
        if(rotate_yaxis) then
           coord_cart(1,nA) = coord_cart(1,nB) - zqd
           coord_cart(2,nA) = coord_cart(2,nB) + yqd
           coord_cart(3,nA) = coord_cart(3,nB) + xqd
        else
           coord_cart(1,nA) = coord_cart(1,nB) + xqd
           coord_cart(2,nA) = coord_cart(2,nB) + yqd
           coord_cart(3,nA) = coord_cart(3,nB) + zqd
        endif
        
     endif  
  enddo    
end subroutine zmat2cart


subroutine zmat2cart_rc(natoms,conn_zmatrix,coord_zmatrix,coord_cart,refcoord_cart)
  use misc_bond_func 
  implicit none
  integer, parameter  ::  dp = selected_real_kind(2*precision(1.0))
  integer, intent(in)  ::  natoms
  integer, dimension(1:3,1:natoms), intent(in)  ::  conn_zmatrix
  real(kind=dp), dimension(1:3,1:natoms), intent(in)  ::  coord_zmatrix
  real(kind=dp), dimension(1:3,1:natoms), intent(out)  ::  coord_cart
  real(kind=dp), dimension(1:3,1:3), intent(in)  ::  refcoord_cart
  integer  ::  nA,nB,nC,nD
  real(kind=dp), dimension(1:3)  ::  cA,cB,cC,cD
  real(kind=dp), dimension(1:3)  ::  vBC,vBA,v1,v2
  real(kind=dp)  ::  dist,angle,dihed
  real(kind=dp)  ::  cosa,sina,cosd,sind
  real(kind=dp)  ::  xa,ya,za,xb,yb,zb,xd,yd,zd
  real(kind=dp)  ::  xpa,ypa,zqa,xpd,ypd,zpd,xqd,yqd,zqd
  real(kind=dp)  ::  cosph,sinph,costh,sinth,coskh,sinkh
  real(kind=dp)  ::  rbc,xyb,yza,tmp1_dp
  logical  ::  rotate_yaxis
  integer  ::  i

  coord_cart = 0.0_dp
  ! Atom 1
  nA = 1
  coord_cart(1,nA) = refcoord_cart(1,nA)
  coord_cart(2,nA) = refcoord_cart(2,nA)
  coord_cart(3,nA) = refcoord_cart(3,nA)
  if (natoms .eq. 1) return
  ! Atom 2
  nA = 2
  nB = conn_zmatrix(1,nA)
  if(nB.ge.nA) then
     write(6,*) 'zmat2cart: error in conn_zmatrix.' ; return
  endif
  if (nB.eq.0) return
  dist = coord_zmatrix(1,nA)
  vBA = refcoord_cart(:,nA)-coord_cart(:,nB)
  vBA = vBA/vector_length(vBA)
  coord_cart(1,nA) = coord_cart(1,nB) + dist*vBA(1)
  coord_cart(2,nA) = coord_cart(2,nB) + dist*vBA(2)
  coord_cart(3,nA) = coord_cart(3,nB) + dist*vBA(3)
  if (natoms .eq. 2) return
  ! Atom 3
  nA = 3
  nB = conn_zmatrix(1,nA)
  nC = conn_zmatrix(2,nA)
  if(nB.ge.nA .or. nC.ge.nA) then
     write(6,*) 'zmat2cart: error in conn_zmatrix.' ; return
  endif
  if (nB.eq.0 .or. nC.eq.0) return
  dist = coord_zmatrix(1,nA)
  angle = coord_zmatrix(2,nA)*deg2rad
  cosa = cos(angle)
  sina = sin(angle)
  vBA = refcoord_cart(:,nA)-coord_cart(:,nB)
  vBC = coord_cart(:,nC)-coord_cart(:,nB)
  v1 = cross_product(vBA,vBC)
  v2 = cross_product(vBC,v1) 
  vBC = vBC/vector_length(vBC)
  v2 = v2/vector_length(v2)
  coord_cart(1,nA) = coord_cart(1,nB) + dist*cosa*vBC(1) + dist*sina*v2(1)
  coord_cart(2,nA) = coord_cart(2,nB) + dist*cosa*vBC(2) + dist*sina*v2(2)
  coord_cart(3,nA) = coord_cart(3,nB) + dist*cosa*vBC(3) + dist*sina*v2(3)
  if (natoms .eq. 3) return
  ! Atoms  = > 4
  do i = 4,natoms
     nA = i
     nB = conn_zmatrix(1,nA)
     nC = conn_zmatrix(2,nA)
     nD = conn_zmatrix(3,nA)
     if(nB.ge.nA .or. nC.ge.nA .or. nD.ge.nA) then
        write(6,*) 'zmat2cart: error in conn_zmatrix.' ; return
     endif
     if (nB.eq.0 .or. nC.eq.0 .or. nD.eq.0) return
     dist = coord_zmatrix(1,nA)
     angle = coord_zmatrix(2,nA)*deg2rad
     dihed = coord_zmatrix(3,nA)*deg2rad
     cB = coord_cart(:,nB)
     cC = coord_cart(:,nC)
     cD = coord_cart(:,nD)
     
     xb = cC(1) - cB(1)
     yb = cC(2) - cB(2)
     zb = cC(3) - cB(3)
     
     rbc = xb**2 + yb**2 + zb**2
     if(rbc.lt.0.0001_dp) then
        write(6,*) 'zmat2cart: error rbc.lt.0001.'; return
     endif
     rbc = 1.0_dp/sqrt(rbc)

     cosa = cos(angle)
     sina = sin(angle)

     if(abs(cosa).ge.0.999999) then
        ! Colinear
        coord_cart(1,nA) = coord_cart(1,nB)+dist*rbc*cosa*xb
        coord_cart(2,nA) = coord_cart(2,nB)+dist*rbc*cosa*yb
        coord_cart(3,nA) = coord_cart(3,nB)+dist*rbc*cosa*zb
     else
        xa = cD(1) - cB(1)
        ya = cD(2) - cB(2)
        za = cD(3) - cB(3)

        sind = -sin(dihed)
        cosd = cos(dihed)
        xd = dist*cosa
        yd = dist*sina*cosd
        zd = dist*sina*sind
        
        xyb = sqrt(xb**2 + yb**2)
        if(xyb .lt. 0.1) then
           ! Rotate about the y-axis
           tmp1_dp=za; za=-xa; xa=tmp1_dp
           tmp1_dp=zb; zb=-xb; xb=tmp1_dp
           xyb = sqrt(xb**2 + yb**2)
           rotate_yaxis = .true.
        else
           rotate_yaxis = .false.
        endif
        
        costh = xb/xyb
        sinth = yb/xyb
        xpa = costh*xa + sinth*ya
        ypa = costh*ya - sinth*xa
        
        sinph = zb*rbc    
        cosph = sqrt(1.0_dp - sinph**2)
        zqa = cosph*za - sinph*xpa
       
        yza = sqrt(ypa**2 + zqa**2)
        if(yza .gt. 1.0e-10_dp) then
           coskh = ypa/yza
           sinkh = zqa/yza
           
           ypd = coskh*yd - sinkh*zd
           zpd = coskh*zd + sinkh*yd
        else
           ! coskh = 1.0
           ! sinkh =  0.0
           ypd = yd
           zpd = zd
        endif
       
        xpd = cosph*xd  - sinph*zpd
        zqd = cosph*zpd + sinph*xd
        xqd = costh*xpd - sinth*ypd
        yqd = costh*ypd + sinth*xpd
 
        if(rotate_yaxis) then
           coord_cart(1,nA) = coord_cart(1,nB) - zqd
           coord_cart(2,nA) = coord_cart(2,nB) + yqd
           coord_cart(3,nA) = coord_cart(3,nB) + xqd
        else
           coord_cart(1,nA) = coord_cart(1,nB) + xqd
           coord_cart(2,nA) = coord_cart(2,nB) + yqd
           coord_cart(3,nA) = coord_cart(3,nB) + zqd
        endif
        
     endif  
  enddo    
end subroutine zmat2cart_rc


