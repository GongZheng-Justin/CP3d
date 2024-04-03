  SumXDir=0.0_RK  ! ux grid
  do k=0,nDistribute
    kd = k+idzc_interp
    prz= RatioZc(k)
    do j=0,nDistribute
      jd  = j+idyc_interp
      pry = RatioYc(j)*prz
      do i=0,nDistribute
        id = i+idxp_interp
        prx= RatioXp(i)
        SumXDir= SumXDir + ux(id,jd,kd)*prx*pry
      enddo
    enddo
  enddo
      
  SumYDir=0.0_RK  ! uy gird
  do k=0,nDistribute
    kd = k+idzc_interp
    prz= RatioZc(k)
    do j=0,nDistribute
      jd = j+idyp_interp
      pry= RatioYp(j)*prz
      do i=0,nDistribute
        id = i+idxc_interp
        prx= RatioXc(i)
        SumYDir= SumYDir + uy(id,jd,kd)*prx*pry
      enddo
    enddo
  enddo
      
  SumZDir=0.0_RK  ! uz grid
  do k=0,nDistribute
    kd = k+idzp_interp
    prz= RatioZp(k)
    do j=0,nDistribute
      jd = j+idyc_interp
      pry= RatioYc(j)*prz
      do i=0,nDistribute
        id = i+idxc_interp
        prx= RatioXc(i)
        SumZDir= SumZDir + uz(id,jd,kd)*prx*pry
      enddo
    enddo
  enddo
