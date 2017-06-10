      Dimension Te(350,350),W(350,250),Tep(350,250),Wp(350,350)
      Dimension x(350),y(350),HTI(350,350),FFI(350,350)	
	real W,Wp,Te,Tep
c========= Geometry Details ==================
	t=.1
	s=.5
	c=.5
	h2=1.

	n1=65
	nc=25
	nt=1
	ns=30
c========= Phiscall Details ==================
      aomega=25.
c	read(*,*)aomega
      rksf=2./t*aomega
	qTe =7.
 	qW =7.
	g1=.9

c=========  R Function Iteraion  ==================
      rx=1.1
	ry1=1.1
	ryc=1.1
	dx=t/2./nt/10
	dy=dx
	do i=1,4000
	rx=rx+.4*(-rx+(1.-(s-t)/2./dx*(1-rx))**(1./(ns)))
	ryc=ryc+.4*(-ryc+(1.-c/2./dy*(1-ryc))**(1./nc))
	ry1=ry1+.4*(-ry1+(1.-1./2./dy*(1-ry1))**(1./n1))
      end do
	write(*,*)" rx , ry ==>>",rx,ry1,ryc 
c=========  Grid Generation  ==================
	dxo=dx
	dyo=dy
	dxoo=dx
	dyoo=dy
	write(*,*)dyo,dyoo
	y(1)=0
      do j=2,n1+1
	y(j)=y(j-1)+dy
	dy=ry1*dy
      end do

	jj=0
	do j=n1+1-1,1,-1
	jj=jj+1
	y(n1+1+jj)=y(n1+1)+(y(n1+1)-y(j))
	end do

	write(*,*)dyo,dyoo
	dy=dyo
	do j=2*n1+1+1,2*n1+1+nc
	y(j)=y(j-1)+dy
	dy=ryc*dy
	end do

	jj=0
	do j=2*n1+1+nc-1,2*n1+1,-1
 	jj=jj+1
      y(2*n1+1+nc+jj)=y(2*n1+1+nc)+(y(2*n1+1+nc)-y(j))
	end do

	ny=2*n1+1+nc+jj
c====================
	dx=t/2./nt 	
      do i=1,nt+1
      x(i)=(i-1)*dx
	end do

	dx=dxo
      do i=nt+2,nt+1+ns
	x(i)=x(i-1)+dx
 	dx=rx*dx
	end do
	
      ii=0
	do i=nt+1+ns-1,1,-1
	ii=ii+1
      x(nt+1+ns+ii)=x(nt+1+ns)+(x(nt+1+ns)-x(i))
	end do
	
	nx=nt+1+ns+ii
c===================
	write(*,*)" nx , ny ==>>",nx,ny

	ic1=1+nt
	ic2=nx-nt
	jc2=2*n1*h2+1
	jc1=2*n1+1
	write(*,*)"Ic1_Ic2_Jc1_Jc2 ",ic1,ic2,jc1,jc2
c	ny=jc1
c	c=0.
c=======================================
	W=0.
	Wp=0.
	Te=0.5
	Tep=0.
c================ Intital ======================
	
	
	do i=2,ny-1
	do j=2,nx-1
	Tep(i,j)=-i*y(j)/1000.
	end do
	end do
	q=1
	qq=0

      i=1
	wbar=1
      iteration=0

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$     WWWWW         WWWWWW       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
1	iteration=iteration + 1
	if (iteration/3000..eq.int(iteration/3000))then
	write(*,20)iteration/1000,ii,jj,Wdif,wbardif,W(ii,jj),fre
	end if
c=========== W ========= W ==================
c=========== W ========= W ==================
	do j = 2, ny-1
	do i = 2, nx-1
	
      if ((i.le.ic1).and.(j.le.jc1))goto 2
	if ((i.ge.ic2).and.(j.le.jc2))goto 2 
      
	dxi=x(i)-x(i-1)
	dxii=x(i+1)-x(i)
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)

      b1=(Wp(i+1,j)*dxi+Wp(i-1,j)*dxii)/(dxi*dxii*.5*(dxi+dxii))
 	b2=(Wp(i,j+1)*dyj+Wp(i,j-1)*dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	b3=+1
	
	a1=-(dxi+dxii)/(dxi*dxii*.5*(dxi+dxii))
	a2=-(dyj+dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	a3=0.
	
     	W(i,j)=Wp(i,j)+g1*(-(b1+b2+b3)/(a1+a2+a3)-Wp(i,j))

2     end do
	end do
c===========  B.C.  W ========================
c==========Right SYMETRI
	i=nx
      do j=1,ny
	W(i,j)=W(i-1,j)
	end do
c==========Left SYMETRI
	i=1
      do j=1,ny
	W(i,j)=W(i+1,j)
	end do
c==========Shroude
 	j=ny
      do i=1,nx
	W(i,j)=0.
	end do
c=========== W bar ========= W bar ==================
	   	wbarnew=.0
      do j = 1, ny-1
   	do i = 1, nx-1
	
c      if ((i.le.ic1).and.(j.le.jc1))goto 3
c	if ((i.ge.ic2).and.(j.le.jc2))goto 3 
      
	dx=x(i+1)-x(i)
	dy=y(j+1)-y(j)
		
     	wbarnew=wbarnew+.25*(W(i,j)+W(i+1,j)+W(i,j+1)+W(i+1,j+1))*dx*dy

3     end do
	end do
	wbarnew=wbarnew/((1+c)*s)
	wbardif=abs((wbarnew-wbar)/wbar)/g1
c=========== fRE ========= fRE ==================
	i=ic1
      BBleft=0
	do j=1,jc1-1
	tx=(W(i+1,j)-W(i,j))/(x(i+1)-x(i))
	BBleft=BBleft+tx*(y(j+1)-y(j))
     	end do
c=============
	BBbase=0.
	j=1
      do i=ic1,ic2
	ty=(W(i,j+1)-W(i,j))/(y(j+1)-Y(j))
   	BBbase=BBbase+ty*(x(i+1)-x(i-1))*.5
  	end do
c=============
	BBright=0
	i=ic2
	do j=1,jc2-1
	tx=(W(i,j)-W(i-1,j))/(x(i-1)-x(i))
	BBright=BBright+tx*(y(j+1)-y(j))
      end do
c=============
	BBtop=0.
	j=ny
      do i=ic1,ic2
	ty=(W(i,j-1)-W(i,j))/(y(j)-Y(j-1))
   	BBtop=BBtop+ty*(x(i+1)-x(i-1))*.5
  	end do

	BBAA=abs(BBtop)+abs(BBright)+abs(BBbase)+abs(BBleft)
	fre=32.*(1.+C)*(S)*BBAA/wbarnew

c ============= ?? RE START  ?? ===================
	do i = 1,nx
     	do j = 1,ny

	if (abs(W(i,j)).ge..1e-9)then
	Wdif=Abs((W(i,j)-Wp(i,j))/W(i,j))/g1
	else
	Wdif=0.
	end if
	ii=i
	jj=j

	If (Wdif-10.**(-qW) .gt.0.0)goto 9

	end do
	end do
	If (Wbardif-10.**(-qW) .gt.0.0)goto 9

   	GoTo 99

c=========  Repalace Pervious =========================
9	do i = 1,nx
	do j = 1,ny
	Wp(i,j)=W(i,j)
	end do
	end do
	wbar=wbarnew
	goto 1

c============Replace W
99      do i = 1,nx
	do j = 1,ny
	W(i,j)=W(i,j)/wbarnew
	end do
	end do
	write(*,*)"Calculating W Finished   W_bar=",wbar,fre
	iteration=0
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$     Te    Te    Te    Te       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	Tab=1
 4	iteration=iteration + 1
      if (iteration/5000..eq.int(iteration/5000))then
	write(*,200)iteration/1000,ii,jj,Tedif,Tabdif,Tab,Te(ii,jj)
	end if
	if (iteration.gt.300000)g1=.5
c========== Te ========== Te ==================
c========== Te ========== Te ==================
	do j = 2, ny-1
	do i = 2, nx-1
	
	dxi=x(i)-x(i-1)
	dxii=x(i+1)-x(i)
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)
	
      b1=(Tep(i+1,j)*dxi+Tep(i-1,j)*dxii)/(dxi*dxii*.5*(dxi+dxii))
 	b2=(Tep(i,j+1)*dyj+Tep(i,j-1)*dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	b3=-W(i,j)/((1+c)*s)
	
	a1=-(dxi+dxii)/(dxi*dxii*.5*(dxi+dxii))
	a2=-(dyj+dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	a3=0.
      if(((i.le.ic1).and.(j.le.jc1)).or.((i.ge.ic2).and.(j.le.jc2)))then
	Te(i,j)=Tep(i,j)+g1*(-(b1+b2)/(a1+a2)-Tep(i,j)) 
	else
	Te(i,j)=Tep(i,j)+g1*(-(b1+b2+b3)/(a1+a2+a3)-Tep(i,j))
	end if
5     end do
	end do

c===========  B.C.  Te  ========================
c==========Left FIN
	i=ic1
	dx=x(ic1+1)-x(ic1)
      do j=2,jc1-1
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)
	b=aomega/(dyj*dyjj*.5*(dyj+dyjj))
	b1=b*Tep(i,j+1)*dyj+b*Tep(i,j-1)*dyjj+Tep(i+1,j)/dx
	b2=b*(dyj+dyjj)+1/dx
	Te(i,j)=Tep(i,j)+(b1/b2-Tep(i,j))
	end do
c==========Right FIN
	i=ic2
	dx=x(ic2)-x(ic2-1)
      do j=2,jc2-1
 	dyj=y(j)-y(j-1)
 	dyjj=y(j+1)-y(j)
	b=aomega/(dyj*dyjj*.5*(dyj+dyjj))
	b1=b*Tep(i,j+1)*dyj+b*Tep(i,j-1)*dyjj+Tep(i-1,j)/dx
	b2=b*(dyj+dyjj)+1/dx
	Te(i,j)=Tep(i,j)+(b1/b2-Tep(i,j))
	end do
c==========Right SYMETRI
	i=nx
      do j=1,ny
	Te(i,j)=Te(i-1,j)
	end do
c==========Left SYMETRI
	i=1
      do j=1,ny
	Te(i,j)=Te(i+1,j)
	end do

c==========Top Left Wall
	j=jc1
      do i=1,ic1
	Te(i,j)=Te(i,j-1)
c	Te(i,j)=(Te(i,j-2)-6*Te(i,j-1)+2*Te(i,j+1))/(+3)

	end do
c==========Top Right Wall
	j=jc2
      do i=ic2,nx
	Te(i,j)=Te(i,j-1)
c	Te(i,j)=(Te(i,j-2)-6*Te(i,j-1)+2*Te(i,j+1))/(+3)
      end do
c==========SHROUD
	j=ny
      do i=1,nx
	Te(i,j)=Te(i,j-1)
	end do
c==========BASE
	j=1
      do i=1,nx
	Te(i,j)=0.
	end do

	Te(ic1,jc1)=(Te(ic1,jc1-1)+Te(ic1-1,jc1))/2.
	Te(ic2,jc2)=(Te(ic2,jc2-1)+Te(ic2+1,jc2))/2.
c=========Tb
   	Tabnew=.0
      do j = 1, ny-1
 	do i = 1, nx-1
	
c      if ((i.le.ic1).and.(j.le.jc1))goto 33
c	if ((i.ge.ic2).and.(j.le.jc2))goto 33 
      
	dx=x(i+1)-x(i)
	dy=y(j+1)-y(j)
	b1=Te(i,j)*W(i,j)+Te(i+1,j)*W(i+1,j)
	b2=Te(i,j+1)*W(i,j+1)+Te(i+1,j+1)*W(i+1,j+1)	
     	Tabnew=Tabnew+.25*(b1+b2)*dx*dy

33    end do
	end do
	Tabnew=Tabnew/((1+c)*s)
	Tabdif=abs((Tabnew-Tab)/Tabnew)
c ============= ?? RE START  ?? ===================
      do i = 1,nx
      do j = 1,ny

	if (abs(Te(i,j)).ge.1e-9) then
      Tedif=Abs((Te(i,j)-Tep(i,j))/Te(i,j))/g1
	else
	Tedif=0.
	end if

	ii=i
	jj=j
	
	If (Tedif-10.**(-qTe) .gt.0.0 )goto 8

	end do
	end do
	If (Tabdif-10.**(-qTe) .gt.0.0 )goto 8

   	GoTo 88
c=========  Repalace Pervious =========================
8	do i = 1,nx
	do j = 1,ny
	Tep(i,j)=Te(i,j)
	end do
	end do
	Tab=Tabnew
	goto 4
c============Replace Fi an Te
88    do i = 1,nx
	do j = 1,ny
	Te(i,j)=Te(i,j)
	end do
	end do
	write(*,*)"Calculating Te Finished"
c===========================================================================
c=============== ENTROPY GENERATION ========================================
c=============== ENTROPY GENERATION ========================================
c===========================================================================

      do i = 2,nx-1
 	do j = 2,ny-1

	dxi=x(i)-x(i-1)
 	dxii=x(i+1)-x(i)
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)

      if ((i.le.ic1).and.(j.le.jc1))goto 11
	if ((i.ge.ic2).and.(j.le.jc2))goto 11 

      tx=(Te(i+1,j)*dxi**2-Te(i-1,j)*dxii**2+Te(i,j)*(dxii**2-dxii**2))
	//(dxi*dxii*(dxi+dxii))
	ty=(Te(i,j+1)*dyj**2-Te(i,j-1)*dyjj**2+Te(i,j)*(dyjj**2-dyjj**2))
	//(dyj*dyjj*(dyj+dyjj))

	wx=( w(i+1,j)*dxi**2- w(i-1,j)*dxii**2+ w(i,j)*(dxii**2-dxii**2))
	//(dxi*dxii*(dxi+dxii))
	wy=( w(i,j+1)*dyj**2- w(i,j-1)*dyjj**2+ w(i,j)*(dyjj**2-dyjj**2))
	//(dyj*dyjj*(dyj+dyjj))

           
	HTI(i,j)=tx**2+ty**2
      FFI(i,j)=(2.*(wx**2)+(wy)**2)

11 	end do
	end do
c========== i=ic1 and i=1
	do j=2,ny-1
      i=ic1
	if (j.gt.jc1)i=1
	dxii=x(i+1)-x(i)
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)
     
      tx=(Te(i+1,j)-Te(i,j))/(dxii)
	ty=(Te(i,j+1)*dyj**2-Te(i,j-1)*dyjj**2+Te(i,j)*(dyjj**2-dyjj**2))
	//(dyj*dyjj*(dyj+dyjj))

	wx=( w(i+1,j)- w(i,j))/(dxii)
	wy=( w(i,j+1)*dyj**2- w(i,j-1)*dyjj**2+ w(i,j)*(dyjj**2-dyjj**2))
	//(dyj*dyjj*(dyj+dyjj))
	
      HTI(i,j)=tx**2+ty**2
      FFI(i,j)=(2.*(wx**2)+(wy)**2)

      end do
c========== i=ic2 and i=nx
	do j=2,ny-1
      i=ic2
	if (j.gt.jc2)i=nx

	dxi=x(i)-x(i-1)
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)
     
      tx=(Te(i,j)-Te(i-1,j))/(dxi)
	ty=(Te(i,j+1)*dyj**2-Te(i,j-1)*dyjj**2+Te(i,j)*(dyjj**2-dyjj**2))
	//(dyj*dyjj*(dyj+dyjj))

	wx=( w(i,j)- w(i-1,j))/(dxi)
	wy=( w(i,j+1)*dyj**2- w(i,j-1)*dyjj**2+ w(i,j)*(dyjj**2-dyjj**2))
	//(dyj*dyjj*(dyj+dyjj))
	
      HTI(i,j)=tx**2+ty**2
      FFI(i,j)=(2.*(wx**2)+(wy)**2)

      end do

c========== j=1
	do i=ic1,ic2
	j=1
	dxi=x(i)-x(i-1)
  	dxii=x(i+1)-x(i)
	dyjj=y(j+1)-y(j)
      tx=(Te(i+1,j)*dxi**2-Te(i-1,j)*dxii**2+Te(i,j)*(dxii**2-dxii**2))
	//(dxi*dxii*(dxi+dxii))
	ty=(Te(i,j+1)-Te(i,j))/(dyjj)

	wx=( w(i+1,j)*dxi**2- w(i-1,j)*dxii**2+ w(i,j)*(dxii**2-dxii**2))
	//(dxi*dxii*(dxi+dxii))
	wy=(w(i,j+1)-w(i,j))/(dyjj)

           
	HTI(i,j)=tx**2+ty**2
      FFI(i,j)=(2.*(wx**2)+(wy)**2)
     
      end do

	do i=2,nx-1
	j=ny

	dxi=x(i)-x(i-1)
  	dxii=x(i+1)-x(i)
	dyj=y(j)-y(j-1)
      tx=(Te(i+1,j)*dxi**2-Te(i-1,j)*dxii**2+Te(i,j)*(dxii**2-dxii**2))
	//(dxi*dxii*(dxi+dxii))
	ty=(Te(i,j)-Te(i,j-1))/(dyj)

	wx=( w(i+1,j)*dxi**2- w(i-1,j)*dxii**2+ w(i,j)*(dxii**2-dxii**2))
	//(dxi*dxii*(dxi+dxii))
	wy=(w(i,j)-w(i,j-1))/(dyj)

           
	HTI(i,j)=tx**2+ty**2
      FFI(i,j)=(2.*(wx**2)+(wy)**2)
     
      end do
c=============  COORNER
	FFI(ic1,1)=(FFI(ic1+1,1)+FFI(ic1,1+1))/2
	FFI(1,ny)=(FFI(1+1,ny)+FFI(1,ny-1))/2
	FFI(ic2,1)=(FFI(ic2-1,1)+FFI(ic2-1,1+1))/2
	FFI(nx,ny)=(FFI(nx-1,ny)+FFI(nx,ny+1))/2

	HTI(ic1,1)=(HTI(ic1+1,1)+HTI(ic1,1+1))/2
	HTI(1,ny)=(HTI(1+1,ny)+HTI(1,ny-1))/2
	HTI(ic2,1)=(HTI(ic2-1,1)+HTI(ic2-1,1+1))/2
	HTI(nx,ny)=(HTI(nx-1,ny)+HTI(nx,ny+1))/2

	FFI(ic1,jc1)=(FFI(ic1-1,jc1)+FFI(ic1,jc1-1))/2
	HTI(ic1,jc1)=(HTI(ic1-1,jc1)+HTI(ic1,jc1-1))/2
	
      FFI(ic2,jc2)=(FFI(ic2+1,jc2)+FFI(ic2,jc2-1))/2
	HTI(ic2,jc2)=(HTI(ic2+1,jc2)+HTI(ic2,jc2-1))/2
c============= Post Process
	SumFFI=0.
 	SumHTI=0.
 	SumEG=0.
      do i=1,nx-1
	do j=1,ny-1
	dx=x(i+1)-x(i)
	dy=y(j+1)-y(j)
	bf=(FFI(i,j)+FFI(i+1,j)+FFI(i,j+1)+FFI(i+1,j+1))
	bh=(HTI(i,j)+HTI(i+1,j)+HTI(i,j+1)+HTI(i+1,j+1))
	SumFFI=SumFFI+.25*bf*dx*dy
	SumHTI=SumHTI+.25*bh*dx*dy
     	end do
	end do
	SumEG=SumFFI+SumHTI


c=========================================================================
c===============  NUSSELT  NUSSELT =======================================	
c===============  NUSSELT  NUSSELT =======================================	
c=========================================================================
	open(20,file="20 Nusselt.txt") 
	open(21,file="21 Teta_Fin.txt") 
c===============left  Fin
       write(20,*)"==left  Fin==" 
	i =ic1
      AAleft=0
	do j=2,jc1
      tx=-(Te(i+1,j)-Te(i,j))/(x(i+1)-x(i))
	b1=tx*(1.)/(Te(i,j)-Tab)
      AAleft=AAleft+b1*(y(j+1)-y(j-1))*.5
     	write(20,*)y(j),b1
 	end do
	AAleft=AAleft/1.
c===============Base
 	write(20,*)"==Base==" 
      j =1
      AAbase1=0
	AAbase2=0
 	do i=ic1+1,ic2-1
      ty1=(Te(i,j+1)-Te(i,j))/(y(j+1)-y(j))
  	dyj=y(j+1)-y(j)
	dyjj=y(j+2)-y(j)
      ty2=(-Te(i,j+2)*dyj**2+Te(i,j+1)*dyjj**2-Te(i,j)*(dyjj**2-dyj**2))
	//(dyj*dyjj*(dyjj-dyj))
	AAbase1=AAbase1+ty1*(x(i+1)-x(i-1))*.5
	AAbase2=AAbase2+ty2*(x(i+1)-x(i-1))*.5
 	end do
	AAbase1=AAbase1/(s)
	AAbase2=AAbase2/(s)
	do i=ic1,ic2
      ty1=(Te(i,j+1)-Te(i,j))/(y(j+1)-y(j))
  	dyj=y(j+1)-y(j)
	dyjj=y(j+2)-y(j)
      ty2=(-Te(i,j+2)*dyj**2+Te(i,j+1)*dyjj**2-Te(i,j)*(dyjj**2-dyj**2))
	//(dyj*dyjj*(dyjj-dyj))
    	write(20,*)x(i),ty2/AAbase2
 	end do
c===============Right Fin   
	write(20,*)"===Right Fin==" 
      i =ic2
      AAright=0
	do j=1,jc2
      tx=-(Te(i-1,j)-Te(i,j))/(x(i)-x(i-1))
 	b1=tx*(1.)/(Te(i,j)-Tab)
      AAright=AAright+b1*(y(j+1)-y(j-1))*.5
	write(20,*)y(j),b1
 	end do
	AAright=AAright/h2
c===============Over All===========   	 
c===============Over All===========   	 
	i=ic1
      BBleft=0
	do j=2,jc1-1
	dxi=x(i)-x(i-1)
	dxii=x(i)-x(i-2)
	tx=-(Te(i+1,j)-Te(i,j))/(x(i+1)-x(i))
	BBleft=BBleft+tx*(y(j+1)-y(j))
     	end do
c=============
	BBbase=0.
	j=1
      do i=ic1,ic2
 	dyj=y(j+1)-y(j)
	dyjj=y(j+2)-y(j)
c      ty=(-Te(i,j+2)*dyj**2+Te(i,j+1)*dyjj**2-Te(i,j)*(dyjj**2-dyj**2))
c	//(dyj*dyjj*(dyjj-dyj))
	ty=(Te(i,j+1)-Te(i,j))/(y(j+1)-Y(j))
   	BBbase=BBbase+ty*(x(i+1)-x(i-1))*.5
  	end do
c=============
	BBright=0
	i=ic2
	do j=2,jc2-1
	dxi=x(i)-x(i-1)
	dxii=x(i)-x(i-2)
c      tx=(-Te(i-2,j)*dxi**2+Te(i-1,j)*dxii**2-Te(i,j)*(dxii**2-dxi**2))
c	//(dxi*dxii*(dxii-dxi))
	tx=-(Te(i-1,j)-Te(i,j))/(x(i)-x(i-1))
	BBright=BBright+tx*(y(j+1)-y(j))
      end do
c=============
  	AAover=(abs(BBright)+abs(BBbase)+abs(BBleft))/(1+s+h2)

	i =ic2
      do j=1,jc1
	write(21,*)y(j),Te(i,j)/Tab
	end do
c=========================================================================
c===============  PLOT  PLOT =============================================
c===============  PLOT  PLOT =============================================
c=========================================================================
1000     write(*,*)"PLOT TO FILE"
  	open(0,file="0 Details.txt") 
      open(1,file="1 W.plt")
  	open(2,file="2 Te.plt")
      open(3,file="3 HTI.plt")
 	open(4,file="4 FFI.plt")

 	write(1,*)"zone", "   i=", ny, "   j=", nx
 	write(2,*)"zone", "   i=", ny, "   j=", nx
 	write(3,*)"zone", "   i=", ny, "   j=", nx
 	write(4,*)"zone", "   i=", ny, "   j=", nx
  
      do i = 1,nx
 	do j = 1,ny
 	write(1,*) x(i),y(j),W(i,j)
 	write(2,*) x(i),y(j),Te(i,j)
 	write(3,*) x(i),y(j),FFI(i,j)
 	write(4,*) x(i),y(j),HTI(i,j) 
 	
       end do
       end do
    
	close(1) 
	close(2) 
	close(1) 
	close(2) 
	
      write(*,*)iteration
c==============Details
	write(0,21)t,rksf
	write(0,22)qTe,qW,g1
	write(0,23)n1,nc,nt,ns 
	write(0,*)"Nx=",nx,"        Ny=",ny
	write(0,*)"Ic1=",ic1,"   Ic2=",ic2,"   Jc1=",jc1,"   Jc2=",jc2
	write(0,*)"=================================="
	write(0,*)"Rx=",rx,"        Ry=",ry1,ryc
	write(0,*)"  C  =",c
	write(0,*)"  S  =",s
	write(0,*)"Omega=",rksf*t/2 
      write(0,*)"=================================="
	write(0,*)" W _ bar=",wbar
	write(0,*)" fRe_fRe=",32./wbar*((1.+c)*(s)/(1.+h2+2.*s))**2
	write(0,*)"T_b  T_b=",Tab
      write(0,*)"=================================="
	write(0,*)"Nusselt   Overall=",AAover/(Tab),1/(1+s+h2)/(Tab)
	write(0,*)"Nusselt     Left=",BBleft/1
	write(0,*)"Nusselt     Base=",BBbase/s/Tab
	write(0,*)"Nusselt    Right=",BBright/h2
	write(0,*)"Fin/Total=",(BBleft+BBright)/(BBleft-(BBbase)+BBright)
 	write(0,*)AAleft,AAbase1,AAbase1,AAright
      write(0,*)"=================================="
	SumEG=SumFFI+SumHTI
	write(0,*)"<FFI> =",SumFFI
	write(0,*)"<HTI> =",SumHTI
	write(0,*)"< EG> =",SumEG

	stop
20	format(i7," >> ",i3,i3,2e10.4," *** ",2e13.6)
200	format(i7," >> ",i3,i3,2e10.4," *** ",2e13.6)
21	format("T=",f7.4,"           K_f/F=",f7.0)
22	format("Error Te=",f4.1,"     Error W=",f4.1,"     S.U.R=",f4.2) 
23    format("nYFin=",i3,"   nYTop=",i3,"     nXFin=",i3,"   nBase=",i3)
c 	write(*,20)iteration,ii,jj,diflanda,Tedif,Te(ii,jj)

	end
