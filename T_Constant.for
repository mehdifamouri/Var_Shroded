      Dimension Te(350,350),W(350,250),Tep(350,250),Wp(350,350)
	/,x(350),y(350),HTI(350,350),ANS(5,5,5),
     /HTIB(350,350),FFIB(350,350),FFI(350,350),BNs(350,350),BeB(350,350)
	CHARACTER(20) plotTe(5),plotw(5),plotFFI(5),plotHTI(5)
	/,plotNS(5),plotBe(5)
	/,dettxt(5),Nutxt(5),Tetafintxt(5)
  	open(99,file="000 All.txt") 
c========= Geometry Details ==================
	t=.01
	s=.5
	c=.5
	h2=1.
	write(*,*)"H2=??"
c	read(*,*)h2
c=========  Grid Generation  ======================================
	n1=35
	nc=15
	nt=3
	ns=15

	dx=t/2./nt*6
	dy=dx
	call Grid(dx,dy,x,y,ry1,ryc,rx,n1,nc,ns,nt,nx,ny,c,s,t)
	call Name(plotTe,plotw,plotFFI,plotHTI,plotNs,plotBe
	/,dettxt,Nutxt,Tetafintxt)
c========= Intital & Error ============================================
    	Br=.1
	ZO=.05
	qTe =7.
 	qW =8.
	g1=.9
	alanda=-.5
	wbar=1.
      aomega=1.
	kk=1
c==================  Z  Z  Preparing  Z  Z==================
c==================  Z  Z  Preparing  Z  Z==================
111   iteration=0
c=================
	ic1=1+nt
	ic2=nx-nt
	jc1=2*n1+1
	do j=2,jc1
	if ((y(j).ge.h2).and.(y(j-1).lt.h2))jc2=j
	end do
c	ny=jc1
c	c=0.
	write(*,*)" H2 , Omega ==>>", h2,aOmega
	write(*,*)" Ic1 .. Jc2 ==>>",ic1,ic2,jc1,jc2
	write(*,*)" nx  ,   ny ==>>",nx,ny
      rksf=2./t*aomega
c=================
	do i=1,ic1
	do j=1,jc1
	w(i,j)=0.
	end do
	end do
c==================	
      do i=ic2,nx
	do j=1,jc2
	w(i,j)=0.
	end do
	end do
c=================
      do i = 1,nx
	do j = 1,ny
	W(i,j)=W(i,j)*wbar
	end do
	end do
c==================
      do i = 1,nx
 	do j = 1,ny
	Te(i,j)=Te(i,j)/alanda
	end do
	end do

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$     WWWWW         WWWWWW       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
1	iteration=iteration+1
	if (iteration/3000..eq.int(iteration/3000))then
	write(*,20)iteration/1000,ii,jj,difwbar,Wdif,fre1,fre2
	end if
	if (iteration.gt.30000)g1=.5
c=========== W ========= W ==================
c=========== W ========= W ==================
	do j = 2, ny-1
	do i = 2, nx-1
	
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
      end do
	end do
c===========  B.C.  W ========================
c==========Right SYMETRI
	i=nx
      do j=1,jc2
	W(i,j)=0.
	end do
	do j=jc2+1,ny-1

	dxi=x(i)-x(i-1)
	dxii=dxi
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)

      b1=(W(i-1,j)*dxi+W(i-1,j)*dxii)/(dxi*dxii*.5*(dxi+dxii))
 	b2=(W(i,j+1)*dyj+W(i,j-1)*dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	b3=1.
	
	a1=-(dxi+dxii)/(dxi*dxii*.5*(dxi+dxii))
	a2=-(dyj+dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	a3=0.
	
     	W(i,j)=Wp(i,j)+1.*(-(b1+b2+b3)/(a1+a2+a3)-Wp(i,j))

c	W(i,j)=W(i-1,j)
	end do
c==========Left SYMETRI
	i=1
      do j=1,jc1
	W(i,j)=0.
	end do
	do j=jc1+1,ny-1

	dxi=x(i+1)-x(i)
	dxii=dxi
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)

      b1=(W(i+1,j)*dxi+W(i+1,j)*dxii)/(dxi*dxii*.5*(dxi+dxii))
 	b2=(W(i,j+1)*dyj+W(i,j-1)*dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	b3=1.
	
	a1=-(dxi+dxii)/(dxi*dxii*.5*(dxi+dxii))
	a2=-(dyj+dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	a3=0.
	
     	W(i,j)=Wp(i,j)+1.*(-(b1+b2+b3)/(a1+a2+a3)-Wp(i,j))

c	W(i,j)=W(i+1,j)
	end do
c==========Shroude
	j=ny
      do i=1,nx
	W(i,j)=0.
	end do
c==========Base
 	j=1
      do i=1,nx
	W(i,j)=0.
	end do

c=========== W bar ========= W bar ==================
c=========== W bar ========= W bar ==================
   	wbarnew=.0
      do j = 1, ny-1
 	do i = 1, nx-1
	dx=x(i+1)-x(i)
	dy=y(j+1)-y(j)
     	wbarnew=wbarnew+.25*(W(i,j)+W(i+1,j)+W(i,j+1)+W(i+1,j+1))*dx*dy
      end do
	end do
	wbarnew=wbarnew/((1.+c)*s)
c=========== fRE ========= fRE ==================
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
	
      fre2=32.*(1.+C)*(S)*BBAA/wbarnew/(1+h2+2*s)**2	
      fre1=32/wbarnew*((1+c)*(s)/(1+h2+2*s))**2
c============= ?? RE START  ?? ===================
c============= ?? RE START  ?? ===================
     	difwbar=abs(wbarnew-wbar)/abs(wbarnew) 
	do i = 1,nx
     	do j = 1,ny
	if (abs(W(i,j)).ge..1e-12)then
	Wdif=Abs((W(i,j)-Wp(i,j))/W(i,j))/g1
	else
	Wdif=0.
	end if
	ii=i
	jj=j
	If (Wdif-10.**(-qW) .gt.0.0)goto 9
	end do
	end do
	If (difwbar-10.**(-qW) .gt.0.0)goto 9
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
99    do i = 1,nx
	do j = 1,ny
	W(i,j)=W(i,j)/wbarnew
	end do
	end do
	write(*,*)"Calculating W Finished   W_bar=",wbar
	iteration=0
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$     Te    Te    Te    Te       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      g1=.9
 4	iteration=iteration + 1
      if (iteration/10000..eq.int(iteration/10000))then
	write(*,200)iteration/1000,ii,jj,Tedif,diflanda,Te(ii,jj),alanda
	end if
c	if (iteration.gt.200000)g1=.7
c	if (iteration.gt.250000)g1=.5
c	if (iteration.gt.300000)g1=.1
c	if (iteration.gt.400000)g1=.05

c========== Te ========== Te ==================
c========== Te ========== Te ==================
	do j = 2, ny-1
	do i = 2, nx-1
	
	dxi=x(i)-x(i-1)
	dxii=x(i+1)-x(i)
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)

	if ((i.eq.ic1).and.(j.le.jc1)) goto 5
	if ((i.eq.ic2).and.(j.le.jc2)) goto 5
	qlanda=Tep(i,j) 
	if ((iteration.lt.10000).and.(kk.eq.1)) qlanda=1.

      b1=(Tep(i+1,j)*dxi+Tep(i-1,j)*dxii)/(dxi*dxii*.5*(dxi+dxii))
 	b2=(Tep(i,j+1)*dyj+Tep(i,j-1)*dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	
      b3=-W(i,j)*alanda*qlanda
	
	a1=-(dxi+dxii)/(dxi*dxii*.5*(dxi+dxii))
	a2=-(dyj+dyjj)/(dyj*dyjj*.5*(dyj+dyjj))
	a3=0
c      if(((i.le.ic1).and.(j.le.jc1)).or.((i.ge.ic2).and.(j.le.jc2)))then
c	Te(i,j)=Tep(i,j)+g1*(-(b1+b2)/(a1+a2)-Tep(i,j)) 
c	Te(i,j)=Te(ic1,j)
c	else
	Te(i,j)=Tep(i,j)+g1*(-(b1+b2+b3)/(a1+a2+a3)-Tep(i,j))
c	end if
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
	b1=b*Te(i,j+1)*dyj+b*Te(i,j-1)*dyjj+Te(i+1,j)/dx
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
	b1=b*Te(i,j+1)*dyj+b*Te(i,j-1)*dyjj+Te(i-1,j)/dx
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
	Te(i,j)=Te(ic1,j-1)
	end do
c==========Top Right Wall
	j=jc2
      do i=ic2,nx
	Te(i,j)=Te(ic2,j-1)
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
c=========== Landa ========= Landa ==================
c=========== Landa ========= Landa ==================
	sumlanda=0.
	if((iteration.lt.10000).and.(kk.eq.1))sumlanda=0.5

	do j = 1, ny-1
 	do i = 1, nx-1

	dx=x(i+1)-x(i)
	dy=y(j+1)-y(j)
	b1=Te(i,j)*W(i,j)+Te(i+1,j)*W(i+1,j)
	b2=Te(i,j+1)*W(i,j+1)+Te(i+1,j+1)*W(i+1,j+1)
	sumlanda=sumlanda+.25*(b2+b1)*dx*dy

6     end do
	end do
	alandanew=((1.+c)*(s))/sumlanda
c ============= ?? RE START  ?? ===================
c ============= ?? RE START  ?? ===================
     	diflanda=abs(alandanew-alanda)/abs(alandanew)/g1 
      do i = 1,nx
      do j = 1,ny
	if (abs(Te(i,j)).ge.1e-12) then
      Tedif=Abs((Te(i,j)-Tep(i,j))/Te(i,j))/g1
	else
	Tedif=0.
	end if
	ii=i
	jj=j
	If (Tedif-10.**(-qTe) .gt.0.0 )goto 8
	end do
	end do
 	if(diflanda-10.**(-qTe).gt.0.0 )goto 8
   	GoTo 88
c=========  Repalace Pervious =========================
8	do i = 1,nx
	do j = 1,ny
	Tep(i,j)=Te(i,j)
	end do
	end do

	if ((iteration.gt.20000).or.(kk.ne.1))then
	alanda=alanda+g1*(-alanda+alandanew)
	end if
	 
	if ((iteration.eq.200000).or.(iteration.eq.300000)
	/.or.(iteration.eq.400000).or.(iteration.eq.500000)
	/.or.(iteration.eq.600000).or.(iteration.eq.700000))then
	write(*,*)"Dou You Want To Continue? [1=Continue]"
      read(*,*)qqq,g1
	if (qqq.ne.1) goto 88
	end if

	goto 4
c============Replace Fi an Te
88    do i = 1,nx
	do j = 1,ny
	Te(i,j)=Te(i,j)*alanda
	end do
	end do
	write(*,*)"Calculating Te Finished   Landa=",alanda
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

      tx=(Te(i+1,j)*dxi**2-Te(i-1,j)*dxii**2+Te(i,j)*(dxii**2-dxi**2))
	//(dxi*dxii*(dxi+dxii))
	ty=(Te(i,j+1)*dyj**2-Te(i,j-1)*dyjj**2+Te(i,j)*(dyjj**2-dyj**2))
	//(dyj*dyjj*(dyj+dyjj))

	wx=( w(i+1,j)*dxi**2- w(i-1,j)*dxii**2+ w(i,j)*(dxii**2-dxi**2))
	//(dxi*dxii*(dxi+dxii))
	wy=( w(i,j+1)*dyj**2- w(i,j-1)*dyjj**2+ w(i,j)*(dyjj**2-dyj**2))
	//(dyj*dyjj*(dyj+dyjj))

           
	HTI(i,j)=tx**2+ty**2+(alanda*Te(i,j))**2
      FFI(i,j)=(wx)**2+(wy)**2


 	end do
	end do
c========== i= and i=1
      i=1
	do j=2,ny-1
	dxii=x(i+1)-x(i)
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)
     
      tx=(Te(i+1,j)-Te(i,j))/(dxii)
	ty=(Te(i,j+1)*dyj**2-Te(i,j-1)*dyjj**2+Te(i,j)*(dyjj**2-dyj**2))
	//(dyj*dyjj*(dyj+dyjj))

	wx=( w(i+1,j)- w(i,j))/(dxii)
	wy=( w(i,j+1)*dyj**2- w(i,j-1)*dyjj**2+ w(i,j)*(dyjj**2-dyj**2))
	//(dyj*dyjj*(dyj+dyjj))
	
	HTI(i,j)=tx**2+ty**2+(alanda*Te(i,j))**2
      FFI(i,j)=(wx)**2+(wy)**2

      end do
c========== i=ic2 and i=nx
      i=nx
	do j=2,ny-1
	dxi=x(i)-x(i-1)
	dyj=y(j)-y(j-1)
	dyjj=y(j+1)-y(j)
     
      tx=(Te(i,j)-Te(i-1,j))/(dxi)
	ty=(Te(i,j+1)*dyj**2-Te(i,j-1)*dyjj**2+Te(i,j)*(dyjj**2-dyj**2))
	//(dyj*dyjj*(dyj+dyjj))

	wx=( w(i,j)- w(i-1,j))/(dxi)
	wy=( w(i,j+1)*dyj**2- w(i,j-1)*dyjj**2+ w(i,j)*(dyjj**2-dyj**2))
	//(dyj*dyjj*(dyj+dyjj))
	
	HTI(i,j)=tx**2+ty**2+(alanda*Te(i,j))**2
      FFI(i,j)=(wx)**2+(wy)**2

      end do

c========== j=1
	do i=2,nx-1
	j=1
	dxi=x(i)-x(i-1)
  	dxii=x(i+1)-x(i)
	dyjj=y(j+1)-y(j)
      tx=(Te(i+1,j)*dxi**2-Te(i-1,j)*dxii**2+Te(i,j)*(dxii**2-dxi**2))
	//(dxi*dxii*(dxi+dxii))
	ty=(Te(i,j+1)-Te(i,j))/(dyjj)

	wx=( w(i+1,j)*dxi**2- w(i-1,j)*dxii**2+ w(i,j)*(dxii**2-dxi**2))
	//(dxi*dxii*(dxi+dxii))
	wy=(w(i,j+1)-w(i,j))/(dyjj)

           
	HTI(i,j)=tx**2+ty**2+(alanda*Te(i,j))**2
      FFI(i,j)=0.
     
      end do

	do i=2,nx-1
	j=ny

	dxi=x(i)-x(i-1)
  	dxii=x(i+1)-x(i)
	dyj=y(j)-y(j-1)

      tx=(Te(i+1,j)*dxi**2-Te(i-1,j)*dxii**2+Te(i,j)*(dxii**2-dxi**2))
  	//(dxi*dxii*(dxi+dxii))
	ty=(Te(i,j)-Te(i,j-1))/(dyj)

  	wx=( w(i+1,j)*dxi**2- w(i-1,j)*dxii**2+ w(i,j)*(dxii**2-dxi**2))
	//(dxi*dxii*(dxi+dxii))
	wy=(w(i,j)-w(i,j-1))/(dyj)

           
	HTI(i,j)=tx**2+ty**2+(alanda*Te(i,j))**2
      FFI(i,j)=0.
     
      end do
c=============  Fins
      do i=1,ic1
	do j=1,jc1
	FFI(i,j)=0.
	end do
	end do

      do i=ic2,nx
	do j=1,jc2
	FFI(i,j)=0.
	end do
	end do
	j=ny
	do i=1,nx
	FFI(i,j)=0.
	end do
c===========symretry
	i=1
	do j=jc1,ny
	FFI(i,j)=FFI(i+1,j)
	HTI(i,j)=HTI(i+1,j)
	end do
	i=nx
 	do j=jc2,ny
 	FFI(i,j)=FFI(i-1,j)
	HTI(i,j)=HTI(i-1,j)
	end do

c=============  COORNER
	FFI(1,1)=0.
	FFI(1,ny)=0.
	FFI(nx,1)=0.
	FFI(nx,ny)=0.

	HTI(1,1)=(HTI(1+1,1)+HTI(1,1+1))/2.
	HTI(1,ny)=(HTI(1+1,ny)+HTI(1,ny-1))/2.
	HTI(nx,1)=(HTI(nx-1,1)+HTI(nx,1+1))/2.
	HTI(nx,ny)=(HTI(nx-1,ny)+HTI(nx,ny+1))/2.

c============= Post Process =================
c============= Post Process =================
	do lll=1,5
	if (lll.eq.1)Br=.02
	if (lll.eq.2)Br=.04
	if (lll.eq.3)Br=.06
	if (lll.eq.4)Br=.08
	if (lll.eq.5)Br=.10

	do kkk=1,5
	if (kkk.eq.1)Zo=.02
	if (kkk.eq.2)Zo=.04
	if (kkk.eq.3)Zo=.06
	if (kkk.eq.4)Zo=.08
	if (kkk.eq.5)Zo=.10
 	HTIB=0.
	FFIB=0.

      do i=1,nx
	do j=1,ny
 	HTIB(i,j)=Zo**2.*HTI(i,j)/(1.+Zo*Te(i,j))**2
	FFIB(i,j)=Br*FFI(i,j)/(1.+Zo*Te(i,j))
	BNs(i,j)=FFIB(i,j)+HTIB(i,j)
	if (FFIB(i,j).eq.0.)then
	BeB(i,j)=1.
	else
      BeB(i,j)=HTIB(i,j)/BNs(i,j)
	end if
      end do
	end do
c=====
	SumFFI=0.
 	SumHTI=0.
 	SumNS=0.
      do i=1,nx-1
	do j=1,ny-1
	dx=x(i+1)-x(i)
	dy=y(j+1)-y(j)
	bf=(FFIB(i,j)+FFIB(i+1,j)+FFIB(i,j+1)+FFIB(i+1,j+1))
	bh=(HTIB(i,j)+HTIB(i+1,j)+HTIB(i,j+1)+HTIB(i+1,j+1))
	SumFFI=SumFFI+.25*bf*dx*dy
	SumHTI=SumHTI+.25*bh*dx*dy
     	end do
	end do
	SumFFI=SumFFI/(s*(1.+c))
	SumHTI=SumHTI/(s*(1.+c))
	SumNs=SumFFI+SumHTI
	ANS(kkk,lll,1)=SumFFI
	ANS(kkk,lll,2)=SumHTI
	ANS(kkk,lll,3)=SumNs
	ANS(kkk,lll,4)=SumHTI/SumNs
	end do
	end do								 
	write(*,*)"Calculating EG Finished   <Ns>=",SumNs
c=========================================================================
c===============  NUSSELT  NUSSELT =======================================	
c===============  NUSSELT  NUSSELT =======================================	
c=========================================================================
	open(20,file=Nutxt(kk)) 
	open(21,file=Tetafintxt(kk)) 
c===============left  Fin
       write(20,*)"==left  Fin==" 
	i =ic1
      AAleft=0
	do j=2,jc1
      tx=(Te(i+1,j)-Te(i,j))/(x(i+1)-x(i))
 	if ((1.-Te(i,j)).ne.0.)b1=tx*(1.)/(1-Te(i,j))
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
      tx=(Te(i-1,j)-Te(i,j))/(x(i)-x(i-1))
 	if ((1.-Te(i,j)).ne.0.)b1=tx*(1.)/(1-Te(i,j))
      AAright=AAright+b1*(y(j+1)-y(j-1))*.5
	write(20,*)y(j),b1
 	end do
	AAright=AAright/h2
c===============Over All===========   	 
c===============Over All===========   	 
	i=ic1
      BBleft=0
	do j=2,jc1
	tx=(Te(i+1,j)-Te(i,j))/(x(i+1)-x(i))
	BBleft=BBleft+tx*(y(j+1)-y(j-1))*.5
     	end do
c=============
	BBbase=0.
	j=1
      do i=ic1,ic2
	ty=(Te(i,j+1)-Te(i,j))/(y(j+1)-Y(j))
   	BBbase=BBbase+ty*(x(i+1)-x(i-1))*.5
  	end do
c=============
	BBright=0
	i=ic2
	do j=2,jc2
	tx=(Te(i-1,j)-Te(i,j))/(x(i)-x(i-1))
	BBright=BBright+tx*(y(j+1)-y(j-1))*.5
      end do
c=============
  	AAover=(abs(BBright)+abs(BBbase)+abs(BBleft))/(1+s+h2)

	i =ic2
      do j=1,jc1
	write(21,*)y(j),Te(i,j)
	end do
	write(*,*)"Calculating Nu Finished   Nu(Overall)=",AAover

c=========================================================================
c===============  PLOT  PLOT =============================================
c===============  PLOT  PLOT =============================================
c=========================================================================
1000  l=l
  	open(0,file=dettxt(kk)) 
      open(1,file=plotTe(kk))
  	open(2,file=plotw(kk))
      open(3,file=plotHTI(kk))
 	open(4,file=plotFFI(kk))
 	open(5,file=plotNs(kk))
 	open(6,file=plotBe(kk))

 	write(1,*)"zone", "   i=", ny, "   j=", nx
 	write(2,*)"zone", "   i=", ny, "   j=", nx
 	write(3,*)"zone", "   i=", ny, "   j=", nx
 	write(4,*)"zone", "   i=", ny, "   j=", nx
	write(5,*)"zone", "   i=", ny, "   j=", nx
	write(6,*)"zone", "   i=", ny, "   j=", nx
  
      do i = 1,nx
 	do j = 1,ny
 	write(1,*) x(i),y(j),Te(i,j)
 	write(2,*) x(i),y(j),W(i,j)
 	write(3,*) x(i),y(j),FFIB(i,j)
 	write(4,*) x(i),y(j),HTIB(i,j) 
	write(5,*) x(i),y(j),BNs(i,j) 
	write(6,*) x(i),y(j),BeB(i,j) 
 	
      end do
      end do

	close(1) 
	close(2) 
	close(3) 
	close(4)
	close(5)
	close(6)
      close(20)
	close(21) 
c==============Details
	write(0,22)qTe,qW,g1
	write(0,23)n1,nc,ns 
	write(0,*)"Nx=",nx,"        Ny=",ny
	write(0,*)"Ic1=",ic1,"   Ic2=",ic2,"   Jc1=",jc1,"   Jc2=",jc2
	write(0,*)"Rx=",rx,"        Ry=",ry1,ryc
	write(0,*)"=================================="
	write(0,*)"  H2 =",h2 
      write(0,*)"  C  =",c
	write(0,*)"  S  =",s
	write(0,*)"Omega=",aomega 
      write(0,*)"=================================="
	write(0,*)" W _ bar=",wbar
	write(0,*)" fRe_fRe=",fre1,fre2
	write(0,*)"   Landa=",alanda
      write(0,*)"=================================="
	write(0,*)"Nusselt   Overall=",AAover
	write(0,*)"Nusselt/Left/Base/Right=",BBleft/1,BBbase/s,BBright/h2
	write(0,*)"Fin/Total=",(BBleft+BBright)/(BBleft+BBbase+BBright)
      write(0,*)"=================================="

	write(0,*)" Br=.02   Omega==>>  .02   .04   .06    .08   .1"
	write(0,24)ANS(1,1,1),ANS(2,1,1),ANS(3,1,1),ANS(4,1,1),ANS(5,1,1)
	write(0,25)ANS(1,1,2),ANS(2,1,2),ANS(3,1,2),ANS(4,1,2),ANS(5,1,2)
	write(0,26)ANS(1,1,3),ANS(2,1,3),ANS(3,1,3),ANS(4,1,3),ANS(5,1,3)
	write(0,27)ANS(1,1,4),ANS(2,1,4),ANS(3,1,4),ANS(4,1,4),ANS(5,1,4)
      write(0,*)
	write(0,*)" Br=.04   Omega==>>  .02   .04   .06    .08   .1"
	write(0,24)ANS(1,2,1),ANS(2,2,1),ANS(3,2,1),ANS(4,2,1),ANS(5,2,1)
	write(0,25)ANS(1,2,2),ANS(2,2,2),ANS(3,2,2),ANS(4,2,2),ANS(5,2,2)
	write(0,26)ANS(1,2,3),ANS(2,2,3),ANS(3,2,3),ANS(4,2,3),ANS(5,2,3)
	write(0,27)ANS(1,2,4),ANS(2,2,4),ANS(3,2,4),ANS(4,2,4),ANS(5,2,4)
      write(0,*)
	write(0,*)" Br=.06   Omega==>>  .02   .04   .06    .08   .1"
	write(0,24)ANS(1,3,1),ANS(2,3,1),ANS(3,3,1),ANS(4,3,1),ANS(5,3,1)
	write(0,25)ANS(1,3,2),ANS(2,3,2),ANS(3,3,2),ANS(4,3,2),ANS(5,3,2)
	write(0,26)ANS(1,3,3),ANS(2,3,3),ANS(3,3,3),ANS(4,3,3),ANS(5,3,3)
	write(0,27)ANS(1,3,4),ANS(2,3,4),ANS(3,3,4),ANS(4,3,4),ANS(5,3,4)
      write(0,*)
	write(0,*)" Br=.08   Omega==>>  .02   .04   .06    .08   .1"
	write(0,24)ANS(1,4,1),ANS(2,4,1),ANS(3,4,1),ANS(4,4,1),ANS(5,4,1)
	write(0,25)ANS(1,4,2),ANS(2,4,2),ANS(3,4,2),ANS(4,4,2),ANS(5,4,2)
	write(0,26)ANS(1,4,3),ANS(2,4,3),ANS(3,4,3),ANS(4,4,3),ANS(5,4,3)
	write(0,27)ANS(1,4,4),ANS(2,4,4),ANS(3,4,4),ANS(4,4,4),ANS(5,4,4)
      write(0,*)
	write(0,*)" Br=.10   Omega==>>  .02   .04   .06    .08   .1"
	write(0,24)ANS(1,5,1),ANS(2,5,1),ANS(3,5,1),ANS(4,5,1),ANS(5,5,1)
	write(0,25)ANS(1,5,2),ANS(2,5,2),ANS(3,5,2),ANS(4,5,2),ANS(5,5,2)
	write(0,26)ANS(1,5,3),ANS(2,5,3),ANS(3,5,3),ANS(4,5,3),ANS(5,5,3)
	write(0,27)ANS(1,5,4),ANS(2,5,4),ANS(3,5,4),ANS(4,5,4),ANS(5,5,4)

      close(0)
	write(*,*)"PLOT TO FILE" 
      write(*,*)"====================================================="
      write(99,*)"=================  ",aOmega,"  ================="
	write(99,*)fre1,AAover
	write(99,24)ANS(1,1,1),ANS(2,1,1),ANS(3,1,1),ANS(4,1,1),ANS(5,1,1)
	write(99,25)ANS(1,1,2),ANS(2,1,2),ANS(3,1,2),ANS(4,1,2),ANS(5,1,2)
	write(99,*)"Br=.04"
	write(99,24)ANS(1,2,1),ANS(2,2,1),ANS(3,2,1),ANS(4,2,1),ANS(5,2,1)
	write(99,25)ANS(1,2,2),ANS(2,2,2),ANS(3,2,2),ANS(4,2,2),ANS(5,2,2)
	write(99,*)"Br=.06"
	write(99,24)ANS(1,3,1),ANS(2,3,1),ANS(3,3,1),ANS(4,3,1),ANS(5,3,1)
	write(99,25)ANS(1,3,2),ANS(2,3,2),ANS(3,3,2),ANS(4,3,2),ANS(5,3,2)
	write(99,*)"Br=.08"
	write(99,24)ANS(1,4,1),ANS(2,4,1),ANS(3,4,1),ANS(4,4,1),ANS(5,4,1)
	write(99,25)ANS(1,4,2),ANS(2,4,2),ANS(3,4,2),ANS(4,4,2),ANS(5,4,2)
	write(99,*)"Br=.10"
	write(99,24)ANS(1,5,1),ANS(2,5,1),ANS(3,5,1),ANS(4,5,1),ANS(5,5,1)
	write(99,25)ANS(1,5,2),ANS(2,5,2),ANS(3,5,2),ANS(4,5,2),ANS(5,5,2)

	Kk=kk+1
	if (kk.eq.2)aOmega=5.
	if (kk.eq.3)aOmega=10.
	if (kk.eq.4)aOmega=25.
	if (kk.eq.5)aOmega=1000.
	if (kk.eq.6)stop
      write(*,*)"Omega=",aOmega
	write(*,*)"====================================================="
	goto 111
	stop

20	format(i7," >> ",i3,i3,2e10.2," *** ",2e10.4)
200	format(i7," >> ",i3,i3,2e10.2," *** ",2e10.3)
21	format("T=",f7.4,"           K_f/F=",f7.0)
22	format("Error Te=",f4.1,"     Error W=",f4.1,"     S.U.R=",f4.2) 
23    format("nYFin=",i3,"   nYTop=",i3,"     nXFin=",i3,"   nBase=",i3)
24    format("<FFI> =",5e14.6)
25    format("<HTI> =",5e14.6)
26    format("  < Ns> =",5e14.6)
27    format("  < Be> =",5e14.6)

	end
c==============================================================
c================	subroutine Name   ===========================
c==============================================================
	subroutine Name(plotTe,plotw,plotFFI,plotHTI,plotNs,plotBe
	/,dettxt,Nutxt,Tetafintxt)

 	CHARACTER(20) plotTe(5),plotw(5),plotFFI(5),plotHTI(5)
	/,dettxt(5),Nutxt(5),Tetafintxt(5),plotNs(5),plotBe(5)

	plotTe(1)='Om=01 Temp.plt'
	plotTe(2)='Om=05 Temp.plt'
	plotTe(3)='Om=10 Temp.plt'
	plotTe(4)='Om=25 Temp.plt'
	plotTe(5)='Om=99 Temp.plt'

	plotw(1)='Om=01 WWW.plt'
 	plotw(2)='Om=05 WWW.plt'
 	plotw(3)='Om=10 WWW.plt'
	plotw(4)='Om=25 WWW.plt'
	plotw(5)='Om=99 WWW.plt'

	plotFFI(1)='Om=01 FFI.plt'
 	plotFFI(2)='Om=05 FFI.plt'
 	plotFFI(3)='Om=10 FFI.plt'
	plotFFI(4)='Om=25 FFI.plt'
	plotFFI(5)='Om=99 FFI.plt'

	plotHTI(1)='Om=01 HTI.plt'
 	plotHTI(2)='Om=05 HTI.plt'
 	plotHTI(3)='Om=10 HTI.plt'
	plotHTI(4)='Om=25 HTI.plt'
	plotHTI(5)='Om=99 HTI.plt'

	plotNs(1)='Om=01 Ns.plt'
 	plotNs(2)='Om=05 Ns.plt'
 	plotNs(3)='Om=10 Ns.plt'
	plotNs(4)='Om=25 Ns.plt'
	plotNs(5)='Om=99 Ns.plt'

	plotBe(1)='Om=01 Be.plt'
 	plotBe(2)='Om=05 Be.plt'
 	plotBe(3)='Om=10 Be.plt'
	plotBe(4)='Om=25 Be.plt'
	plotBe(5)='Om=99 Be.plt'

	dettxt(1)='Om=01 Details.txt'
 	dettxt(2)='Om=05 Details.txt'
 	dettxt(3)='Om=10 Details.txt'
	dettxt(4)='Om=25 Details.txt'
	dettxt(5)='Om=99 Details.txt'

 	Nutxt(1)='Om=01 Nusselt.txt'
 	Nutxt(2)='Om=05 Nusselt.txt'
 	Nutxt(3)='Om=10 Nusselt.txt'
	Nutxt(4)='Om=25 Nusselt.txt'
	Nutxt(5)='Om=99 Nusselt.txt'

	Tetafintxt(1)='Om=01 Teta_fin.txt'
 	Tetafintxt(2)='Om=05 Teta_fin.txt'
 	Tetafintxt(3)='Om=10 Teta_fin.txt'
	Tetafintxt(4)='Om=25 Teta_fin.txt'
	Tetafintxt(5)='Om=99 Teta_fin.txt'
	return
	end
c==============================================================
c==============	subroutine Grid =============================
c==============================================================
	subroutine Grid(dx,dy,x,y,ry1,ryc,rx,n1,nc,ns,nt,nx,ny,c,s,t)

	dimension x(350),y(350)
c=========  R Function Iteraion  ================== 
      rx=1.1
 	ry1=1.1
	ryc=1.1
      do i=1,4000
	rx=rx+.1*(-rx+(1.-(s-t)/2./dx*(1-rx))**(1./(ns)))
	ryc=ryc+.1*(-ryc+(1.-c/2./dy*(1-ryc))**(1./nc))
	ry1=ry1+.1*(-ry1+(1.-1./2./dy*(1-ry1))**(1./(n1)))
      end do
c=====================================================
	dxo=dx
	dyo=dy
	y(1)=0
      do j=2,n1+1
	y(j)=y(j-1)+dy
	dy=ry1*dy
      end do
	write(*,*)y(n1+1)
	jj=0
	do j=n1+1-1,1,-1
	jj=jj+1
	y(n1+1+jj)=y(n1+1)+(y(n1+1)-y(j))
	end do
	write(*,*)y(n1+1+jj)
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
	write(*,*)"  rx  , ry  ==>>",rx,ry1,ryc

	return
	end