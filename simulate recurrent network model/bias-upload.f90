!	simulations for the bias statistics
!	questions to david.hansel@parisdescartes.fr
	implicit none
	integer n,nc,nu,nce,nci
	integer ne,ni,nn,nnn
	real*8 kkee,kkei,kkii,kkie
	real*8 kkee0,kkei0,kkii0,kkie0
	parameter(n=32000)  
	parameter(nu=1) 
        parameter(nc=4000,nce=4000,nci=4000)
        parameter(ne=32000,ni=8000)
	integer i,j,k,iu,it,itmax
	integer itcueon,itcueoff,itrans
	integer i1,i2,i3,npref1,npref2,narch,iarch
real*8, allocatable :: switche(:),switchi(:)
real*8, allocatable :: se(:),si(:)		  
real*8, allocatable :: ue(:),ui(:)	
real*8, allocatable :: semact(:),simact(:)
real*8, allocatable :: seed(:),seem(:)
real*8, allocatable :: siid(:),siim(:)
real*8, allocatable :: thetae(:),thetai(:)
real*8, allocatable :: he(:),hi(:) 
real*8, allocatable :: hee(:),hie(:),hei(:),hii(:)
real*8, allocatable :: heesp(:),hiesp(:),heisp(:),hiisp(:)
real*8, allocatable :: fee2(:),fei2(:),fii2(:),fie2(:)
real*8, allocatable :: sebis(:),sibis(:)
real*8, allocatable :: tspike(:),tspikeold(:)
real*8, allocatable :: tspiki(:),tspikiold(:)
real*8,  allocatable :: sumspe(:),sumsp2e(:),sumspi(:)
real*8,  allocatable :: sumsp2i(:)
real*8,  allocatable :: sumecv2(:),sumicv2(:)
real*8, allocatable :: dte(:),dteold(:),dti(:),dtiold(:)
real*8, allocatable :: fire(:),firi(:)
real*8, allocatable :: cve(:),cvi(:)
real*8, allocatable :: cv2e(:),cv2i(:)
real*8, allocatable :: inputsponte(:),inputsponti(:)
real*8, allocatable :: inputcuee(:),inputcuei(:)
integer, allocatable :: mee(:,:),mei(:,:),mii(:,:),mie(:,:)
integer, allocatable :: meesp(:,:),meisp(:,:),miisp(:,:),miesp(:,:)
integer, allocatable :: ncee(:),ncei(:),ncii(:),ncie(:)
integer, allocatable :: nceesp(:),nceisp(:),nciisp(:),nciesp(:)
integer, allocatable :: nspike(:),nspiki(:)
	real*8 eps,taufac,taurec,bigu  
	real*8 taufacee,taurecee,biguee
	real*8 taufacii,taurecii,biguii
	real*8 gee,gei,gii,gie
	real*8 ieon,ieoff,iion,iioff
	real*8 epsilonee,epsilonie,epsilonei,epsilonii,epsilon
	real*8 sigmaee,sigmaei,sigmaii,sigmaie,sigmainp   
	real*8 gpec
	real*8 tauact,epsact,umact
	real*8 sumg,sigma,sumgee 		
	real*8 xx,yy,zz
	real*8 ggubfs
	real*8 stim,time,threshdec,threshact
	real*8 avfire1,avfire2,avfiri1,avfiri2
	integer ng,ngee		
        integer ngei,ngie,ngii
	real*8 f,x
	real*8 rnd,pi,dpi
	real*8 amplecue1,amplicue1,amplecue2,amplicue2,theta
	real*8 thetacue1,thetacue2,amplion,ampleon
	real*8 inputsponte0,inputsponti0
	real*8 inputcuee0,inputcuei0
	real*8 semarch,seminit
	real*8 timedec,avtimedec
	character*300 nombfile,nombfile0
        character*5 type_ee,type_ei,type_ie,type_ii
	character*15 type_anneal
	character*300 r(100),rarch(100)
        real*8 tee2,tei2,tii2,tie2
	real*8 s1,s2,s3,dtf,dt2
	real*8 exee2,exei2,exii2,exie2
	real*8 expi2
	real*8 dt,eps9
	real*8 ser,set,sir,sit,spr,spt,dtf2
	real*8 snew,uu
	real*8 taue,taui
	integer deltit,idir,anneal
	integer itrial,ntrial,nact1
allocate (switche(ne),switchi(ni))
allocate (se(n),si(n))
allocate (ue(3),ui(3))
allocate (semact(n),simact(n))
allocate (seed(n),seem(n))
allocate (siid(n),siim(n))
allocate (thetae(ne),thetai(ni))
allocate (he(ne),hi(ni))
allocate (hee(ne),hie(ni),hei(ne),hii(ni))
allocate (heesp(ne),hiesp(ni),heisp(ne),hiisp(ni))
allocate (fee2(n),fei2(n),fii2(n),fie2(n))
allocate (sebis(ne),sibis(ni))
allocate (tspike(n),tspikeold(n))
allocate (tspiki(n),tspikiold(n))
allocate (sumspe(n),sumsp2e(n),sumspi(n))
allocate (sumsp2i(n))
allocate (sumecv2(n),sumicv2(n))
allocate (dte(n),dteold(n),dti(n),dtiold(n))
allocate (fire(n),firi(n))
allocate (cve(n),cvi(n))
allocate (cv2e(n),cv2i(n))
allocate (inputsponte(ne),inputsponti(ni))
allocate (inputcuee(ne),inputcuei(ni))
allocate (mee(nce,ne),mei(nce,ni),mii(nci,ni),mie(nci,ne))
allocate (meesp(nce,ne),meisp(nce,ni),miisp(nci,ni),miesp(nci,ne))
allocate (ncee(n),ncei(n),ncii(n),ncie(n))
allocate (nceesp(n),nceisp(n),nciisp(n),nciesp(n))
allocate (nspike(n),nspiki(n))

	eps9=1e-09
	uu=0.0001
	ntrial=500
	narch=10
	open(2,file='bias-upload.ini')
	read(2,*) semarch 
	read(2,*) seminit
	read(2,*) itmax
	read(2,*) dtf
	read(2,*) deltit
	read(2,*) gee
	read(2,*) gei
	read(2,*) gii
	read(2,*) gie
	read(2,*) epsilonee
	read(2,*) epsilonei
	read(2,*) epsilonii
	read(2,*) epsilonie
	read(2,*) sigmaee
	read(2,*) sigmaei
	read(2,*) sigmaii
	read(2,*) sigmaie
	read(2,*) taue
	read(2,*) taui
        read(2,*) tee2
        read(2,*) tei2
        read(2,*) tii2
        read(2,*) tie2
        read(2,*) ser
        read(2,*) set
        read(2,*) sir
        read(2,*) sit
	read(2,*) kkee
	read(2,*) kkei
	read(2,*) kkii
	read(2,*) kkie
	read(2,*) tauact
	read(2,*) itcueon
	read(2,*) itcueoff
	read(2,*) itrans
 	read(2,*) inputsponte0
 	read(2,*) inputsponti0
 	read(2,*) inputcuee0
 	read(2,*) inputcuei0
	read(2,*) threshdec
	read(2,*) nnn
	read(2,90) nombfile0
	close(2)
90	format(a300)

!	open(4,file='conect-choice-'//nombfile0)
 	open(96,file='bias-choice-'//nombfile0)
 	open(98,file='choice-choice-'//nombfile0)
!	open(97,file='orderparameter-choice-'//nombfile0)
!	open(99,file='numrec-choice-'//nombfile0)
!	open(100,file='numrecsparse-choice-'//nombfile0)
!
	kkee0=kkee
	kkei0=kkei
	kkie0=kkie
	kkii0=kkii


        pi=dacos(-1.d0)
        dpi=2.d0*pi
        dtf2=dtf/2
        eps=1/taufac
        epsact=1/tauact
	exee2=exp(-dtf/tee2)
	exei2=exp(-dtf/tei2)
	exii2=exp(-dtf/tii2)
	exie2=exp(-dtf/tie2)

	umact=exp(-dtf/tauact)

	open(12,file='pr-choice-'//nombfile0)
	write(12,*) 'ne,ni=',ne,ni
	write(12,*) 'semarch=',semarch
	write(12,*) 'seminit=',seminit
	write(12,*) 'itmax=',itmax
    	write(12,*) 'ser=',ser
        write(12,*) 'set=',set
        write(12,*) 'sir=',sir
        write(12,*) 'sit=',sit
	write(12,*) 'gee=',gee
	write(12,*) 'gei=',gei
	write(12,*) 'gii=',gii
	write(12,*) 'gie=',gie
	write(12,*) 'epsilonee=',epsilonee
	write(12,*) 'epsilonei=',epsilonei
	write(12,*) 'epsilonii=',epsilonii
	write(12,*) 'epsilonie=',epsilonie
	write(12,*) 'sigmaee =',sigmaee
	write(12,*) 'sigmaei=',sigmaei
	write(12,*) 'sigmaii=',sigmaii
	write(12,*) 'sigmaie=',sigmaie
        write(12,*) 'ne=',ne
	write(12,*) 'ni=',ni
	write(12,*) 'kkee=',kkee0
	write(12,*) 'kkei=',kkei0
	write(12,*) 'kkii=',kkii0
	write(12,*) 'kkie=',kkie0
	write(12,*) 'itrans=',itrans
        write(12,*) 'inputsponte0=',inputsponte0
        write(12,*) 'inputsponti0=',inputsponti0
        write(12,*) 'inputcuee0=',inputcuee0
        write(12,*) 'inputcuei0=',inputcuei0
	write(12,*) 'threshdec=',threshdec
	write(12,*) 'tauact=',tauact
	write(12,*) 'tee2=',tee2
	write(12,*) 'tei2=',tei2
	write(12,*) 'tii2=',tii2
	write(12,*) 'tie2=',tie2
	write(12,*) 'taue=',taue
	write(12,*) 'taui=',taui
	write(12,*) 'nombfile=',nombfile0
	call flush(12)


	gee=gee/sqrt((kkee))/tee2
	gei=gei/sqrt((kkei))/tei2
	gii=gii/sqrt((kkii))/tii2
	gie=gie/sqrt((kkie))/tie2
	do i=1,ne
	thetae(i)=dpi*float(i)/float(ne)
	end do
	do i=1,ni
	thetai(i)=dpi*float(i)/float(ni)
	end do
 	rarch(1)=trim(nombfile0)//'-net01'
 	rarch(2)=trim(nombfile0)//'-net02'
 	rarch(3)=trim(nombfile0)//'-net03'
 	rarch(4)=trim(nombfile0)//'-net04'
 	rarch(5)=trim(nombfile0)//'-net05'
 	rarch(6)=trim(nombfile0)//'-net06'
 	rarch(7)=trim(nombfile0)//'-net07'
 	rarch(8)=trim(nombfile0)//'-net08'
 	rarch(9)=trim(nombfile0)//'-net09'
 	rarch(10)=trim(nombfile0)//'-net10'
	do iarch=1,narch
	avtimedec=0.d0
	nact1=0
	nombfile=rarch(iarch)

  	call crmat(semarch,ne,ne,sigmaee,kkee,nce,'ee',ncee,mee,mei,mie,mii)
        call crmat(semarch,ne,ni,sigmaei,kkei,nce,'ei',ncei,mee,mei,mie,mii)
        call crmat(semarch,ni,ne,sigmaie,kkie,nci,'ie',ncie,mee,mei,mie,mii)
        call crmat(semarch,ni,ni,sigmaii,kkii,nci,'ii',ncii,mee,mei,mie,mii)
  	call crmatsp(semarch,ne,ne,sigma,kkee,nce,'ee',nceesp,meesp,meisp,miesp,miisp)
        call crmatsp(semarch,ne,ni,sigma,kkei,nce,'ei',nceisp,meesp,meisp,miesp,miisp)
        call crmatsp(semarch,ni,ne,sigma,kkie,nci,'ie',nciesp,meesp,meisp,miesp,miisp)
        call crmatsp(semarch,ni,ni,sigma,kkii,nci,'ii',nciisp,meesp,meisp,miesp,miisp)


	inputcuee=dsqrt(kkee)*inputcuee0
	inputcuei=dsqrt(kkii)*inputcuei0
	inputsponte=dsqrt(kkee)*inputsponte0
	inputsponti=dsqrt(kkii)*inputsponti0

 	r(1)=trim(nombfile)//'-r01'
 	r(2)=trim(nombfile)//'-r02'
 	r(3)=trim(nombfile)//'-r03'
 	r(4)=trim(nombfile)//'-r04'
 	r(5)=trim(nombfile)//'-r05'
 	r(6)=trim(nombfile)//'-r06'
 	r(7)=trim(nombfile)//'-r07'
 	r(8)=trim(nombfile)//'-r08'
 	r(9)=trim(nombfile)//'-r09'
 	r(10)=trim(nombfile)//'-r10'
 	r(11)=trim(nombfile)//'-r11'
 	r(12)=trim(nombfile)//'-r12'
 	r(13)=trim(nombfile)//'-r13'
 	r(14)=trim(nombfile)//'-r14'
 	r(15)=trim(nombfile)//'-r15'
 	r(16)=trim(nombfile)//'-r16'
 	r(17)=trim(nombfile)//'-r17'
 	r(18)=trim(nombfile)//'-r18'
 	r(19)=trim(nombfile)//'-r19'
 	r(20)=trim(nombfile)//'-r20'
 	r(21)=trim(nombfile)//'-r21'
 	r(22)=trim(nombfile)//'-r22'
 	r(23)=trim(nombfile)//'-r23'
 	r(24)=trim(nombfile)//'-r24'
 	r(25)=trim(nombfile)//'-r25'
 	r(26)=trim(nombfile)//'-r26'
 	r(27)=trim(nombfile)//'-r27'
 	r(28)=trim(nombfile)//'-r28'
 	r(29)=trim(nombfile)//'-r29'
 	r(30)=trim(nombfile)//'-r30'
	  do itrial=1,ntrial
!       epsilon=0.1*float(itrial-1)+1.0
!       epsilon=2.5
!       epsilonei=epsilon
!       epsilonii=epsilon


	 	fee2=0
	 	fie2=0
	 	fei2=0
	 	fii2=0
	 	he=0 
	 	hi=0 
	 	hee=0 
	 	hei=0 
	 	hii=0 
	 	hie=0 
	 	heesp=0 
	 	heisp=0 
	 	hiisp=0 
	 	hiesp=0 
                nspike=0
                nspiki=0
                tspikeold=0
                tspikiold=0
                sumspe=0
                sumsp2e=0
                sumspi=0
                sumsp2i=0
		sumecv2=0
		sumicv2=0
		dte=0
		dteold=0
		dti=0
		dtiold=0
		semact=0
		simact=0
	do i=1,n
		rnd=ggubfs(seminit)
	 	se(i)=ser+rnd*(set-ser)
		rnd=ggubfs(seminit)
	 	si(i)=sir+rnd*(sit-sir)
! 	seed(i)=1.
! 	siid(i)=1.
! 	seem(i)=biguee
! 	siim(i)=biguii
	end do

!	open(3,file='ux-choice-'//r(itrial))
!       open(26,file='spike-choice-'//r(itrial))
!       open(27,file='spiki-choice-'//r(itrial))
!       open(31,file='fire-choice-'//r(itrial))
!       open(32,file='firi-choice-'//r(itrial))
                nspike=0
                nspiki=0
		stim=0.

	do it=1,itmax
	time=(float(it)+float(itrial-1)*float(itmax))*dtf*0.01d0
	if(it.eq.itcueon) stim=1
	if(it.eq.itcueoff) then
	stim=0
	end if
	do k=1,ne
 	if(semact(k).gt.eps9) semact(k)=semact(k)*umact
	end do  
	do k=1,ni
 	if(simact(k).gt.eps9) simact(k)=simact(k)*umact
	end do  

	hee=hee*exee2
	heesp=heesp*exee2
	do k=1,ne
	if(switche(k).eq.1) then
 	do j=1,ncee(k)
 	nn=mee(j,k)
	hee(nn)=hee(nn)+fee2(k)	
	end do
 	do j=1,nceesp(k)
 	nn=meesp(j,k)
	heesp(nn)=heesp(nn)+fee2(k)	
	end do
	end if
	end do

	hei=hei*exei2
	heisp=heisp*exei2
	do k=1,ni
	if(switchi(k).eq.1) then
 	do j=1,ncei(k)
	nn=mei(j,k)
	hei(nn)=hei(nn)+fei2(k)
	end do
 	do j=1,nceisp(k)
	nn=meisp(j,k)
	heisp(nn)=heisp(nn)+fei2(k)
	end do
	end if
	end do

	he=gee*(hee+heesp*epsilonee)+gei*(hei+heisp*epsilonei)+inputsponte+stim*inputcuee

	hie=hie*exie2
	hiesp=hiesp*exie2
	do k=1,ne
	if(switche(k).eq.1) then
        do j=1,ncie(k)
        nn=mie(j,k)
	hie(nn)=hie(nn)+fie2(k)
	end do
        do j=1,nciesp(k)
        nn=miesp(j,k)
	hiesp(nn)=hiesp(nn)+fie2(k)
	end do
	end if
	end do	     

	hii=hii*exii2
	hiisp=hiisp*exii2
	do k=1,ni
	if(switchi(k).eq.1) then
        do j=1,ncii(k)
	nn=mii(j,k)
	hii(nn)=hii(nn)+fii2(k) 
	end do
        do j=1,nciisp(k)
	nn=miisp(j,k)
	hiisp(nn)=hiisp(nn)+fii2(k) 
	end do
	end if
	end do
	hi=gie*(hie+hiesp*epsilonie)+gii*(hii+hiisp*epsilonii)+inputsponti+stim*inputcuei


      sebis=se/(1+dtf/taue)+dtf*he/(taue+dtf)
      sibis=si/(1+dtf/taui)+dtf*hi/(taui+dtf)



	do k=1,ne
	switche(k)=0
		if(sebis(k).gt.set) then
			switche(k)=1
			dt2=dtf*(sebis(k)-set)/(sebis(k)-se(k))
	snew=ser+(sebis(k)-set)*(1+dtf/taue*(-ser+se(k))/(sebis(k)-se(k)))
                        tspike(k)=it*dtf
                        dte(k)=tspike(k)-tspikeold(k)
        seem(k)=seem(k)*(1.-biguee)*dexp(-dte(k)/taufacee)+biguee
                xx=dexp(-dte(k)/taurecee)
               	seed(k)=seed(k)*(1.-seem(k))*xx+1.-xx
!  		fee2(k)=dexp(-dt2/tee2)*seem(k)*seed(k)
   		fee2(k)=exp(-dt2/tee2)
!  		fie2(k)=exp(-dt2/tie2)*seem(k)*seed(k)
   		fie2(k)=exp(-dt2/tie2)
!	write(26,*) (it)*dtf*0.01,k
!	call flush(26)
			se(k)=snew
 		semact(k)=semact(k)+epsact*exp(-dt2/tauact)
			if(it.gt.itrans .and. it.lt.itcueoff) then
                        	nspike(k)=nspike(k)+1
                        	sumspe(k)=sumspe(k)+dte(k)
                        	sumsp2e(k)=sumsp2e(k)+dte(k)**2
		sumecv2(k)=sumecv2(k)+ 2*abs(dte(k)-dteold(k))/(dte(k)+dteold(k))
			endif
			dteold(k)=dte(k)
                	tspikeold(k)=tspike(k)
		else 
	  		se(k)=sebis(k)
		end if
	end do


	do k=1,ni
	switchi(k)=0
		if(sibis(k).gt.sit) then
		switchi(k)=1
			dt2=dtf*(sibis(k)-sit)/(sibis(k)-si(k))
	snew=sir+(sibis(k)-sit)*(1+dtf/taui*(-sir+si(k))/(sibis(k)-si(k)))
                        tspiki(k)=it*dtf
                        dti(k)=tspiki(k)-tspikiold(k)
       siim(k)=siim(k)*(1.-biguii)*dexp(-dti(k)/taufacii)+biguii
               xx=dexp(-dti(k)/taurecee)
               siid(k)=siid(k)*(1.-siim(k))*xx+1.-xx
!		fei2(k)=exp(-dt2/tei2)*siim(k)*siid(k)
   		fei2(k)=exp(-dt2/tei2)
!  		fii2(k)=exp(-dt2/tii2)*siim(k)*siid(k)
   		fii2(k)=exp(-dt2/tii2)
!	write(27,*) (it)*dtf*0.01,k
 	call flush(27)
			si(k)=snew
 		simact(k)=simact(k)+epsact*exp(-dt2/tauact)

                        tspiki(k)=it*dtf
                        dti(k)=tspiki(k)-tspikiold(k)
			if(it.gt.itrans .and. it.lt.itcueoff) then
                        	nspiki(k)=nspiki(k)+1
                        	sumspi(k)=sumspi(k)+dti(k)
                        	sumsp2i(k)=sumsp2i(k)+dti(k)**2
	sumicv2(k)=sumicv2(k)+ 2*abs(dti(k)-dtiold(k))/(dti(k)+dtiold(k))
			endif
                        tspikiold(k)=tspiki(k)
			dtiold(k)=dti(k)
		else 
	  		si(k)=sibis(k)
		end if
	end do



	if(mod(it,deltit).eq.10) then

	  ue=0.
	  ui=0.
	  do i=1,ne/2
	  	ue(1)=ue(1)+100*semact(i)
	  	ue(2)=ue(2)+100*semact(i+ne/2)
	  end do	
	  do i=1,ni/2
	  	ui(1)=ui(1)+100*simact(i)
	  	ui(2)=ui(2)+100*simact(i+ni/2)
	  end do	
	ue=ue/float(ne/2)
	ui=ui/float(ni/2)
!	write(3,*) sngl(it*dtf*0.01d0),sngl(ue(1)),sngl(ue(2)),sngl(ui(1)),sngl(ui(2))
!call flush(3)
	if(dabs(ue(1)-ue(2))/(ue(1)+ue(2)).gt.threshdec) then
	timedec=float(it)*dtf*0.01
	if(ue(1)-ue(2).gt.0) then
	nact1=nact1+1
	write(98,*) 1,timedec
	else
	write(98,*) 0,timedec
	end if
	avtimedec=avtimedec+timedec
	goto 999
        end if	
	end if
	end do !it
999	continue

	close(3)

        do i=1,ne
                fire(i)=nspike(i)/float(itcueoff-itrans)/dtf
                if(nspike(i).gt.2) then
                        sumspe(i)=sumspe(i)/nspike(i)
                        sumsp2e(i)=sumsp2e(i)/nspike(i)
                        cve(i)=sqrt(sumsp2e(i)-sumspe(i)**2)/sumspe(i)
			cv2e(i)=sumecv2(i)/nspike(i)
                endif
                if((nspike(i).gt.2)) then
! write(31,*) sngl(i*360.d0/ne),sngl(100*fire(i)),sngl(cve(i)),sngl(cv2e(i))
		else
!        write(31,*) sngl(i*360.d0/ne),sngl(100.d0*fire(i)),0,0
		endif
        end do


        do i=1,ni
                firi(i)=nspiki(i)/float(itcueoff-itrans)/dtf
                if(nspiki(i).gt.2) then
                        sumspi(i)=sumspi(i)/nspiki(i)
                        sumsp2i(i)=sumsp2i(i)/nspiki(i)
                        cvi(i)=sqrt(sumsp2i(i)-sumspi(i)**2)/sumspi(i)
			cv2i(i)=sumicv2(i)/nspiki(i)
                endif
                if((nspiki(i).gt.2)) then
!write(32,*) sngl(i*360.d0/ni),sngl(100.d0*firi(i)),sngl(cvi(i)),sngl(cv2i(i))
		else
!	write(32,*) sngl(i*360.d0/ni),sngl(100.d0*firi(i)),0,0
		endif
        end do
	avfire1=0
	avfire2=0
	avfiri1=0
	avfiri2=0
	do i=1,ne/2
	avfire1=avfire1+fire(i)
	avfire2=avfire2+fire(i+ne/2)
	end do
	do i=1,ni/2
	avfiri1=avfiri1+firi(i)
	avfiri2=avfiri2+firi(i+ni/2)
	end do
	avfire1=2.*avfire1/ne
	avfire2=2.*avfire2/ne
	avfiri1=2.*avfiri1/ni
	avfiri2=2.*avfiri2/ni
!	write(97,*) sngl(epsilon),sngl(dabs(avfire1-avfire2)/(avfire1+avfire2)),sngl((avfire1+avfire2)/2.d0*100.),sngl(dabs(avfiri1-avfiri2)/(avfiri1+avfiri2)),sngl((avfiri1+avfiri2)/2.d0*100.)
	write(12,*) 'itrial=',itrial
	call flush(12)
		end do !itrial
	write(96,*) ntrial,float(nact1)/float(ntrial),avtimedec/ntrial
	end do !iarch


	end 

!==================================================================
 subroutine crmat(sem,n1,n2,sigma,k12,nm,typcon,ncmjj,mee,mei,mie,mii)
        implicit none
        integer n,nc,kk,n1,n2,nce,nci,nm,ne,ni
        real*8 k12
        parameter(n=32000)
        parameter(nc=4000,nce=4000,nci=4000)
        parameter(ne=32000,ni=8000)
        integer ng              !ancho de gg (en sitios)
        integer i,j,il,ncm
        real*8  pp,rnd,ncmm,sumg,ran0
        real*8 sigma,pi,dpi,p12,alpha
        real*8 theta1(n),theta2(n),deltatheta,sigmarad
        real*8 ggubfs,sem
        integer mee(nce,ne),mei(nce,ni),mii(nci,ni),mie(nci,ne)
        integer ncmjj(n),ij
        character*2 typcon

        pi=acos(-1.d0)
        dpi=2.d0*acos(-1.d0)
	alpha=0.5/sigma**2

        do i=1,n1
        theta1(i)=dpi*float(i)/float(n1)
        end do
        do i=1,n2
        theta2(i)=dpi*float(i)/float(n2)
        end do

	do j=1,n2
        ncm=0
        do i=1,n1
        deltatheta=theta1(i)-theta2(j)
        p12=dexp(-alpha*(deltatheta)**2)
        do kk=1,4
        p12=p12+dexp(-alpha*(deltatheta+kk*dpi)**2)+dexp(-alpha*(deltatheta-kk*dpi)**2)
        end do
        p12=p12*k12/dsqrt(dpi*sigma**2)/dfloat(n2)
        p12=k12/dfloat(n2)
        rnd=ggubfs(sem)
        if (rnd.lt.p12) then
                          ncm=ncm+1
	if(typcon.eq.'ee') mee(ncm,j)=i
        if(typcon.eq.'ie') mie(ncm,j)=i
        if(typcon.eq.'ii') mii(ncm,j)=i
        if(typcon.eq.'ei') mei(ncm,j)=i
                        if (ncm.eq.nm) then
          write(6,*) 'error en crmat: ncm=nc',i,typcon
  	  stop
                	end if
        end if
                end do
        ncmjj(j)=ncm
!	write(99,*) j,ncm
        end do
        return
        end
!==================================================================
 subroutine crmatsp(sem,n1,n2,sigma,k12,nm,typcon,ncmjj,mee,mei,mie,mii)
  	implicit none
        integer n,nc,kk,n1,n2,nce,nci,nm,ne,ni
        parameter(n=32000)
        parameter(nc=4000,nce=4000,nci=4000)
        parameter(ne=32000,ni=8000)
        integer ng              !ancho de gg (en sitios)
        integer i,j,il,ncm
        real*8 k12,x12
        real*8  pp,rnd,ncmm,sumg,ran0
        real*8 sigma,pi,dpi,p12,alpha
        real*8 theta1(n),theta2(n),deltatheta,sigmarad
        real*8 ggubfs,sem 
	integer mee(nce,ne),mei(nce,ni),mii(nci,ni),mie(nci,ne)
        integer ncmjj(n),ij
        character*2 typcon

        pi=acos(-1.d0)
        dpi=2.d0*acos(-1.d0)
        alpha=0.5/sigma**2

        do i=1,n1
        theta1(i)=dpi*float(i)/float(n1)
        end do
        do i=1,n2
        theta2(i)=dpi*float(i)/float(n2)
        end do

        do j=1,n2
        ncm=0
        do i=1,n1
!       deltatheta=theta1(i)-theta2(j)
!       p12=dexp(-alpha*(deltatheta)**2)
!       do kk=1,4
!       p12=p12+dexp(-alpha*(deltatheta+kk*dpi)**2)+dexp(-alpha*(deltatheta-kk*dpi)**2)
!       end do
!       p12=p12*sqrt(k12)/dsqrt(dpi*sigma**2)/dfloat(n2)
	x12=0
	if(j.lt.n2/2.and.i.gt.n1/2) x12=1
	if(j.gt.n2/2.and.i.lt.n1/2) x12=1
	p12=sqrt(k12)/dfloat(n2/2)*x12
        rnd=ggubfs(sem)
        if (rnd.lt.p12) then
                          ncm=ncm+1
        if(typcon.eq.'ee') mee(ncm,j)=i
        if(typcon.eq.'ie') mie(ncm,j)=i
        if(typcon.eq.'ii') mii(ncm,j)=i
        if(typcon.eq.'ei') mei(ncm,j)=i
                        if (ncm.eq.nm) then
          write(6,*) 'error en crmat: ncm=nc',i,typcon
          stop
                        end if
        end if
                end do
        ncmjj(j)=ncm
!       write(100,*) j,ncm
        end do
        return
	end

        function ggubfs(seed)
        real*8 ggubfs
        real*8 seed
        real*8 d2p31m,d2p31
        data d2p31m /2147483647.d0/
        data d2p31  /2147483711.d0/

        seed = dmod(16807.d0*seed,d2p31m)
        ggubfs=seed/d2p31
        return
        end

