************************************************************
!
! in home directory, only 8 orientation (PLUS background)
!
! 21 Feb, intracellular tuning data
!
!  7 Feb, threshold w/ background rates
!
!  Analysis of Drifting Grating Simulations
!
!  1) Orientation tuning curve (at fixed contrast)
!    a) mean firing rate as a function of angle
!    b) optimal angle and f1/f0 at optimal angle
!    c) CV (circular variance)
!    d) bandwidth
!
!************************************************************
      MODULE parameters
      parameter ( nex=72 ,ney=120, nix=36, niy=60)
      parameter ( nmax = nex*ney+nix*niy )
      END MODULE 
      PROGRAM analyzeData
      USE parameters
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      dimension nlgni(nmax)
      character (LEN=4) fdr,fdr0
      open(1111,file='analInput_hahaha',status='old',
     1    access='sequential',form='formatted')
      do i=1,nmax
          read(1111,*) nlgni(i)
      enddo
      read(1111,*) n
      read(1111,*) neps
      read(1111,*) nr
      read(1111,*) lgnmax
      read(1111,*) lgnmin
      read(1111,*) iperiod
      if (nr.GE.lgnmax)then
          nr = 1
      endif
      do i=1,neps
          read(1111,*) fdr
          read(fdr,*) eps
          if (i.EQ.1) fdr0 = fdr
          write(*,'(A11,F5.1,A1)') 'contrast = ',eps/10.0,'%'
          print *, 'in', fdr, ' and ', fdr0
          call assemble(fdr0,fdr,n,nlgni,lgnmax,nr,lgnmin,iperiod)
      enddo
      close(1111)
      end
      subroutine assemble(fdr0,fdr,n,nlgni,lgnmax,nr,lgnmin,iperiod)
      USE parameters
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      character*50 fdat(33)
      character*2 thetafdr
      character*4 fdr,fdr0
      character*50 opfdr,theta_wise,samples,fsample,spike_wise,fspike
      dimension nlgni(nmax)
      dimension nspikes(nmax,25,65),nread(nmax*25,65)
      real*4 t(2*n+1)
      dimension fexc(nmax*25,65),finh(nmax*25,65)
      dimension fexc3(nmax,25,65),finh3(nmax,25,65)
      dimension cfexc(nmax*25,65), cfinh(nmax*25,65)
      dimension cfexc3(nmax,25,65),cfinh3(nmax,25,65)
      dimension glgn(nmax*25,65),gexc(nmax*25,65),
     1 ginh(nmax*25,65)
      dimension glgn3(nmax,25,65),gexc3(nmax,25,65),
     1 ginh3(nmax,25,65)
      dimension clgn(nmax*25,65),cexc(nmax*25,65),
     1 cinh(nmax*25,65),cpota(nmax*25,65)
      dimension clgn3(nmax,25,65),cexc3(nmax,25,65),
     1 cinh3(nmax,25,65),cpota3(nmax,25,65)
      
      dimension gtot(nmax*25,65),curr(nmax*25,65),
     1 vslave(nmax*25,65),gpota(nmax*25,65)
      dimension gtot3(nmax,25,65),curr3(nmax,25,65),
     1 vslv3(nmax,25,65),gpota3(nmax,25,65)
      dimension vmem(nmax*25,65),vmem3(nmax,25,65)
      dimension currEt(nmax,2*n+1),currEc(nmax,2*n+1),currIc(nmax,2*n+1)
      dimension currEn(nmax,2*n+1),currIn(nmax,2*n+1)
      dimension currEt2(nmax,2*n+1),currEc2(nmax,2*n+1),
     1 currIc2(nmax,2*n+1)
      dimension currEn2(nmax,2*n+1),currIn2(nmax,2*n+1)
      dimension currp(nmax,2*n+1),currp2(nmax,2*n+1)
      dimension gg(nmax),cc(nmax),vv(nmax),gl(nmax),ge(nmax),gi(nmax)
      dimension fe(nmax),fi(nmax),vs(nmax),gp(nmax),gpavg(nmax,2*n+1)
      dimension gavg(nmax,2*n+1),cavg(nmax,2*n+1),vsavg(nmax,2*n+1)
      dimension glavg(nmax,2*n+1),geavg(nmax,2*n+1),giavg(nmax,2*n+1)
      dimension feavg(nmax,2*n+1),fiavg(nmax,2*n+1),vavg(nmax,2*n+1)
      dimension rates(nmax,2*n+1),frate(nmax,2*n), stifr(nmax,2*n)
      dimension cv2(nmax),theta(2*n),pkrate(nmax),cv1(nmax)
      dimension porate(nmax),sc(nmax),sc2(nmax),fmrate(nmax),
     1 oorate(nmax),pirate(nmax),oirate(nmax)
      dimension glgnsc(nmax,2*n+1),glgnsc2(nmax,2*n+1)
      dimension gexcsc(nmax,2*n+1),gexcsc2(nmax,2*n+1)
      dimension vslavesc(nmax,2*n+1),vslavesc2(nmax,2*n+1)
      dimension po(nmax,iperiod),brate(nmax),oo(nmax,iperiod)
      dimension indpo(nmax),indoo(nmax),indpi(nmax),indoi(nmax)
      dimension excite(nmax)
      complex*16 eye,ftheta,f1,f2,ftheta1,vtheta,vtheta1
      equivalence (nspikes,nread),(gpota3,gpota)
      equivalence (gtot3,gtot), (curr3,curr), (vslv3,vslave)
      equivalence (glgn3,glgn), (gexc3,gexc), (ginh3,ginh)
      equivalence (fexc3,fexc), (finh3,finh), (vmem3,vmem)
      equivalence (clgn3,clgn), (cexc3,cexc), (cinh3,cinh)
      equivalence (cfexc3,cfexc), (cfinh3,cfinh), (cpota,cpota3)
      logical excite1d(nmax)
!------------------------------------------------------------
!        prepare input for spikes
!------------------------------------------------------------
        speriod = dble(iperiod)
        if (iperiod.gt.25) then
          print *, 'pls recomplie with larger
     1             iperiod for equivalence objects'
        stop
        endif
        opfdr = 'processed/'//fdr
        theta_wise = 'processed/theta_wise/'//fdr
        spike_wise = 'processed/spike_wise/'//fdr
        samples = 'processed/samples/'//fdr
        call system('mkdir -p processed')
        call system('mkdir -p '//trim(opfdr))
        call system('mkdir -p '//trim(theta_wise))
        call system('mkdir -p '//trim(spike_wise))
        call system('mkdir -p '//trim(samples))
        do i=1,n 
            write(thetafdr,'(I2.2)') i-1
            fdat(i) = fdr//'/'//thetafdr//'/i-and-f.dat1'
        enddo

        write(thetafdr,'(I2.2)') n
        if (trim(fdr).EQ.trim(fdr0)) then 
            fdat(n+1) = fdr0//'/'//thetafdr//'/i-and-f.dat1'
            call system('mv '//trim(fdat(n+1))//' '
     1      //trim(theta_wise)//'/'//thetafdr//'-i-and-f.dat1')
            print *,fdat(n+1),' moved'
        endif
        theta_wise = 'processed/theta_wise/'//fdr0
        fdat(n+1)=trim(theta_wise)//'/'//thetafdr//'-i-and-f.dat1' 
        theta_wise = 'processed/theta_wise/'//fdr
!------------------------------------------------------------
!        read in exc
!------------------------------------------------------------
      open(10,
     1 file=fdr//'/00/i-and-f.dat3',
     2 status='old',form='unformatted',access='direct',recl=8*nmax)
      read(10,rec=5) excite
      close(10)

      do i=1,nmax
        excite1d(i) = .true.
      if (excite(i).lt.0.5D0) excite1d(i) = .false.
      enddo
!------------------------------------------------------------
!  Read in ALL drifting grating spikes
!------------------------------------------------------------
      do i=1,n
      open(10,file=fdat(i),status='old',form='unformatted',
     1       access='direct',recl=4*iperiod*nmax,iostat=ierr)
      read(10,rec=1) (nread(j,i),j=1,nmax*iperiod)
      close(10)
      open(10,file=fdat(i),status='old',form='unformatted',
     1       access='direct',recl=4)
      read(10,rec=iperiod*nmax+1) t(i)
      if (ierr.NE.0.0) then
          print *,fdat(i),' error occurred '
      endif
      close(10)
        enddo

        do i=1,n
      do j=1,nmax*iperiod
        nread(j,i+n)=nread(j,i)
      enddo
      t(i+n) = t(i)
      enddo
      open(10,file=fdat(n+1),status='old',form='unformatted',
     1       access='direct',recl=4*iperiod*nmax,iostat=ierr)
      read(10,rec=1) (nread(j,2*n+1),j=1,nmax*iperiod)
      if (ierr.NE.0.0) then
          print *,fdat(n+1),' error occurred '
      endif
      close(10)
      open(10,file=fdat(n+1),status='old',form='unformatted',
     1       access='direct',recl=4,iostat=ierr)
      read(10,rec=iperiod*nmax+1) t(2*n+1)
      if (ierr.NE.0.0) then
          print *,fdat(n+1),' error occurred '
      endif
      close(10)
!-----------------------------------------------------------
!    Move to processed folder
!-----------------------------------------------------------
            do i=1,n 
        write(thetafdr,'(I2.2)') i-1
        fdat(i) = fdr//'/'//thetafdr//'/i-and-f.dat1'
        call system('mv '//trim(fdat(i))//' '
     1    //trim(theta_wise)//'/'//thetafdr//'-i-and-f.dat1')
             enddo

         if (trim(fdr).EQ.trim(fdr0)) then 
             n = n + 1
         endif
             do i=1,n
         write(thetafdr,'(I2.2)') i-1

         fdat(i) = fdr//'/'//thetafdr//'/gammaOsci.dat'
         call system('mv '//trim(fdat(i))//' '
     1    //trim(theta_wise)//'/'//thetafdr//'-gammaOsci.dat')

         fsample = fdr//'/'//thetafdr//'/samples.dat'
         call system('mv '//trim(fsample)//' '
     1    //trim(samples)//'/'//thetafdr//'-samples.dat')

         fspike = fdr//'/'//thetafdr//'/spikes.dat'
         call system('mv '//trim(fspike)//' '
     1    //trim(spike_wise)//'/'//thetafdr//'-spikes.dat')

         fspike = fdr//'/'//thetafdr//'/onpoi.dat'
         call system('mv '//trim(fspike)//' '
     1    //trim(spike_wise)//'/'//thetafdr//'-onpoi.dat')

         fspike = fdr//'/'//thetafdr//'/offpoi.dat'
         call system('mv '//trim(fspike)//' '
     1    //trim(spike_wise)//'/'//thetafdr//'-offpoi.dat')

         fspike = fdr//'/'//thetafdr//'/enoi.dat'
         call system('mv '//trim(fspike)//' '
     1    //trim(spike_wise)//'/'//thetafdr//'-enoi.dat')

         fspike = fdr//'/'//thetafdr//'/inoi.dat'
         call system('mv '//trim(fspike)//' '
     1    //trim(spike_wise)//'/'//thetafdr//'-inoi.dat')
            enddo
        if (trim(fdr).EQ.trim(fdr0)) then 
            n = n - 1
        endif
!------------------------------------------------------------
!  Compute mean firing rate for each neuron at each angle
!------------------------------------------------------------
      frmaxE = 0.d0
      frmaxI = 0.d0
        do i=1,nmax
      fmrate(i) = 0.d0
      fmax     = -100.d0
        do j=1,2*n
      nsp = 0
        do k=1,iperiod
          nsp = nsp + nspikes(i,k,j)
        enddo
      frate(i,j) = (1.D0 * nsp) / t(j)
      if (fmax.lt.frate(i,j)) fmax = frate(i,j)
      rates(i,j) = frate(i,j)
      fmrate(i)  = fmrate(i) + frate(i,j)
        enddo
      fmrate(i) = fmrate(i) / (2.D0*n)
      pkrate(i) = fmax
      if(pkrate(i).GT.frmaxE.and.excite1d(i))then
              frmaxE = pkrate(i)
      endif
      if(pkrate(i).GT.frmaxI.and.(.not.excite1d(i)))then
              frmaxI = pkrate(i)
      endif
!------------------------------------------------------------
!  background rate j=(2*n+1)
!------------------------------------------------------------
      nsp = 0
      do k=1,iperiod
        nsp = nsp + nspikes(i,k,2*n+1)
      enddo
      brate(i) = (1.D0 * nsp) / t(2*n+1)
      rates(i,2*n+1) = brate(i)
        enddo
      do i=1,nmax
        do j=1,2*n
      stifr(i,j) = frate(i,j)-brate(i)
      if (stifr(i,j).lt.0.0) stifr(i,j) = 0
        enddo
      enddo
!------------------------------------------------------------
!  Compute CV for each neuron w/ tuning curve
!------------------------------------------------------------
      twopi = 8.D0*datan(1.D0)
      eye   = (0.D0,1.D0)
      do j=1,2*n
        theta(j) = (j-1)*twopi/(2.D0*n)
      enddo
      do i=1,nmax
        fmean = 0.0
        fmean1 = 0.0
        ftheta = (0.D0,0.D0)
        ftheta1 = (0.D0,0.D0)
      do j=1,2*n
        fmean = fmean + frate(i,j)
        fmean1 = fmean1 + stifr(i,j)
        ftheta = ftheta + frate(i,j)*cdexp(2.D0*eye*theta(j))
        ftheta1 = ftheta1 + stifr(i,j)*cdexp(2.D0*eye*theta(j))
      enddo
        if (fmean.gt.1.D-8) then
          cv2(i) = 1.D0 - cdabs(ftheta/fmean)
        else
          cv2(i) = 1.D0
        endif
        if (fmean1.gt.1.D-8) then
          cv1(i) = 1.D0 - cdabs(ftheta1/fmean1)
        else
          cv1(i) = 1.D0
        endif
      enddo
      eqcv = 0.D0
      yqcv = 0.D0
      rmincv = 1.0
      rmaxcv = 0.0
      resmincv = 1.0
      resmaxcv = 0.0
      recmincv = 1.0
      recmaxcv = 0.0
      do i=1,nmax
      if (excite1d(i)) then
         eqcv = cv2(i) + eqcv
      else
         yqcv = cv2(i) + yqcv
      endif
         if (cv2(i).gt.rmaxcv) rmaxcv = cv2(i)
         if (cv2(i).lt.rmincv) rmincv = cv2(i)
      if (nr.GE.lgnmax)then
         if (nlgni(i).LE.lgnmin.and.excite1d(i)) then
            if (cv2(i).gt.recmaxcv) recmaxcv = cv2(i)
             if (cv2(i).lt.recmincv) recmincv = cv2(i)
         endif
      else
         if (nlgni(i).LE.nr.and.excite1d(i)) then
            if (cv2(i).gt.recmaxcv) recmaxcv = cv2(i)
             if (cv2(i).lt.recmincv) recmincv = cv2(i)
         endif
      endif
         if (nlgni(i).GE.lgnmax.and.excite1d(i)) then
            if (cv2(i).gt.resmaxcv) resmaxcv = cv2(i)
             if (cv2(i).lt.resmincv) resmincv = cv2(i)
         endif
      enddo
          print *,'computed iCVs:', sngl(yqcv/(nix*niy))
          print *,'computed eCVs:', sngl(eqcv/(nex*ney))
          print *,'lowest CV:',sngl(rmincv),
     1    ' highest CV:',sngl(rmaxcv)
          print *,'c lowest CV:',sngl(recmincv),
     1    ' highest CV:',sngl(recmaxcv)
          print *,'s lowest CV:',sngl(resmincv),
     1    ' highest CV:',sngl(resmaxcv)
!------------------------------------------------------------
!  Find optimal orientation and output cycle-averaged firing rate at PO
!------------------------------------------------------------
      eye   = (0.D0,-1.D0)
      do i=1,nmax
        fmax = -100.0
        indpo(i) = -1
        do j=1,2*n
          if (frate(i,j).gt.fmax) then
            fmax      = frate(i,j)
            indpo(i) = j
          endif
        enddo
        
      if (indpo(i).gt.n) indpo(i) = indpo(i) - n
      indoo(i) = mod(indpo(i)+n/2-1,n)+1

      f0 = 0.D0
      fo0 = 0.D0
      f1 = (0.D0,0.D0)
      f2 = (0.D0,0.D0)
      do k=1,iperiod
       angle   = (k-1)*twopi/speriod
       po(i,k) = nspikes(i,k,indpo(i))*1.D0/
     1    (t(indpo(i))/speriod)
       oo(i,k) = nspikes(i,k,indoo(i))*1.D0/
     1    (t(indoo(i))/speriod)
       f0      = f0 + po(i,k)
       fo0     = fo0 + oo(i,k)
       f1      = f1 + po(i,k)*cdexp(eye*angle)
       f2      = f2 + po(i,k)*cdexp(2.D0*eye*angle)
      enddo
      porate(i) = f0
      oorate(i) = fo0
      if (f0.gt.1D-8) then
        sc(i)   = cdabs(2*f1/f0)
        sc2(i)  = cdabs(2*f2/f0)
      else
        sc(i)   = 0.0
        sc2(i)  = 0.0
      endif
      enddo
      rep = 0.D0
      reo = 0.D0
      cep = 0.D0
      ceo = 0.D0
      sep = 0.D0
      seo = 0.D0
      nne = 0
      nnec = 0
      nnes = 0
      rip = 0.D0
      rio = 0.D0
      nni = 0
      rebrate = 0.D0
      ribrate = 0.D0
      do i=1,nmax 
         if (excite1d(i)) then 
             rep = rep + porate(i)
             reo = reo + oorate(i)
             rebrate = rebrate + brate(i)
             nne = nne + 1
             if (nlgni(i).EQ.lgnmin)then
                 nnec = nnec+1
                 cep = cep + porate(i)
                 ceo = ceo + oorate(i)
             endif
             if (nlgni(i).EQ.lgnmax) then
                 nnes = nnes+1
                 sep = sep + porate(i)
                 seo = seo + oorate(i)
             endif
         else
             rip = rip + porate(i)
             rio = rio + oorate(i) 
             ribrate = ribrate + brate(i)
             nni = nni + 1
         endif
      enddo
      rep = rep/(nne*speriod)
      reo = reo/(nne*speriod)
      cep = cep/(nnec*speriod)
      ceo = ceo/(nnec*speriod)
      sep = sep/(nnes*speriod)
      seo = seo/(nnes*speriod)
      rip = rip/(nni*speriod)
      rio = rio/(nni*speriod)
      rebrate = rebrate/nne
      ribrate = ribrate/nni
      print *,'computed POs'
      print *,'exc rate:',sngl(rep),sngl(reo),
     1    'max:',sngl(frmaxE),
     2    ' brate:',sngl(rebrate),
     3    ' complex:',sngl(cep),sngl(ceo)
      print *,'inh rate:',sngl(rip),sngl(rio),
     1    'max:',sngl(frmaxI),
     2    ' brate:',sngl(ribrate),
     3    ' simple:',sngl(sep),sngl(seo)
!------------------------------------------------------------
!
!  Dump averages to files
!
!------------------------------------------------------------
      open(11,file=trim(opfdr)//'/cv.dat',
     1    status='replace',form='unformatted',
     2    access='direct',recl=8*nmax)
      write(11,rec=1) cv2
      write(11,rec=2) sc
      write(11,rec=3) sc2
      write(11,rec=4) fmrate
      write(11,rec=5) pkrate
      write(11,rec=6) porate
      write(11,rec=7) ((indpo(i)-1)*twopi/(2.D0*n),i=1,nmax)
      write(11,rec=8) excite
      write(11,rec=9) brate
      write(11,rec=10) cv1
      write(11,rec=11) oorate
      close(11)


      open(11,file=trim(opfdr)//'/frate.dat',
     1    status='replace',form='unformatted',
     2    access='direct',recl=8*nmax)
      do i=1,2*n 
        write(11,rec=i) (frate(j,i),j=1,nmax)
      enddo
      do i=1,2*n 
        write(11,rec=2*n+i) (stifr(j,i),j=1,nmax)
      enddo
      close(11)
!------------------------------------------------------------
!  Now onto intracellular stuff 
!------------------------------------------------------------
      do i=1,n 
          write(thetafdr,'(I2.2)') i-1
          fdat(i) = fdr//'/'//thetafdr//'/i-and-f.dat4'
      enddo

      write(thetafdr,'(I2.2)') n
      if (trim(fdr).EQ.trim(fdr0)) then 
        fdat(n+1) = fdr0//'/'//thetafdr//'/i-and-f.dat4'
        call system('mv '//trim(fdat(n+1))//' '
     1    //trim(theta_wise)//'/'//thetafdr//'-i-and-f.dat4')
      endif
      theta_wise = 'processed/theta_wise/'//fdr0
      fdat(n+1)=trim(theta_wise)//'/'//thetafdr//'-i-and-f.dat4' 
!------------------------------------------------------------
!  Read in ALL drifting grating data (cycled averaged data)
!------------------------------------------------------------
      do i=1,n
        open(10,file=fdat(i),status='old',form='unformatted',
     1       access='direct',recl=8*iperiod*nmax,iostat=ierr)
        read(10,rec=1) (glgn(j,i),j=1,nmax*iperiod)
        read(10,rec=2) (gexc(j,i),j=1,nmax*iperiod)
        read(10,rec=3) (ginh(j,i),j=1,nmax*iperiod)
        read(10,rec=4) (gtot(j,i),j=1,nmax*iperiod)
        read(10,rec=5) (curr(j,i),j=1,nmax*iperiod)
        read(10,rec=6) (vslave(j,i),j=1,nmax*iperiod)
        read(10,rec=7) (vmem(j,i),j=1,nmax*iperiod)
        read(10,rec=8) (fexc(j,i),j=1,nmax*iperiod)
        read(10,rec=9) (finh(j,i),j=1,nmax*iperiod)
        read(10,rec=19) (gpota(j,i),j=1,nmax*iperiod)
        read(10,rec=21) (clgn(j,i),j=1,nmax*iperiod)
        read(10,rec=22) (cexc(j,i),j=1,nmax*iperiod)
        read(10,rec=23) (cinh(j,i),j=1,nmax*iperiod)
        read(10,rec=24) (cfexc(j,i),j=1,nmax*iperiod)
        read(10,rec=25) (cfinh(j,i),j=1,nmax*iperiod)
        read(10,rec=26) (cpota(j,i),j=1,nmax*iperiod)
        close(10)
        if (ierr.NE.0.0) then
            print *,fdat(i),' error occurred '
        endif
      enddo
      i=(2*n+1)
        open(10,file=fdat(n+1),status='old',form='unformatted',
     1       access='direct',recl=8*iperiod*nmax,iostat=ierr)
        read(10,rec=1) (glgn(j,i),j=1,nmax*iperiod)
        read(10,rec=2) (gexc(j,i),j=1,nmax*iperiod)
        read(10,rec=3) (ginh(j,i),j=1,nmax*iperiod)
        read(10,rec=4) (gtot(j,i),j=1,nmax*iperiod)
        read(10,rec=5) (curr(j,i),j=1,nmax*iperiod)
        read(10,rec=6) (vslave(j,i),j=1,nmax*iperiod)
        read(10,rec=7) (vmem(j,i),j=1,nmax*iperiod)
        read(10,rec=8) (fexc(j,i),j=1,nmax*iperiod)
        read(10,rec=9) (finh(j,i),j=1,nmax*iperiod)
        read(10,rec=19) (gpota(j,i),j=1,nmax*iperiod)
        read(10,rec=21) (clgn(j,i),j=1,nmax*iperiod)
        read(10,rec=22) (cexc(j,i),j=1,nmax*iperiod)
        read(10,rec=23) (cinh(j,i),j=1,nmax*iperiod)
        read(10,rec=24) (cfexc(j,i),j=1,nmax*iperiod)
        read(10,rec=25) (cfinh(j,i),j=1,nmax*iperiod)
        read(10,rec=26) (cpota(j,i),j=1,nmax*iperiod)
        close(10)
        if (ierr.NE.0.0) then
            print *,fdat(n+1),' error occurred '
        endif
      do i=1,n
      do j=1,nmax*iperiod
        glgn(j,i+n) = glgn(j,i)
        gexc(j,i+n) = gexc(j,i)
        ginh(j,i+n) = ginh(j,i)
        gtot(j,i+n) = gtot(j,i)
        curr(j,i+n) = curr(j,i)
        vslave(j,i+n) = vslave(j,i)
        vmem(j,i+n) = vmem(j,i)
        fexc(j,i+n) = fexc(j,i)
        finh(j,i+n) = finh(j,i)
        gpota(j,i+n) = gpota(j,i)

        clgn(j,i+n) = clgn(j,i)
        cexc(j,i+n) = cexc(j,i)
        cinh(j,i+n) = cinh(j,i)
        cfexc(j,i+n) = cfexc(j,i)
        cfinh(j,i+n) = cfinh(j,i)
        cpota(j,i+n) = cpota(j,i)
      enddo
      enddo
!------------------------------------------------------------
!        output all orientation's averages over cycles 
!------------------------------------------------------------
      open(11,file=trim(opfdr)//'/cycles.dat',
     1    status='replace',form='unformatted',
     2     access='direct',recl=8*nmax*iperiod)
      do i=1,n
      write(11,rec=17*(i-1)+1) (1.D0*nread(j,i), j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+2) (glgn(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+3) (gexc(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+4) (ginh(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+5) (fexc(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+6) (finh(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+7) (vmem(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+8) (gtot(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+9) (curr(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+10) (vslave(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+11) (gpota(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+12) (clgn(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+13) (cexc(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+14) (cinh(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+15) (cfexc(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+16) (cfinh(j,i),j=1,iperiod*nmax)
      write(11,rec=17*(i-1)+17) (cpota(j,i),j=1,iperiod*nmax)
      enddo
      close(11)
!------------------------------------------------------------
      do i=1,nmax
      do k=1,(2*n+1)
      sumg = 0.0
      sumc = 0.0
      sumvs = 0.0
      sumv = 0.0
      sumgl = 0.0
      sumge = 0.0
      sumgi = 0.0
      sumfe = 0.0
      sumfi = 0.0
      sumgp = 0.0
        do j=1,iperiod
      sumg = sumg + gtot3(i,j,k)
      sumc = sumc + curr3(i,j,k)
      sumvs = sumvs + vslv3(i,j,k)
      sumv  = sumv  + vmem3(i,j,k)
      sumgl = sumgl + glgn3(i,j,k)
      sumge = sumge + gexc3(i,j,k)
      sumgi = sumgi + ginh3(i,j,k)
      sumfe = sumfe + fexc3(i,j,k)
      sumfi = sumfi + finh3(i,j,k)
      sumgp = sumgp + gpota3(i,j,k)
        enddo
      gavg(i,k) = sumg/speriod
      cavg(i,k) = sumc/speriod
      vsavg(i,k) = sumvs/speriod
      glavg(i,k) = sumgl/speriod
      geavg(i,k) = sumge/speriod
      giavg(i,k) = sumgi/speriod
      feavg(i,k) = sumfe/speriod
      fiavg(i,k) = sumfi/speriod
      vavg(i,k) = sumv/speriod
      gpavg(i,k) = sumgp/speriod
      enddo
      enddo
      do i=1,nmax
      sumg = 0.0
      sumc = 0.0
      sumvs = 0.0
      sumgl = 0.0
      sumge = 0.0
      sumgi = 0.0
      sumfe = 0.0
      sumfi = 0.0
      sumv = 0.0
      sumgp = 0.0
          do k=1,2*n
      sumg = sumg + gavg(i,k)
      sumc = sumc + cavg(i,k)
      sumvs = sumvs + vsavg(i,k)
      sumgl = sumgl + glavg(i,k)
      sumge = sumge + geavg(i,k)
      sumgi = sumgi + giavg(i,k)
      sumfe = sumfe + feavg(i,k)
      sumfi = sumfi + fiavg(i,k)
      sumv = sumv + vavg(i,k)
      sumgp = sumgp + gpavg(i,k)
          enddo
      gg(i) = sumg/(2.D0*n)
      cc(i) = sumc/(2.D0*n)
      vs(i) = sumvs/(2.D0*n)
      gl(i) = sumgl/(2.D0*n)
      ge(i) = sumge/(2.D0*n)
      gi(i) = sumgi/(2.D0*n)
      fe(i) = sumfe/(2.D0*n)
      fi(i) = sumfi/(2.D0*n)
      vv(i) = sumv/(2.D0*n)
      gp(i) = sumgp/(2.D0*n)
      enddo
!------------------------------------------------------------
!
!  Find F1/F0 of gLGN, gE, Veff of all angles
!
!------------------------------------------------------------
        do j = 1,2*n+1
      do i=1,nmax
        f0 = 0.D0
        f1 = (0.D0,0.D0)
        f2 = (0.D0,0.D0)
        do k=1,iperiod
          angle   = (k-1)*twopi/speriod
          po(i,k) = glgn3(i,k,j)
          f0      = f0 + po(i,k)
          f1      = f1 + po(i,k)*cdexp(eye*angle)
          f2      = f2 + po(i,k)*cdexp(2.D0*eye*angle)
        enddo
        if (f0.gt.1D-8) then
          glgnsc(i,j)   = cdabs(2*f1/f0)
          glgnsc2(i,j)  = cdabs(2*f2/f0)
        else
          glgnsc(i,j)   = 0.0
          glgnsc2(i,j)  = 0.0
        endif
      enddo
        enddo

       	do i = 1,nmax
        fmax = -100.0
        indpi(i) = -1
        do j=1,2*n
          fff = glavg(i,j)*(1+glgnsc(i,j))
          if (fff.gt.fmax) then
            fmax      = fff
            indpi(i) = j
          endif
        enddo
          if (indpi(i).gt.n) indpi(i) = indpi(i) - n
          indoi(i) = mod(indpi(i)+n/2-1,n)+1
        enddo
      open(11,file=trim(opfdr)//'/cv.dat',
     1    status='old',form='unformatted',
     2    access='direct',recl=8*nmax)
      write(11,rec=12) ((indpi(i)-1)*twopi/(2.D0*n),i=1,nmax)
      close(11)

      Print *,'computed F1/F0 of gLGN'

        do j = 1,2*n+1
      do i = 1,nmax
        f0 = 0.D0
        f1 = (0.D0,0.D0)
        f2 = (0.D0,0.D0)
        do k=1,iperiod
          angle   = (k-1)*twopi/speriod
          po(i,k) = gexc3(i,k,j)
          f0      = f0 + po(i,k)
          f1      = f1 + po(i,k)*cdexp(eye*angle)
          f2      = f2 + po(i,k)*cdexp(2.D0*eye*angle)
        enddo
        if (f0.gt.1D-8) then
          gexcsc(i,j)   = cdabs(2*f1/f0)
          gexcsc2(i,j)  = cdabs(2*f2/f0)
        else
          gexcsc(i,j)   = 0.0
          gexcsc2(i,j)  = 0.0
        endif
      enddo
        enddo
      Print *,'computed F1/F0 of gE'

        do j = 1,2*n+1
      do i=1,nmax
        f0 = 0.D0
        f1 = (0.D0,0.D0)
        f2 = (0.D0,0.D0)
        do k=1,iperiod
          angle   = (k-1)*twopi/speriod
          po(i,k) = vslv3(i,k,j)
          f0      = f0 + po(i,k)
          f1      = f1 + po(i,k)*cdexp(eye*angle)
          f2      = f2 + po(i,k)*cdexp(2.D0*eye*angle)
        enddo
        if (f0.gt.1D-8) then
          vslavesc(i,j)   = cdabs(2*f1/f0)
          vslavesc2(i,j)  = cdabs(2*f2/f0)
        else
          vslavesc(i,j)   = 0.0
          vslavesc2(i,j)  = 0.0
        endif
      enddo
        enddo
      Print *,'computed F1/F0 of V_s'

!----------------avgs----------------------------
      open(11,file=trim(opfdr)//'/intra.dat',
     1    status='replace',form='unformatted',
     2    access='direct',recl=8*nmax)
      do i=1,2*n+1
        write(11,rec=i) (gavg(j,i),j=1,nmax)
      enddo
      do i=1,2*n+1
        write(11,rec=2*n+1+i)  (cavg(j,i),j=1,nmax) 
      enddo
      do i=1,2*n+1
          write(11,rec=2*(2*n+1)+i)  (vsavg(j,i),j=1,nmax)
      enddo
      do i=1,2*n+1              
          write(11,rec=3*(2*n+1)+i)  (glavg(j,i),j=1,nmax)  
      enddo                    
      do i=1,2*n+1
          write(11,rec=4*(2*n+1)+i)  (geavg(j,i),j=1,nmax)
      enddo                    
      do i=1,2*n+1             
          write(11,rec=5*(2*n+1)+i)  (giavg(j,i),j=1,nmax)  
      enddo                    
      do i=1,2*n+1             
          write(11,rec=6*(2*n+1)+i)  (feavg(j,i),j=1,nmax)  
      enddo                    
      do i=1,2*n+1             
          write(11,rec=7*(2*n+1)+i)  (fiavg(j,i),j=1,nmax)  
      enddo                    
      do i=1,2*n+1             
          write(11,rec=8*(2*n+1)+i)  (rates(j,i),j=1,nmax)  
      enddo                     
      do i=1,2*n+1              
          write(11,rec=9*(2*n+1)+i)  (vavg(j,i), j=1,nmax)   
      enddo
      do i=1,2*n+1              
          write(11,rec=10*(2*n+1)+i) (glgnsc(j,i),j=1,nmax)   
      enddo
      do i=1,2*n+1              
          write(11,rec=11*(2*n+1)+i) (glgnsc2(j,i),j=1,nmax)   
      enddo
      do i=1,2*n+1              
          write(11,rec=12*(2*n+1)+i) (gexcsc(j,i),j=1,nmax)   
      enddo
      do i=1,2*n+1              
          write(11,rec=13*(2*n+1)+i) (gexcsc2(j,i),j=1,nmax)   
      enddo
      do i=1,2*n+1              
          write(11,rec=14*(2*n+1)+i) (vslavesc(j,i),j=1,nmax)   
      enddo
      do i=1,2*n+1              
          write(11,rec=15*(2*n+1)+i) (vslavesc2(j,i),j=1,nmax)   
      enddo
      do i=1,2*n+1              
          write(11,rec=16*(2*n+1)+i) (gpavg(j,i),j=1,nmax)   
      enddo
      close(11)

      open(11,file=trim(opfdr)//'/intraavg.dat',
     1    status='replace',form='unformatted',
     2    access='direct',recl=8*nmax)
      write(11,rec=1) gg
      write(11,rec=2) cc
      write(11,rec=3) vs
      write(11,rec=4) gl
      write(11,rec=5) ge
      write(11,rec=6) gi
      write(11,rec=7) fe
      write(11,rec=8) fi
      write(11,rec=9) (gavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=10) (cavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=11) (vsavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=12) (glavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=13) (geavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=14) (giavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=15) (feavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=16) (fiavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=17) gp 
      write(11,rec=18) (gpavg(i,(2*n+1)),i=1,nmax)
!------------------------------------------------------------
!
!  Find F1/F0 of vslave
!
!------------------------------------------------------------
      do i=1,nmax
         sc(i)   = vslavesc(i,indpo(i));
         sc2(i)  = vslavesc2(i,indpo(i));
      enddo

      write(11,rec=17) sc
      write(11,rec=18) sc2
      write(11,rec=19) vv
      write(11,rec=20) (vavg(i,(2*n+1)),i=1,nmax)
!------------------------------------------------------------
      close(11)
!------------------------------------------------------------
!        read in second order information 
!------------------------------------------------------------
      do i=1,n
        open(10,file=fdat(i),status='old',form='unformatted',
     1       access='direct',recl=8*iperiod*nmax,iostat=ierr)
        read(10,rec=10) (glgn(j,i),j=1,nmax*iperiod)
        read(10,rec=11) (gexc(j,i),j=1,nmax*iperiod)
        read(10,rec=12) (ginh(j,i),j=1,nmax*iperiod)
        read(10,rec=13) (gtot(j,i),j=1,nmax*iperiod)
        read(10,rec=14) (curr(j,i),j=1,nmax*iperiod)
        read(10,rec=15) (vslave(j,i),j=1,nmax*iperiod)
        read(10,rec=16) (vmem(j,i),j=1,nmax*iperiod)
        read(10,rec=17) (fexc(j,i),j=1,nmax*iperiod)
        read(10,rec=18) (finh(j,i),j=1,nmax*iperiod)
        read(10,rec=20) (gpota(j,i),j=1,nmax*iperiod)
        read(10,rec=27) (clgn(j,i),j=1,nmax*iperiod)
        read(10,rec=28) (cexc(j,i),j=1,nmax*iperiod)
        read(10,rec=29) (cinh(j,i),j=1,nmax*iperiod)
        read(10,rec=30) (cfexc(j,i),j=1,nmax*iperiod)
        read(10,rec=31) (cfinh(j,i),j=1,nmax*iperiod)
        read(10,rec=32) (cpota(j,i),j=1,nmax*iperiod)
        close(10)
        if (ierr.NE.0.0) then
            print *,fdat(i),' error occurred '
        endif
      enddo
      i=(2*n+1)
        open(10,file=fdat(n+1),status='old',form='unformatted',
     1       access='direct',recl=8*iperiod*nmax,iostat=ierr)
        read(10,rec=10) (glgn(j,i),j=1,nmax*iperiod)
        read(10,rec=11) (gexc(j,i),j=1,nmax*iperiod)
        read(10,rec=12) (ginh(j,i),j=1,nmax*iperiod)
        read(10,rec=13) (gtot(j,i),j=1,nmax*iperiod)
        read(10,rec=14) (curr(j,i),j=1,nmax*iperiod)
        read(10,rec=15) (vslave(j,i),j=1,nmax*iperiod)
        read(10,rec=16) (vmem(j,i),j=1,nmax*iperiod)
        read(10,rec=17) (fexc(j,i),j=1,nmax*iperiod)
        read(10,rec=18) (finh(j,i),j=1,nmax*iperiod)
        read(10,rec=20) (gpota(j,i),j=1,nmax*iperiod)
        read(10,rec=27) (clgn(j,i),j=1,nmax*iperiod)
        read(10,rec=28) (cexc(j,i),j=1,nmax*iperiod)
        read(10,rec=29) (cinh(j,i),j=1,nmax*iperiod)
        read(10,rec=30) (cfexc(j,i),j=1,nmax*iperiod)
        read(10,rec=31) (cfinh(j,i),j=1,nmax*iperiod)
        read(10,rec=32) (cpota(j,i),j=1,nmax*iperiod)
        close(10)
        if (ierr.NE.0.0) then
            print *,fdat(n+1),' error occurred '
        endif
      do i=1,n
      do j=1,nmax*iperiod
      glgn(j,i+n) = glgn(j,i)
      gexc(j,i+n) = gexc(j,i)
      ginh(j,i+n) = ginh(j,i)
      gtot(j,i+n) = gtot(j,i)
      curr(j,i+n) = curr(j,i)
      vslave(j,i+n) = vslave(j,i)
      vmem(j,i+n) = vmem(j,i)
      fexc(j,i+n) = fexc(j,i)
      finh(j,i+n) = finh(j,i)
      gpota(j,i+n) = gpota(j,i)

      clgn(j,i+n) = clgn(j,i)
      cexc(j,i+n) = cexc(j,i)
      cinh(j,i+n) = cinh(j,i)
      cfexc(j,i+n) = cfexc(j,i)
      cfinh(j,i+n) = cfinh(j,i)
      cpota(j,i+n) = cpota(j,i)
      enddo
      enddo
!-----------------------------------------------------------
!    Move to processed folder
!--------------------------------------------------------------
      theta_wise = 'processed/theta_wise/'//fdr
      do i=1,n 
          write(thetafdr,'(I2.2)') i-1
          fdat(i) = fdr//'/'//thetafdr//'/i-and-f.dat4'
        call system('mv '//trim(fdat(i))//' '
     1    //trim(theta_wise)//'/'//thetafdr//'-i-and-f.dat4')
      enddo
!------------------------------------------------------------
!        output all orientation 2nd moment over cycles 
!------------------------------------------------------------
      open(11,file=trim(opfdr)//'/cycles.dat',
     1    status='old',form='unformatted',
     2    access='direct',recl=8*nmax*iperiod)
      do i=1,n
      write(11,rec=17*n+(i-1)*16+1) (glgn(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+(i-1)*16+2) (gexc(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+(i-1)*16+3) (ginh(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+(i-1)*16+4) (fexc(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+(i-1)*16+5) (finh(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+(i-1)*16+6) (vmem(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+(i-1)*16+7) (gtot(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+(i-1)*16+8) (curr(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+(i-1)*16+9) (vslave(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+(i-1)*16+10) (gpota(j,i),j=1,nmax*iperiod)
      write(11,rec=17*n+16*(i-1)+11) (clgn(j,i),j=1,iperiod*nmax)
      write(11,rec=17*n+16*(i-1)+12) (cexc(j,i),j=1,iperiod*nmax)
      write(11,rec=17*n+16*(i-1)+13) (cinh(j,i),j=1,iperiod*nmax)
      write(11,rec=17*n+16*(i-1)+14) (cfexc(j,i),j=1,iperiod*nmax)
      write(11,rec=17*n+16*(i-1)+15) (cfinh(j,i),j=1,iperiod*nmax)
      write(11,rec=17*n+16*(i-1)+16) (cpota(j,i),j=1,iperiod*nmax)
      enddo
      close(11)
!------------------------------------------------------------
!     compute average gtot, curr, vslave
!------------------------------------------------------------
      do i=1,nmax
      do k=1,(2*n+1)
        sumg = 0.0
        sumc = 0.0
        sumv = 0.0
        sumgl = 0.0
        sumge = 0.0
        sumgi = 0.0
        sumfe = 0.0
        sumfi = 0.0
        sumvs = 0.0
        sumgp = 0.0
        do j=1,iperiod
      sumg = sumg + gtot3(i,j,k)
      sumc = sumc + curr3(i,j,k)
      sumvs = sumvs + vslv3(i,j,k)
      sumgl = sumgl + glgn3(i,j,k)
      sumge = sumge + gexc3(i,j,k)
      sumgi = sumgi + ginh3(i,j,k)
      sumfe = sumfe + fexc3(i,j,k)
      sumfi = sumfi + finh3(i,j,k)
      sumv = sumv + vmem3(i,j,k)
      sumgp = sumgp + gpota3(i,j,k)
        enddo
      gavg(i,k) = sumg/speriod
      cavg(i,k) = sumc/speriod
      vsavg(i,k) = sumvs/speriod
      glavg(i,k) = sumgl/speriod
      geavg(i,k) = sumge/speriod
      giavg(i,k) = sumgi/speriod
      feavg(i,k) = sumfe/speriod
      fiavg(i,k) = sumfi/speriod
      vavg(i,k) = sumv/speriod
      gpavg(i,k) = sumgp/speriod
      enddo
      enddo
        do i=1,nmax
      sumg = 0.0
      sumc = 0.0
      sumvs = 0.0
      sumgl = 0.0
      sumge = 0.0
      sumgi = 0.0
      sumfe = 0.0
      sumfi = 0.0
      sumv = 0.0
      sump = 0.0
          do k=1,2*n
      sumg = sumg + gavg(i,k)
      sumc = sumc + cavg(i,k)
      sumvs = sumvs + vsavg(i,k)
      sumgl = sumgl + glavg(i,k)
      sumge = sumge + geavg(i,k)
      sumgi = sumgi + giavg(i,k)
      sumfe = sumfe + feavg(i,k)
      sumfi = sumfi + fiavg(i,k)
      sumv = sumv + vavg(i,k)
      sump = sump + gpavg(i,k)
          enddo
      gg(i) = sumg/(2.D0*n)
      cc(i) = sumc/(2.D0*n)
      vs(i) = sumvs/(2.D0*n)
      gl(i) = sumgl/(2.D0*n)
      ge(i) = sumge/(2.D0*n)
      gi(i) = sumgi/(2.D0*n)
      fe(i) = sumfe/(2.D0*n)
      fi(i) = sumfi/(2.D0*n)
      vv(i) = sumv/(2.D0*n)
      gp(i) = sumgp/(2.D0*n)
      enddo

      open(11,file=trim(opfdr)//'/intra2.dat',
     1    status='replace',form='unformatted',
     2    access='direct',recl=8*nmax)
      do i=1,2*n+1
        write(11,rec=i) (gavg(j,i),j=1,nmax)
      enddo
      do i=1,2*n+1
        write(11,rec=2*n+1+i)  (cavg(j,i),j=1,nmax) 
      enddo                     
      do i=1,2*n+1              
          write(11,rec=2*(2*n+1)+i)  (vsavg(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=3*(2*n+1)+i)  (glavg(j,i),j=1,nmax)  
      enddo                    
      do i=1,2*n+1             
          write(11,rec=4*(2*n+1)+i)  (geavg(j,i),j=1,nmax)  
      enddo                    
      do i=1,2*n+1             
          write(11,rec=5*(2*n+1)+i)  (giavg(j,i),j=1,nmax)  
      enddo                    
      do i=1,2*n+1             
          write(11,rec=6*(2*n+1)+i)  (feavg(j,i),j=1,nmax)  
      enddo                    
      do i=1,2*n+1             
          write(11,rec=7*(2*n+1)+i)  (fiavg(j,i),j=1,nmax)  
      enddo                    
      do i=1,2*n+1             
          write(11,rec=8*(2*n+1)+i)  (rates(j,i),j=1,nmax)  
      enddo                     
      do i=1,2*n+1              
          write(11,rec=9*(2*n+1)+i)  (vavg(j,i),j=1,nmax)
      enddo
      do i=1,2*n+1              
          write(11,rec=10*(2*n+1)+i)  (gpavg(j,i),j=1,nmax)
      enddo
      close(11)

      open(11,file=trim(opfdr)//'/intraavg2.dat',
     1    status='replace',form='unformatted',
     2    access='direct',recl=8*nmax)
      write(11,rec=1) gg
      write(11,rec=2) cc
      write(11,rec=3) vs
      write(11,rec=4) gl
      write(11,rec=5) ge
      write(11,rec=6) gi
      write(11,rec=7) fe
      write(11,rec=8) fi
      write(11,rec=9) (gavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=10) (cavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=11) (vsavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=12) (glavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=13) (geavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=14) (giavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=15) (feavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=16) (fiavg(i,(2*n+1)),i=1,nmax)
      write(11,rec=17) gp 
      write(11,rec=18) (gpavg(i,(2*n+1)),i=1,nmax)
      close(11)
!======================================================
!   read in current data merge to one file
!======================================================
      do i=1,n 
          write(thetafdr,'(I2.2)') i-1
          fdat(i) = fdr//'/'//thetafdr//'/i-and-f.dat2'
      enddo

      write(thetafdr,'(I2.2)') n
      if (trim(fdr).EQ.trim(fdr0)) then 
        fdat(n+1) = fdr0//'/'//thetafdr//'/i-and-f.dat2'
        call system('mv '//trim(fdat(n+1))//' '
     1    //trim(theta_wise)//'/'//thetafdr//'-i-and-f.dat2')
      endif
      theta_wise = 'processed/theta_wise/'//fdr0
      fdat(n+1)=trim(theta_wise)//'/'//thetafdr//'-i-and-f.dat2' 
!=============================================================
      do i=1,n
        open(10,file=fdat(i),status='old',form='unformatted',
     1       access='direct',recl=8*nmax,iostat=ierr)
          read(10,rec=1) (currEt(j,i),j=1,nmax)
          read(10,rec=2) (currEc(j,i),j=1,nmax)
          read(10,rec=3) (currIc(j,i),j=1,nmax)
          read(10,rec=4) (currEn(j,i),j=1,nmax)
          read(10,rec=5) (currIn(j,i),j=1,nmax)
          read(10,rec=6) (currEt2(j,i),j=1,nmax)
          read(10,rec=7) (currEc2(j,i),j=1,nmax)
          read(10,rec=8) (currIc2(j,i),j=1,nmax)
          read(10,rec=9) (currEn2(j,i),j=1,nmax)
          read(10,rec=10) (currIn2(j,i),j=1,nmax)
          read(10,rec=11) (currp(j,i),j=1,nmax)
          read(10,rec=12) (currp2(j,i),j=1,nmax)
        close(10)
        if (ierr.NE.0.0) then
            print *,fdat(i),' error occurred '
        endif
      enddo
      i=(2*n+1)
        open(10,file=fdat(n+1),status='old',form='unformatted',
     1       access='direct',recl=8*nmax,iostat=ierr)
          read(10,rec=1) (currEt(j,i),j=1,nmax)
          read(10,rec=2) (currEc(j,i),j=1,nmax)
          read(10,rec=3) (currIc(j,i),j=1,nmax)
          read(10,rec=4) (currEn(j,i),j=1,nmax)
          read(10,rec=5) (currIn(j,i),j=1,nmax)
          read(10,rec=6) (currEt2(j,i),j=1,nmax)
          read(10,rec=7) (currEc2(j,i),j=1,nmax)
          read(10,rec=8) (currIc2(j,i),j=1,nmax)
          read(10,rec=9) (currEn2(j,i),j=1,nmax)
          read(10,rec=10) (currIn2(j,i),j=1,nmax)
          read(10,rec=11) (currp(j,i),j=1,nmax)
          read(10,rec=12) (currp2(j,i),j=1,nmax)
        close(10)
        if (ierr.NE.0.0) then
            print *,fdat(i),' error occurred '
        endif
      do i=1,n
        do j=1,nmax
            currEt(j,i+n) = currEt(j,i)
            currEc(j,i+n) = currEc(j,i)
            currIc(j,i+n) = currIc(j,i)
            currEn(j,i+n) = currEn(j,i)
            currIn(j,i+n) = currIn(j,i)
            currp(j,i+n) = currp(j,i)
            currEt2(j,i+n) = currEt2(j,i)
            currEc2(j,i+n) = currEc2(j,i)
            currIc2(j,i+n) = currIc2(j,i)
            currEn2(j,i+n) = currEn2(j,i)
            currIn2(j,i+n) = currIn2(j,i)
            currp2(j,i+n) = currp2(j,i)
        enddo
      enddo
!======================================================
!    move raw data 
!======================================================
      theta_wise = 'processed/theta_wise/'//fdr
      do i=1,n 
        write(thetafdr,'(I2.2)') i-1
        fdat(i) = fdr//'/'//thetafdr//'/i-and-f.dat2'
        call system('mv '//trim(fdat(i))//' '
     1    //trim(theta_wise)//'/'//thetafdr//'-i-and-f.dat2')
      enddo
!======================================================
!   write to file
!======================================================
      open(11,file=trim(opfdr)//'/curr.dat',
     1    status='replace',form='unformatted',
     2    access='direct',recl=8*nmax,iostat=ierr)
      do i=1,2*n+1
        write(11,rec=i) (currEt(j,i),j=1,nmax)
      enddo
      do i=1,2*n+1
          write(11,rec=1*(2*n+1)+i)  (currEc(j,i),j=1,nmax) 
      enddo                     
      do i=1,2*n+1              
          write(11,rec=2*(2*n+1)+i)  (currIc(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=3*(2*n+1)+i)  (currEn(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=4*(2*n+1)+i)  (currIn(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=5*(2*n+1)+i)  (currEt2(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=6*(2*n+1)+i)  (currEc2(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=7*(2*n+1)+i)  (currIc2(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=8*(2*n+1)+i)  (currEn2(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=9*(2*n+1)+i)  (currIn2(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=10*(2*n+1)+i)  (currp(j,i),j=1,nmax)   
      enddo                     
      do i=1,2*n+1              
          write(11,rec=11*(2*n+1)+i)  (currp2(j,i),j=1,nmax)   
      enddo                     
      if (ierr.NE.0.0) then
        print *,'curr.dat error occurred'
      endif
      close(11)
!============== prefer orthogonal in out preference info===

      open(11,file=trim(opfdr)//'/cv.dat',
     1    status='old',form='unformatted',
     2    access='direct',recl=8*nmax)
      write(11,rec=13) 1.D0*indpo
      write(11,rec=14) 1.D0*indoo
      write(11,rec=15) 1.D0*indpi
      write(11,rec=16) 1.D0*indoi
      close(11)
!======================================================
      return
      END 
!************************************************************
