!************************************************************
!
!  Integrate-and-Fire Model LGN -> V1 
!
!************************************************************
      PROGRAM ompa
      USE parameters
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension pex(ney,nex), pey(ney,nex),pix(niy,nix),piy(niy,nix)
      dimension see(nmax),sei(nmax),sie(nmax),sii(nmax)
      dimension gl(nmax),thetas(nmax),nSub(nmax)
      dimension v(nmax),vnew(nmax),alpha1(nmax),beta1(nmax)
      dimension glo(nlgn), slo1(nlgn), slo2(nlgn), slo3(nlgn)
      dimension glon(nlgn),slon1(nlgn),slon2(nlgn),slon3(nlgn)
      dimension glf(nlgn), slf1(nlgn), slf2(nlgn), slf3(nlgn)
      dimension glfn(nlgn),slfn1(nlgn),slfn2(nlgn),slfn3(nlgn)
      dimension gl2o(nlgn), sl2o1(nlgn), sl2o2(nlgn), sl2o3(nlgn)
      dimension gl2on(nlgn),sl2on1(nlgn),sl2on2(nlgn),sl2on3(nlgn)
      dimension gl2f(nlgn), sl2f1(nlgn), sl2f2(nlgn), sl2f3(nlgn)
      dimension gl2fn(nlgn),sl2fn1(nlgn),sl2fn2(nlgn),sl2fn3(nlgn)
      dimension gx(nmax),sx1(nmax),sx2(nmax),sx3(nmax)
      dimension gy(nmax),sy1(nmax),sy2(nmax),sy3(nmax)
      dimension gn(nmax),sn1(nmax),sn2(nmax),sn3(nmax)
      dimension gm(nmax),sm1(nmax),sm2(nmax),sm3(nmax)
      dimension g2x(nmax),s2x1(nmax),s2x2(nmax),s2x3(nmax)
      dimension g2y(nmax),s2y1(nmax),s2y2(nmax),s2y3(nmax)
      dimension g2n(nmax),s2n1(nmax),s2n2(nmax),s2n3(nmax)
      dimension g2m(nmax),s2m1(nmax),s2m2(nmax),s2m3(nmax)
      dimension gp(nmax),sp1(nmax),sp2(nmax),sp3(nmax)
      dimension currEt(nmax,max_ncaj),currEt2(nmax,max_ncaj)
      dimension currEc(nmax,max_ncaj),currEc2(nmax,max_ncaj)
      dimension currIc(nmax,max_ncaj),currIc2(nmax,max_ncaj)
      dimension currEn(nmax,max_ncaj),currEn2(nmax,max_ncaj)
      dimension currIn(nmax,max_ncaj),currIn2(nmax,max_ncaj)
      dimension currp(nmax,max_ncaj),currp2(nmax,max_ncaj)
      dimension nca(nmax,max_ncaj),sumnca(nmax)
      dimension tcurrEc(nmax),tcurrIc(nmax),tcurrEn(nmax)
      dimension tcurrEt(nmax),tcurrIn(nmax),tcurrp(nmax)
      dimension tcurrEc2(nmax),tcurrIc2(nmax),tcurrEn2(nmax)
      dimension tcurrEt2(nmax),tcurrIn2(nmax),tcurrp2(nmax)
      dimension ge(nmax),se1(nmax),se2(nmax),se3(nmax)
      dimension gf(nmax),sf1(nmax),sf2(nmax),sf3(nmax)
      dimension gi(nmax),si1(nmax),si2(nmax),si3(nmax)
      dimension gj(nmax),sj1(nmax),sj2(nmax),sj3(nmax)
      dimension slist(256,nmax)
      dimension aa_ee(nmax),aa_ei(nmax),aa_ie(nmax),aa_ii(nmax)
      dimension xlgn(nlgn),ylgn(nlgn)
      dimension nonlgn(nmax),noflgn(nmax),nlgni(nmax)
      dimension ionlgn(nmax,nlgn),ioflgn(nmax,nlgn)
      dimension sonlgn(nmax,nlgn),soflgn(nmax,nlgn)
      dimension in(nmax)
      dimension pspon(nlgn),pspoff(nlgn),pspike(nmax)
      dimension irate(nmax,25)
      dimension glgn(nmax,25),gexc(nmax,25),ginh(nmax,25),gtot(nmax,25)
      dimension cond(nmax,25),vslave(nmax,25),gpota(nmax,25)
      dimension glgn2(nmax,25),gexc2(nmax,25),ginh2(nmax,25)
      dimension gtot2(nmax,25),cond2(nmax,25),vslave2(nmax,25)
      dimension vmem(nmax,25),vmem2(nmax,25),gpota2(nmax,25)
      dimension gnexc(nmax,25),gninh(nmax,25)
      dimension gnexc2(nmax,25),gninh2(nmax,25)

      dimension clgn(nmax,25),cexc(nmax,25),cinh(nmax,25)
      dimension cpota(nmax,25),cpota2(nmax,25)
      dimension clgn2(nmax,25),cexc2(nmax,25),cinh2(nmax,25)
      dimension cnexc(nmax,25),cninh(nmax,25)
      dimension cnexc2(nmax,25),cninh2(nmax,25)

      dimension tspike(nmax),ispike(nmax)
      dimension tsonpoi(nlgn),ionpoi(nlgn)
      dimension tsoffpoi(nlgn),ioffpoi(nlgn)
      dimension tsenoi(nmax),ienoi(nmax)
      dimension tsinoi(nmax),iinoi(nmax)
      dimension exc(nmax)
      dimension IpostSyn(nmax,nmax)
      dimension NpostSynE(nmax), NpostSynI(nmax)
      dimension NpreSynE(nmax), NpreSynI(nmax)
      dimension profile(256),Istrength(nmax,nmax)
      dimension gexcitMax(nex*ney),ginhibMax(nex*ney),erate(nex*ney)
      dimension gll(nex*ney),rpdf(6),npdf(7)
      logical excite(nmax),sample(nmax),recordAllspikes
      character(LEN=6), dimension(10) :: label
      character (LEN=25) fence0,fence1
      character (LEN=9 ) finput
      character (LEN=10) flgn
      character (LEN=11) fSprofile 
      character (LEN=30) fsp,fn,fonpoi,foffpoi,fenoi,finoi,fn_id
      character (LEN=30) fspn,fonpoin,foffpoin,fenoin,finoin
      character (LEN=4 ) epsfdr
      character (LEN=2 ) thetafdr
      character (LEN=20) f1,f2,f3,f4,fc
      character (LEN=21) fo
      character (LEN=20) str,fcm
      character (LEN=32) arg
      real*8 lgnsatur,gexcitMax0,gexcitrMax,glrate,unitary
      real*8 ginhibMax0,ginhibrMax,erateEGmax,glgE
      INTEGER*4 OMP_GET_THREAD_NUM,nsp,rk,ntrans
      INTEGER*4 nlgnrateMax,nlgngEmax
      INTEGER*4 nPreEGMax,nPreIGMax
      INTEGER*4 nPreERMax,nPreIRMax,idap
      data iword / 8 /
      external chain
!------------------------------------------------------------
!
!  For each V1 cell ---
!   1) see, sei, sie, sii: intracortical coupling strengths
!   2) gl: summed LGN conductance (summing over all on & off LGN cells
!   3) nonlgn, noflgn: number of On and Off LGN cells sending afferents
!   4) ionlgn, ioflgn: index of On and Off LGN cells sending afferents
!    5) sonlgn, soflgn: strength of On and Off LGN cells sending afferents
!
!  For each LGN cell
!   xlgn, ylgn: coordinates of the center of its "receptive field"
!
!------------------------------------------------------------
      common / PDF / rpdf, npdf
      common / npos / pex, pey, pix, piy
      common / smatrx / see,sei,sie,sii,spotaE,spotaI
      common / lgncnd / gl
      common / lgnmap / sonlgn,soflgn,ionlgn,ioflgn,
     1   nonlgn, noflgn
      common / lgnpos / xlgn,ylgn,jlgn
!------------------------------------------------------------
! ================ p means pre dt value ===================
!
!  glo, glon: On-centered LGN cells, AMPA & NMDA channels
!  glf, glfn: Off-centered LGN cells, AMPA & NMDA channels
!
!------------------------------------------------------------
      common / chaino / glo, slo1, slo2, slo3
      common / chaino2/ glon,slon1,slon2,slon3
      common / chainf / glf, slf1, slf2, slf3
      common / chainf2/ glfn,slfn1,slfn2,slfn3
      common / chain2o / gl2o, sl2o1, sl2o2, sl2o3
      common / chain2o2/ gl2on,sl2on1,sl2on2,sl2on3
      common / chain2f / gl2f, sl2f1, sl2f2, sl2f3
      common / chain2f2/ gl2fn,sl2fn1,sl2fn2,sl2fn3
!------------------------------------------------------------
!
!  ge, gf: intracortical excitation, AMPA & NMDA, resp.
!  gi, gj: intracortical inhibition, GABA_A & B, resp.
!
!------------------------------------------------------------
      common / chaine / ge,se1,se2,se3
      common / chaine2/ gf,sf1,sf2,sf3
      common / chaini / gi,si1,si2,si3
      common / chainj / gj,sj1,sj2,sj3
!------------------------------------------------------------
!
!  gx, gy: "background" excitation, AMPA & NMDA
!  gn, gm: "background" inhibition, GABA_A & B
!
!------------------------------------------------------------
      common / chainx / gx,sx1,sx2,sx3
      common / chainy / gy,sy1,sy2,sy3
      common / chainn / gn,sn1,sn2,sn3
      common / chainm / gm,sm1,sm2,sm3
      common / chain2x / g2x,s2x1,s2x2,s2x3
      common / chain2y / g2y,s2y1,s2y2,s2y3
      common / chain2n / g2n,s2n1,s2n2,s2n3
      common / chain2m / g2m,s2m1,s2m2,s2m3
!------------------------------------------------------------
!
!  gp: potassium calcium activated channel (AHP)
!
!------------------------------------------------------------
      common / chainp / gp,sp1,sp2,sp3
!------------------------------------------------------------
!  Various reversal potentials and time constants
!------------------------------------------------------------
      common / vconst / vthres,vreset,vexcit,vinhib,gleak,gleakI
      common / tconst / tau_e,tau_i,tau0,tau1,tau2,tnrise,tndamp,
     1 tau_a0, tau_a1
!------------------------------------------------------------
!  fnmdat(c) Fraction of NMDA for thalamocortical (& intracortical)
!  fgaba     Fraction of GABA_B
!------------------------------------------------------------
      common /  NMDA  / fnmdatE,fnmdacE,fnmdatI,fnmdacI,
     1   fnmdanE, fnmdanI, fgaba
!------------------------------------------------------------
!  a_xy: spatial kernels
!  excite:  = 1 for excitatory neurons, = 0 for inhibitory
!  nlgni:  No. of LGN afferents
!------------------------------------------------------------
      common / neuron / excite,nlgni,trefE,trefI,nr,
     1   lgnmax,lgnmin,meanlgn
!------------------------------------------------------------
!  Stimulus parameters:
!    omega: temporal frequency
!    (gkx,gky): grating direction vector
!    gphi: grating spatial phase
!
!------------------------------------------------------------
      common / lgnIndivid / pspon,pspoff,gkx,gky,eps,satur,gtheta
      common / lgnSat / xS,yS,linear,rI0
      common /   rhs  / alpha1,beta1
      common /  avgs  / condavg,curravg,condlgn,condexc,condinh,
     1   geavg,giavg,gpavg
      common /  avgsi / vsuma,vsumb,vsuml,vsumen,vsumin,vsume,vsumi,
     1   vsump
      common / curr / currEc,currIc,currEt,currEn,currIn,currp,nca,ncaj
      common / curr2 / currEc2,currIc2,currEt2,currEn2,currIn2,
     1  currp2,sumnca
      common / tcurr / tcurrEc,tcurrIc,tcurrEt,tcurrEn,tcurrIn,
     1  tcurrp
      common / tcurr2 / tcurrEc2,tcurrIc2,tcurrEt2,tcurrEn2,tcurrIn2,
     1  tcurrp2
      common /   isi  / pspike,nsptot,isptot,jsptot,ksptot,msptot
      common / spikes / tspike,ispike,nspike
      common / onpoi / tsonpoi,ionpoi,nonpoi
      common / offpoi / tsoffpoi,ioffpoi,noffpoi
      common / enoi / tsenoi,ienoi,nenoi
      common / inoi / tsinoi,iinoi,ninoi
      common / persp /  irate
      common / ttotal / tfinal,twindow,rstart
      common / connec / synfailEE,synfailEI,synfailIE,synfailII,
     1   neffEE,neffEI,neffIE,neffII
      common / jumpsize / aa_ee,aa_ei,aa_ie,aa_ii,profile,
     1   slist
      common / post / NpostSynE, NpostSynI, IpostSyn,
     1   NpreSynE, NpreSynI,Istrength
      common / poststatics/ nEE, nEI, nIE, nII
      common / prestatics/ npEE, npEI, npIE, npII
      common / lgnGain / tmpphi,conAmp,g0,frtlgn0,treflgn,
     1   gphi, tstart, gfail,cond0
      common / lgnConvol / gk2, taulgn, omega
      common / files / f1,f2,f3,f4,fn,fo,thetafdr,str
      common / noises / ce0_e,ci0_e,frtinh_e,frtexc_e,
     1   ce0_i,ci0_i,frtinh_i,frtexc_i
      common / counters/ i,ii,iii,j,jj,ij,ijj,ncount
      common / stats / glgn,gexc,ginh,gtot,cond,vslave,vmem,gnexc,
     1   gninh,gpota
      common / stats2 / glgn2,gexc2,ginh2,gtot2,cond2,vslave2,
     1   vmem2,gnexc2,gninh2, gpota2
      common / ctats / clgn,cexc,cinh,cnexc,cninh,cpota
      common / ctats2 / clgn2,cexc2,cinh2,cnexc2,cninh2,cpota2
      common / volt / v, vnew
      common / ggIc / gIIc, gEIc
      common / nums / nv1e,nv1s,nv1c,nv1m,nv1i
      common / fspikes / fsp,fonpoi,foffpoi,fenoi,finoi,label
      common / fnspikes / fspn,fonpoin,foffpoin,fenoin,finoin
      common / indRate / erate
      common / EIFconst / DeltaT,vTheta,phiVAS
!$OMP THREADPRIVATE(
!$OMP& /chaino/, /chaino2/, /chainf/, /chainf2/,/ctats/,
!$OMP& /chainx/, /chainy/,  /chainn/, /chainm/,/ctats2/,
!$OMP& /chain2o/,/chain2o2/,/chain2f/,/chain2f2/,
!$OMP& /chain2x/,/chain2y/, /chain2n/,/chain2m/, /chainp/,
!$OMP& /chaine/,/chaine2/,/chaini/,/chainj/,/rhs/,/lgnIndivid/,
!$OMP& /lgncnd/,/spikes/,/isi/,/avgs/,/avgsi/,/files/,/volt/,
!$OMP& /fspikes/,/stats/,/stats2/,/onpoi/,/offpoi/,/enoi/,
!$OMP& /inoi/,/fnspikes/,/persp/,/indRate/)
!------------------------------------------------------------
!
!  Some (not all) INPUT & OUTPUT file declarations
!
!------------------------------------------------------------
      call getarg(1, arg)
      read(arg,'(I1)') neps
      print *, neps
      iword2 = iword/2
      finput = 'adg_INPUT'
      flgn   = 'lgnmap.out'
      fcm = 'conMat.mat'
      fSprofile = 'profile.out'
      fence1 = '-------------------------'
      fence0 = '========================='
!------------------------------------------------------------
!
! Read INPUT file for run parameters
!
!------------------------------------------------------------
      open(unit=13,file=finput,access='sequential',form='formatted')
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*) dt,tfinal,tstep1,tstep2,tstep3,tstep4,iperiod,iseed
      read(13,*)
      read(13,*)
      read(13,*) vthres,vreset,vexcit,vinhib,gleak,DeltaT,vTheta,gleakI
      read(13,*)
      read(13,*)
      read(13,*) tau_e,tau_i,tnrise,tndamp,tau2,treflgn
      read(13,*)
      read(13,*)
      read(13,*) dE,dI,denexc,axnexc,deninh,axninh,gIIc,gEIc
      read(13,*)
      read(13,*)
      read(13,*) ignore,g0,gfail,rI0,tau0,tau1,rall
      read(13,*)
      read(13,*)
      read(13,*) omega,gk,ntheta,gphi,tstart,twindow,taulgn,xS,yS
      read(13,*)
      read(13,*)
      read(13,*)seemax,siemax,ceemax,ciemax,seimax,siimax,ceimax,ciimax,
     1    unitary,nProfile,cr
      read(13,*)
      read(13,*)
      read(13,*) rk, frtlgn0, jlgn, rstart, linear, nsample, trans
      read(13,*)
      read(13,*)
      read(13,*) si4all, se4all, ce4all,trefI,trefE,nr,sigo,retroInh
      read(13,*)
      read(13,*)
      read(13,*) neffII,synfailII,neffIE,synfailIE, readCon,idap
      read(13,*)
      read(13,*)
      read(13,*) neffEI,synfailEI,neffEE,synfailEE, spotaE,spotaI
      read(13,*)
      read(13,*)
      read(13,*) fnmdatE,fnmdacE,fnmdanE,fnmdatI,fnmdacI,fnmdanI,fgaba
      read(13,*)
      read(13,*)
      read(13,*) tau_e4,tau_i4,tnrise4,tndamp4,tau04,tau14,tau24,
     1  tau_a0, tau_a1
      read(13,*)
      read(13,*)
      read(13,*) frtexc_e,ce0_e,frtinh_e,ci0_e,frtexc_i,ce0_i,frtinh_i,
     1  ci0_i
      read(13,*)
      read(13,*)
      read(13,*) SiRatio, SeRatio

      close(13)
      phiVAS = DeltaT*dexp((vreset-vthres)/DeltaT)
      vreset = vreset - phiVAS
      phiVAS = DeltaT*dexp((vreset-vthres)/DeltaT)
      ntrans = abs(trans/dt)
      if (rk.NE.2.and.rk.NE.4) stop 'rk need to be 2 or 4'
      if (rk.eq.4) then
        tau_e = tau_e4
        tau_i = tau_i4
        tnrise = tnrise4
        tndamp = tndamp4
        tau2 = tau24
        tau0 = tau04
        tau1 = tau14
      endif
      if (rall.gt.0.5D0) then
      recordAllspikes = .true.
      else
      recordAllspikes = .false.
      endif
      do i=1,ney
          do j=1,nex
            pex(i,j) = dE/2.D0 + (j-1) * dE
            pey(i,j) = dE/2.D0 + (i-1) * dE
           enddo
      enddo
      do i=1,niy
          do j=1,nix
            pix(i,j) = dI/2.D0 + (j-1) * dI
            piy(i,j) = dI/2.D0 + (i-1) * dI
          enddo
      enddo
      print *,' Inh spans', pix(niy,nix),piy(niy,nix)
      print *,' Exc spans', pex(ney,nex),pey(ney,nex)
      dtheta = atan(1.D0)*4.D0/(ntheta*3)
      if (si4all.GE.0.D0) then
          seimax = si4all
          siimax = si4all
          ciimax = si4all
          ceimax = si4all
      endif
      if (se4all.GE.0.D0) then
          seemax = se4all
          siemax = se4all
      endif
      if (ce4all.GE.0.D0) then
          ceemax = ce4all
          ciemax = ce4all
      endif
      seemax = seemax*SeRatio
      ceemax = ceemax*SeRatio
      siemax = siemax*SeRatio
      ciemax = ciemax*SeRatio
      siimax = siimax*SiRatio
      ciimax = ciimax*SiRatio
      seimax = seimax*SiRatio
      ceimax = ceimax*SiRatio
      if (synfailEI.GT.1.D0) synfailEI = 1.D0 
      if (synfailEE.GT.1.D0) synfailEE = 1.D0 
      if (synfailIE.GT.1.D0) synfailIE = 1.D0 
      if (synfailII.GT.1.D0) synfailII = 1.D0 
!------------------------------------------------------------
!     Read in LGN location and LGN-V1 map from file flgn
!------------------------------------------------------------
      if (rstart.GE.tfinal) stop 'rstart need to be smaller than'//
     1   'tfinal'
      open(14,file=flgn,access='sequential',form='formatted')
      print *,' Reading LGN locations'
      do i=1,9
         read(14,*) 
      enddo
      do i=1,jlgn
         read(14,*) xlgn(i),ylgn(i)
      enddo
      print *,' Reading LGN-V1 map'

!----------------------E/I location can be random
      call e_or_i(excite,exc,iseed)
!-------------------------------------------------
      nv1e = nex*ney
      nv1i = nix*niy
      nv1c = 0
      nv1s = 0
      nv1m = 0
      lgnmin = 1000
      lgnmax = 0
      meanlgn = 0
      do i=1,nmax
      read(14,*) nSub(i), thetas(i)
      nonlgn(i) = 0
      noflgn(i) = 0
      do j=1,nSub(i)
         read(14,*) nonoff
         if (nonoff.gt.0) then
            do k=1,nonoff
                jj = nonlgn(i)+k
                read(14,*) ionlgn(i,jj), sonlgn(i,jj)
            enddo
            nonlgn(i) = nonlgn(i) + nonoff
         else
            do k=1,-nonoff
                jj = noflgn(i)+k
                read(14,*) ioflgn(i,jj), soflgn(i,jj)
                ioflgn(i,jj) = -ioflgn(i,jj)
            enddo
            noflgn(i) = noflgn(i) - nonoff
         endif
      enddo
         nlgni(i) = nonlgn(i) + noflgn(i)
         if(nlgni(i).LT.lgnmin.and.excite(i)) then
            lgnmin=nlgni(i)
         endif
         if(nlgni(i).GT.lgnmax.and.excite(i)) then
            lgnmax=nlgni(i)
         endif
         if (excite(i)) then
            meanlgn = meanlgn + nlgni(i)
         endif
      enddo
      meanlgn = nint(meanlgn / (1.D0*nex*ney))
      if (nr.LT.lgnmin) nr=lgnmin
      do i=1,nmax
         if (nlgni(i).EQ.lgnmin.and.excite(i)) nv1c = nv1c + 1
         if (nlgni(i).EQ.lgnmax.and.excite(i)) nv1s = nv1s + 1
         if (nlgni(i).EQ.meanlgn.and.excite(i))nv1m = nv1m + 1
      enddo
      close(14)
      print *,'lgnmax: ',lgnmax,'lgnmin: ',lgnmin
      print *,'meanlgn: ',meanlgn
!------------------------------------------------------------
!     Read in profile
!------------------------------------------------------------
      open(14,file=fSprofile,access='sequential',form='formatted')
      if (unitary.eq.0.5) then
        do i=1,nProfile
          read(14,*) profile(i)
        enddo
        profile(nProfile) = profile(nProfile)*cr
      endif
      if (dabs(unitary-0.6D0).lt.1e-14)then
        print *,'using different spread of EPSP'
        amaxslist=0.D0
        do j=1,nex*ney
         do i=1,nProfile
           read(14,*) slist(i,j)
      if ((slist(i,j).gt.amaxslist).and.
     1  (dabs(slist(i,j)-1.D0).gt.1e-14))then
          amaxslist=slist(i,j)
      endif
         enddo
        enddo
        do j=nex*ney+1,nmax
            do i=1,nProfile
                slist(i,j) = 1.D0
            enddo
        enddo
      endif
      close(14)
!------------------------------------------------------------
!     Initialize spatial constants and spatial kernels
!------------------------------------------------------------
      HWHM = 2.D0*log(2.D0)
      alee2  =  (denexc*denexc + axnexc*axnexc)/HWHM ! axn: axon
      alei2  =  (denexc*denexc + axninh*axninh)/HWHM ! den: dendrites
      alie2  =  (deninh*deninh + axnexc*axnexc)/HWHM
      alii2  =  (deninh*deninh + axninh*axninh)/HWHM
      twopi  =  8.D0*atan(1.D0)
      dpi    =  0.5D0*twopi
!-------------------------------------------------------------
!  cond0 is normalized to maximum g0
!-------------------------------------------------------------
      satur =  lgnsatur(frtlgn0)
      print *, satur
      cond0  =  g0/satur/tau_e/(1.D0*meanlgn)/gfail
!-------------------------------------------------------------
! linear S between nlgn = 0 & nlgn = 30
!-------------------------------------------------------------

!      if (unitary.EQ.0.5D0) then
!        do i=1,nmax
!            see(i) = 1.D0
!            sei(i) = ceimax
!            sie(i) = ciemax
!            sii(i) = ciimax
!          enddo
!      endif

      do i=1,nmax
          rlambda = (nlgni(i)-lgnmin)/(1.D0*(lgnmax-lgnmin)) 
          if (rlambda.GT.1.D0) rlambda = 1.D0
          see(i) = ceemax + (seemax-ceemax)*rlambda
          sei(i) = ceimax + (seimax-ceimax)*rlambda
          sie(i) = ciemax + (siemax-ciemax)*rlambda
          sii(i) = ciimax + (siimax-ciimax)*rlambda
      enddo
!------------------------------------------------------------
!
!  Now pick 100 neurons at random to keep conductances/potentials
!
!------------------------------------------------------------
      fn_id = 'sample_id.dat' 
      do i=1,nmax
        sample(i) = .false.
      enddo
      open(17,
     1   file=fn_id,status='replace',form='unformatted',
     2   access='direct',recl=iword2)
      if (nsample.LT.nmax) then
        j = 0
        do while (j.LT.nsample)
          i = int(ran2(iseed) * nmax)+1
          if (.not.sample(i)) then
            sample(i) = .true.
              j = j + 1
          endif    
        enddo
        j = 1
        do i = 1,nmax
            if (sample(i)) then
                in(j) = i
                write(17,rec=j) in(j)
                if (j.ge.nsample) then
                    exit
                else
                    j = j + 1
                endif
            endif
        enddo
        if (nsample .eq. 1) then 
            print *,'1 sample neuron: '
        else
            print *,nsample, ' sample neurons : '
        endif
        print *,(in(j),j=1,nsample)
      else
        do i=1,nmax
            in(i)=i
            write(17,rec=i) i
        enddo
        print *,'All neurons are recorded'
      endif
      close(17)
!---------------------Synaptic Failure with fixed connection
      fc = 'cMatrix_summary'
        if (readCon.gt.0.D0) then
      print *,'read in connections'
      call readConnection(fcm,fc,excite)
      print *,'connections read'
        else
      print *,'Form Connections'
      call connection(dE,dI,alee2,alei2,alie2,alii2,thetas,sigo,
     1    retroInh,nlgni,iseed,nr,excite,fc,lgnmax)
      print *,'finished connecting'
        endif
      print *,'--------------------------------------------'
      print *,'    LIF Network of ',ney,' x ',nex,'Exc Neurons'
      print *,'--------------------------------------------'
      print *,' ',niy,' x ',nix,'Inh Neurons'
      print *,'--------------------------------------------'
      print *,'         leak  = ',sngl(gleak)
      print *,'       Vthres  = ',sngl(vthres),'     Ve = ',sngl(vexcit)
      print *,'       Vreset  = ',sngl(vreset),'     Vi = ',sngl(vinhib)
      print *,'--------------------------------------------'
      print *,'        spatial coupling : '
      print *,'        den_e  = ',sngl(denexc),'  den_i = ',sngl(deninh)
      print *,'        axn_e  = ',sngl(axnexc),'  axn_i = ',sngl(axninh)
      print *,'--------------------------------------------'
      print *,' synaptic time constants : '
      print *,'        tau_e  = ',sngl(tau_e), 'tau_0 = ',sngl(tau0)
      print *,'        tau_i  = ',sngl(tau_i), 'tau_1 = ',sngl(tau1)
      print *,'        tau_i2 = ',sngl(tau2)
      if (idap.eq.1)then
      print *,'        tau_a0 = ',sngl(tau_a0), 'tau_a1 = ',sngl(tau_a1)
      endif
      print *,'--------------------------------------------'
      print *,'      synaptic strengths : '
      print *,'          See   = ',sngl(seemax),
     1        '    Sie = ',sngl(siemax)
      print *,'          Sei   = ',sngl(seimax),
     1        '    Sii = ',sngl(siimax)
      print *,'          Cee   = ',sngl(ceemax),
     1        '    Cie = ',sngl(ciemax)
      print *,'          Cei   = ',sngl(ceimax),
     1        '    Cii = ',sngl(ciimax)
      if (idap.eq.1)then
      print *,'    spotaE = ',sngl(spotaE), 'spotaI = ',sngl(spotaI)
      endif
      print *,'--------------------------------------------'
      print *,'        noise parameters : '
      print *,' poisson rate  = ',sngl(frtexc_e),'   str = ',sngl(ce0_e)
      print *,'      (inhib)  = ',sngl(frtinh_e),'   str = ',sngl(ci0_e)
      print *,' poisson rate  = ',sngl(frtexc_i),'   str = ',sngl(ce0_i)
      print *,'      (inhib)  = ',sngl(frtinh_i),'   str = ',sngl(ci0_i)
      print *,'--------------------------------------------'
      print *,'      grating parameters : '
      print *,'          freq  = ',sngl(omega), '      k = ',sngl(gk)
      print *,'       #angle  = ',ntheta,'  phase = ',sngl(gphi)
      print *,'       tstart  = ',sngl(tstart),
     1        'twindow = ',sngl(twindow)
      print *,'--------------------------------------------'
      print *,'output and time-stepping : '
      print *,'       tfinal  = ',sngl(tfinal),'     dt = ',sngl(dt)
      print *,'       tstep1  = ',sngl(tstep1),' tstep2 = ',sngl(tstep2)
      print *,'       tstep3  = ',sngl(tstep3),' tstep4 = ',sngl(tstep4)
      print *,'============================================'
      print *, 'Synaptic Failure EE EI:', sngl(synfailEE),
     1  sngl(synfailEI)
      print *, 'Synaptic Failure IE II:', sngl(synfailIE),
     1  sngl(synfailII)
      if (readCon.lt.0.5D0) then
      print *, 'prescribed clean pre Connections:'
      print *, 'EE EI:', neffEE,', ',neffEI
      print *, 'IE II:', neffIE,', ',neffII
      endif
      print *, 'Pre connection:'
      rneffee = (1.D0*npEE)/nv1e
      rneffei = (1.D0*npEI)/nv1e
      rneffie = (1.D0*npIE)/nv1i
      rneffii = (1.D0*npII)/nv1i
      print *, 'nEE: ', sngl(rneffee), 'nEI: ',
     1  sngl(rneffei)
      print *, 'nIE: ', sngl(rneffie), 'nII: ',
     1  sngl(rneffii)
      if (unitary.eq.0.D0)then
        if (synfailEE.GT.0.D0) then
      print *, 'EE jump(max) = ', ceemax/(rneffee*synfailEE)/tau_e
        endif
        if (synfailEI.GT.0.D0) then
      print *, 'EI jump(max) = ', ceimax/(rneffei*synfailEI)/tau_i
        endif
        if (synfailIE.GT.0.D0) then
      print *, 'IE jump(max) = ', ciemax/(rneffie*synfailIE)/tau_e
        endif
        if (synfailII.GT.0.D0) then
      print *, 'II jump(max) = ', ciimax/(rneffii*synfailII)/tau_i
        endif
      endif    

      if (unitary.gt.0.5D0.and.unitary.lt.1.D0) then
        if (synfailEE.GT.0.D0) then
        if (unitary.eq.0.5D0) then
        print *, 'EE jump profile:'
      do i=2,nProfile
        print *, sngl(profile(i)/tau_e)
      enddo
        endif
        if (dabs(unitary-0.6D0).lt.1e-14)then
            print *, 'max EE jump:'
            print *,sngl(amaxslist/tau_e)
        endif
        endif
        if (synfailEI.GT.0.D0) then
      print *, 'EI jump = ', ceimax/tau_i
        endif
        if (synfailIE.GT.0.D0) then
      print *, 'IE jump = ', ciemax/tau_e
        endif
        if (synfailII.GT.0.D0) then
      print *, 'II jump = ', ciimax/tau_i
        endif
      endif

      if (unitary.eq.1.D0)then
        if (synfailEE.GT.0.D0) then
      print *, 'EE jump(max) = ', ceemax/tau_e
        endif
        if (synfailEI.GT.0.D0) then
      print *, 'EI jump(max) = ', ceimax/tau_i
        endif
        if (synfailIE.GT.0.D0) then
      print *, 'IE jump(max) = ', ciemax/tau_e
        endif
        if (synfailII.GT.0.D0) then
      print *, 'II jump(max) = ', ciimax/tau_i
        endif
      endif
      print *, 'Post connection:'
      rneffee = (1.D0*nEE)/nv1e
      rneffei = (1.D0*nEI)/nv1i
      rneffie = (1.D0*nIE)/nv1e
      rneffii = (1.D0*nII)/nv1i
      print *, 'nEE: ', sngl(rneffee), 'nEI: ',
     1  sngl(rneffei)
      print *, 'nIE: ', sngl(rneffie), 'nII: ',
     1  sngl(rneffii)
      print *, 'cond0 = ',sngl(cond0)
      print *, 'sLGN = ',sngl(cond0*tau_e)
      print *,'--------------------------------------------'
      print *,'  LGN-driving parameters : '
      print *,'  #contrast  = ',neps,'     g0 = ',sngl(g0)
      satur =  lgnsatur(frtlgn0)
      print *,'  spontaneous rate = ',sngl(satur), 
     1  '    gfail = ',sngl(gfail)
      print *,' increase  effective #LGN by ', sngl(rI0)
      print *,'--------------------------------------------'
      count =  0.D0
      do i=1,nmax
      count = count + exc(i)
      enddo
      print *,'No. of Exc/Inh Cells : ',int(count),nmax-int(count)
      print *, 'g0 = ',sngl(g0)
      print *, 'gIIc = ',sngl(gIIc),' gEIc = ',sngl(gEIc)
      print *,'      tauNMDA = ',sngl(tnrise),sngl(tndamp)
      print *,' fNMDA tc/cc/n E= ',sngl(fnmdatE),sngl(fnmdacE),
     1      sngl(fnmdanE)
      print *,' fNMDA tc/cc/n I= ',sngl(fnmdatI),sngl(fnmdacI),
     1      sngl(fnmdanI)
      print *,'  fgaba = ',sngl(fgaba)
      print *,' min lgn input', lgnmin
      print *, '# nv1c = ', nv1c
      print *, '# nv1s = ', nv1s
      print *, '# nv1m = ', nv1m
      ntotal = nint(tfinal/dt)
      nstep1 = nint(tstep1/dt)
      nstep2 = nint(tstep2/dt)
      nstep3 = nint(tstep3/dt)
      nstep4 = nint(tstep4/dt)
      period = 1.D0/omega
      tstep = period/dt/iperiod
      print *, iperiod, 'period parts, ',sngl(tstep),' steps per part'
      nstep0 = nint(tstep)
      nsteperiod = iperiod*nstep0
      steperiod = dble(nsteperiod)
      if (tstart .lt. period) tstart = period 
      do i = 1,nmax
          if (unitary.ge.0.5D0) then
        if (excite(i)) then
            aa_ee(i) = 1.D0/tau_e
            aa_ei(i) = 1.D0/tau_i
        else
              aa_ie(i) = 1.D0/tau_e
              aa_ii(i) = 1.D0/tau_i
        endif
          endif

          if (unitary.eq.0.D0) then
        if (excite(i)) then
            aa_ee(i) = 1.D0/(NpreSynE(i)*synfailEE)/tau_e
            aa_ei(i) = 1.D0/(NpreSynI(i)*synfailEI)/tau_i
        else
              aa_ie(i) = 1.D0/(NpreSynE(i)*synfailIE)/tau_e
              aa_ii(i) = 1.D0/(NpreSynI(i)*synfailII)/tau_i
        endif
          endif
      enddo
      if (synfailEE.LE.0.D0)then
          do i = 1,nmax
              if (excite(i)) then
                aa_ee(i) = 0.D0
              endif
          enddo
      endif
      if (synfailEI.LE.0.D0)then
          do i = 1,nmax
              if (excite(i)) then
                aa_ei(i) = 0.D0
              endif
          enddo
      endif
      if (synfailIE.LE.0.D0)then
          do i = 1,nmax
              if (.not.excite(i)) then
                aa_ie(i) = 0.D0
              endif
          enddo
      endif
      if (synfailII.LE.0.D0)then
          do i = 1,nmax
              if (.not.excite(i)) then
                aa_ii(i) = 0.D0
              endif
          enddo
      endif
!------------------------------------------------------------
!  Input omega in Hertz
!------------------------------------------------------------
      omega  = twopi*omega
      gk2 = (twopi*gk)**2
!!----------------------Determine DRIFTING GRATING time course
      call lgnrf(tmpphi,conAmp)
!============================================================
!    Contrast Loop
!============================================================
!    initialize parallelization
      call OMP_SET_DYNAMIC(.FALSE.)
!------------------------------------------------------------
      do ineps = neps,neps
          epss = 0.125D0*(2**(ineps-1))
          write(epsfdr, '(I4.4)') INT(epss*1000)
          call system('mkdir -p '//TRIM(epsfdr))
          nntheta = ntheta
        call OMP_SET_NUM_THREADS(nntheta)
          print *,' EPSS = ', epss
!============================================================
!    DG Angle Loop
!============================================================
!    parallel do start
!$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE)
!$OMP& SHARED(/smatrx/,/lgnmap/,/lgnpos/,/vconst/,/tconst/,
!$OMP& /NMDA/,/neuron/,/ttotal/,/jumpsize/,/nums/,/connec/,
!$OMP& /post/,/lgnGain/,/noises/,/lgnConvol/,/ggIc/,epss,exc,
!$OMP& epsfdr,gk,twopi,ntheta,iword,ceimax,ceemax,nstep0,dt,
!$OMP& ciimax,ciemax,siimax,siemax,seemax,seimax,nntheta,tist1,
!$OMP& iword2,ntotal,in,axninh,axnexc,deninh,denexc,rk,dtheta,
!$OMP& nsample,recordAllspikes,iperiod,nsteperiod,nc,ntrans,
!$OMP& steperiod,period,ineps,rneffei,rneffie,rneffee,rneffii,
!$OMP& /prestatics/,/EIFconst/,fence0,fence1,unitary,idap)
!$OMP& PRIVATE(gei,gii,gen,gin,t,rctot,vEMean,vIMean,vEMax,
!$OMP& rstot,rmtot,ritot,rtot,vsi,nsp,vpop,gexcitMax,ginhibMax,
!$OMP& /curr/,/curr2/,/tcurr/,/tcurr2/,gexcitMax0,erateMax,
!$OMP& ncycle,/counters/,iid,intheta,ginhibMax0,ginhibrMax,
!$OMP& nlgngEmax,nPreERMax,nPreIRMax,gexcitrMax,erateEGMax,
!$OMP& nlgnrateMax,nPreEGMax,nPreIGMax,glrate,glgE,gll,/PDF/)
!$OMP& FIRSTPRIVATE(iseed)
!------------------------------------------------------------
      do intheta = 1,nntheta
          iid = OMP_GET_THREAD_NUM()
          iseed = iseed + iid
          gtheta = dtheta*(intheta-1)
          eps = epss
          print *, iid ,':',sngl(eps),' angle ',intheta,':',
     1    sngl(conAmp)
!!------------------------------------------------------------
!!  Grating parameters: Drift frequency, spatial frequency
!!------------------------------------------------------------
      gkx    = gk*twopi*cos(gtheta)
      gky    = gk*twopi*sin(gtheta)
!------------------------------------------------------------
!     Initialize output files
!------------------------------------------------------------
      write(thetafdr, '(I2.2)') (intheta-1)
      str = trim(epsfdr)//'/'//trim(thetafdr)    
      call system('mkdir -p '//TRIM(str))
!print *, intheta, ': folder '//str//' created'
!  histogram
      f1     = trim(str)//'/i-and-f.dat1'
!  average current input
      f2     = trim(str)//'/i-and-f.dat2'
!  neuron properties
      f3     = trim(str)//'/i-and-f.dat3' 
!  cycled averages
      f4     = trim(str)//'/i-and-f.dat4'
!  gamma oscillation
      fo     = trim(str)//'/gammaOsci.dat'
!  recording window
      fsp     = trim(str)//'/spikes.dat'
      fonpoi     = trim(str)//'/onpoi.dat'
      foffpoi     = trim(str)//'/offpoi.dat'
      fenoi     = trim(str)//'/enoi.dat'
      finoi     = trim(str)//'/inoi.dat'
!  number of spikes
      fspn     = trim(str)//'/nspikes.dat'
      fonpoin     = trim(str)//'/nonpoi.dat'
      foffpoin     = trim(str)//'/noffpoi.dat'
      fenoin     = trim(str)//'/nenoi.dat'
      finoin     = trim(str)//'/ninoi.dat'

      nsp = 0
      open(intheta+6,
     1  file=fo,status='replace',form='unformatted',
     2  access='direct',recl=iword2)
      close(intheta+6)

      open(intheta+6,
     1  file=fsp,status='replace',form='unformatted',
     2  access='direct',recl=iword2)
      close(intheta+6)

      open(intheta+6,
     1  file=fspn,status='replace',form='unformatted',
     2  access='direct',recl=iword2)
      write(intheta+6,rec=1) nsp
      close(intheta+6)

      open(intheta+6,
     1  file=fonpoi,status='replace',form='unformatted',
     2  access='direct',recl=iword2)
      close(intheta+6)

      open(intheta+6,
     1  file=fonpoin,status='replace',form='unformatted',
     2  access='direct',recl=iword2)
      write(intheta+6,rec=1) nsp
      close(intheta+6)

      open(intheta+6,
     1  file=foffpoi,status='replace',form='unformatted',
     2  access='direct',recl=iword2)
      close(intheta+6)

      open(intheta+6,
     1  file=foffpoin,status='replace',form='unformatted',
     2  access='direct',recl=iword2)
      write(intheta+6,rec=1) nsp
      close(intheta+6)

      open(intheta+6,
     1  file=fenoi,status='replace',form='unformatted',
     2  access='direct',recl=iword2)
      close(intheta+6)

      open(intheta+6,
     1  file=fenoin,status='replace',form='unformatted',
     2  access='direct',recl=iword2*nmax)
      write(intheta+6,rec=1) nsp
      close(intheta+6)

      open(intheta+6,
     1  file=finoi,status='replace',form='unformatted',
     2  access='direct',recl=iword2)
      close(intheta+6)

      open(intheta+6,
     1  file=finoin,status='replace',form='unformatted',
     2  access='direct',recl=iword2*nmax)
      write(intheta+6,rec=1) nsp
      close(intheta+6)

!  output histogram
      open(intheta+6,file=f1,status='replace',form='unformatted',
     1       access='direct',recl=iword2*iperiod*nmax+1)
      close(intheta+6)

      open(intheta+6,file=f2,status='replace',form='unformatted',
     1       access='direct',recl=iword*nmax)
      close(intheta+6)

!   print *,nmax, iword, iword*nmax
      open(intheta+6,file=f3,status='replace',form='unformatted',
     1       access='direct',recl=iword*nmax)
      close(intheta+6)

      nc = nint((tfinal-rstart)/period)
      open(intheta+6,file=f4,status='replace',form='unformatted',
     1       access='direct',recl=iword*nmax*iperiod)
      close(intheta+6)


      ncount = 0
      fn = 'samples.dat'
      fn = trim(epsfdr)//'/'//trim(thetafdr)//'/'//trim(fn)
  
      open(intheta*20000,
     1  file=fn,status='replace',form='unformatted',
     2  access='direct',recl=iword*9) 
      close(intheta*20000)


      open(intheta+6,file=f3,status='old',form='unformatted',
     1       access='direct',recl=iword*nmax)
      write(intheta+6,rec=1) (see(i),i=1,nmax)
      write(intheta+6,rec=2) (sei(i),i=1,nmax)
      write(intheta+6,rec=3) (sie(i),i=1,nmax)
      write(intheta+6,rec=4) (sii(i),i=1,nmax)
      write(intheta+6,rec=5) (exc(i),i=1,nmax)
      write(intheta+6,rec=6) (nlgni(i)*1.D0,i=1,nmax)
      close(intheta+6)

!---------------------------------------Transient noise & lgn 
      do i=1,jlgn
      pspon(i)  = -10000.D0
      pspoff(i) = -10000.D0
      glo(i) = 0.D0
      slo1(i) = 0.D0
      slo2(i) = 0.D0
      slo3(i) = 0.D0
      glon(i) = 0.D0
      slon1(i) = 0.D0
      slon2(i) = 0.D0
      slon3(i) = 0.D0
      glf(i) = 0.D0
      slf1(i) = 0.D0
      slf2(i) = 0.D0
      slf3(i) = 0.D0
      glfn(i) = 0.D0
      slfn1(i) = 0.D0
      slfn2(i) = 0.D0
      slfn3(i) = 0.D0
      enddo
      do i=1,nmax
      v(i) = 0.D0
      vnew(i) = 0.D0
      pspike(i) = -10000.D0
      alpha1(i) = 0.D0
      beta1(i) = 0.D0
      ge(i) = 0.D0
      se1(i) = 0.D0
      se2(i) = 0.D0
      se3(i) = 0.D0
      gf(i) = 0.D0
      sf1(i) = 0.D0
      sf2(i) = 0.D0
      sf3(i) = 0.D0
      gi(i) = 0.D0
      si1(i) = 0.D0
      si2(i) = 0.D0
      si3(i) = 0.D0
      gj(i) = 0.D0
      sj1(i) = 0.D0
      sj2(i) = 0.D0
      sj3(i) = 0.D0
          gx(i) = 0.D0
          sx1(i) = 0.D0
          sx2(i) = 0.D0
          sx3(i) = 0.D0
          gy(i) = 0.D0
          sy1(i) = 0.D0
          sy2(i) = 0.D0
          sy3(i) = 0.D0
          gn(i) = 0.D0
          sn1(i) = 0.D0
          sn2(i) = 0.D0
          sn3(i) = 0.D0
          gm(i) = 0.D0
          sm1(i) = 0.D0
          sm2(i) = 0.D0
          sm3(i) = 0.D0
          gp(i) = 0.D0
          sp1(i) = 0.D0
          sp2(i) = 0.D0
          sp3(i) = 0.D0
      enddo

      do i=1,ntrans
        t = -ntrans*dt+(i-1)*dt
        if (rk.EQ.2) then
            call eif_rk2(t,dt,period,iperiod,iid,idap)
        else
            if (rk.EQ.4) then
                call visual4(t,dt,iseed)
                call rk4(t,dt,period,iperiod,iid)
!!! ONLY for test with RK4 !!!!!!!!!!!!!!!!!!!
!              call rk2(t,dt,period,iperiod)
!!! ONLY for test with RK4 !!!!!!!!!!!!!!!!!!!
            else
                stop 'not implemented yet'
            endif
        endif
        call update(dt,iseed,rk,unitary,idap)
        do j=1,nmax
            v(j) = vnew(j)
        enddo
      enddo
      print *,iid,'transient finished'

      ncaj = 1
      do i=1,nmax
      do j = 1,iperiod
        glgn(i,j) = 0.D0
        gexc(i,j) = 0.D0
        ginh(i,j) = 0.D0
        gnexc(i,j) = 0.D0
        gninh(i,j) = 0.D0
        gtot(i,j) = 0.D0
        cond(i,j) = 0.D0
        vslave(i,j) = 0.D0
        vmem(i,j) = 0.D0
        gpota(i,j) = 0.D0
        glgn2(i,j) = 0.D0
        gexc2(i,j) = 0.D0
        ginh2(i,j) = 0.D0
        gnexc2(i,j) = 0.D0
        gninh2(i,j) = 0.D0
        gtot2(i,j) = 0.D0
        cond2(i,j) = 0.D0
        vslave2(i,j) = 0.D0
        vmem2(i,j) = 0.D0
        gpota2(i,j) = 0.D0
        irate(i,j) = 0
        clgn(i,j) = 0.D0
        cexc(i,j) = 0.D0
        cinh(i,j) = 0.D0
        cnexc(i,j) = 0.D0
        cninh(i,j) = 0.D0
        cpota(i,j) = 0.D0
        clgn2(i,j) = 0.D0
        cexc2(i,j) = 0.D0
        cinh2(i,j) = 0.D0
        cnexc2(i,j) = 0.D0
        cninh2(i,j) = 0.D0
        cpota2(i,j) = 0.D0
      enddo
        nca(i,ncaj) = 0
        currEt(i,ncaj) = 0.D0
        currEc(i,ncaj) = 0.D0
        currIc(i,ncaj) = 0.D0
        currEn(i,ncaj) = 0.D0
        currIn(i,ncaj) = 0.D0
        currp(i,ncaj) =0.D0
        currEt2(i,ncaj) =0.D0
        currEc2(i,ncaj) =0.D0
        currIc2(i,ncaj) =0.D0
        currEn2(i,ncaj) =0.D0
        currIn2(i,ncaj) =0.D0
        currp2(i,ncaj) =0.D0
        tcurrEt(i) =0.D0
        tcurrEc(i) =0.D0
        tcurrIc(i) =0.D0
        tcurrEn(i) =0.D0
        tcurrIn(i) =0.D0
        tcurrp(i) =0.D0
        tcurrEt2(i)=0.D0
        tcurrEc2(i)=0.D0
        tcurrIc2(i)=0.D0
        tcurrEn2(i)=0.D0
        tcurrIn2(i)=0.D0
        tcurrp2(i)=0.D0
        if(excite(i))then
           gexcitMax(i) = 0.D0
           ginhibMax(i) = 0.D0
           erate(i) = 0.D0
           gll(i) = 0.D0
        endif
      enddo
      t = 0.D0
      rtot = 0.D0
      rctot = 0.D0
      rstot = 0.D0
      rmtot = 0.D0
      ritot = 0.D0
      nsptot = 0
      isptot = 0
      jsptot = 0
      ksptot = 0
      msptot = 0
      condavg = 0.D0
      curravg = 0.D0
      condlgn = 0.D0
      condexc = 0.D0
      condinh = 0.D0
      geavg = 0.D0
      giavg = 0.D0
      gpavg = 0.D0
      vsuma = 0.D0
      vsumb = 0.D0
      vsuml = 0.D0
      vsume = 0.D0
      vsumi = 0.D0
      vsumen = 0.D0
      vsumin = 0.D0
      vsump = 0.D0
      vEMax = 0.D0
      erateMax = 0.D0
      gexcitMax0 = 0.D0
      vEMean = 0.D0
      vIMax = 0.D0
      vIMean = 0.D0
      millisec = nint(0.001D0/dt)
      !print *,iid,':','fine init',intheta
      
!-------------Start of Main Time Integration Loop for each dt
     
      open(5000+intheta,file=fonpoi,status='old',form='unformatted',
     1       access='direct',recl=iword2)
      open(5100+intheta,file=fonpoin,status='old',form='unformatted',
     1       access='direct',recl=iword2)

      open(6000+intheta,file=foffpoi,status='old',form='unformatted',
     1       access='direct',recl=iword2)
      open(6100+intheta,file=foffpoin,status='old',form='unformatted',
     1       access='direct',recl=iword2)

      open(7000+intheta,file=fenoi,status='old',form='unformatted',
     1       access='direct',recl=iword2)
      open(7100+intheta,file=fenoin,status='old',form='unformatted',
     1       access='direct',recl=iword2)

      open(8000+intheta,file=finoi,status='old',form='unformatted',
     1       access='direct',recl=iword2)
      open(8100+intheta,file=finoin,status='old',form='unformatted',
     1       access='direct',recl=iword2)
      open(9100+intheta,file=fo,status='old',form='unformatted',
     1       access='direct',recl=iword2*2)
      open(10000+intheta,file=fsp,status='old',form='unformatted',
     1       access='direct',recl=iword2)
      open(11000+intheta,file=fspn,status='old',form='unformatted',
     1       access='direct',recl=iword2)
        open(intheta*20000,
     1        file=fn,status='old',form='unformatted',
     2        access='direct',recl=iword*9)    

      do iii=1,ntotal
!------------------------------------------------------------
      if (rk.EQ.2) then
          call eif_rk2(t,dt,period,iperiod,iid,idap)
      else
          if (rk.EQ.4) then
              call visual4(t,dt,iseed)
              call rk4(t,dt,period,iperiod,iid)
!!! ONLY for test with RK4 !!!!!!!!!!!!!!!!!!!
!              call rk2(t,dt,period,iperiod)
!!! ONLY for test with RK4 !!!!!!!!!!!!!!!!!!!
          else
              stop 'not implemented yet'
          endif
      endif
!-------------------------Generate noise and LGN input spikes 
      call update(dt,iseed,rk,unitary,idap)
!---  -------------------------------------finally advance TIME
      vpop = 0.D0
      do i=1,nmax
         if (isptot.lt.0.or.nsptot.lt.0) then
             print *,'wtf too much spikes',nsptot,isptot
         endif
        v(i) = vnew(i)
        vpop = vpop + v(i)
        if (excite(i)) then
            ssse = see(i)
            sssi = sei(i)
            gggC = gEIc
            if (v(i).gt.vEMax) then
                vEMax = v(i)
            endif
            ginhib  = ((1.D0-fgaba)*gi(i) + fgaba*gj(i))*sssi + gEIc
            gexcit  = ((1.D0-fnmdacE)*ge(i) + fnmdacE*gf(i))*ssse
            vEMean = vEMean + v(i)    
            gexcitMax(i) = gexcitMax(i)+gexcit
            ginhibMax(i) = ginhibMax(i)+ginhib
            gll(i) = gll(i) + gl(i)
        else
          vIMean = vIMean + v(i)    
          if (v(i).gt.vIMax) then
              vIMax = v(i)
          endif
        endif
      enddo
!=============================================================
      if (nint(t/dt).ge.nint((tfinal-twindow)/dt)) then
        if(nspike.gt.0) then
      read(11000+intheta,rec=1) nsp
      do i=1,nspike
        write(10000+intheta,rec=2*(nsp+i)-1) sngl(ispike(i)*1.D0)
        write(10000+intheta,rec=2*(nsp+i)) sngl(t+tspike(i))
      enddo
      
      write(11000+intheta,rec=1) nsp + nspike
        endif
!      print *, vpop

      write(9100+intheta,rec=nint((t-tfinal+twindow)/dt)+1) 
     1  sngl(1.D0*nspike),sngl(vpop)
!-----------------------------------------------------------------
      if (recordAllspikes) then
      read(5100+intheta,rec=1) nsp
      do i=1,nonpoi
        write(5000+intheta,rec=(nsp+i)*2-1) sngl(ionpoi(i)*1.D0)
        write(5000+intheta,rec=(nsp+i)*2) sngl(tsonpoi(i))
      enddo

      write(5100+intheta,rec=1) nsp + nonpoi
!------------------------------------------------------------------
      read(6100+intheta,rec=1) nsp
      do i=1,noffpoi
        write(6000+intheta,rec=(nsp+i)*2-1) sngl(ioffpoi(i)*1.D0)
        write(6000+intheta,rec=(nsp+i)*2) sngl(tsoffpoi(i))
      enddo

      write(6100+intheta,rec=1) nsp + noffpoi
!------------------------------------------------------------------
      read(7100+intheta,rec=1) nsp
      do i=1,nenoi
        write(7000+intheta,rec=(nsp+i)*2-1) sngl(ienoi(i)*1.D0)
        write(7000+intheta,rec=(nsp+i)*2) sngl(tsenoi(i))
      enddo

      write(7100+intheta,rec=1) nsp + nenoi
!-----------------------------------------------------------------
      read(8100+intheta,rec=1) nsp
      do i=1,ninoi
        write(8000+intheta,rec=(nsp+i)*2-1) sngl(iinoi(i)*1.D0)
        write(8000+intheta,rec=(nsp+i)*2) sngl(tsinoi(i))
      enddo      

      write(8100+intheta,rec=1) nsp + ninoi
      endif
      endif
!=============================================================
      t = t + dt
!=============================================================
!------------------------------------Cycle Avg'd Conductances
!                           output file f4
      if((mod(iii,nstep0).eq.0).and.(nint(t/dt).gt.nint(rstart/dt)))then
!------------------------------------------------------------
      ncycle = mod(iii,nsteperiod)/nstep0
      if (ncycle.eq.0) then
        ncycle = iperiod
      endif
      do i=1,nmax
      if ( excite(i) ) then
        gei = ((1.D0-fnmdacE)*ge(i) + fnmdacE*gf(i)) * see(i)
        gii = ((1.D0-fgaba)*gi(i) + fgaba*gj(i)) * sei(i)
        gen = (1.D0-fnmdanE)*gx(i) + fnmdanE*gy(i)
        gin = (1.D0-fgaba)*gn(i) + fgaba*gm(i)
        gps = gp(i)*spotaE
      else
        gei = ((1.D0-fnmdacI)*ge(i) + fnmdacI*gf(i)) * sie(i)
        gii = ((1.D0-fgaba)*gi(i) + fgaba*gj(i)) * sii(i)
        gen = (1.D0-fnmdanI)*gx(i) + fnmdanI*gy(i)
        gin = (1.D0-fgaba)*gn(i) + fgaba*gm(i)
        gps = gp(i)*spotaI
      endif
        ve = v(i)-vexcit
        vi = v(i)-vinhib

!    vs: effective potential
        vsi = -beta1(i)/alpha1(i)
        glgn(i,ncycle) = gl(i) + glgn(i,ncycle)
        gexc(i,ncycle) = gei + gexc(i,ncycle)
        ginh(i,ncycle) = gii + ginh(i,ncycle)
        gtot(i,ncycle) = -alpha1(i) + gtot(i,ncycle)
        cond(i,ncycle) = beta1(i) + cond(i,ncycle)
        vslave(i,ncycle) = vsi + vslave(i,ncycle)
        gnexc(i,ncycle) = gen + gnexc(i,ncycle)
        gninh(i,ncycle) = gin + gninh(i,ncycle)
        gpota(i,ncycle) = gps + gpota(i,ncycle)

        glgn2(i,ncycle) = gl(i)*gl(i) + glgn2(i,ncycle)
        gexc2(i,ncycle) = gei*gei + gexc2(i,ncycle)
        ginh2(i,ncycle) = gii*gii + ginh2(i,ncycle)
        gtot2(i,ncycle) = alpha1(i)*alpha1(i) + gtot2(i,ncycle)
        cond2(i,ncycle) = beta1(i)*beta1(i) + cond2(i,ncycle)
        vslave2(i,ncycle) = vsi*vsi + vslave2(i,ncycle)
        gnexc2(i,ncycle) = gen*gen + gnexc2(i,ncycle)
        gninh2(i,ncycle) = gin*gin + gninh2(i,ncycle)
        gpota2(i,ncycle) = gps*gps + gpota2(i,ncycle)

        vmem(i,ncycle)  = v(i)      + vmem(i,ncycle)
        vmem2(i,ncycle) = v(i)*v(i) + vmem2(i,ncycle)

        clgn(i,ncycle) =  clgn(i,ncycle) - gl(i)*ve
        cexc(i,ncycle) =  cexc(i,ncycle) - gei*ve  
        cinh(i,ncycle) =  cinh(i,ncycle) - gii*vi  
        cnexc(i,ncycle) = cnexc(i,ncycle) - gen*ve 
        cninh(i,ncycle) = cninh(i,ncycle) - gin*vi 
        cpota(i,ncycle) = cpota(i,ncycle) - gps*vi 

        clgn2(i,ncycle) = gl(i)*ve*gl(i)*ve + clgn2(i,ncycle)
        cexc2(i,ncycle) = gei*ve*gei*ve + cexc2(i,ncycle)
        cinh2(i,ncycle) = gii*vi*gii*vi + cinh2(i,ncycle)
        cnexc2(i,ncycle) = gen*ve*gen*ve + cnexc2(i,ncycle)
        cninh2(i,ncycle) = gin*vi*gin*vi + cninh2(i,ncycle)
        cpota2(i,ncycle) = gps*vi*gps*vi + cpota2(i,ncycle)
      enddo
!------------------------------------------------------------
      endif
!------------------------------------------------------------
      if (nint(t/dt).gt.nint(rstart/dt).or.
     1      nint(t/dt).gt.nint((tfinal-twindow)/dt)) then
!-----------------sample neuron samples' stats
!--------------and current for every neuron per millisecond 
      if (mod(iii,millisec).eq.0) then
        j = 1
        do i=1,nmax
          if ( excite(i) ) then
            tref = trefE
            gei = ((1.D0-fnmdacE)*ge(i) + fnmdacE*gf(i)) * see(i)
            gii = ((1.D0-fgaba)*gi(i) + fgaba*gj(i)) * sei(i)
            gen = (1.D0-fnmdanE)*gx(i) + fnmdanE*gy(i)
            gin = (1.D0-fgaba)*gn(i) + fgaba*gm(i)
            gps = gp(i) * spotaE
          else
            tref = trefI
            gei = ((1.D0-fnmdacI)*ge(i) + fnmdacI*gf(i)) * sie(i)
            gii = ((1.D0-fgaba)*gi(i) + fgaba*gj(i)) * sii(i)
            gen = (1.D0-fnmdanI)*gx(i) + fnmdanI*gy(i)
            gin = (1.D0-fgaba)*gn(i) + fgaba*gm(i)
            gps = gp(i) * spotaI
          endif
          if (nint(t/dt).gt.nint(rstart/dt)) then
          if (pspike(i) + tref.lt. t) then
            ve = v(i)-vexcit
            vi = v(i)-vinhib
            currEc(i,ncaj) = currEc(i,ncaj) - gei*ve
            currIc(i,ncaj) = currIc(i,ncaj) - gii*vi
            currEt(i,ncaj) = currEt(i,ncaj) - gl(i)*ve
            currEn(i,ncaj) = currEn(i,ncaj) - gen*ve
            currIn(i,ncaj) = currIn(i,ncaj) - gin*vi
            currp(i,ncaj) = currp(i,ncaj) - gps*vi
            currEc2(i,ncaj) = currEc2(i,ncaj) + (gei*ve)**2
            currIc2(i,ncaj) = currIc2(i,ncaj) + (gii*vi)**2
            currEt2(i,ncaj) = currEt2(i,ncaj) + (gl(i)*ve)**2
            currEn2(i,ncaj) = currEn2(i,ncaj) + (gen*ve)**2
            currIn2(i,ncaj) = currIn2(i,ncaj) + (gin*vi)**2
            currp2(i,ncaj) = currp2(i,ncaj) + (gps*vi)**2
            nca(i,ncaj) = nca(i,ncaj) + 1
          endif
          endif
          if ((nint(t/dt).ge.nint((tfinal-twindow)/dt)).and.
     1       (i.eq.in(j))) then
            write(intheta*20000,rec=ncount*nsample+j) 
     1       v(i),gl(i),gei,gii,gen,gin,
     2       -alpha1(i),beta1(i),gps
            if (j.eq.nsample) then
              ncount = ncount + 1
            endif
            j = j + 1
          endif
        enddo 
      endif
!------------------------------------------------------------
      endif
!------------------------------------------------------------
      if (mod(iii,nsteperiod).eq.0) then
!------------------------------------------------------------
!$OMP CRITICAL
      label(1) = 'Time:'
      label(2) = 'Nspike:'
      label(3) = 'eRate:'
      label(4) = 'iRate:'
      label(5) = 'cRate:'
      label(6) = 'sRate:'
      label(7) = 'mRate:'
      print 1000, iid,label(1),t,label(2),nsptot,
     1        label(3),sngl(nsptot/period/nv1e),
     2        label(4),sngl(isptot/period/nv1i),
     3        label(5),sngl(jsptot/period/nv1c),
     4        label(6),sngl(ksptot/period/nv1s),
     5        label(7),sngl(msptot/period/nv1m)
1000  format (I12,A7,F7.3,A7,I7,A7,F7.2,A7,F7.2,
     1 A7,F7.2,A7,F7.2,A7,F7.2)
      print *,iid
      label(1) = 'E gtot:'
      label(2) = 'Itot:'
      label(3) = 'vMean:'
      label(4) = 'gLGN:'
      label(5) = 'gExc:'
      label(6) = 'gInh:'
      label(7) = 'gEn:'
      label(8) = 'gIn:'
      label(9) = 'vMax:'
      label(10) = 'gp:'
      print 1001,iid,label(1),sngl(condavg/steperiod/nv1e),
     1        label(2),sngl(curravg/steperiod/nv1e),
     2        label(3),sngl(vEMean/steperiod/nv1e),
     3        label(4),sngl(condlgn/steperiod/nv1e),
     4        label(5),sngl(geavg/steperiod/nv1e),
     5        label(6),sngl(giavg/steperiod/nv1e),
     6        label(7),sngl(condexc/steperiod/nv1e),
     7        label(8),sngl(condinh/steperiod/nv1e),
     8        label(9),sngl(vEMax),
     9        label(10),sngl(gpavg/steperiod/nv1e)
      do i=1,nex*ney
          if(erate(i).gt.erateMax)then
            irID = i
            glrate = gll(i)
            nlgnrateMax = nlgni(i)
            erateMax = erate(i)
            gexcitrMax = gexcitMax(i)
            nPreERMax = NpreSynE(i)
            ginhibrMax = ginhibMax(i)
            nPreIRMax = NpreSynI(i)
          endif
          if(gexcitMax(i).gt.gexcitMax0)then
            igID = i
            glgE = gll(i)
            nlgngEmax = nlgni(i)
            erateEGMax = erate(i)
            gexcitMax0 = gexcitMax(i)
            nPreEGMax = NpreSynE(i)
            ginhibMax0 = ginhibMax(i)
            nPreIGMax = NpreSynI(i)
          endif
          gexcitMax(i) = 0.D0
          ginhibMax(i) = 0.D0
          gll(i) = 0.D0
      enddo
      npdf(7)=0
      do i=0,5
            npdf(i+1)=0
            rpdf(i+1) = erateMax/period*dble(i)/6.D0
      enddo
      do i=1,nex*ney
          do j=1,6
            if (erate(i)/period.gt.rpdf(j))then 
                if (j.EQ.6)then
                    npdf(7) = npdf(7) + 1
                else
                    cycle 
                endif
            else
                npdf(j) = npdf(j) + 1
                exit 
            endif
          enddo
          erate(i) = 0.D0
      enddo
1001  format (I12,A7,F7.2,A7,F7.2,A7,F7.2,A7,F7.2,
     1 A7,F7.2,A7,F7.2,A7,F6.2,A7,F6.2,A7,F3.1,A7,F6.2)
      label(1) = 'I gtot:'
      print 1001,iid,label(1),sngl(vsuma/steperiod/nv1i),
     1        label(2),sngl(vsumb/steperiod/nv1i),
     2        label(3),sngl(vIMean/steperiod/nv1i),
     3        label(4),sngl(vsuml/steperiod/nv1i),
     4        label(5),sngl(vsume/steperiod/nv1i),
     5        label(6),sngl(vsumi/steperiod/nv1i),
     6        label(7),sngl(vsumen/steperiod/nv1i),
     7        label(8),sngl(vsumin/steperiod/nv1i),
     8        label(9),sngl(vIMax),
     9        label(10),sngl(vsump/steperiod/nv1i)
      print *,iid
      print 1003,irID,' Rate: ',sngl(erateMax/period),
     1        ' gLGN: ',sngl(glrate/steperiod),
     2        ' nlgn: ',nlgnrateMax,
     3          ' gE: ',sngl(gexcitrMax/steperiod),
     4       ' nPreE: ',nPreERMax,
     5          ' gI: ',sngl(ginhibrMax/steperiod),
     6       ' nPreI: ',nPreIRMax
      print 1003,igID,' Rate: ',sngl(erateEGMax/period),
     1        ' gLGN: ',sngl(glgE/steperiod),
     2        ' nlgn: ',nlgngEmax,
     3          ' gE: ',sngl(gexcitMax0/steperiod),
     4       ' nPreE: ',nPreEGMax,
     5          ' gI: ',sngl(ginhibMax0/steperiod),
     6       ' nPreI: ',nPreIGMax
1003  format (I12,A8,F7.2,A8,F7.2,A8,I7,A8,F7.2,A8,I7,
     1 A8,F7.2,A8,I7)
1004  format (A5,I4,A5,I4,A5,I4,A5,I4,A5,I4,A5,I4,A5,I4,A5)
1005  format (A8,F5.1,A4,F5.1,A4,F5.1,A4,F5.1,A4,F5.1,A4,F5.1,
     1 A4,F5.1)
      print *,fence1
      print 1004,'  |  ',npdf(1),'  |  ',npdf(2),'  |  ',npdf(3),
     1  '  |  ',npdf(4),'  |  ',npdf(5),'  |  ',npdf(6),'  |  ',
     2  npdf(7),'  |  '
      print 1005,'          ',rpdf(1),'    ',rpdf(2),'    ',
     1  rpdf(3),'    ',rpdf(4),'    ',rpdf(5),'    ',rpdf(6),'    ',
     2  sngl(erateMax/period)
      print *,fence1
1002  format (I12,A25,F5.1,A1)
      print 1002,ineps,fence0,sngl(t*100.D0/tfinal),'%'
!$OMP END CRITICAL
      if (t.gt.rstart) rtot = rtot+nsptot/period/nv1e
      if (t.gt.rstart) rctot = rctot+jsptot/period/nv1c
      if (t.gt.rstart) rstot = rstot+ksptot/period/nv1s
      if (t.gt.rstart) rmtot = rmtot+msptot/period/nv1m
      if (t.gt.rstart) ritot = ritot+isptot/period/nv1i
      vEMax = 0.D0
      vEMean = 0.D0
      vIMax = 0.D0
      gexcitMax0 = 0.D0
      erateMax = 0.D0
      vIMean = 0.D0
      isptot = 0
      nsptot = 0
      jsptot = 0
      ksptot = 0
      msptot = 0
      condavg = 0.D0
      curravg = 0.D0
      condlgn = 0.D0
      condexc = 0.D0
      condinh = 0.D0
      geavg = 0.D0
      giavg = 0.D0
      gpavg = 0.D0
      vsuma = 0.D0
      vsumb = 0.D0
      vsuml = 0.D0
      vsume = 0.D0
      vsumi = 0.D0
      vsumen = 0.D0
      vsumin = 0.D0
      vsump = 0.D0
      if (nint(t/dt).gt.nint(rstart/dt)) then
        ncaj = ncaj + 1
        if (ncaj.gt.max_ncaj) then
            print *,'increase record length of current array'
            print *,'ncaj = ', ncaj
        endif
        do i = 1,nmax
           nca(i,ncaj) = 0
        currEc(i,ncaj) = 0.D0
        currIc(i,ncaj) = 0.D0
        currEt(i,ncaj) = 0.D0
        currEn(i,ncaj) = 0.D0
        currIn(i,ncaj) = 0.D0
        currp(i,ncaj) = 0.D0
        currEc2(i,ncaj) =0.D0
        currIc2(i,ncaj) =0.D0
        currEt2(i,ncaj) =0.D0
        currEn2(i,ncaj) =0.D0
        currIn2(i,ncaj) =0.D0
        currp2(i,ncaj) =0.D0
        enddo
      endif
!------------------------------------------------------------
      endif
!------------------------------------------------------------
      enddo
      close(5000+intheta)
      close(5100+intheta)
      close(6000+intheta)
      close(6100+intheta)
      close(7000+intheta)
      close(7100+intheta)
      close(8000+intheta)
      close(8100+intheta)
      close(9100+intheta)
      close(10000+intheta)
      close(11000+intheta)
      close(intheta*20000)
!--------------------------------------End of TimeStepping Loop
!------------------------------------------------Output to f4
      do j=1,iperiod
        do i=1,nmax
            glgn(i,j) = glgn(i,j)/nc
            gexc(i,j) = gexc(i,j)/nc
            ginh(i,j) = ginh(i,j)/nc
            gtot(i,j) = gtot(i,j)/nc
            cond(i,j) = cond(i,j)/nc
            vslave(i,j) = vslave(i,j)/nc
            gnexc(i,j) = gnexc(i,j)/nc
            gninh(i,j) = gninh(i,j)/nc
            vmem(i,j)  = vmem(i,j)/nc
            gpota(i,j)  = gpota(i,j)/nc

            glgn2(i,j) = glgn2(i,j)/nc
            gexc2(i,j) = gexc2(i,j)/nc
            ginh2(i,j) = ginh2(i,j)/nc
            gtot2(i,j) = gtot2(i,j)/nc
            cond2(i,j) = cond2(i,j)/nc
            vslave2(i,j) = vslave2(i,j)/nc
            gnexc2(i,j) = gnexc2(i,j)/nc
            gninh2(i,j) = gninh2(i,j)/nc
            vmem2(i,j) = vmem2(i,j)/nc
            gpota2(i,j)  = gpota2(i,j)/nc

            clgn(i,j) = clgn(i,j)/nc
            cexc(i,j) = cexc(i,j)/nc
            cinh(i,j) = cinh(i,j)/nc
            cnexc(i,j) = cnexc(i,j)/nc
            cninh(i,j) = cninh(i,j)/nc
            cpota(i,j)  = cpota(i,j)/nc

            clgn2(i,j) = clgn2(i,j)/nc
            cexc2(i,j) = cexc2(i,j)/nc
            cinh2(i,j) = cinh2(i,j)/nc
            cnexc2(i,j) = cnexc2(i,j)/nc
            cninh2(i,j) = cninh2(i,j)/nc
            cpota2(i,j)  = cpota2(i,j)/nc
          enddo
      enddo
      open(intheta+6,file=f4,status='old',form='unformatted',
     1        access='direct',recl=iword*nmax*iperiod)
        write(intheta+6,rec=1) glgn
          write(intheta+6,rec=2) gexc
          write(intheta+6,rec=3) ginh
          write(intheta+6,rec=4) gtot
          write(intheta+6,rec=5) cond
          write(intheta+6,rec=6) vslave
          write(intheta+6,rec=7) vmem
          write(intheta+6,rec=8) gnexc
          write(intheta+6,rec=9) gninh
          write(intheta+6,rec=10) glgn2
          write(intheta+6,rec=11) gexc2
          write(intheta+6,rec=12) ginh2
          write(intheta+6,rec=13) gtot2
          write(intheta+6,rec=14) cond2
          write(intheta+6,rec=15) vslave2
          write(intheta+6,rec=16) vmem2
          write(intheta+6,rec=17) gnexc2
          write(intheta+6,rec=18) gninh2
          write(intheta+6,rec=19) gpota
          write(intheta+6,rec=20) gpota2
          write(intheta+6,rec=21) clgn
          write(intheta+6,rec=22) cexc
          write(intheta+6,rec=23) cinh 
          write(intheta+6,rec=24) cnexc
          write(intheta+6,rec=25) cninh
          write(intheta+6,rec=26) cpota 
          write(intheta+6,rec=27) clgn2
          write(intheta+6,rec=28) cexc2
          write(intheta+6,rec=29) cinh2
          write(intheta+6,rec=30) cnexc2
          write(intheta+6,rec=31) cninh2
          write(intheta+6,rec=32) cpota2

      close(intheta+6)
!------------------------------------------------Output to f1
      open(intheta+6,file=f1,status='old',form='unformatted',
     1        access='direct',recl=iword2*iperiod*nmax)
      write(intheta+6,rec=1) irate
      close(intheta+6)
      tist1 = tfinal-rstart
      print *, tist1
      open(intheta+6,file=f1,status='old',form='unformatted',
     1        access='direct',recl=iword2)
      write(intheta+6,rec=iperiod*nmax+1) sngl(tist1)
      close(intheta+6)
!----------------------------------------------Output to f4
      do i=1,nmax 
          sumnca(i) = 0.D0
          do j=1,ncaj-1
              sumnca(i) = sumnca(i)+nca(i,j)
          enddo
      enddo
      !sumnca = sum(nca,2)
      do i = 1,nmax
        do j = 1,ncaj-1
          tcurrEt(i) = tcurrEt(i) + currEt(i,j)/(sumnca(i)*1.D0)
          tcurrEc(i) = tcurrEc(i) + currEc(i,j)/(sumnca(i)*1.D0)
          tcurrIc(i) = tcurrIc(i) + currIc(i,j)/(sumnca(i)*1.D0)
          tcurrEn(i) = tcurrEn(i) + currEn(i,j)/(sumnca(i)*1.D0)
          tcurrIn(i) = tcurrIn(i) + currIn(i,j)/(sumnca(i)*1.D0)
          tcurrp(i) = tcurrp(i) + currp(i,j)/(sumnca(i)*1.D0)
          tcurrEt2(i)= tcurrEt2(i) + currEt2(i,j)/(sumnca(i)*1.D0)
          tcurrEc2(i)= tcurrEc2(i) + currEc2(i,j)/(sumnca(i)*1.D0)
          tcurrIc2(i)= tcurrIc2(i) + currIc2(i,j)/(sumnca(i)*1.D0)
          tcurrEn2(i)= tcurrEn2(i) + currEn2(i,j)/(sumnca(i)*1.D0)
          tcurrIn2(i)= tcurrIn2(i) + currIn2(i,j)/(sumnca(i)*1.D0)
          tcurrp2(i)= tcurrp2(i) + currp2(i,j)/(sumnca(i)*1.D0)
        enddo
      enddo
      open(intheta+6,file=f2,status='old',form='unformatted',
     1       access='direct',recl=iword*nmax)
        write(intheta+6,rec=1) tcurrEt
        write(intheta+6,rec=2) tcurrEc
        write(intheta+6,rec=3) tcurrIc
        write(intheta+6,rec=4) tcurrEn
        write(intheta+6,rec=5) tcurrIn
        write(intheta+6,rec=6) tcurrEt2
        write(intheta+6,rec=7) tcurrEc2
        write(intheta+6,rec=8) tcurrIc2
        write(intheta+6,rec=9) tcurrEn2
        write(intheta+6,rec=10) tcurrIn2
        write(intheta+6,rec=11) tcurrp
        write(intheta+6,rec=12) tcurrp2
      close(intheta+6)

!---------------------------------------------------------
!$OMP ORDERED
      print *,'-------------',sngl(gtheta/twopi*360),'-------------'
      print *,' Exc Rate = ',sngl(rtot/nc)
      print *,' Inh Rate = ',sngl(ritot/nc)
      print *,'    Pure Complex Rate = ',sngl(rctot/nc) 
      print *,'    Pure Simple Rate = ',sngl(rstot/nc)
      print *,'    Middle Mixed Rate = ',sngl(rmtot/nc)
      print *,' # period = ', ncaj-1,' == ', nc
!$OMP END ORDERED
!------------------------------------------------------------
        enddo
!$OMP END PARALLEL DO 
!!------------------------------------------------------------
      print *,'-- Drifting grating of contrast = ',sngl(epss),'ended'
      print *,' ' 
        enddo
!!============================================================
9999  stop
      end
!************************************************************
       !subroutine update(dt,t,iseed)
       subroutine update(dt,iseed,rk,unitary,idap)
       use parameters
!------------------------------------------------------------
!
!  Update chain for excitatory & inhibitory after spikes
!
!------------------------------------------------------------
!* Questions:
!    necessity of disort spikes fisrt
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension ge(nmax),se1(nmax),se2(nmax),se3(nmax)
      dimension gf(nmax),sf1(nmax),sf2(nmax),sf3(nmax)
      dimension gi(nmax),si1(nmax),si2(nmax),si3(nmax)
      dimension gj(nmax),sj1(nmax),sj2(nmax),sj3(nmax)
      dimension gx(nmax),sx1(nmax),sx2(nmax),sx3(nmax)
      dimension gy(nmax),sy1(nmax),sy2(nmax),sy3(nmax)
      dimension gn(nmax),sn1(nmax),sn2(nmax),sn3(nmax)
      dimension gm(nmax),sm1(nmax),sm2(nmax),sm3(nmax)
      dimension gp(nmax),sp1(nmax),sp2(nmax),sp3(nmax)
      dimension tspike(nmax),ispike(nmax),nlgni(nmax)
      dimension gl(nmax),see(nmax),sei(nmax),sie(nmax),sii(nmax)
      dimension aa_ee(nmax),aa_ei(nmax),aa_ie(nmax),aa_ii(nmax)
      dimension beta1(nmax),alpha1(nmax)
      logical excite(nmax)
      dimension IpostSyn(nmax,nmax) 
      dimension NpostSynE(nmax), NpostSynI(nmax)
      dimension NpreSynE(nmax), NpreSynI(nmax)
      dimension profile(256),Istrength(nmax,nmax)
      dimension slist(256,nmax)
      real*8 unitary
      integer rk,iseed
      data iword  / 8 /
      !integer iid, OMP_GET_THREAD_NUM
      common / lgncnd / gl
      common / vconst / vthres,vreset,vexcit,vinhib,gleak,gleakI
      common / smatrx / see,sei,sie,sii,spotaE,spotaI
      common / chaine / ge,se1,se2,se3
      common / chaine2/ gf,sf1,sf2,sf3
      common / chaini / gi,si1,si2,si3
      common / chainj / gj,sj1,sj2,sj3
      common / chainx / gx,sx1,sx2,sx3
      common / chainy / gy,sy1,sy2,sy3
      common / chainn / gn,sn1,sn2,sn3
      common / chainm / gm,sm1,sm2,sm3
      common / chainp / gp,sp1,sp2,sp3
      common /   rhs  / alpha1,beta1
      common /  avgs  / suma,sumb,suml,sumen,sumin,sume,sumi,sump
      common /  avgsi / vsuma,vsumb,vsuml,vsumen,vsumin,vsume,vsumi,
     1 vsump
      common / tconst / tau_e,tau_i,tau0,tau1,tau2,tnrise,tndamp,
     1 tau_a0, tau_a1
      common /  NMDA  / fnmdatE,fnmdacE,fnmdatI,fnmdacI,
     1   fnmdanE, fnmdanI, fgaba
      common / neuron / excite,nlgni,trefE,trefI,nr,
     1   lgnmax,lgnmin,meanlgn
      common / ttotal / tfinal,twindow,rstart
      common / connec / synfailEE,synfailEI,synfailIE,synfailII,
     1   neffEE,neffEI,neffIE,neffII
      common / jumpsize / aa_ee,aa_ei,aa_ie,aa_ii,profile,
     1   slist
      common / post / NpostSynE, NpostSynI, IpostSyn,
     1   NpreSynE, NpreSynI,Istrength
      common / spikes / tspike,ispike,nspike
      common / ggIc / gIIc, gEIc
!$OMP THREADPRIVATE(/spikes/,/chainp/,
!$OMP& /chaine/,/chaine2/,/lgncnd/,
!$OMP& /chaini/,/chainj/,/chainx/,/chainy/,/chainn/,
!$OMP& /chainm/,/rhs/,/avgs/,/avgsi/
!$OMP& )
      iword2 = iword/2
!=============================================================
      if (nspike.gt.0) then
!=============================================================
      do j=1,nspike
      ij = ispike(j)
      dtt = dt - tspike(j)
      if (idap.eq.1) then
          tgp = 0.D0
          tsp1 = 0.D0
          tsp2 = 0.D0
          tsp3 = 1.D0/tau_a0
          if (rk.EQ.2) then
            call decay1_d(tgp,tsp3,tau_a0,tau_a1,dtt)
          else
            if (rk.EQ.4) then
          call decay3_d(tgp,tsp1,tsp2,tsp3,tau_a0,tau_a1,dtt)
            endif
          endif
          s = 1.D0
          gp(ij)  = gp(ij)  + tgp  *s
          sp1(ij) = sp1(ij) + tsp1 *s
          sp2(ij) = sp2(ij) + tsp2 *s
          sp3(ij) = sp3(ij) + tsp3 *s
      endif
!------------------------------------------------------------
        if (synfailEE.GE.0.D0 .and. synfailEE.LT.1.D0.or.
     1     synfailIE.GE.0.D0 .and. synfailIE.LT.1.D0.or.
     2     synfailEI.GE.0.D0 .and. synfailEI.LT.1.D0.or.
     3     synfailII.GE.0.D0 .and. synfailII.LT.1.D0) then
!----------------------------------------------
      !reset temporary condunctance
      if (excite(ij)) then
        tge_e = 0.D0
        tse1_e = 0.D0
        tse2_e = 0.D0
        tse3_e = 1.D0

        tgf_e = 0.D0
        tsf1_e = 0.D0
        tsf2_e = 0.D0
        tsf3_e = 1.D0/tnrise*tau_e
        if (rk.EQ.2) then
          call decay1_d(tge_e,tse3_e,tau_e,tau0,dtt)
          if (fnmdacE.GT.0.D0.or.fnmdacI.GT.0.D0) then
              call decay1_d(tgf_e,tsf3_e,tnrise,tndamp,dtt)
          endif
        else
          if (rk.EQ.4) then
            call decay3_d(tge_e,tse1_e,tse2_e,tse3_e,tau_e,tau0,dtt)
            if (fnmdacE.GT.0.D0.or.fnmdacI.GT.0.D0) then
              call decay3_d(tgf_e,tsf1_e,tsf2_e,
     1                      tsf3_e,tnrise,tndamp,dtt)
            endif
          endif
        endif
        do ijj=1,NpostSynE(ij)+NpostSynI(ij)
            i = IpostSyn(ij,ijj)
            if (unitary .eq. 0.5D0) then
                s = profile(Istrength(ij,ijj))
            else
                if (dabs(unitary-0.6D0).lt.1e-14)then
                    s = slist(Istrength(ij,ijj),i)
                else
                    s = 1.D0;
                endif
            endif
          if (excite(i)) then
            if (ran2(iseed).LT.synfailEE) then
        ge(i)  = ge(i)  + tge_e  * aa_ee(i) *s
        gf(i)  = gf(i)  + tgf_e  * aa_ee(i) *s
        se1(i) = se1(i) + tse1_e * aa_ee(i) *s
        sf1(i) = sf1(i) + tsf1_e * aa_ee(i) *s
        se2(i) = se2(i) + tse2_e * aa_ee(i) *s
        sf2(i) = sf2(i) + tsf2_e * aa_ee(i) *s
        se3(i) = se3(i) + tse3_e * aa_ee(i) *s
        sf3(i) = sf3(i) + tsf3_e * aa_ee(i) *s
            endif
          else
            if (ran2(iseed).LT.synfailIE) then
        ge(i)  = ge(i)  + tge_e  * aa_ie(i) *s
        gf(i)  = gf(i)  + tgf_e  * aa_ie(i) *s
        se1(i) = se1(i) + tse1_e * aa_ie(i) *s
        sf1(i) = sf1(i) + tsf1_e * aa_ie(i) *s
        se2(i) = se2(i) + tse2_e * aa_ie(i) *s
        sf2(i) = sf2(i) + tsf2_e * aa_ie(i) *s
        se3(i) = se3(i) + tse3_e * aa_ie(i) *s
        sf3(i) = sf3(i) + tsf3_e * aa_ie(i) *s
            endif
          endif
        enddo
!----------------------------------------------------------
      else
!----------------------------------------------------------
        tgi_i = 0.D0
        tsi1_i = 0.D0
        tsi2_i = 0.D0
        tsi3_i = 1.D0 

        tgj_i = 0.D0
        tsj1_i = 0.D0
        tsj2_i = 0.D0
        tsj3_i = 1.D0/tau2*tau_i
        if (rk.EQ.2) then
          call decay1_d(tgi_i,tsi3_i,tau_i,tau1,dtt)
          if (fgaba.GT.0.D0) then
            call decay1(tgj_i,tsj3_i,tau2,dtt)
          endif
        else
          if (rk.EQ.4) then
            call decay3_d(tgi_i,tsi1_i,tsi2_i,tsi3_i,tau_i,tau1,dtt)
            if (fgaba.GT.0.D0) then
              call decay3(tgj_i,tsj1_i,tsj2_i,tsj3_i,tau2,dtt)
            endif
          endif
        endif
        do ijj=1,NpostSynE(ij)+NpostSynI(ij)
            i = IpostSyn(ij,ijj)
            if (unitary .eq. 0.5D0) then
                s = profile(Istrength(ij,ijj))
            else
                if (dabs(unitary-0.6D0).lt.1e-14)then
                    s = slist(Istrength(ij,ijj),i)
                else
                    s = 1.D0;
                endif
            endif
          if (excite(i)) then
            if (ran2(iseed).LT.synfailEI) then
        gi(i)  = gi(i)  + tgi_i  * aa_ei(i) *s
        gj(i)  = gj(i)  + tgj_i  * aa_ei(i) *s
        si1(i) = si1(i) + tsi1_i * aa_ei(i) *s
        sj1(i) = sj1(i) + tsj1_i * aa_ei(i) *s
        si2(i) = si2(i) + tsi2_i * aa_ei(i) *s
        sj2(i) = sj2(i) + tsj2_i * aa_ei(i) *s
        si3(i) = si3(i) + tsi3_i * aa_ei(i) *s
        sj3(i) = sj3(i) + tsj3_i * aa_ei(i) *s
            endif
          else
            if (ran2(iseed).LT.synfailII) then
        gi(i)  = gi(i)  + tgi_i  * aa_ii(i) *s
        gj(i)  = gj(i)  + tgj_i  * aa_ii(i) *s
        si1(i) = si1(i) + tsi1_i * aa_ii(i) *s
        sj1(i) = sj1(i) + tsj1_i * aa_ii(i) *s
        si2(i) = si2(i) + tsi2_i * aa_ii(i) *s
        sj2(i) = sj2(i) + tsj2_i * aa_ii(i) *s
        si3(i) = si3(i) + tsi3_i * aa_ii(i) *s
        sj3(i) = sj3(i) + tsj3_i * aa_ii(i) *s
            endif
          endif
        enddo
      endif
!------------ 1 > synfail > 0 ------------------
        else 
        if (synfailEE.eq.1.0.and.synfailIE.eq.1.0
     1      .and.synfailII.eq.1.0.and.synfailEI.eq.1.0) then
!-------------- synfail == 1 -------------------
      if (excite(ij)) then
        tge_e = 0.D0
        tse1_e = 0.D0
        tse2_e = 0.D0
        tse3_e = 1.D0

        tgf_e = 0.D0
        tsf1_e = 0.D0
        tsf2_e = 0.D0
        tsf3_e = 1.D0/tnrise*tau_e
        if (rk.EQ.2) then
        call decay1_d(tge_e,tse3_e,tau_e,tau0,dtt)
        if (fnmdacE.GT.0.D0.or.fnmdacI.gt.0.D0) then
            call decay1_d(tgf_e,tsf3_e,tau_e,tau0,dtt)
        endif
        else
          if (rk.EQ.4) then
            call decay3_d(tge_e,tse1_e,tse2_e,tse3_e,tau_e,tau0,dtt)
            if (fnmdacE.GT.0.D0.or.fnmdacI.gt.0.D0) then
                call decay3_d(tgf_e,tsf1_e,tsf2_e,tsf3_e,tau_e,tau0,dtt)
            endif
          endif
        endif
        do ijj=1,NpostSynE(ij)+NpostSynI(ij)
            i = IpostSyn(ij,ijj)
            if (unitary .eq. 0.5D0) then
                s = profile(Istrength(ij,ijj))
            else
                if (dabs(unitary-0.6D0).lt.1e-14)then
                    s = slist(Istrength(ij,ijj),i)
                else
                    s = 1.D0;
                endif
            endif
          if (excite(i)) then
        ge(i)  = ge(i)  + tge_e * aa_ee(i)  *s
        gf(i)  = gf(i)  + tgf_e * aa_ee(i)  *s
        se1(i) = se1(i) + tse1_e * aa_ee(i) *s
        sf1(i) = sf1(i) + tsf1_e * aa_ee(i) *s
        se2(i) = se2(i) + tse2_e * aa_ee(i) *s
        sf2(i) = sf2(i) + tsf2_e * aa_ee(i) *s
        se3(i) = se3(i) + tse3_e * aa_ee(i) *s
        sf3(i) = sf3(i) + tsf3_e * aa_ee(i) *s
          else
        ge(i)  = ge(i)  + tge_e * aa_ie(i)  *s
        gf(i)  = gf(i)  + tgf_e * aa_ie(i)  *s
        se1(i) = se1(i) + tse1_e * aa_ie(i) *s
        sf1(i) = sf1(i) + tsf1_e * aa_ie(i) *s
        se2(i) = se2(i) + tse2_e * aa_ie(i) *s
        sf2(i) = sf2(i) + tsf2_e * aa_ie(i) *s
        se3(i) = se3(i) + tse3_e * aa_ie(i) *s
        sf3(i) = sf3(i) + tsf3_e * aa_ie(i) *s
          endif
        enddo
      else
        tgi_i = 0.D0
        tsi1_i = 0.D0
        tsi2_i = 0.D0
        tsi3_i = 1.D0

        tgj_i = 0.D0
        tsj1_i = 0.D0
        tsj2_i = 0.D0
        tsj3_i = 1.D0/tau2*tau_i
        if (rk.EQ.2)then
          call decay1_d(tgi_i,tsi3_i,tau_i,tau1,dtt)
          if (fgaba.GT.0.D0) then
          call decay1(tgj_i,tsj3_i,tau2,dtt)
          endif
        else
          if (rk.EQ.4) then
            call decay3_d(tgi_i,tsi1_i,tsi2_i,tsi3_i,tau_i,tau1,dtt)
            if (fgaba.GT.0.D0) then
            call decay3(tgj_i,tsj1_i,tsj2_i,tsj3_i,tau2,dtt)
            endif
          endif
        endif
        do ijj=1,NpostSynE(ij)+NpostSynI(ij)
            i = IpostSyn(ij,ijj)
            if (unitary .eq. 0.5D0) then
                s = profile(Istrength(ij,ijj))
            else
                if (dabs(unitary-0.6D0).lt.1e-14)then
                    s = slist(Istrength(ij,ijj),i)
                else
                    s = 1.D0;
                endif
            endif
          if (excite(i)) then
        gi(i)  = gi(i)  + tgi_i * aa_ei(i)  *s
        gj(i)  = gj(i)  + tgj_i * aa_ei(i)  *s
        si1(i) = si1(i) + tsi1_i * aa_ei(i) *s
        sj1(i) = sj1(i) + tsj1_i * aa_ei(i) *s
        si2(i) = si2(i) + tsi2_i * aa_ei(i) *s
        sj2(i) = sj2(i) + tsj2_i * aa_ei(i) *s
        si3(i) = si3(i) + tsi3_i * aa_ei(i) *s
        sj3(i) = sj3(i) + tsj3_i * aa_ei(i) *s
          else
        gi(i)  = gi(i)  + tgi_i * aa_ii(i)  *s
        gj(i)  = gj(i)  + tgj_i * aa_ii(i)  *s
        si1(i) = si1(i) + tsi1_i * aa_ii(i) *s
        sj1(i) = sj1(i) + tsj1_i * aa_ii(i) *s
        si2(i) = si2(i) + tsi2_i * aa_ii(i) *s
        sj2(i) = sj2(i) + tsj2_i * aa_ii(i) *s
        si3(i) = si3(i) + tsi3_i * aa_ii(i) *s
        sj3(i) = sj3(i) + tsj3_i * aa_ii(i) *s
          endif
        enddo
      endif
!----------------------------------------------
        endif
        endif
!------------------------------------------------------------
      enddo
!======================================================= 
      endif 
!================ default decay are called in rk2 already 
!========================================== update all 
      do i=1,nmax
!=============================== to calc alpha1 and beta1
      if ( excite(i) ) then
      ginhib  = ((1.D0-fgaba)*gi(i) + fgaba*gj(i))*sei(i) + gEIc
      gexcit  = ((1.D0-fnmdacE)*ge(i) + fnmdacE*gf(i))*see(i)
      ginoise   = (1.D0-fgaba)*gn(i) + fgaba*gm(i)
      genoise   = (1.D0-fnmdanE)*gx(i) + fnmdanE*gy(i)
      gps = gp(i) * spotaE
      alpha1(i) = -gleak - gl(i) - ginoise - 
     1    gexcit - ginhib - genoise - gps
      beta1(i)  = (gl(i)+gexcit+genoise)*vexcit + 
     1    (ginoise+ginhib+gps)*vinhib + gleak*vreset
      suma = suma - alpha1(i)
      sumb = sumb + beta1(i)
      suml = suml + gl(i)
      sume = sume + gexcit
      sumi = sumi + ginhib
      sumin = sumin + ginoise
      sumen = sumen + genoise
      sump = sump + gps
      else
      ginhib  = ((1.D0-fgaba)*gi(i) + fgaba*gj(i))*sii(i) + gIIc 
      gexcit  = ((1.D0-fnmdacI)*ge(i) + fnmdacI*gf(i))*sie(i)
      ginoise   = (1.D0-fgaba)*gn(i) + fgaba*gm(i)
      genoise   = (1.D0-fnmdanI)*gx(i) + fnmdanI*gy(i)
      gps = gp(i) * spotaI
      alpha1(i) = -gleakI - gl(i) - ginoise - 
     1    gexcit - ginhib - genoise - gps
      beta1(i)  = (gl(i)+gexcit+genoise)*vexcit + 
     1    (ginoise+ginhib+gps)*vinhib + gleakI*vreset
      vsuma = vsuma - alpha1(i)
      vsumb = vsumb + beta1(i)
      vsuml = vsuml + gl(i)
      vsume = vsume + gexcit
      vsumi = vsumi + ginhib
      vsumin = vsumin + ginoise
      vsumen = vsumen + genoise
      vsump = vsump + gps
      endif
      enddo
!------------------------------------------------------------
      return
      end
!************************************************************
      subroutine e_or_i(excite,exc,iseed)
      use parameters
!------------------------------------------------------------
!  Set up excitatory/inhibitory tag
!  One Quarter of population is inhibitory
!    regular (or random) lattice
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      logical excite(nmax)
      dimension exc(nmax)
!------------------------------------------------------------
      do i=1,nex*ney
        exc(i)  = 1.D0
        excite(i) = .true.
      enddo
      do i=nex*ney+1,nmax
        exc(i) = 0.D0
        excite(i) = .false.
      enddo
!------------------------------------------------------------
      return
      end
!************************************************************
      subroutine visual(t,dt,iseed)
      use parameters
!------------------------------------------------------------
!     Also generate noise (1-1) given firing rates
!
!     Generates LGN spike times using Poisson process 
!    with time-dependent firing rate a Aunction of
!    visual stimulus
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension glo(nlgn), slo1(nlgn), slo2(nlgn), slo3(nlgn)
      dimension glon(nlgn),slon1(nlgn),slon2(nlgn),slon3(nlgn)
      dimension glf(nlgn), slf1(nlgn), slf2(nlgn), slf3(nlgn)
      dimension glfn(nlgn),slfn1(nlgn),slfn2(nlgn),slfn3(nlgn)
      dimension gx(nmax), sx1(nmax), sx2(nmax), sx3(nmax)
      dimension gy(nmax), sy1(nmax), sy2(nmax), sy3(nmax)
      dimension gn(nmax), sn1(nmax), sn2(nmax), sn3(nmax)
      dimension gm(nmax), sm1(nmax), sm2(nmax), sm3(nmax)
      dimension pspon(nlgn),pspoff(nlgn)
      dimension xlgn(nlgn),ylgn(nlgn)
      dimension tsonpoi(nlgn),ionpoi(nlgn)
      dimension tsoffpoi(nlgn),ioffpoi(nlgn)
      dimension tsenoi(nmax),ienoi(nmax)
      dimension tsinoi(nmax),iinoi(nmax)
      dimension nlgni(nmax)
      logical excite(nmax)
      integer iseed
      data iword / 8 /
      common / chaino / glo,slo1,slo2,slo3
      common / chaino2/ glon,slon1,slon2,slon3
      common / chainf / glf,slf1,slf2,slf3
      common / chainf2/ glfn,slfn1,slfn2,slfn3
      common / chainx / gx,sx1,sx2,sx3
      common / chainy / gy,sy1,sy2,sy3
      common / chainn / gn,sn1,sn2,sn3
      common / chainm / gm,sm1,sm2,sm3
      common / lgnIndivid / pspon,pspoff,gkx,gky,eps,satur,gtheta
      common / lgnGain / tmpphi,conAmp,g0,frtlgn0,treflgn,
     1   gphi, tstart, gfail,cond0
      common / lgnpos / xlgn,ylgn,jlgn
      common / noises / ce0_e,ci0_e,frtinh_e,frtexc_e,
     1   ce0_i,ci0_i,frtinh_i,frtexc_i
      common / tconst / tau_e,tau_i,tau0,tau1,tau2,tnrise,tndamp,
     1 tau_a0, tau_a1
      common /  NMDA  / fnmdatE,fnmdacE,fnmdatI,fnmdacI,
     1   fnmdanE, fnmdanI, fgaba
      common / neuron / excite,nlgni,trefE,trefI,nr,
     1   lgnmax,lgnmin,meanlgn
      common / lgnConvol / gk2, taulgn, omega
      common / ttotal / tfinal,twindow,rstart
      common / onpoi / tsonpoi,ionpoi,nonpoi
      common / offpoi / tsoffpoi,ioffpoi,noffpoi
      common / enoi / tsenoi,ienoi,nenoi
      common / inoi / tsinoi,iinoi,ninoi
      real*8 lgnsatur
      !INTEGER*4 nonsp,noffsp,nensp,ninsp
!$OMP THREADPRIVATE(/chaino/,/chaino2/,/chainf/,/chainf2/,
!$OMP& /chainx/,/chainy/,/chainn/,/chainm/,/lgnIndivid/,
!$OMP& /onpoi/,/offpoi/,/enoi/,/inoi/)
      iword2 = iword/2
      !omegat = omega*(t-tstart)
      omegat = omega*t

      nonpoi = 0
      noffpoi = 0
      ninoi = 0
      nenoi = 0

      do i=1,jlgn
       if (t.lt.0.D0) frate = frtlgn0
       rannum = -1.D0
       if (t+dt.gt.pspon(i)+treflgn) then
!-------------Find firing rate as function of visual stimulus
! For CONTRAST REVERSAL, need only amplitude of sinusoid
!============================================================
      if ((t.ge.0.D0).and.(t.lt.tstart)) then ! slowly ramp up
           !frate = frtlgn0 + (t/tstart)**2*eps*
           !frate = frtlgn0 + eps*
           frate = frtlgn0 + t/tstart*eps*
     1        conAmp*sin(omegat+tmpphi-gphi+xlgn(i)*gkx+ylgn(i)*gky)
      endif
      if (t.ge.tstart) then
           frate = frtlgn0 + eps*
     1        conAmp*sin(omegat+tmpphi-gphi+xlgn(i)*gkx+ylgn(i)*gky)
      endif
!------------------------Compute spikes only if non-zero rate
      if (frate.gt.0d0) then
!------------------------------------------------------------
! Contrast Saturation (200 levels, linear interpolate)
!------------------------------------------------------------
      satur =  lgnsatur(frate)
      frate = gfail*satur
!------------------------------------------------------------
      rannum = ran2(iseed)
      endif
       endif
!===========================================================
       dtt = rannum/frate
       if (rannum.lt.dt*frate.and.rannum.ge.0.D0
     1      .and.t+dtt.gt.pspon(i)+treflgn) then
!===========================================================
       pspon(i) = t + dtt
!=============================================================
      if (t.gt.tfinal-twindow) then
        nonpoi = nonpoi + 1
        tsonpoi(nonpoi) = t + dtt
        ionpoi(nonpoi) = i
      endif
!------------------------------------
       !lgnspikeo(i) = lgnspikeo(i)+1
!--------- t + dt -- dt

!--------- before spike | 0 --*--> dtt -----> dt
!--------- ampa exc

       call decay1_d(glo(i),slo3(i),tau_e,tau0,dtt)
       slo3(i) = slo3(i) + cond0
!--------- nmda exc 
       if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
         call decay1_d(glon(i),slon3(i),tnrise,tndamp,dtt)
         slon3(i) = slon3(i) + cond0*tau_e/tnrise
       endif 

!--------- after spike | 0 -----> dtt --*--> dt
       dtt = dt - dtt
       call decay1_d(glo(i),slo3(i),tau_e,tau0,dtt)
       if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
         call decay1_d(glon(i),slon3(i),tnrise,tndamp,dtt)
       endif
!--------- add by daiwei--------------------------------
       else
       dtt = dt
       call decay1_d(glo(i),slo3(i),tau_e,tau0,dtt)
       if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
         call decay1_d(glon(i),slon3(i),tnrise,tndamp,dtt)
       endif
!===========================================================
       endif
!===========================================================
      enddo
!----------------------------------Now off-centered LGN cells

      do i=1,jlgn
        if (t.lt.0.D0) frate = frtlgn0
        rannum = -1.D0
!===========================================================
      if (t+dt.gt.pspoff(i)+treflgn) then
!===========================================================
      if ((t.ge.0.D0).and.(t.lt.tstart)) then
         !frate = frtlgn0 - (t/tstart)**2*eps
         frate = frtlgn0 - t/tstart*eps*
     1   conAmp*sin(omegat+tmpphi-gphi+xlgn(i)*gkx+ylgn(i)*gky)
      endif
      if (t.ge.tstart) then
         frate = frtlgn0 - eps*
     1   conAmp*sin(omegat+tmpphi-gphi+xlgn(i)*gkx+ylgn(i)*gky)
      endif
!------------------------------------------------------------
        if (frate.gt.0.D0) then
!------------------------------------------------------------
! Contrast Saturation
!------------------------------------------------------------
      satur = lgnsatur(frate)
      frate = gfail*satur
      rannum = ran2(iseed)
        endif    
      endif
!===========================================================
      dtt = rannum/frate
      if (rannum.lt.dt*frate.and.rannum.ge.0.D0
     1   .and.t+dtt.gt.pspoff(i)+treflgn) then
!===========================================================
      !lgnspikef(i) = lgnspikef(i)+1
      pspoff(i) = t + dtt
!=============================================================
      if (t.gt.tfinal-twindow) then
        noffpoi = noffpoi + 1
        tsoffpoi(noffpoi) = t+dtt
        ioffpoi(noffpoi) = i
      endif
!------------------------------------
      call decay1_d(glf(i),slf3(i),tau_e,tau0,dtt)
      slf3(i) = slf3(i) + cond0
      if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
        call decay1_d(glfn(i),slfn3(i),tnrise,tndamp,dtt)
          slfn3(i) = slfn3(i) + cond0*tau_e/tnrise
      endif

      dtt = dt - dtt
      call decay1_d(glf(i),slf3(i),tau_e,tau0,dtt)
      if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
        call decay1_d(glfn(i),slfn3(i),tnrise,tndamp,dtt)
      endif
      else
      dtt = dt
      call decay1_d(glf(i),slf3(i),tau_e,tau0,dtt)
      if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
        call decay1_d(glfn(i),slfn3(i),tnrise,tndamp,dtt)
      endif
      endif
      enddo
!======================================== Exc noise units
      if (frtexc_e.gt.0.D0.and.ce0_e.gt.0.D0
     1  .or.frtexc_i.gt.0.D0.and.ce0_i.gt.0.D0)then
      do i=1,nmax
        rannum = ran2(iseed)
      if (excite(i)) then
          frtexc = frtexc_e
          ce0 = ce0_e
      else
          frtexc = frtexc_i
          ce0 = ce0_i
      endif
      if (rannum.lt.dt*frtexc.and.rannum.ge.0.D0) then
          dtt = rannum/frtexc
!------------------------------------------------------------
!  kick each input layer cell at tspike w/ delta fcn of 
!     strength Ce, Ci
!=============================================================
      if (t.gt.tfinal-twindow) then
        nenoi = nenoi + 1
        tsenoi(nenoi) = t + dtt
        ienoi(nenoi) = i
      endif
!------------------------------------
      call decay1_d(gx(i),sx3(i),tau_e,tau0,dtt)
      sx3(i) = sx3(i) + ce0/tau_e
        if (fnmdanE.GT.0.D0.or.fnmdanI.GT.0.D0) then
          call decay1_d(gy(i),sy3(i),tnrise,tndamp,dtt)
          sy3(i) = sy3(i) + ce0/tnrise
        endif
      dtt = dt - dtt
      call decay1_d(gx(i),sx3(i),tau_e,tau0,dtt)
        if (fnmdanE.GT.0.D0.or.fnmdanI.GT.0.D0) then
          call decay1_d(gy(i),sy3(i),tnrise,tndamp,dtt)
        endif
      else
      dtt = dt
      call decay1_d(gx(i),sx3(i),tau_e,tau0,dtt)
        if (fnmdanE.GT.0.D0.or.fnmdanI.GT.0.D0) then
          call decay1_d(gy(i),sy3(i),tnrise,tndamp,dtt)
        endif
      endif
      enddo
      endif
!======================================= Inh noise units
      if (frtinh_e.gt.0.D0.and.ci0_e.gt.0.D0
     1 .or.frtinh_i.gt.0.D0.and.ci0_i.gt.0.D0)then
      do i=1,nmax
          rannum = ran2(iseed)
      if (excite(i)) then
          frtinh = frtinh_e
          ci0 = ci0_e
      else
          frtinh = frtinh_i
          ci0 = ci0_i
      endif
      if (rannum.lt.dt*frtinh.and.rannum.ge.0.D0) then
          dtt = rannum/frtinh
!=============================================================
      if (t.gt.tfinal-twindow) then
        ninoi = ninoi + 1
        tsinoi(ninoi) = t + dtt
        iinoi(ninoi) = i
      endif
!------------------------------------
      call decay1_d(gn(i),sn3(i),tau_i,tau1,dtt)
      sn3(i) = sn3(i) + ci0/tau_i
      if (fgaba.GT.0.D0) then
        call decay1(gm(i),sm3(i),tau2,dtt)
        sm3(i) = sm3(i) + ci0/tau2
      endif

      dtt = dt - dtt 
      call decay1_d(gn(i),sn3(i),tau_i,tau1,dtt)
        if (fgaba.GT.0.D0) then
          call decay1(gm(i),sm3(i),tau2,dtt)
        endif
      else
      dtt = dt
      call decay1_d(gn(i),sn3(i),tau_i,tau1,dtt)
        if (fgaba.GT.0.D0) then
          call decay1(gm(i),sm3(i),tau2,dtt)
        endif
      endif
      enddo
      endif
!------------------------------------------------------------
      return
      end
!************************************************************
      subroutine chain1(ge,se,tau,deltat,nn)
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension ge(nn),se(nn)
      te  = deltat/tau
      ete = exp(-te)
      do i=1,nn
        ge(i) = (ge(i) + se(i)*te) * ete
        se(i) = se(i) * ete
      enddo
      return
      end
      subroutine chain1_d(ge,se,trise,tdamp,deltat,nn)
!------------------------------------------------------------
!  Oct 2000: Updating difference of expon'tial w/o intracortical spikes
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension ge(nn),se(nn)
!------------------------------------------------------------
      tr  = deltat/trise
      etr = exp(-tr)

      td  = deltat/tdamp
      etd = exp(-td)

      const = trise/(tdamp - trise) * (etd - etr)

      do i=1,nn
        ge(i)  =  ge(i) * etd + const * se(i)
        se(i) = se(i) * etr
      enddo
!------------------------------------------------------------
      return
      end
!------------------------------------------------------------
      subroutine chain3(ge,se1,se2,se3,tau,deltat,nn)
!------------------------------------------------------------
!  Oct 2000: Updating chain w/o intracortical spikes
!    PSP handled in rk2/rk4 routine
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension ge(nn),se1(nn),se2(nn),se3(nn)
!------------------------------------------------------------
      te  = deltat/tau
      te2 = te*te/2.D0
      te3 = te*te2/3.D0
      ete = exp(-te)

      do i=1,nn
!      do i=1,nmax
        ge(i)  = (ge(i) + se1(i)*te + se2(i)*te2
     1         +  se3(i)*te3)*ete
        se1(i) = (se1(i) + se2(i)*te + se3(i)*te2)*ete
        se2(i) = (se2(i) + se3(i)*te)*ete
        se3(i) = (se3(i)            )*ete
      enddo
!------------------------------------------------------------
      return
      end
!------------------------------------------------------------
      subroutine chain3_d(ge,se1,se2,se3,trise,tdamp,deltat,nn)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension ge(nn),se1(nn),se2(nn),se3(nn)
      td = deltat/tdamp
      tr = deltat/trise
      t2 = deltat**2
      etr = exp(-tr)
      etd = exp(-td)
      tdamp2 = tdamp**2
      trise2 = trise**2
      tdiff = tdamp-trise
      tmulti = tdamp*trise
      con1 = trise/tdiff * (etd - etr)
      con2 = (etd*tmulti - etr*(deltat*tdiff+tmulti))/tdiff**2

      do i=1,nn
      ge(i) = se3(i)*(2.D0*tmulti**2*etd -  
     1      (t2*(tdamp2 + trise2) - 
     2      (t2-tmulti-deltat*tdiff)*2*tmulti)*etr)/
     3      (2.D0*trise*tdiff**3)+se1(i)*con1+se2(i)*con2+ge(i)*etd

      se1(i) = (se3(i) * 0.5D0*tr*tr + se2(i) * tr + se1(i)) * etr
      se2(i) = (se3(i) * tr  + se2(i)) * etr
      se3(i) = se3(i) * etr 
      enddo
      return
      end
!************************************************************
!    single neuron version of chain and chain2
!************************************************************
      subroutine decay1(ge,se,tau,deltat)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      te  = deltat/tau
      ete = exp(-te)
      ge  = (ge + se * te) *ete
      se  = se * ete
      return
      end
      subroutine decay1_d(ge,se,trise,tdamp,deltat)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      !print *,'before ge ', ge,'se ',se
      tr  = deltat/trise
      etr = exp(-tr)

      td  = deltat/tdamp
      etd = exp(-td)

      const = trise/(tdamp - trise) * (etd - etr)

        ge  =  ge * etd + const * se
        se  =  se * etr
      !print *,'after ge ', ge,'se ',se
      return
      end
      subroutine decay3(ge,se1,se2,se3,tau,deltat)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      te  = deltat/tau
      te2 = te*te/2.D0
      te3 = te*te2/3.D0
      ete = exp(-te)

        ge  = (ge + se1*te + se2*te2
     1         +  se3*te3)*ete
        se1 = (se1 + se2*te + se3*te2)*ete
        se2 = (se2 + se3*te)*ete
        se3 = se3 * ete
      return
      end
      subroutine decay3_d(ge,se1,se2,se3,trise,tdamp,deltat)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 

      td = deltat/tdamp
      tr = deltat/trise
      t2 = deltat**2
      etr = exp(-tr)
      etd = exp(-td)
      tdamp2 = tdamp**2
      trise2 = trise**2
      tdiff = tdamp-trise
      tmulti = tdamp*trise
      con1 = trise/tdiff * (etd - etr)
      con2 = (etd*tmulti - etr*(deltat*tdiff+tmulti))/tdiff**2

      ge = se3*(2.D0*tmulti**2*etd -  
     1     (t2*(tdamp2 + trise2) - 
     2     (t2-tmulti-deltat*tdiff)*2*tmulti)*etr)/
     3     (2.D0*trise*tdiff**3)+se1*con1+se2*con2+ge*etd

      se1 = (se3 * 0.5D0*tr*tr + se2 * tr + se1) * etr
      se2 = (se3 * tr  + se2) * etr
      se3 = se3 * etr 
      return
      end
!************************************************************
      subroutine eif_rk2(t,dt,period,iperiod,iid,idap)
      use parameters
!------------------------------------------------------------
!  Oct 2000
!  uses chain to evaluate PSP (including LGN input & Noise)
!
!  Modified RK2 to solve
!    d v / dt = alpha(t) v + beta(t)
!  when v = vthres it becomes vreset
!------------------------------------------------------------
!  April 2014
!* Questions: 
!     No spikes, i.e., no conductance rise from lgn, noise when estimate spike time
!     Whatif large conductance lead to vnew after spike still larger than vthres
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      dimension v(nmax),vnew(nmax),gl(nmax)
      dimension alpha1(nmax),beta1(nmax)
      dimension alpha0(nmax),beta0(nmax)
      dimension glo(nlgn),slo1(nlgn),slo2(nlgn),slo3(nlgn)
      dimension glon(nlgn),slon1(nlgn),slon2(nlgn),slon3(nlgn)
      dimension glf(nlgn),slf1(nlgn),slf2(nlgn),slf3(nlgn)
      dimension glfn(nlgn),slfn1(nlgn),slfn2(nlgn),slfn3(nlgn)
      dimension ge(nmax),se1(nmax),se2(nmax),se3(nmax)
      dimension gf(nmax),sf1(nmax),sf2(nmax),sf3(nmax)
      dimension gi(nmax),si1(nmax),si2(nmax),si3(nmax)
      dimension gj(nmax),sj1(nmax),sj2(nmax),sj3(nmax)
      dimension gx(nmax),sx1(nmax),sx2(nmax),sx3(nmax)
      dimension gy(nmax),sy1(nmax),sy2(nmax),sy3(nmax)
      dimension gn(nmax),sn1(nmax),sn2(nmax),sn3(nmax)
      dimension gm(nmax),sm1(nmax),sm2(nmax),sm3(nmax)
      dimension gp(nmax),sp1(nmax),sp2(nmax),sp3(nmax)
      dimension ge0(nmax),se30(nmax)
      dimension gf0(nmax),sf30(nmax)
      dimension gi0(nmax),si30(nmax)
      dimension gj0(nmax),sj30(nmax)
      dimension gx0(nmax),sx30(nmax)
      dimension gy0(nmax),sy30(nmax)
      dimension gn0(nmax),sn30(nmax)
      dimension gm0(nmax),sm30(nmax)
      dimension gp0(nmax),sp30(nmax)
      dimension glnmda(nmax),glampa(nmax)
      dimension glnmda2(nmax),glampa2(nmax)
      dimension nonlgn(nmax),noflgn(nmax)
      dimension ionlgn(nmax,nlgn),ioflgn(nmax,nlgn)
      dimension sonlgn(nmax,nlgn),soflgn(nmax,nlgn)
      dimension tspike(nmax),ispike(nmax),nlgni(nmax),pspike(nmax)
      dimension irate(nmax,25),erate(nex*ney)
      dimension see(nmax),sei(nmax),sie(nmax),sii(nmax)
      logical excite(nmax)
      !integer OMP_GET_THREAD_NUM
!----THREADPRIVATE----------------------------------------------
      common / indRate /erate
      common / spikes / tspike,ispike,nspike
      common / persp /  irate
      common / lgncnd / gl
      common /   rhs  / alpha1,beta1
      common / chaino / glo,slo1,slo2,slo3
      common / chaino2/ glon,slon1,slon2,slon3
      common / chainf / glf,slf1,slf2,slf3
      common / chainf2/ glfn,slfn1,slfn2,slfn3
      common / chaine / ge,se1,se2,se3
      common / chaine2/ gf,sf1,sf2,sf3
      common / chaini / gi,si1,si2,si3
      common / chainj / gj,sj1,sj2,sj3
      common / chainx / gx,sx1,sx2,sx3
      common / chainy / gy,sy1,sy2,sy3
      common / chainn / gn,sn1,sn2,sn3
      common / chainm / gm,sm1,sm2,sm3
      common / chainp / gp,sp1,sp2,sp3
      common /   isi  / pspike,nsptot,isptot,jsptot,ksptot,msptot
      common /  volt  / v, vnew
!----SHARED-----------------------------------------------
      common / smatrx / see,sei,sie,sii,spotaE,spotaI
      common / ttotal / tfinal,twindow,rstart
      common / neuron / excite,nlgni,trefE,trefI,nr,
     1   lgnmax,lgnmin,meanlgn
      common / lgnmap / sonlgn,soflgn,ionlgn,ioflgn,
     1   nonlgn, noflgn
      common / ggIc / gIIc, gEIc
      common / vconst / vthres,vreset,vexcit,vinhib,gleak,gleakI
      common / tconst / tau_e,tau_i,tau0,tau1,tau2,tnrise,tndamp,
     1 tau_a0, tau_a1
      common /  NMDA  / fnmdatE,fnmdacE,fnmdatI,fnmdacI,
     1   fnmdanE, fnmdanI, fgaba
      common / prestatics/ npEE, npEI, npIE, npII
      common / EIFconst / DeltaT,vTheta,phiVAS
!$OMP THREADPRIVATE(/chaino/,/chaino2/,/chainf/,/chainf2/,/chainp/,
!$OMP& /chaine/,/chaine2/,/chaini/,/chainj/,/chainx/,/persp/,
!$OMP& /chainy/,/chainn/,/chainm/,/rhs/,/isi/,/lgncnd/,/volt/,
!$OMP& /spikes/,/indRate/)
!------------------------------------------------------------
!  First update total conductance & total current
!------------------------------------------------------------
      do i=1,nmax
        alpha0(i) = alpha1(i)
        beta0(i) = beta1(i)
        glnmda(i) = 0.D0
        glampa(i) = 0.D0
        glnmda2(i) = 0.D0
        glampa2(i) = 0.D0
      do j=1,nonlgn(i)
        glnmda(i) = glnmda(i) + glon(ionlgn(i,j))*sonlgn(i,j)
        glampa(i) = glampa(i) + glo(ionlgn(i,j))*sonlgn(i,j)
        glnmda2(i) = glnmda2(i) + slon3(ionlgn(i,j))*sonlgn(i,j)
        glampa2(i) = glampa2(i) + slo3(ionlgn(i,j))*sonlgn(i,j)
      enddo
      do j=1,noflgn(i)
        glnmda(i) = glnmda(i) + glfn(ioflgn(i,j))*soflgn(i,j)
        glampa(i) = glampa(i) + glf(ioflgn(i,j))*soflgn(i,j)
        glnmda2(i) = glnmda2(i) + slfn3(ioflgn(i,j))*soflgn(i,j)
        glampa2(i) = glampa2(i) + slf3(ioflgn(i,j))*soflgn(i,j)
      enddo
        ge0(i) = ge(i)
        se30(i) = se3(i)
        gf0(i) = gf(i)
        sf30(i) = sf3(i)
        gi0(i) = gi(i)
        si30(i) = si3(i)
        gj0(i) = gj(i)
        sj30(i) = sj3(i)
        gx0(i) = gx(i)
        sx30(i) = sx3(i)
        gy0(i) = gy(i)
        sy30(i) = sy3(i)
        gn0(i) = gn(i)
        sn30(i) = sn3(i)
        gm0(i) = gm(i)
        sm30(i) = sm3(i)
        gp0(i) = gp(i)
        sp30(i) = sp3(i)
      enddo
! thalamic and noises
      call visual(t,dt,iseed)
! ampa cortical
      call chain1_d(ge,se3,tau_e,tau0,dt,nmax)
! nmda cortical 
        if(fnmdacE.GT.0.D0.or.fnmdacI.GT.0.D0) then
      call chain1_d(gf,sf3,tnrise,tndamp,dt,nmax)
        endif

! gabaA cortical
      call chain1_d(gi,si3,tau_i,tau1,dt,nmax)
! gaba2 cortical
        if (fgaba.GT.0.D0)then
      call chain1(gj,sj3,tau2,dt,nmax)
        endif
! potassium 
        if (idap.eq.1) then
      call chain1_d(gp,sp3,tau_a0,tau_a1,dt,nmax)
        endif
!!----- add up conductance to calculate alpha1 and beta1
      do i=1,nmax
      gl(i) = 0.D0
      if ( excite(i) ) then
      do j=1,nonlgn(i)
      gl(i) = gl(i) + (1-fnmdatE)*glo(ionlgn(i,j))*sonlgn(i,j)
     1        + fnmdatE * glon(ionlgn(i,j))*sonlgn(i,j)
      enddo
      do j=1,noflgn(i)
      gl(i) = gl(i) + (1-fnmdatE)*glf(ioflgn(i,j))*soflgn(i,j)
     1        + fnmdatE * glfn(ioflgn(i,j))*soflgn(i,j)
      enddo

      ginhib  = ((1-fgaba)*gi(i) + fgaba*gj(i))*sei(i) + gEIc
      gexcit  = ((1-fnmdacE)*ge(i) + fnmdacE*gf(i))*see(i)
      ginoise   = (1-fgaba)*gn(i) + fgaba*gm(i)
      genoise   = (1-fnmdanE)*gx(i) + fnmdanE*gy(i)
      gps = gp(i)*spotaE
      gleakk = gleak
      else
      do j=1,nonlgn(i)
      gl(i) = gl(i) + (1-fnmdatI)*glo(ionlgn(i,j))*sonlgn(i,j)
     1        + fnmdatI * glon(ionlgn(i,j))*sonlgn(i,j)
      enddo
      do j=1,noflgn(i)
      gl(i) = gl(i) + (1-fnmdatI)*glf(ioflgn(i,j))*soflgn(i,j)
     1        + fnmdatI * glfn(ioflgn(i,j))*soflgn(i,j)
      enddo
      ginhib  = ((1-fgaba)*gi(i) + fgaba*gj(i))*sii(i) + gIIc
      gexcit  = ((1-fnmdacI)*ge(i) + fnmdacI*gf(i))*sie(i)
      ginoise   = (1-fgaba)*gn(i) + fgaba*gm(i)
      genoise   = (1-fnmdanI)*gx(i) + fnmdanI*gy(i)
      gps = gp(i)*spotaI
      gleakk = gleakI
      endif
      alpha1(i) = -gleakk - gl(i) - ginoise - gexcit - 
     1   ginhib - genoise - gps
      beta1(i)  = (gl(i)+gexcit+genoise)*vexcit + 
     1   (ginoise+ginhib+gps)*vinhib + gleakk*vreset
      enddo
!----- rk2 begins
      dt2 = 0.5D0*dt
      nspike = 0
      do i=1,nmax
        a1 = alpha1(i)
        b10 = beta1(i)
        vn = v(i)
        dt0 = dt
!-------come back from refractory -daiwei----------------------
        if ( excite(i) ) then
          tref = trefE
        else
          tref = trefI
        endif

      tBack = pspike(i)+tref-t
      if (tBack.lt.dt.and.tBack.gt.0.D0)then
 
! ampa thalamic
      call decay1_d(glampa(i),glampa2(i),tau_e,tau0,tBack)
! nmda thalamic
        if(fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
      call decay1_d(glnmda(i),glnmda2(i),tnrise,tndamp,tBack)
        endif

! ampa cortical
      call decay1_d(ge0(i),se30(i),tau_e,tau0,tBack)
! nmda cortical 
        if(fnmdacE.GT.0.D0.or.fnmdacI.GT.0.D0) then
      call decay1_d(gf0(i),sf30(i),tnrise,tndamp,tBack)
        endif
! ampa noise 
      call decay1_d(gx0(i),sx30(i),tau_e,tau0,tBack)
! nmda noise 
        if(fnmdacE.GT.0.D0.or.fnmdacI.GT.0.D0) then
      call decay1_d(gy0(i),sy30(i),tnrise,tndamp,tBack)
        endif

! gabaA cortical
      call decay1_d(gi0(i),si30(i),tau_i,tau1,tBack)
! gaba2 cortical
        if (fgaba.GT.0.D0)then
      call decay1(gj0(i),sj30(i),tau2,tBack)
        endif
! gabaA noise 
      call decay1_d(gn0(i),sn30(i),tau_i,tau1,tBack)
! gaba2 noise 
        if (fgaba.GT.0.D0)then
      call decay1(gm0(i),sm30(i),tau2,tBack)
        endif
! potassium
        if (idap.eq.1) then
      call decay1_d(gp0(i),sp30(i),tau_a0,tau_a1,tBack)
        endif
!---
      gl0      = 0.D0
      if ( excite(i) ) then
      gl0 = gl0 + (1-fnmdatE)*glampa(i) + fnmdatE*glnmda(i)
      ginhib  = ((1-fgaba)*gi0(i) + fgaba*gj0(i))*sei(i) + gEIc
      gexcit  = ((1-fnmdacE)*ge0(i) + fnmdacE*gf0(i))*see(i)
      ginoise = (1-fgaba)*gn0(i) + fgaba*gm0(i)
      genoise = (1-fnmdanE)*gx0(i) + fnmdanE*gy0(i)
      gps = gp0(i)*spotaE
      gleakk = gleak;
      else
      gl0 = gl0 + (1-fnmdatI)*glampa(i) + fnmdatI*glnmda(i)
      ginhib  = ((1-fgaba)*gi0(i) + fgaba*gj0(i))*sii(i) + gIIc
      gexcit  = ((1-fnmdacI)*ge0(i) + fnmdacI*gf0(i))*sie(i)
      ginoise = (1-fgaba)*gn0(i) + fgaba*gm0(i)
      genoise = (1-fnmdanI)*gx0(i) + fnmdanI*gy0(i)
      gps = gp0(i)*spotaI
      gleakk = gleakI;
      endif
      a0 = -gleakk - gl0 - ginoise - gexcit - 
     1  ginhib - genoise - gps
      b0 = (gl0+gexcit+genoise)*vexcit + 
     1  (ginoise+ginhib+gps)*vinhib + gleakk*(vreset+phiVAS)

        vn  = vreset
        dt0 = dt - tBack
        fk1 = a0*vn + b0
        v1 = vn + dt0*fk1
        phiV = gleakk*DeltaT*dexp((v1-vthres)/DeltaT)
        b1 = b10 + phiV
        fk2 = a1*v1 + b1
        vnew(i) = vn + dt0/2.D0*(fk1+fk2)
      else
        if (tBack.gt.dt) then
          vnew(i) = vreset
!------ else normal rk2 ----------------------------------
        else
          tBack = 0.D0
          a0 = alpha0(i)
          phiV = gleakk*DeltaT*dexp((vn-vthres)/DeltaT)
          b0 = beta0(i) + phiV
          fk1 = a0*vn + b0 
          v1 = vn + dt*fk1
          phiV = gleakk*DeltaT*dexp((v1-vthres)/DeltaT)
          b1 = b10 + phiV
          fk2 = a1*v1 + b1
          vnew(i) = vn + dt2*(fk1+fk2)
        endif
      endif
!$OMP CRITICAL
      if (v(i).lt.vinhib.and.vnew(i).lt.v(i)) then
          print *, 'wtf', iid, i
          if (excite(i)) then
             ginhib  = ((1-fgaba)*gi(i) + fgaba*gj(i))*sei(i) + gEIc
             gexcit  = ((1-fnmdacE)*ge(i) + fnmdacE*gf(i))*see(i)
             print *, npEI*isptot/(nix*niy),' Inh spikes -> 1 Exc'
             print *, npEE*nsptot/(nex*ney),' Exc spikes -> 1 Exc'
          else
             ginhib  = ((1-fgaba)*gi(i) + fgaba*gj(i))*sii(i) + gIIc
             gexcit  = ((1-fnmdacI)*ge(i) + fnmdacI*gf(i))*sie(i)
             print *, npEI*isptot/(nix*niy),' Inh spikes -> 1 Inh'
             print *, npEE*nsptot/(nex*ney),' Exc spikes -> 1 Inh'
          endif
          print *, gl(i),gexcit,ginhib
          print *, v(i), vnew(i)
          print *, a0, b0
          print *, a1, b1
          print *, fk1, fk2
          print *, gps
          stop
      endif
!$OMP END CRITICAL
!---------------------------------------------if neuron fires
      if (vnew(i).gt.vTheta) then
!------------------------------------------------------------
!  Mike's modified RK2
!    1. estimate spike time (interpolate linearly)
!    2. calculate new init cond for the reset (extrapolate)
!    3. calculate vnew after spike (retake rk2 step)
!------------------------------------------------------------
      dtsp = tBack + (vTheta-vn)/(vnew(i)-vn)*dt0
!---------------------------------------------------------
!  Refractory Period (Wei, 4/28/2014)
!--------------------------------------------------------
      pspike(i) = t+dtsp
      if (pspike(i) + tref.lt. t+dt) then
 
! ampa thalamic
      call decay1_d(glampa(i),glampa2(i),tau_e,tau0,dtsp-tBack)
! nmda thalamic
        if(fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
      call decay1_d(glnmda(i),glnmda2(i),tnrise,tndamp,dtsp-tBack)
        endif

! ampa cortical
      call decay1_d(ge0(i),se30(i),tau_e,tau0,dtsp-tBack)
! nmda cortical 
        if(fnmdacE.GT.0.D0.or.fnmdacI.GT.0.D0) then
      call decay1_d(gf0(i),sf30(i),tnrise,tndamp,dtsp-tBack)
        endif
! ampa noise 
      call decay1_d(gx0(i),sx30(i),tau_e,tau0,dtsp-tBack)
! nmda noise 
        if(fnmdacE.GT.0.D0.or.fnmdacI.GT.0.D0) then
      call decay1_d(gy0(i),sy30(i),tnrise,tndamp,dtsp-tBack)
        endif

! gabaA cortical
      call decay1_d(gi0(i),si30(i),tau_i,tau1,dtsp-tBack)
! gaba2 cortical
        if (fgaba.GT.0.D0)then
      call decay1(gj0(i),sj30(i),tau2,dtsp-tBack)
        endif
! gabaA noise 
      call decay1_d(gn0(i),sn30(i),tau_i,tau1,dtsp-tBack)
! gaba2 noise 
        if (fgaba.GT.0.D0)then
      call decay1(gm0(i),sm30(i),tau2,dtsp-tBack)
        endif
! potassium
        if (idap.eq.1) then
      call decay1_d(gp0(i),sp30(i),tau_a0,tau_a1,dtsp-tBack)
        endif
!---

      gl0      = 0.D0
      if ( excite(i) ) then
      gl0 = gl0 + (1-fnmdatE)*glampa(i) + fnmdatE*glnmda(i)
      ginhib  = ((1-fgaba)*gi0(i) + fgaba*gj0(i))*sei(i) + gEIc
      gexcit  = ((1-fnmdacE)*ge0(i) + fnmdacE*gf0(i))*see(i)
      ginoise = (1-fgaba)*gn0(i) + fgaba*gm0(i)
      genoise = (1-fnmdanE)*gx0(i) + fnmdanE*gy0(i)
      gps = gp0(i)*spotaE
      gleakk = gleak
      else
      gl0 = gl0 + (1-fnmdatI)*glampa(i) + fnmdatI*glnmda(i)
      ginhib  = ((1-fgaba)*gi0(i) + fgaba*gj0(i))*sii(i) + gIIc
      gexcit  = ((1-fnmdacI)*ge0(i) + fnmdacI*gf0(i))*sie(i)
      ginoise = (1-fgaba)*gn0(i) + fgaba*gm0(i)
      genoise = (1-fnmdanI)*gx0(i) + fnmdanI*gy0(i)
      gps = gp0(i)*spotaI
      gleakk = gleakI
      endif
      a0 = -gleakk - gl0 - ginoise - gexcit - 
     1  ginhib - genoise - gps
      b0 = (gl0+gexcit+genoise)*vexcit + 
     1  (ginoise+ginhib+gps)*vinhib + gleakk*(vreset+phiVAS)

        vn  = vreset
        dt0 = dt - dtsp
        fk1 = a0*vn + b0
        v1 = vn + dt0*fk1
        phiV = gleakk*DeltaT*dexp((v1-vthres)/DeltaT)
        b1 = b10 + phiV
        fk2 = a1*v1 + b1
        vnew(i) = vn + dt0/2.D0*(fk1+fk2)
      else
        vnew(i) = vreset
      endif
!  debug daiwei
      if (vnew(i).gt.vTheta) print *,'2 spikes per dt! i:',
     1         i,t+dtsp,'v0:',v(i),'vnew:',vnew(i)
      nspike = nspike + 1
      tspike(nspike) = dtsp
      ispike(nspike) = i    
      if ( .not. excite(i) ) then
          isptot = isptot + 1
      else 
          erate(i) = erate(i) + 1.D0
          nsptot = nsptot + 1
      endif
      if (nlgni(i).EQ.lgnmin.and.excite(i)) jsptot = jsptot + 1
      if (nlgni(i).GE.lgnmax.and.excite(i)) ksptot = ksptot + 1
      if (nlgni(i).EQ.meanlgn.and.excite(i)) msptot = msptot + 1
!------------------------------------------------------------
!  Construct histogram for spike rate: ith neuron nspikes at nirate th bin 
!------------------------------------------------------------
      if (t.ge.rstart) then
        tirate = dmod(t+dtsp,period)
        nirate = floor(tirate*1.D0*iperiod/period) + 1
        irate(i,nirate) = irate(i,nirate) + 1
      endif
!------------------------------------------endif neuron fires
      endif
!------------------------------------------------------------
      enddo
!------------------------------------------------------------
      return
      end
!************************************************************
      subroutine rk4(t,dt,period,iperiod,iid)
      use parameters
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      dimension v(nmax),vnew(nmax),gl(nmax)
      dimension alpha0(nmax),alpha1(nmax),beta0(nmax),beta1(nmax)
      dimension alpha2(nmax),beta2(nmax)
      dimension glo(nlgn), slo1(nlgn), slo2(nlgn), slo3(nlgn)
      dimension glon(nlgn),slon1(nlgn),slon2(nlgn),slon3(nlgn)
      dimension glf(nlgn), slf1(nlgn), slf2(nlgn), slf3(nlgn)
      dimension glfn(nlgn),slfn1(nlgn),slfn2(nlgn),slfn3(nlgn)
      dimension gx(nmax), sx1(nmax), sx2(nmax), sx3(nmax)
      dimension gy(nmax), sy1(nmax), sy2(nmax), sy3(nmax)
      dimension gn(nmax), sn1(nmax), sn2(nmax), sn3(nmax)
      dimension gm(nmax), sm1(nmax), sm2(nmax), sm3(nmax)
      dimension gl2o(nlgn), sl2o1(nlgn), sl2o2(nlgn), sl2o3(nlgn)
      dimension gl2on(nlgn),sl2on1(nlgn),sl2on2(nlgn),sl2on3(nlgn)
      dimension gl2f(nlgn), sl2f1(nlgn), sl2f2(nlgn), sl2f3(nlgn)
      dimension gl2fn(nlgn),sl2fn1(nlgn),sl2fn2(nlgn),sl2fn3(nlgn)
      dimension g2x(nmax), s2x1(nmax), s2x2(nmax), s2x3(nmax)
      dimension g2y(nmax), s2y1(nmax), s2y2(nmax), s2y3(nmax)
      dimension g2n(nmax), s2n1(nmax), s2n2(nmax), s2n3(nmax)
      dimension g2m(nmax), s2m1(nmax), s2m2(nmax), s2m3(nmax)
      dimension ge(nmax),se1(nmax),se2(nmax),se3(nmax)
      dimension gf(nmax),sf1(nmax),sf2(nmax),sf3(nmax)
      dimension gi(nmax),si1(nmax),si2(nmax),si3(nmax)
      dimension gj(nmax),sj1(nmax),sj2(nmax),sj3(nmax)
      dimension nonlgn(nmax),noflgn(nmax)
      dimension ionlgn(nmax,nlgn),ioflgn(nmax,nlgn)
      dimension sonlgn(nmax,nlgn),soflgn(nmax,nlgn)
      dimension tspike(nmax),ispike(nmax),nlgni(nmax),pspike(nmax)
      dimension irate(nmax,25)
      dimension see(nmax),sei(nmax),sie(nmax),sii(nmax)
      dimension erate(nex*ney)
      logical excite(nmax)
      !integer OMP_GET_THREAD_NUM
!----THREADPRIVATE----------------------------------------------
      common / indRate / erate
      common / spikes / tspike,ispike,nspike
      common / persp /  irate
      common / lgncnd / gl
      common /   rhs  / alpha1,beta1
      common / chaino / glo, slo1, slo2, slo3
      common / chaino2/ glon,slon1,slon2,slon3
      common / chainf / glf, slf1, slf2, slf3
      common / chainf2/ glfn,slfn1,slfn2,slfn3
      common / chain2o / gl2o, sl2o1, sl2o2, sl2o3
      common / chain2o2/ gl2on,sl2on1,sl2on2,sl2on3
      common / chain2f / gl2f, sl2f1, sl2f2, sl2f3
      common / chain2f2/ gl2fn,sl2fn1,sl2fn2,sl2fn3
      common / chainx / gx,sx1,sx2,sx3
      common / chainy / gy,sy1,sy2,sy3
      common / chainn / gn,sn1,sn2,sn3
      common / chainm / gm,sm1,sm2,sm3
      common / chain2x / g2x,s2x1,s2x2,s2x3
      common / chain2y / g2y,s2y1,s2y2,s2y3
      common / chain2n / g2n,s2n1,s2n2,s2n3
      common / chain2m / g2m,s2m1,s2m2,s2m3
      common / chaine / ge,se1,se2,se3
      common / chaine2/ gf,sf1,sf2,sf3
      common / chaini / gi,si1,si2,si3
      common / chainj / gj,sj1,sj2,sj3
      common /   isi  / pspike,nsptot,isptot,jsptot,ksptot,msptot
      common /  volt  / v, vnew
!----SHARED-----------------------------------------------
      common / smatrx / see,sei,sie,sii,spotaE,spotaI
      common / ttotal / tfinal,twindow,rstart
      common / neuron / excite,nlgni,trefE,trefI,nr,
     1   lgnmax,lgnmin,meanlgn
      common / lgnmap / sonlgn,soflgn,ionlgn,ioflgn,
     1   nonlgn, noflgn
      common / ggIc / gIIc, gEIc
      common / vconst / vthres,vreset,vexcit,vinhib,gleak,gleakI
      common / tconst / tau_e,tau_i,tau0,tau1,tau2,tnrise,tndamp,
     1 tau_a0, tau_a1
      common /  NMDA  / fnmdatE,fnmdacE,fnmdatI,fnmdacI,
     1   fnmdanE, fnmdanI, fgaba
      common / prestatics/ npEE, npEI, npIE, npII
!$OMP THREADPRIVATE(/indRate/,
!$OMP& /chaino/, /chaino2/, /chainf/, /chainf2/,
!$OMP& /chainx/, /chainy/,  /chainn/, /chainm/,
!$OMP& /chain2o/,/chain2o2/,/chain2f/,/chain2f2/,
!$OMP& /chain2x/,/chain2y/, /chain2n/,/chain2m/,
!$OMP& /chaine/,/chaine2/,/chaini/,/chainj/,/volt/,
!$OMP& /spikes/,/persp/,/rhs/,/isi/,/lgncnd/)

!===== t ======================
      itmax = 100
      dt2 = dt/2.D0
      dt6 = dt/6.D0
      do i=1,nmax
        alpha0(i) = alpha1(i)
        beta0(i)  = beta1(i)
      enddo

!===== t+0.D5dt ======================
! ampa cortical
      call chain3_d(ge,se1,se2,se3,tau_e,tau0,dt2,nmax)
! nmda cortical 
        if(fnmdacE.GT.0.D0.or.fnmdacI.GT.0.D0) then
      call chain3_d(gf,sf1,sf2,sf3,tnrise,tndamp,dt2,nmax)
        endif

! gabaA cortical
      call chain3_d(gi,si1,si2,si3,tau_i,tau1,dt2,nmax)
! gaba2 cortical
        if (fgaba.GT.0.D0)then
      call chain3(gj,sj1,sj2,sj3,tau2,dt2,nmax)
        endif

!!----- add up conductance to calculate alpha2 and beta2
      do i=1,nmax
      gl(i)      = 0.D0
      if ( excite(i) ) then
      do j=1,nonlgn(i)
      gl(i) = gl(i) + (1.D0-fnmdatE)*gl2o(ionlgn(i,j))*sonlgn(i,j)
     1        + fnmdatE * gl2on(ionlgn(i,j))*sonlgn(i,j)
      enddo
      do j=1,noflgn(i)
      gl(i) = gl(i) + (1.D0-fnmdatE)*gl2f(ioflgn(i,j))*soflgn(i,j)
     1        + fnmdatE * gl2fn(ioflgn(i,j))*soflgn(i,j)
      enddo

      ginhib  = ((1.D0-fgaba)*gi(i) + fgaba*gj(i))*sei(i) + gEIc
      gexcit  = ((1.D0-fnmdacE)*ge(i) + fnmdacE*gf(i))*see(i)
      ginoise   = (1.D0-fgaba)*g2n(i) + fgaba*g2m(i)
      genoise   = (1.D0-fnmdanE)*g2x(i) + fnmdanE*g2y(i)
      alpha2(i) = -gleak - gl(i) - ginoise - gexcit - 
     1   ginhib - genoise
      beta2(i)  = (gl(i)+gexcit+genoise)*vexcit + 
     1   (ginoise+ginhib)*vinhib + gleak*vreset
      else
      do j=1,nonlgn(i)
      gl(i) = gl(i) + (1.D0-fnmdatI)*gl2o(ionlgn(i,j))*sonlgn(i,j)
     1        + fnmdatI * gl2on(ionlgn(i,j))*sonlgn(i,j)
      enddo
      do j=1,noflgn(i)
      gl(i) = gl(i) + (1.D0-fnmdatI)*gl2f(ioflgn(i,j))*soflgn(i,j)
     1        + fnmdatI * gl2fn(ioflgn(i,j))*soflgn(i,j)
      enddo
      ginhib  = ((1.D0-fgaba)*gi(i) + fgaba*gj(i))*sii(i) + gIIc
      gexcit  = ((1.D0-fnmdacI)*ge(i) + fnmdacI*gf(i))*sie(i)
      ginoise   = (1.D0-fgaba)*g2n(i) + fgaba*g2m(i)
      genoise   = (1.D0-fnmdanI)*g2x(i) + fnmdanI*g2y(i)
      alpha2(i) = -gleakI - gl(i) - ginoise - gexcit - 
     1  ginhib - genoise
      beta2(i)  = (gl(i)+gexcit+genoise)*vexcit + 
     1  (ginoise+ginhib)*vinhib + gleakI*vreset
      endif
      enddo

!===== t+dt ======================
! ampa cortical
      call chain3_d(ge,se1,se2,se3,tau_e,tau0,dt2,nmax)
! nmda cortical 
        if(fnmdacE.GT.0.D0.or.fnmdacI.GT.0.D0) then
      call chain3_d(gf,sf1,sf2,sf3,tnrise,tndamp,dt2,nmax)
        endif

! gabaA cortical
      call chain3_d(gi,si1,si2,si3,tau_i,tau1,dt2,nmax)
! gaba2 cortical
        if (fgaba.GT.0.D0)then
      call chain3(gj,sj1,sj2,sj3,tau2,dt2,nmax)
        endif

      do i=1,nmax
      gl(i)      = 0.D0
      if ( excite(i) ) then
      do j=1,nonlgn(i)
      gl(i) = gl(i) + (1.D0-fnmdatE)*glo(ionlgn(i,j))*sonlgn(i,j)
     1        + fnmdatE * glon(ionlgn(i,j))*sonlgn(i,j)
      enddo
      do j=1,noflgn(i)
      gl(i) = gl(i) + (1.D0-fnmdatE)*glf(ioflgn(i,j))*soflgn(i,j)
     1        + fnmdatE * glfn(ioflgn(i,j))*soflgn(i,j)
      enddo

      ginhib  = ((1.D0-fgaba)*gi(i) + fgaba*gj(i))*sei(i) + gEIc
      gexcit  = ((1.D0-fnmdacE)*ge(i) + fnmdacE*gf(i))*see(i)
      ginoise   = (1.D0-fgaba)*gn(i) + fgaba*gm(i)
      genoise   = (1.D0-fnmdanE)*gx(i) + fnmdanE*gy(i)
      alpha1(i) = -gleak - gl(i) - ginoise - gexcit - 
     1   ginhib - genoise
      beta1(i)  = (gl(i)+gexcit+genoise)*vexcit + 
     1   (ginoise+ginhib)*vinhib + gleak*vreset
      else
      do j=1,nonlgn(i)
      gl(i) = gl(i) + (1.D0-fnmdatI)*glo(ionlgn(i,j))*sonlgn(i,j)
     1        + fnmdatI * glon(ionlgn(i,j))*sonlgn(i,j)
      enddo
      do j=1,noflgn(i)
      gl(i) = gl(i) + (1.D0-fnmdatI)*glf(ioflgn(i,j))*soflgn(i,j)
     1        + fnmdatI * glfn(ioflgn(i,j))*soflgn(i,j)
      enddo
      ginhib  = ((1.D0-fgaba)*gi(i) + fgaba*gj(i))*sii(i) + gIIc
      gexcit  = ((1.D0-fnmdacI)*ge(i) + fnmdacI*gf(i))*sie(i)
      ginoise   = (1.D0-fgaba)*gn(i) + fgaba*gm(i)
      genoise   = (1.D0-fnmdanI)*gx(i) + fnmdanI*gy(i)
      alpha1(i) = -gleakI - gl(i) - ginoise - gexcit - 
     1  ginhib - genoise
      beta1(i)  = (gl(i)+gexcit+genoise)*vexcit + 
     1  (ginoise+ginhib)*vinhib + gleakI*vreset
      endif
      enddo
!=========== rk4 begins ===============
      nspike = 0
      do i=1,nmax 
        a0 = alpha0(i)
        b0 = beta0(i)
        a1 = alpha1(i)
        b1 = beta1(i)
        ahlf = alpha2(i)
        bhlf = beta2(i)

        if ( excite(i) ) then
          tref = trefE
        else
          tref = trefI
        endif
        dtsp = pspike(i)+tref-t
        if (dtsp.lt.dt.and.dtsp.gt.0.D0)then
          dt22 = dt*dt
          ratio = (dtsp/dt)
          ratio2 = ratio*ratio
          ratio1_2 = (1.D0-ratio)**2
          covn = ratio1_2*(1.D0+2.D0*ratio+dtsp*a0)
          covn1 = ratio2*(3.D0-2.D0*ratio+(dtsp-dt)*a1)
          coeff = dtsp*ratio1_2*b0+ratio2*(dtsp-dt)*b1
          a2dt = ahlf*dt
          a2dt2 = a2dt*a2dt
          factor = 2.D0+2.D0*a2dt+a2dt2
          cfk = 1.D0/24.D0*dt*(b0*(2.D0*factor+a1*dt*a2dt2)+
     1          2.D0*(2.D0*b1+bhlf*(8.D0+2.D0*a2dt+a1*dt*(2.D0+a2dt))))
          covn2 = 1.D0/12.D0*((a0+a1)*dt*factor+1.D0/2.D0*dt22*a2dt2
     1              + 12.D0 + 8.D0*a2dt2+2.D0*a2dt2)
          covnratio = covn/covn2
          vnew(i) = (vreset-coeff+covnratio*cfk)/(covnratio+covn1)
        else
          if (dtsp.ge.dt) then
            vnew(i) = vreset
          else
            fk1 = a0*v(i) + b0
            fk2 = ahlf*(v(i)+0.5D0*fk1*dt) + bhlf
            fk3 = ahlf*(v(i)+0.5D0*fk2*dt) + bhlf
            fk4 = a1*(v(i)+fk3*dt) + b1
            vnew(i) = v(i) + dt6*(fk1+2.D0*(fk2+fk3)+fk4)
!$OMP CRITICAL
      if (v(i).lt.vinhib.and.vnew(i).lt.v(i)) then
          print *, 'wtf', iid+1, i
          if (excite(i)) then
             ginhib  = ((1-fgaba)*gi(i) + fgaba*gj(i))*sei(i) + gEIc
             gexcit  = ((1-fnmdacE)*ge(i) + fnmdacE*gf(i))*see(i)
             print *, npEI*isptot/(nix*niy),' Inh spikes -> 1 Exc'
             print *, npEE*nsptot/(nex*ney),' Exc spikes -> 1 Exc'
          else
             ginhib  = ((1-fgaba)*gi(i) + fgaba*gj(i))*sii(i) + gIIc
             gexcit  = ((1-fnmdacI)*ge(i) + fnmdacI*gf(i))*sie(i)
          print *, npEI*isptot/(nix*niy),' Inh spikes',npEI,' -> 1 Inh'
          print *, npEE*nsptot/(nex*ney),' Exc spikes',npEE,' -> 1 Inh'
          endif
          print *, gl(i),gexcit,ginhib
          print *, v(i), vnew(i)
          print *, a0,b0
          print *, ahlf,bhlf
          print *, a1,b1
          print *, fk1,fk2,fk3,fk4
          stop
      endif
!$OMP END CRITICAL
          endif
        endif
        if (vnew(i).gt.vthres) then
        
        v0 = v(i) 
        v1 = vnew(i)
        vt0 = a0*v0+b0
        vt1 = a1*v1+b1
        dtsp = dt*(vthres-v0)/(v1-v0)
        vvv = herm_int(dtsp,dt,v0,v1,vt0,vt1)
        vvv0 = vvv
        vvt = hermt_int(dtsp,dt,v0,v1,vt0,vt1)
        vvt0 = vvt
        dtsp0 = dtsp
        slope = -1.D0
        iflag = 0
        iSlopeFlag = 0
        diff = vvv/vvt
        if (dabs(vvv/vvt).le.1.D-08) then
            slope = vvt
        else    
          do it = 1,itmax
            dtsp = dtsp - diff
            vvv = herm_int(dtsp,dt,v0,v1,vt0,vt1)
            vvt = hermt_int(dtsp,dt,v0,v1,vt0,vt1)
            diff = (vvv-vthres)/vvt
            if(dabs(diff/dt).le.1.D-08.and.
     1            dtsp.ge.0.D0.and.dtsp.le.dt)then
              slope = vvt 
              iSlopeFlag = 1
              if (slope.GT.0.D0) then
                iflag = 1
              endif
              exit
            endif
          enddo
        endif
        if (iflag.EQ.0)then
          if (vt1.GT.0.D0.and.vt0.GT.0.D0.and.iSlopeFlag.eq.1) then
!======= 3-root case: find from 0 ===
            dtsp = 0.D0
            vvv = herm_int(dtsp,dt,v0,v1,vt0,vt1)
            vvt = hermt_int(dtsp,dt,v0,v1,vt0,vt1)
            dtsp0 = dtsp
            vvv0 = vvv
            vvt0 = vvt
            call newton(iflag,dtsp,vvv,vvt,dt,v0,v1,vt0,vt1,itmax)
            if (iflag.EQ.0) then
!$OMP CRITICAL
              print *,'3-root case failed, increase itmax'
              print *, 'v0', v0, 'v1', v1
              print *, 'ge',ge(i),'gi',gi(i),'gl',gl(i)
              print *, 'gen',gx(i),'gin',gn(i)
              print *, 'before',dtsp0,vvv0, vvt0
              print *, 'after',dtsp,vvv,vvt
              stop
!$OMP END CRITICAL
            endif
          else
!======= find tspike in 1-root cases using bisect_newton =====
            vvv = vvv0
            vvt = vvt0
            dtsp = dtsp0
            if (vvv0.gt.0.D0) then
              call bisect_newton(iflag,dtsp,dt,v0,v1,vt0,vt1,
     1                            0.D0,dtsp,vvv,vvt,itmax)
            else
              call bisect_newton(iflag,dtsp,dt,v0,v1,vt0,vt1,
     1                            dtsp,dt,vvv,vvt,itmax)
            endif
            if (iflag.EQ.0)then
!$OMP CRITICAL
              print *,'1-root case failed'
              print *, 'v0', v0, 'v1', v1
              print *, 'ge',ge(i),'gi',gi(i),'gl',gl(i)
              print *, 'gen',gx(i),'gin',gn(i)
              print *, 'before',dtsp0,vvv0, vvt0
              print *, 'after',dtsp,vvv,vvt
              print *, nsptot,' Exc spikes'
              stop
!$OMP END CRITICAL
            endif
          endif
        endif
!======= finish tspike estimation =====
        pspike(i) = t+dtsp
        if (pspike(i) + tref.lt. t+dt) then
          dtspTref = dtsp+tref
          dt22 = dt*dt
          ratio = (dtspTref/dt)
          ratio2 = ratio*ratio
          ratio1_2 = (1.D0-ratio)**2
          covn = ratio1_2*(1.D0+2.D0*ratio+dtspTref*a0)
          covn1 = ratio2*(3.D0-2.D0*ratio+(dtspTref-dt)*a1)
          coeff = dtspTref*ratio1_2*b0+ratio2*(dtspTref-dt)*b1
          a2dt = ahlf*dt
          a2dt2 = a2dt*a2dt
          factor = 2.D0+2.D0*a2dt+a2dt2
          cfk = 1.D0/24.D0*dt*(b0*(2.D0*factor+a1*dt*a2dt2)+
     1          2.D0*(2.D0*b1+bhlf*(8.D0+2.D0*a2dt+a1*dt*(2.D0+a2dt))))
          covn2 = 1.D0/12.D0*((a0+a1)*dt*factor+1.D0/2.D0*dt22*a2dt2
     1              + 12.D0 + 8.D0*a2dt2+2.D0*a2dt2)
          covnratio = covn/covn2
          vnew(i) = (vreset-coeff+covnratio*cfk)/(covnratio+covn1)
        else
          vnew(i) = vreset
        endif
        if (vnew(i).gt.vthres) print *,'2 spikes per dt! i:',
     1           i,t+dtsp,'v0:',v(i),'vnew:',vnew(i)
        nspike = nspike + 1
        tspike(nspike) = dtsp
        ispike(nspike) = i

        if ( .not. excite(i) ) then
            isptot = isptot + 1
        else 
            erate(i) = erate(i) + 1.D0
            nsptot = nsptot + 1
        endif
        if (nlgni(i).EQ.lgnmin.and.excite(i)) jsptot = jsptot + 1
        if (nlgni(i).GE.lgnmax.and.excite(i)) ksptot = ksptot + 1
        if (nlgni(i).EQ.meanlgn.and.excite(i)) msptot = msptot + 1
!------------------------------------------------------------
!  Construct histogram for spike rate: ith neuron nspikes at nirate th bin 
!------------------------------------------------------------
        if (t.ge.rstart) then
          tirate = dmod(t+dtsp,period)
          nirate = floor(tirate*1.D0*iperiod/period) + 1
          irate(i,nirate) = irate(i,nirate) + 1
        endif
!---------------------------------------neuron firing end
        endif
!      if (vnew(i).lt.-10.D0) then
!!$OMP CRITICAL
!        print *, i, 'th neuron'
!        print *, v(i), vnew(i)
!        print *, '          a0,              a1'
!        print *, sngl(a0),sngl(a1)
!        print *, '          b0,              b1'
!        print *, sngl(b0),sngl(b1)
!        print *, '          fk1,            fk2'
!        print *, sngl(fk1),sngl(fk2)
!        print *, '            t,             dt'
!        print *, sngl(t),sngl(dtsp)
!        stop
!!$OMP END CRITICAL
!      endif
        enddo
        return
      end
!************************************************************
      subroutine visual4(t,dt,iseed)
      use parameters
!------------------------------------------------------------
!     Also generate noise (1-1) given firing rates
!
!     Generates LGN spike times using Poisson process 
!    with time-dependent firing rate a Aunction of
!    visual stimulus
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension glo(nlgn), slo1(nlgn), slo2(nlgn), slo3(nlgn)
      dimension glon(nlgn),slon1(nlgn),slon2(nlgn),slon3(nlgn)
      dimension glf(nlgn), slf1(nlgn), slf2(nlgn), slf3(nlgn)
      dimension glfn(nlgn),slfn1(nlgn),slfn2(nlgn),slfn3(nlgn)
      dimension gx(nmax), sx1(nmax), sx2(nmax), sx3(nmax)
      dimension gy(nmax), sy1(nmax), sy2(nmax), sy3(nmax)
      dimension gn(nmax), sn1(nmax), sn2(nmax), sn3(nmax)
      dimension gm(nmax), sm1(nmax), sm2(nmax), sm3(nmax)
      dimension gl2o(nlgn), sl2o1(nlgn), sl2o2(nlgn), sl2o3(nlgn)
      dimension gl2on(nlgn),sl2on1(nlgn),sl2on2(nlgn),sl2on3(nlgn)
      dimension gl2f(nlgn), sl2f1(nlgn), sl2f2(nlgn), sl2f3(nlgn)
      dimension gl2fn(nlgn),sl2fn1(nlgn),sl2fn2(nlgn),sl2fn3(nlgn)
      dimension g2x(nmax), s2x1(nmax), s2x2(nmax), s2x3(nmax)
      dimension g2y(nmax), s2y1(nmax), s2y2(nmax), s2y3(nmax)
      dimension g2n(nmax), s2n1(nmax), s2n2(nmax), s2n3(nmax)
      dimension g2m(nmax), s2m1(nmax), s2m2(nmax), s2m3(nmax)
      dimension pspon(nlgn),pspoff(nlgn)
      dimension xlgn(nlgn),ylgn(nlgn)
      dimension tsonpoi(nlgn),ionpoi(nlgn)
      dimension tsoffpoi(nlgn),ioffpoi(nlgn)
      dimension tsenoi(nmax),ienoi(nmax)
      dimension tsinoi(nmax),iinoi(nmax)
      dimension nlgni(nmax)
      logical excite(nmax)
      integer iseed
      data iword / 8 /
      common / chaino / glo, slo1, slo2, slo3
      common / chaino2/ glon,slon1,slon2,slon3
      common / chainf / glf, slf1, slf2, slf3
      common / chainf2/ glfn,slfn1,slfn2,slfn3
      common / chain2o / gl2o, sl2o1, sl2o2, sl2o3
      common / chain2o2/ gl2on,sl2on1,sl2on2,sl2on3
      common / chain2f / gl2f, sl2f1, sl2f2, sl2f3
      common / chain2f2/ gl2fn,sl2fn1,sl2fn2,sl2fn3
      common / chainx / gx,sx1,sx2,sx3
      common / chainy / gy,sy1,sy2,sy3
      common / chainn / gn,sn1,sn2,sn3
      common / chainm / gm,sm1,sm2,sm3
      common / chain2x / g2x,s2x1,s2x2,s2x3
      common / chain2y / g2y,s2y1,s2y2,s2y3
      common / chain2n / g2n,s2n1,s2n2,s2n3
      common / chain2m / g2m,s2m1,s2m2,s2m3
      common / lgnIndivid / pspon,pspoff,gkx,gky,eps,satur,gtheta
      common / lgnGain / tmpphi,conAmp,g0,frtlgn0,treflgn,
     1   gphi, tstart, gfail,cond0
      common / lgnpos / xlgn,ylgn,jlgn
      common / noises / ce0_e,ci0_e,frtinh_e,frtexc_e,
     1   ce0_i,ci0_i,frtinh_i,frtexc_i
      common / tconst / tau_e,tau_i,tau0,tau1,tau2,tnrise,tndamp,
     1 tau_a0, tau_a1
      common /  NMDA  / fnmdatE,fnmdacE,fnmdatI,fnmdacI,
     1   fnmdanE, fnmdanI, fgaba
      common / neuron / excite,nlgni,trefE,trefI,nr,
     1   lgnmax,lgnmin,meanlgn
      common / lgnConvol / gk2, taulgn, omega
      common / ttotal / tfinal,twindow,rstart
      common / onpoi / tsonpoi,ionpoi,nonpoi
      common / offpoi / tsoffpoi,ioffpoi,noffpoi
      common / enoi / tsenoi,ienoi,nenoi
      common / inoi / tsinoi,iinoi,ninoi
      real*8 lgnsatur
      !INTEGER*4 nonsp,noffsp,nensp,ninsp
!$OMP THREADPRIVATE(
!$OMP& /chaino/, /chaino2/, /chainf/, /chainf2/,
!$OMP& /chainx/, /chainy/,  /chainn/, /chainm/,
!$OMP& /chain2o/,/chain2o2/,/chain2f/,/chain2f2/,
!$OMP& /chain2x/,/chain2y/, /chain2n/,/chain2m/,
!$OMP& /lgnIndivid/,/onpoi/,/offpoi/,/enoi/,/inoi/)
      iword2 = iword/2
      omegat = omega*(t-tstart)
      dt2 = dt/2.D0

      nonpoi = 0
      noffpoi = 0
      ninoi = 0
      nenoi = 0

      do i=1,jlgn
          gl2o(i)  = glo(i)
          sl2o1(i) = slo1(i)
          sl2o2(i) = slo2(i)
          sl2o3(i) = slo3(i)
          gl2on(i) = glon(i)
          sl2on1(i)= slon1(i)
          sl2on2(i)= slon2(i)
          sl2on3(i)= slon3(i)
       if (t.lt.0.D0) frate = frtlgn0
       rannum = -1.D0
       if (t+dt.gt.pspon(i)+treflgn) then
!-------------Find firing rate as function of visual stimulus
! For CONTRAST REVERSAL, need only amplitude of sinusoid
!============================================================
      if ((t.ge.0.D0).and.(t.lt.tstart)) then ! slowly ramp up
           !frate = frtlgn0 + (t/tstart)**2*eps*
           !frate = frtlgn0 + t/tstart*eps*
           frate = frtlgn0 + eps*
     1        conAmp*sin(omegat+tmpphi-gphi+xlgn(i)*gkx+ylgn(i)*gky)
      endif
      if (t.ge.tstart) then
           frate = frtlgn0 + eps*
     1        conAmp*sin(omegat+tmpphi-gphi+xlgn(i)*gkx+ylgn(i)*gky)
      endif
!------------------------Compute spikes only if non-zero rate
      if (frate.gt.0D0) then
!------------------------------------------------------------
! Contrast Saturation (200 levels, linear interpolate)
!------------------------------------------------------------
      satur =  lgnsatur(frate)
      frate = gfail*satur
!------------------------------------------------------------
      rannum = ran2(iseed)
      endif
       endif
!===========================================================
       dtt = rannum/frate
       if (rannum.lt.dt*frate.and.rannum.ge.0.D0
     1      .and.t+dtt.gt.pspon(i)+treflgn) then
!===========================================================
       pspon(i) = t + dtt
!=============================================================
      if (t.gt.tfinal-twindow) then
        nonpoi = nonpoi + 1
        tsonpoi(nonpoi) = t + dtt
        ionpoi(nonpoi) = i
      endif
!------------------------------------
       !lgnspikeo(i) = lgnspikeo(i)+1
!--------- t + dt -- dt

!--------- before spike | 0 --*--> dtt -----> dt
!--------- ampa exc
      call decay3_d(glo(i),slo1(i),slo2(i),slo3(i),tau_e,tau0,dtt)
      slo3(i) = slo3(i) + cond0
!--------- nmda exc 
      if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
        call decay3_d(glon(i),slon1(i),slon2(i),
     1                 slon3(i),tnrise,tndamp,dtt)
        slon3(i) = slon3(i) + cond0*tau_e/tnrise
      endif 
!-----------rk4 conductantce at dt2-----------
!---------------------------------------
      if (dtt.lt.dt2) then
        gl2o(i)  = glo(i)
        sl2o1(i) = slo1(i)
        sl2o2(i) = slo2(i)
        sl2o3(i) = slo3(i)
        gl2on(i) = glon(i)
        sl2on1(i)= slon1(i)
        sl2on2(i)= slon2(i)
        sl2on3(i)= slon3(i)

        dtt = dt2 - dtt
        call decay3_d(gl2o(i),sl2o1(i),sl2o2(i),sl2o3(i),tau_e,tau0,dtt)
        if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
          call decay3_d(gl2on(i),sl2on1(i),sl2on2(i),
     1                   sl2on3(i),tnrise,tndamp,dtt)
        endif
        ! reset dtt for full-length conductance at dt
        dtt = dt2 - dtt
      else
        call decay3_d(gl2o(i),sl2o1(i),sl2o2(i),sl2o3(i),tau_e,tau0,dt2)
        if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
          call decay3_d(gl2on(i),sl2on1(i),sl2on2(i),
     1                   sl2on3(i),tnrise,tndamp,dt2)
        endif
      endif
!---------------------------------------
!--------- after spike | 0 -----> dtt --*--> dt
      dtt = dt - dtt
      call decay3_d(glo(i),slo1(i),slo2(i),slo3(i),tau_e,tau0,dtt)
      if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
        call decay3_d(glon(i),slon1(i),slon2(i),
     1                 slon3(i),tnrise,tndamp,dtt)
      endif
!--------- if no spike just decay--------------------------------
      else
        call decay3_d(glo(i),slo1(i),slo2(i),slo3(i),tau_e,tau0,dt)
        if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
          call decay3_d(glon(i),slon1(i),slon2(i),
     1                   slon3(i),tnrise,tndamp,dt)
        endif
        call decay3_d(gl2o(i),sl2o1(i),sl2o2(i),sl2o3(i),tau_e,tau0,dt2)
        if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
          call decay3_d(gl2on(i),sl2on1(i),sl2on2(i),
     1                   sl2on3(i),tnrise,tndamp,dt2)
        endif
!===========================================================
       endif
!===========================================================
      enddo
!----------------------------------Now off-centered LGN cells

      do i=1,jlgn
          gl2f(i)  = glf(i)
          gl2fn(i) = glfn(i)
          sl2f1(i) = slf1(i)
          sl2f2(i) = slf2(i)
          sl2f3(i) = slf3(i)
          sl2fn1(i)= slfn1(i)
          sl2fn2(i)= slfn2(i)
          sl2fn3(i)= slfn3(i)
        if (t.lt.0.D0) frate = frtlgn0
        rannum = -1.D0
!===========================================================
      if (t+dt.gt.pspoff(i)+treflgn) then
!===========================================================
      if ((t.ge.0.D0).and.(t.lt.tstart)) then
         !frate = frtlgn0 - (t/tstart)**2*eps
         frate = frtlgn0 - t/tstart*eps*
     1   conAmp*sin(omegat+tmpphi-gphi+xlgn(i)*gkx+ylgn(i)*gky)
      endif
      if (t.ge.tstart) then
         frate = frtlgn0 - eps*
     1   conAmp*sin(omegat+tmpphi-gphi+xlgn(i)*gkx+ylgn(i)*gky)
      endif
!------------------------------------------------------------
        if (frate.gt.0.D0) then
!------------------------------------------------------------
! Contrast Saturation
!------------------------------------------------------------
      satur = lgnsatur(frate)
      frate = gfail*satur
      rannum = ran2(iseed)
        endif    
      endif
!===========================================================
      dtt = rannum/frate
      if (rannum.lt.dt*frate.and.rannum.ge.0.D0
     1  .and.t+dtt.gt.pspoff(i)+treflgn) then
!===========================================================
      !lgnspikef(i) = lgnspikef(i)+1
      pspoff(i) = t + dtt
!=============================================================
      if (t.gt.tfinal-twindow) then
        noffpoi = noffpoi + 1
        tsoffpoi(noffpoi) = t+dtt
        ioffpoi(noffpoi) = i
      endif
!------------------------------------
      call decay3_d(glf(i),slf1(i),slf2(i),slf3(i),tau_e,tau0,dtt)
      slf3(i) = slf3(i) + cond0
      if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
          call decay3_d(glfn(i),slfn1(i),slfn2(i),
     1                  slfn3(i),tnrise,tndamp,dtt)
          slfn3(i) = slfn3(i) + cond0/tnrise*tau_e
      endif
      if (dtt.lt.dt2) then
        gl2f(i)  = glf(i)
        gl2fn(i) = glfn(i)
        sl2f1(i) = slf1(i)
        sl2f2(i) = slf2(i)
        sl2f3(i) = slf3(i)
        sl2fn1(i)= slfn1(i)
        sl2fn2(i)= slfn2(i)
        sl2fn3(i)= slfn3(i)

        dtt = dt2 - dtt
        call decay3_d(gl2f(i),sl2f1(i),sl2f2(i),sl2f3(i),tau_e,tau0,dtt)
        if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
          call decay3_d(gl2fn(i),sl2fn1(i),sl2fn2(i),
     1                  sl2fn3(i),tnrise,tndamp,dtt)
        endif
        dtt = dt2 - dtt
      else
        call decay3_d(gl2f(i),sl2f1(i),sl2f2(i),sl2f3(i),tau_e,tau0,dt2)
        if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
          call decay3_d(gl2fn(i),sl2fn1(i),sl2fn2(i),
     1                  sl2fn3(i),tnrise,tndamp,dt2)
        endif
      endif

      dtt = dt - dtt
      call decay3_d(glf(i),slf1(i),slf2(i),slf3(i),tau_e,tau0,dtt)
      if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
        call decay3_d(glfn(i),slfn1(i),slfn2(i),
     1                  slfn3(i),tnrise,tndamp,dtt)
      endif
!--------- if no spike just decay--------------------------------
      else
        call decay3_d(glf(i),slf1(i),slf2(i),slf3(i),tau_e,tau0,dt)
        if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
          call decay3_d(glfn(i),slfn1(i),slfn2(i),
     1                  slfn3(i),tnrise,tndamp,dt)
        endif
        call decay3_d(gl2f(i),sl2f1(i),sl2f2(i),sl2f3(i),tau_e,tau0,dt2)
        if (fnmdatE.GT.0.D0.or.fnmdatI.GT.0.D0) then
          call decay3_d(gl2fn(i),sl2fn1(i),sl2fn2(i),
     1                  sl2fn3(i),tnrise,tndamp,dt2)
        endif
      endif
      enddo
!======================================== Exc noise units
      if (frtexc_e.gt.0.D0.and.ce0_e.gt.0.D0
     1  .or.frtexc_i.gt.0.D0.and.ce0_i.gt.0.D0)then
      do i=1,nmax
        g2x(i) = gx(i)
        s2x1(i)= sx1(i)
        s2x2(i)= sx2(i)
        s2x3(i)= sx3(i)
        g2y(i) = gy(i)
        s2y1(i)= sy1(i)
        s2y2(i)= sy2(i)
        s2y3(i)= sy3(i)
        if (excite(i)) then
            frtexc = frtexc_e
            ce0 = ce0_e
        else
            frtexc = frtexc_i
            ce0 = ce0_i
        endif
        rannum = ran2(iseed)
      if (rannum.lt.dt*frtexc.and.rannum.ge.0.D0) then
        dtt = rannum/frtexc
!------------------------------------------------------------
!  kick each input layer cell at tspike w/ delta fcn of 
!     strength Ce, Ci
!=============================================================
        if (t.gt.tfinal-twindow) then
          nenoi = nenoi + 1
          tsenoi(nenoi) = t + dtt
          ienoi(nenoi) = i
        endif
!------------------------------------
        call decay3_d(gx(i),sx1(i),sx2(i),sx3(i),tau_e,tau0,dtt)
        sx3(i) = sx3(i) + ce0/tau_e
        if (fnmdanE.GT.0.D0.or.fnmdanI.GT.0.D0) then
          call decay3_d(gy(i),sy1(i),sy2(i),sy3(i),tnrise,tndamp,dtt)
          sy3(i) = sy3(i) + ce0/tnrise
        endif
        if (dtt.lt.dt2) then
          g2x(i) = gx(i)
          s2x1(i)= sx1(i)
          s2x2(i)= sx2(i)
          s2x3(i)= sx3(i)
          g2y(i) = gy(i)
          s2y1(i)= sy1(i)
          s2y2(i)= sy2(i)
          s2y3(i)= sy3(i)
          
          dtt = dt2 - dtt
          call decay3_d(g2x(i),s2x1(i),s2x2(i),s2x3(i),tau_e,tau0,dtt)
          if (fnmdanE.GT.0.D0.or.fnmdanI.GT.0.D0) then
            call decay3_d(g2y(i),s2y1(i),s2y2(i),
     1                  s2y3(i),tnrise,tndamp,dtt)
          endif
          dtt = dt2 - dtt
        else
          call decay3_d(g2x(i),s2x1(i),s2x2(i),s2x3(i),tau_e,tau0,dt2)
          if (fnmdanE.GT.0.D0.or.fnmdanI.GT.0.D0) then
            call decay3_d(g2y(i),s2y1(i),s2y2(i),
     1                  s2y3(i),tnrise,tndamp,dt2)
          endif
        endif
        dtt = dt - dtt
        call decay3_d(gx(i),sx1(i),sx2(i),sx3(i),tau_e,tau0,dtt)
        if (fnmdanE.GT.0.D0.or.fnmdanI.GT.0.D0) then
          call decay3_d(gy(i),sy1(i),sy2(i),sy3(i),tnrise,tndamp,dtt)
        endif
      else
        call decay3_d(gx(i),sx1(i),sx2(i),sx3(i),tau_e,tau0,dt)
        if (fnmdanE.GT.0.D0.or.fnmdanI.GT.0.D0) then
          call decay3_d(gy(i),sy1(i),sy2(i),sy3(i),tnrise,tndamp,dt)
        endif
        call decay3_d(g2x(i),s2x1(i),s2x2(i),s2x3(i),tau_e,tau0,dt2)
        if (fnmdanE.GT.0.D0.or.fnmdanI.GT.0.D0) then
          call decay3_d(g2y(i),s2y1(i),s2y2(i),
     1                  s2y3(i),tnrise,tndamp,dt2)
        endif
      endif
      enddo
      endif
!======================================= Inh noise units
      if (frtinh_e.gt.0.D0.and.ci0_e.gt.0.D0
     1 .or.frtinh_i.gt.0.D0.and.ci0_i.gt.0.D0)then
      do i=1,nmax
         g2n(i) = gn(i)
         s2n1(i)= sn1(i)
         s2n2(i)= sn2(i)
         s2n3(i)= sn3(i)
         g2m(i) = gm(i)
         s2m1(i)= sm1(i)
         s2m2(i)= sm2(i)
         s2m3(i)= sm3(i)
         if (excite(i)) then
             frtinh = frtinh_e
             ci0 = ci0_e
         else
             frtinh = frtinh_i
             ci0 = ci0_i
         endif
         rannum = ran2(iseed)
      if (rannum.lt.dt*frtinh.and.rannum.ge.0.D0) then
         dtt = rannum/frtinh
!=============================================================
        if (t.gt.tfinal-twindow) then
          ninoi = ninoi + 1
          tsinoi(ninoi) = t + dtt
          iinoi(ninoi) = i
        endif
!------------------------------------
        call decay3_d(gn(i),sn1(i),sn2(i),sn3(i),tau_i,tau1,dtt)
        sn3(i) = sn3(i) + ci0/tau_i
        if (fgaba.GT.0.D0) then
          call decay3(gm(i),sm1(i),sm2(i),sm3(i),tau2,dtt)
          sm3(i) = sm3(i) + ci0/tau2
        endif
        
        if (dtt.lt.dt2) then
          g2n(i) = gn(i)
          s2n1(i)= sn1(i)
          s2n2(i)= sn2(i)
          s2n3(i)= sn3(i)
          g2m(i) = gm(i)
          s2m1(i)= sm1(i)
          s2m2(i)= sm2(i)
          s2m3(i)= sm3(i)

          dtt = dt2 - dtt
          call decay3_d(g2n(i),s2n1(i),s2n2(i),s2n3(i),tau_i,tau1,dtt)
          if (fgaba.GT.0.D0) then
            call decay3(g2m(i),s2m1(i),s2m2(i),s2m3(i),tau2,dtt)
          endif
          dtt = dt2 - dtt
        else
          call decay3_d(g2n(i),s2n1(i),s2n2(i),s2n3(i),tau_i,tau1,dt2)
          if (fgaba.GT.0.D0) then
            call decay3(g2m(i),s2m1(i),s2m2(i),s2m3(i),tau2,dt2)
          endif
        endif

        dtt = dt - dtt 
        call decay3_d(gn(i),sn1(i),sn2(i),sn3(i),tau_i,tau1,dtt)
        if (fgaba.GT.0.D0) then
          call decay3(gm(i),sm1(i),sm2(i),sm3(i),tau2,dtt)
        endif
      else
         call decay3_d(gn(i),sn1(i),sn2(i),sn3(i),tau_i,tau1,dt)
         if (fgaba.GT.0.D0) then
           call decay3(gm(i),sm1(i),sm2(i),sm3(i),tau2,dt)
         endif
         call decay3_d(g2n(i),s2n1(i),s2n2(i),s2n3(i),tau_i,tau1,dt2)
         if (fgaba.GT.0.D0) then
           call decay3(g2m(i),s2m1(i),s2m2(i),s2m3(i),tau2,dt2)
         endif
      endif
      enddo
      endif
!------------------------------------------------------------
      return
      end
!************************************************************
      subroutine readConnection(fcm,fc,excite)
      use parameters
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension NpostSynE(nmax), NpostSynI(nmax)
      dimension NpreSynE(nmax), NpreSynI(nmax)
      dimension IpostSyn(nmax,nmax)
      dimension Istrength(nmax,nmax)
      integer*1 IcMatrix(nmax,nmax) ! (i,j) i->j
      logical excite(nmax)
      common / post    / NpostSynE, NpostSynI, IpostSyn,
     1    NpreSynE, NpreSynI, Istrength
      common / poststatics/ nEE, nEI, nIE, nII
      common / prestatics/ npEE, npEI, npIE, npII
      character (LEN=20) fc,fcm

      do i=1,nmax
          NpostSynE(i) = 0
          NpostSynI(i) = 0
          NpreSynE(i) = 0
          NpreSynI(i) = 0
      enddo

      open(16,file=fcm,status='old',form='unformatted',
     1       access='direct',recl=nmax)
      do i=1,nmax
        read(16,rec=i) (IcMatrix(j,i),j=1,nmax)
      enddo
      close(16)
      
      do i=1,nmax
        do j=1,nex*ney
            if (IcMatrix(i,j).GT.0) then
                if (excite(i)) then
                    NpreSynE(j) = NpreSynE(j)+1
                else
                    NpreSynI(j) = NpreSynI(j)+1
                endif
                NpostSynE(i) = NpostSynE(i)+1
                IpostSyn(i,NpostSynE(i)+NpostSynI(i)) = j 
        Istrength(i,NpostSynE(i)+NpostSynI(i))=IcMatrix(i,j)
            endif
        enddo
        do j=nex*ney+1,nmax
            if (IcMatrix(i,j).GT.0) then
                if (excite(i)) then
                    NpreSynE(j) = NpreSynE(j)+1
                else
                    NpreSynI(j) = NpreSynI(j)+1
                endif
                NpostSynI(i) = NpostSynI(i)+1
                IpostSyn(i,NpostSynE(i)+NpostSynI(i)) = j 
        Istrength(i,NpostSynE(i)+NpostSynI(i))=IcMatrix(i,j)
            endif
        enddo
      enddo
      npEE = 0 
      npEI = 0 
      npIE = 0 
      npII = 0 
      nmaxpEE = 0
      nmaxpEI = 0
      nminpEE = nmax
      nminpEI = nmax
      do i=1,nmax
        if (excite(i)) then
            npEE = npEE + NpreSynE(i)
            npEI = npEI + NpreSynI(i)
            if(NpreSynE(i).gt.nmaxpEE)then
                nmaxpEE = NpreSynE(i)
            endif
            if(NpreSynI(i).gt.nmaxpEI)then
                nmaxpEI = NpreSynI(i)
            endif
            if(NpreSynE(i).lt.nminpEE)then
                nminpEE = NpreSynE(i)
            endif
            if(NpreSynI(i).lt.nminpEI)then
                nminpEI = NpreSynI(i)
            endif
        else
            npIE = npIE + NpreSynE(i)
            npII = npII + NpreSynI(i)
        endif
      enddo
      print *,'nPreEE',nminpEE,'=>',nmaxpEE
      print *,'nPreEI',nminpEI,'=>',nmaxpEI
      nEE = 0 
      nEI = 0 
      nIE = 0 
      nII = 0 
      do i=1,nmax
        if (excite(i)) then
            nEE = nEE + NpostSynE(i)
            nIE = nIE + NpostSynI(i)
        else
            nEI = nEI + NpostSynE(i)
            nII = nII + NpostSynI(i)
        endif
      enddo

      open(16,file=fc,status='replace',form='unformatted',
     1       access='direct',recl=4*nmax)
      write(16,rec=1) (NpostSynE(i),i=1,nmax)
      close(16)
      open(16,file=fc,status='old',form='unformatted',
     1       access='direct',recl=4*nmax)
      write(16,rec=2) (NpostSynI(i),i=1,nmax)
      close(16)
      open(16,file=fc,status='old',form='unformatted',
     1       access='direct',recl=4*nmax)
      write(16,rec=3) (NpreSynE(i),i=1,nmax)
      close(16)
      open(16,file=fc,status='old',form='unformatted',
     1       access='direct',recl=4*nmax)
      write(16,rec=4) (NpreSynI(i),i=1,nmax)
      close(16)
    
      open(16,file=fc,status='old',form='unformatted',
     1       access='direct',recl=nmax)
      do i=17,nmax+16
        write(16,rec=i) (IcMatrix(i-16,j),j=1,nmax)
      enddo
      close(16)
      end
      function thetaGauss(dtheta,sigo,A)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
        thetaGauss = A*exp(-(dtheta/sigo)**2/2.D0)
      end
      REAL*8 function prob(A,x2,da2)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
        prob = A*exp(-0.5D0*x2/da2)
        return 
      end
      subroutine connection(dE,dI,aee,aei,aie,aii,thetas,sigo,
     1  retroInh,nlgni,iseed,nr,excite,fc,lgnmax)
      use parameters
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      dimension pex(ney,nex), pey(ney,nex),pix(niy,nix),piy(niy,nix)
      dimension NpostSynE(nmax), NpostSynI(nmax),nlgni(nmax)
      dimension NpreSynE(nmax), NpreSynI(nmax)
      dimension IpostSyn(nmax,nmax)
      integer*1 IcMatrix(nmax,nmax) ! (i,j) i->j
      dimension Istrength(nmax,nmax)
      dimension thetas(nmax)
      logical excite(nmax)
      CHARACTER (LEN=20) fc
      common / connec / synfailEE,synfailEI,synfailIE,synfailII,
     1  neffEE,neffEI,neffIE,neffII
      common / post    / NpostSynE, NpostSynI, IpostSyn,
     1    NpreSynE, NpreSynI, Istrength
      common / poststatics/ nEE, nEI, nIE, nII
      common / prestatics/ npEE, npEI, npIE, npII
      common / npos / pex, pey, pix, piy
      neff = (neffEE+neffEI+neffII+neffIE)/4.D0
      if (neff.GT.0.and.neff.LT.nmax)then
            pi = 4.D0*atan(1.D0)
            sigo = sigo/180.D0*pi
            A = pi/(erf(pi/(2.D0*sqrt(2.D0)*sigo))*
     1    sqrt(2.D0*pi)*sigo)
      de2 = dE**2
      dI2 = dI**2
      do i=1,nmax
          if (thetas(i).GT.pi) thetas(i) = thetas(i) - pi
          NpostSynE(i) = 0
          NpostSynI(i) = 0
          NpreSynE(i) = 0
          NpreSynI(i) = 0
          do j=1,nmax
              IcMatrix(i,j) = 0
          enddo
      enddo

      do ii=1,ney*nex
        iy = mod(ii-1,ney) + 1
        ix = (ii-1)/ney+1
!====================================
!        print *,'E <- E presynaptic',ix,iy
!===================================================
        rlx = (nex-1.D0)/2.D0 * dE
          rly = (ney-1.D0)/2.D0 * dE
        do j=1,ney
            do k=1,nex
            jj = ney*(k-1)+j
            if (jj.NE.ii) then
            dx = pex(iy,ix)-pex(j,k)
            dy = pey(iy,ix)-pey(j,k)
        if (dabs(dx).gt. rlx) then
            dx = 2*mod(pex(iy,ix)+rlx,2*rlx)-pex(j,k)
     1           - pex(iy,ix)
        endif
        if (dabs(dy).gt. rly) then
            dy = 2.D0*mod(pey(iy,ix)+rly,2.D0*rly)-pey(j,k)
     1           - pey(iy,ix)
        endif
            distance2 = dx**2 + dy**2
              prefix = de2*neffEE/(2.D0*pi*aee)
              aa = prob(prefix,distance2,aee)
        if (nlgni(jj).GT.nr) then        !    oriented
              dtheta = dabs(thetas(ii) - thetas(jj))
              if (dtheta.GT.pi/2.D0) dtheta = pi - dtheta
              aa = aa*thetaGauss(dtheta,sigo,A)
        endif
!====================================================
              if (aa.GT.1.D0.or.ran2(iseed).LE.aa) then
            NpostSynE(jj) = NpostSynE(jj)+1
            NpreSynE(ii) = NpreSynE(ii)+1
            IpostSyn(jj,NpostSynE(jj)+NpostSynI(jj))=ii
            IcMatrix(jj,ii) = 1
              endif
            endif
            enddo
        enddo
!===================================
!            print *,'E <- I presynaptic'
!===================================================
        rlx = (nix-1.D0)/2.D0 * dI
        rly = (niy-1.D0)/2.D0 * dI
        do j=1,niy
            do k=1,nix
            jj = ney*nex + niy*(k-1)+j
            dx = pex(iy,ix)-pix(j,k)
            dy = pey(iy,ix)-piy(j,k)
        if (dabs(dx).gt. rlx) then
            dx = 2.D0*mod(pex(iy,ix)+rlx,2.D0*rlx)-pix(j,k)
     1           - pix(iy,ix)
        endif
        if (dabs(dy).gt. rly) then
            dy = 2.D0*mod(pey(iy,ix)+rly,2.D0*rly)-piy(j,k)
     1           - piy(iy,ix)
        endif
            distance2 = dx**2 + dy**2
              prefix = di2*neffEI/(2.D0*pi*aei)
              aa = prob(prefix,distance2,aei)
              if (aa.GT.1.D0.or.ran2(iseed).LE.aa) then
            NpostSynE(jj) = NpostSynE(jj)+1
            NpreSynI(ii) = NpreSynI(ii)+1
            IpostSyn(jj,NpostSynE(jj)+NpostSynI(jj))=ii
            IcMatrix(jj,ii) = 1
            endif
!===================================================
            enddo
        enddo
      enddo
!===================================================
!            Inh
!===================================================
      do ii=ney*nex+1,nmax
        iy = mod((ii-nex*ney)-1,niy) + 1
        ix = (ii-nex*ney-1)/niy+1
!===================================================
!          print *,'I<-E presynaptic'
!===================================================
        rlx = (nex-1.D0)/2.D0 * dE
          rly = (ney-1.D0)/2.D0 * dE
        do j=1,ney
            do k=1,nex
            jj = ney*(k-1)+j
            dx = pix(iy,ix)-pex(j,k)
            dy = piy(iy,ix)-pey(j,k)
        if (dabs(dx).gt. rlx) then
            dx = 2.D0*mod(pix(iy,ix)+rlx,2.D0*rlx)-pex(j,k)
     1           - pex(iy,ix)
        endif
        if (dabs(dy).gt. rly) then
            dy = 2.D0*mod(piy(iy,ix)+rly,2.D0*rly)-pey(j,k)
     1           - pey(iy,ix)
        endif
            distance2 = dx**2 + dy**2
            prefix = de2*neffIE/(2.D0*pi*aie)
            aa = prob(prefix,distance2,aie)
            if (aa.GT.1.D0.or.ran2(iseed).LE.aa) then
              NpostSynI(jj) = NpostSynI(jj)+1
              NpreSynE(ii) = NpreSynE(ii)+1
              IpostSyn(jj,NpostSynE(jj)+NpostSynI(jj))=ii
              IcMatrix(jj,ii) = 1
            endif
            enddo
        enddo
!===================================================
!         print *,'I<-I presynaptic'
!===================================================
        rlx = (nix-1.D0)/2.D0 * dI
        rly = (niy-1.D0)/2.D0 * dI
        do j=1,niy
            do k=1,nix
            jj = nex*ney + niy*(k-1)+j
        if (ii.NE.jj) then
            dx = pix(iy,ix)-pix(j,k)
            dy = piy(iy,ix)-piy(j,k)
        if (dabs(dx).gt. rlx) then
            dx = 2*mod(pix(iy,ix)+rlx,2*rlx)-pix(j,k)
     1           - pix(iy,ix)
        endif
        if (dabs(dy).gt. rly) then
            dy = 2*mod(piy(iy,ix)+rly,2*rly)-piy(j,k)
     1           - piy(iy,ix)
        endif
            distance2 = dx**2 + dy**2
            prefix = di2*neffII/(2*aii*pi)
            aa = prob(prefix,distance2,aii)
            if (aa.GT.1.D0.or.ran2(iseed).LE.aa) then
              NpostSynI(jj) = NpostSynI(jj)+1
              NpreSynI(ii) = NpreSynI(ii)+1
              IpostSyn(jj,NpostSynE(jj)+NpostSynI(jj))=ii
              IcMatrix(jj,ii) = 1
            endif
        endif
            enddo
        enddo
      enddo
!=====================retro==============================
!            I->E postSynaptic
!========================================================
      do ii = nex*ney+1,nmax
          do jj = 1,nex*ney
      if (IcMatrix(jj,ii).EQ.1.and.IcMatrix(ii,jj).NE.1) then
          if (retroInh.GE.1.D0.or.
     1    (retroInh.GT.0.D0.and.ran2(iseed).LE.retroInh)) then
            NpostSynE(ii) = NpostSynE(ii)+1
            NpreSynI(jj) = NpreSynI(jj)+1
            IpostSyn(ii,NpostSynE(ii)+NpostSynI(ii))=jj
            IcMatrix(ii,jj) = 1
          endif
      endif
          enddo
      enddo
!===============================================
!        all to all
!===============================================
        else
      do i=1,nmax
          NpostSynE(i) = nmax*(4.D0/5.D0)
          NpostSynI(i) = nmax/5
          do j=1,nmax
            IpostSyn(i,j)=j
          enddo
      enddo
        endif
      npEE = 0 
      npEI = 0 
      npIE = 0 
      npII = 0 
      do i=1,nmax
        if (excite(i)) then
            npEE = npEE + NpreSynE(i)
            npEI = npEI + NpreSynI(i)
        else
            npIE = npIE + NpreSynE(i)
            npII = npII + NpreSynI(i)
        endif
      enddo

      nEE = 0 
      nEI = 0 
      nIE = 0 
      nII = 0 
      k = 0
      ddmean = 0.D0
      do i=1,nmax
        if (excite(i)) then
            nEE = nEE + NpostSynE(i)
            nIE = nIE + NpostSynI(i)
        else
            nEI = nEI + NpostSynE(i)
            nII = nII + NpostSynI(i)
        endif
        if (nlgni(i).GT.nr) then
            k = k + 1
            dmean = 0.D0
            do j = 1,NpostSynE(i)+NpostSynI(i)
          dtheta = dabs(thetas(i) - thetas(IpostSyn(i,j)))
          if (dtheta.GT.pi/2.D0) dtheta = pi-dtheta
            dmean = dmean + dtheta
            enddo
            ddmean = ddmean + 
     1        dmean/(1.D0*NpostSynE(i)+NpostSynI(i))
        endif
      enddo
      if (nr.LT.lgnmax) then
        print *,'mean dtheta', ddmean/(1.D0*k)
      endif
    
      open(16,file=fc,status='replace',form='unformatted',
     1       access='direct',recl=4*nmax)
      write(16,rec=1) (NpostSynE(i),i=1,nmax)
      close(16)
      open(16,file=fc,status='old',form='unformatted',
     1       access='direct',recl=4*nmax)
      write(16,rec=2) (NpostSynI(i),i=1,nmax)
      close(16)
      open(16,file=fc,status='old',form='unformatted',
     1       access='direct',recl=4*nmax)
      write(16,rec=3) (NpreSynE(i),i=1,nmax)
      close(16)
      open(16,file=fc,status='old',form='unformatted',
     1       access='direct',recl=4*nmax)
      write(16,rec=4) (NpreSynI(i),i=1,nmax)
      close(16)
    
      do i=5,nmax+4
      open(16,file=fc,status='old',form='unformatted',
     1       access='direct',recl=4*nmax)
      write(16,rec=i) (IcMatrix(i-4,j),j=1,nmax)
      close(16)
      enddo
    
      return
      end
!************************************************************
      subroutine lgnrf(tmpphi,conAmp)
!  Map each LGN cell to its receptive field :
!    location & spatiotemporal filter parameters
!
!  Determine each cell's firing rate given visual stimulus :
!
!      f(t)   = f0 + \int_0^t ds \int dy G(t - s) A(x - y) I(y,s)
!      firing rate at time t for LGN cell centered at x
!
!    where
!
!    f0     = background rate
!       G(t)   = t^5[exp(-t/t0)/t0^6 - exp(-t/t1)/t1^6]
!    A(y)   = ampa/siga/pi exp(-y^2/siga^2) 
!        - ampb/sigb/pi exp(-y^2/sigb^2) 
!         (overall minus one factor for off-center LGN cells)
!    I(y,s) = I_0 [ 1 + eps sin(omega s) cos( k y - phi )]
!        (for contrast reversal)
!------------ DG (Dai Wei) --------------------------------------------
!    I(y,s) = I_0 [1 + eps sin(omega s - k dotProduct y) ] 
!        (for drifting grating, expand sin in code, calculate the position part here)
!* Questions:
!    don't know how to perform this convolution
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N)
      common / lgnConvol / gk2, taulgn, omega
!------------------------------------------------------------
!     xlgn, ylgn: center of each LGN RF
!     siga, sigb: sigma of each Gaussian (RF = diff of 2 Gaussians)
!     ampa, ampb: amplitude of each Gaussian 
!     t0,   t1  : LGN kernel time constants
!     frtlgn    : background firing rate
!     ampson, ampsof: amplitude of (co)sine(phi - kx -ky) for on/off cell
!
!     9 Jan 2001: specialize to DG/CR for t >> t0
!    f(t) is a rectified sinusoid
!------------------------------------------------------------
!------------------------------------------------------------
!  Hardwire each LGN cell
!------------------------------------------------------------
      tau0 = taulgn
      tau1 = 5.D0/3.D0*tau0
      ot0  = omega*tau0
      ot1  = omega*tau1
      Gs = 240.D0*((3.D0 - 10.D0*ot0*ot0 + 3.D0*(ot0**4))*ot0/
     1      (1.D0+ot0*ot0)**6-(3.D0-10.D0*ot1*ot1+3.D0*(ot1**4))*ot1/
     2      (1.D0+ot1*ot1)**6)
      Gc = 120.D0*((-1.D0 + 15.D0*ot0*ot0 - 15.D0*(ot0**4) + ot0**6)
     1      * ot0 / (1.D0+ot0*ot0)**6
     2      - (-1.D0 + 15.D0*ot1*ot1 - 15.D0*(ot1**4) + ot1**6)
     3      * ot1 / (1.D0+ot1*ot1)**6)
      tmpphi = atan2(Gs,Gc)
      Aa = 14.88D0
      Ab = Aa*0.97D0
      rc = 5.61D0
      rs = 16.98D0*rc/5.61D0
      siga = rc/sqrt(2.D0)*4.D0*atan(1.D0)/180.D0
      sigb = rs/sqrt(2.D0)*4.D0*atan(1.D0)/180.D0
      A = Aa*exp(-gk2*siga**2/2.D0) - Ab*exp(-gk2*sigb**2/2.D0)
!------------------------------------------------------------
      DG_A = sqrt(Gc**2+Gs**2)
      conAmp = DG_A*A
      print *,'DG Amplitude = ', sngl(conAmp)
      return
      end
!************************************************************
      REAL*8 FUNCTION lgnsatur(x)
!------------------------------------------------------------
!
! Melinda's LGN contrast saturation function
!
!------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
!------------------------------------------------------------
      common / lgnSat / xS, yS, linear, rI0
      rr = xS*x
        if (linear.NE.1) then
      c2 = 0.09827D0
      c3 = -0.001589D0

      satmax = 60.D0
      rjoin=30.D0
      d1 = 40.D0
      pi = 4.D0*atan(1.D0)
      d4 = satmax-d1*pi/2.D0
      zz = tan((rjoin - d4)/d1) + pi
      d2 = (zz*zz+1)/d1
      d3 = zz - rjoin * d2
      if (rr.lt.41.D0) then
        lgnsatur = c3*(rr**3) + c2*(rr**2)
      else
        lgnsatur = d1 * atan(d2*rr + d3) + d4
      endif
        else
            lgnsatur = rr
        endif
      lgnsatur = lgnsatur*yS*rI0

!------------------------------------------------------------
      return
      end
!************************************************************
      double precision function herm_int(s,dt,v0,v1,vt0,vt1)
!------------------------------------------------------------
!  Mike's hermite polynomial interpolant
!------------------------------------------------------------
      implicit none
      double precision s,dt,v0,v1,vt0,vt1,a,b,c,d
      double precision s2,s3,dt22,dt3
!------------------------------------------------------------
      dt22=dt*dt
      dt3=dt22*dt
      a=v0
      b=vt0
      c=(3.D0*(v1-v0)-dt*(2.D0*vt0+vt1))/dt22
      d=(-2.D0*(v1-v0)+dt*(vt0+vt1))/dt3
      s2=s*s
      s3=s2*s
      herm_int = a+b*s+c*s2+d*s3
!------------------------------------------------------------
      return
      end
!************************************************************
      double precision function hermt_int(s,dt,v0,v1,vt0,vt1)
!------------------------------------------------------------
!  Mike's hermite (derivative) polynomial interpolant
!------------------------------------------------------------
      implicit none
      double precision s,dt,v0,v1,vt0,vt1,b,c,d
      double precision s2,dt22,dt3
!------------------------------------------------------------
      dt22=dt*dt
      dt3=dt22*dt
      b=vt0
      c=(3.D0*(v1-v0)-dt*(2.D0*vt0+vt1))/dt22
      d=(-2.D0*(v1-v0)+dt*(vt0+vt1))/dt3
      s2=s*s
      hermt_int = b+2.D0*c*s+3.D0*d*s2
!------------------------------------------------------------
      return
      end
!************************************************************
      subroutine newton(iflag,dtsp,vvv,vvt,dt,v0,v1,vt0,vt1,itmax)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      common / vconst / vthres,vreset,vexcit,vinhib,gleak,gleakI
        diff = (vvv-vthres)/vvt
        do it = 1,itmax
          dtsp = dtsp - diff
          if(dabs(diff/dt).le.1.D-08.and.
     1          dtsp.ge.0.D0.and.dtsp.le.dt)then
            iflag = 1
            return
          endif
          vvv = herm_int(dtsp,dt,v0,v1,vt0,vt1)
          vvt = hermt_int(dtsp,dt,v0,v1,vt0,vt1)
          diff = (vvv-vthres)/vvt
        enddo
        return
      end
      subroutine bisect_newton(iflag,dtsp,dt,v0,v1,vt0,vt1,
     1                          ta,tb,vvv,vvt,itmax)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-N) 
      common / vconst / vthres,vreset,vexcit,vinhib,gleak
      attempt_diff = (vvv-vthres)/vvt
      attempt = dtsp - attempt_diff
      do it=1,itmax
          if (dabs(attempt_diff*vvt).lt.dabs(vvv)*2.D0
     1      .or.(attempt-ta)*(tb-attempt).lt.0.D0) then
! deny newton's attempt, bisect
            diff = (tb-ta)*0.5D0
            dtsp = ta + diff
          else
! accept
            diff = attempt_diff
            dtsp = dtsp - diff
          endif
          if (dabs(diff/dt).le.1.D-08) then 
            iflag = 1
            return
          endif
          vvv = herm_int(dtsp,dt,v0,v1,vt0,vt1)
          vvt = hermt_int(dtsp,dt,v0,v1,vt0,vt1)
          attempt_diff = (vvv-vthres)/vvt
          attempt = dtsp - attempt_diff
          if (vvv.gt.0.D0) then
            tb = dtsp
          else
            ta = dtsp
          endif
      enddo 
      return
      end
      FUNCTION ran2(idum)
      INTEGER*4 idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2D-7,RNMX=1.-EPS)
      INTEGER*4 idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
!$OMP THREADPRIVATE(iv,iy,idum2) 
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
