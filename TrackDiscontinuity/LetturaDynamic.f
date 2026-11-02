      program LetturaDynamics

      Implicit none

      Logical linp

      Character*12 titolo,testo
      Character*100 linea
      Character(len=100) :: filename
 
      Integer*4 marcatori,nelem,npoin,nmats,nprop,ndofn,ndofn_T,
     *          nnode,ndime,ngaus,nstre,nevab,nsvab,nevab_T,nsvab_T,
     *          nvfix,nvfix_T,nTout,harmonic,idra,p,phaseC,PhaseCc,
     *          Grav,celleT,celleP,nelemp_max,nsn_start,nsn_iter,
     *          NnewSN

      Double Precision tfine,tsteps,Sl,vel_in,Lx_min,
     *           cont_y,ocean_y,xsub,vh2o,dist_marc,ht_max,
     *           PressMin,PressMax,TempMin,TempMax,
     *                cf_fix,xin,yin,dist
	                                       						   
       Parameter(marcatori=455489,										! Numero totale di marcatori 
     *          nelem=6813,												! Numero totale di Elementi
     *		npoin=13846,												! Numero totale di punti nodali
     *		nmats=1,													! Numero totale di materiali
     *		nprop=3,													! Numero di proprietà per ogni materiale
     *		ndofn=2,													! Numero di gradi di libertà  MOMENTO
     *		ndofn_T=1,												! Numero di gradi di libertà TEMPERATURA
     *		nnode=6,													! Numero di nodi per elemento
     *		ndime=2,													! Dimensioni del problema
     *		ngaus=7,													! Numero di punti di Gauss 
     *		nstre=3,													! Numero di componenti dello stress
     *          nevab=ndofn*nnode,
     *		nsvab=ndofn*npoin,
     *		nevab_T=ndofn_T*nnode,
     *          nsvab_T=ndofn_T*npoin,
     *          nvfix=633,												! Numero di punti a velocità fissata - 584 Attiva - 557 Grav
     *		nvfix_T=1376,												! Numero di punti a Temperatura fissata
     *          tfine=0.06D0,												! Tempo massimo del modello (Myr)
     *          tsteps=0.05D0,												! Time step (Myr)
     *          nTout=1,
     *          Sl=80.D03,												! Spessore termico iniziale della litosfera
     *          vel_in=5.D0,												! Velocita di ingresso dei marcatori crostali: ATTUALMENTE SIMMETRICA!!!!
     *          Lx_min=-276.0D0,			      								! Limite sinistro modello
     *          cont_y=-30.D0,												! Y della crosta continentale in km!
     *          ocean_y=-10.D0,											! Y della crosta oceanica in km!
     *          xsub=-100.D0,												! X del trench
     *          vh2o=20.D0,												! velocità dei fluidi in cm/yr
     *          dist_marc=0.5D0,											! Distanza tra i marcatori
     *          harmonic=0,												! Average type
     *          idra=1,													! Hydration yes/no
     *          p=10,													! Pressure quadrature points
     *          phaseC=1,												! Cambiamenti fase Mantello
     *          phaseCc=1,												! Cambiamenti fase Crosta
     *          Grav=0,													! Fase Gravitazionale/Attiva
     *          ht_max=2.5D0,											! Altezza massima della topografia
     *          celleT=100,												! Punti per la Temperatura di Perplex
     *          celleP=100,												! Punti per la Pressione di Perplex
     *          PressMin=0.D0,											! Pressione minima del file Perplex
     *          PressMax=7.D0,											! Pressione massima del file Perplex
     *          TempMin=300.D0,											! Temperatura minima del file Perplex
     *          TempMax=1600.D0,
     *          cf_fix=0.50D0, ! Valore fissato di coupling     
     *          nelemp_max=48,
     *          nsn_start=100,
     *          nsn_iter=100,
     *          dist=5.)											! Temperatura massima del file Perplex

      Integer*4 ipoin,ielem,idime,inode,ivfix,idofn,
     *          ievab,idofn_T,imarc,nmarc,kidra,tout,
     *          jdime,jevab,igaus,kgasp,matno,
     *          nrows,ncols,iadentro,jsvab,neqns,
     *          ivfix_T,nodpt,isvab,Nx,Ny,it,tst,
     *          numel,imats,numat,iprop,i,j,nodptl,d1,d2

      Double Precision Dens0,Dens_C,Dens_Cm,Dens_O,Dens_Om,Dens_a,
     *                 Dens_ou,Dens_aOu,Dens_S,Dens_aS,Dens_CL,
     *                 Hm,Hc_C,Hc_Cm,Hc_O,Hc_Om,Hc_Ou,H_E,Densita_E,
     *                 Hair,Hc_S,Viscosita_E,
     *                 Km,Kc_C,Kc_Cm,Kc_O,Kc_Om,Kc_a,Kc_w,Kc_Ou,
     *                 Ts,Tb,Temp0_D,Temp0_V,D_E,
     *                 Visc0,Visc0_m,Visc0_oc,Visc0_cc,Visc0_a,
     *		       Visc0_s,Dens_w,Visc0_ocu,Capacita_a,expans0,
     *                 Lo,Ea_m,Ea_oc,Ea_cc,Ea_ocu,Diff0,Diff_a,Capacita,
     *                 gravita,pivoti,teta,Ra,Rc,Ra_a,Rc_a,Aexp,
     *                 timeo,tfin,tstep,Uo,Ho,Ko,K_E,coordx,coordy,
     *                 Area1,Area2,Area3,E1,E2,Emax,
     *                 x,y,ExxD,EyyD,ExyD,Ee,Tr_V,Tr_D,Ts_N,
     *                 SxxD,SyyD,SxyD,Vr,ViscY,I2,dissi,Had,Hd,
     *                 produzione,conduttivita,adiabatic,        
     *		       Uxsomma,Uxmedia,Uysomma,Uymedia,
     *                 hx,hy,kappax,kappay,alfax,alfay,Pex,Pey,
     *                 Ax,Ay,h,Pe,alfa,dissipazione,Di,
     *                 Temper,dvolu,s,t,
     *                 U1,U2,time,djacb,densita,viscosita

      Integer*4 nofix(nvfix),ifpre(nvfix,ndofn),lnods(nelem,nnode),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T),
     *          Numelem(marcatori),C0(nelem),nelemp(npoin),
     *          lelemp(npoin,nelemp_max),
c     *          lpoins(npoin,nelemp_max),
     *          lpoins(npoin,npoin),
     *          nlpoins(npoin)

      Double Precision, Dimension(nelem) :: D0P,DCP,DOP,DOUP,DSP,
     *                  DensCm,DensOm,Hc,Kc,DcT,Xe,Ye,Te,xmin,xmax,
     *                  DensitaC,DDEns,DCLP,Visc,composiz,Area,
     *                  ymin,ymax,Po,Pelem,Eps

      Double Precision, Dimension(nnode) :: shape,V1loc,V2loc

      Double Precision, Dimension(3) :: posigp,weiigp

      Double Precision, Dimension(nsvab) :: V1,V2

      Double Precision, Dimension(nevab,nevab) :: estif_T_L,estif_L

      Double Precision, Dimension(nevab_T,nevab_T) :: damp,estif_T,
     *                  estif_T_V

      Double Precision, Dimension(ndime,ndime) :: xjacm,xjaci

      Double Precision, Dimension(ndime,nnode) :: elcod,cartd,
     *                  cartd2_x,cartd2_y

      Double Precision, Dimension(nelem,nnode) :: x_pol,y_pol

      Double Precision, Dimension(celleP,celleT) :: DensPerpUC,
     *                  DensPerpCL,DensPerpS,DensPerpO,TPerp,PPerp

      Double Precision Temp(npoin),weigp(ngaus),coord(npoin,ndime),
     *                 presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T),
     *                 posgp(ndime,ngaus),
     *                 estif_V(ndime*nnode,ndime*nnode),
     *                 eload_T(nelem,nevab_T), 
     *                 eload_V(nelem,nnode*ndime),
     *                 stiff(nelem,nevab,nevab),
     *                 dampi(nelem,nevab_T,nevab_T),
     *                 stiff_T(nelem,nevab_T,nevab_T)

      Integer*4, dimension (:), allocatable :: iffixed,iffixed_T,M,
     *           Pkc,Pko,Pkcm,Pkom,Pka,Pkw,Pkou,PKcs,Pkcl

      Double Precision, Dimension(:), allocatable :: eloadg_V,
     *                  eloadg_T,noto,fixed,fixed_T,Vdisp,Tdisp,
     *                  aslod,aslod_T,Xm,Ym,Tm,Pm,Pk,
     *                  Cc,Co,Ccm,Ca,Cou,Cw,Com,Ccs,Ccl

      Double Precision, Dimension(:,:), allocatable :: estifg_V,
     *                  estifg_T_V,deriv,deriv2_s,deriv2_t,dampe,
     *                  astif,astif_T,adampi_T,dampg

c______________________________________________________________________________________ 
c      Variabili per Slip Nodes

      Integer*4 tip(nelem),nsn,snf(npoin)

      Double Precision, Dimension(nelem,nnode) :: cfl

      Double Precision, Dimension(npoin) :: cf,sn


        
      nsn=0

c_______________________________________________________________________________________ 	

c    1-  LETTURA DEI PARAMETRI E FILES DI INPUT
c________________________________________________________________________________________ 
      
      call parametri(Dens0,Dens_C,Dens_Cm,Dens_O,Dens_Om,
     *                  Dens_a,Dens_w,Dens_Ou,Dens_S, 
     *                  Hm,Hc_C,Hc_Cm,Hc_O,Hc_Om,Hair,Hc_S,Hc_Ou,
     *                  Km,Kc_C,Kc_Cm,Kc_O,Kc_Om,Kc_a,Kc_w,Kc_Ou,
     *                  Ts,Tb,Temp0_D,Temp0_V,Tr_V,Tr_D,Ts_N,
     *                  Visc0,Visc0_m,Visc0_oc,Visc0_cc,Visc0_a,
     *	        	      Visc0_s,Visc0_ocu,Dens_CL,
     *                  Lo,Ea_m,Ea_oc,Ea_ocu,Ea_cc,Diff0,Diff_a,
     *                  Capacita,Capacita_a,expans0,
     *                  gravita,pivoti,teta,Ra,Ra_a,Rc,Rc_a,Aexp,
     *                  timeo,Uo,Ho,Ko,Di)        
   
c---------------------------------------------------------------------------
c     LETTURA DATI DI DEFINIZIONE DEL PROBLEMA DAL FILE INPUT

      inquire(file='Input',exist=linp)
      if(.not.linp) print *,'ERRORE - il file Input non esiste!'
      
      open(unit=15, file='Input', status='old')

c      print *,'Nodes =',npoin
c      print *,'Elements =',nelem
c      print *,'Markers =',marcatori
c      print *,'Fixed velocities =',nvfix
c      print *,'Fixed temperatures =',nvfix_T
c      print *,'Initial lithospheric thickness (km) =',Sl/1000.
c      print *,'Continental moho depth (km) =',cont_y
c      print *,'Oceanic moho depth (km) =',ocean_y
c      print *,'Time step (Ma)=',tsteps
c      print *,'Final time (Ma)=',tfine
c      print *,'Plate velocity (cm/y) =',vel_in
	  
c      if(idra.eq.1) then
c        print *,'Hydrated mantle wedge'
c      elseif(idra.eq.0) then
c        print *,'Dry mantle wedge'
c      endif

c      if(harmonic.eq.1) then
c        print *,'Geometric average'
c      elseif(harmonic.eq.0) then
c        print *,'Arithmetic average'
c      endif            
                                          
      read(15,1009) testo

      do i=1,npoin
        read(15,*) ipoin,(coord(ipoin,idime),idime=1,ndime),							! Lettura coordinate Nodali
     *             sn(ipoin),cf(ipoin)
        coord(ipoin,1)=coord(ipoin,1)/700.   										! Coordinate normalizzate
        coord(ipoin,2)=coord(ipoin,2)/700.      
      enddo                                                

      read(15,1009) testo

      Xe=0.D0
      Ye=0.D0
      do ielem=1,nelem
       read(15,*) numel,tip(numel),(lnods(numel,inode),
     *    inode=1,nnode)  ! Lettura: Numero Elemento, Tipo di Elemento, Connessioni nodali per ogni elemento
                											! Calcolo proprietà geometriche dell'elemento


       do inode=1,nnode
        if(tip(ielem).eq.0 .or. tip(ielem).eq.-1) then
         cfl(ielem,inode)=1.0D0
        elseif(tip(ielem).eq.+1) then
         cfl(ielem,inode)=cf(lnods(ielem,inode))
       endif
        Xe(ielem)=Xe(ielem)+(coord(lnods(ielem,inode),1))*700.
        Ye(ielem)=Ye(ielem)+(coord(lnods(ielem,inode),2))*700.
	     x_pol(ielem,inode)=(coord(lnods(ielem,inode),1))*700.   						! Coordinate nodali in km
	     y_pol(ielem,inode)=(coord(lnods(ielem,inode),2))*700.
       enddo
       Xe(ielem)=Xe(ielem)/6												! Baricentro coordinate normalizzato
       Ye(ielem)=Ye(ielem)/6
       Area(ielem)=dabs((1.*(x_pol(ielem,2)*y_pol(ielem,3)-							! Area elemento in km^2
     *                    x_pol(ielem,3)*y_pol(ielem,2))-
     *                    x_pol(ielem,1)*(y_pol(ielem,3)-
     *                    y_pol(ielem,2))+y_pol(ielem,1)*
     *                    (x_pol(ielem,3)-x_pol(ielem,2))))/2
       Po(ielem)=(Area(ielem)/(dist_marc**2))           								! Marcatori massimi dell'elemento
       xmin(ielem)=MIN(x_pol(ielem,1),x_pol(ielem,2),x_pol(ielem,3))
       xmax(ielem)=MAX(x_pol(ielem,1),x_pol(ielem,2),x_pol(ielem,3))        
       ymin(ielem)=MIN(y_pol(ielem,1),y_pol(ielem,2),y_pol(ielem,3))
       ymax(ielem)=MAX(y_pol(ielem,1),y_pol(ielem,2),y_pol(ielem,3))
      enddo			
c -------------------------------------------------------------------------------
c         CREAZIONE CONNESSIONI NELEMP NPOINS
c      call creazione_lelemp_lpoins(nelem,npoin,nnode,ndime,
c     *           nelemp_max,coord,lnods,nelemp,lelemp,
c     *           nlpoins,lpoins)

c      call creazione_lelemp_lpoins_NV(nelem,npoin,nnode,ndime,
c     *           nelemp_max,coord,lnods,nelemp,lelemp,
c     *           nlpoins,lpoins)

c      call creazione_lelemp_lpoins_NV2(nelem,npoin,nnode,ndime,
c     *           nelemp_max,coord,lnods,nelemp,lelemp,
c     *           nlpoins,lpoins)

      call creazione_lelemp_lpoins_NV3(nelem,npoin,nnode,ndime,
     *           nelemp_max,coord,lnods,nelemp,lelemp,
     *           nlpoins,lpoins,dist/700.)

        
c -------------------------------------------------------------------------------
c        NORMALIZATION OF FINAL TIME, TIME-STEP e SPESSORE LITOSFERA
	   
	tfin=tfine/timeo    
	tstep=tsteps/timeo          
c_________________________________________________________________________________________ 
c     APERTURA FILES SCRITTURA FAGLIA

      open(unit=28,file='Fault_nodes_1.dat',status='unknown')
      close(unit=28,status='delete')
      open(unit=28,file='Fault_nodes_1.dat',status='new')
      close(unit=28)

      open(unit=30,file='Fault_nodes_2.dat',status='unknown')
      close(unit=30,status='delete')
      open(unit=30,file='Fault_nodes_2.dat',status='new')
      close(unit=30)

c________________________________________________________________________________________ 
c    2-  INIZIO DEL CICLO TEMPORALE
c_________________________________________________________________________________________ 
      print *,''
      print *,'INIZIA IL TIME STEP TEMPORALE'
 
107   time=0.0D0

      Tout=nTout+1

      do 1 it=1,INT(tfin/tstep)           

      print *,"Iterazione = ",it

      tst=(it-1)

      print *,'Time=',time*timeo,'Ma'

c________________________________________________________________________________________ 
c
c    3-  LETTURA RISULTATI      
c_________________________________________________________________________________________ 
      open(unit=28,file='Fault_nodes_1.dat',access='append',
     *   status='old')

      open(unit=30,file='Fault_nodes_2.dat',access='append',
     *   status='old')

      if(it.eq.1.or.it.eq.Tout) then

c      write(filename,'("I2.Nodi.",I3.3,"Ma.t")') tst
c      open(unit=22,file=filename,status='unknown')

      open(unit=22,file='I2.Nodi.tmp',status='unknown')


c      write(filename,'("Vel_Temp.",I3.3,"Ma.t")') tst
C      open(unit=16,file=filename,status='old')

      open(unit=16,file='Vel_Temp.tmp',status='old')
 
cc      read(16,*) testo
c      read(16,*) testo
c      read(16,*) testo
c      read(16,*) testo
c      read(16,*) testo
      read(16,*) testo
      read(16,*) testo
      
c      print *,testo
      
      do i=1,npoin
c           print *,i
           read(16,*) ipoin,coord(ipoin,1),coord(ipoin,2),
     *	                     V1(ipoin),V2(ipoin),Temp(ipoin)
     
            V1(ipoin)=V1(ipoin)/Uo
            V2(ipoin)=V2(ipoin)/Uo
c         read(16,33333) ipoin,coord(ipoin,1),coord(ipoin,2),
c     *	               (cf(ipoin)*V1(ipoin))*Uo,
c     *                 (cf(ipoin)*V2(ipoin))*Uo,Temp(ipoin)

      end do

      close(16)
     
c      write(filename,'("Strain_Stress_Hd.",I3.3,"Ma.t")') tst    
C      open(unit=24,file=filename,status='old')

      open(unit=24,file='Strain_Stress.tmp',status='old')
 
c      print*,filename  
      read(24,*) testo
      read(24,*) testo

c      print *,testo
      
      do i=1,nelem
         read(24,999) ielem,x,y,
     *            ExxD,EyyD,ExyD,Ee,
     *            SxxD,SyyD,SxyD,
     *            dissi,
     *            Hc(ielem),		!Scrive la Dissipazione e produzione radiogenica in W/m^3
     *            Visc(ielem),
     *            Vr,ViscY,I2,Had
         E1=(ExxD+EyyD)/2.D0+dsqrt((((ExxD-EyyD)**2)/4.D0)+ExyD**2)
         E2=(ExxD+EyyD)/2.D0-dsqrt((((ExxD-EyyD)**2)/4.D0)+ExyD**2)
         Emax=0.5D0*(E1-E2)
c         E1=(SxxD+SyyD)/2.D0+dsqrt((((SxxD-SyyD)**2)/4.D0)+SxyD**2)
c         E2=(SxxD+SyyD)/2.D0-dsqrt((((SxxD-SyyD)**2)/4.D0)+SxyD**2)
c         Emax=0.5D0*(E1-E2)

         Eps(ielem)=Emax
      end do

      close(24)
      end if
c________________________________________________________________________________________ 
c
c    4-  DYNAMIC SLIP NODE
c_________________________________________________________________________________________ 	

c      call dynamic_sn_time(nelem,npoin,V1,V2,Uo,Lo,
c     *           nsn_start,nsn_iter,cf_fix,  
c     *           nelemp_max,nnode,ndime,
c     *           nelemp,lelemp,nlpoins,lpoins,
c     *           coord,lnods,Xe,Ye,time*timeo,
c     *           -275.D0,700.D0,-700.D0,0.D0,xin,yin,
c     *           Eps,nsn,snf,it,tip,cfl,sn,cf)

c      call  dynamic_sn_istantaneo(nelem,npoin,V1,V2,Uo,Lo,
c     *           nsn_start,nsn_iter,cf_fix,  
c     *           nelemp_max,nnode,ndime,
c     *           nelemp,lelemp,nlpoins,lpoins,
c     *           coord,lnods,Xe,Ye,time*timeo,
c     *           -275.D0,700.D0,-700.D0,0.D0,xin,yin,
c     *           Eps,it,tip,nsn,snf,cfl,sn,cf,NnewSN)

      call  dynamic_sn_istantaneo_nv(nelem,npoin,V1,V2,Uo,Lo,
     *           nsn_start,nsn_iter,cf_fix,  
     *           nelemp_max,nnode,ndime,
     *           nelemp,lelemp,nlpoins,lpoins,
     *           coord,lnods,Xe,Ye,time*timeo,
     *           -275.D0,700.D0,-700.D0,0.D0,xin,yin,
     *           Eps,it,tip,nsn,snf,cfl,sn,cf,NnewSN)                        

      close(28)
      close(30)

      time=time+tstep

      if(it.eq.Tout) then
       	                  
       Tout=Tout+nTout    

      endif
            
      
      !SE TIME < Tfin 	PASSA AL TIME STEP SUCCESSIVO
      !ALTRIMENTI	INTERROMPI IL PROGRAMMA
       
      if(it.eq.INT(tfin/tstep)) print *,' *** PROGRAM END ***' 
      if(it.eq.INT(tfin/tstep)) stop


      print *,'END TIME LOOP' 

1      continue

C   FINE DEL N CICLO TEMPORALE E PASSAGGIO AL N+1

c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
C  !!!!!!!!!!!!!!!!!!!!!!!!     PROGRAMMA CONCLUSO       !!!!!!!!!!!!!!!!!!!!!!!!!

c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
c_________________________________________________________________________________________  
c______________________________________________________________________________________ 

c    10-  FORMATI DEI FILES
c_______________________________________________________________________________________  

999       format(i10,f10.3,f10.3,15(2x,g10.4))           
1009  format(a12)
22222 format(4x,4HNODO,3x,7H   X   ,3x,7H   Y   ,3x,
     *       7H   U   ,3x,7H   V   ,3x,7H  Temp ,3x,
     *       7H  Visc ,3x,7H  Dens /)    
33333 format(3x,I5,1x,f8.3,1x,f8.3,1x,g10.3,
     * 1x,g10.3,2x,f8.3,2x,g9.3,2x,f8.3)
20000 format(1x,6HTempo=,f5.1,2HMa)
22223 format(4x,9HMARCATORE,3x,7H   X   ,3x,7H   Y   ,3x,
     *       7H   M   /)
22224 format(4x,8HELEMENTO,3x,7H   Xe  ,3x,7H   Ye  ,3x,
     *       7H   C   ,9HDens(T,C),7HDens(C),11Hlog10(Visc),2HKe,2HHe/)
22225 format(10H      ELEM,6H    Xe,10H        Ye,
     *       12H         Exx,12H         Eyy,12H         Exy,
     *       12H        Eeff,
     *       12H         Sxx,12H         Syy,12H         Sxy,
     *       12H          Hd,12H          Hc,
     *       16H      Visc Usata,12H  Visc Reale,9H    ViscY/)
12345 format(i7,7(1x,f10.2),1x,g10.2,2(1x,f10.2),i5)  
22226 format(6H    Xe,6H    Ye,/)
 
      stop
      end

c====================================================================================

c_______________________________________________________________________________________  
c  	  SUBROUTINES


      subroutine parametri(Dens0,Dens_C,Dens_Cm,Dens_O,Dens_Om,
     *                     Dens_a,Dens_w,Dens_Ou,Dens_S,        
     *                     Hm,Hc_C,Hc_Cm,Hc_O,Hc_Om,Hair,Hc_S,Hc_Ou,
     *                     Km,Kc_C,Kc_Cm,Kc_O,Kc_Om,Kc_a,Kc_w,Kc_Ou,
     *                     Ts,Tb,Temp0_D,Temp0_V,Tr_V,Tr_D,Ts_N,
     *                     Visc0,Visc0_m,Visc0_oc,Visc0_cc,Visc0_a,
     *			           Visc0_s,Visc0_ocu,Dens_LC,
     *                     Lo,Ea_m,Ea_oc,Ea_ocu,Ea_cc,Diff0,Diff_a,
     *			           Capacita,
     *                     Capacita_a,expans0,
     *                     gravita,pivoti,teta,Ra,Ra_a,Rc,Rc_a,Aexp,
     *                     timeo,
     *                     Uo,Ho,Ko,Di)       
      
       implicit none
       
       Double Precision Dens0,Dens_C,Dens_Cm,Dens_O,Dens_Om,Dens_a,
     *                     Dens_w,Hm,Hc_C,Hc_Cm,Hc_O,Hc_Om,Hair,Hc_S,
     *                     Km,Kc_C,Kc_Cm,Kc_O,Kc_Om,Kc_a,Kc_w,
     *                     Ts,Tb,Temp0_D,Temp0_V,
     *                     Visc0,Visc0_m,Visc0_oc,Visc0_cc,Visc0_a,
     *			   Visc0_s,Dens_S,
     *                     Lo,Ea_m,Ea_oc,Ea_cc,Diff0,Diff_a,Capacita,
     *                     Capacita_a,expans0,
     *                     gravita,pivoti,teta,Ra,Ra_a,Rc,Rc_a,Aexp,
     *                     timeo,Dens_LC,
     *                     Uo,Ho,Ko,Tr_V,Tr_D,Ts_N,Di,
     *		           Dens_Ou,Hc_Ou,Kc_Ou,Visc0_ocu,Ea_ocu
    
c----------------------------------------------------------------------------------------- 
c     Parametri dimensionali
      
      Dens0=3200.D0  !Densità di riferimento del mantello alla Temperatura Temp0_D (kg/m^3) 
      Dens_C=2640.   !Densità della crosta CONTINENTALE - 2/3 gneis + 1/3 granite (kg/m^3)
      Dens_LC=2640.   !Densità della crosta CONTINENTALE - 2/3 gneis + 1/3 granite (kg/m^3) 
      Dens_S=2640.   !Densità dei SEDIMENTI - 2/3 gneis + 1/3 granite (kg/m^3) 
      Dens_Cm=3200.  !Densità della mantello CONTINENTALE  - 100% rhyolite  (kg/m^3) 
      Dens_O=2961.   !Densità della crosta OCEANICA - 7%basalto+16%dolerite+77%gabbro(kg/m^3)  
      Dens_Om=3200.  ! Densità della mantello OCEANICA  - 100% basalt MORB (kg/m^3) 
      Dens_a=1.18    ! Densità dell'aria
      Dens_w=1000.   ! Densità dell'acqua
      Dens_Ou=2961.  ! Densità della crosta OCEANICA superiore - 7% basalto + 16% dolerite + 77% gabbro	(kg/m^3) 
      
      Hm=0.002D-06/Dens0    !Produzione di calore mantello per unita' di massa (W/kg)      
      Hc_C=2.5D-06/Dens_C   !Produzione di calore crosta CONTINENTALE per unita' di massa      
      Hc_Cm=0.002D-06/Dens0  !Produzione di calore mantello CONTINENTALE per unita'di m.     
      Hc_O=0.4D-06/Dens_O   !Produzione di calore crosta OCEANICA per unita' di massa 
      Hc_Om=0.002D-06/Dens0 !Produzione di calore mantello OCEANICA  per unita' di massa
      Hair=0.                 !Produzione di calore aria per unita' di massa 
      Hc_S=2D-06/Dens_C     !Produzione di calore sedimenti per unita' di massa
      Hc_Ou=0.4D-06/Dens_Ou  ! Produzione di calore crosta OCEANICA superiore per unita' di massa (W/kg)        
      
      Km=4.15D0			! Conduttività termica mantello 		(W/m°K)      
      Kc_C=3.03D0		! Conduttività termica crosta CONTINENTALE 	(W/m°K)      
      Kc_Cm=4.15D0		! Conduttività termica mantello CONTINENTALE (W/m°K)      
      Kc_O=2.10D0		! Conduttività termica crosta OCEANICA 		(W/m°K)      
      Kc_Om=4.15D0		! Conduttività termica mantello OCEANICA  	(W/m°K)
      Kc_a=0.026                ! Conduttivita termica aria
      Kc_w=0.6                  ! Conduttivita termica acqua
      Kc_Ou=2.10D0		! Conduttività termica crosta OCEANICA superiore(W/m°K)      
     
      Ts=300.D0			! Temperatura alla superficie 			(°K)      
      Tb=2300.D0		! Temperatura alla base 			(°K)     
      Temp0_D=1600.D0		! Temperatura di riferimento per la Densità 	(°K)
      Temp0_V=1600.D0		! Temperatura di riferimento per la Viscosità 	(°K)

      Lo=700.*1.D03		! Profondità della cella (m) 
     
      Capacita=1.D03		! Capacità termica 		(J/kg°K)
      Capacita_a=1.005           !Capacita termica aria 
      
      expans0=3.5D-05		! Coefficiente di espansione termica  di mantello (°K^-1)
     
      gravita=9.81D0 		! Accelerazione di gravità 	(m/s^2)
         
      pivoti=1.D03		! Fattore di Penalty per la pressione
          
      teta=2./3.		! Galerkin: teta=2/3;  Crank Nicolson: teta=1/2
      
      Diff0=1.D-06		! Diffusività termica 	(m^2/s)
      Diff_a=2.4D-05            ! Diffusività termica aria (m^2/s)

c   PARAMETRI REOLOGICI

      Visc0=0.5D21	  !Viscosità di riferimento alla Temperatura Temp0_V     (Pa s^-1)   
      Visc0_m= 5.012098297276852D20	!Mantle 
      Visc0_oc=1.612817973471232D22	!Oceanic crust
      Visc0_cc=3.474008984368633D21	!Continental crust
      Visc0_a=1.0D19		                        !Air end Water
      Visc0_s=3.474008984368633D20   	!Sediments
      Visc0_ocu=1.612817973471232D19	!Oceanic upper crust 
      
      Ea_m=130.2052785923754D03	        !Mantle activation energy (J mol^-1)     
      Ea_oc=103.1914893617021D03    	!Ocenic crust activation energy
      Ea_cc=38.43750000000000D03    	!Continental crust activation energy
      Ea_ocu=103.1914893617021D03		!Ocenic upper crust
      
c---------------------------------------------------------------------------------------- 
c     Parametri derivati adimensionali   
         
      Ra=(Dens0*Gravita*expans0*(Tb-Ts)*(Lo**3))/(Diff0*Visc0) !Numero di Rayleigth "termico"  

c     Ra_a=(Dens_a*Gravita*expans0*(Tb-Ts)*(Lo**3))/(Diff_a*Visc0_a)!Numero di Rayleigh aria 
  
      Rc=(Lo**3)*Gravita/(Diff0*Visc0) 	      !Numero di Rayleigth composizionale          
c      Rc_a=(Lo**3)*Gravita/(Diff_a*Visc0_a)  !Numero di Rayleigth composizionale aria
      
      Aexp=-Dlog(100.D0)*(Tb-Ts)/(Ts-Temp0_D) !Parametro per la dipendenza della viscosità
                                              !dalla temperatura

      Tr_V=(Temp0_V-Ts)/(Tb-Ts)	! Temperatura (Normalizzata) di riferimento per la Viscosità
      Tr_D=(Temp0_D-Ts)/(Tb-Ts)	! Temperatura (NOrmalizzata) di riferimento per la densita
      Ts_N=0.0	! Temperatura superficiale (Normalizzata)

      Di=expans0*Gravita*Lo/Capacita

c---------------------------------------------------------------------------------------- 
c     Fattori di adimensionalizzazione
           
      timeo=(Lo**2)/(Diff0*(1.D06*365*24*60*60)) !Fattore di Normalizzazione per il Tempo
     
      Uo=(Diff0*(365*24*60*60))/(Lo*(1.D-02)) !Fattore di Normalizzazione per la Velocita   
      
      Ho=Capacita*Diff0*(Tb-Ts)/(Lo**2)	!Fattore di Normalizzazione per la produzione di calore
      
      Ko=Dens0*Capacita*Diff0 !Fattore di Normalizzazione per la conduttività termica
      
      return
      end

c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      
      subroutine creazione_lelemp_lpoins_NV(nelem,npoin,nnode,ndime,
     *           nelemp_max,coord,lnods,nelemp,lelemp,
     *           nlpoins,lpoins)
      
      implicit none

      integer*4 npoin,ipoin,i,j,ielem,nelem,idime,ndime,nnode,
     *            inode,k,nelemp_max
      
      Integer*4 lnods(nelem,nnode),numel,nelemp(npoin)
        
      Double precision coord(npoin,ndime)

      integer*4 lelemp(npoin,nelemp_max),lpoins(npoin,nelemp_max),
     *              nlpoins(npoin)
        
        
      character*12 testo

      open(unit=21,file='Corrispondenze_Nodi2Elementi')

      open(unit=234,file='Connessioni_Dirette_Nodi')
  
      lelemp=0
      lpoins=0
 

          
      do ipoin=1,npoin

        k=0
        do ielem=1,nelem
            
          do inode=1,nnode
            if(ipoin.eq.lnods(ielem,inode)) then
              k=k+1
              lelemp(ipoin,k)=ielem
            endif 
          
          enddo
        
          nelemp(ipoin)=k        
        enddo
        
        
1     write(21,*) ipoin,nelemp(ipoin),
     *       (lelemp(ipoin,i),i=1,nelemp(ipoin))
      
      enddo 


        do ipoin=1,npoin
        
     
          k=0
          do i=1,nelemp(ipoin)


            do inode = 1,nnode            
              if(ipoin.eq.lnods(lelemp(ipoin,i),inode)) then
                if(inode.eq.1) then
                    
                    
        do j=1,k
              if(lnods(lelemp(ipoin,i),4)
     *                  .eq.lpoins(ipoin,j)) then
                    go to 100
                  end if
                end do
                k=k+1
                lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
100             continue

        do j=1,k
              if(lnods(lelemp(ipoin,i),5)
     *                  .eq.lpoins(ipoin,j)) then
                    go to 106
                  end if
                end do
                k=k+1
                lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)
106             continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do


                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)                            
                        go to 200
                        




                    else if(inode.eq.2) then
                    
                        do j=1,k
      if(lnods(lelemp(ipoin,i),4).eq.lpoins(ipoin,j)) then
                                go to 101
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
101                     continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                go to 107
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)
107                     continue

                        do j=1,k
      if(lnods(lelemp(ipoin,i),5).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)                            
                        go to 200
                    

                    else if(inode.eq.3) then
                    
                    
                        do j=1,k
      if(lnods(lelemp(ipoin,i),5).eq.lpoins(ipoin,j)) then
                                go to 102
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)
102                     continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),4).eq.lpoins(ipoin,j)) then
                                go to 108
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
108                     continue


                        do j=1,k    
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do    
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)                            
                        go to 200
                        

                    else if(inode.eq.4) then
                    
                    
                        do j=1,k
      if(lnods(lelemp(ipoin,i),1).eq.lpoins(ipoin,j)) then
                                go to 103
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),1)
103                     continue
  
  
                        do j=1,k
      if(lnods(lelemp(ipoin,i),3).eq.lpoins(ipoin,j)) then
                                go to 109
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),3)
109                     continue
  
  
                            
                        do j=1,k    
      if(lnods(lelemp(ipoin,i),2).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),2)  
                                        
                        go to 200
                        

                    else if(inode.eq.5) then
                    
                    
                        do j=1,k
      if(lnods(lelemp(ipoin,i),2).eq.lpoins(ipoin,j)) then
                                go to 104
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),2)
104                     continue

                        do j=1,k
      if(lnods(lelemp(ipoin,i),1).eq.lpoins(ipoin,j)) then
                                go to 110
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),1)
110                    continue



                            
                        do j=1,k    
      if(lnods(lelemp(ipoin,i),3).eq.lpoins(ipoin,j)) then
                                 go to 200
                            end if
                        end do

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),3)                            
                        go to 200
                        

                    else if(inode.eq.6) then

                        do j=1,k
      if(lnods(lelemp(ipoin,i),1).eq.lpoins(ipoin,j)) then
                                go to 105
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),1)
105                     continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),2).eq.lpoins(ipoin,j)) then
                                go to 111
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),2)
111                     continue

                            
                        do j=1,k    
      if(lnods(lelemp(ipoin,i),3).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do    
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),3)                            
                        go to 200
                        
                    end if
                end if

            end do
200         continue


        end do
        nlpoins(ipoin)=k
        enddo
        
        
        do ipoin=1,npoin
          write(234,*) ipoin,nlpoins(ipoin),
     *       (lpoins(ipoin,k),k=1,nlpoins(ipoin))
        
        enddo
      
      close(unit=21)
      close(unit=234)
      
      return
      end
      
      
c--------------------------------------------------------------------    
      
      subroutine creazione_lelemp_lpoins_NV2(nelem,npoin,nnode,ndime,
     *           nelemp_max,coord,lnods,nelemp,lelemp,
     *           nlpoins,lpoins)
      
      implicit none

      integer*4 npoin,ipoin,i,j,ielem,nelem,idime,ndime,nnode,
     *            inode,k,nelemp_max
      
      Integer*4 lnods(nelem,nnode),numel,nelemp(npoin)
        
      Double precision coord(npoin,ndime)

      integer*4 lelemp(npoin,nelemp_max),lpoins(npoin,nelemp_max),
     *              nlpoins(npoin)
        
        
      character*12 testo

      open(unit=21,file='Corrispondenze_Nodi2Elementi')

      open(unit=234,file='Connessioni_Dirette_Nodi')
  
      lelemp=0
      lpoins=0
 

          
      do ipoin=1,npoin

        k=0
        do ielem=1,nelem
            
          do inode=1,nnode
            if(ipoin.eq.lnods(ielem,inode)) then
              k=k+1
              lelemp(ipoin,k)=ielem
            endif 
          
          enddo
        
          nelemp(ipoin)=k        
        enddo
        
        
1     write(21,*) ipoin,nelemp(ipoin),
     *       (lelemp(ipoin,i),i=1,nelemp(ipoin))
      
      enddo 


        do ipoin=1,npoin
        
     
          k=0
          do i=1,nelemp(ipoin)


            do inode = 1,nnode            
              if(ipoin.eq.lnods(lelemp(ipoin,i),inode)) then
                if(inode.eq.1) then
                    
                    
        do j=1,k
              if(lnods(lelemp(ipoin,i),4)
     *                  .eq.lpoins(ipoin,j)) then
                    go to 100
                  end if
                end do
                k=k+1
                lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
100             continue

        do j=1,k
              if(lnods(lelemp(ipoin,i),5)
     *                  .eq.lpoins(ipoin,j)) then
                    go to 106
                  end if
                end do
                k=k+1
                lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)
106             continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do


                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)                            
                        go to 200
                        




                    else if(inode.eq.2) then
                    
                        do j=1,k
      if(lnods(lelemp(ipoin,i),4).eq.lpoins(ipoin,j)) then
                                go to 101
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
101                     continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                go to 107
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)
107                     continue

                        do j=1,k
      if(lnods(lelemp(ipoin,i),5).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)                            
                        go to 200
                    

                    else if(inode.eq.3) then
                    
                    
                        do j=1,k
      if(lnods(lelemp(ipoin,i),5).eq.lpoins(ipoin,j)) then
                                go to 102
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)
102                     continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),4).eq.lpoins(ipoin,j)) then
                                go to 108
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
108                     continue


                        do j=1,k    
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do    
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)                            
                        go to 200
                        
c---------------Nodi intermedi------------------------------

                    else if(inode.eq.4) then
                    
                    
                        do j=1,k
      if(lnods(lelemp(ipoin,i),1).eq.lpoins(ipoin,j)) then
                                go to 103
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),1)
103                     continue
  
  
                         do j=1,k    
      if(lnods(lelemp(ipoin,i),2).eq.lpoins(ipoin,j)) then
                                go to 112
                            end if
                        end do

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),2)  
112                    continue
 
                        do j=1,k
      if(lnods(lelemp(ipoin,i),3).eq.lpoins(ipoin,j)) then
                                go to 109
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),3)
109                     continue
  
 
 
                         do j=1,k
      if(lnods(lelemp(ipoin,i),5).eq.lpoins(ipoin,j)) then
                                go to 113
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)
113                     continue
                            
                        do j=1,k    
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)  
                                        
                        go to 200
                        

                    else if(inode.eq.5) then
                    
                    
                        do j=1,k
      if(lnods(lelemp(ipoin,i),2).eq.lpoins(ipoin,j)) then
                                go to 104
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),2)
104                     continue

                        do j=1,k
      if(lnods(lelemp(ipoin,i),1).eq.lpoins(ipoin,j)) then
                                go to 110
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),1)
110                    continue

                        do j=1,k
      if(lnods(lelemp(ipoin,i),3).eq.lpoins(ipoin,j)) then
                                go to 114
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),3)
114                    continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),4).eq.lpoins(ipoin,j)) then
                                go to 115
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
115                    continue

                            
                        do j=1,k    
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                 go to 200
                            end if
                        end do

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)                            
                        go to 200
                        

                    else if(inode.eq.6) then

                        do j=1,k
      if(lnods(lelemp(ipoin,i),1).eq.lpoins(ipoin,j)) then
                                go to 105
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),1)
105                     continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),2).eq.lpoins(ipoin,j)) then
                                go to 111
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),2)
111                     continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),3).eq.lpoins(ipoin,j)) then
                                go to 116
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),3)
116                     continue


                        do j=1,k
      if(lnods(lelemp(ipoin,i),4).eq.lpoins(ipoin,j)) then
                                go to 117
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
117                     continue

                            
                        do j=1,k    
      if(lnods(lelemp(ipoin,i),5).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do    
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)                            
                        go to 200
                        
                    end if
                end if

            end do
200         continue


        end do
        nlpoins(ipoin)=k
        enddo
        
        
        do ipoin=1,npoin
          write(234,*) ipoin,nlpoins(ipoin),
     *       (lpoins(ipoin,k),k=1,nlpoins(ipoin))
        
        enddo
      
      close(unit=21)
      close(unit=234)
      
      return
      end
      
      



c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine creazione_lelemp_lpoins(nelem,npoin,nnode,ndime,
     *           nelemp_max,coord,lnods,nelemp,lelemp,
     *           nlpoins,lpoins)
      
      implicit none

      integer*4 npoin,ipoin,i,j,ielem,nelem,idime,ndime,nnode,
     *            inode,k,nelemp_max
      
      Integer*4 lnods(nelem,nnode),numel,nelemp(npoin)
        
      Double precision coord(npoin,ndime)

      integer*4 lelemp(npoin,nelemp_max),lpoins(npoin,nelemp_max),
     *              nlpoins(npoin)
        
        
      character*12 testo

      open(unit=21,file='Corrispondenze_Nodi2Elementi')

      open(unit=234,file='Connessioni_Dirette_Nodi')
  
      lelemp=0
      lpoins=0
 
          
      do ipoin=1,npoin

        k=0
        do ielem=1,nelem
            
          do inode=1,nnode
            if(ipoin.eq.lnods(ielem,inode)) then
              k=k+1
              lelemp(ipoin,k)=ielem
            endif 
          
          enddo
        
          nelemp(ipoin)=k        
        enddo
        
        
1     write(21,*) ipoin,nelemp(ipoin),
     *       (lelemp(ipoin,i),i=1,nelemp(ipoin))
      
      enddo 

        do ipoin=1,npoin
          k=0
          do i=1,nelemp(ipoin)
            do inode = 1,nnode            
              if(ipoin.eq.lnods(lelemp(ipoin,i),inode)) then
                if(inode.eq.1) then
                    
                    
        do j=1,k
              if(lnods(lelemp(ipoin,i),4)
     *                  .eq.lpoins(ipoin,j)) then
                    go to 100
                  end if
                end do
                k=k+1
                lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
100             continue

                        do j=1,k
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do


                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)                            
                        go to 200
                        

                    else if(inode.eq.2) then
                        do j=1,k
      if(lnods(lelemp(ipoin,i),4).eq.lpoins(ipoin,j)) then
                                go to 101
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),4)
101                     continue

                        do j=1,k
      if(lnods(lelemp(ipoin,i),5).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)                            
                        go to 200
                    

                    else if(inode.eq.3) then
                        do j=1,k
      if(lnods(lelemp(ipoin,i),5).eq.lpoins(ipoin,j)) then
                                go to 102
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),5)
102                     continue

                        do j=1,k    
      if(lnods(lelemp(ipoin,i),6).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do    
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),6)                            
                        go to 200
                        

                    else if(inode.eq.4) then
                        do j=1,k
      if(lnods(lelemp(ipoin,i),1).eq.lpoins(ipoin,j)) then
                                go to 103
                            end if
                        end do    

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),1)
103                     continue
                            
                        do j=1,k    
      if(lnods(lelemp(ipoin,i),2).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),2)  
                                        
                        go to 200
                        

                    else if(inode.eq.5) then
                        do j=1,k
      if(lnods(lelemp(ipoin,i),2).eq.lpoins(ipoin,j)) then
                                go to 104
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),2)
104                     continue
                            
                        do j=1,k    
      if(lnods(lelemp(ipoin,i),3).eq.lpoins(ipoin,j)) then
                                 go to 200
                            end if
                        end do

                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),3)                            
                        go to 200
                        

                    else if(inode.eq.6) then

                        do j=1,k
      if(lnods(lelemp(ipoin,i),1).eq.lpoins(ipoin,j)) then
                                go to 105
                            end if
                        end do
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),1)
105                     continue
                            
                        do j=1,k    
      if(lnods(lelemp(ipoin,i),3).eq.lpoins(ipoin,j)) then
                                go to 200
                            end if
                        end do    
                        k=k+1
                        lpoins(ipoin,k)=lnods(lelemp(ipoin,i),3)                            
                        go to 200
                        
                    end if
                end if

            end do
200         continue
        end do
        nlpoins(ipoin)=k
        enddo
        
        
        do ipoin=1,npoin
          write(234,*) ipoin,nlpoins(ipoin),
     *       (lpoins(ipoin,k),k=1,nlpoins(ipoin))
        
        enddo
      
      close(unit=21)
      close(unit=234)
      
      return
      end



c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine interpolazione_elementi2nodi(nelem,npoin,nnode,
     *           ndime,coord,lnods,sn,cf,val,valN)

      implicit none
      integer*4 npoin,ipoin,i,ielem,nelem,idime,ndime,nnode,inode,k
          
      Integer*4 lnods(nelem,nnode),tip,numel
        
      Double precision coord(npoin,ndime),sn(npoin),cf(npoin),
     *            val(nelem),valE(nelem,nnode),valN(npoin)        
        
!      print*,'#--Inizio interpolazione'

!      print*,'#--Salvo valori degli elementi nei nodi'
      do ielem=1,nelem
        valE(ielem,1)=val(ielem)
        valE(ielem,2)=val(ielem)
        valE(ielem,3)=val(ielem)
        valE(ielem,4)=val(ielem)
        valE(ielem,5)=val(ielem)
        valE(ielem,6)=val(ielem)
      enddo

!      print*,'#--Somma di tutti i contributi ai nodi'
      do ipoin=1,npoin
        valN(ipoin)=0.0D0
        k=0
        do ielem=1,nelem
          do inode=1,nnode
            if(ipoin.eq.lnods(ielem,inode)) then          
              k=k+1
              valN(ipoin)=valN(ipoin)+valE(ielem,inode)
c               if(ipoin.eq.825) print *,ipoin,ielem,inode,k
            endif
          enddo
        enddo
        valN(ipoin)=valN(ipoin)/k
        

        write(22,*) ipoin,coord(ipoin,1)*700.D0,
     *                    coord(ipoin,2)*700.D0,
     *                    valN(ipoin) 

      enddo


!      print*,'#--Chiudo file'
      close(unit=22)

!      print*,'#--Fine interpolazione'
      return
      end


c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine dynamic_sn_istantaneo(nelem,npoin,V1,V2,Uo,Lo,
     *           nsn_start,nsn_iter,cf_fix,  
     *           nelemp_max,nnode,ndime,
     *           nelemp,lelemp,nlpoins,lpoins,
     *           coord,lnods,Xe,Ye,time,
     *           Lx_min,Lx_max,Ly_min,Ly_max,xin,yin,
     *           I2,it,tip,nsn,snf,cfl,sn,cf,NnewSN)
      
      Implicit none 
      
      integer*4 nelem,npoin,nnode,ndime,nelemp_max,nsn,ipoin,i,
     *          ielem,idime,numel,nsn_pot,it,nsn_start,nsn_iter,
     *          NnewSN

      Double precision Lx_min,Lx_max,Ly_min,Ly_max,yieldp,yield,Uo,
     *                 Lo,Lx1,Lx2,Ly1,Ly2,time,xin,yin,cf_fix

      parameter(Lx1=-100.D0/700.D0,
     *          Lx2=100.D0/700.D0,
     *          Ly1=-100.D0/700.D0,
     *          Ly2= 0.D0/700.D0,
c     *          yield=450.0D06,!su effective stress
     *          yieldp=1.0D-13 !su eff stress
     *          )     

      integer*4 lelemp(npoin,nelemp_max),lpoins(npoin,nelemp_max),
     *          nlpoins(npoin),lnods(nelem,nnode),nelemp(npoin),
     *          sn_pot(npoin),tip(nelem),snf(npoin)  

      Double precision coord(npoin,ndime),sn(npoin),cf(npoin),
     *              I2(nelem),I2n(npoin),n,x,y,v(npoin),u(npoin),
     *              Xe(nelem),Ye(nelem),V1(npoin),V2(npoin),
     *              I2_pot(npoin),cfl(nelem,nnode)

      yield=yieldp!/((Uo*1.D-02/(365*24*60*60))/Lo)

c      print*,''
c      print*,"#########################################"
       print*,"### Inizio calcolo dynamic slip nodes ###"
c      print*,"#########################################"

      NnewSN=0
      cfl=1.D0

c --- interpolo la grandezza di interesse da elem a punti
!      print*,"# Interpolo da elementi a nodi"
      call interpolazione_elementi2nodi(nelem,npoin,nnode,
     *           ndime,coord,lnods,sn,cf,I2,I2n)

      print *, "I2 massimo: ",MAXVAL(I2N,npoin)
c --- attivazione potenziali
!      print*,"# Preparo attivazione sn potenziali"
      sn_pot=0 !vettore contenente il N dei nodi che potrebbero diventare slip node
      nsn_pot=0
      I2_pot=0
 
c     -- escludo quelli delle faglie precedenti  
!      print*,"# Escludo i nodi delle faglie precedenti"
      do i=1,nsn
        ipoin=snf(i)
        sn_pot(ipoin)=2
      end do     

      

c    --- considero solo i rimanenti e all'interno di un rettangolo
!      print*,"# Cerco sn potenziali all'interno di un rettangolo"
      do ipoin=1,npoin
        if(coord(ipoin,1).gt.Lx1.and.coord(ipoin,1).lt.Lx2.
     *   and.coord(ipoin,2).gt.Ly1.and.coord(ipoin,2).lt.Ly2) then
c        print *, ipoin, I2n(ipoin)
          if(I2n(ipoin).ge.yield.and.sn_pot(ipoin).eq.0) then
            sn_pot(ipoin)=1
            I2_pot(ipoin)=I2n(ipoin)
c          print*,ipoin
c          else 
c          print*,0
          end if
        endif
      enddo 

!      print*,'sn potenziali:',nsn_pot
 
      
      do i=1,nsn
        ipoin=snf(i)
 !       print*,ipoin,snf(i)
        do ielem=1,nelemp(ipoin)
          if(ipoin.eq.
     *          lnods(lelemp(ipoin,ielem),4)) then
!            print*,"4"
            sn_pot(lnods(lelemp(ipoin,ielem),3))=0
            I2_pot(lnods(lelemp(ipoin,ielem),3))=0
            sn_pot(lnods(lelemp(ipoin,ielem),5))=0
            I2_pot(lnods(lelemp(ipoin,ielem),5))=0
            sn_pot(lnods(lelemp(ipoin,ielem),6))=0
            I2_pot(lnods(lelemp(ipoin,ielem),6))=0
          else if(ipoin.eq.
     *         lnods(lelemp(ipoin,ielem),5)) then
!           print*,"5"
            sn_pot(lnods(lelemp(ipoin,ielem),1))=0
            I2_pot(lnods(lelemp(ipoin,ielem),1))=0
            sn_pot(lnods(lelemp(ipoin,ielem),4))=0
            I2_pot(lnods(lelemp(ipoin,ielem),4))=0
            sn_pot(lnods(lelemp(ipoin,ielem),6))=0
            I2_pot(lnods(lelemp(ipoin,ielem),6))=0
          else if(ipoin.eq.
     *            lnods(lelemp(ipoin,ielem),6)) then
!            print*,"6"
            sn_pot(lnods(lelemp(ipoin,ielem),2))=0
            I2_pot(lnods(lelemp(ipoin,ielem),2))=0
            sn_pot(lnods(lelemp(ipoin,ielem),4))=0
            I2_pot(lnods(lelemp(ipoin,ielem),4))=0
            sn_pot(lnods(lelemp(ipoin,ielem),5))=0
            I2_pot(lnods(lelemp(ipoin,ielem),5))=0
          end if
!         print*,"end ",ipoin
        end do

      end do 

      do ipoin=1,npoin
        if (sn_pot(ipoin).eq.1) then
          nsn_pot=nsn_pot+1
        endif
      end do

      print*,"Potenziali: ",nsn_pot

c --- generazione faglia

      do ipoin=1,npoin
        U(ipoin)=V1(ipoin)*Uo
        V(ipoin)=V2(ipoin)*Uo
      end do

      if(nsn_pot.gt.0)  then
!        print*,"# Inizio calcolo faglie"
!        print*,'SONO IN IF x FAULT - potenziali:',nsn_pot  
        

        call fault_doppia_insta(sn,nsn_pot,sn_pot,I2_pot,
     *                 U,V,nsn_start,nsn_iter,cf_fix,nelem,
     *                 npoin,nnode,ndime,lnods,coord,nelemp_max,
     *                nelemp,lelemp,nlpoins,lpoins,xe,ye,time,
     *                 Lx_min,Lx_max,Ly_min,Ly_max,7.D05,it,
     *                 xin,yin,nsn,snf,tip,cfl,cf,NnewSN)
!        print*,'fault ok'

      end if

      nsn=0
      snf=0
      cfl=1.D0
      cf=1.0D0
      tip=1

      return
      end
 
c-----------------------------------------      
      subroutine dynamic_sn_istantaneo_nv(nelem,npoin,V1,V2,Uo,Lo,
     *           nsn_start,nsn_iter,cf_fix,  
     *           nelemp_max,nnode,ndime,
     *           nelemp,lelemp,nlpoins,lpoins,
     *           coord,lnods,Xe,Ye,time,
     *           Lx_min,Lx_max,Ly_min,Ly_max,xin,yin,
     *           I2,it,tip,nsn,snf,cfl,sn,cf,NnewSN)
      
      Implicit none
      
      integer*4 nelem,npoin,nnode,ndime,nelemp_max,nsn,ipoin,i,
     *          ielem,idime,numel,nsn_pot,it,nsn_start,nsn_iter,
     *          NnewSN

      Double precision Lx_min,Lx_max,Ly_min,Ly_max,yieldp,yield,Uo,
     *                 Lo,Lx1,Lx2,Ly1,Ly2,time,xin,yin,cf_fix

      parameter(Lx1=-90.D0/700.D0,
     *          Lx2=+100.D0/700.D0,
     *          Ly1=-100.D0/700.D0,
     *          Ly2= 0.D0/700.D0,
c     *          yield=450.0D06,!su effective stress
     *          yieldp=1.0D-13 !su max shear strain rate
c     *          yieldp=0.1D-13 !su max shear stress
     *          )     

      integer*4 lelemp(npoin,nelemp_max),lpoins(npoin,nelemp_max),
     *          nlpoins(npoin),lnods(nelem,nnode),nelemp(npoin),
     *          sn_pot(npoin),tip(nelem),snf(npoin)  

      Double precision coord(npoin,ndime),sn(npoin),cf(npoin),
     *              I2(nelem),I2n(npoin),n,x,y,v(npoin),u(npoin),
     *              Xe(nelem),Ye(nelem),V1(npoin),V2(npoin),
     *              I2_pot(npoin),cfl(nelem,nnode)

      yield=yieldp!/((Uo*1.D-02/(365*24*60*60))/Lo)

c      print*,''
c      print*,"#########################################"
       print*,"### Inizio calcolo dynamic slip nodes ###"
c      print*,"#########################################"

      NnewSN=0
      cfl=1.D0

c --- interpolo la grandezza di interesse da elem a punti
!      print*,"# Interpolo da elementi a nodi"
      call interpolazione_elementi2nodi(nelem,npoin,nnode,
     *           ndime,coord,lnods,sn,cf,I2,I2n)

      print *, "I2 massimo: ",MAXVAL(I2N,npoin)
c --- attivazione potenziali
!      print*,"# Preparo attivazione sn potenziali"
      sn_pot=0 !vettore contenente il N dei nodi che potrebbero diventare slip node
      nsn_pot=0
      I2_pot=0.D0
 
c     -- escludo quelli delle faglie precedenti  
!      print*,"# Escludo i nodi delle faglie precedenti"
      do i=1,nsn
        ipoin=snf(i)
        sn_pot(ipoin)=2
      end do     

      

c    --- considero solo i rimanenti e all'interno di un rettangolo
!      print*,"# Cerco sn potenziali all'interno di un rettangolo"
      do ipoin=1,npoin
        if(coord(ipoin,1).gt.Lx1.and.coord(ipoin,1).lt.Lx2.
     *   and.coord(ipoin,2).gt.Ly1.and.coord(ipoin,2).lt.Ly2) then
c        print *, ipoin, I2n(ipoin)
          if(I2n(ipoin).ge.yield.and.sn_pot(ipoin).eq.0) then
            sn_pot(ipoin)=1
            I2_pot(ipoin)=I2n(ipoin)
c          print*,ipoin
c          else 
c          print*,0
          end if
        endif
      enddo 

!      print*,'sn potenziali:',nsn_pot
 
      
c      do i=1,nsn
c        ipoin=snf(i)
 !       print*,ipoin,snf(i)
c        do ielem=1,nelemp(ipoin)
c          if(ipoin.eq.
c     *          lnods(lelemp(ipoin,ielem),4)) then
!            print*,"4"
c            sn_pot(lnods(lelemp(ipoin,ielem),3))=0
c            I2_pot(lnods(lelemp(ipoin,ielem),3))=0
c            sn_pot(lnods(lelemp(ipoin,ielem),5))=0
c            I2_pot(lnods(lelemp(ipoin,ielem),5))=0
c            sn_pot(lnods(lelemp(ipoin,ielem),6))=0
c            I2_pot(lnods(lelemp(ipoin,ielem),6))=0
c          else if(ipoin.eq.
c     *         lnods(lelemp(ipoin,ielem),5)) then
!           print*,"5"
c            sn_pot(lnods(lelemp(ipoin,ielem),1))=0
c            I2_pot(lnods(lelemp(ipoin,ielem),1))=0
c            sn_pot(lnods(lelemp(ipoin,ielem),4))=0
c            I2_pot(lnods(lelemp(ipoin,ielem),4))=0
c            sn_pot(lnods(lelemp(ipoin,ielem),6))=0
c            I2_pot(lnods(lelemp(ipoin,ielem),6))=0
c          else if(ipoin.eq.
c     *            lnods(lelemp(ipoin,ielem),6)) then
!            print*,"6"
c            sn_pot(lnods(lelemp(ipoin,ielem),2))=0
c            I2_pot(lnods(lelemp(ipoin,ielem),2))=0
c            sn_pot(lnods(lelemp(ipoin,ielem),4))=0
c            I2_pot(lnods(lelemp(ipoin,ielem),4))=0
c            sn_pot(lnods(lelemp(ipoin,ielem),5))=0
c            I2_pot(lnods(lelemp(ipoin,ielem),5))=0
c          end if
!         print*,"end ",ipoin
c        end do
c
c      end do 

      do ipoin=1,npoin
        if (sn_pot(ipoin).eq.1) then
          nsn_pot=nsn_pot+1
        endif
      end do

      print*,"Potenziali: ",nsn_pot

c --- generazione faglia

      do ipoin=1,npoin
        U(ipoin)=V1(ipoin)*Uo
        V(ipoin)=V2(ipoin)*Uo
      end do

      if(nsn_pot.gt.0)  then
!        print*,"# Inizio calcolo faglie"
!        print*,'SONO IN IF x FAULT - potenziali:',nsn_pot  
        

        call fault_doppia_insta_nv(sn,nsn_pot,sn_pot,I2_pot,
     *                 U,V,nsn_start,nsn_iter,cf_fix,nelem,
     *                 npoin,nnode,ndime,lnods,coord,nelemp_max,
     *                nelemp,lelemp,nlpoins,lpoins,xe,ye,time,
     *                 Lx_min,Lx_max,Ly_min,Ly_max,7.D05,it,
     *                 xin,yin,nsn,snf,tip,cfl,cf,NnewSN)
!        print*,'fault ok'

      end if

      nsn=0
      snf=0
      cfl=1.D0
      cf=1.0D0
      tip=1

      return
      end
      
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine dynamic_sn_time(nelem,npoin,V1,V2,Uo,Lo,
     *           nsn_start,nsn_iter,cf_fix,  
     *           nelemp_max,nnode,ndime,
     *           nelemp,lelemp,nlpoins,lpoins,
     *           coord,lnods,Xe,Ye,time,
     *           Lx_min,Lx_max,Ly_min,Ly_max,xin,yin,
     *           I2,nsn,snf,it,tip,cfl,sn,cf)
      
      Implicit none
      
      integer*4 nelem,npoin,nnode,ndime,nelemp_max,nsn,ipoin,i,
     *          ielem,idime,numel,nsn_pot,it,nsn_start,nsn_iter

      Double precision Lx_min,Lx_max,Ly_min,Ly_max,yieldp,yield,Uo,
     *                 Lo,Lx1,Lx2,Ly1,Ly2,time,xin,yin,cf_fix

      parameter(Lx1=-100.D0/700.D0,
     *          Lx2=+100.D0/700.D0,
     *          Ly1=-100.D0/700.D0,
     *          Ly2= 0.D0/700.D0,
c     *          yield=450.0D06,!su effective stress
     *          yieldp=0.8D-13 !su eff stress
     *          )     

      integer*4 lelemp(npoin,nelemp_max),lpoins(npoin,nelemp_max),
     *          nlpoins(npoin),lnods(nelem,nnode),nelemp(npoin),
     *          sn_pot(npoin),tip(nelem),snf(npoin)  

      Double precision coord(npoin,ndime),sn(npoin),cf(npoin),
     *              I2(nelem),I2n(npoin),n,x,y,v(npoin),u(npoin),
     *              Xe(nelem),Ye(nelem),V1(npoin),V2(npoin),
     *              I2_pot(npoin),cfl(nelem,nnode)

      yield=yieldp!/((Uo*1.D-02/(365*24*60*60))/Lo)

c      print*,''
c      print*,"#########################################"
       print*,"### Inizio calcolo dynamic slip nodes ###"
c      print*,"#########################################"



c --- interpolo la grandezza di interesse da elem a punti
!      print*,"# Interpolo da elementi a nodi"
      call interpolazione_elementi2nodi(nelem,npoin,nnode,
     *           ndime,coord,lnods,sn,cf,I2,I2n)

      print *, "I2 massimo: ",MAXVAL(I2N,npoin)
c --- attivazione potenziali
!      print*,"# Preparo attivazione sn potenziali"
      sn_pot=0 !vettore contenente il N dei nodi che potrebbero diventare slip node
      nsn_pot=0
      I2_pot=0
 
c     -- escludo quelli delle faglie precedenti  
!      print*,"# Escludo i nodi delle faglie precedenti"
      do i=1,nsn
        ipoin=snf(i)
        sn_pot(ipoin)=2
      end do     

c    --- considero solo i rimanenti e all'interno di un rettangolo
!      print*,"# Cerco sn potenziali all'interno di un rettangolo"
      do ipoin=1,npoin
        if(coord(ipoin,1).gt.Lx1.and.coord(ipoin,1).lt.Lx2.
     *   and.coord(ipoin,2).gt.Ly1.and.coord(ipoin,2).lt.Ly2) then
c        print *, ipoin, I2n(ipoin)
          if(I2n(ipoin).ge.yield.and.sn_pot(ipoin).eq.0) then
            sn_pot(ipoin)=1
            I2_pot(ipoin)=I2n(ipoin)
c          print*,ipoin
c          else 
c          print*,0
          end if
        endif
      enddo 

!      print*,'sn potenziali:',nsn_pot
 
      
      do i=1,nsn
        ipoin=snf(i)
 !       print*,ipoin,snf(i)
        do ielem=1,nelemp(ipoin)
          if(ipoin.eq.
     *          lnods(lelemp(ipoin,ielem),4)) then
!            print*,"4"
            sn_pot(lnods(lelemp(ipoin,ielem),3))=0
            I2_pot(lnods(lelemp(ipoin,ielem),3))=0
            sn_pot(lnods(lelemp(ipoin,ielem),5))=0
            I2_pot(lnods(lelemp(ipoin,ielem),5))=0
            sn_pot(lnods(lelemp(ipoin,ielem),6))=0
            I2_pot(lnods(lelemp(ipoin,ielem),6))=0
          else if(ipoin.eq.
     *         lnods(lelemp(ipoin,ielem),5)) then
!           print*,"5"
            sn_pot(lnods(lelemp(ipoin,ielem),1))=0
            I2_pot(lnods(lelemp(ipoin,ielem),1))=0
            sn_pot(lnods(lelemp(ipoin,ielem),4))=0
            I2_pot(lnods(lelemp(ipoin,ielem),4))=0
            sn_pot(lnods(lelemp(ipoin,ielem),6))=0
            I2_pot(lnods(lelemp(ipoin,ielem),6))=0
          else if(ipoin.eq.
     *            lnods(lelemp(ipoin,ielem),6)) then
!            print*,"6"
            sn_pot(lnods(lelemp(ipoin,ielem),2))=0
            I2_pot(lnods(lelemp(ipoin,ielem),2))=0
            sn_pot(lnods(lelemp(ipoin,ielem),4))=0
            I2_pot(lnods(lelemp(ipoin,ielem),4))=0
            sn_pot(lnods(lelemp(ipoin,ielem),5))=0
            I2_pot(lnods(lelemp(ipoin,ielem),5))=0
          end if
!         print*,"end ",ipoin
        end do

      end do 

      do ipoin=1,npoin
        if (sn_pot(ipoin).eq.1) then
          nsn_pot=nsn_pot+1
        endif
      end do

      print*,"Potenziali: ",nsn_pot

c --- generazione faglia

      do ipoin=1,npoin
        U(ipoin)=V1(ipoin)*Uo
        V(ipoin)=V2(ipoin)*Uo
      end do

      if(nsn_pot.gt.0)  then
!        print*,"# Inizio calcolo faglie"
!        print*,'SONO IN IF x FAULT - potenziali:',nsn_pot  
        

        call fault_doppia(sn,nsn_pot,sn_pot,I2_pot,
     *                 U,V,nsn_start,nsn_iter,cf_fix,nelem,
     *                 npoin,nnode,ndime,lnods,coord,nelemp_max,
     *                nelemp,lelemp,nlpoins,lpoins,xe,ye,time,
     *                 Lx_min,Lx_max,Ly_min,Ly_max,it,
     *                 xin,yin,nsn,snf,tip,cfl,cf)
!        print*,'fault ok'

      end if

      return
      end

c===================================================================

      subroutine fault_doppia(sn,nsn_pot,sn_pot,I2_pot,
     *                 U,V,nsn_start,nsn_iter,cf_fix,nelem,
     *                 npoin,nnode,ndime,lnods,coord,nelemp_max,
     *                 nelemp,lelemp,nlpoins,lpoins,xe,ye,time,
     *                 Lx_min,Lx_max,Ly_min,Ly_max,it,
     *                 xin,yin,nsn,snf,tip,cfl,cf)
       
      implicit none
       
      Integer*4   npoin,nelem,ipoin,ielem,idime,ndime,inode,
     *            nelemp_max,n,nsn_start,nsn_iter,it,
     *            nnode,lnods(nelem,nnode),nsn_pot,sn_pot(npoin),
     *            nelemp(npoin),lelemp(npoin,nelemp_max),ielemp,
     *            nlpoins(npoin),lpoins(npoin,nelemp_max),j,
     *            i,k,ss,np_pol,nsn,cont,
     *            tip(nelem),divsx,divdx,start,new_nsn
     
      Integer*4 snf(npoin),snf2(nsn_iter+1)
       
      Double precision Lx_min,Lx_max,Ly_min,Ly_max,vsx,vdx,vtot,
     *                 xin,yin,xfin,yfin,cf_fix,I2_0(nelemp_max)

      Double precision coord(npoin,ndime),I2_pot(npoin),I2_max,
     *                 v(npoin),u(npoin),xe(nelem),ye(nelem),
     *                 c(npoin),time,cfl(nelem,nnode),sn(npoin),
     *                 cf(npoin)

      Double precision,Dimension(:),allocatable:: x_pol, y_pol

      Character(len=100) :: filename
      
    
      I2_max=MAXVAL(I2_pot,nsn_pot)
      start=nsn+1
      new_nsn=0
      
      if(it.eq.1) then 
      
c      print*,' '
c      print*, '#######################'
c      print *,"------- NUCLEO --------"
        do ipoin=1,npoin       
         if(I2_pot(ipoin).eq.I2_max.and.
     *     sn_pot(ipoin).eq.1) then   
                          
            snf(start)=ipoin    
            sn_pot(snf(start))=2 
            sn(snf(start))=1 
            cf(snf(start))=cf_fix
            new_nsn=new_nsn+1      
            write(28,*) snf(start),coord(snf(start),1)*700.,
     *                         coord(snf(start),2)*700.,time
     
            write(30,*) snf(start),coord(snf(start),1)*700.,
     *                         coord(snf(start),2)*700.,time
     
c            print*, snf(start),coord(snf(start),1)*700.,
c     *                         coord(snf(start),2)*700.,time,
c     *                         I2_pot(snf(start))
            do ielem=1,nelemp(ipoin)
              if(ipoin.eq.
     *            lnods(lelemp(ipoin,ielem),4)) then
                sn_pot(lnods(lelemp(ipoin,ielem),3))=0
                sn_pot(lnods(lelemp(ipoin,ielem),5))=0
                sn_pot(lnods(lelemp(ipoin,ielem),6))=0
              else if(ipoin.eq.
     *            lnods(lelemp(ipoin,ielem),5)) then
                sn_pot(lnods(lelemp(ipoin,ielem),1))=0
                sn_pot(lnods(lelemp(ipoin,ielem),4))=0
                sn_pot(lnods(lelemp(ipoin,ielem),6))=0
              else if(ipoin.eq.
     *            lnods(lelemp(ipoin,ielem),6)) then
                sn_pot(lnods(lelemp(ipoin,ielem),2))=0
                sn_pot(lnods(lelemp(ipoin,ielem),4))=0
                sn_pot(lnods(lelemp(ipoin,ielem),5))=0
              end if
            end do
c          print*, start,snf(start),coord(snf(start),1)*700.,
c     *                         coord(snf(start),2)*700.
            exit              
          endif 
        enddo      

      else 
        start=start-1
      endif

c ---- propagazione faglia in coda
 
c      print *,"------- PRIMO RAMO FAGLIA --------"
c      print*,"start:",start
      do k=start,npoin
        if (it.eq.1.and.new_nsn.eq.nsn_start) then
          exit
        else if (it.ne.1.and.new_nsn.eq.nsn_iter) then
          exit
        end if
    
      
c        print*,"K=",k
        I2_0=0.0D0
        ss=0
        n=nlpoins(snf(k))
c        print*,"snf(k)=",snf(k)  
c        print*,"n=",n
        do i=1,n   
          if(sn_pot(lpoins(snf(k),i)).eq.1) then
            ss=1
            I2_0(i)=I2_pot(lpoins(snf(k),i))    
          endif
        end do

        if(ss.eq.0) then
c          print*,"ss=0"
          exit
        else
          I2_max=MAXVAL(I2_0,nelemp_max)
          do i=1,n
            if(I2_0(i).eq.I2_max .and.
     *           sn_pot(lpoins(snf(k),i)).eq.1) then                    
              snf(k+1)=lpoins(snf(k),i)
              sn(snf(k+1))=1
              cf(snf(k+1))=cf_fix
              sn_pot(snf(k+1))=2 
              new_nsn=new_nsn+1
c              print*, k+1,snf(k+1),coord(snf(k+1),1)*700.,
c     *                         coord(snf(k+1),2)*700.,
c     *                         I2_pot(snf(k+1))
              write(28,*) snf(k+1),coord(snf(k+1),1)*700.,
     *                         coord(snf(k+1),2)*700.,time



              do ielem=1,nelemp(snf(k+1))

c                print*,ielem,snf(k+1),
c     *           (lnods(lelemp(snf(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf(k+1),ielem),inode)),
c     *           inode=1,6)
                
                if(snf(k+1).eq.
     *           lnods(lelemp(snf(k+1),ielem),4)) then
                  sn_pot(lnods(lelemp(snf(k+1),ielem),3))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),5))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),6))=0
                else if(snf(k+1).eq.
     *             lnods(lelemp(snf(k+1),ielem),5)) then
                  sn_pot(lnods(lelemp(snf(k+1),ielem),1))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),4))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),6))=0
                else if(snf(k+1).eq.
     *             lnods(lelemp(snf(k+1),ielem),6)) then
                  sn_pot(lnods(lelemp(snf(k+1),ielem),2))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),4))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),5))=0
                end if
c                print*,ielem,snf(k+1),
c     *           (lnods(lelemp(snf(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf(k+1),ielem),inode)),
c     *           inode=1,6)
              end do

              
              exit              
            endif   
          enddo    
        endif


      enddo    


      nsn=nsn+new_nsn



c ---- propagazione faglia in testa
      new_nsn=0
      start=snf(1)
      snf2(1)=start

c      print *,"------- SECONDO RAMO FAGLIA --------"
c      print*,"start:",start

       write(28,1000)
1000   format(1H>)

      do k=1,nsn_iter
        if (it.eq.1) then
c          print*,"it:",it
          exit
       end if
c        print*,"K=",k
        I2_0=0.0D0
        ss=0
        n=nlpoins(snf2(k))
c        print*,"snf2(k)=",snf2(k)  
c        print*,"n=",n
        do i=1,n
          if(sn_pot(lpoins(snf2(k),i)).eq.1) then
            ss=1
            I2_0(i)=I2_pot(lpoins(snf2(k),i))    
          endif
        end do

        if(ss.eq.0) then
c          print*,"ss=0"
          exit
        else
          I2_max=MAXVAL(I2_0,nelemp_max)
          do i=1,n
            if(I2_0(i).eq.I2_max .and.
     *           sn_pot(lpoins(snf2(k),i)).eq.1) then                    
              snf2(k+1)=lpoins(snf2(k),i)
              sn(snf2(k+1))=1
              cf(snf2(k+1))=cf_fix
              sn_pot(snf2(k+1))=2 
              new_nsn=new_nsn+1
c              print*, k+1,snf2(k+1),coord(snf2(k+1),1)*700.,
c     *                         coord(snf2(k+1),2)*700.,
c     *                         I2_pot(snf2(k+1))
              write(29,*) snf2(k+1),coord(snf2(k+1),1)*700.,
     *                         coord(snf2(k+1),2)*700.,time

              do ielem=1,nelemp(snf2(k+1))

c                print*,ielem,snf2(k+1),
c     *           (lnods(lelemp(snf2(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf2(k+1),ielem),inode)),
c     *           inode=1,6)
                
                if(snf2(k+1).eq.
     *           lnods(lelemp(snf2(k+1),ielem),4)) then
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),3))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),5))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),6))=0
                else if(snf2(k+1).eq.
     *             lnods(lelemp(snf2(k+1),ielem),5)) then
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),1))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),4))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),6))=0
                else if(snf2(k+1).eq.
     *             lnods(lelemp(snf2(k+1),ielem),6)) then
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),2))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),4))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),5))=0
                end if
c                print*,ielem,snf2(k+1),
c     *           (lnods(lelemp(snf2(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf2(k+1),ielem),inode)),
c     *           inode=1,6)
              end do

              
              exit              
            endif   
          enddo    
        endif


      enddo   

c ----- aggiungo il ramo ribaltato in testa a snf (shift) -----

      do i=1,nsn
        j=nsn+1-i
        snf(j+new_nsn)=snf(j)
      end do

      do i=1,new_nsn
        j=new_nsn+2-i
        snf(i)=snf2(j)
      end do
c      print *,"chiuso"

      nsn=nsn+new_nsn

c ----- creazione del poligono

      !print*,' '
      !print*,'------ CREAZIONE POLIGONO ------'
      !print*,"nsn,new_nsn:",nsn,new_nsn
      !print*,' '
      
      
      np_pol=nsn+5

      allocate(x_pol(np_pol))
      allocate(y_pol(np_pol))
      
c      if(it.eq.1) then
        xin=(3*coord(snf(1),1)-2*coord(snf(2),1))*700.0D0
        yin=(3*coord(snf(1),2)-2*coord(snf(2),2))*700.0D0
c      end if

      xfin=(3*coord(snf(nsn),1)-2*coord(snf(nsn-1),1))*700.0D0
      yfin=(3*coord(snf(nsn),2)-2*coord(snf(nsn-1),2))*700.0D0

      x_pol(1)=xin
      y_pol(1)=yin

      do i=1,nsn
        x_pol(i+1)=coord(snf(i),1)*700.
        y_pol(i+1)=coord(snf(i),2)*700.
      end do

      x_pol(nsn+2)=xfin
      y_pol(nsn+2)=yfin
      
      if (xin.ge.xfin.and.
     *    yin.ge.yfin) then
        x_pol(nsn+3)=Lx_min
        y_pol(nsn+3)=Ly_min
        x_pol(nsn+4)=Lx_max
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_max
        y_pol(nsn+5)=Ly_max
        
       else if (xin.ge.xfin.and.
     *          yin.lt.yfin) then
        x_pol(nsn+3)=Lx_min
        y_pol(nsn+3)=Ly_max
        x_pol(nsn+4)=Lx_min
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_max
        y_pol(nsn+5)=Ly_min

      else if (xin.lt.xfin.and.
     *         yin.ge.yfin) then
        x_pol(nsn+3)=Lx_max
        y_pol(nsn+3)=Ly_min
        x_pol(nsn+4)=Lx_min
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_min
        y_pol(nsn+5)=Ly_max

      else if (xin.lt.xfin.and.
     *         yin.lt.yfin) then
        x_pol(nsn+3)=Lx_max
        y_pol(nsn+3)=Ly_max
        x_pol(nsn+4)=Lx_max
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_min
        y_pol(nsn+5)=Ly_min   
      end if

      do i=1,np_pol
        write(34,*)x_pol(i),y_pol(i),time
c        print*,x_pol(i),y_pol(i),time
      end do


c ------ definizione destro-sinistro



      !print*,' '
      !print *,'------- DEFINISCO DX-SX -------'
      tip=0
      do i=1,nsn
        do ielemp=1,nelemp(snf(i))
          ielem=lelemp(snf(i),ielemp)
          call in_out (ielem,xe(ielem),ye(ielem),np_pol,x_pol,y_pol,
     *                 tip(ielem))
          if (tip(ielem).eq.0) tip(ielem)=-1
             write(29,*)  ielem,xe(ielem),ye(ielem),tip(ielem),
     *        snf(i),time
c             print *, i, snf(i), ielem, xe(ielem), ye(ielem), tip(ielem)
        end do
      end do

      deallocate(x_pol)
      deallocate(y_pol)

c      close(unit=29)



c ------ coupling

      do ipoin=1,npoin
        c(ipoin)=0.0D0
      end do
      !      print*,' '
      !print*,"----- Calcolo del coupling -----"
c      write(30,*)"CALCOLO DEL COUPLING"
      !print*,"i    ipoin     vsx      vdx          c"     
c      write(30,*)"i    ipoin     vsx      vdx          c"

      do i=1,nsn
c        print*,i
        vsx=0
        vdx=0
        divsx=0
        divdx=0
        do ielemp=1,nelemp(snf(i))
          ielem=lelemp(snf(i),ielemp)
          do inode=1,nnode
            ipoin=lnods(ielem,inode)
            vtot=sqrt(U(ipoin)**2+V(ipoin)**2)
            if(tip(ielem).eq.-1)then
              vsx=vsx+vtot
              divsx=divsx+1
              cfl(ielem,inode)=1.0D0
            else if(tip(ielem).eq.1)then
              vdx=vdx+vtot
              divdx=divdx+1
              cfl(ielem,inode)=cf_fix
            end if
          end do
        end do

      !        print*,divsx,divdx

        if(divsx.ne.0.and.divdx.ne.0) then
          vsx=vsx/divsx
          vdx=vdx/divdx
          c(snf(i))=vdx/vsx
        else 
          c(snf(i))=0.5
        endif 


      !        print*,i,snf(i),vsx,vdx,c(snf(i))
c        write(30,*) i,snf(i),vsx,vdx,c(snf(i))
      end do
      
c      close(unit=28) 

c      close(unit=29)

c      close(unit=30)     
      
      !print*,' '
      !print *, "----- FINE CICLO -----"

      print*, "nsn: ", nsn
c      do i=1,nsn
c        print*,snf(i)
c      end do

c      print*,"=============="
c      print*,""

      !print*,' '
      !print*,'-------------------------'
      !print*,' '
      !print*,' '


c      do j=1,nsn
c        print*,j,nsn,snf(j)
c      end do

c      print *,'uscita'

      return
      end


c=================================

      subroutine fault_doppia_insta_nv(sn,nsn_pot,sn_pot,I2_pot,
     *                 U,V,nsn_start,nsn_iter,cf_fix,nelem,
     *                 npoin,nnode,ndime,lnods,coord,nelemp_max,
     *                 nelemp,lelemp,nlpoins,lpoins,xe,ye,time,
     *                 Lx_min,Lx_max,Ly_min,Ly_max,Lo,it,
     *                 xin,yin,nsn,snf,tip,cfl,cf,w)
       
       
      implicit none
       
      Integer*4   npoin,nelem,ipoin,ielem,idime,ndime,inode,
     *            nelemp_max,n,nsn_start,nsn_iter,it,
     *            nnode,lnods(nelem,nnode),nsn_pot,sn_pot(npoin),
     *            nelemp(npoin),lelemp(npoin,nelemp_max),ielemp,
     *            nlpoins(npoin),lpoins(npoin,npoin),j,
     *            i,k,ss,np_pol,nsn,cont,
     *            tip(nelem),divsx,divdx,start,new_nsn
     
      Integer*4 snf(npoin),snf2(nsn_iter+1),w
       
      Double precision Lx_min,Lx_max,Ly_min,Ly_max,vsx,vdx,vtot,
     *                 xin,yin,xfin,yfin,cf_fix,vsxp,vdxp,Lo

      Double precision coord(npoin,ndime),I2_pot(npoin),I2_max,
     *                 v(npoin),u(npoin),xe(nelem),ye(nelem),
     *                 c(npoin),time,cfl(nelem,nnode),sn(npoin),
     *                 I2_0(npoin),cf(npoin),tetaSlip(npoin)

      Double precision,Dimension(:),allocatable:: x_pol, y_pol

      Character(len=100) :: filename
      
      I2_max=MAXVAL(I2_pot,nsn_pot)
      start=nsn+1
      new_nsn=0
      
!      if(it.eq.1) then 
      if(nsn.eq.0) then 
      
c      !print*,' '
c      !print*, '#######################'
      print*,"------- NUCLEO --------"
c         print*,I2_max
        w=0
        do ipoin=1,npoin       
         !!print*,ipoin,sn_pot(ipoin),'NUCLEO'
         if(I2_pot(ipoin).eq.I2_max.and.
     *     sn_pot(ipoin).eq.1) then   
            snf(start)=ipoin    
            sn_pot(snf(start))=2 
            sn(snf(start))=1 
            cf(snf(start))=cf_fix
            new_nsn=new_nsn+1 
            w=w+1     
           print*,ipoin,I2_max
           
            write(28,*) snf(start),coord(snf(start),1)*Lo/1.D03,
     *                    coord(snf(start),2)*Lo/1.D03,it,w
            write(30,*) snf(start),coord(snf(start),1)*Lo/1.D03,
     *                    coord(snf(start),2)*Lo/1.D03,it,w
     
            !print*, snf(start),coord(snf(start),1)*Lo/1.D03,
c     *                         coord(snf(start),2)*Lo/1.D03,time,
c     *                         I2_pot(snf(start))
c            do ielem=1,nelemp(ipoin)
c              if(ipoin.eq.
c     *            lnods(lelemp(ipoin,ielem),4)) then
c                sn_pot(lnods(lelemp(ipoin,ielem),3))=0
c                sn_pot(lnods(lelemp(ipoin,ielem),5))=0
c                sn_pot(lnods(lelemp(ipoin,ielem),6))=0
c              else if(ipoin.eq.
c     *            lnods(lelemp(ipoin,ielem),5)) then
c                sn_pot(lnods(lelemp(ipoin,ielem),1))=0
c                sn_pot(lnods(lelemp(ipoin,ielem),4))=0
c                sn_pot(lnods(lelemp(ipoin,ielem),6))=0
c              else if(ipoin.eq.
c     *            lnods(lelemp(ipoin,ielem),6)) then
c                sn_pot(lnods(lelemp(ipoin,ielem),2))=0
c                sn_pot(lnods(lelemp(ipoin,ielem),4))=0
c                sn_pot(lnods(lelemp(ipoin,ielem),5))=0
c              end if
c            end do
          !print*, 'start',snf(start),coord(snf(start),1)*Lo/1.D03,
c     *                         coord(snf(start),2)*Lo/1.D03
            exit              
          endif 
        enddo      

      else 
        start=start-1
      endif

c ---- propagazione faglia in coda
 
      print*,"------- PRIMO RAMO FAGLIA --------"
      !print*,"start:",start
      do k=start,npoin
!        if (it.eq.1.and.new_nsn.eq.nsn_start) then
        if (nsn.eq.0.and.new_nsn.eq.nsn_start) then
          exit
!        else if (it.ne.1.and.new_nsn.eq.nsn_iter) then
         else if (nsn.eq.0.and.new_nsn.eq.nsn_iter) then
         exit
        end if
    
      
c        !print*,"K=",k
        I2_0=0.0D0
        ss=0
        n=nlpoins(snf(k))
        !print*,"snf(k)=",snf(k)  
c        !print*,"n=",n
        do i=1,n  
          if(sn_pot(lpoins(snf(k),i)).eq.1) then
            ss=1
            I2_0(i)=I2_pot(lpoins(snf(k),i))
            print*,"I2_0:",lpoins(snf(k),i),I2_0(i)   
          endif
        end do

        if(ss.eq.0) then
          !print*,"ss=0"
          exit
        else
c          print *,I2_0,n
          
          I2_max=MAXVAL(I2_0,npoin)
          print*, "I2_max=", I2_max
          do i=1,n
            if(I2_0(i).eq.I2_max.and.
     *           sn_pot(lpoins(snf(k),i)).eq.1) then                    
              snf(k+1)=lpoins(snf(k),i)
              sn(snf(k+1))=1
              cf(snf(k+1))=cf_fix
              sn_pot(snf(k+1))=2 
              new_nsn=new_nsn+1
               w=w+1   
              !print*, k+1,snf(k+1),coord(snf(k+1),1)*Lo/1.D03,
c     *                         coord(snf(k+1),2)*Lo/1.D03,
c     *                         I2_pot(snf(k+1))
              write(28,*) snf(k+1),coord(snf(k+1),1)*Lo/1.D03,
     *                         coord(snf(k+1),2)*Lo/1.D03,it,w

c              do ielem=1,nelemp(snf(k+1))

                !print*,ielem,snf(k+1),
c     *           (lnods(lelemp(snf(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf(k+1),ielem),inode)),
c     *           inode=1,6)
                
c                if(snf(k+1).eq.
c     *           lnods(lelemp(snf(k+1),ielem),4)) then
c                  sn_pot(lnods(lelemp(snf(k+1),ielem),3))=0
c                  sn_pot(lnods(lelemp(snf(k+1),ielem),5))=0
c                  sn_pot(lnods(lelemp(snf(k+1),ielem),6))=0
c                else if(snf(k+1).eq.
c     *             lnods(lelemp(snf(k+1),ielem),5)) then
c                  sn_pot(lnods(lelemp(snf(k+1),ielem),1))=0
c                  sn_pot(lnods(lelemp(snf(k+1),ielem),4))=0
c                  sn_pot(lnods(lelemp(snf(k+1),ielem),6))=0
c                else if(snf(k+1).eq.
c     *             lnods(lelemp(snf(k+1),ielem),6)) then
c                  sn_pot(lnods(lelemp(snf(k+1),ielem),2))=0
c                  sn_pot(lnods(lelemp(snf(k+1),ielem),4))=0
c                  sn_pot(lnods(lelemp(snf(k+1),ielem),5))=0
c                end if
c                !print*,ielem,snf(k+1),
c     *           (lnods(lelemp(snf(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf(k+1),ielem),inode)),
c     *           inode=1,6)
c              end do

              
              exit              
            endif   
          enddo    
        endif


      enddo    


      nsn=nsn+new_nsn



c ---- propagazione faglia in testa
      new_nsn=0
      start=snf(1)
      snf2(1)=start

      print*,"------- SECONDO RAMO FAGLIA --------"
      !print*,"start:",start
       
      do k=1,nsn_iter
c        if (it.eq.1) then
c          !print*,"it:",it
c          exit
c       end if

c        !print*,"K=",k
        I2_0=0.0D0
        ss=0
        n=nlpoins(snf2(k))
c        !print*,"snf2(k)=",snf2(k)  
c        !print*,"n=",n
        do i=1,n
          if(sn_pot(lpoins(snf2(k),i)).eq.1) then
            ss=1
            I2_0(i)=I2_pot(lpoins(snf2(k),i))
            
            print*,"I2_0:",lpoins(snf2(k),i),I2_0(i)   
          endif
        end do

        if(ss.eq.0) then
          !print*,"ss=0"
          exit
        else
          I2_max=MAXVAL(I2_0,npoin)
          print *, "I2_max=", I2_max
          do i=1,n
            if(I2_0(i).eq.I2_max .and.
     *           sn_pot(lpoins(snf2(k),i)).eq.1) then                    
              snf2(k+1)=lpoins(snf2(k),i)
              sn(snf2(k+1))=1
              cf(snf2(k+1))=cf_fix
              sn_pot(snf2(k+1))=2 
              new_nsn=new_nsn+1
               w=w+1   
c              !print*, k+1,snf2(k+1),coord(snf2(k+1),1)*Lo/1.D03,
c     *                         coord(snf2(k+1),2)*Lo/1.D03,
c     *                         I2_pot(snf2(k+1))
              write(30,*) snf2(k+1),coord(snf2(k+1),1)*Lo/1.D03,
     *                         coord(snf2(k+1),2)*Lo/1.D03,it,w

c              do ielem=1,nelemp(snf2(k+1))

c                !print*,ielem,snf2(k+1),
c     *           (lnods(lelemp(snf2(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf2(k+1),ielem),inode)),
c     *           inode=1,6)
                
c                if(snf2(k+1).eq.
c     *           lnods(lelemp(snf2(k+1),ielem),4)) then
c                  sn_pot(lnods(lelemp(snf2(k+1),ielem),3))=0
c                  sn_pot(lnods(lelemp(snf2(k+1),ielem),5))=0
c                  sn_pot(lnods(lelemp(snf2(k+1),ielem),6))=0
c                else if(snf2(k+1).eq.
c     *             lnods(lelemp(snf2(k+1),ielem),5)) then
c                  sn_pot(lnods(lelemp(snf2(k+1),ielem),1))=0
c                  sn_pot(lnods(lelemp(snf2(k+1),ielem),4))=0
c                  sn_pot(lnods(lelemp(snf2(k+1),ielem),6))=0
c                else if(snf2(k+1).eq.
c     *             lnods(lelemp(snf2(k+1),ielem),6)) then
c                  sn_pot(lnods(lelemp(snf2(k+1),ielem),2))=0
c                  sn_pot(lnods(lelemp(snf2(k+1),ielem),4))=0
c                  sn_pot(lnods(lelemp(snf2(k+1),ielem),5))=0
c                end if
c                !print*,ielem,snf2(k+1),
c     *           (lnods(lelemp(snf2(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf2(k+1),ielem),inode)),
c     *           inode=1,6)
c              end do

              
              exit              
            endif   
          enddo    
        endif


      enddo   

c ----- aggiungo il ramo ribaltato in testa a snf (shift) -----

      do i=1,nsn
        j=nsn+1-i
        snf(j+new_nsn)=snf(j)
      end do

      do i=1,new_nsn
        j=new_nsn+2-i
        snf(i)=snf2(j)
      end do
c      !!print*,"chiuso"

      nsn=nsn+new_nsn

c ----- creazione del poligono

      !!print*,' '
      !print*,'------ CREAZIONE POLIGONO ------'
      !print*,"nsn,new_nsn:",nsn,new_nsn
      !!print*,' '
            
      np_pol=nsn+5

        !!print*,'nsn',nsn,'np_pol',np_pol
        allocate(x_pol(np_pol))
        allocate(y_pol(np_pol))
      

c      if(it.eq.1) then
        xin=(3*coord(snf(1),1)-2*coord(snf(2),1))*Lo/1.D03
        yin=(3*coord(snf(1),2)-2*coord(snf(2),2))*Lo/1.D03
c      end if

         !!print*,nsn,snf(nsn),'prima'
      xfin=(3*coord(snf(nsn),1)-2*coord(snf(nsn-1),1))*Lo/1.D03
      yfin=(3*coord(snf(nsn),2)-2*coord(snf(nsn-1),2))*Lo/1.D03
          !!print*,nsn,'dopo'

      x_pol(1)=xin
      y_pol(1)=yin

      do i=1,nsn
        x_pol(i+1)=coord(snf(i),1)*Lo/1.D03
        y_pol(i+1)=coord(snf(i),2)*Lo/1.D03
      end do

      x_pol(nsn+2)=xfin
      y_pol(nsn+2)=yfin
      
      if (xin.ge.xfin.and.
     *    yin.ge.yfin) then
        x_pol(nsn+3)=Lx_min
        y_pol(nsn+3)=Ly_min
        x_pol(nsn+4)=Lx_max
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_max
        y_pol(nsn+5)=Ly_max
        
       else if (xin.ge.xfin.and.
     *          yin.lt.yfin) then
        x_pol(nsn+3)=Lx_min
        y_pol(nsn+3)=Ly_max
        x_pol(nsn+4)=Lx_min
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_max
        y_pol(nsn+5)=Ly_min

      else if (xin.lt.xfin.and.
     *         yin.ge.yfin) then
        x_pol(nsn+3)=Lx_max
        y_pol(nsn+3)=Ly_min
        x_pol(nsn+4)=Lx_min
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_min
        y_pol(nsn+5)=Ly_max

      else if (xin.lt.xfin.and.
     *         yin.lt.yfin) then
        x_pol(nsn+3)=Lx_max
        y_pol(nsn+3)=Ly_max
        x_pol(nsn+4)=Lx_max
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_min
        y_pol(nsn+5)=Ly_min   
      end if

      do i=1,np_pol
        write(34,*)x_pol(i),y_pol(i),time
        !print*,x_pol(i),y_pol(i),time
      end do


c ------ definizione destro-sinistro



      !!print*,' '
      !!print*,'------- DEFINISCO DX-SX -------'
      tip=0
      do i=1,nsn
        do ielemp=1,nelemp(snf(i))
          ielem=lelemp(snf(i),ielemp)
          call in_out (ielem,xe(ielem),ye(ielem),np_pol,x_pol,y_pol,
     *                 tip(ielem))
          if (tip(ielem).eq.0) tip(ielem)=-1
             write(29,*)  ielem,xe(ielem),ye(ielem),tip(ielem),
     *        snf(i),time
c             !!print*, i, snf(i), ielem, xe(ielem), ye(ielem), tip(ielem)
        end do
      end do

      deallocate(x_pol)
      deallocate(y_pol)

c      close(unit=29)



       !!print*,'CALCOLO TETASLIP'
       
      do i=2,nsn      
        if(coord(snf(i),1).eq.coord(snf(i-1),1)) then
          tetaSlip(snf(i))=90. !codice x slip verticale  
        else if(coord(snf(i),2).eq.coord(snf(i-1),2)) then      
          tetaSlip(snf(i))=0.   !codice x slip orizzontale   
        else 
         tetaSlip(snf(i))=datan((coord(snf(i),2)-coord(snf(i-1),2))/
     *                   (coord(snf(i),1)-coord(snf(i-1),1)))
        endif
     
c        !print*,i,snf(i),tetaSlip(snf(i))*180.D0/3.14D0
      end do






c ------ coupling

      do ipoin=1,npoin
        c(ipoin)=0.0D0
      end do
      !      !print*,' '
c      !print*,"----- Calcolo del coupling -----"
c      write(30,*)"CALCOLO DEL COUPLING"
      !!print*,"i    ipoin     vsx      vdx          c"     
c      write(30,*)"i    ipoin     vsx      vdx          c"

      do i=1,nsn
c        !print*,i
        vsx=0
        vdx=0
        divsx=0
        divdx=0
        do ielemp=1,nelemp(snf(i))
          ielem=lelemp(snf(i),ielemp)
          do inode=1,nnode
            ipoin=lnods(ielem,inode)
            vtot=sqrt(U(ipoin)**2+V(ipoin)**2)
            if(tip(ielem).eq.-1)then
              vsx=vsx+vtot
              divsx=divsx+1
            else if(tip(ielem).eq.1.and.sn(ipoin).eq.1)then
              vdx=vdx+vtot
              divdx=divdx+1
            end if
          end do
        end do
        
        vsx=vsx/divsx
        vdx=vdx/divdx
         
        vsxp=vsx*dcos(tetaslip(snf(i)))+vsx*dsin(tetaslip(snf(i)))
        
        vdxp=vdx*dcos(tetaslip(snf(i)))+vdx*dsin(tetaslip(snf(i)))
        
        if(vsxp*vdxp.lt.0.) c(snf(i))=-1.
        if(vsxp*vdxp.gt.0.) c(snf(i))=+1.
        
        c(snf(i))=+1.       !da togliere quando aggiusteremo
        
         cf(snf(i))=cf(snf(i))*c(snf(i))


        do ielemp=1,nelemp(snf(i))
          ielem=lelemp(snf(i),ielemp)
          do inode=1,nnode
            ipoin=lnods(ielem,inode)
            if(tip(ielem).eq.-1)then
              cfl(ielem,inode)=1.0D0
            else if(tip(ielem).eq.1.and.sn(ipoin).eq.1)then
              cfl(ielem,inode)=cf(snf(i))
            end if
          end do
        end do








      !        !print*,divsx,divdx

c        if(divsx.ne.0.and.divdx.ne.0) then
c          vsx=vsx/divsx
c          vdx=vdx/divdx
c           c(snf(i))=vdx/vsx
c        else 
c          c(snf(i))=0.5
c        endif 




      !        !print*,i,snf(i),vsx,vdx,c(snf(i))
c        write(30,*) i,snf(i),vsx,vdx,c(snf(i))

      end do
      
c      close(unit=28) 

c      close(unit=29)

c      close(unit=30)     
      
      !!print*,' '
c      !!print*, "----- FINE CICLO -----"

      print*, "nsn: ", nsn
c      do i=1,nsn
c        !print*,snf(i)
c      end do

c      !print*,"=============="
c      !print*,""

      !!print*,' '
      !!print*,'-------------------------'
      !!print*,' '
      !!print*,' '

c      do j=1,nsn
c        !print*,j,nsn,snf(j)
c      end do

c      !!print*,'uscita'

      return
      end

C-----------------------------------------

      subroutine fault_doppia_insta(sn,nsn_pot,sn_pot,I2_pot,
     *                 U,V,nsn_start,nsn_iter,cf_fix,nelem,
     *                 npoin,nnode,ndime,lnods,coord,nelemp_max,
     *                 nelemp,lelemp,nlpoins,lpoins,xe,ye,time,
     *                 Lx_min,Lx_max,Ly_min,Ly_max,Lo,it,
     *                 xin,yin,nsn,snf,tip,cfl,cf,w)
       
       
      implicit none
       
      Integer*4   npoin,nelem,ipoin,ielem,idime,ndime,inode,
     *            nelemp_max,n,nsn_start,nsn_iter,it,
     *            nnode,lnods(nelem,nnode),nsn_pot,sn_pot(npoin),
     *            nelemp(npoin),lelemp(npoin,nelemp_max),ielemp,
     *            nlpoins(npoin),lpoins(npoin,nelemp_max),j,
     *            i,k,ss,np_pol,nsn,cont,
     *            tip(nelem),divsx,divdx,start,new_nsn
     
      Integer*4 snf(npoin),snf2(nsn_iter+1),w
       
      Double precision Lx_min,Lx_max,Ly_min,Ly_max,vsx,vdx,vtot,
     *                 xin,yin,xfin,yfin,cf_fix,vsxp,vdxp,Lo

      Double precision coord(npoin,ndime),I2_pot(npoin),I2_max,
     *                 v(npoin),u(npoin),xe(nelem),ye(nelem),
     *                 c(npoin),time,cfl(nelem,nnode),sn(npoin),
     *                 I2_0(nelemp_max),cf(npoin),tetaSlip(npoin)

      Double precision,Dimension(:),allocatable:: x_pol, y_pol

      Character(len=100) :: filename
      
      I2_max=MAXVAL(I2_pot,nsn_pot)
      start=nsn+1
      new_nsn=0
      
!      if(it.eq.1) then 
      if(nsn.eq.0) then 
      
c      !print*,' '
c      !print*, '#######################'
      !!print*,"------- NUCLEO --------"
        w=0
        do ipoin=1,npoin       
         !!print*,ipoin,sn_pot(ipoin),'NUCLEO'
         if(I2_pot(ipoin).eq.I2_max.and.
     *     sn_pot(ipoin).eq.1) then   
            snf(start)=ipoin    
            sn_pot(snf(start))=2 
            sn(snf(start))=1 
            cf(snf(start))=cf_fix
            new_nsn=new_nsn+1 
            w=w+1     
            
            write(28,*) snf(start),coord(snf(start),1)*Lo/1.D03,
     *                    coord(snf(start),2)*Lo/1.D03,it,w

            write(30,*) snf(start),coord(snf(start),1)*Lo/1.D03,
     *                    coord(snf(start),2)*Lo/1.D03,it,w
     
            !print*, snf(start),coord(snf(start),1)*Lo/1.D03,
c     *                         coord(snf(start),2)*Lo/1.D03,time,
c     *                         I2_pot(snf(start))
            do ielem=1,nelemp(ipoin)
              if(ipoin.eq.
     *            lnods(lelemp(ipoin,ielem),4)) then
                sn_pot(lnods(lelemp(ipoin,ielem),3))=0
                sn_pot(lnods(lelemp(ipoin,ielem),5))=0
                sn_pot(lnods(lelemp(ipoin,ielem),6))=0
              else if(ipoin.eq.
     *            lnods(lelemp(ipoin,ielem),5)) then
                sn_pot(lnods(lelemp(ipoin,ielem),1))=0
                sn_pot(lnods(lelemp(ipoin,ielem),4))=0
                sn_pot(lnods(lelemp(ipoin,ielem),6))=0
              else if(ipoin.eq.
     *            lnods(lelemp(ipoin,ielem),6)) then
                sn_pot(lnods(lelemp(ipoin,ielem),2))=0
                sn_pot(lnods(lelemp(ipoin,ielem),4))=0
                sn_pot(lnods(lelemp(ipoin,ielem),5))=0
              end if
            end do
          !print*, 'start',snf(start),coord(snf(start),1)*Lo/1.D03,
c     *                         coord(snf(start),2)*Lo/1.D03
            exit              
          endif 
        enddo      

      else 
        start=start-1
      endif

c ---- propagazione faglia in coda
 
      !!print*,"------- PRIMO RAMO FAGLIA --------"
      !print*,"start:",start
      do k=start,npoin
!        if (it.eq.1.and.new_nsn.eq.nsn_start) then
        if (nsn.eq.0.and.new_nsn.eq.nsn_start) then
          exit
!        else if (it.ne.1.and.new_nsn.eq.nsn_iter) then
         else if (nsn.eq.0.and.new_nsn.eq.nsn_iter) then
         exit
        end if
    
      
c        !print*,"K=",k
        I2_0=0.0D0
        ss=0
        n=nlpoins(snf(k))
        !print*,"snf(k)=",snf(k)  
c        !print*,"n=",n
        do i=1,n   
          if(sn_pot(lpoins(snf(k),i)).eq.1) then
            ss=1
            I2_0(i)=I2_pot(lpoins(snf(k),i))    
          endif
        end do

        if(ss.eq.0) then
          !print*,"ss=0"
          exit
        else
          I2_max=MAXVAL(I2_0,nelemp_max)
          do i=1,n
            if(I2_0(i).eq.I2_max.and.
     *           sn_pot(lpoins(snf(k),i)).eq.1) then                    
              snf(k+1)=lpoins(snf(k),i)
              sn(snf(k+1))=1
              cf(snf(k+1))=cf_fix
              sn_pot(snf(k+1))=2 
              new_nsn=new_nsn+1
               w=w+1   
              !print*, k+1,snf(k+1),coord(snf(k+1),1)*Lo/1.D03,
c     *                         coord(snf(k+1),2)*Lo/1.D03,
c     *                         I2_pot(snf(k+1))
              write(28,*) snf(k+1),coord(snf(k+1),1)*Lo/1.D03,
     *                         coord(snf(k+1),2)*Lo/1.D03,it,w

              do ielem=1,nelemp(snf(k+1))

                !print*,ielem,snf(k+1),
c     *           (lnods(lelemp(snf(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf(k+1),ielem),inode)),
c     *           inode=1,6)
                
                if(snf(k+1).eq.
     *           lnods(lelemp(snf(k+1),ielem),4)) then
                  sn_pot(lnods(lelemp(snf(k+1),ielem),3))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),5))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),6))=0
                else if(snf(k+1).eq.
     *             lnods(lelemp(snf(k+1),ielem),5)) then
                  sn_pot(lnods(lelemp(snf(k+1),ielem),1))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),4))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),6))=0
                else if(snf(k+1).eq.
     *             lnods(lelemp(snf(k+1),ielem),6)) then
                  sn_pot(lnods(lelemp(snf(k+1),ielem),2))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),4))=0
                  sn_pot(lnods(lelemp(snf(k+1),ielem),5))=0
                end if
c                !print*,ielem,snf(k+1),
c     *           (lnods(lelemp(snf(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf(k+1),ielem),inode)),
c     *           inode=1,6)
              end do

              
              exit              
            endif   
          enddo    
        endif


      enddo    


      nsn=nsn+new_nsn



c ---- propagazione faglia in testa
      new_nsn=0
      start=snf(1)
      snf2(1)=start
        
      !!print*,"------- SECONDO RAMO FAGLIA --------"
      !print*,"start:",start
       write(28,1000)
1000   format(1H>)

        
      do k=1,nsn_iter
c        if (it.eq.1) then
c          !print*,"it:",it
c          exit
c       end if

c        !print*,"K=",k
        I2_0=0.0D0
        ss=0
        n=nlpoins(snf2(k))
c        !print*,"snf2(k)=",snf2(k)  
c        !print*,"n=",n
        do i=1,n
          if(sn_pot(lpoins(snf2(k),i)).eq.1) then
            ss=1
            I2_0(i)=I2_pot(lpoins(snf2(k),i))    
          endif
        end do

        if(ss.eq.0) then
          !print*,"ss=0"
          exit
        else
          I2_max=MAXVAL(I2_0,nelemp_max)
          do i=1,n
            if(I2_0(i).eq.I2_max .and.
     *           sn_pot(lpoins(snf2(k),i)).eq.1) then                    
              snf2(k+1)=lpoins(snf2(k),i)
              sn(snf2(k+1))=1
              cf(snf2(k+1))=cf_fix
              sn_pot(snf2(k+1))=2 
              new_nsn=new_nsn+1
               w=w+1   
c              !print*, k+1,snf2(k+1),coord(snf2(k+1),1)*Lo/1.D03,
c     *                         coord(snf2(k+1),2)*Lo/1.D03,
c     *                         I2_pot(snf2(k+1))
              write(30,*) snf2(k+1),coord(snf2(k+1),1)*Lo/1.D03,
     *                         coord(snf2(k+1),2)*Lo/1.D03,it,w

              do ielem=1,nelemp(snf2(k+1))

c                !print*,ielem,snf2(k+1),
c     *           (lnods(lelemp(snf2(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf2(k+1),ielem),inode)),
c     *           inode=1,6)
                
                if(snf2(k+1).eq.
     *           lnods(lelemp(snf2(k+1),ielem),4)) then
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),3))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),5))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),6))=0
                else if(snf2(k+1).eq.
     *             lnods(lelemp(snf2(k+1),ielem),5)) then
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),1))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),4))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),6))=0
                else if(snf2(k+1).eq.
     *             lnods(lelemp(snf2(k+1),ielem),6)) then
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),2))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),4))=0
                  sn_pot(lnods(lelemp(snf2(k+1),ielem),5))=0
                end if
c                !print*,ielem,snf2(k+1),
c     *           (lnods(lelemp(snf2(k+1),ielem),inode),
c     *           sn_pot(lnods(lelemp(snf2(k+1),ielem),inode)),
c     *           inode=1,6)
              end do

              
              exit              
            endif   
          enddo    
        endif


      enddo   

c ----- aggiungo il ramo ribaltato in testa a snf (shift) -----

      do i=1,nsn
        j=nsn+1-i
        snf(j+new_nsn)=snf(j)
      end do

      do i=1,new_nsn
        j=new_nsn+2-i
        snf(i)=snf2(j)
      end do
c      !!print*,"chiuso"

      nsn=nsn+new_nsn

c ----- creazione del poligono

      !!print*,' '
      !print*,'------ CREAZIONE POLIGONO ------'
      !print*,"nsn,new_nsn:",nsn,new_nsn
      !!print*,' '
            
      np_pol=nsn+5

        !!print*,'nsn',nsn,'np_pol',np_pol
        allocate(x_pol(np_pol))
        allocate(y_pol(np_pol))
      

c      if(it.eq.1) then
        xin=(3*coord(snf(1),1)-2*coord(snf(2),1))*Lo/1.D03
        yin=(3*coord(snf(1),2)-2*coord(snf(2),2))*Lo/1.D03
c      end if

         !!print*,nsn,snf(nsn),'prima'
      xfin=(3*coord(snf(nsn),1)-2*coord(snf(nsn-1),1))*Lo/1.D03
      yfin=(3*coord(snf(nsn),2)-2*coord(snf(nsn-1),2))*Lo/1.D03
          !!print*,nsn,'dopo'

      x_pol(1)=xin
      y_pol(1)=yin

      do i=1,nsn
        x_pol(i+1)=coord(snf(i),1)*Lo/1.D03
        y_pol(i+1)=coord(snf(i),2)*Lo/1.D03
      end do

      x_pol(nsn+2)=xfin
      y_pol(nsn+2)=yfin
      
      if (xin.ge.xfin.and.
     *    yin.ge.yfin) then
        x_pol(nsn+3)=Lx_min
        y_pol(nsn+3)=Ly_min
        x_pol(nsn+4)=Lx_max
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_max
        y_pol(nsn+5)=Ly_max
        
       else if (xin.ge.xfin.and.
     *          yin.lt.yfin) then
        x_pol(nsn+3)=Lx_min
        y_pol(nsn+3)=Ly_max
        x_pol(nsn+4)=Lx_min
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_max
        y_pol(nsn+5)=Ly_min

      else if (xin.lt.xfin.and.
     *         yin.ge.yfin) then
        x_pol(nsn+3)=Lx_max
        y_pol(nsn+3)=Ly_min
        x_pol(nsn+4)=Lx_min
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_min
        y_pol(nsn+5)=Ly_max

      else if (xin.lt.xfin.and.
     *         yin.lt.yfin) then
        x_pol(nsn+3)=Lx_max
        y_pol(nsn+3)=Ly_max
        x_pol(nsn+4)=Lx_max
        y_pol(nsn+4)=Ly_min
        x_pol(nsn+5)=Lx_min
        y_pol(nsn+5)=Ly_min   
      end if

      do i=1,np_pol
        write(34,*)x_pol(i),y_pol(i),time
        !print*,x_pol(i),y_pol(i),time
      end do


c ------ definizione destro-sinistro



      !!print*,' '
      !!print*,'------- DEFINISCO DX-SX -------'
      tip=0
      do i=1,nsn
        do ielemp=1,nelemp(snf(i))
          ielem=lelemp(snf(i),ielemp)
          call in_out (ielem,xe(ielem),ye(ielem),np_pol,x_pol,y_pol,
     *                 tip(ielem))
          if (tip(ielem).eq.0) tip(ielem)=-1
             write(29,*)  ielem,xe(ielem),ye(ielem),tip(ielem),
     *        snf(i),time
c             !!print*, i, snf(i), ielem, xe(ielem), ye(ielem), tip(ielem)
        end do
      end do

      deallocate(x_pol)
      deallocate(y_pol)

c      close(unit=29)



       !!print*,'CALCOLO TETASLIP'
       
      do i=2,nsn      
        if(coord(snf(i),1).eq.coord(snf(i-1),1)) then
          tetaSlip(snf(i))=90. !codice x slip verticale  
        else if(coord(snf(i),2).eq.coord(snf(i-1),2)) then      
          tetaSlip(snf(i))=0.   !codice x slip orizzontale   
        else 
         tetaSlip(snf(i))=datan((coord(snf(i),2)-coord(snf(i-1),2))/
     *                   (coord(snf(i),1)-coord(snf(i-1),1)))
        endif
     
c        !print*,i,snf(i),tetaSlip(snf(i))*180.D0/3.14D0
      end do






c ------ coupling

      do ipoin=1,npoin
        c(ipoin)=0.0D0
      end do
      !      !print*,' '
c      !print*,"----- Calcolo del coupling -----"
c      write(30,*)"CALCOLO DEL COUPLING"
      !!print*,"i    ipoin     vsx      vdx          c"     
c      write(30,*)"i    ipoin     vsx      vdx          c"

      do i=1,nsn
c        !print*,i
        vsx=0
        vdx=0
        divsx=0
        divdx=0
        do ielemp=1,nelemp(snf(i))
          ielem=lelemp(snf(i),ielemp)
          do inode=1,nnode
            ipoin=lnods(ielem,inode)
            vtot=sqrt(U(ipoin)**2+V(ipoin)**2)
            if(tip(ielem).eq.-1)then
              vsx=vsx+vtot
              divsx=divsx+1
            else if(tip(ielem).eq.1.and.sn(ipoin).eq.1)then
              vdx=vdx+vtot
              divdx=divdx+1
            end if
          end do
        end do
        
        vsx=vsx/divsx
        vdx=vdx/divdx
         
        vsxp=vsx*dcos(tetaslip(snf(i)))+vsx*dsin(tetaslip(snf(i)))
        
        vdxp=vdx*dcos(tetaslip(snf(i)))+vdx*dsin(tetaslip(snf(i)))
        
        if(vsxp*vdxp.lt.0.) c(snf(i))=-1.
        if(vsxp*vdxp.gt.0.) c(snf(i))=+1.
        
        c(snf(i))=+1.       !da togliere quando aggiusteremo
        
         cf(snf(i))=cf(snf(i))*c(snf(i))


        do ielemp=1,nelemp(snf(i))
          ielem=lelemp(snf(i),ielemp)
          do inode=1,nnode
            ipoin=lnods(ielem,inode)
            if(tip(ielem).eq.-1)then
              cfl(ielem,inode)=1.0D0
            else if(tip(ielem).eq.1.and.sn(ipoin).eq.1)then
              cfl(ielem,inode)=cf(snf(i))
            end if
          end do
        end do








      !        !print*,divsx,divdx

c        if(divsx.ne.0.and.divdx.ne.0) then
c          vsx=vsx/divsx
c          vdx=vdx/divdx
c           c(snf(i))=vdx/vsx
c        else 
c          c(snf(i))=0.5
c        endif 




      !        !print*,i,snf(i),vsx,vdx,c(snf(i))
c        write(30,*) i,snf(i),vsx,vdx,c(snf(i))

      end do
      
c      close(unit=28) 

c      close(unit=29)

c      close(unit=30)     
      
      !!print*,' '
c      !!print*, "----- FINE CICLO -----"

      print*, "nsn: ", nsn
c      do i=1,nsn
c        !print*,snf(i)
c      end do

c      !print*,"=============="
c      !print*,""

      !!print*,' '
      !!print*,'-------------------------'
      !!print*,' '
      !!print*,' '

c      do j=1,nsn
c        !print*,j,nsn,snf(j)
c      end do

c      !!print*,'uscita'

      return
      end

c----------------------------------------------------------------------------
       SUBROUTINE IN_OUT (ielem,x,y,np_pol,x_pol,y_pol,iadentro)

c      Determina se un punto si trova dentro un poligono
c      -------------------------------------------------
       
       implicit integer (i-n)
       implicit double precision (a-h,o-z)

       integer np_pol
       dimension x_pol(np_pol)
       dimension y_pol(np_pol)
       
	circulacion=0
        pol_long=0
        parcial_circ=0
	irecuerdo=0	   

        do 10 i=1,np_pol
            x0=x_pol(i)
            y0=y_pol(i)
	    
	    if (irecuerdo.eq.1) then
		x0=xf
		y0=yf
	    endif	

            if (i.lt.np_pol) then
                 xf=x_pol(i+1)
                 yf=y_pol(i+1)
            else
                 xf=x_pol(1)
                 yf=y_pol(1)
            endif

	   irecuerdo=0
	   if (xf.eq.x0.and.yf.eq.y0) then
		irecuerdo=1
		if (i.lt.np_pol-1) then
			xf=(x_pol(i+2)+x0)/2
			yf=(y_pol(i+2)+y0)/2
		else
			if (i.eq.np_pol) then
				xf=(x_pol(2)+x0)/2
				yf=(y_pol(2)+x0)/2
			endif
			if (i.eq.np_pol-1) then
				xf=(x_pol(1)+x0)/2
				yf=(y_pol(1)+x0)/2
			endif
		endif
	   endif       
            
            D2=(x0-x)**2+(y0-y)**2
            E2=(xf-x)**2+(yf-y)**2
	    rL2=(xf-x0)**2+(yf-y0)**2

c           Distanza fra il punto (x,y) e il primo punto (x0,y0) del segmento
	    D=sqrt(D2)

c           Distanza fra il punto (x,y) e il secondo punto (xf,yf) del segmento
	    E=sqrt(E2)
	   
c           Distanza fra i due punti estremi (x0,y0) e (xf,yf) del segmento 	    
            rL=sqrt(rL2)
	    
c           lunghezza parziale del poligono	    
	    pol_long=pol_long+rL
            
c           Prodotto vettoriale fra i vettori (x,y)-(x0,y0) e (x,y)-(xf,yf),
c           diviso la lunghezza rL

	    rK=((x0-x)*(yf-y0)-(y0-y)*(xf-x0))/rL
	    
c           Prodotto scalare fra i vettori (x,y)-(x0,y0) e (x,y)-(xf,yf),
c           diviso la lunghezza rL, moltiplicato per 2
   
             b=2*((y0-y)*(yf-y0)+(x0-x)*(xf-x0))/rL

c           Punto (x,y) coincidente con il primo punto del segmento (x0,y0)

            if(D.eq.0.0) then 
	      iadentro=1	    
	      go to 100
	    
c          Punto (x,y) coincidente con il secondo punto del segmento (xf,yf)

           elseif(E.eq.0.0) then 
	     iadentro=1
	     go to 100

c          punto sul segmento (x0,y0)-(xf,yf)
           	  
           elseif(rK.eq.0.0) then
	     DE_rL=(D+E)-rL
	     if(DE_rL.le.1.D-05) DE_rL=0.0
	     if(DE_rL.eq.0.0) then 
	       iadentro=1
	       go to 100	     	  
	     endif

           else 
	   
            disc=b*b-4*D2
           
	    if (disc.ge.0) then
                 sq=sqrt(disc)
            else
                 sq=sqrt(-disc)
            endif

c           OJO, RESTITUIR LA LINEA SIGUIENTE:
            if (sq.eq.0) sq=0.001

            if (disc.ge.0) then
                 parcial_circ=rK/sq*(log((2*rL+b-sq)/(2*rL+b+sq))-
     &           log((b-sq)/(b+sq)))
            else
                 parcial_circ=rK/sq*(atan((2*rL+b)/sq)-atan(b/sq))
            endif
            circulacion=circulacion+parcial_circ
           
       endif
       	       	    
10     continue

       if (abs(circulacion).ge.3.14159-.1) then
       
        iadentro=1
       else
         iadentro=0
       endif

100    return

	end

 
 
 
 
c--------------------------------------------------------------------    
      
      subroutine creazione_lelemp_lpoins_NV3(nelem,npoin,nnode,ndime,
     *           nelemp_max,coord,lnods,nelemp,lelemp,
     *           nlpoins,lpoins,dist)
      
      implicit none

      integer*4 npoin,ipoin,jpoin,i,j,ielem,nelem,idime,ndime,nnode,
     *            inode,k,nelemp_max
      
      Integer*4 lnods(nelem,nnode),numel,nelemp(npoin)
        
      Double precision coord(npoin,ndime),dist,dij

      integer*4 lelemp(npoin,nelemp_max),lpoins(npoin,npoin),
     *              nlpoins(npoin)
        
        
      character*12 testo

      open(unit=21,file='Corrispondenze_Nodi2Elementi',status='unknown')


      open(unit=234,file='Connessioni_Dirette_Nodi',status='unknown')
      open(unit=235,file='Connessioni_Dirette_Nodi.xy',status='unknown')
  
      lelemp=0
      lpoins=0
 
        
      do ipoin=1,npoin

        k=0
        do ielem=1,nelem
            
          do inode=1,nnode
            if(ipoin.eq.lnods(ielem,inode)) then
              k=k+1
              lelemp(ipoin,k)=ielem
            endif 
          
          enddo
        
          nelemp(ipoin)=k        
        enddo
        
        
1     write(21,*) ipoin,nelemp(ipoin),
     *       (lelemp(ipoin,i),i=1,nelemp(ipoin))
      
      enddo 



        do ipoin=1,npoin
c           print *,ipoin           
          k=0
          
          do jpoin=1,npoin
          
c            print *,ipoin,jpoin
            
            if(jpoin.ne.ipoin) then 
            
            dij=dsqrt((coord(ipoin,1)-coord(jpoin,1))**2+
     *                (coord(ipoin,2)-coord(jpoin,2))**2)
          
            if(dij.le.dist) then
            k=k+1
            lpoins(ipoin,k)=jpoin          
            endif
            
            endif
          enddo
          
          nlpoins(ipoin)=k
        enddo
        
        
        do ipoin=1,npoin
          write(234,*) ipoin,nlpoins(ipoin),
     *       (lpoins(ipoin,k),k=1,nlpoins(ipoin))
           
           do k=1,nlpoins(ipoin)
             write(235,*) ipoin,lpoins(ipoin,k),
     *                    coord(lpoins(ipoin,k),1)*700.,
     *                    coord(lpoins(ipoin,k),2)*700.
           enddo
        enddo
      
      close(unit=21)
      close(unit=234)
      close(unit=235)
      
      return
      end
      
      

 
 
 
 
 
 