      program Test_Regrid

      Implicit none

      Logical linp

      Character*12 testo
 
      Integer*4 marcatori,nmats,nprop,ndofn,ndofn_T,nelem,npoin,
     *          nsvab,nsvab_T,new_nelem,new_npoin,
     *          nnode,ndime,ngaus,nstre,nevab,nevab_T,
     *          nvfix,nvfix_T,nTout,harmonic,idra,p,phaseC,PhaseCc,
     *          Grav,celleT,celleP,nelemp_max,nsn_start,nsn_iter

      Double Precision tfine,tsteps,Sl,vel_in,Lx_min,
     *           Lx_max,Ly_min,Ly_max,dist,
     *           cont_y,ocean_y,xsub,vh2o,dist_marc,ht_max,
     *           PressMin,PressMax,TempMin,TempMax,cf_fix
                                                        
       Parameter(marcatori=455489,                ! Numero totale di marcatori 
!     *          nelem=6813,                    ! Numero totale di Elementi
!     *    npoin=13846,                      ! Numero totale di punti nodali
     *    nmats=1,                          ! Numero totale di materiali
     *    nprop=3,                          ! Numero di proprietà per ogni materiale
     *    ndofn=2,                          ! Numero di gradi di libertà  MOMENTO
     *    ndofn_T=1,                        ! Numero di gradi di libertà TEMPERATURA
     *    nnode=6,                          ! Numero di nodi per elemento
     *    ndime=2,                          ! Dimensioni del problema
     *    ngaus=7,                          ! Numero di punti di Gauss 
     *    nstre=3,                          ! Numero di componenti dello stress
     *      nelemp_max=15,                     ! Massimo numero di elementi a cui un nodo può appartenere
     *      nsn_start=3,                      ! Numero di slip nodes creati alla prima iterazione
     *      nsn_iter=1,                       ! Numero di slip nodes creati in ogni ramo ad ogni iterazione
     *          nevab=ndofn*nnode,
!     *    nsvab=ndofn*npoin,
     *    nevab_T=ndofn_T*nnode,
!     *          nsvab_T=ndofn_T*npoin,
     *          nvfix=1,                        ! Numero di punti a velocità fissata - 584 Attiva - 557 Grav
     *    nvfix_T=1,                        ! Numero di punti a Temperatura fissata
     *          tfine=55.D0,                      ! Tempo massimo del modello (Myr)
     *          tsteps=0.05D0,                    ! Time step (Myr)
     *          nTout=10,
     *          Sl=80.0D03,                        ! Spessore termico iniziale della litosfera
     *          vel_in=5.D0,                      ! Velocita di ingresso dei marcatori crostali: ATTUALMENTE SIMMETRICA!!!!
     *          Lx_min=-276.0D0,                      ! Limite sinistro modello
     *          Lx_max=700.0D0,                       ! Limite destro modello
     *          Ly_min=-700.0D0,                      ! Limite inferiore modello
     *          Ly_max=0.0D0,                        ! Limite superiore modello
     *          cont_y=-30.0D0,                    ! Y della crosta continentale in km!
     *          ocean_y=-10.0D0,                    ! Y della crosta oceanica in km!
     *          xsub=-100.0D0,                      ! X del trench
     *          vh2o=20.0D0,                      ! velocità dei fluidi in cm/yr
     *          dist_marc=0.5D0,                    ! Distanza tra i marcatori
     *          harmonic=0,                      ! Average type
     *          idra=1,                          ! Hydration yes/no
     *          p=10,                          ! Pressure quadrature points
     *          phaseC=1,                        ! Cambiamenti fase Mantello
     *          phaseCc=1,                        ! Cambiamenti fase Crosta
     *          Grav=0,                          ! Fase Gravitazionale/Attiva
     *          ht_max=2.5D0,                      ! Altezza massima della topografia
     *          celleT=100,                      ! Punti per la Temperatura di Perplex
     *          celleP=100,                      ! Punti per la Pressione di Perplex
     *          PressMin=0.0D0,                    ! Pressione minima del file Perplex
     *          PressMax=7.D0,                    ! Pressione massima del file Perplex
     *          TempMin=300.0D0,                    ! Temperatura minima del file Perplex
     *          TempMax=1600.0D0,                    ! Temperatura massima del file Perplex
     *          dist=5,                             ! Distanza entro cui nodi sono vicini
     *          cf_fix=1)                         ! Coupling factor

      Integer*4 ipoin,ielem,idime,inode,ivfix,idofn,
     *          ivfix_T,
     *          numel,i

      Double Precision Uo

      Integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T)
    
      Double Precision presc(nvfix,ndofn), presc_T(nvfix_T,ndofn_T)

      

C --------------------- DA ALLOCARE 

c     nelem
      Double Precision, Dimension(:), allocatable ::ymin,ymax,Po,Xe,Ye,
     *                  xmin,xmax,Area                       
      
     
      
      !npoin
      Integer*4, Dimension(:), allocatable :: nelemp,nlpoins

      Integer*4, Dimension(:,:), allocatable :: lnods, !nelem,nnode
     *          lelemp, !npoin,nelemp_max
     *          lpoins !npoin,npoin

      Double Precision, Dimension(:,:), allocatable ::
     *                 coord!(npoin,ndime)
      

c______________________________________________________________________________________ 
c      Variabili per Slip Nodes

      
      Integer*4 nsn,sn2
      
      Integer*4, Dimension(:), allocatable :: tip !nelem
      
      Integer*4, Dimension(:), allocatable :: snf !npoin
      
      Double Precision, Dimension(:,:), allocatable :: cfl !nelem,nnode

      Double Precision, Dimension(:), allocatable :: cf,sn !npoin

      nsn=0
c      
c______________________________________________________________________________________
c       Allocazione variabili

      nelem=32
      npoin=81
      nsvab=ndofn*npoin
      nsvab_T=ndofn_T*npoin

      
      
      allocate(lnods(nelem,nnode))
    
      allocate(coord(npoin,ndime))

      allocate(tip(nelem))
      allocate(snf(npoin))
      allocate(cfl(nelem,nnode))
      allocate(cf(npoin))
      allocate(sn(npoin)) 
      
      
c---------------------------------------------------------------------------
c      COPIA FILE DI INPUT


      inquire(file='Input_2D_iso',exist=linp)
      !inquire(file='Input',exist=linp)
      !inquire(file='Input_30-51',exist=linp)
      if(.not.linp) print *,'ERRORE - il file Input non esiste!'
      
      open(unit=15, file='Input_2D_iso', status='old')
      !open(unit=15, file='Input', status='old')
      !open(unit=15, file='Input_30-51', status='old')
      read(15,1009) testo
      write(14,1009) testo
      
      do i=1,npoin
        read(15,1003) ipoin,(coord(ipoin,idime),idime=1,ndime),              ! Lettura coordinate Nodali
     *             sn(ipoin),cf(ipoin)
        write(14,1003) ipoin,(coord(ipoin,idime),idime=1,ndime),              ! Lettura coordinate Nodali
     *             sn(ipoin),cf(ipoin)
      enddo     
      
      read(15,1009) testo
      write(14,1009) testo

      do ielem=1,nelem
        read(15,1005) numel,tip(ielem),(lnods(numel,inode),                  ! Lettura: Numero Elemento, Tipo di Elemento, Connessioni nodali per ogni elemento
     *            inode=1,nnode)                        ! Calcolo proprietà geometriche dell'elemento
        write(14,1005) numel,tip(ielem),(lnods(numel,inode),                  ! Lettura: Numero Elemento, Tipo di Elemento, Connessioni nodali per ogni elemento
     *            inode=1,nnode)
      enddo      
    
           
      read(15,1009) testo                         
      write(14,1009) testo

      do ivfix=1,nvfix 
        read(15,*) nofix(ivfix),(ifpre(ivfix,idofn), idofn=1,ndofn),
     *               (presc(ivfix,idofn), idofn=1,ndofn)
        write(14,*)nofix(ivfix),(ifpre(ivfix,idofn), idofn=1,ndofn),
     *               (presc(ivfix,idofn), idofn=1,ndofn)
      enddo  

         
      read(15,1009) testo
      write(14,1009) testo
      do ivfix_T=1,nvfix_T
        read(15,*) nofix_T(ivfix_T),ifpre_T(ivfix_T,1),
     *      presc_T(ivfix_T,1)                           ! Temperatura normalizzata
        write(14,*) nofix_T(ivfix_T),ifpre_T(ivfix_T,1),
     *      presc_T(ivfix_T,1)    
    
      enddo        

      close(15)
      close(14)

c---------------------------------------------------------------------------
c         LETTURA FILE DI INPUT

      open(unit=14, file='fort.14', status='old')
      read(14,1009) testo

      do i=1,npoin
        read(14,1003) ipoin,(coord(ipoin,idime),idime=1,ndime),              ! Lettura coordinate Nodali
     *             sn(ipoin),cf(ipoin)
        coord(ipoin,1)=coord(ipoin,1)/700.                       ! Coordinate normalizzate
        coord(ipoin,2)=coord(ipoin,2)/700.    
      enddo     
      
      read(14,1009) testo

      do ielem=1,nelem        
        read(14,1005) numel,tip(numel),(lnods(numel,inode),                  ! Lettura: Numero Elemento, Tipo di Elemento, Connessioni nodali per ogni elemento
     *            inode=1,nnode)                      ! Calcolo proprietà geometriche dell'elemento
        do inode=1,nnode
          if(tip(ielem).eq.0.or.tip(ielem).eq.-1) then
            cfl(ielem,inode)=1.0D0
          elseif(tip(ielem).eq.+1) then
            cfl(ielem,inode)=cf(lnods(ielem,inode))
          endif
        enddo      
      enddo

      
c -------------------------------------------------------------------------------
c         LETTURA VELOCITA FISSATE       
            
      read(14,1009) testo                         

      do ivfix=1,nvfix 
        read(14,*) nofix(ivfix),(ifpre(ivfix,idofn), idofn=1,ndofn),
     *               (presc(ivfix,idofn), idofn=1,ndofn)
        presc(ivfix,1)=presc(ivfix,1)/Uo                          ! Velocita normalizzate
        presc(ivfix,2)=presc(ivfix,2)/Uo
      enddo  

c -------------------------------------------------------------------------------
c         LETTURA TEMPERATURE FISSATE 

      if(nvfix_T.lt.0.or.nvfix_T.gt.npoin) write(6,*) 'nvfix_T errato' 
      if(nvfix_T.lt.0.or.nvfix_T.gt.npoin) stop            
      read(14,1009) testo
      do ivfix_T=1,nvfix_T
        read(14,*) nofix_T(ivfix_T),ifpre_T(ivfix_T,1),
     *      presc_T(ivfix_T,1)                           ! Temperatura normalizzata
      enddo        

      close(unit=14)


C===============================================
C================== REGRID =====================

      new_nelem=nelem
      new_npoin=npoin

      ! definire la faglia qui
      snf=0
      nsn=2 ! Numero di slip nodes da creare -> devo mettere un snf(i) per ognuno di questi


      snf(1)=32 ! definisco gli slip nodes in ordine, tramite la numerazione della griglia
      snf(2)=33
      !snf(3)=23
      !snf(4)=32
      !snf(5)=50
      !snf(6)=60
      !snf(7)=70
      !snf(8)=81

       

      do i=1,nsn
        write(7,*) snf(i),coord(snf(i),1)*700.0D0,
     *   coord(snf(i),2)*700.0D0
      end do
      close(7)

      sn2=snf(1)
      do i=1,npoin
        print*,"IT:",i,snf(i),snf(i+1),sn2

        if(snf(i).eq.sn2) then
          sn2=snf(i+1)
        end if  
          if(snf(i+1).eq.0) exit
          call RegridFaglia(snf(i),snf(i+1),snf,nsn,
     *     new_npoin,new_nelem,ndime,nnode,nelemp_max,
     *           cf_fix,nvfix,nofix,ifpre,ndofn,presc,Uo,
     *    nvfix_T,nofix_T,ifpre_T,presc_T)
        
      end do
c     Faccio girare di nuovo per riempire eventuali buchi

      print*,""
      print*, "       --- SECONDO GIRO - controllo buchi ---"
      print*,""

      sn2=snf(1)
      do i=1,npoin
        print*,"IT:",i,snf(i),snf(i+1),sn2

        if(snf(i).eq.sn2) then
          sn2=snf(i+1)
        end if
          if(snf(i+1).eq.0) exit
          call RegridFaglia(snf(i),snf(i+1),snf,nsn,
     *     new_npoin,new_nelem,ndime,nnode,nelemp_max,
     *           cf_fix,nvfix,nofix,ifpre,ndofn,presc,Uo,
     *    nvfix_T,nofix_T,ifpre_T,presc_T)
       
      end do

      print*,"----------------"
      print*,"Fine regrid:"
      print*,"nelem",nelem
      print*,"npoin",npoin
      print*,"new_nelem",new_nelem
      print*,"new_npoin",new_npoin
      print*,"nsn",nsn
      print*,"nnode",nnode


C================================================
C================================================

      !print*,""
      !print*,"Deallocazione"

      deallocate(lnods)
      !print*,"a"

      
      !deallocate(nlpoins)
      !deallocate(lpoins)
      !deallocate(lelemp)
      !deallocate(nelemp)
      !print*,"b"

      deallocate(coord)

      deallocate(tip)
      deallocate(cfl)
      deallocate(cf)
      deallocate(sn)

      nelem=new_nelem
      npoin=new_npoin

      allocate(Xe(nelem))
      allocate(Ye(nelem))
      allocate(xmin(nelem))
      allocate(xmax(nelem))
      allocate(Area(nelem))
      allocate(ymin(nelem))
      allocate(ymax(nelem))
      allocate(Po(nelem))
      
      allocate(lnods(nelem,nnode))
    
      allocate(coord(npoin,ndime))

      allocate(tip(nelem))
      allocate(cfl(nelem,nnode))
      allocate(cf(npoin))
      allocate(sn(npoin)) 

      allocate(nelemp(npoin))
      allocate(nlpoins(npoin))
      allocate(lelemp(npoin,nelemp_max))
      allocate(lpoins(npoin,npoin))
      
      !print*,"Ricalcolo variabili"
      open(unit=14, file='fort.14', status='old')

      read(14,1009) testo
      
      do i=1,npoin
        
        read(14,*) ipoin,(coord(ipoin,idime),idime=1,ndime),          ! Lettura coordinate Nodali
     *             sn(ipoin),cf(ipoin)

        coord(ipoin,1)=coord(ipoin,1)/700.                       ! Coordinate normalizzate
        coord(ipoin,2)=coord(ipoin,2)/700.    
      enddo     
      
      read(14,1009) testo

      Xe=0.0D0
      Ye=0.0D0
      do ielem=1,nelem      
        read(14,*) numel,tip(numel),(lnods(numel,inode),                  ! Lettura: Numero Elemento, Tipo di Elemento, Connessioni nodali per ogni elemento
     *            inode=1,nnode)                      ! Calcolo proprietà geometriche dell'elemento
        do inode=1,nnode
          if(tip(ielem).eq.0.or.tip(ielem).eq.-1) then
            cfl(ielem,inode)=1.0D0
          elseif(tip(ielem).eq.+1) then
            cfl(ielem,inode)=cf(lnods(ielem,inode))
          endif
        enddo
      enddo      

c -------------------------------------------------------------------------------
c         LETTURA VELOCITA FISSATE       
            
      read(14,1009) testo                         

      do ivfix=1,nvfix 
        read(14,*) nofix(ivfix),(ifpre(ivfix,idofn), idofn=1,ndofn),
     *               (presc(ivfix,idofn), idofn=1,ndofn)
        presc(ivfix,1)=presc(ivfix,1)/Uo                          ! Velocita normalizzate
        presc(ivfix,2)=presc(ivfix,2)/Uo
      enddo  

c -------------------------------------------------------------------------------
c         LETTURA TEMPERATURE FISSATE 

      if(nvfix_T.lt.0.or.nvfix_T.gt.npoin) write(6,*) 'nvfix_T errato' 
      if(nvfix_T.lt.0.or.nvfix_T.gt.npoin) stop            
      read(14,1009) testo
      do ivfix_T=1,nvfix_T
        read(14,*) nofix_T(ivfix_T),ifpre_T(ivfix_T,1),
     *      presc_T(ivfix_T,1)                           ! Temperatura normalizzata
      enddo        

      close(14)
            
  
c -------------------------------------------------------------------------------
c        CREAZIONE LELEMP LPOINS

      call neighboring_gridbounded(nelem,npoin,nnode,ndime,
     *           nelemp_max,coord,lnods,nelemp,lelemp,
     *           nlpoins,lpoins) 

      print*,"Nuovo snf"
      do i=1,nsn
        print*,snf(i)
        write(8,*) snf(i),coord(snf(i),1)*700.0D0,
     *   coord(snf(i),2)*700.0D0
      end do

      print*,"----------------"

      close(8)
      print*,"END PROGRAM"
1009  format(a12)
1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)     
1005  format(10(1x,i10),5x,6(i2)) 

      end program

C============================================================
C============================================================

      subroutine RegridFaglia(sn1,sn2,snf,nsn,
     *           npoin,nelem,ndime,nnode,nelemp_max,
     *           cf_fix,nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
     

C     RICALCOLO DELLA MESH QUADRATICA ATTORNO ALLA FAGLIA

      implicit none
      
      character*100 linea
      Integer*4 i,j,k,nnode,npoin,nelem,ielem,ndime,
     *          nsn,
     *          sn1,sn2,nelemp_max
           
      Integer*4 control

      Integer*4  ee(nelem),tt(nelem), !nelem
     *              nn(npoin),snf(npoin),nelemp(npoin) !npoin

      integer*4 lnods(nelem,nnode),lelemp(npoin,nelemp_max)


      double precision xx(npoin),yy(npoin),sn(npoin),  !npoin
     *                 cf(npoin),cf_fix

      double precision coord(npoin,ndime),
     *                 asci(nnode,nelem),
     *                 ordi(nnode,nelem)

      integer*4 nvfix,nvfix_t,ndofn,ndofn_T

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T),
     *          lpoins(npoin,npoin),
     *          nlpoins(npoin)


      double precision presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

c     Apertura file di input

      print*,"--INIZIO REGRID--"
      Open(unit=1,file='fort.14',status='old')
      read(1,*) linea 

c     Legge le coordinate nodali e le scrive nel file "Nodi"
      
      Open(unit=2,file='Nodi.tmp',status='unknown')
      close(unit=2,status='delete')
      Open(unit=2,file='Nodi.tmp',status='new')   

      !print*,"npoin",npoin
      !print*,"nelem",nelem

      do i=1,npoin
        read(1,*) nn(i),xx(i),yy(i),sn(i),cf(i)
        write(2,1003) nn(i),xx(i),yy(i),sn(i),cf(i)     
        coord(i,1)=xx(i)
        coord(i,2)=yy(i)  
      enddo
       
      read(1,*)  linea 
             
c     Legge le connessioni nodali c(k,j) e il tipo di materiale t(j)
       
      Open(unit=3,file='Elementi.tmp',status='unknown')
      close(unit=3,status='delete')
      Open(unit=3,file='Elementi.tmp',status='new')
      do j=1,nelem
        read(1,*) ee(j),tt(j),(lnods(j,k),k=1,nnode) !elementi, tipo materiale, connessioni nodali
      enddo
      
      close(unit=1) 

      
c.....Per ogni nodo di connessione sostituisce il numero del nodo 
c     con le coordinate nodali

      do j=1,nelem
        do k=1,nnode
          do i=1,npoin
            if(lnods(j,k).eq.nn(i)) then
              asci(k,j)=xx(i)
              ordi(k,j)=yy(i)
            endif
          end do
        end do        
      end do 

c---- connessioni nodi-elementi

c      call RicalcoloLelemp(nelem,npoin,nnode,ndime,
c     *           nelemp_max,lnods,nelemp,lelemp)
      control=0

      call neighboring_gridbounded(nelem,npoin,nnode,ndime,
     *           nelemp_max,coord,lnods,nelemp,lelemp,
     *           nlpoins,lpoins)


c     verifico a che elementi appartengono i due slip nodes
      do i = 1,nelemp(sn1)
        do j = 1,nelemp(sn2)
c       se appartengono allo stesso elemento 3 casi - subroutine 1 elem
          if(lelemp(sn1,i).eq.lelemp(sn2,j)) then
            ielem=lelemp(sn1,i)
c         taglia un solo elemento, quindi 4 casi
c         - corner_side
c         - side_side
c         - corner_corner
c         - no_cut
            print *,"Stesso elemento:", ielem,lelemp(sn2,j)

            call taglio1elem(sn1,sn2,snf,ielem,nelem,npoin,nsn,
     *          nnode,asci,ordi,ee,tt,nelemp_max,nelemp,ndime,
     *           lelemp,lnods,cf_fix,nvfix,nofix,ifpre,ndofn,
     *          presc,Uo,
     *          nvfix_T,nofix_T,ifpre_T,presc_T)
              close(2)
              close(3)

            control=1

            go to 200
          end if
        end do
      end do
200   continue

      if(control.ne.1) then
        print *,"Elementi diversi"
        call taglio2elem(coord,asci,ordi,
     *   nelem,npoin,ndime,sn1,sn2,nelemp_max,nelemp,lelemp,
     * nnode,lnods,snf,nsn,ee,tt,cf_fix,
     * nvfix,nofix,
     *     ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
        close(2)
        close(3)
      end if


c         - intersezione sn1 coincide con intersezione sn2 - subroutine 2 elem
c         - intersezione sn1 non coincide con intersezione sn2
c           --> la faglia taglia elementi tra ielem e jelem = kelem1, kelem2, kelem3...
c           trovo il kelem che confina con ielem
c           trovo il kelem2 che confina con kelem1
c           finché kelemn = jelem
c           per ogni coppia richiamo subroutine 2 elem


1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)     

      end

c--------------------------------------
      subroutine taglio1elem(sn1,sn2,snf,ielem,nelem,npoin,nsn,
     *           nnode,asci,ordi,ee,tt,nelemp_max,nelemp,ndime,
     *           lelemp,lnods,cf_fix,nvfix,nofix,ifpre,ndofn,
     *          presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
      implicit none     
      integer*4 nelem,ielem,nnode,inode,node1,node2,node,
     *        npoin,sn1,sn2,k,i,nsn,nelemp_max,new_nelem

      Integer*4 ee(nelem),tt(nelem),
     *          lnods(nelem,nnode),snf(npoin),
     *          lelemp(npoin,nelemp_max),nelemp(npoin)

      double precision x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,
     *                 cf_fix 

      double precision asci(nnode,nelem),ordi(nnode,nelem)

      integer*4 ndime,nvfix,nvfix_t,ndofn,ndofn_T

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T)

      double precision presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

      new_nelem=nelem
      do inode=1,nnode
        if(lnods(ielem,inode).eq.sn1) node1=inode
        if(lnods(ielem,inode).eq.sn2) node2=inode
      end do
      node=MIN(node1,node2)
      node2=MAX(node1,node2)
      node1=node

      print*,"node1,2:",node1,node2

      x1=asci(1,ielem)
      y1=ordi(1,ielem) 
    
      x2=asci(2,ielem)
      y2=ordi(2,ielem)  
    
      x3=asci(3,ielem)
      y3=ordi(3,ielem)
      
      x4=asci(4,ielem)
      y4=ordi(4,ielem)  
    
      x5=asci(5,ielem)
      y5=ordi(5,ielem)  
    
      x6=asci(6,ielem)
      y6=ordi(6,ielem)

      if(node1.eq.1.and.node2.eq.5.or. 
     *     node1.eq.3.and.node2.eq.4.or.
     *     node1.eq.2.and.node2.eq.6) then
        print *,"  Caso corner-side",ielem
        call corner_side(ielem,node1,node2,sn2,snf,nsn,nelem,
     *    new_nelem,nnode,npoin,nelemp_max,lelemp,nelemp,asci,ordi,
     *      ndime, ee,tt,cf_fix,lnods,
     *    nvfix,nofix,ifpre,ndofn,
     *          presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      else if(node1.eq.1.and.node2.eq.3.or. 
     *          node1.eq.1.and.node2.eq.2.or.
     *          node1.eq.2.and.node2.eq.3) then
        print *,"  Caso corner-corner",ielem
        call corner_corner(ielem,node1,node2,sn2,snf,
     *   nsn,nelem,new_nelem,nnode,npoin,ndime,ee,tt,
     *   lnods,nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      else if(node1.eq.4.and.node2.eq.6.or. 
     *          node1.eq.4.and.node2.eq.5.or.
     *          node1.eq.5.and.node2.eq.6) then
      print *,"  Caso side-side",ielem
    
        call side_side(node1,node2,ielem,nelem,new_nelem,npoin,
     *   lnods,nnode,nelemp,lelemp,nelemp_max,snf,ndime,asci,ordi,
     * nsn,sn2,ee,tt,cf_fix,nvfix,nofix,ifpre,ndofn,
     *          presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
      else
        print *,"  Caso nodi adiacenti sullo stesso bordo",ielem
        do i=1,nelem
           write(3,1005) ee(i),tt(i),(lnods(i,k),k=1,nnode)
        end do
      end if

      close(2)
      close(3)

      nelem=new_nelem
      print*,""
 
1005  format(10(1x,i10),5x,6(i2))    

      end

c--------------------------------------
      subroutine taglio2elem(coord,asci,ordi,
     *   nelem,npoin,ndime,sn1,sn2,nelemp_max,nelemp,lelemp,
     * nnode,lnods,snf,nsn,ee,tt,cf_fix,nvfix,nofix,
     *     ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
      
      implicit none

      integer*4 nelem,npoin,ndime,sn1,sn2,i,j,nelemp_max,
     *      nnode,sn_new,control,inode,nsn,ielem,node1,ielemp,
     *      node2,node3,new_nelem,a,b

      integer*4 nelemp(npoin),lelemp(npoin,nelemp_max),
     * lnods(nelem,nnode),snf(npoin),ee(nelem),
     * tt(nelem),nordi(nnode)

      double precision coord(npoin,ndime),asci(nnode,nelem),
     * ordi(nnode,nelem)

      double precision xmin,xmax,ymin,ymax,m_r,q_r,m_e,q_e,
     *                xe_min,xe_max,ye_min,ye_max,xi,yi,
     *               cf_fix,dl

      integer*4 nvfix,nvfix_t,ndofn,ndofn_T

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T)

      double precision presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

      print*,"Taglio 2 elem"
      
      new_nelem=nelem
      dl=1.25D0
      dl=2.5D0
       
c         creo la retta tra i due sn

      if(coord(sn1,1).ne.coord(sn2,1)) then
        m_r=(coord(sn1,2)-coord(sn2,2))/
     *      (coord(sn1,1)-coord(sn2,1))
        q_r=coord(sn1,2) - m_r*coord(sn1,1)
      else
        m_r=0
        q_r=0
      end if

c         creo la finestra entro cui può esserci intersezione
      xmin=MIN(coord(sn1,1),coord(sn2,1))
      xmax=MAX(coord(sn1,1),coord(sn2,1))
      ymin=MIN(coord(sn1,2),coord(sn2,2))
      ymax=MAX(coord(sn1,2),coord(sn2,2))
      !print*,xmin,xmax,ymin,ymax        
c         vado  a vedere tutti gli elementi a cui il primo  sn appartiene
      do i=1,nelemp(sn1)
c           per ogni elemento valuto la retta passante per ogni lato 
c           (4=retta da nodo 1 a 2, 5= retta tra 2 e 3,  6=retta tra 3 e 1)
c           e l'intersezione con la retta passante per gli sn
        ielemp=lelemp(sn1,i)
        do j=4,6
          if(j.eq.4) then
            a=1
            b=2
          else if (j.eq.5) then
            a=2
            b=3
          else if(j.eq.6) then
            a=3
            b=1
          end if

          if(asci(a,ielemp).ne.asci(b,ielemp)) then
            m_e=(ordi(a,ielemp)-ordi(b,ielemp))/
     *          (asci(a,ielemp)-asci(b,ielemp))
            q_e=ordi(a,ielemp) - m_e*asci( a,ielemp)
          else
            m_e=0
            q_e=0
          end if

          dl=sqrt((asci(a,ielemp)-asci(b,ielemp))**2 +
     *            (ordi(a,ielemp)-ordi(b,ielemp))**2 )*0.13D0

          print*,"------------------"
          print*,"RETTE X INTERSEZIONE"
          print*,"elemento, nodo",ielemp,j
          print*,"Retta faglia",m_r,q_r
          print*,"Retta bordo",m_e,q_e
          print*,"Punto iniziale",coord(sn1,1),coord(sn1,2)
          print*,"Zona di aggancio", dl
          print*,""

          xe_min=MIN(asci(a,ielemp),asci(b,ielemp))
          xe_max=MAX(asci(a,ielemp),asci(b,ielemp))
          ye_min=MIN(ordi(a,ielemp),ordi(b,ielemp))
          ye_max=MAX(ordi(a,ielemp),ordi(b,ielemp))

          if((m_r.eq.m_e).and.(q_e.eq.q_r))then
            !controllo nodi
            if(asci(a,ielemp).le.xmax.and.
     *         asci(a,ielemp).ge.xmin.and.
     *         ordi(a,ielemp).le.ymax.and.
     *         ordi(a,ielemp).ge.ymin) then
              xi=asci(a,ielemp)
              yi=ordi(a,ielemp)
              print*,"- rette coincidenti - intersezione in a",xi,yi
            else if(asci(b,ielemp).le.xmax.and.
     *         asci(b,ielemp).ge.xmin.and.
     *         ordi(b,ielemp).le.ymax.and.
     *         ordi(b,ielemp).ge.ymin) then
              xi=asci(b,ielemp)
              yi=ordi(b,ielemp)
              print*,"- rette coincidenti- intersezione in b",xi,yi
            end if
          else if(coord(sn1,1).eq.coord(sn2,1)) then            
            xi=coord(sn1,1)
            yi=m_e*xi+q_e
            print*," - faglia verticale",xi,yi
          else if(asci(a,ielemp).eq.asci(b,ielemp)) then
            xi=asci(a,ielemp)
            yi=m_r*xi+q_r
            print*," - bordo verticale",xi,yi
          else
            xi=(q_r-q_e)/(m_e-m_r)
            yi=m_e*xi+q_e
            print*," - rette incidenti",xi,yi
          end if

          if(xi.le.xmax) print*,"1- xi.le.xmax"
          if(xi.ge.xmin) print*,"2- xi.ge.xmin"
          if(xi.le.xe_max) print*,"3- xi.le.xe_max"
          if(xi.ge.xe_min) print*,"4- xi.ge.xe_min"
          if(yi.le.ymax) print*,"5- yi.le.ymax"
          if(yi.ge.ymin) print*,"6- yi.ge.ymin"
          if(yi.le.ye_max) print*,"7- yi.le.ye_max"
          if(yi.ge.ye_min) print*,"8- yi.ge.ye_min"
          if(.not.(xi.eq.coord(sn1,1).and.yi.eq.coord(sn1,2)))
     *     print*,"9-"
    
          if (xi.le.xmax.and.xi.ge.xmin.and.
     *          yi.le.ymax.and.yi.ge.ymin.and.
     *          xi.le.xe_max.and.xi.ge.xe_min.and.
     *          yi.le.ye_max.and.yi.ge.ye_min.and.
     *         .not.
     *           (xi.le.coord(sn1,1)+dl/3.and.
     *            xi.ge.coord(sn1,1)-dl/3.and.
     *            yi.le.coord(sn1,2)+dl/3.and.
     *            yi.ge.coord(sn1,2)-dl/3))then
              ielem=ielemp
c               ho trovato l'intersezione: ha coordinate xi yi e si trova sul bordo dell'elemento i
              go to 100 !esco dal check degli elementi a cui sn1 appartiene
c               verifico se il punto di intersezione coincide con un nodo preesistente
          end if
        end do
      end do

      print*,"! Something about the intersection went wrong"

100   continue
      print*,"(xi,yi)",xi,yi,a,j,b
      do inode=1,nnode
        if(asci(inode,ielem).le.xi+dl.and.
     *     asci(inode,ielem).ge.xi-dl.and.
     *     ordi(inode,ielem).le.yi+dl.and.
     *     ordi(inode,ielem).ge.yi-dl.and.
     *     lnods(ielem,inode).ne.sn1) then
c           coincide con un nodo preesistente, che diventa slip node
          print*,"coincide con nodo preesistente:",lnods(ielem,inode)
          control=1
          sn_new=lnods(ielem,inode)
          call shift(sn2,sn_new,snf,npoin,nsn)

c             richiamo la subroutine 1 elem per sn_new e sn1
          call taglio1elem(sn1,sn_new,snf,ielem,new_nelem,npoin,nsn,
     *           nnode,asci,ordi,ee,tt,nelemp_max,nelemp,ndime,
     *           lelemp,lnods,cf_fix,nvfix,nofix,ifpre,ndofn,
     *          presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
         go to 101
        end if
      end do


        !print*,"ielem",ielem,sn1
      if(control.ne.1) then
        print*,"Intersezione non coincide con un nodo preesistente"
c             non coincide con un nodo, xi yi diventa uno slip node

        do inode=1,nnode
          if(lnods(ielem,inode).eq.sn1) node1=inode          
        end do

        xmin=MIN(asci(a,ielem),asci(j,ielem))
        xmax=MAX(asci(a,ielem),asci(j,ielem))
        ymin=MIN(ordi(a,ielem),ordi(j,ielem))
        ymax=MAX(ordi(a,ielem),ordi(j,ielem))
        if(xi.le.xmax.and.xi.ge.xmin.and.
     *       yi.le.ymax.and.yi.ge.ymin) then
          node2=a
          node3=j
         else
          node2=j
          node3=b
         
        end if

200     continue 
      
        if (node1.eq.1.or.
     *       node1.eq.2.or.
     *       node1.eq.3) then
          !print*,"call corner-bordo"
          call corner_bordo(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)


        else
          if(node1.eq.4)then
            if(node2.eq.2.or.node2.eq.6) 
     *       call side_bordo_basso(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
            if(node2.eq.5.or.node2.eq.3) 
     *        call side_bordo_alto(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
      
          else if(node1.eq.5)then
            if(node2.eq.3.or.node2.eq.4)  
     *       call side_bordo_basso(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
            if(node2.eq.6.or.node2.eq.1)  
     *        call side_bordo_alto(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
          else if(node1.eq.6)then
            if(node2.eq.1.or.node2.eq.5)  
     *       call side_bordo_basso(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
            if(node2.eq.4.or.node2.eq.2) 
     *        call side_bordo_alto(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
          end if
        end if
      end if

101   continue

      nelem=new_nelem
       

      end

c--------------------------------------
      subroutine corner_side(ielem,node1,node2,sn2,snf,nsn,nelem,
     *     new_nelem,nnode,npoin,nelemp_max,lelemp,nelemp,asci,ordi,
     *       ndime,ee,tt,cf_fix,lnods,nvfix,
     *      nofix,
     *     ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      implicit none

      integer*4 i,j,k,ielem,jelem,nelem,nnode,npoin,
     *          node1,node2,nelemp_max,sn2,nsn,new_nelem

      integer*4 lnods(nelem,nnode),nelemp(npoin),
     *          lelemp(npoin,nelemp_max),snf(npoin),
     *          ee(nelem),tt(nelem)

      double precision x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,
     *                 k1x,k1y,k2x,k2y,k3x,k3y,cf_fix
     
      double precision asci(nnode,nelem),ordi(nnode,nelem)
      
      integer*4 ndime,nvfix,nvfix_t,ndofn,ndofn_T

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T)


      double precision presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

      print*,"npoin pre corner-side",npoin
      do j=1,nelemp(lnods(ielem,node2))
        if (lelemp(lnods(ielem,node2),j).ne.ielem) then
          jelem=lelemp(lnods(ielem,node2),j)
        end if
      end do
      

      x1=asci(1,ielem)
      y1=ordi(1,ielem) 
    
      x2=asci(2,ielem)
      y2=ordi(2,ielem)  
    
      x3=asci(3,ielem)
      y3=ordi(3,ielem)
      
      x4=asci(4,ielem)
      y4=ordi(4,ielem)  
    
      x5=asci(5,ielem)
      y5=ordi(5,ielem)  
    
      x6=asci(6,ielem)
      y6=ordi(6,ielem)
     
      ! taglio tra i nodi 1-5
      if(node1.eq.1) then      
        print*,"  Taglio 1-5"      
        k1x=(x2+x5)/2
        k1y=(y2+y5)/2
      
        k2x=(x5+x3)/2
        k2y=(y5+y3)/2
      
        k3x=(x1+x5)/2
        k3y=(y1+y5)/2   
                  
        write(2,1003) npoin+1,k1x,k1y,0.0D0,1.0D0
        write(2,1003) npoin+2,k2x,k2y,0.0D0,1.0D0
        write(2,1003) npoin+3,k3x,k3y,1.0D0,cf_fix


        !k3 diventa slip node
        !sn(npoin+3)=1

        call shift(sn2,npoin+3,snf,npoin,nsn)

        write(3,1005) ielem,1,
     *          lnods(ielem,2),lnods(ielem,5),lnods(ielem,1),
     *                    npoin+1,npoin+3,lnods(ielem,4)
        write(3,1005) new_nelem+1,1,
     *    lnods(ielem,5),lnods(ielem,3),lnods(ielem,1),
     *                      npoin+2,lnods(ielem,6),npoin+3
            
      ! taglio tra i nodi 2-6
      else if(node1.eq.2) then      
        print*,"  Taglio 2-6"       
        k1x=(x3+x6)/2
        k1y=(y3+y6)/2
      
        k2x=(x6+x1)/2
        k2y=(y6+y1)/2
      
        k3x=(x2+x6)/2
        k3y=(y2+y6)/2
                  

        write(2,1003) npoin+1,k1x,k1y,0.0D0,1.0D0
        write(2,1003) npoin+2,k2x,k2y,0.0D0,1.0D0
        write(2,1003) npoin+3,k3x,k3y,1.0D0,cf_fix 
    
        
        !k3 diventa slip node
        !sn(npoin+3)=1
        call shift(sn2,npoin+3,snf,npoin,nsn)

        write(3,1005) ielem,1,
     *            lnods(ielem,3),lnods(ielem,6),lnods(ielem,2),
     *                   npoin+1,npoin+3,lnods(ielem,5)
        write(3,1005) new_nelem+1,1,
     *         lnods(ielem,6),lnods(ielem,1),lnods(ielem,2),
     *          npoin+2,lnods(ielem,4),npoin+3


      !taglio tra i nodi 3-4
      else if(node1.eq.3) then      
        print*,"  Taglio 3-4"        
        k1x=(x1+x4)/2
        k1y=(y1+y4)/2
      
        k2x=(x4+x2)/2
        k2y=(y4+y2)/2
      
        k3x=(x3+x4)/2
        k3y=(y3+y4)/2
                  
        !aggiorno file  nodi
         write(2,1003) npoin+1,k1x,k1y,0.0D0,1.0D0
         write(2,1003) npoin+2,k2x,k2y,0.0D0,1.0D0
         write(2,1003) npoin+3,k3x,k3y,1.0D0,cf_fix
    
        !k3 diventa slip node
        !sn(npoin+3)=1
        call shift(sn2,npoin+3,snf,npoin,nsn)

        !aggiorno elementi
        write(3,1005) ielem,1,
     *       lnods(ielem,1),lnods(ielem,4),lnods(ielem,3),
     *        npoin+1,npoin+3,lnods(ielem,6)
        write(3,1005) new_nelem+1,1,
     *        lnods(ielem,4),lnods(ielem,2),lnods(ielem,3),
     *                    npoin+2,lnods(ielem,5),npoin+3
      end if    
    
      new_nelem=new_nelem+1
      npoin=npoin+3

c      call shapef_6(ndime,nnode,k1x,k1y,shape,deriv,
c     *             deriv2_s,deriv2_t)
c      do i=1,nnode
c        V1_0_new(npoin+1)=V1_0_new(npoin+1)+shape*V1_0(lnods(ielem,i))
c        V2_0_new(npoin+1)=V2_0_new(npoin+1)+shape*V2_0(lnods(ielem,i))
c        T_new(npoin+1)=T_new(npoin+1)+shape*T(lnods(ielem,i))
c      end do

c      call shapef_6(ndime,nnode,k2x,k2y,shape,deriv,
c     *             deriv2_s,deriv2_t)
c      do i=1,nnode
c        V1_0_new(npoin+2)=V1_0_new(npoin+2)+shape*V1_0(lnods(ielem,i))
c        V2_0_new(npoin+2)=V2_0_new(npoin+2)+shape*V2_0(lnods(ielem,i))
c        T_new(npoin+2)=T_new(npoin+2)+shape*T(lnods(ielem,i))
c      end do

c      call shapef_6(ndime,nnode,k3x,k3y,shape,deriv,
c     *             deriv2_s,deriv2_t)
c      do i=1,nnode
c        V1_0_new(npoin+3)=V1_0_new(npoin+3)+shape*V1_0(lnods(ielem,i))
c        V2_0_new(npoin+3)=V2_0_new(npoin+3)+shape*V2_0(lnods(ielem,i))
c        T_new(npoin+3)=T_new(npoin+3)+shape*T(lnods(ielem,i))
c      end do


      print*,"npoin_corner",npoin

      if(jelem.ne.0) call Raccordo_corner(jelem,
     *  lnods(ielem,node2),npoin-1,npoin-2,
     *          nelem,nnode,npoin,lnods,
     *           asci,ordi,new_nelem)
      

      do i=1,nelem    
        if(i.ne.ielem.and.i.ne.jelem) then
          write(3,1005) ee(i),tt(i),(lnods(i,k),k=1,nnode)   
        end if
      end do

      close(2)
      close(3)

      
      print*,""

      call AssemblaInput(new_nelem,npoin,ndime,nnode,
     *   nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)   
1005  format(10(1x,i10),5x,6(i2)) 
      
      end

c--------------------------------------
      subroutine side_side(node1,node2,ielem,nelem,new_nelem,npoin,
     *   lnods,nnode,nelemp,lelemp,nelemp_max,snf,ndime,asci,ordi,
     * nsn,sn2,ee,tt,cf_fix,nvfix,
     *      nofix,
     *     ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      implicit none 

      integer*4 i,k,ielem,nelem,jelem1,jelem2,jelem3,
     *          npoin,nnode,nelemp_max,node1,node2,
     *          sn2,nsn,new_nelem

      integer*4 lnods(nelem,nnode),lelemp(npoin,nelemp_max),
     *          nelemp(npoin),snf(npoin),
     *          ee(nelem),tt(nelem)

      double precision x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,
     *       k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,k5x,k5y,k6x,k6y,
     *       k7x,k7y,k8x,k8y,k9x,k9y,cf_fix
       
      double precision asci(nnode,nelem),ordi(nnode,nelem)

      integer*4 ndime,nvfix,nvfix_t,ndofn,ndofn_T

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T)


      double precision presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

      !print*,"  in side-side"
      x1=asci(1,ielem)
      y1=ordi(1,ielem) 
    
      x2=asci(2,ielem)
      y2=ordi(2,ielem)  
    
      x3=asci(3,ielem)
      y3=ordi(3,ielem)
      
      x4=asci(4,ielem)
      y4=ordi(4,ielem)  
    
      x5=asci(5,ielem)
      y5=ordi(5,ielem)  
    
      x6=asci(6,ielem)
      y6=ordi(6,ielem)

      
      jelem1=0
      jelem2=0
      jelem3=0

      do i=1,nelemp(lnods(ielem,4))
        if (lelemp(lnods(ielem,4),i).ne.ielem) then
          jelem1=lelemp(lnods(ielem,4),i)
        end if
      end do

      do i=1,nelemp(lnods(ielem,5))
        if (lelemp(lnods(ielem,5),i).ne.ielem) then
          jelem2=lelemp(lnods(ielem,5),i)
        end if
      end do

      do i=1,nelemp(lnods(ielem,6))
        if (lelemp(lnods(ielem,6),i).ne.ielem) then
          jelem3=lelemp(lnods(ielem,6),i)
        end if
      end do
c      CASO A   


      k1x=(x1+x4)/2
      k1y=(y1+y4)/2
    
      k2x=(x4+x2)/2
      k2y=(y4+y2)/2
    
      k3x=(x2+x5)/2
      k3y=(y2+y5)/2
      
      k4x=(x5+x3)/2
      k4y=(y5+y3)/2
    
      k5x=(x3+x6)/2
      k5y=(y3+y6)/2
    
      k6x=(x6+x1)/2
      k6y=(y6+y1)/2      

      k7x=(x6+x4)/2
      k7y=(y6+y4)/2      
      
      k8x=(x5+x4)/2
      k8y=(y5+y4)/2      
      
      k9x=(x5+x6)/2
      k9y=(y5+y6)/2      
              

      !k3 diventa slip node
      if(node1.eq.4.and.node2.eq.5) then
        print*,"  Taglio 4-5"
        write(2,1003) npoin+1,k1x,k1y,0.0D0,1.0D0
        write(2,1003) npoin+2,k2x,k2y,0.0D0,1.0D0
        write(2,1003) npoin+3,k3x,k3y,0.0D0,1.0D0
                   
        write(2,1003) npoin+4,k4x,k4y,0.0D0,1.0D0
        write(2,1003) npoin+5,k5x,k5y,0.0D0,1.0D0
        write(2,1003) npoin+6,k6x,k6y,0.0D0,1.0D0
         
        write(2,1003) npoin+7,k7x,k7y,0.0D0,1.0D0
        write(2,1003) npoin+8,k8x,k8y,1.0D0,cf_fix
        write(2,1003) npoin+9,k9x,k9y,0.0D0,1.0D0

        !sn(npoin+8)=1
        call shift(sn2,npoin+8,snf,npoin,nsn)
      else if(node1.eq.4.and.node2.eq.6) then
        print*,"  Taglio 4-6"
        write(2,1003) npoin+1,k1x,k1y,0.0D0,1.0D0
        write(2,1003) npoin+2,k2x,k2y,0.0D0,1.0D0
        write(2,1003) npoin+3,k3x,k3y,0.0D0,1.0D0
                   
        write(2,1003) npoin+4,k4x,k4y,0.0D0,1.0D0
        write(2,1003) npoin+5,k5x,k5y,0.0D0,1.0D0
        write(2,1003) npoin+6,k6x,k6y,0.0D0,1.0D0
          
        write(2,1003) npoin+7,k7x,k7y,1.0D0,cf_fix
        write(2,1003) npoin+8,k8x,k8y,0.0D0,1.0D0
        write(2,1003) npoin+9,k9x,k9y,0.0D0,1.0D0

        !sn(npoin+7)=1
        call shift(sn2,npoin+7,snf,npoin,nsn)
      else if(node1.eq.5.and.node2.eq.6) then
        print*,"  Taglio 5-6"
        write(2,1003) npoin+1,k1x,k1y,0.0D0,1.0D0 !npoin-8
        write(2,1003) npoin+2,k2x,k2y,0.0D0,1.0D0 !npoin-7
        write(2,1003) npoin+3,k3x,k3y,0.0D0,1.0D0 !npoin-6
                    
        write(2,1003) npoin+4,k4x,k4y,0.0D0,1.0D0 !npoin-5
        write(2,1003) npoin+5,k5x,k5y,0.0D0,1.0D0 !npoin-4
        write(2,1003) npoin+6,k6x,k6y,0.0D0,1.0D0 !npoin-3
         
        write(2,1003) npoin+7,k7x,k7y,0.0D0,1.0D0 !npoin-2
        write(2,1003) npoin+8,k8x,k8y,0.0D0,1.0D0 !npoin-1
        write(2,1003) npoin+9,k9x,k9y,1.0D0,cf_fix !npoin

        !sn(npoin+9)=1
        call shift(sn2,npoin+9,snf,npoin,nsn)
      end if

      write(3,1005) ielem,1,
     *                    lnods(ielem,6),lnods(ielem,5),lnods(ielem,3),
     *                          npoin+9,npoin+4,npoin+5
      write(3,1005) new_nelem+1,1,
     *           lnods(ielem,6),lnods(ielem,4),lnods(ielem,5),
     *                           npoin+7,npoin+8,npoin+9
      write(3,1005) new_nelem+2,1,
     *           lnods(ielem,1),lnods(ielem,4),lnods(ielem,6),
     *                          npoin+1,npoin+7,npoin+6
      write(3,1005) new_nelem+3,1,
     *          lnods(ielem,4),lnods(ielem,2),lnods(ielem,5),
     *                           npoin+2,npoin+3,npoin+8
      
       new_nelem=new_nelem+3
      npoin=npoin+9


      if(jelem1.ne.0)
     *  call Raccordo_corner(jelem1,lnods(ielem,4),npoin-7,npoin-8, !k2,k1
     *           nelem,nnode,npoin,lnods,
     *           asci,ordi,new_nelem)
      if(jelem2.ne.0)
     *  call Raccordo_corner(jelem2,lnods(ielem,5),npoin-6,npoin-7, !k4,k3
     *           nelem,nnode,npoin,lnods,
     *           asci,ordi,new_nelem)

      if(jelem3.ne.0)
     *  call Raccordo_corner(jelem3,lnods(ielem,6),npoin-5,npoin-6, !k6,k5
     *           nelem,nnode,npoin,lnods,
     *           asci,ordi,new_nelem)

      do i=1,nelem    
        if(i.ne.ielem.and.i.ne.jelem1.and.i.ne.jelem2.and.
     *   i.ne.jelem3) then
             write(3,1005) ee(i),tt(i),(lnods(i,k),k=1,nnode)   
        end if
      end do

      close(2)
      close(3)

      call AssemblaInput(new_nelem,npoin,ndime,nnode,
     *   nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      nelem=new_nelem
1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)     
1005  format(10(1x,i10),5x,6(i2)) 

      end

c--------------------------------------
      subroutine corner_corner(ielem,node1,node2,sn2,snf,
     *   nsn,nelem,new_nelem,nnode,npoin,ndime,ee,tt,
     *   lnods,nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
      
      implicit none

      integer*4 node1,node2,ielem,npoin,nsn,nelem,k,i,nnode,
     *          sn2,new_nelem

      integer*4 lnods(nelem,nnode),snf(npoin),ee(nelem),tt(nelem)

      !double precision sn(npoin),cf_fix,cf(npoin)

      integer*4 ndime,nvfix,nvfix_t,ndofn,ndofn_T

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T)

      double precision presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

      if(node1.eq.1.and.node2.eq.2) then
        !sn(lnods(ielem,4))=1
        !cf(lnods(ielem,4))=cf_fix
        call shift(sn2,lnods(ielem,4),snf,npoin,nsn)
      else if(node1.eq.2.and.node2.eq.3) then
        !sn(lnods(ielem,5))=1
        !cf(lnods(ielem,5))=cf_fix
        call shift(sn2,lnods(ielem,5),snf,npoin,nsn)
      else if(node1.eq.1.and.node2.eq.3) then
        !sn(lnods(ielem,6))=1
        !cf(lnods(ielem,6))=cf_fix
        call shift(sn2,lnods(ielem,6),snf,npoin,nsn)
      endif

      do i=1,nelem      
         write(3,1005) ee(i),tt(i),(lnods(i,k),k=1,nnode)   
      end do

      close(3)
      close(2)

       call AssemblaInput(new_nelem,npoin,ndime,nnode,
     *   nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

c1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)    
1005  format(10(1x,i10),5x,6(i2)) 

      end
c--------------------------------------
      subroutine corner_bordo(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      implicit none 

      integer*4 i,k,ielem,nelem,jelem1,jelem2,
     *          npoin,nnode,nelemp_max,node1,node2,
     *          nsn,new_nelem,sn_new,a,b,c,d,e,f,sn2,new_npoin,
     *          k1,k2,k4

      integer*4 lnods(nelem,nnode),lelemp(npoin,nelemp_max),
     *          nelemp(npoin),snf(npoin),
     *          ee(nelem),tt(nelem)

      double precision xa,xb,xc,xd,xe,xf,ya,yb,yc,yd,ye,yf,
     *       k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,k5x,k5y,
     *       cf_fix,xi,yi
       
      double precision asci(nnode,nelem),ordi(nnode,nelem)

      integer*4 ndime,nvfix,nvfix_t,ndofn,ndofn_T

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T)


      double precision presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

      print*,"  Corner-bordo", ielem
       
      new_npoin=npoin

      print*,"node1-node2",node1,node2


      if(node1.eq.1.and.node2.eq.2) then
        print*,"  1-2"
        a=1
        b=4
        c=2
        d=5
        e=3
        f=6

      else if(node1.eq.1.and.node2.eq.5) then
        print*,"  1-5"
        a=1
        b=6
        c=3
        d=5
        e=2
        f=4

      else if(node1.eq.2.and.node2.eq.3) then
        print*,"  2-3"
        a=2
        b=5
        c=3
        d=6
        e=1
        f=4

      else if(node1.eq.2.and.node2.eq.6) then
        print*,"  2-6"
        a=2
        b=4
        c=1
        d=6
        e=3
        f=5
      
      else if(node1.eq.3.and.node2.eq.1) then
        print*,"  3-1"
        a=3
        b=6
        c=1
        d=4
        e=2
        f=5

      else if(node1.eq.3.and.node2.eq.4) then
        print*,"  3-4"
        a=3
        b=5
        c=2
        d=4
        e=1
        f=6
      end if
      
      xa=asci(a,ielem)
      ya=ordi(a,ielem) 
    
      xb=asci(b,ielem)
      yb=ordi(b,ielem)  
    
      xc=asci(c,ielem)
      yc=ordi(c,ielem)
      
      xd=asci(d,ielem)
      yd=ordi(d,ielem)  
    
      xe=asci(e,ielem)
      ye=ordi(e,ielem)  
    
      xf=asci(f,ielem)
      yf=ordi(f,ielem)

      jelem1=0
      jelem2=0

      do i=1,nelemp(lnods(ielem,d))
        if (lelemp(lnods(ielem,d),i).ne.ielem) then
          jelem1=lelemp(lnods(ielem,d),i)
        end if
      end do

 
      
      k1x=(xi+xc)/2
      k1y=(yi+yc)/2
    
      k2x=(xi+xd)/2
      k2y=(yi+yd)/2
    
      k3x=(xi+xa)/2
      k3y=(yi+ya)/2
      
      k4x=(xe+xd)/2
      k4y=(ye+yd)/2
    
      k5x=(xa+xd)/2
      k5y=(ya+yd)/2
      
      !print*,"(x1,y1)",new_npoin+1,k1x,k1y
      !print*,"(x2,y2)",new_npoin+2,k2x,k2y
      !print*,"(x3,y3)",new_npoin+3,k3x,k3y
      !print*,"(x4,y4)",new_npoin+4,k4x,k4y
      !print*,"(x5,y5)",new_npoin+5,k5x,k5y
      !print*,"(x6,y6)",new_npoin+6,k6x,k6y
      !print*,"(x7,y7)",new_npoin+7,k7x,k7y
      !print*,"(x8,y8)",new_npoin+8,k8x,k8y
      !print*,"(xi,yi)",new_npoin+9,xi,yi
          
      write(2,1003)  new_npoin+1,k1x,k1y,0.0D0,1.0D0
      write(2,1003)  new_npoin+2,k2x,k2y,0.0D0,1.0D0
      write(2,1003)  new_npoin+3,k3x,k3y,1.0D0,cf_fix
                  
      write(2,1003)  new_npoin+4,k4x,k4y,0.0D0,1.0D0
      write(2,1003)  new_npoin+5,k5x,k5y,0.0D0,1.0D0
      write(2,1003)  new_npoin+6,xi,yi,1.0D0,cf_fix

      print 1003,  new_npoin+1,k1x,k1y,0.0D0,1.0D0
      print 1003,  new_npoin+2,k2x,k2y,0.0D0,1.0D0
      print 1003,  new_npoin+3,k3x,k3y,1.0D0,cf_fix
                  
      print 1003,  new_npoin+4,k4x,k4y,0.0D0,1.0D0
      print 1003,  new_npoin+5,k5x,k5y,0.0D0,1.0D0
      print 1003,  new_npoin+6,xi,yi,1.0D0,cf_fix

      k1=new_npoin+1
      k2=new_npoin+2
      k4=new_npoin+4

      sn_new= new_npoin+6
      call shift(sn2,sn_new,snf, new_npoin,nsn)
      call shift(sn_new, new_npoin+3,snf, new_npoin,nsn)


      if (node2.eq.2.or.node2.eq.3.or.node2.eq.1) then
        write(3,1005) ielem,1,
     *                  lnods(ielem,d),lnods(ielem,e),lnods(ielem,a),
     *                        new_npoin+4,lnods(ielem,f),new_npoin+5
        write(3,1005) new_nelem+1,1,
     *           sn_new,lnods(ielem,a),lnods(ielem,c),
     *                     new_npoin+3, lnods(ielem,b), new_npoin+1
        write(3,1005) new_nelem+2,1,
     *           lnods(ielem,d),lnods(ielem,a),sn_new,
     *                      new_npoin+5, new_npoin+3, new_npoin+2

        print 1005, ielem,1,
     *                  lnods(ielem,d),lnods(ielem,e),lnods(ielem,a),
     *                        new_npoin+4,lnods(ielem,f),new_npoin+5
       print 1005, new_nelem+1,1,
     *           sn_new,lnods(ielem,a),lnods(ielem,c),
     *                     new_npoin+3, lnods(ielem,b), new_npoin+1
        print 1005, new_nelem+2,1,
     *           lnods(ielem,d),lnods(ielem,a),sn_new,
     *                      new_npoin+5, new_npoin+3, new_npoin+2
      
      else if (node2.eq.5.or.node2.eq.6.or.node2.eq.4) then
        write(3,1005) ielem,1,
     *                  lnods(ielem,d),lnods(ielem,a),lnods(ielem,e),
     *                        new_npoin+5,lnods(ielem,f),new_npoin+4
        write(3,1005) new_nelem+1,1,
     *           sn_new,lnods(ielem,c),lnods(ielem,a),
     *                     new_npoin+1, lnods(ielem,b), new_npoin+3
        write(3,1005) new_nelem+2,1,
     *           lnods(ielem,d),sn_new,lnods(ielem,a),
     *                      new_npoin+2, new_npoin+3, new_npoin+5

        print 1005, ielem,1,
     *                  lnods(ielem,d),lnods(ielem,a),lnods(ielem,e),
     *                        new_npoin+4,lnods(ielem,f),new_npoin+5
        print 1005, new_nelem+1,1,
     *           sn_new,lnods(ielem,c),lnods(ielem,a),
     *                     new_npoin+1, lnods(ielem,b), new_npoin+3
        print 1005, new_nelem+2,1,
     *           lnods(ielem,d),sn_new,lnods(ielem,a),
     *                      new_npoin+2, new_npoin+3, new_npoin+5
      
      end if

      new_nelem=new_nelem+2
      new_npoin= new_npoin+6
       
      if(jelem1.ne.0) 
     *   call Raccordo_bordo(jelem1,jelem2,lnods(ielem,d),
     *      lnods(ielem,c),xi,yi,sn_new,k1,k2,k4,
     *      nelem,nnode,npoin,lnods,nelemp,
     *           lelemp,nelemp_max,asci,ordi,new_nelem,new_npoin)
 

      do i=1,nelem    
        if(i.ne.ielem.and.i.ne.jelem1.and.
     *     i.ne.jelem2) then
            write(3,1005) ee(i),tt(i),(lnods(i,k),k=1,nnode)   
        end if
      end do

      close(2)
      close(3)

      call AssemblaInput(new_nelem,new_npoin,ndime,nnode,
     *   nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      nelem=new_nelem
      npoin=new_npoin

1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)     
1005  format(10(1x,i10),5x,6(i2)) 

      end

c--------------------------------------
      subroutine side_bordo_basso(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      implicit none 

      integer*4 i,k,ielem,nelem,jelem1,jelem2,jelem3,
     *          npoin,nnode,nelemp_max,node1,node2,
     *          nsn,new_nelem,sn_new,a,b,c,d,e,f,sn2,new_npoin,
     *          k1,k2,k4,k7,k8

      integer*4 lnods(nelem,nnode),lelemp(npoin,nelemp_max),
     *          nelemp(npoin),snf(npoin),
     *          ee(nelem),tt(nelem)

      double precision xa,xb,xc,xd,xe,xf,ya,yb,yc,yd,ye,yf,
     *       k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,k5x,k5y,k6x,k6y,
     *       k7x,k7y,k8x,k8y,cf_fix,xi,yi
       
      double precision asci(nnode,nelem),ordi(nnode,nelem)

      integer*4 ndime,nvfix,nvfix_t,ndofn,ndofn_T

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T)


      double precision presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

      print*,"  Side-bordo basso",ielem
       
      new_npoin=npoin
      if(node1.eq.4.and.node2.eq.2) then
        print*,"  4-2"
        a=4
        b=2
        c=5
        d=3
        e=6
        f=1

      else if(node1.eq.4.and.node2.eq.6) then
        print*,"  4-6"
        a=4
        b=1
        c=6
        d=3
        e=5
        f=2

      else if(node1.eq.5.and.node2.eq.4) then
        print*,"  5-4"
        a=5
        b=2
        c=4
        d=1
        e=6
        f=3

      else if(node1.eq.5.and.node2.eq.3) then
        print*,"  5-3"
        a=5
        b=3
        c=6
        d=1
        e=4
        f=2
      
      else if(node1.eq.6.and.node2.eq.5) then
        print*,"  6-5"
        a=6
        b=3
        c=5
        d=2
        e=4
        f=1

      else if(node1.eq.6.and.node2.eq.1) then
        print*,"  6-1"
        a=6
        b=1
        c=4
        d=2
        e=5
        f=3
      end if
      
      xa=asci(a,ielem)
      ya=ordi(a,ielem) 
    
      xb=asci(b,ielem)
      yb=ordi(b,ielem)  
    
      xc=asci(c,ielem)
      yc=ordi(c,ielem)
      
      xd=asci(d,ielem)
      yd=ordi(d,ielem)  
    
      xe=asci(e,ielem)
      ye=ordi(e,ielem)  
    
      xf=asci(f,ielem)
      yf=ordi(f,ielem)

      jelem1=0
      jelem2=0
      jelem3=0

      do i=1,nelemp(lnods(ielem,a))
        if (lelemp(lnods(ielem,a),i).ne.ielem) then
          jelem1=lelemp(lnods(ielem,a),i)
        end if
      end do

      do i=1,nelemp(lnods(ielem,c))
        if (lelemp(lnods(ielem,c),i).ne.ielem) then
          jelem2=lelemp(lnods(ielem,c),i)
        end if
      end do
 
      
      k1x=(xi+xb)/2
      k1y=(yi+yb)/2
    
      k2x=(xi+xc)/2
      k2y=(yi+yc)/2
    
      k3x=(xi+xa)/2
      k3y=(yi+ya)/2
      
      k4x=(xc+xd)/2
      k4y=(yc+yd)/2
    
      k5x=(xc+xf)/2
      k5y=(yc+yf)/2
    
      k6x=(xc+xa)/2
      k6y=(yc+ya)/2      

      k7x=(xf+xa)/2
      k7y=(yf+ya)/2      
      
      k8x=(xa+xb)/2
      k8y=(ya+yb)/2 
      
      !print*,"(x1,y1)",new_npoin+1,k1x,k1y
      !print*,"(x2,y2)",new_npoin+2,k2x,k2y
      !print*,"(x3,y3)",new_npoin+3,k3x,k3y
      !print*,"(x4,y4)",new_npoin+4,k4x,k4y
      !print*,"(x5,y5)",new_npoin+5,k5x,k5y
      !print*,"(x6,y6)",new_npoin+6,k6x,k6y
      !print*,"(x7,y7)",new_npoin+7,k7x,k7y
      !print*,"(x8,y8)",new_npoin+8,k8x,k8y
      !print*,"(xi,yi)",new_npoin+9,xi,yi
          
      write(2,1003)  new_npoin+1,k1x,k1y,0.0D0,1.0D0
      write(2,1003)  new_npoin+2,k2x,k2y,0.0D0,1.0D0
      write(2,1003)  new_npoin+3,k3x,k3y,1.0D0,cf_fix
                  
      write(2,1003)  new_npoin+4,k4x,k4y,0.0D0,1.0D0
      write(2,1003)  new_npoin+5,k5x,k5y,0.0D0,1.0D0
      write(2,1003)  new_npoin+6,k6x,k6y,0.0D0,1.0D0
        
      write(2,1003)  new_npoin+7,k7x,k7y,0.0D0,1.0D0
      write(2,1003)  new_npoin+8,k8x,k8y,0.0D0,1.0D0
      write(2,1003)  new_npoin+9,xi,yi,1.0D0,cf_fix

      print 1003,  new_npoin+1,k1x,k1y,0.0D0,1.0D0
      print 1003,  new_npoin+2,k2x,k2y,0.0D0,1.0D0
      print 1003,  new_npoin+3,k3x,k3y,1.0D0,cf_fix
                  
      print 1003,  new_npoin+4,k4x,k4y,0.0D0,1.0D0
      print 1003,  new_npoin+5,k5x,k5y,0.0D0,1.0D0
      print 1003,  new_npoin+6,k6x,k6y,0.0D0,1.0D0
        
      print 1003,  new_npoin+7,k7x,k7y,0.0D0,1.0D0
      print 1003,  new_npoin+8,k8x,k8y,0.0D0,1.0D0
      print 1003,  new_npoin+9,xi,yi,1.0D0,cf_fix

      k1=new_npoin+1
      k2=new_npoin+2
      k4=new_npoin+4
      k7=new_npoin+7
      k8=new_npoin+8

      sn_new= new_npoin+9
      call shift(sn2,sn_new,snf, new_npoin,nsn)
      call shift(sn_new, new_npoin+3,snf, new_npoin,nsn)


      if(node2.eq.4.or.node2.eq.5.or.node2.eq.6) then
        write(3,1005) ielem,1,
     *                  lnods(ielem,c),lnods(ielem,d),lnods(ielem,f),
     *                        new_npoin+4,lnods(ielem,e),new_npoin+5
        write(3,1005) new_nelem+1,1,
     *           lnods(ielem,a),lnods(ielem,c),lnods(ielem,f),
     *                      new_npoin+6, new_npoin+5, new_npoin+7
        write(3,1005) new_nelem+2,1,
     *           sn_new,lnods(ielem,a),lnods(ielem,b),
     *                     new_npoin+3, new_npoin+8, new_npoin+1
        write(3,1005) new_nelem+3,1,
     *          sn_new,lnods(ielem,c),lnods(ielem,a),
     *                       new_npoin+2, new_npoin+6, new_npoin+3

        print 1005, ielem,1,
     *                  lnods(ielem,c),lnods(ielem,d),lnods(ielem,f),
     *                        new_npoin+4,lnods(ielem,e),new_npoin+5
        print 1005, new_nelem+1,1,
     *            lnods(ielem,a),lnods(ielem,c),lnods(ielem,f),
     *                      new_npoin+6, new_npoin+5, new_npoin+7
        print 1005, new_nelem+2,1,
     *           sn_new,lnods(ielem,a),lnods(ielem,b),
     *                     new_npoin+3, new_npoin+8, new_npoin+1
        print 1005, new_nelem+3,1,
     *          sn_new,lnods(ielem,c),lnods(ielem,a),
     *                       new_npoin+2, new_npoin+6, new_npoin+3
      else if(node2.eq.1.or.node2.eq.2.or.node2.eq.3) then
        write(3,1005) ielem,1,
     *                  lnods(ielem,c),lnods(ielem,f),lnods(ielem,d),
     *                        new_npoin+5,lnods(ielem,e),new_npoin+4
        write(3,1005) new_nelem+1,1,
     *           lnods(ielem,a),lnods(ielem,f),lnods(ielem,c),
     *                      new_npoin+7, new_npoin+5, new_npoin+6
        write(3,1005) new_nelem+2,1,
     *           sn_new,lnods(ielem,b),lnods(ielem,a),
     *                     new_npoin+1, new_npoin+8, new_npoin+3
        write(3,1005) new_nelem+3,1,
     *          sn_new,lnods(ielem,a),lnods(ielem,c),
     *                       new_npoin+3, new_npoin+6, new_npoin+2

        print 1005, ielem,1,
     *                  lnods(ielem,c),lnods(ielem,f),lnods(ielem,d),
     *                        new_npoin+5,lnods(ielem,e),new_npoin+4
        print 1005, new_nelem+1,1,
     *            lnods(ielem,a),lnods(ielem,f),lnods(ielem,c),
     *                      new_npoin+7, new_npoin+5, new_npoin+6
        print 1005, new_nelem+2,1,
     *           sn_new,lnods(ielem,b),lnods(ielem,a),
     *                     new_npoin+1, new_npoin+8, new_npoin+3
        print 1005, new_nelem+3,1,
     *          sn_new,lnods(ielem,a),lnods(ielem,c),
     *                       new_npoin+3, new_npoin+6, new_npoin+2
      end if
      new_nelem=new_nelem+3
      new_npoin= new_npoin+9
       
      if(jelem2.ne.0) then
        call Raccordo_bordo(jelem2,jelem3,lnods(ielem,c),
     *      lnods(ielem,b),xi,yi,sn_new,k1,k2,
     *           k4,nelem,nnode,npoin,lnods,nelemp,
     *           lelemp,nelemp_max,asci,ordi,new_nelem,new_npoin)

      end if
       

      if(jelem1.ne.0) then
        if(node2.eq.4.or.node2.eq.5.or.node2.eq.6) then
          call Raccordo_corner(jelem1,lnods(ielem,a),k7, 
     *           k8,nelem,nnode,new_npoin,lnods,
     *           asci,ordi,new_nelem)
        else
          call Raccordo_corner(jelem1,lnods(ielem,a),k8, 
     *           k7,nelem,nnode,new_npoin,lnods,
     *           asci,ordi,new_nelem)
        end if
      end if


      do i=1,nelem    
        if(i.ne.ielem.and.i.ne.jelem1.and.
     *     i.ne.jelem2.and.i.ne.jelem3) then
            write(3,1005) ee(i),tt(i),(lnods(i,k),k=1,nnode)   
        end if
      end do

      close(2)
      close(3)

      call AssemblaInput(new_nelem,new_npoin,ndime,nnode,
     *   nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      nelem=new_nelem
      npoin=new_npoin
1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)     
1005  format(10(1x,i10),5x,6(i2)) 

      end

c--------------------------------------
      subroutine side_bordo_alto(sn2,node1,node2,xi,yi,ielem,
     *   nelem,new_nelem,npoin,lnods,nnode,nelemp,lelemp,
     *   nelemp_max,snf,ndime,asci,
     *   ordi,nsn,ee,tt,cf_fix,nvfix,sn_new,
     *   nofix,
     *   ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      implicit none 

      integer*4 i,k,ielem,nelem,jelem1,jelem2,jelem3,
     *          npoin,nnode,nelemp_max,node1,node2,
     *          nsn,new_nelem,sn_new,a,b,c,d,e,f,sn2,new_npoin,
     *          k1,k2,k4,k7,k8

      integer*4 lnods(nelem,nnode),lelemp(npoin,nelemp_max),
     *          nelemp(npoin),snf(npoin),
     *          ee(nelem),tt(nelem)

      double precision xa,xb,xc,xd,xe,xf,ya,yb,yc,yd,ye,yf,
     *       k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,k5x,k5y,k6x,k6y,
     *       k7x,k7y,k8x,k8y,cf_fix,xi,yi
       
      double precision asci(nnode,nelem),ordi(nnode,nelem)

      integer*4 ndime,nvfix,nvfix_t,ndofn,ndofn_T

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T)


      double precision presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

      print*,"  Side-bordo alto",ielem
       
      new_npoin=npoin
      if(node1.eq.4.and.node2.eq.5) then
        print*,"  4-5"
        a=4
        b=2
        c=5
        d=3
        e=6
        f=1

      else if(node1.eq.4.and.node2.eq.3) then
        print*,"  4-3"
        a=4
        b=1
        c=6
        d=3
        e=5
        f=2

      else if(node1.eq.5.and.node2.eq.1) then
        print*,"  5-1"
        a=5
        b=2
        c=4
        d=1
        e=6
        f=3

      else if(node1.eq.5.and.node2.eq.6) then
        print*,"  5-6"
        a=5
        b=3
        c=6
        d=1
        e=4
        f=2
      
      else if(node1.eq.6.and.node2.eq.2) then
        print*,"  6-2"
        a=6
        b=3
        c=5
        d=2
        e=4
        f=1

      else if(node1.eq.6.and.node2.eq.4) then
        print*,"  6-4"
        a=6
        b=1
        c=4
        d=2
        e=5
        f=3
      end if
      
      xa=asci(a,ielem)
      ya=ordi(a,ielem) 
    
      xb=asci(b,ielem)
      yb=ordi(b,ielem)  
    
      xc=asci(c,ielem)
      yc=ordi(c,ielem)
      
      xd=asci(d,ielem)
      yd=ordi(d,ielem)  
    
      xe=asci(e,ielem)
      ye=ordi(e,ielem)  
    
      xf=asci(f,ielem)
      yf=ordi(f,ielem)

      jelem1=0
      jelem2=0
      jelem3=0


      do i=1,nelemp(lnods(ielem,a))
        if (lelemp(lnods(ielem,a),i).ne.ielem) then
          jelem1=lelemp(lnods(ielem,a),i)
        end if
      end do

      do i=1,nelemp(lnods(ielem,c))
        if (lelemp(lnods(ielem,c),i).ne.ielem) then
          jelem2=lelemp(lnods(ielem,c),i)
        end if
      end do
 
      
      k1x=(xi+xd)/2
      k1y=(yi+yd)/2
    
      k2x=(xi+xc)/2
      k2y=(yi+yc)/2
    
      k3x=(xi+xa)/2
      k3y=(yi+ya)/2
      
      k4x=(xc+xb)/2
      k4y=(yc+yb)/2
    
      k5x=(xi+xf)/2
      k5y=(yi+yf)/2
    
      k6x=(xc+xa)/2
      k6y=(yc+ya)/2      

      k7x=(xf+xa)/2
      k7y=(yf+ya)/2      
      
      k8x=(xa+xb)/2
      k8y=(ya+yb)/2 
      
      !print*,"(x1,y1)",new_npoin+1,k1x,k1y
      !print*,"(x2,y2)",new_npoin+2,k2x,k2y
      !print*,"(x3,y3)",new_npoin+3,k3x,k3y
      !print*,"(x4,y4)",new_npoin+4,k4x,k4y
      !print*,"(x5,y5)",new_npoin+5,k5x,k5y
      !print*,"(x6,y6)",new_npoin+6,k6x,k6y
      !print*,"(x7,y7)",new_npoin+7,k7x,k7y
      !print*,"(x8,y8)",new_npoin+8,k8x,k8y
      !print*,"(xi,yi)",new_npoin+9,xi,yi
          
      write(2,1003)  new_npoin+1,k1x,k1y,0.0D0,1.0D0
      write(2,1003)  new_npoin+2,k2x,k2y,0.0D0,1.0D0
      write(2,1003)  new_npoin+3,k3x,k3y,1.0D0,cf_fix
                  
      write(2,1003)  new_npoin+4,k4x,k4y,0.0D0,1.0D0
      write(2,1003)  new_npoin+5,k5x,k5y,0.0D0,1.0D0
      write(2,1003)  new_npoin+6,k6x,k6y,0.0D0,1.0D0
        
      write(2,1003)  new_npoin+7,k7x,k7y,0.0D0,1.0D0
      write(2,1003)  new_npoin+8,k8x,k8y,0.0D0,1.0D0
      write(2,1003)  new_npoin+9,xi,yi,1.0D0,cf_fix

      print 1003,  new_npoin+1,k1x,k1y,0.0D0,1.0D0
      print 1003,  new_npoin+2,k2x,k2y,0.0D0,1.0D0
      print 1003,  new_npoin+3,k3x,k3y,1.0D0,cf_fix
                  
      print 1003,  new_npoin+4,k4x,k4y,0.0D0,1.0D0
      print 1003,  new_npoin+5,k5x,k5y,0.0D0,1.0D0
      print 1003,  new_npoin+6,k6x,k6y,0.0D0,1.0D0
        
      print 1003,  new_npoin+7,k7x,k7y,0.0D0,1.0D0
      print 1003,  new_npoin+8,k8x,k8y,0.0D0,1.0D0
      print 1003,  new_npoin+9,xi,yi,1.0D0,cf_fix

      k1=new_npoin+1
      k2=new_npoin+2
      k4=new_npoin+4
      k7=new_npoin+7
      k8=new_npoin+8

      sn_new= new_npoin+9
      call shift(sn2,sn_new,snf, new_npoin,nsn)
      call shift(sn_new, new_npoin+3,snf, new_npoin,nsn)

      if(node2.eq.1.or.node2.eq.2.or.node2.eq.3) then
        write(3,1005) ielem,1,
     *                  sn_new,lnods(ielem,d),lnods(ielem,f),
     *                        new_npoin+1,lnods(ielem,e),new_npoin+5
        write(3,1005) new_nelem+1,1,
     *           lnods(ielem,a),sn_new,lnods(ielem,f),
     *                      new_npoin+3,new_npoin+5, new_npoin+7
        write(3,1005) new_nelem+2,1,
     *          lnods(ielem,a),lnods(ielem,b),lnods(ielem,c),
     *                     new_npoin+8, new_npoin+4, new_npoin+6
        write(3,1005) new_nelem+3,1,
     *          sn_new,lnods(ielem,a),lnods(ielem,c),
     *                       new_npoin+3, new_npoin+6, new_npoin+2

        print 1005, ielem,1,
     *                  sn_new,lnods(ielem,d),lnods(ielem,f),
     *                        new_npoin+1,lnods(ielem,e),new_npoin+5
        print 1005, new_nelem+1,1,
     *           lnods(ielem,a),sn_new,lnods(ielem,f),
     *                      new_npoin+3,new_npoin+5, new_npoin+7
        print 1005, new_nelem+2,1,
     *          lnods(ielem,a),lnods(ielem,b),lnods(ielem,c),
     *                     new_npoin+8, new_npoin+4, new_npoin+6
        print 1005, new_nelem+3,1,
     *          sn_new,lnods(ielem,a),lnods(ielem,c),
     *                       new_npoin+3, new_npoin+6, new_npoin+2
      else if(node2.eq.4.or.node2.eq.5.or.node2.eq.6) then
        write(3,1005) ielem,1,
     *                  sn_new,lnods(ielem,f),lnods(ielem,d),
     *                        new_npoin+5,lnods(ielem,e),new_npoin+1
        write(3,1005) new_nelem+1,1,
     *           lnods(ielem,a),lnods(ielem,f),sn_new,
     *                      new_npoin+7,new_npoin+5, new_npoin+3
        write(3,1005) new_nelem+2,1,
     *          lnods(ielem,a),lnods(ielem,c),lnods(ielem,b),
     *                     new_npoin+6, new_npoin+4, new_npoin+8
        write(3,1005) new_nelem+3,1,
     *          sn_new,lnods(ielem,c),lnods(ielem,a),
     *                       new_npoin+2, new_npoin+6, new_npoin+3

        print 1005, ielem,1,
     *                  sn_new,lnods(ielem,f),lnods(ielem,d),
     *                        new_npoin+5,lnods(ielem,e),new_npoin+1
        print 1005, new_nelem+1,1,
     *           lnods(ielem,a),lnods(ielem,f),sn_new,
     *                      new_npoin+7,new_npoin+5, new_npoin+3
        print 1005, new_nelem+2,1,
     *          lnods(ielem,a),lnods(ielem,c),lnods(ielem,b),
     *                     new_npoin+6, new_npoin+4, new_npoin+8
        print 1005, new_nelem+3,1,
     *          sn_new,lnods(ielem,c),lnods(ielem,a),
     *                       new_npoin+2, new_npoin+6, new_npoin+3
      end if
      new_nelem=new_nelem+3
      new_npoin= new_npoin+9
       
      if(jelem2.ne.0) 
     *   call Raccordo_bordo(jelem2,jelem3,lnods(ielem,c),
     *      lnods(ielem,d),xi,yi,sn_new,k1,k2,
     *           k4,nelem,nnode,npoin,lnods,nelemp,
     *           lelemp,nelemp_max,asci,ordi,new_nelem,new_npoin)

       
      
       

      if(jelem1.ne.0) then
        if(node2.eq.4.or.node2.eq.5.or.node2.eq.6) then
          call Raccordo_corner(jelem1,lnods(ielem,a),k8, 
     *           k7,nelem,nnode, new_npoin,lnods,
     *           asci,ordi,new_nelem)
        else
          call Raccordo_corner(jelem1,lnods(ielem,a),k8, 
     *           k7,nelem,nnode,new_npoin,lnods,
     *           asci,ordi,new_nelem)
        end if
      end if

      do i=1,nelem    
        if(i.ne.ielem.and.i.ne.jelem1.and.
     *     i.ne.jelem2.and.i.ne.jelem3) then
            write(3,1005) ee(i),tt(i),(lnods(i,k),k=1,nnode)   
        end if
      end do

      close(2)
      close(3)

      call AssemblaInput(new_nelem,new_npoin,ndime,nnode,
     *   nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)

      nelem=new_nelem
      npoin=new_npoin
1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)     
1005  format(10(1x,i10),5x,6(i2)) 

      end      


c--------------------------------------
      subroutine Raccordo_corner(ielem,poin,k1,k2,
     *           nelem,nnode,npoin,lnods,
     *           asci,ordi,new_nelem)
      implicit none
      integer*4 ielem,nelem,nnode,npoin,poin,
     *          node,k1,k2,inode,new_nelem

      integer*4 lnods(nelem,nnode)

      double precision x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,
     *                 k3x,k3y

      double precision asci(nnode,nelem),ordi(nnode,nelem)
c     CASO B  del refinement quadratico
      
      print*,""
      print*,"    Raccordo corner",ielem,poin
       

      x1=asci(1,ielem)
      y1=ordi(1,ielem) 
    
      x2=asci(2,ielem)
      y2=ordi(2,ielem)  
    
      x3=asci(3,ielem)
      y3=ordi(3,ielem)
      
      x4=asci(4,ielem)
      y4=ordi(4,ielem)  
    
      x5=asci(5,ielem)
      y5=ordi(5,ielem)  
    
      x6=asci(6,ielem)
      y6=ordi(6,ielem)

      do inode=1,nnode
        if(lnods(ielem,inode).eq.poin) node=inode
      end do
      
      !nodi 3-4
      if(node.eq.3.or.node.eq.4) then      
              print*,"    3 o 4"
        k3x=(x3+x4)/2
        k3y=(y3+y4)/2
                   
        !aggiorno file  nodi
        write(2,1003) npoin+1,k3x,k3y,0.0D0,1.0D0

        print 1003, npoin+1,k3x,k3y,0.0D0,1.0D0
        !aggiorno elementi
        write(3,1005) ielem,1,
     *         lnods(ielem,1),lnods(ielem,4),lnods(ielem,3),
     *          k1,npoin+1,lnods(ielem,6)
        write(3,1005) new_nelem+1,1,
     *            lnods(ielem,4),lnods(ielem,2),lnods(ielem,3),
     *           k2,lnods(ielem,5),npoin+1

        print 1005, ielem,1,
     *         lnods(ielem,1),lnods(ielem,4),lnods(ielem,3),
     *          k1,npoin+1,lnods(ielem,6)
        print 1005, new_nelem+1,1,
     *            lnods(ielem,4),lnods(ielem,2),lnods(ielem,3),
     *           k2,lnods(ielem,5),npoin+1
       print*,"npoin",npoin
      new_nelem=new_nelem+1
        npoin=npoin+1  
        print*,"new_npoin",npoin

      !nodi 1-5
      else if(node.eq.1.or.node.eq.5) then     
        print*,"    1 o 5"
        k3x=(x1+x5)/2
        k3y=(y1+y5)/2
                   
         write(2,1003) npoin+1,k3x,k3y,0.0D0,1.0D0    
         print 1003, npoin+1,k3x,k3y,0.0D0,1.0D0
    
         write(3,1005) ielem,1,
     *          lnods(ielem,2),lnods(ielem,5),lnods(ielem,1),
     *            k1,npoin+1,lnods(ielem,4)
         write(3,1005) new_nelem+1,1,
     *             lnods(ielem,5),lnods(ielem,3),lnods(ielem,1),
     *              k2,lnods(ielem,6),npoin+1
     
         print 1005, ielem,1,
     *          lnods(ielem,2),lnods(ielem,5),lnods(ielem,1),
     *            k1,npoin+1,lnods(ielem,4)
         print 1005, new_nelem+1,1,
     *             lnods(ielem,5),lnods(ielem,3),lnods(ielem,1),
     *              k2,lnods(ielem,6),npoin+1
        new_nelem=new_nelem+1
        npoin=npoin+1  
      !nodi 2-6
      else if(node.eq.2.or.node.eq.6) then      
        print*,"    2 o 6"
        k3x=(x2+x6)/2
        k3y=(y2+y6)/2
        
         write(2,1003) npoin+1,k3x,k3y,0.0D0,1.0D0
         print 1003, npoin+1,k3x,k3y,0.0D0,1.0D0

         write(3,1005) ielem,1,
     *       lnods(ielem,3),lnods(ielem,6),lnods(ielem,2),
     *                          k1,npoin+1,lnods(ielem,5)
         write(3,1005) new_nelem+1,1,
     *            lnods(ielem,6),lnods(ielem,1),lnods(ielem,2),
     *                     k2,lnods(ielem,4),npoin+1

        print 1005, ielem,1,
     *       lnods(ielem,3),lnods(ielem,6),lnods(ielem,2),
     *                          k1,npoin+1,lnods(ielem,5)
         print 1005, new_nelem+1,1,
     *            lnods(ielem,6),lnods(ielem,1),lnods(ielem,2),
     *                     k2,lnods(ielem,4),npoin+1
        new_nelem=new_nelem+1
        npoin=npoin+1
      endif

1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)     
1005  format(10(1x,i10),5x,6(i2))       
      end

c----------------------------------------

      subroutine Raccordo_bordo(ielem,jelem,poin1,poin2,xi,yi,sn_new,
     *           k1,k2,k4,nelem,nnode,npoin,lnods,nelemp,lelemp,
     *           nelemp_max,asci,ordi,new_nelem,new_npoin)
      implicit none 

      integer*4 i,k1,k2,k4,ielem,nelem,jelem,inode,poin1,poin2,
     *          npoin,nnode,nelemp_max,node1,node2,new_nelem,
     *          sn_new,a,b,c,d,e,f,new_npoin

      integer*4 lnods(nelem,nnode),lelemp(npoin,nelemp_max),
     *          nelemp(npoin),snf(npoin),
     *          ee(nelem),tt(nelem)

      double precision xa,xb,xc,xd,xe,xf,ya,yb,yc,yd,ye,yf,
     *       k3x,k3y,k5x,k5y,k6x,k6y,
     *       k7x,k7y,k8x,k8y,xi,yi
       
      double precision asci(nnode,nelem),ordi(nnode,nelem)

      integer*4 ndime

      double precision Uo

      print*,"    Raccordo bordo", ielem
       
      print*,"poin1,poin2",poin1,poin2

c      do i=1,nelem
c        print*,i,(lnods(i,inode),inode=1,6)
c      end do

      do inode=1,nnode
        if(lnods(ielem,inode).eq.poin1) node1=inode
        if(lnods(ielem,inode).eq.poin2) node2=inode
      end do

      print*,"node1,node2",node1,node2



      if(node1.eq.4.and.node2.eq.1) then
        print*,"    4-1"
        a=1
        b=4
        c=2
        d=5
        e=3
        f=6

      else if(node1.eq.4.and.node2.eq.2) then
        print*,"    4-2"
        a=2
        b=4
        c=1
        d=6
        e=3
        f=5

      else if(node1.eq.5.and.node2.eq.2) then
        print*,"    5-2"
        a=2
        b=5
        c=3
        d=6
        e=1
        f=4

      else if(node1.eq.5.and.node2.eq.3) then
        print*,"    5-3"
        a=3
        b=5
        c=2
        d=4
        e=1
        f=6
      
      else if(node1.eq.6.and.node2.eq.1) then
        print*,"    6-1"
        a=1
        b=6
        c=3
        d=5
        e=2
        f=4

      else if(node1.eq.6.and.node2.eq.3) then
        print*,"    6-3"
        a=3
        b=6
        c=1
        d=4
        e=2
        f=5
      end if
      
      xa=asci(a,ielem)
      ya=ordi(a,ielem) 
    
      xb=asci(b,ielem)
      yb=ordi(b,ielem)  
    
      xc=asci(c,ielem)
      yc=ordi(c,ielem)
      
      xd=asci(d,ielem)
      yd=ordi(d,ielem)  
    
      xe=asci(e,ielem)
      ye=ordi(e,ielem)  
    
      xf=asci(f,ielem)
      yf=ordi(f,ielem)


      do i=1,nelemp(lnods(ielem,f))
        if (lelemp(lnods(ielem,f),i).ne.ielem) then
          jelem=lelemp(lnods(ielem,f),i)
        end if
      end do

      print*,"jelem",jelem
      

      k3x=(xi+xf)/2
      k3y=(yi+yf)/2
    
      k5x=(xb+xe)/2
      k5y=(yb+ye)/2
    
      k6x=(xb+xf)/2
      k6y=(yb+yf)/2      

      k7x=(xe+xf)/2
      k7y=(ye+yf)/2      
      
      k8x=(xf+xa)/2
      k8y=(yf+ya)/2      
          
      
      !print*,"(x3,y3)",new_npoin+1,k3x,k3y
      !print*,"(x5,y5)",new_npoin+2,k5x,k5y
      !print*,"(x6,y6)",new_npoin+3,k6x,k6y
      !print*,"(x7,y7)",new_npoin+4,k7x,k7y
      !print*,"(x8,y8)",new_npoin+5,k8x,k8y
      !print*,"(xi,yi)",sn_new,xi,yi
          

      write(2,1003) new_npoin+1,k3x,k3y,0.0D0,1.0D0
                  
      write(2,1003) new_npoin+2,k5x,k5y,0.0D0,1.0D0
      write(2,1003) new_npoin+3,k6x,k6y,0.0D0,1.0D0
        
      write(2,1003) new_npoin+4,k7x,k7y,0.0D0,1.0D0
      write(2,1003) new_npoin+5,k8x,k8y,0.0D0,1.0D0


      print 1003, new_npoin+1,k3x,k3y,0.0D0,1.0D0
                  
      print 1003, new_npoin+2,k5x,k5y,0.0D0,1.0D0
      print 1003, new_npoin+3,k6x,k6y,0.0D0,1.0D0
        
      print 1003, new_npoin+4,k7x,k7y,0.0D0,1.0D0
      print 1003, new_npoin+5,k8x,k8y,0.0D0,1.0D0
      

      print*,"a",lnods(ielem,a)
      print*,"b",lnods(ielem,b)
      print*,"c",lnods(ielem,c)
      print*,"d",lnods(ielem,d)
      print*,"e",lnods(ielem,e)
      print*,"f",lnods(ielem,f)


      if((node1.eq.6.and.node2.eq.1).or.
     *   (node1.eq.4.and.node2.eq.2).or.
     *   (node1.eq.5.and.node2.eq.3)) then
        write(3,1005) ielem,1,
     *             lnods(ielem,c),lnods(ielem,b),lnods(ielem,e),
     *                        k4,new_npoin+2,lnods(ielem,d)
        write(3,1005) new_nelem+1,1,
     *           lnods(ielem,b),lnods(ielem,f),lnods(ielem,e),
     *                     new_npoin+3,new_npoin+4,new_npoin+2
        write(3,1005) new_nelem+2,1,
     *           sn_new,lnods(ielem,f),lnods(ielem,b),
     *                      new_npoin+1,new_npoin+3,k2
        write(3,1005) new_nelem+3,1,
     *          sn_new,lnods(ielem,a),lnods(ielem,f),
     *                         k1,new_npoin+5,new_npoin+1

        print 1005, ielem,1,
     *             lnods(ielem,c),lnods(ielem,b),lnods(ielem,e),
     *                        k4,new_npoin+2,lnods(ielem,d)
        print 1005, new_nelem+1,1,
     *           lnods(ielem,b),lnods(ielem,f),lnods(ielem,e),
     *                     new_npoin+3,new_npoin+4,new_npoin+2
        print 1005, new_nelem+2,1,
     *           sn_new,lnods(ielem,f),lnods(ielem,b),
     *                      new_npoin+1,new_npoin+3,k2
        print 1005, new_nelem+3,1,
     *          sn_new,lnods(ielem,a),lnods(ielem,f),
     *                         k1,new_npoin+5,new_npoin+1
      else
        write(3,1005) ielem,1,
     *             lnods(ielem,c),lnods(ielem,e),lnods(ielem,b),
     *                        lnods(ielem,d),new_npoin+2,k4
        write(3,1005) new_nelem+1,1,
     *           lnods(ielem,b),lnods(ielem,e),lnods(ielem,f),
     *                     new_npoin+2,new_npoin+4,new_npoin+3
        write(3,1005) new_nelem+2,1,
     *           sn_new,lnods(ielem,b),lnods(ielem,f),
     *                      k2,new_npoin+3,new_npoin+1
        write(3,1005) new_nelem+3,1,
     *          sn_new,lnods(ielem,f),lnods(ielem,a),
     *                         new_npoin+1,new_npoin+5,k1

        print 1005, ielem,1,
     *             lnods(ielem,c),lnods(ielem,e),lnods(ielem,b),
     *                        lnods(ielem,d),new_npoin+2,k4
        print 1005, new_nelem+1,1,
     *           lnods(ielem,b),lnods(ielem,e),lnods(ielem,f),
     *                     new_npoin+2,new_npoin+4,new_npoin+3
        print 1005, new_nelem+2,1,
     *           sn_new,lnods(ielem,b),lnods(ielem,f),
     *                      k2,new_npoin+3,new_npoin+1
        print 1005, new_nelem+3,1,
     *          sn_new,lnods(ielem,f),lnods(ielem,a),
     *                         new_npoin+1,new_npoin+5,k1
      end if
      new_nelem=new_nelem+3
      new_npoin=new_npoin+5

      !print*,lnods(ielem,a)

      if(jelem.ne.0) then
        if((node1.eq.6.and.node2.eq.1).or.
     *   (node1.eq.4.and.node2.eq.2).or.
     *   (node1.eq.5.and.node2.eq.3)) then
          call Raccordo_corner(jelem,lnods(ielem,f),new_npoin-1,
     *          new_npoin,nelem,nnode,new_npoin,lnods,
     *           asci,ordi,new_nelem)

        else
          call Raccordo_corner(jelem,lnods(ielem,f),new_npoin,
     *          new_npoin-1,nelem,nnode,new_npoin,lnods,
     *           asci,ordi,new_nelem)
        end if
      end if
       

!      close(2)
!      close(3)

1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)    
1005  format(10(1x,i10),5x,6(i2)) 

            
      end
c--------------------------------------
      subroutine shift(sn2,sn_new,snf,npoin,nsn)
      implicit none

      integer*4 sn2,sn_new,npoin,k,i,nsn,j
      integer*4 snf(npoin)


      do i=1,nsn
        if(snf(i).eq.sn2) k=i
      end do

      do i=1,nsn+1-k
        j=nsn+1-i
        snf(j+1)=snf(j)
      end do
      snf(k)=sn_new
      nsn=nsn+1

      end

c--------------------------------------

      subroutine AssemblaInput(nelem,npoin,ndime,nnode,
     *   nvfix,nofix,ifpre,ndofn,presc,Uo,
     *   nvfix_T,nofix_T,ifpre_T,presc_T)
      implicit none

      integer*4 nelem,npoin,ndime,nvfix,nvfix_t,ndofn,ndofn_T,
     *          nnode,ivfix,ivfix_T,ipoin,inode,ielem,idofn,i

      integer*4 nofix(nvfix),ifpre(nvfix,ndofn),
     *          nofix_T(nvfix_T),ifpre_T(nvfix_T,ndofn_T),
     *         lnods(nelem,nnode),tip(nelem)


      double precision coord(npoin,ndime),sn(npoin),cf(npoin),
     *                 presc(nvfix,ndofn),presc_T(nvfix_T,ndofn_T)

      double precision Uo

      print*,""
      print*,"-- Assemblo nuovo input"
      
      print*,"nelem",nelem
      !print*,"new_nelem",new_nelem
      print*,"npoin",npoin

      !sn=0.0D0
      !cf=0.0D0

c     Apertura file di input
      Open(unit=1,file='fort.14',status='unknown')
      close(unit=1,status='delete')
      Open(unit=1,file='fort.14',status='new')  

c     Legge le coordinate nodali e le scrive nel file "Nodi"
      
      Open(unit=2,file='Nodi.tmp',status='old')
      Open(unit=3,file='Elementi.tmp',status='old')

      write(1,1009) "aaaaaaaaaaaa"
      
c      print*,nelem
c      print*,npoin
c      print*,ndime

      do i=1,npoin
        read(2,1003) ipoin,coord(ipoin,1),coord(ipoin,2),
     *   sn(ipoin),cf(ipoin)
!         print*,ipoin
      end do

      do ipoin=1,npoin
        write(1,1003) ipoin,coord(ipoin,1),coord(ipoin,2), ! Lettura coordinate Nodali
     *             sn(ipoin),cf(ipoin)
      enddo

      write(1,1009) "aaaaaaaaaaaa"

      do i=1,nelem! Calcolo proprietà geometriche dell'elemento
        read(3,1005) ielem,tip(ielem),(lnods(ielem,inode),
     *            inode=1,nnode)
      end do
      
      do ielem=1,nelem
        write(1,1005) ielem,tip(ielem),(lnods(ielem,inode), ! Lettura: Numero Elemento, Tipo di Elemento, Connessioni nodali per ogni elemento
     *            inode=1,nnode)
      enddo
    
                       
      write(1,1009) "aaaaaaaaaaaa"

      do ivfix=1,nvfix 
        write(1,*) nofix(ivfix),(ifpre(ivfix,idofn), idofn=1,ndofn),
     *               (presc(ivfix,idofn)*Uo, idofn=1,ndofn)
      enddo  


      write(1,1009) "aaaaaaaaaaaa"
      do ivfix_T=1,nvfix_T  ! Temperatura normalizzata
        write(1,*) nofix_T(ivfix_T),ifpre_T(ivfix_T,1),
     *      presc_T(ivfix_T,1)    
    
      enddo        

      !print*,"assemblaggio finito"
      close(1)
      close(2)
      close(3)

c      call neighboring_gridbounded(nelem,npoin,nnode,ndime,
c     *           nelemp_max,coord,lnods,nelemp,lelemp,
c     *           nlpoins,lpoins) 

1009  format(a12)
1003  format(i10,2(1x,f10.5),3x,f3.1,1x,f6.2)     
1005  format(10(1x,i10),5x,6(i2)) 

      end



c--------------------------------------------------------------

      subroutine neighboring_gridbounded(nelem,npoin,nnode,ndime,
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


c--------------------------------------------------------------

      subroutine neighboring_free(nelem,npoin,nnode,ndime,
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
 
c       nlpoins,lpoins
        
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



c--------------------------------------------------------------
      subroutine RicalcoloLelemp(nelem,npoin,nnode,
     *           nelemp_max,lnods,nelemp,lelemp)
      
      implicit none

      integer*4 npoin,ipoin,i,ielem,nelem,nnode,
     *            inode,k,nelemp_max
      
      Integer*4 lnods(nelem,nnode),nelemp(npoin)
        

      integer*4 lelemp(npoin,nelemp_max),lpoins(npoin,nelemp_max)
        
        

      open(unit=21,file='Corrispondenze_Nodi2Elementi',status="unknown")
      close(unit=21,status='delete')
      open(unit=21,file='Corrispondenze_Nodi2Elementi',status="new")
  
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
        
        
      write(21,*) ipoin,nelemp(ipoin),
     *       (lelemp(ipoin,i),i=1,nelemp(ipoin))
      
        close(unit=21)
      enddo 

      end
