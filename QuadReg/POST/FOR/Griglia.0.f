c      Griglia.f

c      Anna Maria Marotta
c     20/07/2004

c      Programma per leggere i dati e plottare la Griglia

C     PROGRAMMA PER LEGGERE IL FILE INPUT E GRAFICARE LA MESH

      implicit none
      
      Integer*4 i,j,k,nnode,nodi,elementi
      
      parameter(nnode=6,
     *          nodi=81,
     *          elementi=32)
         
      double precision xx(nodi),yy(nodi)      
      double precision asci(nnode,elementi),ordi(nnode,elementi)
      double precision asci_M(elementi),ordi_M(elementi)
      
      Integer*4 cc(nnode,elementi),ee(elementi),tt(elementi),nn(nodi)
      character*100 testo
                             
c     Apertura file di input
      Open(unit=1,file='Input.tmp',status='old')
      
c     Legge le coordinate nodali e le scrive nel file "Nodi"
      
      Open(unit=2,file='Nodi.tmp',status='unknown')
      close(unit=2,status='delete')
      Open(unit=2,file='Nodi.tmp',status='new')
      
      read(1,*) testo

      do i=1,nodi
        read(1,*) nn(i),xx(i),yy(i)
        write(2,*) nn(i),xx(i),yy(i)
      enddo
      close(unit=2)
       
c     Legge le connessioni nodali c(k,j) e il tipo di materiale t(j)
      
       read(1,*) testo
       do j=1,elementi
         read(1,*) ee(j),tt(j),(cc(k,j),k=1,nnode)
       enddo
       
       close(unit=1) 
       
c.....Per ogni nodo di connessione sostinuisce il numero del nodo 
c     con le coordinate nodali

       do 300 j=1,elementi
         do 400 k=1,nnode
           do 500 i=1,nodi
             if(cc(k,j).eq.nn(i)) then
               asci(k,j)=xx(i)
               ordi(k,j)=yy(i)               
             endif
500         continue
400       continue        
300     continue         

      Open(unit=11,file='Elem-1.tmp', status='unknown')
      Open(unit=12,file='Elem-2.tmp', status='unknown')
      Open(unit=13,file='Elem-3.tmp', status='unknown')
      Open(unit=14,file='Elem-4.tmp', status='unknown')
      Open(unit=15,file='Elem-5.tmp', status='unknown')
            
      do 901 j=1,elementi
         if(tt(j).eq.1) go to 111
         if(tt(j).eq.2) go to 112
	 if(tt(j).eq.3) go to 113
         if(tt(j).eq.4) go to 114
	 if(tt(j).eq.5) go to 115
	 

111      write(11,1000) 
         write(11,*) asci(1,j),ordi(1,j)
         write(11,*) asci(2,j),ordi(2,j)
         write(11,*) asci(3,j),ordi(3,j)
         go to 901
	 
112      write(12,1000) 
         write(12,*) asci(1,j),ordi(1,j)
         write(12,*) asci(2,j),ordi(2,j)
         write(12,*) asci(3,j),ordi(3,j) 
         go to 901
	 
113      write(13,1000) 
         write(13,*) asci(1,j),ordi(1,j)
         write(13,*) asci(2,j),ordi(2,j)
         write(13,*) asci(3,j),ordi(3,j) 
         go to 901
	 
114      write(14,1000) 
         write(14,*) asci(1,j),ordi(1,j)
         write(14,*) asci(2,j),ordi(2,j)
         write(14,*) asci(3,j),ordi(3,j) 
         go to 901
	 
115      write(15,1000) 
         write(15,*) asci(1,j),ordi(1,j)
         write(15,*) asci(2,j),ordi(2,j)
         write(15,*) asci(3,j),ordi(3,j) 
         go to 901
 	 
901     continue
       close(10)
       close(8)
      
c     Posiziona il numero di ogni elemento al centro dell'elemento

      open(unit=10,file='N_elem.tmp',status='unknown')
      close(unit=10,status='delete')
      open(unit=10,file='N_elem.tmp',status='new')
      rewind(10)
      
       do j=1,elementi
         asci_M(j)=(asci(1,j)+asci(2,j)+asci(3,j))/3
         ordi_M(j)=(ordi(1,j)+ordi(2,j)+ordi(3,j))/3
         write(10,*) j,asci_M(j),ordi_M(j)
       enddo        
       close(10)


c.....formats
1000   format('>')
1001   format(2x,4HELEM,1x,4HPROP,1x,18HCONNESSIONI NODALI)
1002  format(3x,4HNODO,5x,1HX,9x,1HY)
      end
