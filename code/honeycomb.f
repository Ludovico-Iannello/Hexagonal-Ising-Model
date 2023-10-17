c programma per la simulazione del modello di ising bidimensionale
      program ising
      parameter(nlatt=20,nvol=2*nlatt**2,nbeta=20,measures=20000)
      integer iflag,i_decorrel,term
      real magn,ene,chi,cv,b,dm,de,dchi,dcv,db
      real magn_arr,ene_arr,var_ene, var_magn
      real spin,beta,extfield
      real t1,t2
      common/arr/magn_arr(measures),ene_arr(measures)
      common/lattice/spin(nvol)     !! matrice con le variabili di spin
                                        !! che rappresenta lo stato del sistema
      common/various/ beta,extfield       !! parametri della simulazione

      call second(t1)

      call ranstart                 !! inizializza random number generator
                                    !! per chiamare un numero random usare
                                    !! la funzione ran2()

cc apertura file sul quale scrivere le misure della magnetizzazione
      open(2,file='data_L20_histo_bassa.dat',status='new')
      open(4,file='ising_L20_histo_bassa.dat',status='new')

cc def parametri della simulazione
      iflag=1
      extfield=0
      i_decorrel=2
      term=6000

cc simulazione al variare di beta
      do ibeta = 1,nbeta
        write(*,*) ibeta
        magn=0
        ene=0
        cv=0
        chi=0
        var_ene=0
        var_magn=0


        beta=0.5+((0.8-0.5)/nbeta)*ibeta
        if (abs(beta-0.65).lt.0.1) then
          i_decorrel=40
        else
          i_decorrel=20
        end if
        
        write(*,*) beta
        CALL geometry()                     !!inizializza condizioni al bordo
        CALL initialize_lattice(iflag)      !!inizializza configurazione iniziale


        do idec = 1,i_decorrel*term
          call update_metropolis()       !!per termalizzare
        enddo

        do iter = 1,measures

CC   AGGIORNAMENTO CONFIGURAZIONE: i_decorrel spazzate di tutto il reticolo

          do idec = 1,i_decorrel
             call update_metropolis()
          enddo

CC   MISURA DELLE OSSERVABILI FISICHE
          call magnetization(xmagn)
          call energy(xene)
CC   SCRIVO LE MISURE SU FILE PER POI EFFETTUARE L'ANALISI
          write(2,*) xmagn,xene

CC abs predno il val ass di magnetizzzione
          xmagn=abs(xmagn)



          magn_arr(iter)=xmagn          !!def array
          ene_arr(iter)=xene
       enddo !! measures

       magn=sum(magn_arr)/float(measures)
       ene=sum(ene_arr)/float(measures)               !!medie valori

       do i=1,measures
         var_magn=(magn_arr(i)-magn)**2+var_magn      !!calcolo varianze
         var_ene=(ene_arr(i)-ene)**2+var_ene
       enddo

CC     calcolo cv e suscettività dalle var
       cv=beta**2*var_ene*nvol/float(measures)
       chi=beta*var_magn*nvol/float(measures)

CC calcolo errori e binder
       call bootstrap(ene_arr,de,dcv)
       call bootstrap(magn_arr,dm,dchi)
       call binder(magn_arr,b,db)

       write(4,*) beta,magn,dm,ene,de,chi,dchi*beta,cv,dcv*beta**2,b,db
       call second(t2)
       write(*,*) 'il tempo impiegato è ', t2-t1
      enddo

CC=========TERMINE SIMULAZIONE MONTE-CARLO===========

CC SALVO CONFIGURAZIONE E STATO GEN. RANDOM PER
CC POTER EVENTUALMENTE RIPARTIRE
c      open(3,file='lattice',status='unknown')
c      write(3,*) spin
c      close(3)
      close(2)
      close(4)
      call ranfinish

      call second(t2)
      write(*,*) 'il tempo impiegato è ', (t2-t1)/60

      stop
      end program

CC      INIZIO SUBROUTINES

c*****************************************************************
      subroutine geometry()
c*****************************************************************
c per ogni coordinata definisco il passo in avanti o indietro
c con le opportune condizioni al bordo
c=================================================================
      parameter (nlatt=20,nvol = 2*nlatt**2)
      common/move/ npp(nvol+2*nlatt),nmm(nvol+2*nlatt)

      !! le funzioni npp ed nmm sono costruite come dei vettori
      !! di interi, in modo da non essere ricalcolate ogni volta
      !! e rendere tutto piu` efficiente, prendono in input una
      !! coordinata e restituiscono la coordinata in avanti o
      !! indietro, tenendo conto delle opportune condizioni

      do i = 1,nvol
         npp(i+2*nlatt-1) = i
         nmm(i+1) = i
      enddo

      do i=1, 2*nlatt-1
        npp(i)=nvol-2*nlatt+1+i
        nmm(i+nvol+1)=i
      enddo


      npp(nvol+2*nlatt) = 1             !! RIAGGIUSTO IL CALCOLO AI BORDI PER TENERE
      nmm(1) = nvol             !! CONTO DELLE CONDIZIONI AL BORDO PERIODICHE

      return
      end

c=================================================================
      subroutine initialize_lattice(iflag)
c ASSEGNO LA CONFIGURAZIONE DI PARTENZA DELLA CATENA DI MARKOV
      parameter (nlatt=20, nvol = 2*nlatt**2)
      common/lattice/ spin(nvol)

CC  PARTENZA A FREDDO (tutti gli spin a 1 come se fosse T = 0)
      if (iflag.eq.0) then
         do i = 1,nvol                  !! loop su tutti i siti
            spin(i) = 1.0             !! del reticolo
          enddo
CC  ... A CALDO ... (spin random, come se fosse T = infinito)
      elseif (iflag.eq.1) then
         do i = 1,nvol                  !! loop su tutti i siti del reticolo
           x = ran2()                    !! ran2() random fra 0 e 1
           spin(i) = 1.0
           if (x.lt.0.5) spin(i) = -1.0
         enddo
CC  ... O DA DOVE ERO RIMASTO L'ULTIMA VOLTA
      else
         open(9,file='lattice',status='old')
         read(9,*) spin
         close(9)
      endif

      return
      end
c=================================================================

c*****************************************************************
      subroutine magnetization(xmagn)
c calcolo magnetizzazione media del reticolo
c*****************************************************************
      parameter (nlatt=20,nvol = 2*nlatt**2)
      common/lattice/ spin(nvol)

      xmagn = 0.0                       !! inizializzo xmagn a zero
      do i = 1,nvol                    !! faccio il loop su tutto il reticolo
        xmagn = xmagn + spin(i)        !! e sommo tutti i valori del campo
       enddo

      xmagn = xmagn/float(nvol)         !! normalizzo dividendo per il volume

      return
      end
c=================================================================

c*****************************************************************
      subroutine energy(xene)
C energia media (= 0 per configurazione ordinata e campo esterno 0)
c=================================================================
      parameter (nlatt=20, nvol = 2*nlatt**2)
      common/lattice/ spin(nvol)
      common/move/ npp(nvol+2*nlatt),nmm(nvol+2*nlatt)
      common/various/ beta,extfield

      xene = 0.0                             !! inizializzo variabile
      do i = 1,nvol                        !! inizio il loop sul reticolo
        if (modulo(i,2)==0) then
          force = spin(npp(i+2*nlatt)) + spin(nmm(i)) +     !! somma dei 3 primi vicini
     $            spin(npp(i))
        else
          force = spin(npp(i+2*nlatt)) + spin(nmm(i)) +     !! somma dei 3 primi vicini
     $            spin(nmm(i+2*nlatt))
        endif
        xene = xene -  0.5*force*spin(i)  !! 1/2 per conteggio giusto
c        xene = xene - extfield*spin(i)       !! contributo campo esterno
      enddo

      xene = xene/float(nvol)            !! normalizzo -> densita` di energia

      return
      end
c=================================================================

      subroutine update_metropolis()
cc  faccio aggiornamenti locali delle variabili di spin con metropolis
cc  la variabile di spin di prova e` sempre quella opposta a quella attuale
      parameter (nlatt=20, nvol = 2*nlatt**2)
      common/lattice/ spin(nvol)
      common/move/ npp(nvol+2*nlatt),nmm(nvol+2*nlatt)
      common/various/ beta,extfield

      do ivol = 1,nvol            !! loop su tutti i siti

        i = int(ran2()*nvol + 1.0)  !! scelgo a caso un sito del reticolo

        if (modulo(i,2)==0) then
          force = spin(npp(i+2*nlatt)) + spin(nmm(i)) +     !! somma dei 3 primi vicini
     $            spin(npp(i))
        else
          force = spin(npp(i+2*nlatt)) + spin(nmm(i)) +     !! somma dei 3 primi vicini
     $            spin(nmm(i+2*nlatt))
        endif

        force = beta*(force + extfield)     !! campo esterno e divido per
                                                !! (kT): ottengo quello che
                                                !! moltiplicato per lo spin in
                                                !! (i,j) mi da la parte di
                                                !! energia che dipende solo
                                                !! dallo spin in (i,j), da cui
                                                !! deriva il rapporto di prob.

        phi =  spin(i)                  !! phi = valore attuale dello spin.
                                        !! Il tentativo e` sempre invertire lo
                                        !! lo spin. Calcolo il rapporto p_rat
        p_rat = exp(-2.0*phi*force)     !! di prob. fra il caso invertito e non

        x = ran2()                            !! METRO-TEST! x = random (0,1)
                                              !! x < p_rat verifica anche caso
        if (x.lt.p_rat) spin(i) = -phi        !! p_rat > 1 -> se si accetto


      enddo

      return
      end

CC---  bootstrap con correlazione ----------------------------------------------
c----  k numero simulazioni, p lunghezza della catena da tagliare---------------

      subroutine bootstrap (x,dx,dx1)
        integer j,q,i,N,l,k,p,it,nlatt,nvol
        parameter (nlatt=20,k=5000,N=20000,p=100,nvol=2*nlatt**2)
        real dx,dx1,var,media_x, media_dx
        real sim_arr(k),varsim_arr(k)       !!media e varianza simulazioni
        real a(N),x(N)        !!a array provvisorio

        dx=0
        dx1=0

        do j=1,k
          var=0
          do i=1,N,p
            l=int((N-p)*rand()+1)
            do q=1,p !punti successivi della catena che taglio
              a(i+q-1)=x(l+q)
            end do
          end do
          sim_arr(j)=sum(a)/float(N)
          do it=1,N
            var=(a(it)-sim_arr(j))**2+var      !!calcolo varianze
          enddo
          varsim_arr(j)=var*nvol/float(N)   !!chi/heat manca una potenza di beta
        end do

        media_x=sum(sim_arr)/float(k)
        media_dx=sum(varsim_arr)/float(k)

        do i = 1,k
          dx=(sim_arr(i)-media_x)**2+dx
          dx1=(varsim_arr(i)-media_dx)**2+dx1
        end do

        dx=sqrt(dx/((k-1)))
        dx1=sqrt(dx1/((k-1)))
        return
      end
c============================================================================
c  binder
c============================================================================
      subroutine binder (x,b,db)
       integer j,q,i,N,l,k,p,nlatt,nvol
       parameter (nlatt=20,k=5000,N=20000,p=100,nvol=2*nlatt**2)
       real db,media_b,b,y4,y2
       real sim_arr(k)       !!media e varianza simulazioni
       real x(N)        !!a array provvisorio

       db=0
       b=N*sum(x**4)/sum(x**2)**2

       do j=1,k
         y4=0
         y2=0
         do i=1,N,p
           l=int((N-p)*rand()+1)
           do q=1,p !punti successivi della catena che taglio
             y4=x(l+q)**4+y4
             y2=x(l+q)**2+y2
           end do
         end do
         sim_arr(j)=y4*N/(y2**2)
       end do
       media_b=sum(sim_arr)/float(k)

       do i = 1,k
         db=(sim_arr(i)-media_b)**2+db
       end do
       db=sqrt(db/((k-1)))
       return
      end

c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c=============================================================================
      subroutine ranstart
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1
      close(23)
 118  continue

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end
c=============================================================================
