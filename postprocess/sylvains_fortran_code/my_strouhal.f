      program strouhal
c
      parameter (ntg=450000)
      parameter (ncut=32768)
      parameter (iover=16384)
      parameter (lomb=int(ntg/ncut))
      parameter (nt=lomb*ncut)
      parameter (mx=ncut+2)
      parameter (mw=3*mx/2)
      parameter (mxf=2*ncut)
      parameter (nspec=2*lomb-1)
      
 
c
      character*3 car
      character*40 filename
c
      common /wrk/    work(mw)
      common /wrkx1/  worky1(mx)
      common /ftx/    ifaxy(20)
      common /exx/    trigsy(mxf)
      common /exx2/   trigsy2(mxf)

      real(4),dimension(ntg,303) :: tb4
      real(4),dimension(ntg) :: signal_G1
      real(4),dimension(nt) :: signal_G
      real(4),dimension(mx) :: signal
      real(4),dimension(ncut) :: hann
      real(4),dimension(ncut,nspec) :: spec_f,signal1    



      dt=0.010
c      ixsonde=60

      open(775,file='MOY2440FFT',form='unformatted',
     2         status='unknown')
      read(775) tb4
      close(775)

      print *,'read file done'



      do ixsonde=1,10!1,303

      print *,'ixsonde',ixsonde

      do i=1,ntg
         signal_G(i)=0.
         signal_G1(i)=0.
      enddo
      do i=1,ntg
         signal_G(i)=tb4(i,ixsonde)
         signal_G1(i)=tb4(i,ixsonde)
      enddo
c      open(775,file='signali1.dat',form='formatted',status='unknown')      
c      do i=1,nt
c      read(775,*) dy,signal_G(i)
c      enddo
c      close(775)

      open(775,file='signal_SF.dat',form='formatted',status='unknown')      
      do i=1,ntg
      write(775,*) (i-1)*dt,signal_G1(i)
      enddo
      close(775)

         umm=0.
         do i=1,nt
            umm=umm+signal_G(i)
         enddo
         umm=umm/float(nt)
         write(*,*) 'Average value for the signal :', umm
         do i=1,nt
            signal_G(i)=signal_G(i)-umm
         enddo



c     POUR TRACER LE SIGNAL TOTAL POUR VERIF
c      stop

c     mise en oeuvre du periodogramme avec fenetre de hanning
      do kk=1,nspec
         print *,'nspec',kk
c     decoupage de la grande sonde en 'LOMB' sondes de taille 'NCUT'    

         do i=1,ncut
            signal(i)=0.
            signal(i)=signal_G((kk-1)*ncut/2+i) !pour du 50%
         enddo
c         print *,kk,(kk-1)*ncut/2+1,(kk-1)*ncut/2+ncut
c ***********************************
c     fenetrage (hanning)
c ***********************************
c
         call shanning(hann,ncut)
c
         do i=1,ncut
            signal(i)=hann(i)*signal(i)
         enddo
c         write(*,*)'Fenetrage OK'
         
         call fftfax (ncut,ifaxy,trigsy)
         call fft999 (signal,work,trigsy,ifaxy,1,mx,ncut,1,-1)
c
         do i=1,ncut
            spec_f(i,kk)=signal(i)
         enddo
c         print *,'FFT done!'
      enddo


c     On aditionne le tout
      do kk=1,nspec
         do i=1,ncut
            signal(i)=0.
         enddo
         do i=1,ncut,2
            signal(i)=(spec_f(i,kk)*spec_f(i,kk)+
     1        spec_f(i+1,kk)*spec_f(i+1,kk))
         enddo
         do i=1,ncut/2
            signal1(i,kk)=0.
            signal1(i,kk)=signal(2*i-1)
         enddo
      enddo
      do kk=2,nspec
      do i=1,ncut/2
         signal1(i,1)=signal1(i,1)+signal1(i,kk)
      enddo
      enddo

      print *,'NSPEC',nspec

      open(20,file='test22.dat',
     1     form='formatted',status='unknown')
      ir=1
      do i=3,ncut/2
         freq=float(ir)/(dt*ncut)
         write(20,*) ir,freq,signal1(i,1)/nspec
         ir=ir+1
      enddo
      close(20)
  990 format('P_SP',I4.4)

      
       jj=(ixsonde-1)*1.2664/303*1000.
       print *,jj,jj,'JJJJJ'
      write(filename, 990) jj
!      call system('gnuplot courbe.gnu')
!      call system('gnuplot courbe_jpeg.gnu')
!      call system('mv test.jpeg '//filename(1:7)//'.jpeg')
!      call system('mv test.ps '//filename(1:7)//'.ps') 
      call system('mv test22.dat '//filename(1:8)//'_1.dat')

      call system('mv signal_SF.dat '//filename(1:8)//'_2.dat')
      enddo

      end
c
c ****************************************************
c
      subroutine shanning(hann,nxb)
c
      dimension hann(nxb)
c      
      dpisn=2.*acos(-1.)/float(nxb)
      coeff=sqrt(2./3.)
      do  i=1,nxb
         hann(i)=coeff*(1.-cos(dpisn*float(i)))
      enddo
      return
      end
c
c
c***********************************************************************
c
      subroutine fct99 (data,work,trig1,trig2,a2,ifax,
     1                  inc,jump,n,m,isign)
c
c***********************************************************************
c
c  (**  Does fast cosine transforms from real to cosine space.   **)
c
c       Based on the algorithm outlined in 
c       "The Fast Fourier Transform Algorithm
c  Programming Considerations in the Calculation of Sine, Cosine
c  and Laplace Transforms" by Cooley, Lewis and Welch
c  (Sound Vib. (1970) 12, 315-337)
c
c  The idea is to reduce the cosine coefficients into a form
c  which can be used by FFT991.  This version is vectorized.  **)
c
c  (**  INPUT :
c       DATA : Field to be transformed
c       WORK : Work array to rearrange DATA for FFT991 and
c              to vectorize.  Length = N+2 + M if not vectorized,
c              = (N+3)*M if vectorized.
c       TRIG1 : Coefficients for FFT991.
c       TRIG2 : Coefficients used after FFT991 to rearrange into
c               final result.
c       IFAX : Integer coefficients for Temperton version of
c              FFT991.
c       INC : Increment between elements of a single transform.
c       JUMP : Jump between elements of different transforms.
c       N : Length of transform.
c       M : Number of transforms to be done at once.
c       ISIGN : -1 : Real to cosines.
c                1 : Cosines to real.     **)
c
      dimension trig1(1),trig2(1),ifax(1),work(1),data(1)
      dimension a2(1)
c
      jump1 = n + 2
c
      if (isign.eq.-1) then
      fac1 = 1./n
      fac2=.25     
      else
      fac1=1.
      fac2=.5
      endif
c
c
c  (**  Rearrange data for FFT999  **)
      do 30 j=1,m
      jj = (j-1)*jump+1
      work(1) = data(jj)*fac1
      work(n+1) = data(jj+n*inc)*fac1
c      work(2) = 0.
      work(n+2) = 0.
c
      sum = 0.0
      do 2 i=2,n,2
         sum=sum + data((i-1)*inc+jj)
 2    continue
      a2(j)= 2.*fac1*sum
c
      do 20 i = 3,n,2
      work(i) = fac1*data(jj+(i-1)*inc)
      work(i+1) = fac1*(data(jj+(i-2)*inc)-data(jj+i*inc))
   20 continue
      do 25 i = 1,jump1
   25 data(jj+(i-1)*inc) = work(i)
c
   30 continue
c
      call FFT999 (data,work,trig1,ifax,inc,jump,n,m,1)
c
c  (**  Rearrange output into final format  **)
c
      do 40 j = 1,m
         jj= (j-1)*jump+1
         do 35 i = 1,jump1
   35       work(i) = data(jj+(i-1)*inc)
         do 36 i=2,n
            ii=(i-1)*inc+jj
            data(ii)=fac2*((work(i)+work(n+2-i)) -
     1                     .5*(work(i)-work(n+2-i))/trig2(i))
   36    continue
      data(jj)= 2.*fac2*(work(1)+a2(j))
      data(n*inc+jj)= 2.*fac2*(work(1)-a2(j))
      if(data(jj).ne.0.0.or.data(jj+n*inc).ne.0.0) then
c     print *,data(jj),data(n*inc+jj),j
      endif
   40 continue
      return
      end
c
c***********************************************************************
c
      subroutine fst99 (data,work,trig1,trig2,ifax,inc,jump,n,m,isign)
c
c***********************************************************************
c
c  (**  Does fast sine transforms from real to sine space.   **)
c
c       Based on the algorithm outlined in :
c       "The Fast Fourier Transform Algorithm:
c  Programming Considerations in the Calculation of Sine, Cosine
c  and Laplace Transforms" by Cooley, Lewis and Welch.
c  (Sound Vib. (1970) 12, 315-337).
c
c  The idea is to reduce the sine coefficients into a form
c  which can be used by FFT991.  This version is vectorized.  **)
c
c  (**  INPUT :
c       DATA : Field to be transformed.
c       WORK : Work array to rearrange DATA for FFT991 and
c              to vectorize.  Length = N+2 + M if not vectorized,
c              = (N+3)*M if vectorized
c       TRIG1 : Coefficients for FFT991.
c       TRIG2 : Coefficients used after FFT991 to rearrange into
c               final result.
c       IFAX : Integer coefficients for Temperton version of
c              FFT991.
c       INC : Increment between elements of a single transform.
c       JUMP : Jump between elements of different transforms.
c       N : Length of transform.
c       M : Number of transforms to be done at once.
c       ISIGN : -1 : Real to sines.
c                1 : Sines to real.     **)
c
      dimension trig1(1),trig2(1),ifax(1),work(1),data(1)
c
      jump1 = n + 2
c
      if (isign.eq.-1) then
      fac1 = 1./n
      fac2=.25     
      else
      fac1=1.
      fac2=.5
      endif
c
c  (**  Rearrange data for FFT999  **)
c
      do 30 j=1,m
        jj = (j-1)*jump+1
        work(1) = -data(jj+inc)*fac1*2.
        work(n+1) = data(jj+(n-1)*inc)*fac1*2.
c        work(2) = 0.
        work(n+2) = 0.
c
        do 20 i = 3,n,2
          work(i) = fac1*(data(jj+(i-2)*inc)-data(jj+i*inc))
          work(i+1) = -fac1*data(jj+(i-1)*inc)
   20   continue
        do 25 i = 1,jump1
   25     data(jj+(i-1)*inc) = work(i)
   30 continue
c
      call FFT999 (data,work,trig1,ifax,inc,jump,n,m,1)
c
c  (**  Rearrange output into final format  **)
      do 40 j = 1,m
         jj= (j-1)*jump+1
         do 35 i = 1,jump1
   35       work(i) = data(jj+(i-1)*inc)
c
         data(jj) = 0.
         data(jj+n*inc) = 0.
c
        do 40 i=2,n
           ii=(i-1)*inc+jj
           data(ii)=fac2*((work(i)-work(n+2-i)) -
     1               .5*(work(i)+work(n+2-i))/trig2(i))
   40   continue
c
      return
      end
c
c***********************************************************************
c
      subroutine fft99 (a,work,trigs,ifax,inc,jump,n,lot,isign)
c
c***********************************************************************
* c06-summation of series                                       b6.1/3 *
*                                                                      *
*                                                               fft99  *
*                                                               fft991 *
*                                                                      *
*                                                                      *
* subprogram       subroutine    fft99                                 *
*                                fft991                                *
*                                                                      *
* purpose          perform multiple fast fourier transforms            *
*                                                                      *
*                                                                      *
* version          cyber                         cray-1                *
*                                                                      *
*                  jan 1979 original             jan 1979 original     *
*                                                                      *
* usage                                                                *
*                  call fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)   *
*                  call fft991(a,work,trigs,ifax,inc,jump,n,m,isign)   *
*                                                                      *
* arguments        1.dimension                                         *
*                       a(idim),work((n+1)*m),trigs(3*n/2),ifax(10)    *
*                       work is a work array                           *
*                                                                      *
*                  2.input                                             *
*                      a - an array containing the input data or       *
*                          coefficient vectors.                        *
*                         this array is overwritten by the results.    *
*                      trigs and ifax - arrays set up by fftrig and fax*
*                                     - see writeup of fftrig and fax  *
*                      inc - the word increment between successive     *
*                           elements of each data or coefficient vector*
*                           e.g. inc=1 for consecutively stored data.  *
*                      jump - the word increment between the first     *
*                            elements of successive data or coefficient*
*                            vectors.                                  *
*                      n - the length of each transform. (see note x)  *
*                      m - the number of transforms to be done         *
*                          simultaneously.                             *
*                      isign - +1 for a transform from fourier         *
*                              coefficients to data values.            *
*                              -1 for a transform from data values     *
*                              to fourier coefficients.                *
*                                                                      *
*                  3.output                                            *
*                      a - contains either the coefficients or the     *
*                          data values,depending on isign.             *
*                          in each case n independent quantities       *
*                          occupy n+2 words.   the coefficients are    *
*                          stored as successive pairs of real and      *
*                          imaginary parts -                           *
*                          a(k),b(k) , k=0,1,...n/2                    *
*                          b(0) and b(n/2) are stored although they    *
*                          must be 0.                                  *
*                      for fft99 the data is stored with explicit      *
*                          periodicity -                               *
*                          x(n-1),x(0),x(1),....x(n-1),x(0)            *
*                      for fft991 the data appears as -                *
*                          x(0),x(1),x(2),......x(n-1),0,0             *
*                                                                      *
* notes            1. on cray-1, arrange data so that jump is not a    *
*                     multiple of 8 (to avoid memory bank conflicts)   *
*                                                                      *
* write up         computer bulletin b6.6/1                            *
*                                                                      *
* entry points        fft99,fft991                                     *
*                                                                      *
* common blocks    none                                                *
*                                                                      *
* i/o              none                                                *
*                                                                      *
* precision        single                                              *
*                                                                      *
* other routines   fft99a,fft99b,vpassm          (cy)                  *
*       required   cal99,cpass                   (cr)                  *
*                                                                      *
*                                                                      *
* 7/80                      fft99-1                                    *
*                                                                      *
************************************************************************
*                                                                      *
* c06-summation of series                                       b6.1/3 *
*                                                                      *
*                                                               fft99  *
*                                                               fft991 *
*                                                                      *
* access (object)  cyber:                                              *
*                           attach,eclib.                              *
*                           ldset(lib=eclib)                           *
*                  cray 1:                                             *
*                           ldr(lib=eclib...)                          *
*                                                                      *
* access (source)           attach,oldpl,eclibpl                       *
*                                                                      *
*                  cyber :         %define cyber                       *
*                  cray:           %define cray                        *
*                                  %c    fft99,fft991                  *
*                                                                      *
* language         fortran                                             *
*                  but cray implementation of pass is in cal           *
*                                                                      *
* specialist       clive temperton                                     *
*                                                                      *
* history          written by c.temperton      jan     1979            *
*                                                                      *
* algorithm        the algorithm is the self-sorting (temperton)       *
*                  version of the fast fourier transform               *
*                                                                      *
* references       ecmwf technical report no.3                         *
*                  ecmwf internal report no.21 -   c.temperton         *
*                                                                      *
* object size               fft991  fft99  (octal words)               *
*                  cyber:    2665    2676                              *
*                  cray :    1250    1260                              *
*                                                                      *
*                                                                      *
* accuracy                                                             *
*                                                                      *
* timing           vectorization is on vectors of length m.      (cr)  *
*                  hence timing is strongly dependent on m.            *
*                  time per transform on cray-1 (microseconds)         *
*                  n    m=4    m=16    m=64                            *
*                 64     46      17      10                            *
*                128     81      33      21                            *
*                180    150      58      37                            *
*                192    149      58      36                            *
*                240    192      76      49                            *
*                256    191      76      49                            *
*                288    219      89      58                            *
*                300    253     102      68                            *
*                320    248     101      66                            *
*                360    286     118      79                            *
*               1024    898     359     238                            *
*                                                                      *
* portability      standard fortran                                    *
*                  standard cal  (cr)                                  *
*                                                                      *
* system routines  none                                                *
*        required                                                      *
*                                                                      *
* 7/80                      fft99-1                                    *
*                                                                      *
************************************************************************
c
c     subroutine 'fft99' - multiple fast real periodic transform
c     corresponding to old scalar routine fft9
c     procedure used to convert to half-length complex transform
c     is given by cooley, lewis ' welch (j. sound vib., vol. 12
c     (1970), 315-337)
c
c     a is the array containing input ' output data
c     work is an area of size (n+1)*lot
c     trigs is a previously prepared list of trig function values
c     ifax is a previously prepared list of factors of n/2
c     inc is the increment within each data "vector"
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isign = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c     ordering of coefficients:
c         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
c         where b(0)=b(n/2)=0; (n+2) locations required
c
c     ordering of data:
c         x(n-1),x(0),x(1),x(2),...,x(n),x(0)
c         i.e. explicit cyclic continuity; (n+2) locations required
c
c     vectorization is achieved on cray by doing the transforms in
c     parallel
c
c     *** n.b. n is assumed to be an even number
c
c     definition of transforms:
c     -------------------------
c
c     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c
c     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
c               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
c
c
c     line following next is not subroutine header(only comment)
c
c     subroutine fft99(a,work,trigs,ifax,inc,jump,n,lot,isign)
c
c     end
      dimension a(n),work(n),trigs(n),ifax(1)
c
      nfax=ifax(1)
      if(nfax.le.0) go to 99
      nx=n+1
      nh=n/2

      ink=inc+inc
      if (isign.eq.+1) go to 30
c
c     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=inc+1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
*VDIR NODEP
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
c
      igo=60
      go to 40
c
c     preprocessing (isign=+1)
c     ------------------------
c
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     complex transform
c     -----------------
c
   40 continue
      ia=inc+1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c
      if (isign.eq.-1) go to 130
c
c     if necessary, transfer data from work area
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=ia
      do 100 l=1,lot
      i=ibase
      j=jbase
*VDIR NODEP
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     fill in cyclic boundary points
  110 continue
      ia=1
      ib=n*inc+1
*VDIR NODEP
      do 120 l=1,lot
      a(ia)=a(ib)
      a(ib+inc)=a(ia+inc)
      ia=ia+jump
      ib=ib+jump
  120 continue
      go to 140
c
c     postprocessing (isign=-1):
c     --------------------------
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
      return
c  ** error exit   ifax(1) le 0 **
99    print *,'fft99 called but factors not supplied'
c     call abort
      stop 'sur ABORT'
      end
      subroutine fft99a(a,work,trigs,inc,jump,n,lot)
c
c     subroutine fft99a - preprocessing step for fft99, isign=+1
c     (spectral to gridpoint transform)
c
c     line following next is not subroutine header(only comment)
c
c     subroutine fft99a(a,work,trigs,inc,jump,n,lot)
c     end
      dimension a(n),work(n),trigs(n)
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) ' a(n/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
*VDIR NODEP
      do 10 l=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   10 continue
c
c     remaining wavenumbers
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1
c
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
*VDIR NODEP
      do 20 l=1,lot
      work(ja)=(a(ia)+a(ib))-
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(jb)=(a(ia)+a(ib))+
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+
     *    (a(ia+inc)-a(ib+inc))
      work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))-
     *    (a(ia+inc)-a(ib+inc))
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   20 continue
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
*VDIR NODEP
      do 40 l=1,lot
      work(ja)=2.0*a(ia)
      work(ja+1)=-2.0*a(ia+inc)
      ia=ia+jump
      ja=ja+nx
   40 continue
c
   50 continue
      return
      end
      subroutine fft99b(work,a,trigs,inc,jump,n,lot)
c
c     subroutine fft99b - postprocessing step for fft99, isign=-1
c     (gridpoint to spectral transform)
c
c
c     line follwing next is not subroutine header(only comment)
c
c     subroutine fft99b(work,a,trigs,inc,jump,n,lot)
c     end
      dimension work(n),a(n),trigs(n)
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) ' a(n/2)
      scale=1.0/float(n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
*VDIR NODEP
      do 10 l=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc)=0.0
      a(jb+inc)=0.0
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   10 continue
c
c     remaining wavenumbers
      scale=0.5*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1
c
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
*VDIR NODEP
      do 20 l=1,lot
      a(ja)=scale*((work(ia)+work(ib))
     *   +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(jb)=scale*((work(ia)+work(ib))
     *   -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(ja+inc)=scale*((c*(work(ia)-work(ib))
     *    -s*(work(ia+1)+work(ib+1)))
     *    +(work(ib+1)-work(ia+1)))
      a(jb+inc)=scale*((c*(work(ia)-work(ib))
     *    -s*(work(ia+1)+work(ib+1)))
     *    -(work(ib+1)-work(ia+1)))
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   20 continue
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
      scale=2.0*scale
*VDIR NODEP
      do 40 l=1,lot
      a(ja)=scale*work(ia)
      a(ja+inc)=-scale*work(ia+1)
      ia=ia+nx
      ja=ja+jump
   40 continue
c
   50 continue
      return
      end
c
c***********************************************************************
c
      subroutine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
c
c***********************************************************************
c
c     subroutine 'fft991' - multiple real/half-complex periodic
c     fast fourier transform
c
************************************************************************
*                                                                      *
* c06-summation of series                                       b6.1/3 *
*                                                                      *
*                                                               fft99  *
*                                                               fft991 *
*                                                                      *
*                                                                      *
* subprogram       subroutine    fft99                                 *
*                                fft991                                *
*                                                                      *
* purpose          perform multiple fast fourier transforms            *
*                                                                      *
*                                                                      *
* version          cyber                         cray-1                *
*                                                                      *
*                  jan 1979 original             jan 1979 original     *
*                                                                      *
* usage                                                                *
*                  call fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)   *
*                  call fft991(a,work,trigs,ifax,inc,jump,n,m,isign)   *
*                                                                      *
* arguments        1.dimension                                         *
*                       a(idim),work((n+1)*m),trigs(3*n/2),ifax(10)    *
*                       work is a work array                           *
*                                                                      *
*                  2.input                                             *
*                      a - an array containing the input data or       *
*                          coefficient vectors.                        *
*                         this array is overwritten by the results.    *
*                      trigs and ifax - arrays set up by fftrig and fax*
*                                     - see writeup of fftrig and fax  *
*                      inc - the word increment between successive     *
*                           elements of each data or coefficient vector*
*                           e.g. inc=1 for consecutively stored data.  *
*                      jump - the word increment between the first     *
*                            elements of successive data or coefficient*
*                            vectors.                                  *
*                      n - the length of each transform. (see note x)  *
*                      m - the number of transforms to be done         *
*                          simultaneously.                             *
*                      isign - +1 for a transform from fourier         *
*                              coefficients to data values.            *
*                              -1 for a transform from data values     *
*                              to fourier coefficients.                *
*                                                                      *
*                  3.output                                            *
*                      a - contains either the coefficients or the     *
*                          data values,depending on isign.             *
*                          in each case n independent quantities       *
*                          occupy n+2 words.   the coefficients are    *
*                          stored as successive pairs of real and      *
*                          imaginary parts -                           *
*                          a(k),b(k) , k=0,1,...n/2                    *
*                          b(0) and b(n/2) are stored although they    *
*                          must be 0.                                  *
*                      for fft99 the data is stored with explicit      *
*                          periodicity -                               *
*                          x(n-1),x(0),x(1),....x(n-1),x(0)            *
*                      for fft991 the data appears as -                *
*                          x(0),x(1),x(2),......x(n-1),0,0             *
*                                                                      *
* notes            1. on cray-1, arrange data so that jump is not a    *
*                     multiple of 8 (to avoid memory bank conflicts)   *
*                                                                      *
* write up         computer bulletin b6.6/1                            *
*                                                                      *
* entry points        fft99,fft991                                     *
*                                                                      *
* common blocks    none                                                *
*                                                                      *
* i/o              none                                                *
*                                                                      *
* precision        single                                              *
*                                                                      *
* other routines   fft99a,fft99b,vpassm          (cy)                  *
*       required   cal99,cpass                   (cr)                  *
*                                                                      *
*                                                                      *
* 7/80                      fft99-1                                    *
*                                                                      *
************************************************************************
*                                                                      *
* c06-summation of series                                       b6.1/3 *
*                                                                      *
*                                                               fft99  *
*                                                               fft991 *
*                                                                      *
* access (object)  cyber:                                              *
*                           attach,eclib.                              *
*                           ldset(lib=eclib)                           *
*                  cray 1:                                             *
*                           ldr(lib=eclib...)                          *
*                                                                      *
* access (source)           attach,oldpl,eclibpl                       *
*                                                                      *
*                  cyber :         %define cyber                       *
*                  cray:           %define cray                        *
*                                  %c    fft99,fft991                  *
*                                                                      *
* language         fortran                                             *
*                  but cray implementation of pass is in cal           *
*                                                                      *
* specialist       clive temperton                                     *
*                                                                      *
* history          written by c.temperton      jan     1979            *
*                                                                      *
* algorithm        the algorithm is the self-sorting (temperton)       *
*                  version of the fast fourier transform               *
*                                                                      *
* references       ecmwf technical report no.3                         *
*                  ecmwf internal report no.21 -   c.temperton         *
*                                                                      *
* object size               fft991  fft99  (octal words)               *
*                  cyber:    2665    2676                              *
*                  cray :    1250    1260                              *
*                                                                      *
*                                                                      *
* accuracy                                                             *
*                                                                      *
* timing           vectorization is on vectors of length m.      (cr)  *
*                  hence timing is strongly dependent on m.            *
*                  time per transform on cray-1 (microseconds)         *
*                  n    m=4    m=16    m=64                            *
*                 64     46      17      10                            *
*                128     81      33      21                            *
*                180    150      58      37                            *
*                192    149      58      36                            *
*                240    192      76      49                            *
*                256    191      76      49                            *
*                288    219      89      58                            *
*                300    253     102      68                            *
*                320    248     101      66                            *
*                360    286     118      79                            *
*               1024    898     359     238                            *
*                                                                      *
* portability      standard fortran                                    *
*                  standard cal  (cr)                                  *
*                                                                      *
* system routines  none                                                *
*        required                                                      *
*                                                                      *
* 7/80                      fft99-1                                    *
*                                                                      *
************************************************************************
c
c     same as fft99 except that ordering of data corresponds to
c     that in mrfft2
c
c     procedure used to convert to half-length complex transform
c     is given by cooley, lewis ' welch (j. sound vib., vol. 12
c     (1970), 315-337)
c
c     a is the array containing input ' output data
c     work is an area of size (n+1)*lot
c     trigs is a previously prepared list of trig function values
c     ifax is a previously prepared list of factors of n/2
c     inc is the increment within each data "vector"
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isign = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c     ordering of coefficients:
c         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
c         where b(0)=b(n/2)=0; (n+2) locations required
c
c     ordering of data:
c         x(0),x(1),x(2),...,x(n-1)
c
c     vectorization is achieved on cray by doing the transforms in
c     parallel
c
c     *** n.b. n is assumed to be an even number
c
c     definition of transforms:
c     -------------------------
c
c     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c
c     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
c               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
c
c     subroutine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
c     end
      dimension a(n),work(n),trigs(n),ifax(1)
c
      nfax=ifax(1)
      if(nfax.le.0) go to 99
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
c
c     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
*VDIR NODEP
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
c
      igo=60
      go to 40
c
c     preprocessing (isign=+1)
c     ------------------------
c
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     complex transform
c     -----------------
c
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c
      if (isign.eq.-1) go to 130
c
c     if necessary, transfer data from work area
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
*VDIR NODEP
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     fill in zeros at end
  110 continue
      ib=n*inc+1
*VDIR NODEP
      do 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 continue
      go to 140
c
c     postprocessing (isign=-1):
c     --------------------------
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
      return
c   **  error     ifax(1) le 0  **
99    print *,' fft991 called but factors not supplied '
c      print *,'ifax = ',ifax(1)
c     call abort
      stop 'sur ABORT'
      end
      subroutine vpassm(a,b,c,d,trigs,
     *  inc1,inc2,inc3,inc4,lot,n,ifac,la)
c
c     subroutine 'vpassm' - multiple version of 'vpassa'
c     performs one pass through data
c     as part of multiple complex fft routine
c     a is first real input vector
c     b is first imaginary input vector
c     c is first real output vector
c     d is first imaginary output vector
c     trigs is precalculated table of sines ' cosines
c     inc1 is addressing increment for a and b
c     inc2 is addressing increment for c and d
c     inc3 is addressing increment between a's & b's
c     inc4 is addressing increment between c's & d's
c     lot is the number of vectors
c     n is length of vectors
c     ifac is current factor of n
c     la is product of previous factors
c
c     subroutine vpassm(a,b,c,d,trigs,
c     *   inc1,inc2,inc3,inc4,lot,n,ifac,la)
c
c     end
      dimension a(n),b(n),c(n),d(n),trigs(n)
      data sin36/0.587785252292473/,cos36/0.809016994374947/,
     *     sin72/0.951056516295154/,cos72/0.309016994374947/,
     *     sin60/0.866025403784437/
c
c      print *,'Dans vpassm, lot=',lot
      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
c     check factors are correct - ensure non-negative
      if (igo.le.0) goto 998
      if (igo.gt.4) go to 999
      go to (10,50,90,130),igo
c
c     coding for factor 2
c
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 l=1,la
      i=ibase
      j=jbase
*VDIR NODEP
      do 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 l=1,la
      i=ibase
      j=jbase
*VDIR NODEP
      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
      return
c
c     coding for factor 3
c
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 l=1,la
      i=ibase
      j=jbase
*VDIR NODEP
      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 l=1,la
      i=ibase
      j=jbase
*VDIR NODEP
      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=
     *   c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *  -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=
     *   s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *  +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=
     *   c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *  -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=
     *   s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *  +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
      return
c
c     coding for factor 4
c
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 l=1,la
      i=ibase
      j=jbase
*VDIR NODEP
      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 l=1,la
      i=ibase
      j=jbase
*VDIR NODEP
      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=
     *    c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=
     *    s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=
     *    c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=
     *    s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=
     *    c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=
     *    s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
      return
c
c     coding for factor 5
c
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 l=1,la
      i=ibase
      j=jbase
*VDIR NODEP
      DO 135 IJK = 1, LOT
        C(1+J+(IJK-1)*INC4) = A(1+I+(IJK-1)*INC3) + (A(IB+I+(IJK-1)*
     1      INC3)+A(IE+I+(IJK-1)*INC3)) + (A(IC+I+(IJK-1)*INC3)+A(ID+I+(
     2      IJK-1)*INC3))
        D(1+J+(IJK-1)*INC4) = B(1+I+(IJK-1)*INC3) + (B(IB+I+(IJK-1)*
     1      INC3)+B(IE+I+(IJK-1)*INC3)) + (B(IC+I+(IJK-1)*INC3)+B(ID+I+(
     2      IJK-1)*INC3))
        C(JB+J+(IJK-1)*INC4) = (A(1+I+(IJK-1)*INC3)+COS72*(A(IB+I+(IJK-
     1      1)*INC3)+A(IE+I+(IJK-1)*INC3))-COS36*(A(IC+I+(IJK-1)*INC3)+A
     2      (ID+I+(IJK-1)*INC3))) - (SIN72*(B(IB+I+(IJK-1)*INC3)-B(IE+I+
     3      (IJK-1)*INC3))+SIN36*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-1)*
     4      INC3)))
        C(JE+J+(IJK-1)*INC4) = (A(1+I+(IJK-1)*INC3)+COS72*(A(IB+I+(IJK-
     1      1)*INC3)+A(IE+I+(IJK-1)*INC3))-COS36*(A(IC+I+(IJK-1)*INC3)+A
     2      (ID+I+(IJK-1)*INC3))) + (SIN72*(B(IB+I+(IJK-1)*INC3)-B(IE+I+
     3      (IJK-1)*INC3))+SIN36*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-1)*
     4      INC3)))
        D(JB+J+(IJK-1)*INC4) = (B(1+I+(IJK-1)*INC3)+COS72*(B(IB+I+(IJK-
     1      1)*INC3)+B(IE+I+(IJK-1)*INC3))-COS36*(B(IC+I+(IJK-1)*INC3)+B
     2      (ID+I+(IJK-1)*INC3))) + (SIN72*(A(IB+I+(IJK-1)*INC3)-A(IE+I+
     3      (IJK-1)*INC3))+SIN36*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     4      INC3)))
        D(JE+J+(IJK-1)*INC4) = (B(1+I+(IJK-1)*INC3)+COS72*(B(IB+I+(IJK-
     1      1)*INC3)+B(IE+I+(IJK-1)*INC3))-COS36*(B(IC+I+(IJK-1)*INC3)+B
     2      (ID+I+(IJK-1)*INC3))) - (SIN72*(A(IB+I+(IJK-1)*INC3)-A(IE+I+
     3      (IJK-1)*INC3))+SIN36*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     4      INC3)))
        C(JC+J+(IJK-1)*INC4) = (A(1+I+(IJK-1)*INC3)-COS36*(A(IB+I+(IJK-
     1      1)*INC3)+A(IE+I+(IJK-1)*INC3))+COS72*(A(IC+I+(IJK-1)*INC3)+A
     2      (ID+I+(IJK-1)*INC3))) - (SIN36*(B(IB+I+(IJK-1)*INC3)-B(IE+I+
     3      (IJK-1)*INC3))-SIN72*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-1)*
     4      INC3)))
        C(JD+J+(IJK-1)*INC4) = (A(1+I+(IJK-1)*INC3)-COS36*(A(IB+I+(IJK-
     1      1)*INC3)+A(IE+I+(IJK-1)*INC3))+COS72*(A(IC+I+(IJK-1)*INC3)+A
     2      (ID+I+(IJK-1)*INC3))) + (SIN36*(B(IB+I+(IJK-1)*INC3)-B(IE+I+
     3      (IJK-1)*INC3))-SIN72*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-1)*
     4      INC3)))
        D(JC+J+(IJK-1)*INC4) = (B(1+I+(IJK-1)*INC3)-COS36*(B(IB+I+(IJK-
     1      1)*INC3)+B(IE+I+(IJK-1)*INC3))+COS72*(B(IC+I+(IJK-1)*INC3)+B
     2      (ID+I+(IJK-1)*INC3))) + (SIN36*(A(IB+I+(IJK-1)*INC3)-A(IE+I+
     3      (IJK-1)*INC3))-SIN72*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     4      INC3)))
        D(JD+J+(IJK-1)*INC4) = (B(1+I+(IJK-1)*INC3)-COS36*(B(IB+I+(IJK-
     1      1)*INC3)+B(IE+I+(IJK-1)*INC3))+COS72*(B(IC+I+(IJK-1)*INC3)+B
     2      (ID+I+(IJK-1)*INC3))) - (SIN36*(A(IB+I+(IJK-1)*INC3)-A(IE+I+
     3      (IJK-1)*INC3))-SIN72*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     4      INC3)))
  135 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 l=1,la
      i=ibase
      j=jbase
*VDIR NODEP
      DO 145 IJK = 1, LOT
         C(1+J+(IJK-1)*INC4) = A(1+I+(IJK-1)*INC3) + (A(IB+I+(IJK-1)*
     1      INC3)+A(IE+I+(IJK-1)*INC3)) + (A(IC+I+(IJK-1)*INC3)+A(ID+I+(
     2      IJK-1)*INC3))
         D(1+J+(IJK-1)*INC4) = B(1+I+(IJK-1)*INC3) + (B(IB+I+(IJK-1)*
     1      INC3)+B(IE+I+(IJK-1)*INC3)) + (B(IC+I+(IJK-1)*INC3)+B(ID+I+(
     2      IJK-1)*INC3))
         C(JB+J+(IJK-1)*INC4) = C1*((A(1+I+(IJK-1)*INC3)+COS72*(A(IB+I+(
     1      IJK-1)*INC3)+A(IE+I+(IJK-1)*INC3))-COS36*(A(IC+I+(IJK-1)*
     2      INC3)+A(ID+I+(IJK-1)*INC3)))-(SIN72*(B(IB+I+(IJK-1)*INC3)-B(
     3      IE+I+(IJK-1)*INC3))+SIN36*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-
     4      1)*INC3)))) - S1*((B(1+I+(IJK-1)*INC3)+COS72*(B(IB+I+(IJK-1)
     5      *INC3)+B(IE+I+(IJK-1)*INC3))-COS36*(B(IC+I+(IJK-1)*INC3)+B(
     6      ID+I+(IJK-1)*INC3)))+(SIN72*(A(IB+I+(IJK-1)*INC3)-A(IE+I+(
     7      IJK-1)*INC3))+SIN36*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     8      INC3))))
         D(JB+J+(IJK-1)*INC4) = S1*((A(1+I+(IJK-1)*INC3)+COS72*(A(IB+I+(
     1      IJK-1)*INC3)+A(IE+I+(IJK-1)*INC3))-COS36*(A(IC+I+(IJK-1)*
     2      INC3)+A(ID+I+(IJK-1)*INC3)))-(SIN72*(B(IB+I+(IJK-1)*INC3)-B(
     3      IE+I+(IJK-1)*INC3))+SIN36*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-
     4      1)*INC3)))) + C1*((B(1+I+(IJK-1)*INC3)+COS72*(B(IB+I+(IJK-1)
     5      *INC3)+B(IE+I+(IJK-1)*INC3))-COS36*(B(IC+I+(IJK-1)*INC3)+B(
     6      ID+I+(IJK-1)*INC3)))+(SIN72*(A(IB+I+(IJK-1)*INC3)-A(IE+I+(
     7      IJK-1)*INC3))+SIN36*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     8      INC3))))
         C(JE+J+(IJK-1)*INC4) = C4*((A(1+I+(IJK-1)*INC3)+COS72*(A(IB+I+(
     1      IJK-1)*INC3)+A(IE+I+(IJK-1)*INC3))-COS36*(A(IC+I+(IJK-1)*
     2      INC3)+A(ID+I+(IJK-1)*INC3)))+(SIN72*(B(IB+I+(IJK-1)*INC3)-B(
     3      IE+I+(IJK-1)*INC3))+SIN36*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-
     4      1)*INC3)))) - S4*((B(1+I+(IJK-1)*INC3)+COS72*(B(IB+I+(IJK-1)
     5      *INC3)+B(IE+I+(IJK-1)*INC3))-COS36*(B(IC+I+(IJK-1)*INC3)+B(
     6      ID+I+(IJK-1)*INC3)))-(SIN72*(A(IB+I+(IJK-1)*INC3)-A(IE+I+(
     7      IJK-1)*INC3))+SIN36*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     8      INC3))))
         D(JE+J+(IJK-1)*INC4) = S4*((A(1+I+(IJK-1)*INC3)+COS72*(A(IB+I+(
     1      IJK-1)*INC3)+A(IE+I+(IJK-1)*INC3))-COS36*(A(IC+I+(IJK-1)*
     2      INC3)+A(ID+I+(IJK-1)*INC3)))+(SIN72*(B(IB+I+(IJK-1)*INC3)-B(
     3      IE+I+(IJK-1)*INC3))+SIN36*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-
     4      1)*INC3)))) + C4*((B(1+I+(IJK-1)*INC3)+COS72*(B(IB+I+(IJK-1)
     5      *INC3)+B(IE+I+(IJK-1)*INC3))-COS36*(B(IC+I+(IJK-1)*INC3)+B(
     6      ID+I+(IJK-1)*INC3)))-(SIN72*(A(IB+I+(IJK-1)*INC3)-A(IE+I+(
     7      IJK-1)*INC3))+SIN36*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     8      INC3))))
         C(JC+J+(IJK-1)*INC4) = C2*((A(1+I+(IJK-1)*INC3)-COS36*(A(IB+I+(
     1      IJK-1)*INC3)+A(IE+I+(IJK-1)*INC3))+COS72*(A(IC+I+(IJK-1)*
     2      INC3)+A(ID+I+(IJK-1)*INC3)))-(SIN36*(B(IB+I+(IJK-1)*INC3)-B(
     3      IE+I+(IJK-1)*INC3))-SIN72*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-
     4      1)*INC3)))) - S2*((B(1+I+(IJK-1)*INC3)-COS36*(B(IB+I+(IJK-1)
     5      *INC3)+B(IE+I+(IJK-1)*INC3))+COS72*(B(IC+I+(IJK-1)*INC3)+B(
     6      ID+I+(IJK-1)*INC3)))+(SIN36*(A(IB+I+(IJK-1)*INC3)-A(IE+I+(
     7      IJK-1)*INC3))-SIN72*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     8      INC3))))
         D(JC+J+(IJK-1)*INC4) = S2*((A(1+I+(IJK-1)*INC3)-COS36*(A(IB+I+(
     1      IJK-1)*INC3)+A(IE+I+(IJK-1)*INC3))+COS72*(A(IC+I+(IJK-1)*
     2      INC3)+A(ID+I+(IJK-1)*INC3)))-(SIN36*(B(IB+I+(IJK-1)*INC3)-B(
     3      IE+I+(IJK-1)*INC3))-SIN72*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-
     4      1)*INC3)))) + C2*((B(1+I+(IJK-1)*INC3)-COS36*(B(IB+I+(IJK-1)
     5      *INC3)+B(IE+I+(IJK-1)*INC3))+COS72*(B(IC+I+(IJK-1)*INC3)+B(
     6      ID+I+(IJK-1)*INC3)))+(SIN36*(A(IB+I+(IJK-1)*INC3)-A(IE+I+(
     7      IJK-1)*INC3))-SIN72*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     8      INC3))))
         C(JD+J+(IJK-1)*INC4) = C3*((A(1+I+(IJK-1)*INC3)-COS36*(A(IB+I+(
     1      IJK-1)*INC3)+A(IE+I+(IJK-1)*INC3))+COS72*(A(IC+I+(IJK-1)*
     2      INC3)+A(ID+I+(IJK-1)*INC3)))+(SIN36*(B(IB+I+(IJK-1)*INC3)-B(
     3      IE+I+(IJK-1)*INC3))-SIN72*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-
     4      1)*INC3)))) - S3*((B(1+I+(IJK-1)*INC3)-COS36*(B(IB+I+(IJK-1)
     5      *INC3)+B(IE+I+(IJK-1)*INC3))+COS72*(B(IC+I+(IJK-1)*INC3)+B(
     6      ID+I+(IJK-1)*INC3)))-(SIN36*(A(IB+I+(IJK-1)*INC3)-A(IE+I+(
     7      IJK-1)*INC3))-SIN72*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     8      INC3))))
         D(JD+J+(IJK-1)*INC4) = S3*((A(1+I+(IJK-1)*INC3)-COS36*(A(IB+I+(
     1      IJK-1)*INC3)+A(IE+I+(IJK-1)*INC3))+COS72*(A(IC+I+(IJK-1)*
     2      INC3)+A(ID+I+(IJK-1)*INC3)))+(SIN36*(B(IB+I+(IJK-1)*INC3)-B(
     3      IE+I+(IJK-1)*INC3))-SIN72*(B(IC+I+(IJK-1)*INC3)-B(ID+I+(IJK-
     4      1)*INC3)))) + C3*((B(1+I+(IJK-1)*INC3)-COS36*(B(IB+I+(IJK-1)
     5      *INC3)+B(IE+I+(IJK-1)*INC3))+COS72*(B(IC+I+(IJK-1)*INC3)+B(
     6      ID+I+(IJK-1)*INC3)))-(SIN36*(A(IB+I+(IJK-1)*INC3)-A(IE+I+(
     7      IJK-1)*INC3))-SIN72*(A(IC+I+(IJK-1)*INC3)-A(ID+I+(IJK-1)*
     8      INC3))))
  145 CONTINUE
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue
      return
c  ** error - factor less than 1  not allowed **
998   print *,' fft99: factors are incorrect '
c     call abort
      stop 'sur ABORT'
c  ** error - factor higher than 5 not allowed **
999   print *,' fft99: factors higher than 5 are not supported '
c     call abort
      stop 'sur ABORT'
      end
      subroutine fax(ifax,n,mode)
c***********************************************************************
c                                                                      *
c c06-summatiom
c c06-summation of series                                      b6.1/3  *
c                                                                      *
c                                                              fftrig  *
c                                                              fax     *
c                                                                      *
c                                                                      *
c sup
c subprogram       subroutine   fftrig                                 *
c                               fax                                    *
c                                                                      *
c purpose          setup routines for fft packages                     *
c                                                                      *
c                                                                      *
c version          cyber                         cray-1                *
c                                                                      *
c                  jan 1979 original             jan 1979 original     *
c                                                                      *
c usage                                                                *
c                  call fftrig(trigs,n,3)                              *
c                  call fax   (ifax ,n,3)                              *
c                                                                      *
c arguments        1.dimension                                         *
c                       trigs(dimension 3*n/2 - add 1 if n/2 is odd)   *
c                       ifax(10)                                       *
c                                                                      *
c                  2.input                                             *
c                      n - the lenght of the transforms to be performed*
c                          n must be even.                             *
c                          the number of words of ifax used increases  *
c                          logarithmically with n.                     *
c                          ifax(10) suffices for practical purposes.   *
c                          (transforms of lenght at least 10000)       *
c                                                                      *
c                  3.output                                            *
c                      trigs - fftrig returns an array of trigonometric*
c                              function values subsequently used by    *
c                              fft routines.                           *
c                      ifax  - fax factorizes n/2 into a product of    *
c                              4"s and 2"s and higher prime numbers.   *
c                              ifax(1) contains the number of factors. *
c                              and the factors themselves are stored   *
c                              in ascending order in ifax(2),ifax(3).. *
c                              if fax is called with n odd ,ifax(1)    *
c                              is set to -99(error condition) and no   *
c                              factorization is done.                  *
c                                                                      *
c write up         none                                                *
c                                                                      *
c entry points           fftrig,  fax                                  *
c                                                                      *
c common blocks    none                                                *
c i/o              none                                                *
c precision        single                                              *
c other routines   none                                                *
c       required                                                       *
c 7/80                     fftrig-1                                    *
c                                                                      *
c***********************************************************************
c                                                                      *
c co6-summation of series                                       b6.1/3 *
c                                                                      *
c                                                              fftrig  *
c                                                              fax     *
c                                                                      *
c acsses (object)  cyber:                                              *
c                           attach,eclib.                              *
c                           ldset(lib=eclib)                           *
c                  cray 1:                                             *
c                           ldr(lib=eclib...)                          *
c                                                                      *
c access (source)           attach,oldpl,eclibpl                       *
c                                                                      *
c                  cyber :         %define cyber                       *
c                  cray:           %define cray                        *
c                                  %c   fftrig,   fax                  *
c                                                                      *
c language         fortran                                             *
c                                                                      *
c specialist       clive temperton                                     *
c                                                                      *
c history          written by c.temperton      jan     1979            *
c                                                                      *
c algorithm                                                            *
c references                                                           *
c                                                                      *
c object size               fftrig  fax  (octal words)                 *
c                  cyber:     145   127                                *
c                  cray :     221   157                                *
c                                                                      *
c                                                                      *
c accuracy                                                             *
c                                                                      *
c timing                                                               *
c                                                                      *
c portability      standard fortran                                    *
c                                                                      *
c system routines  none                                                *
c        required                                                      *
c
c 7/80                      fftrig-2                                   *
c
c***********************************************************************
c     end
      dimension ifax(10)
      nn=n
      if (iabs(mode).eq.1) go to 10
      if (iabs(mode).eq.8) go to 10
      nn=n/2
      if ((nn+nn).eq.n)  go to 10
      ifax(1)=-99
      return
   10 k=1
c     test for factors of 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
c     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
c     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40

c     now find remaining factors
   50 l=5
      inc=2
c     inc alternatively takes on values 2 and 4
   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 l=l+inc
      inc=6-inc
      go to 60
   80 ifax(1)=k-1
c     ifax(1) contains number of factors
      nfax=ifax(1)
c     sort factors into ascending order
      if (nfax.eq.1) go to 110
      do 100 ii=2,nfax
      istop=nfax+2-ii
      do 90 i=2,istop
      if (ifax(i+1).ge.ifax(i)) go to 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 continue
  100 continue
  110 continue
      return
      end
      subroutine fftrig (trigs,n,mode)
      dimension trigs(n)
c    fftrig returns an array of trigonometric function values
c           subsequently used by    f f t    routines
c    see comments in routine    f a x
c     end
      pi=2.0*asin(1.0)
      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/float(nn)
      l=nn+nn
      do 10 i=1,l,2

      angle=0.5*float(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      if (imode.eq.1) return
      if (imode.eq.8) return
      del=0.5*del
      nh=(nn+1)/2
      l=nh+nh
      la=nn+nn
      do 20 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(la+i)=cos(angle)
      trigs(la+i+1)=sin(angle)
   20 continue
      if (imode.le.3) return
      del=0.5*del
      la=la+nn
      if (mode.eq.5) go to 40
      do 30 i=2,nn
      angle=float(i-1)*del
      trigs(la+i)=2.0*sin(angle)
   30 continue
      return
   40 continue
      del=0.5*del
      do 50 i=2,n
      angle=float(i-1)*del
      trigs(la+i)=sin(angle)
   50 continue
      return
      end
************************************************************************
      subroutine fftcc(z,w,ex,ifax,inc,jump,n,nft,isign)
* fft complex-->complex using temperton's real-->complex fft
* same arguments as fft991
* n must be even
      if(isign.eq.1)then
      call trccc(z,w,inc,jump,n,nft,isign)
      call fft991(z,w,ex,ifax,inc,jump,n,nft,isign)
      return
      else
      call fft991(z,w,ex,ifax,inc,jump,n,nft,isign)
      call trccc(z,w,inc,jump,n,nft,isign)
      return
      endif
      end
************************************************************************
      subroutine trccc(f,w,inc,jump,n,nft,isign)
*
* this routine is used to compute complex to complex transforms
*       using temperton's routine fft991
*
*
*
* after ft along x1, one has
*
*    sum   f(x)cos(k1x1)      sum   f(x)sin(k1x1)
*     x1                       x1
*
*--------------------------------------------------------------
* after ft along x2 using fft991 one has
*
*     sum  f(x)cos(k1x1)sin(k2x2)   sum f(x)sin(k1x1)sin(k2x2)
*      x                             x
*
*     sum  f(x)cos(k1x1)cos(k2x2)   sum f(x)sin(k1x1)cos(k2x2)
*      x                             x
*
*---------------------------------------------------------------
* what we want is
*
*     sum  f(x)( cos(k1x1)cos(k2x2)-sin(k1x1)sin(k2x2) )
*      x
*
*     sum  f(x)( cos(k1x1)sin(k2x2)+sin(k1x1)cos(k2x2) )
*      x
*
*---------------------------------------------------------------
      real f(1),w(2,0:1)
      nh=n/2
      nfth=nft/2
      i1=1
      incd=inc*2
      if(isign.eq.-1)then
        do 1 ift=1,nfth
        i2=i1+1
        i1p=i1+inc
        i2p=i2+inc
*VDIR NODEP
        do 3 j=0,nh-1
        w(1,j)=f(i1+incd*j)-f(i2p+incd*j)
        w(2,j)=f(i2+incd*j)+f(i1p+incd*j)
        w(1,n-j)=f(i1+incd*j)+f(i2p+incd*j)
        w(2,n-j)=f(i2+incd*j)-f(i1p+incd*j)
  3     continue
        w(1,nh)=f(i1+incd*nh)+f(i2p+incd*nh)
        w(2,nh)=f(i2+incd*nh)-f(i1p+incd*nh)
*VDIR NODEP
        do 5 j=0,n-1
        f(i1+inc*j)=w(1,j)
        f(i2+inc*j)=w(2,j)
  5     continue
        i1=i1+jump*2
  1     continue
      else
        do 2 ift=1,nfth
        i2=i1+1
        i1n=i1+n*inc
        i2n=i2+n*inc
*VDIR NODEP
        do 4 j=1,nh-1
        w(1,2*j)=.5*(f(i1+inc*j)+f(i1n-inc*j))
        w(2,2*j)=.5*(f(i2+inc*j)+f(i2n-inc*j))
        w(1,2*j+1)=.5*(f(i2+inc*j)-f(i2n-inc*j))
        w(2,2*j+1)=.5*(f(i1n-inc*j)-f(i1+inc*j))
  4     continue
        w(1,0)=f(i1)
        w(2,0)=f(i2)
        w(1,1)=0.
        w(2,1)=0.
        w(1,n)=f(i1n)
        w(2,n)=f(i2n)
        w(1,n+1)=0.
        w(2,n+1)=0.
*VDIR NODEP
        do 6 j=0,n-1
        f(i1+inc*j)=w(1,j)
        f(i2+inc*j)=w(2,j)
  6     continue
        i1=i1+jump*2
  2     continue
      endif
      return
      end
************************************************************************
      subroutine set99(ex,ifax,n)
      call fftrig(ex,n,3)
      call fax(ifax,n,3)
      return
      end
      subroutine fftfax(n,ifax,ex)
      call set99(ex,ifax,n)
      return
      end
      subroutine cftfax(n,ifax,ex)
      call set99(ex,ifax,n)
      return
      end
      subroutine fft999(a,w,ex,ifax,inc,jump,n,nft,isign) 
      dimension a(1),w(1),ex(1),ifax(1)
      call fft991(a,w,ex,ifax,inc,jump,n,nft,isign)
      return
      end
      subroutine cfft999(a,b,w,ex,ifax,inc,jump,n,nft,isign)
      dimension a(1),b(1),w(1),ex(1),ifax(1)
      call fftcc(a,w,ex,ifax,inc,jump/2,n,nft*2,isign)
      return
      end 
c
c***********************************************************************
c
      subroutine fctfax (n,ifax,trigs1,trigs2)
c
c***********************************************************************
c
c  (**    Determine coefficients for FCT99.  **)
c
      dimension trigs1(1),trigs2(1),ifax(1)
      data pi /3.14159265358979323846/
c
c  (**    Get coefficients for fft991.  **)
c  (**    call fftfax if Temperton fft991 used.  **)
c
      call fftfax (n,ifax,trigs1)
c  (**    Get coefficients for rearranging data after fft991.  **)
      do 10 i = 1,n
   10 trigs2(i) = sin((i-1)*pi/n)
      return
      end
c
c******************************************************************
c
      subroutine numcar (num,car)
c
c******************************************************************
c
      implicit real(a-h,o-z)
c
      character*3 car
c
      if (num .ge.100) then
         write (car,1) num
 1       format (i3)
      else
      if (num. ge. 10) then
         write (car,2) num
 2       format ('0',i2)
      else
         write (car,3) num
 3       format ('00',i1)
      endif
      endif
c
      return
      end
