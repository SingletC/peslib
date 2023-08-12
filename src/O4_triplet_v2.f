C   System:                     O4
C   Functional form:            permutation-invariant polynomials
C   Common name:                O4 triplet(adiabatic ground state)
C   Number of derivatives:      1
C   Number of bodies:           4
C   Number of electronic surfaces: 1
C   Interface: Section-2
C
C   Reference - "Potential energy surface of triplet O4," Yuliya 
C   Paukku, Zoltan Varga, and Donald G. Truhlar, to be published 
C   (Feb. 09, 2018)
C
C   Notes:    PES of triplet O4 with special emphasis on
C                  O2 + O2 -> O2 + O + O
C            - Fit is based on a dataset of 10180 points
C            - Mixed-exponential-Gaussian (MEG)
C               variables are applied to describe
C               the long-range interactions
C            - D3 dispersion correction is implemented  
C               for pairwise interactions
C
C     O1--O2
C
C     O3--O4
C
C   Input: X(4),Y(4),Z(4)               in unit of bohr
C   Output: E                           in unit of hartree
C   Output: dEdX(4),dEdY(4),dEdZ(4)     hartree/bohr
C
      subroutine pot(X,Y,Z,E,dEdX,dEdY,dEdZ)

      implicit double precision (a-h,o-z)

C
C Conversion factors 
C Cconv: bohr to Angstrom 
C        1 bohr = 0.52917721092 angstrom
C Econv: kcal/mol to hartree 
C        1 kcal/mol = 0.159360144 * 10**-2 hartree
C Gconv: kcal/(mol*Angstrom) to hartree/bohr
C        1 kcal mol^-1 angstrom^-1 = 0.843297564 * 10**-3 hartree/bohr
C
      double precision Cconv
      double precision Econv
      double precision Gconv
      parameter(Cconv=0.52917721092d0)
      parameter(Econv=0.159360144d-2)
      parameter(Gconv=0.843297564d-3)
C
C Reference energy of infinitely separated O2 + O2 in hartrees 
C  (taken from calculations)
C
      double precision Eref
      parameter(Eref=-299.89131442d0)

      integer i
      double precision E,V
      double precision X(4),Y(4),Z(4),dEdX(4),dEdY(4),dEdZ(4)
      double precision Xcart(12),dVdX(12)

C Convert to local variables
      do i=1,4
        Xcart(3*i-2)=X(i)*Cconv
        Xcart(3*i-1)=Y(i)*Cconv
        Xcart(3*i)=Z(i)*Cconv
      enddo

      call o4pes(Xcart,v,dVdX,1)

C Convert local output to the ones ANT wants

      E = v*Econv + Eref
      do i=1,4
        dEdX(i)=dVdX(3*i-2)*Gconv
        dEdY(i)=dVdX(3*i-1)*Gconv
        dEdZ(i)=dVdX(3*i)*Gconv
      enddo

      end

C
C Local variables used in the O4 PES subroutine
C input coordinate matrix: X(12)in Ang
C                          o1: X(1),X(2),X(3)
C                          O2: X(4),X(5),X(6)
C                          O3: X(7),X(8),X(9)
C                          O4: X(10),X(11),X(12)
C input flag: igrad     igrad=0 energy-only calculation
C                       igrad=1 energy + gradient
C output potential energy:      v    in kcal/mol
C output gradient:              dVdX in kcal/(mol*Ang)
C
      subroutine o4pes(X,v,dVdX,igrad)
***********************************************************************
* Subroutine to calculate the potential energy V and gradient dVdX
* for given Cartesian coordinates X(12)  
* R:            Interatomic bond distance (6)
* V:            Calculated potential energy
* dVdX:         The derivative of V w.r.t. X, dim(12)
* dVdR:         The derivative of V w.r.t. R, dim(6) 
* dPdR:         The derivative of basis functions w.r.t. R
*               dim(6*306)
* dMdR:         The derivative of monomials w.r.t. R
*               dim(6*112)
* dRdX:         The derivative of R w.r.t. X, dim(6*12)
***********************************************************************
      
      implicit double precision (a-h,o-z) 

      integer i,igrad,j,nob,k
      double precision V
      double precision dVdX(12),X(12) 

C      common /coord/    R(6)
C      common /epot/     rMs(6),rM(0:111),P(0:305),C(430),B(430)
C      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:305),dVdR(6),
C     $                  dRdX(6,12),dBdR(6,276)
C      common /msprmt/   a,ab,ra,rb

C Read Cartesian coordinate from input file
      call coord_convt(X)

      if (igrad .le. 1) then
C Call subroutine Evv to evaluate potential energy V
        call evv(V)

        if (igrad .eq. 1) then
C Call EvdVdX to evaluate the derivatives of V w.r.t. X
          call evdvdx(X,dVdX)
        endif
      else
        write (*,*) 'Only igrad = 0, 1 is allowed!'
      endif

      end 
      subroutine coord_convt(X)
***********************************************************************
*  Program to calculate the six interatomic distance 
*  by reading XYZ coordinate
***********************************************************************
      implicit double precision (a-h,o-z) 

      integer i
      double precision X(12)

      common /coord/    R(6)
      
***********************************************************************
*  Now, calculate the interatomic distances
*  r1 = r(O1O2)    r2 = r(O1O3)
*  r3 = r(O1O4)    r4 = r(O2O3)
*  r5 = r(O2O4)    r6 = r(O3O4)       
***********************************************************************

      R(1)=Sqrt((X(4)-X(1))**2 + (X(5)-X(2))**2 + (X(6)-X(3))**2)
      R(2)=Sqrt((X(7)-X(1))**2 + (X(8)-X(2))**2 + (X(9)-X(3))**2)
      R(3)=Sqrt((X(10)-X(1))**2 + (X(11)-X(2))**2 + (X(12)-X(3))**2)
      R(4)=Sqrt((X(4)-X(7))**2 + (X(5)-X(8))**2 + (X(6)-X(9))**2)
      R(5)=Sqrt((X(4)-X(10))**2 + (X(5)-X(11))**2 + (X(6)-X(12))**2)
      R(6)=Sqrt((X(7)-X(10))**2 + (X(8)-X(11))**2 + (X(9)-X(12))**2)
 
      return

      end
      subroutine EvV(V)
***********************************************************************
* Subroutine to evaluate V for given R 
* V(R) = C*P
* C:            Coefficients, stored in 'dim.inc' 
* P:            Basis functions evaluated for given R
* rMs:          rMs(6), six mixed exponetial-Gaussian terms (MEG)
* a:            Nonlinear parameters in Morse terms (Angstrom)
* ab:           Nonlinear parameters in Gauss terms (Angstrom^2)
* re:           Equilibrium bond length(Angstrom)
* nop:          number of points
* nom:          number of monomials
* nob:          number of basis functions (polynomials)
* rM(0:111):    Array to store monomials
* P(0:305):     Array to store polynomials
* B(1:276):     Array to store basis functions
***********************************************************************
      
      implicit double precision (a-h,o-z)

      integer i,j,k
      double precision dist,dv2dr,V,V2,disp,dispdr(6)

      common /coord/    R(6)
      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)

C Calculate the six MEG terms for each point
      call evmorse

C Calculate the monomials for each point by using six MEG terms
      call evmono

C Calculate the polynomials (basis functions) by using monomials
      call evpoly 

C Calculate the basis functions by removing unconnected and 2-body terms
      call evbas

C Initialized v to be 2De 2*120.243 kcal/mol
      v=240.486d0
C Evaluate 2-body interactions
      do i=1,6
        dist=r(i)
        call ev2gm2(dist,v2,dv2dr,4,0)
        v=v+v2
      enddo

C Add D3 dispersion correction
        call d3disp(r,disp,dispdr,0)
        v=v+disp


C Evaluate V by taken the product of C and Basis function array
      do i=1,430
        v=v + c(i)*b(i)
      enddo

C      Write(*,9999) V 
C 9999 Format('The potential energy is ',F20.14,' kcal/mol')

      return

      end 
      subroutine EvdVdX(X,dVdX)
***********************************************************************
* Subroutine to evaluate dRdX for given R and X 
* R:            R(6), 6 bond lengths
* X:            X(12), 12 Cartesian coordinates
* rM(0:111):    Array to store monomials
* P(0:305):     Array to store polynomials
* dVdX:         dVdX(12), derivatives of V w.r.t. Cartesian coordinates 
* dVdR:         dVdR(6), derivatives of V w.r.t.6 bond length
* dRdX:         dRdX(6,12), derivatives of R(6) w.r.t. 12  
*               Cartesian coordinates
***********************************************************************
      
      implicit double precision (a-h,o-z)

      integer i,j
      double precision dVdX(12),X(12)

      common /coord/    R(6)
      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:465),dVdR(6),
     $                  dRdX(6,12),dBdR(6,430)

C Initialize dVdX
      do i=1,12
        dVdX(i)=0.0d0
      enddo

C Call EvdVdR to evaluate dVdR(6)
      Call evdvdr

C Call EvdRdX to evaluate dRdX(6,12)
      Call evdrdx(X)  

C Calculate dVdX by using chain rule: dV/dXi=(dV/dRj)*(dRj/dXi), j=1 to
C 6
      do i=1,12
        do j=1,6
          dVdX(i)=dVdX(i) + dVdR(j)*dRdX(j,i)
        enddo
      enddo

C      write(*,*) 'The 12 dVdX are:'
C      Write(*,9999) (dVdX(i),i=1,12) 
C 9999 Format(1x,3F15.8)

      return
      end 
      subroutine EvMorse
      
      implicit double precision (a-h,o-z)

      integer i

      common /coord/    R(6)
      common /msprmt/   a,ab,ra,rb
      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)
      
C mixed exponential-Gaussian terms = exp(-(r-re)/a-(r-re)^2/ab)
C re:   equlibrium bond length(1.208A)
C a:    nonlinear parameter optimized for diatomic molecules, unit
C Anstrom
C ab:   nonlinear parameter optimized for diatomic molecules, unit
C Anstrom^2
      do i=1,6
         rms(i)=Exp(-(r(i)-ra)/a-((r(i)-rb)**2.0d0)/ab)
      enddo

      end 
      subroutine EvMono
***********************************************************************
*  The subroutine reads six MEG variables(X) and calculates the
*  monomials (M), the polynomial is built from the monomials.
*  For A4 with max. degree 10, the number of monomials is nom.
***********************************************************************

      implicit double precision (a-h,o-z)

      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)

      rm(0) = 1.0d0
      rm(1) = rms(6)
      rm(2) = rms(5)
      rm(3) = rms(4)
      rm(4) = rms(3)
      rm(5) = rms(2)
      rm(6) = rms(1)
      rm(7) = rm(3)*rm(4)
      rm(8) = rm(2)*rm(5)
      rm(9) = rm(1)*rm(6)
      rm(10) = rm(1)*rm(2)
      rm(11) = rm(1)*rm(3)
      rm(12) = rm(2)*rm(3)
      rm(13) = rm(1)*rm(4)
      rm(14) = rm(2)*rm(4)
      rm(15) = rm(1)*rm(5)
      rm(16) = rm(3)*rm(5)
      rm(17) = rm(4)*rm(5)
      rm(18) = rm(2)*rm(6)
      rm(19) = rm(3)*rm(6)
      rm(20) = rm(4)*rm(6)
      rm(21) = rm(5)*rm(6)
      rm(22) = rm(1)*rm(7)
      rm(23) = rm(2)*rm(7)
      rm(24) = rm(1)*rm(8)
      rm(25) = rm(2)*rm(16)
      rm(26) = rm(2)*rm(17)
      rm(27) = rm(3)*rm(17)
      rm(28) = rm(1)*rm(18)
      rm(29) = rm(1)*rm(19)
      rm(30) = rm(1)*rm(20)
      rm(31) = rm(3)*rm(20)
      rm(32) = rm(1)*rm(21)
      rm(33) = rm(2)*rm(21)
      rm(34) = rm(1)*rm(12)
      rm(35) = rm(1)*rm(17)
      rm(36) = rm(2)*rm(20)
      rm(37) = rm(3)*rm(21)
      rm(38) = rm(1)*rm(14)
      rm(39) = rm(1)*rm(16)
      rm(40) = rm(2)*rm(19)
      rm(41) = rm(4)*rm(21)
      rm(42) = rm(2)*rm(27)
      rm(43) = rm(1)*rm(31)
      rm(44) = rm(1)*rm(33)
      rm(45) = rm(1)*rm(23)
      rm(46) = rm(1)*rm(25)
      rm(47) = rm(1)*rm(26)
      rm(48) = rm(1)*rm(27)
      rm(49) = rm(1)*rm(40)
      rm(50) = rm(1)*rm(36)
      rm(51) = rm(2)*rm(31)
      rm(52) = rm(1)*rm(37)
      rm(53) = rm(2)*rm(37)
      rm(54) = rm(1)*rm(41)
      rm(55) = rm(2)*rm(41)
      rm(56) = rm(3)*rm(41)
      rm(57) = rm(1)*rm(42)
      rm(58) = rm(1)*rm(51)
      rm(59) = rm(1)*rm(53)
      rm(60) = rm(1)*rm(55)
      rm(61) = rm(1)*rm(56)
      rm(62) = rm(2)*rm(56)
      rm(63) = rm(1)*rm(62)
      rm(64) = rm(2)*rm(57)
      rm(65) = rm(3)*rm(57)
      rm(66) = rm(4)*rm(57)
      rm(67) = rm(5)*rm(57)
      rm(68) = rm(1)*rm(58)
      rm(69) = rm(3)*rm(58)
      rm(70) = rm(4)*rm(58)
      rm(71) = rm(1)*rm(59)
      rm(72) = rm(2)*rm(59)
      rm(73) = rm(1)*rm(60)
      rm(74) = rm(2)*rm(60)
      rm(75) = rm(1)*rm(61)
      rm(76) = rm(2)*rm(62)
      rm(77) = rm(3)*rm(61)
      rm(78) = rm(3)*rm(62)
      rm(79) = rm(4)*rm(61)
      rm(80) = rm(4)*rm(62)
      rm(81) = rm(5)*rm(59)
      rm(82) = rm(5)*rm(60)
      rm(83) = rm(5)*rm(62)
      rm(84) = rm(6)*rm(58)
      rm(85) = rm(6)*rm(59)
      rm(86) = rm(6)*rm(60)
      rm(87) = rm(6)*rm(61)
      rm(88) = rm(2)*rm(64)
      rm(89) = rm(3)*rm(65)
      rm(90) = rm(4)*rm(66)
      rm(91) = rm(5)*rm(67)
      rm(92) = rm(1)*rm(68)
      rm(93) = rm(3)*rm(69)
      rm(94) = rm(4)*rm(70)
      rm(95) = rm(1)*rm(71)
      rm(96) = rm(2)*rm(72)
      rm(97) = rm(1)*rm(73)
      rm(98) = rm(2)*rm(74)
      rm(99) = rm(1)*rm(75)
      rm(100) = rm(2)*rm(76)
      rm(101) = rm(3)*rm(77)
      rm(102) = rm(3)*rm(78)
      rm(103) = rm(4)*rm(79)
      rm(104) = rm(4)*rm(80)
      rm(105) = rm(5)*rm(81)
      rm(106) = rm(5)*rm(82)
      rm(107) = rm(5)*rm(83)
      rm(108) = rm(6)*rm(84)
      rm(109) = rm(6)*rm(85)
      rm(110) = rm(6)*rm(86)
      rm(111) = rm(6)*rm(87)


      return

      end subroutine EvMono
      subroutine EvPoly
***********************************************************************
*  The subroutine reads monomials (m) and calculates the
*  permutation invariant polynomials (p).
*  For A4 with max. degree 10, the number of polynomials is nob.
***********************************************************************

      implicit double precision (a-h,o-z)

      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)

      p(0) = rm(0)
      p(1) = rm(1) + rm(2) + rm(3) + rm(4) + rm(5) + rm(6)
      p(2) = rm(7) + rm(8) + rm(9)
      p(3) = rm(10) + rm(11) + rm(12) + rm(13) + rm(14) + rm(15) +
     $       rm(16) + rm(17) + rm(18) + rm(19) + rm(20) + rm(21)
      p(4) = p(1)*p(1) - p(3) - p(2) - p(3) - p(2)
      p(5) = rm(22) + rm(23) + rm(24) + rm(25) + rm(26) + rm(27) +
     $       rm(28) + rm(29) + rm(30) + rm(31) + rm(32) + rm(33)
      p(6) = rm(34) + rm(35) + rm(36) + rm(37)
      p(7) = rm(38) + rm(39) + rm(40) + rm(41)
      p(8) = p(1)*p(2) - p(5)
      p(9) = p(1)*p(3) - p(6) - p(7) - p(5) - p(6)
     $      - p(7) - p(5) - p(6) - p(7)
      p(10) = p(1)*p(4) - p(9) - p(8)
      p(11) = rm(42) + rm(43) + rm(44)
      p(12) = rm(45) + rm(46) + rm(47) + rm(48) + rm(49) + rm(50) +
     $       rm(51) + rm(52) + rm(53) + rm(54) + rm(55) + rm(56)
      p(13) = p(2)*p(3) - p(12)
      p(14) = p(1)*p(5) - p(12) - p(11) - p(13) - p(12)
     $      - p(11) - p(11) - p(11)
      p(15) = p(1)*p(6) - p(12)
      p(16) = p(1)*p(7) - p(12)
      p(17) = p(2)*p(2) - p(11) - p(11)
      p(18) = p(3)*p(3) - p(12) - p(11) - p(15) - p(16)
     $      - p(14) - p(12) - p(11) - p(15) - p(16) - p(14)
     $      - p(12) - p(11) - p(12) - p(11)
      p(19) = p(2)*p(4) - p(14)
      p(20) = p(3)*p(4) - p(15) - p(16) - p(13)
      p(21) = p(1)*p(10) - p(20) - p(19)
      p(22) = rm(57) + rm(58) + rm(59) + rm(60) + rm(61) + rm(62)
      p(23) = p(1)*p(11) - p(22)
      p(24) = p(2)*p(6)
      p(25) = p(2)*p(7)
      p(26) = p(1)*p(12) - p(22) - p(24) - p(25) - p(22)
     $      - p(22) - p(22)
      p(27) = p(2)*p(5) - p(22) - p(23) - p(22)
      p(28) = p(3)*p(5) - p(22) - p(26) - p(24) - p(25)
     $      - p(23) - p(22) - p(24) - p(25) - p(23) - p(22)
     $      - p(22)
      p(29) = p(3)*p(6) - p(22) - p(26) - p(22)
      p(30) = p(3)*p(7) - p(22) - p(26) - p(22)
      p(31) = p(2)*p(9) - p(26) - p(28)
      p(32) = p(1)*p(14) - p(26) - p(23) - p(28)
      p(33) = p(4)*p(6) - p(25)
      p(34) = p(4)*p(7) - p(24)
      p(35) = p(1)*p(17) - p(27)
      p(36) = p(1)*p(18) - p(29) - p(30) - p(28)
      p(37) = p(2)*p(10) - p(32)
      p(38) = p(3)*p(10) - p(33) - p(34) - p(31)
      p(39) = p(1)*p(21) - p(38) - p(37)
      p(40) = rm(63)
      p(41) = rm(64) + rm(65) + rm(66) + rm(67) + rm(68) + rm(69) +
     $       rm(70) + rm(71) + rm(72) + rm(73) + rm(74) + rm(75) +
     $       rm(76) + rm(77) + rm(78) + rm(79) + rm(80) + rm(81) +
     $       rm(82) + rm(83) + rm(84) + rm(85) + rm(86) + rm(87)
      p(42) = p(1)*p(22) - p(40) - p(41) - p(40) - p(40)
     $      - p(40) - p(40) - p(40)
      p(43) = p(2)*p(11) - p(40) - p(40) - p(40)
      p(44) = p(2)*p(12) - p(41)
      p(45) = p(3)*p(11) - p(41)
      p(46) = p(5)*p(6) - p(41)
      p(47) = p(5)*p(7) - p(41)
      p(48) = p(6)*p(7) - p(40) - p(40) - p(40) - p(40)
      p(49) = p(4)*p(11) - p(42)
      p(50) = p(2)*p(15) - p(46)
      p(51) = p(2)*p(16) - p(47)
      p(52) = p(4)*p(12) - p(41) - p(50) - p(51)
      p(53) = p(2)*p(14) - p(42) - p(49) - p(42)
      p(54) = p(6)*p(6) - p(42) - p(42)
      p(55) = p(7)*p(7) - p(42) - p(42)
      p(56) = p(3)*p(17) - p(44)
      p(57) = p(2)*p(18) - p(48)
      p(58) = p(3)*p(14) - p(41) - p(52) - p(46) - p(47)
     $      - p(45) - p(45)
      p(59) = p(6)*p(9) - p(41) - p(52) - p(47)
      p(60) = p(7)*p(9) - p(41) - p(52) - p(46)
      p(61) = p(2)*p(20) - p(52) - p(58)
      p(62) = p(1)*p(32) - p(52) - p(49) - p(58)
      p(63) = p(6)*p(10) - p(51)
      p(64) = p(7)*p(10) - p(50)
      p(65) = p(2)*p(17) - p(43)
      p(66) = p(3)*p(18) - p(46) - p(47) - p(45) - p(59)
     $      - p(60) - p(58)
      p(67) = p(2)*p(19) - p(49)
      p(68) = p(1)*p(36) - p(59) - p(60) - p(58) - p(57)
     $      - p(66) - p(66)
      p(69) = p(2)*p(21) - p(62)
      p(70) = p(3)*p(21) - p(63) - p(64) - p(61)
      p(71) = p(1)*p(39) - p(70) - p(69)
      p(72) = p(40)*p(1)
      p(73) = p(2)*p(22) - p(72)
      p(74) = p(6)*p(11)
      p(75) = p(7)*p(11)
      p(76) = p(3)*p(22) - p(72) - p(74) - p(75) - p(72)
     $      - p(72) - p(72)
      p(77) = rm(88) + rm(89) + rm(90) + rm(91) + rm(92) + rm(93) +
     $       rm(94) + rm(95) + rm(96) + rm(97) + rm(98) + rm(99) +
     $       rm(100) + rm(101) + rm(102) + rm(103) + rm(104) + rm(105) +
     $       rm(106) + rm(107) + rm(108) + rm(109) + rm(110) + rm(111)
      p(78) = p(1)*p(42) - p(72) - p(76)
      p(79) = p(5)*p(11) - p(72) - p(73) - p(72)
      p(80) = p(2)*p(26) - p(76) - p(77)
      p(81) = p(6)*p(12) - p(72) - p(76) - p(72)
      p(82) = p(7)*p(12) - p(72) - p(76) - p(72)
      p(83) = p(8)*p(11) - p(72)
      p(84) = p(6)*p(17)
      p(85) = p(7)*p(17)
      p(86) = p(9)*p(11) - p(76) - p(77)
      p(87) = p(2)*p(29) - p(81)
      p(88) = p(6)*p(14) - p(75) - p(75)
      p(89) = p(2)*p(30) - p(82)
      p(90) = p(7)*p(14) - p(74) - p(74)
      p(91) = p(1)*p(48) - p(76) - p(81) - p(82)
      p(92) = p(10)*p(11) - p(78)
      p(93) = p(2)*p(33) - p(88)
      p(94) = p(2)*p(34) - p(90)
      p(95) = p(10)*p(12) - p(77) - p(93) - p(94)
      p(96) = p(2)*p(27) - p(73) - p(79)
      p(97) = p(2)*p(28) - p(76) - p(86)
      p(98) = p(1)*p(53) - p(80) - p(79) - p(97)
      p(99) = p(1)*p(54) - p(81)
      p(100) = p(1)*p(55) - p(82)
      p(101) = p(5)*p(18) - p(76) - p(91) - p(87) - p(89)
     $      - p(86)
      p(102) = p(6)*p(18) - p(74) - p(90)
      p(103) = p(7)*p(18) - p(75) - p(88)
      p(104) = p(2)*p(31) - p(77) - p(86)
      p(105) = p(2)*p(36) - p(91) - p(101)
      p(106) = p(3)*p(32) - p(77) - p(95) - p(88) - p(90)
     $      - p(86)
      p(107) = p(4)*p(29) - p(82) - p(80) - p(99)
      p(108) = p(4)*p(30) - p(81) - p(80) - p(100)
      p(109) = p(2)*p(38) - p(95) - p(106)
      p(110) = p(1)*p(62) - p(95) - p(92) - p(106)
      p(111) = p(6)*p(21) - p(94)
      p(112) = p(7)*p(21) - p(93)
      p(113) = p(1)*p(65) - p(96)
      p(114) = p(1)*p(66) - p(102) - p(103) - p(101)
      p(115) = p(2)*p(37) - p(92)
      p(116) = p(10)*p(18) - p(99) - p(100) - p(97)
      p(117) = p(2)*p(39) - p(110)
      p(118) = p(3)*p(39) - p(111) - p(112) - p(109)
      p(119) = p(1)*p(71) - p(118) - p(117)
      p(120) = p(40)*p(2)
      p(121) = p(40)*p(3)
      p(122) = p(40)*p(4)
      p(123) = p(11)*p(12) - p(121)
      p(124) = p(2)*p(42) - p(122)
      p(125) = p(6)*p(22) - p(121)
      p(126) = p(7)*p(22) - p(121)
      p(127) = p(2)*p(41) - p(121) - p(123) - p(121)
      p(128) = p(6)*p(23) - p(123)
      p(129) = p(7)*p(23) - p(123)
      p(130) = p(3)*p(41) - p(122) - p(121) - p(120) - p(125)
     $      - p(128) - p(126) - p(129) - p(124) - p(123) - p(122)
     $      - p(121) - p(120) - p(125) - p(126) - p(124) - p(123)
     $      - p(122) - p(121) - p(120) - p(122) - p(121) - p(120)
     $      - p(120) - p(120) - p(120) - p(120)
      p(131) = p(3)*p(42) - p(121) - p(125) - p(126) - p(121)
      p(132) = p(4)*p(41) - p(121) - p(131) - p(128) - p(129)
     $      - p(127) - p(121)
      p(133) = p(1)*p(78) - p(122) - p(131)
      p(134) = p(11)*p(11) - p(120) - p(120)
      p(135) = p(2)*p(48) - p(130)
      p(136) = p(11)*p(17) - p(120)
      p(137) = p(2)*p(44) - p(123)
      p(138) = p(2)*p(45) - p(121)
      p(139) = p(11)*p(14) - p(122) - p(124) - p(122)
      p(140) = p(6)*p(27) - p(127)
      p(141) = p(2)*p(54)
      p(142) = p(7)*p(27) - p(127)
      p(143) = p(2)*p(52) - p(131) - p(132)
      p(144) = p(1)*p(81) - p(125) - p(141) - p(135) - p(125)
      p(145) = p(2)*p(55)
      p(146) = p(1)*p(82) - p(126) - p(135) - p(145) - p(126)
      p(147) = p(11)*p(18) - p(130)
      p(148) = p(6)*p(28) - p(129) - p(123) - p(143)
      p(149) = p(7)*p(28) - p(128) - p(123) - p(143)
      p(150) = p(6)*p(30) - p(121) - p(146)
      p(151) = p(11)*p(19) - p(122)
      p(152) = p(2)*p(50) - p(128)
      p(153) = p(2)*p(51) - p(129)
      p(154) = p(11)*p(20) - p(131) - p(132)
      p(155) = p(2)*p(59) - p(144) - p(148)
      p(156) = p(6)*p(32) - p(129)
      p(157) = p(2)*p(60) - p(146) - p(149)
      p(158) = p(7)*p(32) - p(128)
      p(159) = p(6)*p(34) - p(122) - p(145) - p(122)
      p(160) = p(11)*p(21) - p(133)
      p(161) = p(2)*p(63) - p(156)
      p(162) = p(2)*p(64) - p(158)
      p(163) = p(12)*p(21) - p(132) - p(161) - p(162)
      p(164) = p(2)*p(53) - p(124) - p(139)
      p(165) = p(2)*p(58) - p(131) - p(154)
      p(166) = p(3)*p(54) - p(125) - p(144)
      p(167) = p(3)*p(55) - p(126) - p(146)
      p(168) = p(3)*p(65) - p(137)
      p(169) = p(17)*p(18) - p(135)
      p(170) = p(1)*p(98) - p(143) - p(139) - p(165)
      p(171) = p(4)*p(54) - p(135)
      p(172) = p(4)*p(55) - p(135)
      p(173) = p(2)*p(66) - p(150)
      p(174) = p(1)*p(101) - p(148) - p(149) - p(147) - p(173)
     $      - p(165) - p(147)
      p(175) = p(1)*p(102) - p(150) - p(148) - p(166)
      p(176) = p(1)*p(103) - p(150) - p(149) - p(167)
      p(177) = p(2)*p(61) - p(132) - p(154)
      p(178) = p(2)*p(68) - p(159) - p(174)
      p(179) = p(3)*p(62) - p(132) - p(163) - p(156) - p(158)
     $      - p(154)
      p(180) = p(6)*p(38) - p(132) - p(163) - p(157)
      p(181) = p(7)*p(38) - p(132) - p(163) - p(155)
      p(182) = p(2)*p(70) - p(163) - p(179)
      p(183) = p(1)*p(110) - p(163) - p(160) - p(179)
      p(184) = p(6)*p(39) - p(162)
      p(185) = p(7)*p(39) - p(161)
      p(186) = p(2)*p(65) - p(136)
      p(187) = p(3)*p(66) - p(148) - p(149) - p(147) - p(175)
     $      - p(176) - p(174)
      p(188) = p(2)*p(67) - p(151)
      p(189) = p(4)*p(66) - p(166) - p(167) - p(165)
      p(190) = p(2)*p(69) - p(160)
      p(191) = p(18)*p(21) - p(171) - p(172) - p(169)
      p(192) = p(2)*p(71) - p(183)
      p(193) = p(3)*p(71) - p(184) - p(185) - p(182)
      p(194) = p(1)*p(119) - p(193) - p(192)
      p(195) = p(40)*p(5)
      p(196) = p(40)*p(6)
      p(197) = p(40)*p(7)
      p(198) = p(40)*p(8)
      p(199) = p(40)*p(9)
      p(200) = p(40)*p(10)
      p(201) = p(11)*p(22) - p(195)
      p(202) = p(12)*p(22) - p(196) - p(197) - p(195) - p(196)
     $      - p(197) - p(195) - p(196) - p(197)
      p(203) = p(17)*p(22) - p(198)
      p(204) = p(6)*p(43)
      p(205) = p(7)*p(43)
      p(206) = p(11)*p(26) - p(199) - p(202)
      p(207) = p(2)*p(76) - p(199) - p(202)
      p(208) = p(2)*p(78) - p(200)
      p(209) = p(6)*p(41) - p(199) - p(195) - p(202) - p(195)
      p(210) = p(6)*p(42) - p(197) - p(197) - p(197)
      p(211) = p(7)*p(41) - p(199) - p(195) - p(202) - p(195)
      p(212) = p(7)*p(42) - p(196) - p(196) - p(196)
      p(213) = p(11)*p(29) - p(209)
      p(214) = p(11)*p(30) - p(211)
      p(215) = p(18)*p(22) - p(199) - p(213) - p(214)
      p(216) = p(2)*p(77) - p(199) - p(206)
      p(217) = p(6)*p(49) - p(205)
      p(218) = p(7)*p(49) - p(204)
      p(219) = p(3)*p(77) - p(200) - p(199) - p(198) - p(209)
     $      - p(217) - p(211) - p(218) - p(207) - p(204) - p(205)
     $      - p(200) - p(199) - p(198) - p(200) - p(198) - p(200)
     $      - p(198)
      p(220) = p(3)*p(78) - p(199) - p(210) - p(212)
      p(221) = p(10)*p(41) - p(199) - p(220) - p(217) - p(218)
     $      - p(216)
      p(222) = p(1)*p(133) - p(200) - p(220)
      p(223) = p(11)*p(27) - p(195) - p(203)
      p(224) = p(1)*p(134) - p(201)
      p(225) = p(2)*p(81) - p(209)
      p(226) = p(2)*p(80) - p(202) - p(206)
      p(227) = p(2)*p(82) - p(211)
      p(228) = p(2)*p(91) - p(215) - p(219)
      p(229) = p(11)*p(28) - p(199) - p(207)
      p(230) = p(6)*p(53) - p(205)
      p(231) = p(5)*p(54) - p(209)
      p(232) = p(7)*p(53) - p(204)
      p(233) = p(7)*p(54) - p(196)
      p(234) = p(5)*p(55) - p(211)
      p(235) = p(6)*p(55) - p(197)
      p(236) = p(11)*p(35) - p(198)
      p(237) = p(6)*p(65)
      p(238) = p(7)*p(65)
      p(239) = p(2)*p(86) - p(199) - p(229)
      p(240) = p(1)*p(139) - p(206) - p(229) - p(224)
      p(241) = p(17)*p(29) - p(225)
      p(242) = p(2)*p(99) - p(231)
      p(243) = p(17)*p(30) - p(227)
      p(244) = p(2)*p(95) - p(220) - p(221)
      p(245) = p(4)*p(81) - p(202) - p(242) - p(227)
      p(246) = p(2)*p(100) - p(234)
      p(247) = p(4)*p(82) - p(202) - p(225) - p(246)
      p(248) = p(11)*p(36) - p(215) - p(219)
      p(249) = p(2)*p(102) - p(233)
      p(250) = p(6)*p(58) - p(206) - p(214) - p(244) - p(214)
      p(251) = p(2)*p(103) - p(235)
      p(252) = p(7)*p(58) - p(213) - p(206) - p(244) - p(213)
      p(253) = p(1)*p(150) - p(215) - p(233) - p(235)
      p(254) = p(11)*p(37) - p(200)
      p(255) = p(2)*p(93) - p(217)
      p(256) = p(2)*p(94) - p(218)
      p(257) = p(11)*p(38) - p(220) - p(221)
      p(258) = p(2)*p(107) - p(245) - p(250)
      p(259) = p(6)*p(62) - p(218)
      p(260) = p(2)*p(108) - p(247) - p(252)
      p(261) = p(7)*p(62) - p(217)
      p(262) = p(6)*p(64) - p(200) - p(246) - p(200)
      p(263) = p(11)*p(39) - p(222)
      p(264) = p(2)*p(111) - p(259)
      p(265) = p(2)*p(112) - p(261)
      p(266) = p(12)*p(39) - p(221) - p(264) - p(265)
      p(267) = p(2)*p(98) - p(208) - p(240)
      p(268) = p(6)*p(54) - p(210)
      p(269) = p(7)*p(55) - p(212)
      p(270) = p(2)*p(96) - p(203) - p(223)
      p(271) = p(2)*p(97) - p(207) - p(229)
      p(272) = p(2)*p(101) - p(215) - p(248)
      p(273) = p(2)*p(106) - p(220) - p(257)
      p(274) = p(6)*p(59) - p(220) - p(215) - p(209)
      p(275) = p(7)*p(60) - p(220) - p(215) - p(211)
      p(276) = p(5)*p(66) - p(215) - p(253) - p(249) - p(251)
     $      - p(248)
      p(277) = p(6)*p(66) - p(213) - p(252)
      p(278) = p(7)*p(66) - p(214) - p(250)
      p(279) = p(2)*p(104) - p(216) - p(239)
      p(280) = p(2)*p(105) - p(219) - p(248)
      p(281) = p(1)*p(170) - p(244) - p(240) - p(273)
      p(282) = p(10)*p(54) - p(227)
      p(283) = p(10)*p(55) - p(225)
      p(284) = p(2)*p(114) - p(253) - p(276)
      p(285) = p(1)*p(174) - p(250) - p(252) - p(248) - p(276)
     $      - p(273)
      p(286) = p(6)*p(68) - p(217) - p(261) - p(251)
      p(287) = p(7)*p(68) - p(218) - p(259) - p(249)
      p(288) = p(2)*p(109) - p(221) - p(257)
      p(289) = p(2)*p(116) - p(262) - p(285)
      p(290) = p(3)*p(110) - p(221) - p(266) - p(259) - p(261)
     $      - p(257)
      p(291) = p(6)*p(70) - p(221) - p(266) - p(260)
      p(292) = p(7)*p(70) - p(221) - p(266) - p(258)
      p(293) = p(2)*p(118) - p(266) - p(290)
      p(294) = p(1)*p(183) - p(266) - p(263) - p(290)
      p(295) = p(6)*p(71) - p(265)
      p(296) = p(7)*p(71) - p(264)
      p(297) = p(1)*p(186) - p(270)
      p(298) = p(1)*p(187) - p(277) - p(278) - p(276)
      p(299) = p(2)*p(115) - p(254)
      p(300) = p(1)*p(189) - p(286) - p(287) - p(285) - p(284)
     $      - p(298)
      p(301) = p(2)*p(117) - p(263)
      p(302) = p(18)*p(39) - p(282) - p(283) - p(280)
      p(303) = p(2)*p(119) - p(294)
      p(304) = p(3)*p(119) - p(295) - p(296) - p(293)
      p(305) = p(1)*p(194) - p(304) - p(303)
      p(306) = p(40)*p(11)
      p(307) = p(40)*p(12)
      p(308) = p(40)*p(17)
      p(309) = p(40)*p(13)
      p(310) = p(40)*p(14)
      p(311) = p(40)*p(15)
      p(312) = p(40)*p(16)
      p(313) = p(40)*p(18)
      p(314) = p(40)*p(19)
      p(315) = p(40)*p(20)
      p(316) = p(40)*p(21)
      p(317) = p(11)*p(42) - p(310)
      p(318) = p(11)*p(44) - p(307)
      p(319) = p(12)*p(43) - p(309) - p(318)
      p(320) = p(17)*p(42) - p(314)
      p(321) = p(2)*p(125) - p(311)
      p(322) = p(2)*p(126) - p(312)
      p(323) = p(11)*p(48) - p(313)
      p(324) = p(12)*p(42) - p(311) - p(312) - p(307) - p(307)
      p(325) = p(6)*p(79) - p(319)
      p(326) = p(11)*p(54)
      p(327) = p(7)*p(79) - p(319)
      p(328) = p(2)*p(131) - p(315) - p(324)
      p(329) = p(22)*p(29) - p(313) - p(311) - p(326) - p(311)
      p(330) = p(11)*p(55)
      p(331) = p(22)*p(30) - p(313) - p(312) - p(330) - p(312)
      p(332) = p(2)*p(127) - p(309) - p(319)
      p(333) = p(6)*p(83) - p(318)
      p(334) = p(7)*p(83) - p(318)
      p(335) = p(11)*p(52) - p(315) - p(324)
      p(336) = p(2)*p(130) - p(313) - p(323) - p(313)
      p(337) = p(2)*p(133) - p(316)
      p(338) = p(6)*p(77) - p(315) - p(309) - p(322)
      p(339) = p(6)*p(78) - p(312)
      p(340) = p(7)*p(77) - p(315) - p(309) - p(321)
      p(341) = p(7)*p(78) - p(311)
      p(342) = p(11)*p(59) - p(329) - p(338)
      p(343) = p(11)*p(60) - p(331) - p(340)
      p(344) = p(3)*p(130) - p(315) - p(311) - p(312) - p(309)
     $      - p(329) - p(338) - p(331) - p(340) - p(328) - p(321)
     $      - p(322) - p(311) - p(312)
      p(345) = p(18)*p(42) - p(310) - p(326) - p(330) - p(310)
      p(346) = p(2)*p(132) - p(315) - p(335)
      p(347) = p(6)*p(92) - p(334)
      p(348) = p(7)*p(92) - p(333)
      p(349) = p(3)*p(132) - p(316) - p(315) - p(314) - p(338)
     $      - p(347) - p(340) - p(348) - p(336) - p(333) - p(334)
     $      - p(316) - p(315) - p(314) - p(316) - p(314) - p(316)
     $      - p(314)
      p(350) = p(3)*p(133) - p(315) - p(339) - p(341)
      p(351) = p(4)*p(132) - p(315) - p(344) - p(342) - p(343)
     $      - p(332)
      p(352) = p(1)*p(222) - p(316) - p(350)
      p(353) = p(2)*p(134) - p(306)
      p(354) = p(2)*p(135) - p(323)
      p(355) = p(3)*p(134) - p(319)
      p(356) = p(11)*p(53) - p(310) - p(320)
      p(357) = p(2)*p(144) - p(329) - p(338)
      p(358) = p(2)*p(143) - p(324) - p(335)
      p(359) = p(6)*p(81) - p(311) - p(324)
      p(360) = p(2)*p(146) - p(331) - p(340)
      p(361) = p(2)*p(150) - p(344)
      p(362) = p(7)*p(82) - p(312) - p(324)
      p(363) = p(11)*p(65) - p(308)
      p(364) = p(2)*p(137) - p(318)
      p(365) = p(17)*p(45) - p(307)
      p(366) = p(4)*p(134) - p(317)
      p(367) = p(6)*p(96) - p(332)
      p(368) = p(17)*p(54)
      p(369) = p(7)*p(96) - p(332)
      p(370) = p(17)*p(55)
      p(371) = p(2)*p(159) - p(345) - p(349)
      p(372) = p(2)*p(147) - p(313)
      p(373) = p(11)*p(58) - p(315) - p(328)
      p(374) = p(2)*p(148) - p(329) - p(342)
      p(375) = p(6)*p(98) - p(327)
      p(376) = p(2)*p(166) - p(359)
      p(377) = p(14)*p(54) - p(323)
      p(378) = p(2)*p(149) - p(331) - p(343)
      p(379) = p(7)*p(98) - p(325)
      p(380) = p(7)*p(99) - p(311) - p(359)
      p(381) = p(2)*p(167) - p(362)
      p(382) = p(14)*p(55) - p(323)
      p(383) = p(6)*p(100) - p(312) - p(362)
      p(384) = p(11)*p(66) - p(344)
      p(385) = p(6)*p(101) - p(325) - p(343) - p(379)
      p(386) = p(7)*p(101) - p(342) - p(327) - p(375)
      p(387) = p(6)*p(103) - p(313) - p(382)
      p(388) = p(11)*p(67) - p(314)
      p(389) = p(2)*p(152) - p(333)
      p(390) = p(2)*p(153) - p(334)
      p(391) = p(2)*p(154) - p(315) - p(373)
      p(392) = p(1)*p(240) - p(335) - p(373) - p(366)
      p(393) = p(2)*p(155) - p(338) - p(342)
      p(394) = p(2)*p(171) - p(377)
      p(395) = p(2)*p(157) - p(340) - p(343)
      p(396) = p(2)*p(163) - p(350) - p(351)
      p(397) = p(6)*p(95) - p(315) - p(350) - p(340)
      p(398) = p(2)*p(172) - p(382)
      p(399) = p(7)*p(95) - p(315) - p(350) - p(338)
      p(400) = p(11)*p(68) - p(345) - p(349)
      p(401) = p(2)*p(175) - p(380) - p(385)
      p(402) = p(6)*p(106) - p(335) - p(343) - p(396)
      p(403) = p(2)*p(176) - p(383) - p(386)
      p(404) = p(7)*p(106) - p(342) - p(335) - p(396)
      p(405) = p(4)*p(150) - p(328) - p(359) - p(362)
      p(406) = p(11)*p(69) - p(316)
      p(407) = p(2)*p(161) - p(347)
      p(408) = p(2)*p(162) - p(348)
      p(409) = p(11)*p(70) - p(350) - p(351)
      p(410) = p(2)*p(180) - p(397) - p(402)
      p(411) = p(6)*p(110) - p(348)
      p(412) = p(2)*p(181) - p(399) - p(404)
      p(413) = p(7)*p(110) - p(347)
      p(414) = p(6)*p(112) - p(316) - p(398) - p(316)
      p(415) = p(11)*p(71) - p(352)
      p(416) = p(2)*p(184) - p(411)
      p(417) = p(2)*p(185) - p(413)
      p(418) = p(12)*p(71) - p(351) - p(416) - p(417)
      p(419) = p(2)*p(164) - p(320) - p(356)
      p(420) = p(2)*p(165) - p(328) - p(373)
      p(421) = p(2)*p(170) - p(337) - p(392)
      p(422) = p(1)*p(268) - p(359)
      p(423) = p(1)*p(269) - p(362)
      p(424) = p(2)*p(174) - p(345) - p(400)
      p(425) = p(6)*p(102) - p(345) - p(326)
      p(426) = p(7)*p(103) - p(345) - p(330)
      p(427) = p(3)*p(186) - p(364)
      p(428) = p(18)*p(65) - p(354)
      p(429) = p(17)*p(66) - p(361)
      p(430) = p(2)*p(179) - p(350) - p(409)
      p(431) = p(4)*p(166) - p(361) - p(357) - p(422)
      p(432) = p(4)*p(167) - p(361) - p(360) - p(423)
      p(433) = p(2)*p(187) - p(387)
      p(434) = p(14)*p(66) - p(328) - p(405) - p(376) - p(381)
     $      - p(373)
      p(435) = p(1)*p(277) - p(387) - p(385) - p(425)
      p(436) = p(1)*p(278) - p(387) - p(386) - p(426)
      p(437) = p(2)*p(177) - p(346) - p(391)
      p(438) = p(2)*p(178) - p(349) - p(400)
      p(439) = p(1)*p(281) - p(396) - p(392) - p(430)
      p(440) = p(21)*p(54) - p(370)
      p(441) = p(21)*p(55) - p(368)
      p(442) = p(2)*p(189) - p(405) - p(434)
      p(443) = p(1)*p(285) - p(402) - p(404) - p(400) - p(434)
     $      - p(430)
      p(444) = p(6)*p(116) - p(347) - p(413) - p(403)
      p(445) = p(7)*p(116) - p(348) - p(411) - p(401)
      p(446) = p(2)*p(182) - p(351) - p(409)
      p(447) = p(2)*p(191) - p(414) - p(443)
      p(448) = p(3)*p(183) - p(351) - p(418) - p(411) - p(413)
     $      - p(409)
      p(449) = p(6)*p(118) - p(351) - p(418) - p(412)
      p(450) = p(7)*p(118) - p(351) - p(418) - p(410)
      p(451) = p(2)*p(193) - p(418) - p(448)
      p(452) = p(1)*p(294) - p(418) - p(415) - p(448)
      p(453) = p(6)*p(119) - p(417)
      p(454) = p(7)*p(119) - p(416)
      p(455) = p(2)*p(186) - p(363)
      p(456) = p(3)*p(187) - p(385) - p(386) - p(384) - p(435)
     $      - p(436) - p(434)
      p(457) = p(2)*p(188) - p(388)
      p(458) = p(4)*p(187) - p(425) - p(426) - p(424)
      p(459) = p(2)*p(190) - p(406)
      p(460) = p(21)*p(66) - p(422) - p(423) - p(420)
      p(461) = p(2)*p(192) - p(415)
      p(462) = p(18)*p(71) - p(440) - p(441) - p(438)
      p(463) = p(2)*p(194) - p(452)
      p(464) = p(3)*p(194) - p(453) - p(454) - p(451)
      p(465) = p(1)*p(305) - p(464) - p(463)

      return

      end subroutine EvPoly
      subroutine evbas
***********************************************************************
*  This subroutine eliminates the 2-body terms in Bowman's approach.
***********************************************************************

      implicit double precision (a-h,o-z)

      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)
      
      integer i
      double precision b1(466) 

C Pass P(0:465) to BM1(1:466)
      do i=1,466
        b1(i)=p(i-1)
      enddo


C Remove unconnected terms and 2-body terms and pass to B(1:430)
      b(1)=b1(4)

      do i=2,4
        b(i)=b1(i+4)
      enddo

      b(5)=b1(10)

      do i=6,11
        b(i)=b1(i+6)
      enddo

      b(12)=b1(19)
      b(13)=b1(21)

      do i=14,26
        b(i)=b1(i+9)
      enddo

      b(27)=b1(37)
      b(28)=b1(39)

      do i=29,53
        b(i)=b1(i+12)
      enddo

      b(54)=b1(67)
      b(55)=b1(69)
      b(56)=b1(71)

      do i=57,97
        b(i)=b1(i+16)
      enddo

      b(98)=b1(115)
      b(99)=b1(117)
      b(100)=b1(119)
      
      do i=101,166
        b(i)=b1(i+20)
      enddo

      b(167)=b1(188)
      b(168)=b1(190)
      b(169)=b1(192)
      b(170)=b1(194)

      do i=171,272
        b(i)=b1(i+25)
      enddo

      b(273)=b1(299)
      b(274)=b1(301)
      b(275)=b1(303)
      b(276)=b1(305)

      do i=277,425
        b(i)=b1(i+30)
      enddo

      b(426)=b1(457)
      b(427)=b1(459)
      b(428)=b1(461)
      b(429)=b1(463)
      b(430)=b1(465)

      return

      end 
      subroutine ev2gm2(r,v,grad,imol,igrad) 
***********************************************************************
*
* Compute the diatomic potential of ground-state triplet O2
*
* References: J. Chem. Phys. 132, 074307 (2010)
*
* Input:  r      interatomic distance in Angstrom
* Output: V      potential in kcal/mol
*         grad   gradient (kcal/mol)/Angstrom
*
***********************************************************************
      implicit none
      integer imol
      double precision  r
      double precision v, grad
C Parameters of analytical even-tempered Gaussian expansions for the
C ground
C state potential energy curve of O2 CBS+SR+SO+CV. Units: alpha in
C Anstromgs^-2, beta=dimensionless, a_k in milihartree.
      double precision :: alpha,beta,a(0:7)
      integer :: k, igrad
! Original parameters
!      alpha = 0.785d0
!      beta = 1.307d0
!      a(0) = -2388.5641690d0
!      a(1) = 18086.977116d0
!      a(2) = -71760.197585d0
!      a(3) = 154738.09175d0
!      a(4) = -215074.85646d0
!      a(5) = 214799.54567d0
!      a(6) = -148395.42850d0
!      a(7) = 73310.781453d0

! Modified parameters for D3(BJ)
       alpha = 9.439784362354936d-1
       beta =  1.262242998506810d0
       a(0) = -1.488979427684798d3
       a(1) =  1.881435846488955d4
       a(2) = -1.053475425838226d5
       a(3) =  2.755135591229064d5
       a(4) = -4.277588997761775d5
       a(5) =  4.404104009614092d5
       a(6) = -2.946204062950765d5
       a(7) =  1.176861219078620d5

      v=0.d0
      do k=0,7
       v= v + a(k)*dexp(-alpha*beta**k*r**2)
      enddo
C From milihartree to kcal/mol
      v=v*627.509523475149d-3
C Compute the gradient if needed
      if (igrad.eq.1) then
       grad=0.d0
       do k=0,7
         grad=grad-2.d0*a(k)*alpha*beta**k*r*dexp(-alpha*beta**k*r**2)
       enddo
C Convert from milihartree/A to i(kcal/mol)/A
         grad=grad*627.509523475149d-3
      endif
      return
      end

      subroutine EvdVdR
***********************************************************************
* Subroutine to evaluate dVdR for given R 
* dVdR = dV2dR + C*dBdR
* C:            Coefficients, stored in 'dim.inc' 
* P:            Basis functions evaluated for given R
* M:            Monomials evaluated for given R
* dV2dR:        Gradient of 2-body interactions
* dMsdR:        dMsdR(6,6), 6 MEG terms w.r.t. 6 bond lengths
* dMdR:         dMdR(6,nom), nom monomials w.r.t. 6 bond lengths
* dPdR:         dPdR(6,nob), nop polynomial basis functions 
*               w.r.t. 6 bond length
* nom:          number of monomials
* nob:          number of basis functions (polynomials)
* M(nom):       Array to store monomials
* P(nob):       Array to store polynomials
***********************************************************************
      
      implicit double precision (a-h,o-z)
      
      integer i,j
      double precision dist,v2,dv2dr,disp,dispdr(6)

      common /coord/    R(6)
      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:465),dVdR(6),
     $                  dRdX(6,12),dBdR(6,430)

C Initialize dVdR(6)
      do i=1,6
        dVdR(i)=0.0d0
      enddo

C Add dV2dR(i) to dVdR
      do i=1,6
        dist=R(i)
        call ev2gm2(dist,v2,dv2dr,4,1)
        dVdR(i)=dv2dr
      enddo

C Add numerical gradient of D3 dispersion correction
       call d3disp(R,disp,dispdr,1)
      do i=1,6
        dVdR(i)= dVdR(i) + dispdr(i)
      enddo

C Calculate dMEG/dr(6,6) for given R(6)
      call evdmsdr

C Calculate the monomials for each point by using six MEG terms
      call evdmdr

C Calculate the polynomials by using monomials
      call evdpdr 

C Remove 2-body interactions and unconnected terms from polynomials
      call evdbdr

C Evaluate dVdR(6) by taken the product of C(j) and dPdR(i,j)
      do i=1,6      
        do j=1,430
         dVdR(i)=dVdR(i) + c(j)*dBdR(i,j)
        enddo
      enddo

      return
      end 
      subroutine EvdRdX(X)
***********************************************************************
* Subroutine to evaluate dRdX for given R and X 
* R:            R(6), 6 bond lengths
* X:            X(12), 12 Cartesian coordinates
* 
* dMdR:         dMdR(6,nom), nom monomials w.r.t.6 bond length
* dPdR:         dPdR(6,nob), nop polynomial basis functions 
*               w.r.t. 6 bond length
* M(nom):       Array to store monomials
* P(nob):       Array to store polynomials
***********************************************************************
      
      implicit double precision (a-h,o-z)

      integer i,j
      double precision X(12)

      common /coord/    R(6)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:465),dVdR(6),
     $                  dRdX(6,12),dBdR(6,430)

C Initialize dRdX(6,12)
      do i=1,6
        do j=1,12
          dRdX(i,j)=0.0d0
        enddo
      enddo

C Start to calculate the non-zero dRdX
C dr1dx
      dRdX(1,1)=(x(1)-x(4))/r(1)
      dRdX(1,2)=(x(2)-x(5))/r(1)
      dRdX(1,3)=(x(3)-x(6))/r(1)
      dRdX(1,4)=-dRdX(1,1)
      dRdX(1,5)=-dRdX(1,2)
      dRdX(1,6)=-dRdX(1,3)

C dr2dx
      dRdX(2,1)=(x(1)-x(7))/r(2)
      dRdX(2,2)=(x(2)-x(8))/r(2)
      dRdX(2,3)=(x(3)-x(9))/r(2)
      dRdX(2,7)=-dRdX(2,1)
      dRdX(2,8)=-dRdX(2,2)
      dRdX(2,9)=-dRdX(2,3)

C dr3dx
      dRdX(3,1)=(x(1)-x(10))/r(3)
      dRdX(3,2)=(x(2)-x(11))/r(3)
      dRdX(3,3)=(x(3)-x(12))/r(3)
      dRdX(3,10)=-dRdX(3,1)
      dRdX(3,11)=-dRdX(3,2)
      dRdX(3,12)=-dRdX(3,3)

C dr4dx
      dRdX(4,4)=(x(4)-x(7))/r(4)
      dRdX(4,5)=(x(5)-x(8))/r(4)
      dRdX(4,6)=(x(6)-x(9))/r(4)
      dRdX(4,7)=-dRdX(4,4)
      dRdX(4,8)=-dRdX(4,5)
      dRdX(4,9)=-dRdX(4,6)

C dr5dx
      dRdX(5,4)=(x(4)-x(10))/r(5)
      dRdX(5,5)=(x(5)-x(11))/r(5)
      dRdX(5,6)=(x(6)-x(12))/r(5)
      dRdX(5,10)=-dRdX(5,4)
      dRdX(5,11)=-dRdX(5,5)
      dRdX(5,12)=-dRdX(5,6)

C dr6dx
      dRdX(6,7)=(x(7)-x(10))/r(6)
      dRdX(6,8)=(x(8)-x(11))/r(6)
      dRdX(6,9)=(x(9)-x(12))/r(6)
      dRdX(6,10)=-dRdX(6,7)
      dRdX(6,11)=-dRdX(6,8)
      dRdX(6,12)=-dRdX(6,9)
C Finish the calculation of non-zero dRdX

      return

      end 
      subroutine EvdMsdR
***********************************************************************
* Subroutine to evaluate the derivatives of MEG term X
* w.r.t. interatomic distance R(6)
* dmsdR:        Local variables, dirm(6,6)
* a:            Nonlinear pamameter(Angstrom)
* ab:       Nonlinear pamameter(Angstrom^2)
* re:           equilibrium bond length(Angstrom)
***********************************************************************
      
      implicit double precision (a-h,o-z)

      integer i,j

      common /coord/    R(6)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:465),dVdR(6),
     $                  dRdX(6,12),dBdR(6,430)
      common /msprmt/   a,ab,ra,rb

C Initialize dmsdr
      do i=1,6
        do j=1,6
          dmsdr(i,j)=0.0d0
        enddo
      enddo
C
C MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
C dmsdr(i,j)=0  i!=j
C
      do i=1,6
         dmsdr(i,i)=(-2.0d0*(r(i)-rb)/ab-1/a)*
     $Exp(-(r(i)-ra)/a-((r(i)-rb)**2.0d0)/ab)
                 
      enddo 

      return

      end 
      subroutine EvdMdR
***********************************************************************
*  The subroutine reads M(nom) and dMSdR(6,6) and calculates the
*  dMdR(6,nom), the polynomial is built from the monomials.
*  For A4 with max. degree 10, the number of monomials is nom.
***********************************************************************

      implicit double precision (a-h,o-z)

      integer i

      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:465),dVdR(6),
     $                  dRdX(6,12),dBdR(6,430)

      do i=1,6
        dmdr(i,0) = 0.0d0
        dmdr(i,1) = dmsdr(i,6)
        dmdr(i,2) = dmsdr(i,5)
        dmdr(i,3) = dmsdr(i,4)
        dmdr(i,4) = dmsdr(i,3)
        dmdr(i,5) = dmsdr(i,2)
        dmdr(i,6) = dmsdr(i,1)
        dmdr(i,7) = dmdr(i,3)*rm(4) + rm(3)*dmdr(i,4)
        dmdr(i,8) = dmdr(i,2)*rm(5) + rm(2)*dmdr(i,5)
        dmdr(i,9) = dmdr(i,1)*rm(6) + rm(1)*dmdr(i,6)
        dmdr(i,10) = dmdr(i,1)*rm(2) + rm(1)*dmdr(i,2)
        dmdr(i,11) = dmdr(i,1)*rm(3) + rm(1)*dmdr(i,3)
        dmdr(i,12) = dmdr(i,2)*rm(3) + rm(2)*dmdr(i,3)
        dmdr(i,13) = dmdr(i,1)*rm(4) + rm(1)*dmdr(i,4)
        dmdr(i,14) = dmdr(i,2)*rm(4) + rm(2)*dmdr(i,4)
        dmdr(i,15) = dmdr(i,1)*rm(5) + rm(1)*dmdr(i,5)
        dmdr(i,16) = dmdr(i,3)*rm(5) + rm(3)*dmdr(i,5)
        dmdr(i,17) = dmdr(i,4)*rm(5) + rm(4)*dmdr(i,5)
        dmdr(i,18) = dmdr(i,2)*rm(6) + rm(2)*dmdr(i,6)
        dmdr(i,19) = dmdr(i,3)*rm(6) + rm(3)*dmdr(i,6)
        dmdr(i,20) = dmdr(i,4)*rm(6) + rm(4)*dmdr(i,6)
        dmdr(i,21) = dmdr(i,5)*rm(6) + rm(5)*dmdr(i,6)
        dmdr(i,22) = dmdr(i,1)*rm(7) + rm(1)*dmdr(i,7)
        dmdr(i,23) = dmdr(i,2)*rm(7) + rm(2)*dmdr(i,7)
        dmdr(i,24) = dmdr(i,1)*rm(8) + rm(1)*dmdr(i,8)
        dmdr(i,25) = dmdr(i,2)*rm(16) + rm(2)*dmdr(i,16)
        dmdr(i,26) = dmdr(i,2)*rm(17) + rm(2)*dmdr(i,17)
        dmdr(i,27) = dmdr(i,3)*rm(17) + rm(3)*dmdr(i,17)
        dmdr(i,28) = dmdr(i,1)*rm(18) + rm(1)*dmdr(i,18)
        dmdr(i,29) = dmdr(i,1)*rm(19) + rm(1)*dmdr(i,19)
        dmdr(i,30) = dmdr(i,1)*rm(20) + rm(1)*dmdr(i,20)
        dmdr(i,31) = dmdr(i,3)*rm(20) + rm(3)*dmdr(i,20)
        dmdr(i,32) = dmdr(i,1)*rm(21) + rm(1)*dmdr(i,21)
        dmdr(i,33) = dmdr(i,2)*rm(21) + rm(2)*dmdr(i,21)
        dmdr(i,34) = dmdr(i,1)*rm(12) + rm(1)*dmdr(i,12)
        dmdr(i,35) = dmdr(i,1)*rm(17) + rm(1)*dmdr(i,17)
        dmdr(i,36) = dmdr(i,2)*rm(20) + rm(2)*dmdr(i,20)
        dmdr(i,37) = dmdr(i,3)*rm(21) + rm(3)*dmdr(i,21)
        dmdr(i,38) = dmdr(i,1)*rm(14) + rm(1)*dmdr(i,14)
        dmdr(i,39) = dmdr(i,1)*rm(16) + rm(1)*dmdr(i,16)
        dmdr(i,40) = dmdr(i,2)*rm(19) + rm(2)*dmdr(i,19)
        dmdr(i,41) = dmdr(i,4)*rm(21) + rm(4)*dmdr(i,21)
        dmdr(i,42) = dmdr(i,2)*rm(27) + rm(2)*dmdr(i,27)
        dmdr(i,43) = dmdr(i,1)*rm(31) + rm(1)*dmdr(i,31)
        dmdr(i,44) = dmdr(i,1)*rm(33) + rm(1)*dmdr(i,33)
        dmdr(i,45) = dmdr(i,1)*rm(23) + rm(1)*dmdr(i,23)
        dmdr(i,46) = dmdr(i,1)*rm(25) + rm(1)*dmdr(i,25)
        dmdr(i,47) = dmdr(i,1)*rm(26) + rm(1)*dmdr(i,26)
        dmdr(i,48) = dmdr(i,1)*rm(27) + rm(1)*dmdr(i,27)
        dmdr(i,49) = dmdr(i,1)*rm(40) + rm(1)*dmdr(i,40)
        dmdr(i,50) = dmdr(i,1)*rm(36) + rm(1)*dmdr(i,36)
        dmdr(i,51) = dmdr(i,2)*rm(31) + rm(2)*dmdr(i,31)
        dmdr(i,52) = dmdr(i,1)*rm(37) + rm(1)*dmdr(i,37)
        dmdr(i,53) = dmdr(i,2)*rm(37) + rm(2)*dmdr(i,37)
        dmdr(i,54) = dmdr(i,1)*rm(41) + rm(1)*dmdr(i,41)
        dmdr(i,55) = dmdr(i,2)*rm(41) + rm(2)*dmdr(i,41)
        dmdr(i,56) = dmdr(i,3)*rm(41) + rm(3)*dmdr(i,41)
        dmdr(i,57) = dmdr(i,1)*rm(42) + rm(1)*dmdr(i,42)
        dmdr(i,58) = dmdr(i,1)*rm(51) + rm(1)*dmdr(i,51)
        dmdr(i,59) = dmdr(i,1)*rm(53) + rm(1)*dmdr(i,53)
        dmdr(i,60) = dmdr(i,1)*rm(55) + rm(1)*dmdr(i,55)
        dmdr(i,61) = dmdr(i,1)*rm(56) + rm(1)*dmdr(i,56)
        dmdr(i,62) = dmdr(i,2)*rm(56) + rm(2)*dmdr(i,56)
        dmdr(i,63) = dmdr(i,1)*rm(62) + rm(1)*dmdr(i,62)
        dmdr(i,64) = dmdr(i,2)*rm(57) + rm(2)*dmdr(i,57)
        dmdr(i,65) = dmdr(i,3)*rm(57) + rm(3)*dmdr(i,57)
        dmdr(i,66) = dmdr(i,4)*rm(57) + rm(4)*dmdr(i,57)
        dmdr(i,67) = dmdr(i,5)*rm(57) + rm(5)*dmdr(i,57)
        dmdr(i,68) = dmdr(i,1)*rm(58) + rm(1)*dmdr(i,58)
        dmdr(i,69) = dmdr(i,3)*rm(58) + rm(3)*dmdr(i,58)
        dmdr(i,70) = dmdr(i,4)*rm(58) + rm(4)*dmdr(i,58)
        dmdr(i,71) = dmdr(i,1)*rm(59) + rm(1)*dmdr(i,59)
        dmdr(i,72) = dmdr(i,2)*rm(59) + rm(2)*dmdr(i,59)
        dmdr(i,73) = dmdr(i,1)*rm(60) + rm(1)*dmdr(i,60)
        dmdr(i,74) = dmdr(i,2)*rm(60) + rm(2)*dmdr(i,60)
        dmdr(i,75) = dmdr(i,1)*rm(61) + rm(1)*dmdr(i,61)
        dmdr(i,76) = dmdr(i,2)*rm(62) + rm(2)*dmdr(i,62)
        dmdr(i,77) = dmdr(i,3)*rm(61) + rm(3)*dmdr(i,61)
        dmdr(i,78) = dmdr(i,3)*rm(62) + rm(3)*dmdr(i,62)
        dmdr(i,79) = dmdr(i,4)*rm(61) + rm(4)*dmdr(i,61)
        dmdr(i,80) = dmdr(i,4)*rm(62) + rm(4)*dmdr(i,62)
        dmdr(i,81) = dmdr(i,5)*rm(59) + rm(5)*dmdr(i,59)
        dmdr(i,82) = dmdr(i,5)*rm(60) + rm(5)*dmdr(i,60)
        dmdr(i,83) = dmdr(i,5)*rm(62) + rm(5)*dmdr(i,62)
        dmdr(i,84) = dmdr(i,6)*rm(58) + rm(6)*dmdr(i,58)
        dmdr(i,85) = dmdr(i,6)*rm(59) + rm(6)*dmdr(i,59)
        dmdr(i,86) = dmdr(i,6)*rm(60) + rm(6)*dmdr(i,60)
        dmdr(i,87) = dmdr(i,6)*rm(61) + rm(6)*dmdr(i,61)
        dmdr(i,88) = dmdr(i,2)*rm(64) + rm(2)*dmdr(i,64)
        dmdr(i,89) = dmdr(i,3)*rm(65) + rm(3)*dmdr(i,65)
        dmdr(i,90) = dmdr(i,4)*rm(66) + rm(4)*dmdr(i,66)
        dmdr(i,91) = dmdr(i,5)*rm(67) + rm(5)*dmdr(i,67)
        dmdr(i,92) = dmdr(i,1)*rm(68) + rm(1)*dmdr(i,68)
        dmdr(i,93) = dmdr(i,3)*rm(69) + rm(3)*dmdr(i,69)
        dmdr(i,94) = dmdr(i,4)*rm(70) + rm(4)*dmdr(i,70)
        dmdr(i,95) = dmdr(i,1)*rm(71) + rm(1)*dmdr(i,71)
        dmdr(i,96) = dmdr(i,2)*rm(72) + rm(2)*dmdr(i,72)
        dmdr(i,97) = dmdr(i,1)*rm(73) + rm(1)*dmdr(i,73)
        dmdr(i,98) = dmdr(i,2)*rm(74) + rm(2)*dmdr(i,74)
        dmdr(i,99) = dmdr(i,1)*rm(75) + rm(1)*dmdr(i,75)
        dmdr(i,100) = dmdr(i,2)*rm(76) + rm(2)*dmdr(i,76)
        dmdr(i,101) = dmdr(i,3)*rm(77) + rm(3)*dmdr(i,77)
        dmdr(i,102) = dmdr(i,3)*rm(78) + rm(3)*dmdr(i,78)
        dmdr(i,103) = dmdr(i,4)*rm(79) + rm(4)*dmdr(i,79)
        dmdr(i,104) = dmdr(i,4)*rm(80) + rm(4)*dmdr(i,80)
        dmdr(i,105) = dmdr(i,5)*rm(81) + rm(5)*dmdr(i,81)
        dmdr(i,106) = dmdr(i,5)*rm(82) + rm(5)*dmdr(i,82)
        dmdr(i,107) = dmdr(i,5)*rm(83) + rm(5)*dmdr(i,83)
        dmdr(i,108) = dmdr(i,6)*rm(84) + rm(6)*dmdr(i,84)
        dmdr(i,109) = dmdr(i,6)*rm(85) + rm(6)*dmdr(i,85)
        dmdr(i,110) = dmdr(i,6)*rm(86) + rm(6)*dmdr(i,86)
        dmdr(i,111) = dmdr(i,6)*rm(87) + rm(6)*dmdr(i,87)
      enddo

      return

      end subroutine EvdMdR
      subroutine EvdPdr
***********************************************************************
*  The subroutine reads monomials (m) and calculates the
*  permutation invariant polynomials (p)
*  For A4 with max. degree 10, the number of polynomials is nob.
***********************************************************************

      implicit double precision (a-h,o-z)

      integer i

      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:465),dVdR(6),
     $                  dRdX(6,12),dBdR(6,430)

      do i=1,6
      dpdr(i,0) = dmdr(i,0)
      dpdr(i,1) = dmdr(i,1) + dmdr(i,2) + dmdr(i,3) + dmdr(i,4) +
     $       dmdr(i,5) + dmdr(i,6)
      dpdr(i,2) = dmdr(i,7) + dmdr(i,8) + dmdr(i,9)
      dpdr(i,3) = dmdr(i,10) + dmdr(i,11) + dmdr(i,12) + dmdr(i,13) +
     $       dmdr(i,14) + dmdr(i,15) + dmdr(i,16) + dmdr(i,17) +
     $       dmdr(i,18) + dmdr(i,19) + dmdr(i,20) + dmdr(i,21)
      dpdr(i,4) = dpdr(i,1)*p(1) + p(1)*dpdr(i,1)
     $      - dpdr(i,3) - dpdr(i,2) - dpdr(i,3) - dpdr(i,2)
      dpdr(i,5) = dmdr(i,22) + dmdr(i,23) + dmdr(i,24) + dmdr(i,25) +
     $       dmdr(i,26) + dmdr(i,27) + dmdr(i,28) + dmdr(i,29) +
     $       dmdr(i,30) + dmdr(i,31) + dmdr(i,32) + dmdr(i,33)
      dpdr(i,6) = dmdr(i,34) + dmdr(i,35) + dmdr(i,36) + dmdr(i,37)
      dpdr(i,7) = dmdr(i,38) + dmdr(i,39) + dmdr(i,40) + dmdr(i,41)
      dpdr(i,8) = dpdr(i,1)*p(2) + p(1)*dpdr(i,2)
     $      - dpdr(i,5)
      dpdr(i,9) = dpdr(i,1)*p(3) + p(1)*dpdr(i,3)
     $      - dpdr(i,6) - dpdr(i,7) - dpdr(i,5) - dpdr(i,6)
     $      - dpdr(i,7) - dpdr(i,5) - dpdr(i,6) - dpdr(i,7)
      dpdr(i,10) = dpdr(i,1)*p(4) + p(1)*dpdr(i,4)
     $      - dpdr(i,9) - dpdr(i,8)
      dpdr(i,11) = dmdr(i,42) + dmdr(i,43) + dmdr(i,44)
      dpdr(i,12) = dmdr(i,45) + dmdr(i,46) + dmdr(i,47) + dmdr(i,48) +
     $       dmdr(i,49) + dmdr(i,50) + dmdr(i,51) + dmdr(i,52) +
     $       dmdr(i,53) + dmdr(i,54) + dmdr(i,55) + dmdr(i,56)
      dpdr(i,13) = dpdr(i,2)*p(3) + p(2)*dpdr(i,3)
     $      - dpdr(i,12)
      dpdr(i,14) = dpdr(i,1)*p(5) + p(1)*dpdr(i,5)
     $      - dpdr(i,12) - dpdr(i,11) - dpdr(i,13) - dpdr(i,12)
     $      - dpdr(i,11) - dpdr(i,11) - dpdr(i,11)
      dpdr(i,15) = dpdr(i,1)*p(6) + p(1)*dpdr(i,6)
     $      - dpdr(i,12)
      dpdr(i,16) = dpdr(i,1)*p(7) + p(1)*dpdr(i,7)
     $      - dpdr(i,12)
      dpdr(i,17) = dpdr(i,2)*p(2) + p(2)*dpdr(i,2)
     $      - dpdr(i,11) - dpdr(i,11)
      dpdr(i,18) = dpdr(i,3)*p(3) + p(3)*dpdr(i,3)
     $      - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) - dpdr(i,16)
     $      - dpdr(i,14) - dpdr(i,12) - dpdr(i,11) - dpdr(i,15)
     $      - dpdr(i,16) - dpdr(i,14) - dpdr(i,12) - dpdr(i,11)
     $      - dpdr(i,12) - dpdr(i,11)
      dpdr(i,19) = dpdr(i,2)*p(4) + p(2)*dpdr(i,4)
     $      - dpdr(i,14)
      dpdr(i,20) = dpdr(i,3)*p(4) + p(3)*dpdr(i,4)
     $      - dpdr(i,15) - dpdr(i,16) - dpdr(i,13)
      dpdr(i,21) = dpdr(i,1)*p(10) + p(1)*dpdr(i,10)
     $      - dpdr(i,20) - dpdr(i,19)
      dpdr(i,22) = dmdr(i,57) + dmdr(i,58) + dmdr(i,59) + dmdr(i,60) +
     $       dmdr(i,61) + dmdr(i,62)
      dpdr(i,23) = dpdr(i,1)*p(11) + p(1)*dpdr(i,11)
     $      - dpdr(i,22)
      dpdr(i,24) = dpdr(i,2)*p(6) + p(2)*dpdr(i,6)
      dpdr(i,25) = dpdr(i,2)*p(7) + p(2)*dpdr(i,7)
      dpdr(i,26) = dpdr(i,1)*p(12) + p(1)*dpdr(i,12)
     $      - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) - dpdr(i,22)
     $      - dpdr(i,22) - dpdr(i,22)
      dpdr(i,27) = dpdr(i,2)*p(5) + p(2)*dpdr(i,5)
     $      - dpdr(i,22) - dpdr(i,23) - dpdr(i,22)
      dpdr(i,28) = dpdr(i,3)*p(5) + p(3)*dpdr(i,5)
     $      - dpdr(i,22) - dpdr(i,26) - dpdr(i,24) - dpdr(i,25)
     $      - dpdr(i,23) - dpdr(i,22) - dpdr(i,24) - dpdr(i,25)
     $      - dpdr(i,23) - dpdr(i,22) - dpdr(i,22)
      dpdr(i,29) = dpdr(i,3)*p(6) + p(3)*dpdr(i,6)
     $      - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
      dpdr(i,30) = dpdr(i,3)*p(7) + p(3)*dpdr(i,7)
     $      - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
      dpdr(i,31) = dpdr(i,2)*p(9) + p(2)*dpdr(i,9)
     $      - dpdr(i,26) - dpdr(i,28)
      dpdr(i,32) = dpdr(i,1)*p(14) + p(1)*dpdr(i,14)
     $      - dpdr(i,26) - dpdr(i,23) - dpdr(i,28)
      dpdr(i,33) = dpdr(i,4)*p(6) + p(4)*dpdr(i,6)
     $      - dpdr(i,25)
      dpdr(i,34) = dpdr(i,4)*p(7) + p(4)*dpdr(i,7)
     $      - dpdr(i,24)
      dpdr(i,35) = dpdr(i,1)*p(17) + p(1)*dpdr(i,17)
     $      - dpdr(i,27)
      dpdr(i,36) = dpdr(i,1)*p(18) + p(1)*dpdr(i,18)
     $      - dpdr(i,29) - dpdr(i,30) - dpdr(i,28)
      dpdr(i,37) = dpdr(i,2)*p(10) + p(2)*dpdr(i,10)
     $      - dpdr(i,32)
      dpdr(i,38) = dpdr(i,3)*p(10) + p(3)*dpdr(i,10)
     $      - dpdr(i,33) - dpdr(i,34) - dpdr(i,31)
      dpdr(i,39) = dpdr(i,1)*p(21) + p(1)*dpdr(i,21)
     $      - dpdr(i,38) - dpdr(i,37)
      dpdr(i,40) = dmdr(i,63)
      dpdr(i,41) = dmdr(i,64) + dmdr(i,65) + dmdr(i,66) + dmdr(i,67) +
     $       dmdr(i,68) + dmdr(i,69) + dmdr(i,70) + dmdr(i,71) +
     $       dmdr(i,72) + dmdr(i,73) + dmdr(i,74) + dmdr(i,75) +
     $       dmdr(i,76) + dmdr(i,77) + dmdr(i,78) + dmdr(i,79) +
     $       dmdr(i,80) + dmdr(i,81) + dmdr(i,82) + dmdr(i,83) +
     $       dmdr(i,84) + dmdr(i,85) + dmdr(i,86) + dmdr(i,87)
      dpdr(i,42) = dpdr(i,1)*p(22) + p(1)*dpdr(i,22)
     $      - dpdr(i,40) - dpdr(i,41) - dpdr(i,40) - dpdr(i,40)
     $      - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
      dpdr(i,43) = dpdr(i,2)*p(11) + p(2)*dpdr(i,11)
     $      - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
      dpdr(i,44) = dpdr(i,2)*p(12) + p(2)*dpdr(i,12)
     $      - dpdr(i,41)
      dpdr(i,45) = dpdr(i,3)*p(11) + p(3)*dpdr(i,11)
     $      - dpdr(i,41)
      dpdr(i,46) = dpdr(i,5)*p(6) + p(5)*dpdr(i,6)
     $      - dpdr(i,41)
      dpdr(i,47) = dpdr(i,5)*p(7) + p(5)*dpdr(i,7)
     $      - dpdr(i,41)
      dpdr(i,48) = dpdr(i,6)*p(7) + p(6)*dpdr(i,7)
     $      - dpdr(i,40) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
      dpdr(i,49) = dpdr(i,4)*p(11) + p(4)*dpdr(i,11)
     $      - dpdr(i,42)
      dpdr(i,50) = dpdr(i,2)*p(15) + p(2)*dpdr(i,15)
     $      - dpdr(i,46)
      dpdr(i,51) = dpdr(i,2)*p(16) + p(2)*dpdr(i,16)
     $      - dpdr(i,47)
      dpdr(i,52) = dpdr(i,4)*p(12) + p(4)*dpdr(i,12)
     $      - dpdr(i,41) - dpdr(i,50) - dpdr(i,51)
      dpdr(i,53) = dpdr(i,2)*p(14) + p(2)*dpdr(i,14)
     $      - dpdr(i,42) - dpdr(i,49) - dpdr(i,42)
      dpdr(i,54) = dpdr(i,6)*p(6) + p(6)*dpdr(i,6)
     $      - dpdr(i,42) - dpdr(i,42)
      dpdr(i,55) = dpdr(i,7)*p(7) + p(7)*dpdr(i,7)
     $      - dpdr(i,42) - dpdr(i,42)
      dpdr(i,56) = dpdr(i,3)*p(17) + p(3)*dpdr(i,17)
     $      - dpdr(i,44)
      dpdr(i,57) = dpdr(i,2)*p(18) + p(2)*dpdr(i,18)
     $      - dpdr(i,48)
      dpdr(i,58) = dpdr(i,3)*p(14) + p(3)*dpdr(i,14)
     $      - dpdr(i,41) - dpdr(i,52) - dpdr(i,46) - dpdr(i,47)
     $      - dpdr(i,45) - dpdr(i,45)
      dpdr(i,59) = dpdr(i,6)*p(9) + p(6)*dpdr(i,9)
     $      - dpdr(i,41) - dpdr(i,52) - dpdr(i,47)
      dpdr(i,60) = dpdr(i,7)*p(9) + p(7)*dpdr(i,9)
     $      - dpdr(i,41) - dpdr(i,52) - dpdr(i,46)
      dpdr(i,61) = dpdr(i,2)*p(20) + p(2)*dpdr(i,20)
     $      - dpdr(i,52) - dpdr(i,58)
      dpdr(i,62) = dpdr(i,1)*p(32) + p(1)*dpdr(i,32)
     $      - dpdr(i,52) - dpdr(i,49) - dpdr(i,58)
      dpdr(i,63) = dpdr(i,6)*p(10) + p(6)*dpdr(i,10)
     $      - dpdr(i,51)
      dpdr(i,64) = dpdr(i,7)*p(10) + p(7)*dpdr(i,10)
     $      - dpdr(i,50)
      dpdr(i,65) = dpdr(i,2)*p(17) + p(2)*dpdr(i,17)
     $      - dpdr(i,43)
      dpdr(i,66) = dpdr(i,3)*p(18) + p(3)*dpdr(i,18)
     $      - dpdr(i,46) - dpdr(i,47) - dpdr(i,45) - dpdr(i,59)
     $      - dpdr(i,60) - dpdr(i,58)
      dpdr(i,67) = dpdr(i,2)*p(19) + p(2)*dpdr(i,19)
     $      - dpdr(i,49)
      dpdr(i,68) = dpdr(i,1)*p(36) + p(1)*dpdr(i,36)
     $      - dpdr(i,59) - dpdr(i,60) - dpdr(i,58) - dpdr(i,57)
     $      - dpdr(i,66) - dpdr(i,66)
      dpdr(i,69) = dpdr(i,2)*p(21) + p(2)*dpdr(i,21)
     $      - dpdr(i,62)
      dpdr(i,70) = dpdr(i,3)*p(21) + p(3)*dpdr(i,21)
     $      - dpdr(i,63) - dpdr(i,64) - dpdr(i,61)
      dpdr(i,71) = dpdr(i,1)*p(39) + p(1)*dpdr(i,39)
     $      - dpdr(i,70) - dpdr(i,69)
      dpdr(i,72) = dpdr(i,40)*p(1) + p(40)*dpdr(i,1)
      dpdr(i,73) = dpdr(i,2)*p(22) + p(2)*dpdr(i,22)
     $      - dpdr(i,72)
      dpdr(i,74) = dpdr(i,6)*p(11) + p(6)*dpdr(i,11)
      dpdr(i,75) = dpdr(i,7)*p(11) + p(7)*dpdr(i,11)
      dpdr(i,76) = dpdr(i,3)*p(22) + p(3)*dpdr(i,22)
     $      - dpdr(i,72) - dpdr(i,74) - dpdr(i,75) - dpdr(i,72)
     $      - dpdr(i,72) - dpdr(i,72)
      dpdr(i,77) = dmdr(i,88) + dmdr(i,89) + dmdr(i,90) + dmdr(i,91) +
     $       dmdr(i,92) + dmdr(i,93) + dmdr(i,94) + dmdr(i,95) +
     $       dmdr(i,96) + dmdr(i,97) + dmdr(i,98) + dmdr(i,99) +
     $       dmdr(i,100) + dmdr(i,101) + dmdr(i,102) + dmdr(i,103) +
     $       dmdr(i,104) + dmdr(i,105) + dmdr(i,106) + dmdr(i,107) +
     $       dmdr(i,108) + dmdr(i,109) + dmdr(i,110) + dmdr(i,111)
      dpdr(i,78) = dpdr(i,1)*p(42) + p(1)*dpdr(i,42)
     $      - dpdr(i,72) - dpdr(i,76)
      dpdr(i,79) = dpdr(i,5)*p(11) + p(5)*dpdr(i,11)
     $      - dpdr(i,72) - dpdr(i,73) - dpdr(i,72)
      dpdr(i,80) = dpdr(i,2)*p(26) + p(2)*dpdr(i,26)
     $      - dpdr(i,76) - dpdr(i,77)
      dpdr(i,81) = dpdr(i,6)*p(12) + p(6)*dpdr(i,12)
     $      - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
      dpdr(i,82) = dpdr(i,7)*p(12) + p(7)*dpdr(i,12)
     $      - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
      dpdr(i,83) = dpdr(i,8)*p(11) + p(8)*dpdr(i,11)
     $      - dpdr(i,72)
      dpdr(i,84) = dpdr(i,6)*p(17) + p(6)*dpdr(i,17)
      dpdr(i,85) = dpdr(i,7)*p(17) + p(7)*dpdr(i,17)
      dpdr(i,86) = dpdr(i,9)*p(11) + p(9)*dpdr(i,11)
     $      - dpdr(i,76) - dpdr(i,77)
      dpdr(i,87) = dpdr(i,2)*p(29) + p(2)*dpdr(i,29)
     $      - dpdr(i,81)
      dpdr(i,88) = dpdr(i,6)*p(14) + p(6)*dpdr(i,14)
     $      - dpdr(i,75) - dpdr(i,75)
      dpdr(i,89) = dpdr(i,2)*p(30) + p(2)*dpdr(i,30)
     $      - dpdr(i,82)
      dpdr(i,90) = dpdr(i,7)*p(14) + p(7)*dpdr(i,14)
     $      - dpdr(i,74) - dpdr(i,74)
      dpdr(i,91) = dpdr(i,1)*p(48) + p(1)*dpdr(i,48)
     $      - dpdr(i,76) - dpdr(i,81) - dpdr(i,82)
      dpdr(i,92) = dpdr(i,10)*p(11) + p(10)*dpdr(i,11)
     $      - dpdr(i,78)
      dpdr(i,93) = dpdr(i,2)*p(33) + p(2)*dpdr(i,33)
     $      - dpdr(i,88)
      dpdr(i,94) = dpdr(i,2)*p(34) + p(2)*dpdr(i,34)
     $      - dpdr(i,90)
      dpdr(i,95) = dpdr(i,10)*p(12) + p(10)*dpdr(i,12)
     $      - dpdr(i,77) - dpdr(i,93) - dpdr(i,94)
      dpdr(i,96) = dpdr(i,2)*p(27) + p(2)*dpdr(i,27)
     $      - dpdr(i,73) - dpdr(i,79)
      dpdr(i,97) = dpdr(i,2)*p(28) + p(2)*dpdr(i,28)
     $      - dpdr(i,76) - dpdr(i,86)
      dpdr(i,98) = dpdr(i,1)*p(53) + p(1)*dpdr(i,53)
     $      - dpdr(i,80) - dpdr(i,79) - dpdr(i,97)
      dpdr(i,99) = dpdr(i,1)*p(54) + p(1)*dpdr(i,54)
     $      - dpdr(i,81)
      dpdr(i,100) = dpdr(i,1)*p(55) + p(1)*dpdr(i,55)
     $      - dpdr(i,82)
      dpdr(i,101) = dpdr(i,5)*p(18) + p(5)*dpdr(i,18)
     $      - dpdr(i,76) - dpdr(i,91) - dpdr(i,87) - dpdr(i,89)
     $      - dpdr(i,86)
      dpdr(i,102) = dpdr(i,6)*p(18) + p(6)*dpdr(i,18)
     $      - dpdr(i,74) - dpdr(i,90)
      dpdr(i,103) = dpdr(i,7)*p(18) + p(7)*dpdr(i,18)
     $      - dpdr(i,75) - dpdr(i,88)
      dpdr(i,104) = dpdr(i,2)*p(31) + p(2)*dpdr(i,31)
     $      - dpdr(i,77) - dpdr(i,86)
      dpdr(i,105) = dpdr(i,2)*p(36) + p(2)*dpdr(i,36)
     $      - dpdr(i,91) - dpdr(i,101)
      dpdr(i,106) = dpdr(i,3)*p(32) + p(3)*dpdr(i,32)
     $      - dpdr(i,77) - dpdr(i,95) - dpdr(i,88) - dpdr(i,90)
     $      - dpdr(i,86)
      dpdr(i,107) = dpdr(i,4)*p(29) + p(4)*dpdr(i,29)
     $      - dpdr(i,82) - dpdr(i,80) - dpdr(i,99)
      dpdr(i,108) = dpdr(i,4)*p(30) + p(4)*dpdr(i,30)
     $      - dpdr(i,81) - dpdr(i,80) - dpdr(i,100)
      dpdr(i,109) = dpdr(i,2)*p(38) + p(2)*dpdr(i,38)
     $      - dpdr(i,95) - dpdr(i,106)
      dpdr(i,110) = dpdr(i,1)*p(62) + p(1)*dpdr(i,62)
     $      - dpdr(i,95) - dpdr(i,92) - dpdr(i,106)
      dpdr(i,111) = dpdr(i,6)*p(21) + p(6)*dpdr(i,21)
     $      - dpdr(i,94)
      dpdr(i,112) = dpdr(i,7)*p(21) + p(7)*dpdr(i,21)
     $      - dpdr(i,93)
      dpdr(i,113) = dpdr(i,1)*p(65) + p(1)*dpdr(i,65)
     $      - dpdr(i,96)
      dpdr(i,114) = dpdr(i,1)*p(66) + p(1)*dpdr(i,66)
     $      - dpdr(i,102) - dpdr(i,103) - dpdr(i,101)
      dpdr(i,115) = dpdr(i,2)*p(37) + p(2)*dpdr(i,37)
     $      - dpdr(i,92)
      dpdr(i,116) = dpdr(i,10)*p(18) + p(10)*dpdr(i,18)
     $      - dpdr(i,99) - dpdr(i,100) - dpdr(i,97)
      dpdr(i,117) = dpdr(i,2)*p(39) + p(2)*dpdr(i,39)
     $      - dpdr(i,110)
      dpdr(i,118) = dpdr(i,3)*p(39) + p(3)*dpdr(i,39)
     $      - dpdr(i,111) - dpdr(i,112) - dpdr(i,109)
      dpdr(i,119) = dpdr(i,1)*p(71) + p(1)*dpdr(i,71)
     $      - dpdr(i,118) - dpdr(i,117)
      dpdr(i,120) = dpdr(i,40)*p(2) + p(40)*dpdr(i,2)
      dpdr(i,121) = dpdr(i,40)*p(3) + p(40)*dpdr(i,3)
      dpdr(i,122) = dpdr(i,40)*p(4) + p(40)*dpdr(i,4)
      dpdr(i,123) = dpdr(i,11)*p(12) + p(11)*dpdr(i,12)
     $      - dpdr(i,121)
      dpdr(i,124) = dpdr(i,2)*p(42) + p(2)*dpdr(i,42)
     $      - dpdr(i,122)
      dpdr(i,125) = dpdr(i,6)*p(22) + p(6)*dpdr(i,22)
     $      - dpdr(i,121)
      dpdr(i,126) = dpdr(i,7)*p(22) + p(7)*dpdr(i,22)
     $      - dpdr(i,121)
      dpdr(i,127) = dpdr(i,2)*p(41) + p(2)*dpdr(i,41)
     $      - dpdr(i,121) - dpdr(i,123) - dpdr(i,121)
      dpdr(i,128) = dpdr(i,6)*p(23) + p(6)*dpdr(i,23)
     $      - dpdr(i,123)
      dpdr(i,129) = dpdr(i,7)*p(23) + p(7)*dpdr(i,23)
     $      - dpdr(i,123)
      dpdr(i,130) = dpdr(i,3)*p(41) + p(3)*dpdr(i,41)
     $      - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,125)
     $      - dpdr(i,128) - dpdr(i,126) - dpdr(i,129) - dpdr(i,124)
     $      - dpdr(i,123) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120)
     $      - dpdr(i,125) - dpdr(i,126) - dpdr(i,124) - dpdr(i,123)
     $      - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,122)
     $      - dpdr(i,121) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120)
     $      - dpdr(i,120) - dpdr(i,120)
      dpdr(i,131) = dpdr(i,3)*p(42) + p(3)*dpdr(i,42)
     $      - dpdr(i,121) - dpdr(i,125) - dpdr(i,126) - dpdr(i,121)
      dpdr(i,132) = dpdr(i,4)*p(41) + p(4)*dpdr(i,41)
     $      - dpdr(i,121) - dpdr(i,131) - dpdr(i,128) - dpdr(i,129)
     $      - dpdr(i,127) - dpdr(i,121)
      dpdr(i,133) = dpdr(i,1)*p(78) + p(1)*dpdr(i,78)
     $      - dpdr(i,122) - dpdr(i,131)
      dpdr(i,134) = dpdr(i,11)*p(11) + p(11)*dpdr(i,11)
     $      - dpdr(i,120) - dpdr(i,120)
      dpdr(i,135) = dpdr(i,2)*p(48) + p(2)*dpdr(i,48)
     $      - dpdr(i,130)
      dpdr(i,136) = dpdr(i,11)*p(17) + p(11)*dpdr(i,17)
     $      - dpdr(i,120)
      dpdr(i,137) = dpdr(i,2)*p(44) + p(2)*dpdr(i,44)
     $      - dpdr(i,123)
      dpdr(i,138) = dpdr(i,2)*p(45) + p(2)*dpdr(i,45)
     $      - dpdr(i,121)
      dpdr(i,139) = dpdr(i,11)*p(14) + p(11)*dpdr(i,14)
     $      - dpdr(i,122) - dpdr(i,124) - dpdr(i,122)
      dpdr(i,140) = dpdr(i,6)*p(27) + p(6)*dpdr(i,27)
     $      - dpdr(i,127)
      dpdr(i,141) = dpdr(i,2)*p(54) + p(2)*dpdr(i,54)
      dpdr(i,142) = dpdr(i,7)*p(27) + p(7)*dpdr(i,27)
     $      - dpdr(i,127)
      dpdr(i,143) = dpdr(i,2)*p(52) + p(2)*dpdr(i,52)
     $      - dpdr(i,131) - dpdr(i,132)
      dpdr(i,144) = dpdr(i,1)*p(81) + p(1)*dpdr(i,81)
     $      - dpdr(i,125) - dpdr(i,141) - dpdr(i,135) - dpdr(i,125)
      dpdr(i,145) = dpdr(i,2)*p(55) + p(2)*dpdr(i,55)
      dpdr(i,146) = dpdr(i,1)*p(82) + p(1)*dpdr(i,82)
     $      - dpdr(i,126) - dpdr(i,135) - dpdr(i,145) - dpdr(i,126)
      dpdr(i,147) = dpdr(i,11)*p(18) + p(11)*dpdr(i,18)
     $      - dpdr(i,130)
      dpdr(i,148) = dpdr(i,6)*p(28) + p(6)*dpdr(i,28)
     $      - dpdr(i,129) - dpdr(i,123) - dpdr(i,143)
      dpdr(i,149) = dpdr(i,7)*p(28) + p(7)*dpdr(i,28)
     $      - dpdr(i,128) - dpdr(i,123) - dpdr(i,143)
      dpdr(i,150) = dpdr(i,6)*p(30) + p(6)*dpdr(i,30)
     $      - dpdr(i,121) - dpdr(i,146)
      dpdr(i,151) = dpdr(i,11)*p(19) + p(11)*dpdr(i,19)
     $      - dpdr(i,122)
      dpdr(i,152) = dpdr(i,2)*p(50) + p(2)*dpdr(i,50)
     $      - dpdr(i,128)
      dpdr(i,153) = dpdr(i,2)*p(51) + p(2)*dpdr(i,51)
     $      - dpdr(i,129)
      dpdr(i,154) = dpdr(i,11)*p(20) + p(11)*dpdr(i,20)
     $      - dpdr(i,131) - dpdr(i,132)
      dpdr(i,155) = dpdr(i,2)*p(59) + p(2)*dpdr(i,59)
     $      - dpdr(i,144) - dpdr(i,148)
      dpdr(i,156) = dpdr(i,6)*p(32) + p(6)*dpdr(i,32)
     $      - dpdr(i,129)
      dpdr(i,157) = dpdr(i,2)*p(60) + p(2)*dpdr(i,60)
     $      - dpdr(i,146) - dpdr(i,149)
      dpdr(i,158) = dpdr(i,7)*p(32) + p(7)*dpdr(i,32)
     $      - dpdr(i,128)
      dpdr(i,159) = dpdr(i,6)*p(34) + p(6)*dpdr(i,34)
     $      - dpdr(i,122) - dpdr(i,145) - dpdr(i,122)
      dpdr(i,160) = dpdr(i,11)*p(21) + p(11)*dpdr(i,21)
     $      - dpdr(i,133)
      dpdr(i,161) = dpdr(i,2)*p(63) + p(2)*dpdr(i,63)
     $      - dpdr(i,156)
      dpdr(i,162) = dpdr(i,2)*p(64) + p(2)*dpdr(i,64)
     $      - dpdr(i,158)
      dpdr(i,163) = dpdr(i,12)*p(21) + p(12)*dpdr(i,21)
     $      - dpdr(i,132) - dpdr(i,161) - dpdr(i,162)
      dpdr(i,164) = dpdr(i,2)*p(53) + p(2)*dpdr(i,53)
     $      - dpdr(i,124) - dpdr(i,139)
      dpdr(i,165) = dpdr(i,2)*p(58) + p(2)*dpdr(i,58)
     $      - dpdr(i,131) - dpdr(i,154)
      dpdr(i,166) = dpdr(i,3)*p(54) + p(3)*dpdr(i,54)
     $      - dpdr(i,125) - dpdr(i,144)
      dpdr(i,167) = dpdr(i,3)*p(55) + p(3)*dpdr(i,55)
     $      - dpdr(i,126) - dpdr(i,146)
      dpdr(i,168) = dpdr(i,3)*p(65) + p(3)*dpdr(i,65)
     $      - dpdr(i,137)
      dpdr(i,169) = dpdr(i,17)*p(18) + p(17)*dpdr(i,18)
     $      - dpdr(i,135)
      dpdr(i,170) = dpdr(i,1)*p(98) + p(1)*dpdr(i,98)
     $      - dpdr(i,143) - dpdr(i,139) - dpdr(i,165)
      dpdr(i,171) = dpdr(i,4)*p(54) + p(4)*dpdr(i,54)
     $      - dpdr(i,135)
      dpdr(i,172) = dpdr(i,4)*p(55) + p(4)*dpdr(i,55)
     $      - dpdr(i,135)
      dpdr(i,173) = dpdr(i,2)*p(66) + p(2)*dpdr(i,66)
     $      - dpdr(i,150)
      dpdr(i,174) = dpdr(i,1)*p(101) + p(1)*dpdr(i,101)
     $      - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,173)
     $      - dpdr(i,165) - dpdr(i,147)
      dpdr(i,175) = dpdr(i,1)*p(102) + p(1)*dpdr(i,102)
     $      - dpdr(i,150) - dpdr(i,148) - dpdr(i,166)
      dpdr(i,176) = dpdr(i,1)*p(103) + p(1)*dpdr(i,103)
     $      - dpdr(i,150) - dpdr(i,149) - dpdr(i,167)
      dpdr(i,177) = dpdr(i,2)*p(61) + p(2)*dpdr(i,61)
     $      - dpdr(i,132) - dpdr(i,154)
      dpdr(i,178) = dpdr(i,2)*p(68) + p(2)*dpdr(i,68)
     $      - dpdr(i,159) - dpdr(i,174)
      dpdr(i,179) = dpdr(i,3)*p(62) + p(3)*dpdr(i,62)
     $      - dpdr(i,132) - dpdr(i,163) - dpdr(i,156) - dpdr(i,158)
     $      - dpdr(i,154)
      dpdr(i,180) = dpdr(i,6)*p(38) + p(6)*dpdr(i,38)
     $      - dpdr(i,132) - dpdr(i,163) - dpdr(i,157)
      dpdr(i,181) = dpdr(i,7)*p(38) + p(7)*dpdr(i,38)
     $      - dpdr(i,132) - dpdr(i,163) - dpdr(i,155)
      dpdr(i,182) = dpdr(i,2)*p(70) + p(2)*dpdr(i,70)
     $      - dpdr(i,163) - dpdr(i,179)
      dpdr(i,183) = dpdr(i,1)*p(110) + p(1)*dpdr(i,110)
     $      - dpdr(i,163) - dpdr(i,160) - dpdr(i,179)
      dpdr(i,184) = dpdr(i,6)*p(39) + p(6)*dpdr(i,39)
     $      - dpdr(i,162)
      dpdr(i,185) = dpdr(i,7)*p(39) + p(7)*dpdr(i,39)
     $      - dpdr(i,161)
      dpdr(i,186) = dpdr(i,2)*p(65) + p(2)*dpdr(i,65)
     $      - dpdr(i,136)
      dpdr(i,187) = dpdr(i,3)*p(66) + p(3)*dpdr(i,66)
     $      - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,175)
     $      - dpdr(i,176) - dpdr(i,174)
      dpdr(i,188) = dpdr(i,2)*p(67) + p(2)*dpdr(i,67)
     $      - dpdr(i,151)
      dpdr(i,189) = dpdr(i,4)*p(66) + p(4)*dpdr(i,66)
     $      - dpdr(i,166) - dpdr(i,167) - dpdr(i,165)
      dpdr(i,190) = dpdr(i,2)*p(69) + p(2)*dpdr(i,69)
     $      - dpdr(i,160)
      dpdr(i,191) = dpdr(i,18)*p(21) + p(18)*dpdr(i,21)
     $      - dpdr(i,171) - dpdr(i,172) - dpdr(i,169)
      dpdr(i,192) = dpdr(i,2)*p(71) + p(2)*dpdr(i,71)
     $      - dpdr(i,183)
      dpdr(i,193) = dpdr(i,3)*p(71) + p(3)*dpdr(i,71)
     $      - dpdr(i,184) - dpdr(i,185) - dpdr(i,182)
      dpdr(i,194) = dpdr(i,1)*p(119) + p(1)*dpdr(i,119)
     $      - dpdr(i,193) - dpdr(i,192)
      dpdr(i,195) = dpdr(i,40)*p(5) + p(40)*dpdr(i,5)
      dpdr(i,196) = dpdr(i,40)*p(6) + p(40)*dpdr(i,6)
      dpdr(i,197) = dpdr(i,40)*p(7) + p(40)*dpdr(i,7)
      dpdr(i,198) = dpdr(i,40)*p(8) + p(40)*dpdr(i,8)
      dpdr(i,199) = dpdr(i,40)*p(9) + p(40)*dpdr(i,9)
      dpdr(i,200) = dpdr(i,40)*p(10) + p(40)*dpdr(i,10)
      dpdr(i,201) = dpdr(i,11)*p(22) + p(11)*dpdr(i,22)
     $      - dpdr(i,195)
      dpdr(i,202) = dpdr(i,12)*p(22) + p(12)*dpdr(i,22)
     $      - dpdr(i,196) - dpdr(i,197) - dpdr(i,195) - dpdr(i,196)
     $      - dpdr(i,197) - dpdr(i,195) - dpdr(i,196) - dpdr(i,197)
      dpdr(i,203) = dpdr(i,17)*p(22) + p(17)*dpdr(i,22)
     $      - dpdr(i,198)
      dpdr(i,204) = dpdr(i,6)*p(43) + p(6)*dpdr(i,43)
      dpdr(i,205) = dpdr(i,7)*p(43) + p(7)*dpdr(i,43)
      dpdr(i,206) = dpdr(i,11)*p(26) + p(11)*dpdr(i,26)
     $      - dpdr(i,199) - dpdr(i,202)
      dpdr(i,207) = dpdr(i,2)*p(76) + p(2)*dpdr(i,76)
     $      - dpdr(i,199) - dpdr(i,202)
      dpdr(i,208) = dpdr(i,2)*p(78) + p(2)*dpdr(i,78)
     $      - dpdr(i,200)
      dpdr(i,209) = dpdr(i,6)*p(41) + p(6)*dpdr(i,41)
     $      - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
      dpdr(i,210) = dpdr(i,6)*p(42) + p(6)*dpdr(i,42)
     $      - dpdr(i,197) - dpdr(i,197) - dpdr(i,197)
      dpdr(i,211) = dpdr(i,7)*p(41) + p(7)*dpdr(i,41)
     $      - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
      dpdr(i,212) = dpdr(i,7)*p(42) + p(7)*dpdr(i,42)
     $      - dpdr(i,196) - dpdr(i,196) - dpdr(i,196)
      dpdr(i,213) = dpdr(i,11)*p(29) + p(11)*dpdr(i,29)
     $      - dpdr(i,209)
      dpdr(i,214) = dpdr(i,11)*p(30) + p(11)*dpdr(i,30)
     $      - dpdr(i,211)
      dpdr(i,215) = dpdr(i,18)*p(22) + p(18)*dpdr(i,22)
     $      - dpdr(i,199) - dpdr(i,213) - dpdr(i,214)
      dpdr(i,216) = dpdr(i,2)*p(77) + p(2)*dpdr(i,77)
     $      - dpdr(i,199) - dpdr(i,206)
      dpdr(i,217) = dpdr(i,6)*p(49) + p(6)*dpdr(i,49)
     $      - dpdr(i,205)
      dpdr(i,218) = dpdr(i,7)*p(49) + p(7)*dpdr(i,49)
     $      - dpdr(i,204)
      dpdr(i,219) = dpdr(i,3)*p(77) + p(3)*dpdr(i,77)
     $      - dpdr(i,200) - dpdr(i,199) - dpdr(i,198) - dpdr(i,209)
     $      - dpdr(i,217) - dpdr(i,211) - dpdr(i,218) - dpdr(i,207)
     $      - dpdr(i,204) - dpdr(i,205) - dpdr(i,200) - dpdr(i,199)
     $      - dpdr(i,198) - dpdr(i,200) - dpdr(i,198) - dpdr(i,200)
     $      - dpdr(i,198)
      dpdr(i,220) = dpdr(i,3)*p(78) + p(3)*dpdr(i,78)
     $      - dpdr(i,199) - dpdr(i,210) - dpdr(i,212)
      dpdr(i,221) = dpdr(i,10)*p(41) + p(10)*dpdr(i,41)
     $      - dpdr(i,199) - dpdr(i,220) - dpdr(i,217) - dpdr(i,218)
     $      - dpdr(i,216)
      dpdr(i,222) = dpdr(i,1)*p(133) + p(1)*dpdr(i,133)
     $      - dpdr(i,200) - dpdr(i,220)
      dpdr(i,223) = dpdr(i,11)*p(27) + p(11)*dpdr(i,27)
     $      - dpdr(i,195) - dpdr(i,203)
      dpdr(i,224) = dpdr(i,1)*p(134) + p(1)*dpdr(i,134)
     $      - dpdr(i,201)
      dpdr(i,225) = dpdr(i,2)*p(81) + p(2)*dpdr(i,81)
     $      - dpdr(i,209)
      dpdr(i,226) = dpdr(i,2)*p(80) + p(2)*dpdr(i,80)
     $      - dpdr(i,202) - dpdr(i,206)
      dpdr(i,227) = dpdr(i,2)*p(82) + p(2)*dpdr(i,82)
     $      - dpdr(i,211)
      dpdr(i,228) = dpdr(i,2)*p(91) + p(2)*dpdr(i,91)
     $      - dpdr(i,215) - dpdr(i,219)
      dpdr(i,229) = dpdr(i,11)*p(28) + p(11)*dpdr(i,28)
     $      - dpdr(i,199) - dpdr(i,207)
      dpdr(i,230) = dpdr(i,6)*p(53) + p(6)*dpdr(i,53)
     $      - dpdr(i,205)
      dpdr(i,231) = dpdr(i,5)*p(54) + p(5)*dpdr(i,54)
     $      - dpdr(i,209)
      dpdr(i,232) = dpdr(i,7)*p(53) + p(7)*dpdr(i,53)
     $      - dpdr(i,204)
      dpdr(i,233) = dpdr(i,7)*p(54) + p(7)*dpdr(i,54)
     $      - dpdr(i,196)
      dpdr(i,234) = dpdr(i,5)*p(55) + p(5)*dpdr(i,55)
     $      - dpdr(i,211)
      dpdr(i,235) = dpdr(i,6)*p(55) + p(6)*dpdr(i,55)
     $      - dpdr(i,197)
      dpdr(i,236) = dpdr(i,11)*p(35) + p(11)*dpdr(i,35)
     $      - dpdr(i,198)
      dpdr(i,237) = dpdr(i,6)*p(65) + p(6)*dpdr(i,65)
      dpdr(i,238) = dpdr(i,7)*p(65) + p(7)*dpdr(i,65)
      dpdr(i,239) = dpdr(i,2)*p(86) + p(2)*dpdr(i,86)
     $      - dpdr(i,199) - dpdr(i,229)
      dpdr(i,240) = dpdr(i,1)*p(139) + p(1)*dpdr(i,139)
     $      - dpdr(i,206) - dpdr(i,229) - dpdr(i,224)
      dpdr(i,241) = dpdr(i,17)*p(29) + p(17)*dpdr(i,29)
     $      - dpdr(i,225)
      dpdr(i,242) = dpdr(i,2)*p(99) + p(2)*dpdr(i,99)
     $      - dpdr(i,231)
      dpdr(i,243) = dpdr(i,17)*p(30) + p(17)*dpdr(i,30)
     $      - dpdr(i,227)
      dpdr(i,244) = dpdr(i,2)*p(95) + p(2)*dpdr(i,95)
     $      - dpdr(i,220) - dpdr(i,221)
      dpdr(i,245) = dpdr(i,4)*p(81) + p(4)*dpdr(i,81)
     $      - dpdr(i,202) - dpdr(i,242) - dpdr(i,227)
      dpdr(i,246) = dpdr(i,2)*p(100) + p(2)*dpdr(i,100)
     $      - dpdr(i,234)
      dpdr(i,247) = dpdr(i,4)*p(82) + p(4)*dpdr(i,82)
     $      - dpdr(i,202) - dpdr(i,225) - dpdr(i,246)
      dpdr(i,248) = dpdr(i,11)*p(36) + p(11)*dpdr(i,36)
     $      - dpdr(i,215) - dpdr(i,219)
      dpdr(i,249) = dpdr(i,2)*p(102) + p(2)*dpdr(i,102)
     $      - dpdr(i,233)
      dpdr(i,250) = dpdr(i,6)*p(58) + p(6)*dpdr(i,58)
     $      - dpdr(i,206) - dpdr(i,214) - dpdr(i,244) - dpdr(i,214)
      dpdr(i,251) = dpdr(i,2)*p(103) + p(2)*dpdr(i,103)
     $      - dpdr(i,235)
      dpdr(i,252) = dpdr(i,7)*p(58) + p(7)*dpdr(i,58)
     $      - dpdr(i,213) - dpdr(i,206) - dpdr(i,244) - dpdr(i,213)
      dpdr(i,253) = dpdr(i,1)*p(150) + p(1)*dpdr(i,150)
     $      - dpdr(i,215) - dpdr(i,233) - dpdr(i,235)
      dpdr(i,254) = dpdr(i,11)*p(37) + p(11)*dpdr(i,37)
     $      - dpdr(i,200)
      dpdr(i,255) = dpdr(i,2)*p(93) + p(2)*dpdr(i,93)
     $      - dpdr(i,217)
      dpdr(i,256) = dpdr(i,2)*p(94) + p(2)*dpdr(i,94)
     $      - dpdr(i,218)
      dpdr(i,257) = dpdr(i,11)*p(38) + p(11)*dpdr(i,38)
     $      - dpdr(i,220) - dpdr(i,221)
      dpdr(i,258) = dpdr(i,2)*p(107) + p(2)*dpdr(i,107)
     $      - dpdr(i,245) - dpdr(i,250)
      dpdr(i,259) = dpdr(i,6)*p(62) + p(6)*dpdr(i,62)
     $      - dpdr(i,218)
      dpdr(i,260) = dpdr(i,2)*p(108) + p(2)*dpdr(i,108)
     $      - dpdr(i,247) - dpdr(i,252)
      dpdr(i,261) = dpdr(i,7)*p(62) + p(7)*dpdr(i,62)
     $      - dpdr(i,217)
      dpdr(i,262) = dpdr(i,6)*p(64) + p(6)*dpdr(i,64)
     $      - dpdr(i,200) - dpdr(i,246) - dpdr(i,200)
      dpdr(i,263) = dpdr(i,11)*p(39) + p(11)*dpdr(i,39)
     $      - dpdr(i,222)
      dpdr(i,264) = dpdr(i,2)*p(111) + p(2)*dpdr(i,111)
     $      - dpdr(i,259)
      dpdr(i,265) = dpdr(i,2)*p(112) + p(2)*dpdr(i,112)
     $      - dpdr(i,261)
      dpdr(i,266) = dpdr(i,12)*p(39) + p(12)*dpdr(i,39)
     $      - dpdr(i,221) - dpdr(i,264) - dpdr(i,265)
      dpdr(i,267) = dpdr(i,2)*p(98) + p(2)*dpdr(i,98)
     $      - dpdr(i,208) - dpdr(i,240)
      dpdr(i,268) = dpdr(i,6)*p(54) + p(6)*dpdr(i,54)
     $      - dpdr(i,210)
      dpdr(i,269) = dpdr(i,7)*p(55) + p(7)*dpdr(i,55)
     $      - dpdr(i,212)
      dpdr(i,270) = dpdr(i,2)*p(96) + p(2)*dpdr(i,96)
     $      - dpdr(i,203) - dpdr(i,223)
      dpdr(i,271) = dpdr(i,2)*p(97) + p(2)*dpdr(i,97)
     $      - dpdr(i,207) - dpdr(i,229)
      dpdr(i,272) = dpdr(i,2)*p(101) + p(2)*dpdr(i,101)
     $      - dpdr(i,215) - dpdr(i,248)
      dpdr(i,273) = dpdr(i,2)*p(106) + p(2)*dpdr(i,106)
     $      - dpdr(i,220) - dpdr(i,257)
      dpdr(i,274) = dpdr(i,6)*p(59) + p(6)*dpdr(i,59)
     $      - dpdr(i,220) - dpdr(i,215) - dpdr(i,209)
      dpdr(i,275) = dpdr(i,7)*p(60) + p(7)*dpdr(i,60)
     $      - dpdr(i,220) - dpdr(i,215) - dpdr(i,211)
      dpdr(i,276) = dpdr(i,5)*p(66) + p(5)*dpdr(i,66)
     $      - dpdr(i,215) - dpdr(i,253) - dpdr(i,249) - dpdr(i,251)
     $      - dpdr(i,248)
      dpdr(i,277) = dpdr(i,6)*p(66) + p(6)*dpdr(i,66)
     $      - dpdr(i,213) - dpdr(i,252)
      dpdr(i,278) = dpdr(i,7)*p(66) + p(7)*dpdr(i,66)
     $      - dpdr(i,214) - dpdr(i,250)
      dpdr(i,279) = dpdr(i,2)*p(104) + p(2)*dpdr(i,104)
     $      - dpdr(i,216) - dpdr(i,239)
      dpdr(i,280) = dpdr(i,2)*p(105) + p(2)*dpdr(i,105)
     $      - dpdr(i,219) - dpdr(i,248)
      dpdr(i,281) = dpdr(i,1)*p(170) + p(1)*dpdr(i,170)
     $      - dpdr(i,244) - dpdr(i,240) - dpdr(i,273)
      dpdr(i,282) = dpdr(i,10)*p(54) + p(10)*dpdr(i,54)
     $      - dpdr(i,227)
      dpdr(i,283) = dpdr(i,10)*p(55) + p(10)*dpdr(i,55)
     $      - dpdr(i,225)
      dpdr(i,284) = dpdr(i,2)*p(114) + p(2)*dpdr(i,114)
     $      - dpdr(i,253) - dpdr(i,276)
      dpdr(i,285) = dpdr(i,1)*p(174) + p(1)*dpdr(i,174)
     $      - dpdr(i,250) - dpdr(i,252) - dpdr(i,248) - dpdr(i,276)
     $      - dpdr(i,273)
      dpdr(i,286) = dpdr(i,6)*p(68) + p(6)*dpdr(i,68)
     $      - dpdr(i,217) - dpdr(i,261) - dpdr(i,251)
      dpdr(i,287) = dpdr(i,7)*p(68) + p(7)*dpdr(i,68)
     $      - dpdr(i,218) - dpdr(i,259) - dpdr(i,249)
      dpdr(i,288) = dpdr(i,2)*p(109) + p(2)*dpdr(i,109)
     $      - dpdr(i,221) - dpdr(i,257)
      dpdr(i,289) = dpdr(i,2)*p(116) + p(2)*dpdr(i,116)
     $      - dpdr(i,262) - dpdr(i,285)
      dpdr(i,290) = dpdr(i,3)*p(110) + p(3)*dpdr(i,110)
     $      - dpdr(i,221) - dpdr(i,266) - dpdr(i,259) - dpdr(i,261)
     $      - dpdr(i,257)
      dpdr(i,291) = dpdr(i,6)*p(70) + p(6)*dpdr(i,70)
     $      - dpdr(i,221) - dpdr(i,266) - dpdr(i,260)
      dpdr(i,292) = dpdr(i,7)*p(70) + p(7)*dpdr(i,70)
     $      - dpdr(i,221) - dpdr(i,266) - dpdr(i,258)
      dpdr(i,293) = dpdr(i,2)*p(118) + p(2)*dpdr(i,118)
     $      - dpdr(i,266) - dpdr(i,290)
      dpdr(i,294) = dpdr(i,1)*p(183) + p(1)*dpdr(i,183)
     $      - dpdr(i,266) - dpdr(i,263) - dpdr(i,290)
      dpdr(i,295) = dpdr(i,6)*p(71) + p(6)*dpdr(i,71)
     $      - dpdr(i,265)
      dpdr(i,296) = dpdr(i,7)*p(71) + p(7)*dpdr(i,71)
     $      - dpdr(i,264)
      dpdr(i,297) = dpdr(i,1)*p(186) + p(1)*dpdr(i,186)
     $      - dpdr(i,270)
      dpdr(i,298) = dpdr(i,1)*p(187) + p(1)*dpdr(i,187)
     $      - dpdr(i,277) - dpdr(i,278) - dpdr(i,276)
      dpdr(i,299) = dpdr(i,2)*p(115) + p(2)*dpdr(i,115)
     $      - dpdr(i,254)
      dpdr(i,300) = dpdr(i,1)*p(189) + p(1)*dpdr(i,189)
     $      - dpdr(i,286) - dpdr(i,287) - dpdr(i,285) - dpdr(i,284)
     $      - dpdr(i,298)
      dpdr(i,301) = dpdr(i,2)*p(117) + p(2)*dpdr(i,117)
     $      - dpdr(i,263)
      dpdr(i,302) = dpdr(i,18)*p(39) + p(18)*dpdr(i,39)
     $      - dpdr(i,282) - dpdr(i,283) - dpdr(i,280)
      dpdr(i,303) = dpdr(i,2)*p(119) + p(2)*dpdr(i,119)
     $      - dpdr(i,294)
      dpdr(i,304) = dpdr(i,3)*p(119) + p(3)*dpdr(i,119)
     $      - dpdr(i,295) - dpdr(i,296) - dpdr(i,293)
      dpdr(i,305) = dpdr(i,1)*p(194) + p(1)*dpdr(i,194)
     $      - dpdr(i,304) - dpdr(i,303)
      dpdr(i,306) = dpdr(i,40)*p(11) + p(40)*dpdr(i,11)  
      dpdr(i,307) = dpdr(i,40)*p(12) + p(40)*dpdr(i,12)  
      dpdr(i,308) = dpdr(i,40)*p(17) + p(40)*dpdr(i,17)  
      dpdr(i,309) = dpdr(i,40)*p(13) + p(40)*dpdr(i,13)  
      dpdr(i,310) = dpdr(i,40)*p(14) + p(40)*dpdr(i,14)  
      dpdr(i,311) = dpdr(i,40)*p(15) + p(40)*dpdr(i,15)  
      dpdr(i,312) = dpdr(i,40)*p(16) + p(40)*dpdr(i,16)  
      dpdr(i,313) = dpdr(i,40)*p(18) + p(40)*dpdr(i,18)  
      dpdr(i,314) = dpdr(i,40)*p(19) + p(40)*dpdr(i,19)  
      dpdr(i,315) = dpdr(i,40)*p(20) + p(40)*dpdr(i,20)  
      dpdr(i,316) = dpdr(i,40)*p(21) + p(40)*dpdr(i,21)  
      dpdr(i,317) = dpdr(i,11)*p(42) + p(11)*dpdr(i,42) - dpdr(i,310)
      dpdr(i,318) = dpdr(i,11)*p(44) + p(11)*dpdr(i,44) - dpdr(i,307)
      dpdr(i,319) = dpdr(i,12)*p(43) + p(12)*dpdr(i,43) - dpdr(i,309)
     $ - dpdr(i,318)
      dpdr(i,320) = dpdr(i,17)*p(42) + p(17)*dpdr(i,42) - dpdr(i,314)
      dpdr(i,321) = dpdr(i,2)*p(125) + p(2)*dpdr(i,125) - dpdr(i,311)
      dpdr(i,322) = dpdr(i,2)*p(126) + p(2)*dpdr(i,126) - dpdr(i,312)
      dpdr(i,323) = dpdr(i,11)*p(48) + p(11)*dpdr(i,48) - dpdr(i,313)
      dpdr(i,324) = dpdr(i,12)*p(42) + p(12)*dpdr(i,42) - dpdr(i,311)
     $ - dpdr(i,312) - dpdr(i,307) - dpdr(i,307)
      dpdr(i,325) = dpdr(i,6)*p(79)  + p(6)*dpdr(i,79) - dpdr(i,319)
      dpdr(i,326) = dpdr(i,11)*p(54) + p(11)*dpdr(i,54)  
      dpdr(i,327) = dpdr(i,7)*p(79)  + p(7)*dpdr(i,79) - dpdr(i,319)
      dpdr(i,328) = dpdr(i,2)*p(131) + p(2)*dpdr(i,131) - dpdr(i,315)
     $ - dpdr(i,324)
      dpdr(i,329) = dpdr(i,22)*p(29) + p(22)*dpdr(i,29) - dpdr(i,313)
     $ - dpdr(i,311) - dpdr(i,326) - dpdr(i,311)
      dpdr(i,330) = dpdr(i,11)*p(55) + p(11)*dpdr(i,55)  
      dpdr(i,331) = dpdr(i,22)*p(30) + p(22)*dpdr(i,30) - dpdr(i,313)
     $ - dpdr(i,312) - dpdr(i,330) - dpdr(i,312)
      dpdr(i,332) = dpdr(i,2)*p(127) + p(2)*dpdr(i,127) - dpdr(i,309)
     $ - dpdr(i,319)
      dpdr(i,333) = dpdr(i,6)*p(83)  + p(6)*dpdr(i,83) - dpdr(i,318)
      dpdr(i,334) = dpdr(i,7)*p(83)  + p(7)*dpdr(i,83) - dpdr(i,318)
      dpdr(i,335) = dpdr(i,11)*p(52) + p(11)*dpdr(i,52) - dpdr(i,315) 
     $ - dpdr(i,324)
      dpdr(i,336) = dpdr(i,2)*p(130) + p(2)*dpdr(i,130) - dpdr(i,313) 
     $ - dpdr(i,323) - dpdr(i,313)
      dpdr(i,337) = dpdr(i,2)*p(133) + p(2)*dpdr(i,133) - dpdr(i,316)
      dpdr(i,338) = dpdr(i,6)*p(77)  + p(6)*dpdr(i,77)  - dpdr(i,315)
     $ - dpdr(i,309) - dpdr(i,322)
      dpdr(i,339) = dpdr(i,6)*p(78)  + p(6)*dpdr(i,78) - dpdr(i,312)
      dpdr(i,340) = dpdr(i,7)*p(77)  + p(7)*dpdr(i,77) - dpdr(i,315)
     $ - dpdr(i,309) - dpdr(i,321)
      dpdr(i,341) = dpdr(i,7)*p(78)  + p(7)*dpdr(i,78) - dpdr(i,311)
      dpdr(i,342) = dpdr(i,11)*p(59) + p(11)*dpdr(i,59) - dpdr(i,329)
     $ - dpdr(i,338)
      dpdr(i,343) = dpdr(i,11)*p(60) + p(11)*dpdr(i,60) - dpdr(i,331)
     $ - dpdr(i,340)
      dpdr(i,344) = dpdr(i,3)*p(130) + p(3)*dpdr(i,130) - dpdr(i,315)
     $ - dpdr(i,311) - dpdr(i,312) - dpdr(i,309) - dpdr(i,329) 
     $ - dpdr(i,338) - dpdr(i,331) - dpdr(i,340) - dpdr(i,328) 
     $ - dpdr(i,321) - dpdr(i,322) - dpdr(i,311) - dpdr(i,312)
      dpdr(i,345) = dpdr(i,18)*p(42) + p(18)*dpdr(i,42) - dpdr(i,310)
     $ - dpdr(i,326) - dpdr(i,330) - dpdr(i,310)
      dpdr(i,346) = dpdr(i,2)*p(132) + p(2)*dpdr(i,132) - dpdr(i,315)
     $ - dpdr(i,335)
      dpdr(i,347) = dpdr(i,6)*p(92)  + p(6)*dpdr(i,92) - dpdr(i,334)
      dpdr(i,348) = dpdr(i,7)*p(92)  + p(7)*dpdr(i,92) - dpdr(i,333)
      dpdr(i,349) = dpdr(i,3)*p(132) + p(3)*dpdr(i,132) - dpdr(i,316)
     $ - dpdr(i,315) - dpdr(i,314) - dpdr(i,338) - dpdr(i,347)
     $ - dpdr(i,340) - dpdr(i,348) - dpdr(i,336) - dpdr(i,333)
     $ - dpdr(i,334) - dpdr(i,316) - dpdr(i,315) - dpdr(i,314)
     $ - dpdr(i,316) - dpdr(i,314) - dpdr(i,316) - dpdr(i,314)
      dpdr(i,350) = dpdr(i,3)*p(133) + p(3)*dpdr(i,133) - dpdr(i,315)
     $ - dpdr(i,339) - dpdr(i,341)
      dpdr(i,351) = dpdr(i,4)*p(132) + p(4)*dpdr(i,132) - dpdr(i,315) 
     $ - dpdr(i,344) - dpdr(i,342) - dpdr(i,343) - dpdr(i,332)
      dpdr(i,352) = dpdr(i,1)*p(222) + p(1)*dpdr(i,222) - dpdr(i,316)
     $ - dpdr(i,350)
      dpdr(i,353) = dpdr(i,2)*p(134) + p(2)*dpdr(i,134) - dpdr(i,306)
      dpdr(i,354) = dpdr(i,2)*p(135) + p(2)*dpdr(i,135) - dpdr(i,323)
      dpdr(i,355) = dpdr(i,3)*p(134) + p(3)*dpdr(i,134) - dpdr(i,319)
      dpdr(i,356) = dpdr(i,11)*p(53) + p(11)*dpdr(i,53) - dpdr(i,310)
     $ - dpdr(i,320)
      dpdr(i,357) = dpdr(i,2)*p(144) + p(2)*dpdr(i,144) - dpdr(i,329)
     $ - dpdr(i,338)
      dpdr(i,358) = dpdr(i,2)*p(143) + p(2)*dpdr(i,143) - dpdr(i,324)
     $ - dpdr(i,335)
      dpdr(i,359) = dpdr(i,6)*p(81)  + p(6)*dpdr(i,81) - dpdr(i,311)
     $ - dpdr(i,324)
      dpdr(i,360) = dpdr(i,2)*p(146) + p(2)*dpdr(i,146) - dpdr(i,331)
     $ - dpdr(i,340)
      dpdr(i,361) = dpdr(i,2)*p(150) + p(2)*dpdr(i,150) - dpdr(i,344)
      dpdr(i,362) = dpdr(i,7)*p(82)  + p(7)*dpdr(i,82) - dpdr(i,312)
     $ - dpdr(i,324)
      dpdr(i,363) = dpdr(i,11)*p(65) + p(11)*dpdr(i,65) - dpdr(i,308)
      dpdr(i,364) = dpdr(i,2)*p(137) + p(2)*dpdr(i,137) - dpdr(i,318)
      dpdr(i,365) = dpdr(i,17)*p(45) + p(17)*dpdr(i,45) - dpdr(i,307)
      dpdr(i,366) = dpdr(i,4)*p(134) + p(4)*dpdr(i,134) - dpdr(i,317)
      dpdr(i,367) = dpdr(i,6)*p(96)  + p(6)*dpdr(i,96) - dpdr(i,332)
      dpdr(i,368) = dpdr(i,17)*p(54) + p(17)*dpdr(i,54)  
      dpdr(i,369) = dpdr(i,7)*p(96)  + p(7)*dpdr(i,96) - dpdr(i,332)
      dpdr(i,370) = dpdr(i,17)*p(55) + p(17)*dpdr(i,55)  
      dpdr(i,371) = dpdr(i,2)*p(159) + p(2)*dpdr(i,159) - dpdr(i,345)
     $ - dpdr(i,349)
      dpdr(i,372) = dpdr(i,2)*p(147) + p(2)*dpdr(i,147) - dpdr(i,313)
      dpdr(i,373) = dpdr(i,11)*p(58) + p(11)*dpdr(i,58) - dpdr(i,315)
     $ - dpdr(i,328)
      dpdr(i,374) = dpdr(i,2)*p(148) + p(2)*dpdr(i,148) - dpdr(i,329)
     $ - dpdr(i,342)
      dpdr(i,375) = dpdr(i,6)*p(98)  + p(6)*dpdr(i,98) - dpdr(i,327)
      dpdr(i,376) = dpdr(i,2)*p(166) + p(2)*dpdr(i,166) - dpdr(i,359)
      dpdr(i,377) = dpdr(i,14)*p(54) + p(14)*dpdr(i,54) - dpdr(i,323)
      dpdr(i,378) = dpdr(i,2)*p(149) + p(2)*dpdr(i,149) - dpdr(i,331)
     $ - dpdr(i,343)
      dpdr(i,379) = dpdr(i,7)*p(98)  + p(7)*dpdr(i,98) - dpdr(i,325)
      dpdr(i,380) = dpdr(i,7)*p(99)  + p(7)*dpdr(i,99) - dpdr(i,311)
     $ - dpdr(i,359)
      dpdr(i,381) = dpdr(i,2)*p(167) + p(2)*dpdr(i,167) - dpdr(i,362)
      dpdr(i,382) = dpdr(i,14)*p(55) + p(14)*dpdr(i,55) - dpdr(i,323)
      dpdr(i,383) = dpdr(i,6)*p(100) + p(6)*dpdr(i,100) - dpdr(i,312)
     $ - dpdr(i,362)
      dpdr(i,384) = dpdr(i,11)*p(66) + p(11)*dpdr(i,66) - dpdr(i,344)
      dpdr(i,385) = dpdr(i,6)*p(101) + p(6)*dpdr(i,101) - dpdr(i,325)
     $ - dpdr(i,343) - dpdr(i,379)
      dpdr(i,386) = dpdr(i,7)*p(101) + p(7)*dpdr(i,101) - dpdr(i,342)
     $ - dpdr(i,327) - dpdr(i,375)
      dpdr(i,387) = dpdr(i,6)*p(103) + p(6)*dpdr(i,103) - dpdr(i,313)
     $ - dpdr(i,382)
      dpdr(i,388) = dpdr(i,11)*p(67) + p(11)*dpdr(i,67) - dpdr(i,314)
      dpdr(i,389) = dpdr(i,2)*p(152) + p(2)*dpdr(i,152) - dpdr(i,333)
      dpdr(i,390) = dpdr(i,2)*p(153) + p(2)*dpdr(i,153) - dpdr(i,334)
      dpdr(i,391) = dpdr(i,2)*p(154) + p(2)*dpdr(i,154) - dpdr(i,315)
     $ - dpdr(i,373)
      dpdr(i,392) = dpdr(i,1)*p(240) + p(1)*dpdr(i,240) - dpdr(i,335)
     $ - dpdr(i,373) - dpdr(i,366)
      dpdr(i,393) = dpdr(i,2)*p(155) + p(2)*dpdr(i,155) - dpdr(i,338)
     $ - dpdr(i,342)
      dpdr(i,394) = dpdr(i,2)*p(171) + p(2)*dpdr(i,171) - dpdr(i,377)
      dpdr(i,395) = dpdr(i,2)*p(157) + p(2)*dpdr(i,157) - dpdr(i,340)
     $ - dpdr(i,343)
      dpdr(i,396) = dpdr(i,2)*p(163) + p(2)*dpdr(i,163) - dpdr(i,350)
     $ - dpdr(i,351)
      dpdr(i,397) = dpdr(i,6)*p(95)  + p(6)*dpdr(i,95) - dpdr(i,315)
     $ - dpdr(i,350) - dpdr(i,340)
      dpdr(i,398) = dpdr(i,2)*p(172) + p(2)*dpdr(i,172) - dpdr(i,382)
      dpdr(i,399) = dpdr(i,7)*p(95)  + p(7)*dpdr(i,95) - dpdr(i,315)
     $ - dpdr(i,350) - dpdr(i,338)
      dpdr(i,400) = dpdr(i,11)*p(68) + p(11)*dpdr(i,68) - dpdr(i,345)
     $ - dpdr(i,349)
      dpdr(i,401) = dpdr(i,2)*p(175) + p(2)*dpdr(i,175) - dpdr(i,380)
     $ - dpdr(i,385)
      dpdr(i,402) = dpdr(i,6)*p(106) + p(6)*dpdr(i,106) - dpdr(i,335)
     $ - dpdr(i,343) - dpdr(i,396)
      dpdr(i,403) = dpdr(i,2)*p(176) + p(2)*dpdr(i,176) - dpdr(i,383)
     $ - dpdr(i,386)
      dpdr(i,404) = dpdr(i,7)*p(106) + p(7)*dpdr(i,106) - dpdr(i,342)
     $ - dpdr(i,335) - dpdr(i,396)
      dpdr(i,405) = dpdr(i,4)*p(150) + p(4)*dpdr(i,150) - dpdr(i,328)
     $ - dpdr(i,359) - dpdr(i,362)
      dpdr(i,406) = dpdr(i,11)*p(69) + p(11)*dpdr(i,69) - dpdr(i,316)
      dpdr(i,407) = dpdr(i,2)*p(161) + p(2)*dpdr(i,161) - dpdr(i,347)
      dpdr(i,408) = dpdr(i,2)*p(162) + p(2)*dpdr(i,162) - dpdr(i,348)
      dpdr(i,409) = dpdr(i,11)*p(70) + p(11)*dpdr(i,70) - dpdr(i,350)
     $ - dpdr(i,351)
      dpdr(i,410) = dpdr(i,2)*p(180) + p(2)*dpdr(i,180) - dpdr(i,397)
     $ - dpdr(i,402)
      dpdr(i,411) = dpdr(i,6)*p(110) + p(6)*dpdr(i,110) - dpdr(i,348)
      dpdr(i,412) = dpdr(i,2)*p(181) + p(2)*dpdr(i,181) - dpdr(i,399)
     $ - dpdr(i,404)
      dpdr(i,413) = dpdr(i,7)*p(110) + p(7)*dpdr(i,110) - dpdr(i,347)
      dpdr(i,414) = dpdr(i,6)*p(112) + p(6)*dpdr(i,112) - dpdr(i,316)
     $ - dpdr(i,398) - dpdr(i,316)
      dpdr(i,415) = dpdr(i,11)*p(71) + p(11)*dpdr(i,71) - dpdr(i,352)
      dpdr(i,416) = dpdr(i,2)*p(184) + p(2)*dpdr(i,184) - dpdr(i,411)
      dpdr(i,417) = dpdr(i,2)*p(185) + p(2)*dpdr(i,185) - dpdr(i,413)
      dpdr(i,418) = dpdr(i,12)*p(71) + p(12)*dpdr(i,71) - dpdr(i,351)
     $ - dpdr(i,416) - dpdr(i,417)
      dpdr(i,419) = dpdr(i,2)*p(164) + p(2)*dpdr(i,164) - dpdr(i,320)
     $ - dpdr(i,356)
      dpdr(i,420) = dpdr(i,2)*p(165) + p(2)*dpdr(i,165) - dpdr(i,328)
     $ - dpdr(i,373)
      dpdr(i,421) = dpdr(i,2)*p(170) + p(2)*dpdr(i,170) - dpdr(i,337)
     $ - dpdr(i,392)
      dpdr(i,422) = dpdr(i,1)*p(268) + p(1)*dpdr(i,268) - dpdr(i,359)
      dpdr(i,423) = dpdr(i,1)*p(269) + p(1)*dpdr(i,269) - dpdr(i,362)
      dpdr(i,424) = dpdr(i,2)*p(174) + p(2)*dpdr(i,174) - dpdr(i,345)
     $ - dpdr(i,400)
      dpdr(i,425) = dpdr(i,6)*p(102) + p(6)*dpdr(i,102) - dpdr(i,345)
     $ - dpdr(i,326)
      dpdr(i,426) = dpdr(i,7)*p(103) + p(7)*dpdr(i,103) - dpdr(i,345)
     $ - dpdr(i,330)
      dpdr(i,427) = dpdr(i,3)*p(186) + p(3)*dpdr(i,186) - dpdr(i,364)
      dpdr(i,428) = dpdr(i,18)*p(65) + p(18)*dpdr(i,65) - dpdr(i,354)
      dpdr(i,429) = dpdr(i,17)*p(66) + p(17)*dpdr(i,66) - dpdr(i,361)
      dpdr(i,430) = dpdr(i,2)*p(179) + p(2)*dpdr(i,179) - dpdr(i,350)
     $ - dpdr(i,409)
      dpdr(i,431) = dpdr(i,4)*p(166) + p(4)*dpdr(i,166) - dpdr(i,361)
     $ - dpdr(i,357) - dpdr(i,422)
      dpdr(i,432) = dpdr(i,4)*p(167) + p(4)*dpdr(i,167) - dpdr(i,361)
     $ - dpdr(i,360) - dpdr(i,423)
      dpdr(i,433) = dpdr(i,2)*p(187) + p(2)*dpdr(i,187) - dpdr(i,387)
      dpdr(i,434) = dpdr(i,14)*p(66) + p(14)*dpdr(i,66) - dpdr(i,328)
     $ - dpdr(i,405) - dpdr(i,376) - dpdr(i,381) - dpdr(i,373)
      dpdr(i,435) = dpdr(i,1)*p(277) + p(1)*dpdr(i,277) - dpdr(i,387)
     $ - dpdr(i,385) - dpdr(i,425)
      dpdr(i,436) = dpdr(i,1)*p(278) + p(1)*dpdr(i,278) - dpdr(i,387)
     $ - dpdr(i,386) - dpdr(i,426)
      dpdr(i,437) = dpdr(i,2)*p(177) + p(2)*dpdr(i,177) - dpdr(i,346)
     $ - dpdr(i,391)
      dpdr(i,438) = dpdr(i,2)*p(178) + p(2)*dpdr(i,178) - dpdr(i,349)
     $ - dpdr(i,400)
      dpdr(i,439) = dpdr(i,1)*p(281) + p(1)*dpdr(i,281) - dpdr(i,396)
     $ - dpdr(i,392) - dpdr(i,430)
      dpdr(i,440) = dpdr(i,21)*p(54) + p(21)*dpdr(i,54) - dpdr(i,370)
      dpdr(i,441) = dpdr(i,21)*p(55) + p(21)*dpdr(i,55) - dpdr(i,368)
      dpdr(i,442) = dpdr(i,2)*p(189) + p(2)*dpdr(i,189) - dpdr(i,405)
     $ - dpdr(i,434)
      dpdr(i,443) = dpdr(i,1)*p(285) + p(1)*dpdr(i,285) - dpdr(i,402)
     $ - dpdr(i,404) - dpdr(i,400) - dpdr(i,434)  - dpdr(i,430)
      dpdr(i,444) = dpdr(i,6)*p(116) + p(6)*dpdr(i,116) - dpdr(i,347)
     $ - dpdr(i,413) - dpdr(i,403)
      dpdr(i,445) = dpdr(i,7)*p(116) + p(7)*dpdr(i,116) - dpdr(i,348)
     $ - dpdr(i,411) - dpdr(i,401)
      dpdr(i,446) = dpdr(i,2)*p(182) + p(2)*dpdr(i,182) - dpdr(i,351)
     $ - dpdr(i,409)
      dpdr(i,447) = dpdr(i,2)*p(191) + p(2)*dpdr(i,191) - dpdr(i,414)
     $ - dpdr(i,443)
      dpdr(i,448) = dpdr(i,3)*p(183) + p(3)*dpdr(i,183) - dpdr(i,351)
     $ - dpdr(i,418) - dpdr(i,411) - dpdr(i,413) - dpdr(i,409)
      dpdr(i,449) = dpdr(i,6)*p(118) + p(6)*dpdr(i,118) - dpdr(i,351)
     $ - dpdr(i,418) - dpdr(i,412)
      dpdr(i,450) = dpdr(i,7)*p(118) + p(7)*dpdr(i,118) - dpdr(i,351)
     $ - dpdr(i,418) - dpdr(i,410)
      dpdr(i,451) = dpdr(i,2)*p(193) + p(2)*dpdr(i,193) - dpdr(i,418)
     $ - dpdr(i,448)
      dpdr(i,452) = dpdr(i,1)*p(294) + p(1)*dpdr(i,294) - dpdr(i,418)
     $ - dpdr(i,415) - dpdr(i,448)
      dpdr(i,453) = dpdr(i,6)*p(119) + p(6)*dpdr(i,119) - dpdr(i,417)
      dpdr(i,454) = dpdr(i,7)*p(119) + p(7)*dpdr(i,119) - dpdr(i,416)
      dpdr(i,455) = dpdr(i,2)*p(186) + p(2)*dpdr(i,186) - dpdr(i,363)
      dpdr(i,456) = dpdr(i,3)*p(187) + p(3)*dpdr(i,187) - dpdr(i,385)
     $ - dpdr(i,386) - dpdr(i,384) - dpdr(i,435) - dpdr(i,436)
     $ - dpdr(i,434)
      dpdr(i,457) = dpdr(i,2)*p(188) + p(2)*dpdr(i,188) - dpdr(i,388)
      dpdr(i,458) = dpdr(i,4)*p(187) + p(4)*dpdr(i,187) - dpdr(i,425)
     $ - dpdr(i,426) - dpdr(i,424)
      dpdr(i,459) = dpdr(i,2)*p(190) + p(2)*dpdr(i,190) - dpdr(i,406)
      dpdr(i,460) = dpdr(i,21)*p(66) + p(21)*dpdr(i,66) - dpdr(i,422)
     $ - dpdr(i,423) - dpdr(i,420)
      dpdr(i,461) = dpdr(i,2)*p(192) + p(2)*dpdr(i,192) - dpdr(i,415)
      dpdr(i,462) = dpdr(i,18)*p(71) + p(18)*dpdr(i,71) - dpdr(i,440)
     $ - dpdr(i,441) - dpdr(i,438)
      dpdr(i,463) = dpdr(i,2)*p(194) + p(2)*dpdr(i,194) - dpdr(i,452)
      dpdr(i,464) = dpdr(i,3)*p(194) + p(3)*dpdr(i,194) - dpdr(i,453)
     $ - dpdr(i,454) - dpdr(i,451)
      dpdr(i,465) = dpdr(i,1)*p(305) + p(1)*dpdr(i,305) - dpdr(i,464)
     $ - dpdr(i,463)
         
      enddo

      return

      end subroutine EvdPdR
      subroutine evdbdr
***********************************************************************
*  This subroutine eliminates the 2-body terms in Bowman's approach.
***********************************************************************

      implicit double precision (a-h,o-z)

      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:465),dVdR(6),
     $                  dRdX(6,12),dBdR(6,430)
      
      integer i
      double precision db1dr(6,466) 

C Pass P(0:465) to BM1(1:466)
      do j=1,6
      do i=1,466
        db1dr(j,i)=dpdr(j,i-1)
      enddo
      enddo

C Remove unconnected terms and 2-body terms and pass to B(1:430)
      do j=1,6
      dbdr(j,1)=db1dr(j,4)

      do i=2,4
        dbdr(j,i)=db1dr(j,i+4)
      enddo

      dbdr(j,5)=db1dr(j,10)

      do i=6,11
        dbdr(j,i)=db1dr(j,i+6)
      enddo

      dbdr(j,12)=db1dr(j,19)
      dbdr(j,13)=db1dr(j,21)

      do i=14,26
        dbdr(j,i)=db1dr(j,i+9)
      enddo

      dbdr(j,27)=db1dr(j,37)
      dbdr(j,28)=db1dr(j,39)

      do i=29,53
        dbdr(j,i)=db1dr(j,i+12)
      enddo

      dbdr(j,54)=db1dr(j,67)
      dbdr(j,55)=db1dr(j,69)
      dbdr(j,56)=db1dr(j,71)

      do i=57,97
        dbdr(j,i)=db1dr(j,i+16)
      enddo

      dbdr(j,98)=db1dr(j,115)
      dbdr(j,99)=db1dr(j,117)
      dbdr(j,100)=db1dr(j,119)
      
      do i=101,166
        dbdr(j,i)=db1dr(j,i+20)
      enddo

      dbdr(j,167)=db1dr(j,188)
      dbdr(j,168)=db1dr(j,190)
      dbdr(j,169)=db1dr(j,192)
      dbdr(j,170)=db1dr(j,194)

      do i=171,272
        dbdr(j,i)=db1dr(j,i+25)
      enddo

      dbdr(j,273)=db1dr(j,299)
      dbdr(j,274)=db1dr(j,301)
      dbdr(j,275)=db1dr(j,303)
      dbdr(j,276)=db1dr(j,305)

      do i=277,425
        dbdr(j,i)=db1dr(j,i+30)
      enddo

      dbdr(j,426)=db1dr(j,457)
      dbdr(j,427)=db1dr(j,459)
      dbdr(j,428)=db1dr(j,461)
      dbdr(j,429)=db1dr(j,463)
      dbdr(j,430)=db1dr(j,465) 
      enddo

      return

      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Dispersion correction based on Grimme's D3(BJ) calculation for
C diatomic pairs
C
C Several subroutines of DFTD3 V3.1 Rev 1 by Grimme were merged into 
C subroutine edisp and they have been heavily modified to calculate
C only dispersion energy corrections that are needed.
C
C S. Grimme, J. Antony, S. Ehrlich and H. Krieg
C J. Chem. Phys, 132 (2010), 154104
C and 
C S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011),
C 1456-1465
C
C The C6 values are fixed.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine d3disp(dist,disp,dispdr,igrad)

      double precision cn(4),s6,s8,rs6,rs8
      double precision dist(6), e6(6), e8(6), disp, dispdr(6), c6(6)
      double precision e6dr(6),e8dr(6)
      integer iz(4), mxc(94), i, j, igrad
      double precision c6ab(94,94,5,5,3)
      double precision r2r4(94)
      double precision autoang,autokcal

      autoang =0.52917726d0
      autokcal=627.509541d0

! Generalized parameters for BJ damping
      s6= 1.0d0
      s8= 2.0d0
      rs6= 0.5299d0
      rs8= 2.20d0

      do i=1,6
      dist(i)=dist(i)/autoang
      enddo

C iz for O4 system
      iz(1)=8
      iz(2)=8
      iz(3)=8
      iz(4)=8
C C6 for O4 system
      c6(1)=12.8d0
      c6(2)=12.8d0
      c6(3)=12.8d0
      c6(4)=12.8d0
      c6(5)=12.8d0
      c6(6)=12.8d0

C Calculate dispersion correction
      call edisp(94,5,4,dist,iz,mxc,
     .     rs6,rs8,e6,e8,e6dr,e8dr,c6,0)

      disp = 0.0d0

      do i=1,6
      disp =disp + (-s6*e6(i)-s8*e8(i))*autokcal
      enddo

      if (igrad .eq. 1) then
      call edisp(94,5,4,dist,iz,mxc,
     .     rs6,rs8,e6,e8,e6dr,e8dr,c6,1)

      dispdr(:) = 0.0d0

      do i=1,6
      dispdr(i) =dispdr(i) + (-s6*e6dr(i)-s8*e8dr(i))*autokcal/autoang
      enddo
      endif

      do i=1,6
      dist(i)=dist(i)*autoang
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute energy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine edisp(max_elem,maxc,n,dist,iz,mxc,
     .           rs6,rs8,e6,e8,e6dr,e8dr,c6a,igrad)

      implicit none  
      integer n,iz(4),max_elem,maxc,mxc(max_elem) 
      double precision dist(6),r2r4(max_elem),r0ab(max_elem,max_elem)
      double precision rs6,rs8,rcov(max_elem)
      double precision c6ab(max_elem,max_elem,maxc,maxc,3)
      double precision e6(6), e8(6), c6a(6), e6dr(6), e8dr(6)
       
      integer iat,jat,igrad
      double precision r,tmp,c6,c8,a1,a2
      double precision damp6,damp8
      double precision cn(n)                             
      double precision r2ab(n*n),cc6ab(n*n),dmp(n*n)
      integer step

      e6(:) =0.0d0
      e8(:) =0.0d0

      e6dr(:) =0.0d0
      e8dr(:) =0.0d0

      a1=rs6
      a2=rs8 

!  r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
      r2r4(1:94)=(/
     . 2.00734898,  1.56637132,  5.01986934,  3.85379032,  3.64446594,
     . 3.10492822,  2.71175247,  2.59361680,  2.38825250,  2.21522516,
     . 6.58585536,  5.46295967,  5.65216669,  4.88284902,  4.29727576,
     . 4.04108902,  3.72932356,  3.44677275,  7.97762753,  7.07623947,
     . 6.60844053,  6.28791364,  6.07728703,  5.54643096,  5.80491167,
     . 5.58415602,  5.41374528,  5.28497229,  5.22592821,  5.09817141,
     . 6.12149689,  5.54083734,  5.06696878,  4.87005108,  4.59089647,
     . 4.31176304,  9.55461698,  8.67396077,  7.97210197,  7.43439917,
     . 6.58711862,  6.19536215,  6.01517290,  5.81623410,  5.65710424,
     . 5.52640661,  5.44263305,  5.58285373,  7.02081898,  6.46815523,
     . 5.98089120,  5.81686657,  5.53321815,  5.25477007, 11.02204549,
     .10.15679528,  9.35167836,  9.06926079,  8.97241155,  8.90092807,
     . 8.85984840,  8.81736827,  8.79317710,  7.89969626,  8.80588454,
     . 8.42439218,  8.54289262,  8.47583370,  8.45090888,  8.47339339,
     . 7.83525634,  8.20702843,  7.70559063,  7.32755997,  7.03887381,
     . 6.68978720,  6.05450052,  5.88752022,  5.70661499,  5.78450695,
     . 7.79780729,  7.26443867,  6.78151984,  6.67883169,  6.39024318,
     . 6.09527958, 11.79156076, 11.10997644,  9.51377795,  8.67197068,
     . 8.77140725,  8.65402716,  8.53923501,  8.85024712 /)

! these new data are scaled with k2=4./3. and converted to a_0 via
! autoang=0.52917726d0
      rcov(1:94)=(/
     . 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865,
     . 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527,
     . 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820,
     . 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730,
     . 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923,
     . 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188,
     . 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349,
     . 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216,
     . 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,
     . 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967,
     . 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625,
     . 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657,
     . 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833,
     . 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098,
     . 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878,
     . 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790,
     . 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584,
     . 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,
     . 3.82984466, 3.85504098, 3.88023730, 3.90543362 /)

C DFT-D3
      step=0
      do iat=1,n-1
         do jat=iat+1,n
         step=step+1
         r=dist(step)
         c6=c6a(step)
c r2r4 stored in main as sqrt
         c8 =3.0d0*c6*r2r4(iz(iat))*r2r4(iz(jat))

c energy for BJ damping
          tmp=sqrt(c8/c6)              
          e6(step)= c6/(r**6+(a1*tmp+a2)**6)
          e8(step)= c8/(r**8+(a1*tmp+a2)**8)
C calculate gradients
         if (igrad .eq. 1) then
c grad for BJ damping
          e6dr(step)=c6*(-6*r**5)/(r**6+(a1*tmp+a2)**6)**2
          e8dr(step)=c8*(-8*r**7)/(r**8+(a1*tmp+a2)**8)**2
         endif
         enddo
      enddo

      end subroutine edisp

C Begin
       block data prmt

      implicit double precision (a-h,o-z)

      common /coord/    R(6)
      common /epot/     rMs(6),rM(0:111),P(0:465),C(430),B(430)
      common /gradt/    dMsdR(6,6),dMdR(6,0:111),dPdR(6,0:465),dVdR(6),
     $                  dRdX(6,12),dBdR(6,430)
      common /msprmt/   a,ab,ra,rb

C Nonlinear parameters of a(1.0 Ang) and ab (1.5 Ang^2) and re (1.208
C Ang)
      data a    /0.9d0/
      data ab    /4.5d0/
      data ra   /1.3d0/
      data rb   /1.1d0/

C Linear parameters optimized by the weighted-least square fitting
      data C(   1)   /   0.434720431798D+02 /
      data C(   2)   /  -0.554668579283D+03 /
      data C(   3)   /   0.257445901353D+04 /
      data C(   4)   /  -0.213807688288D+04 /
      data C(   5)   /  -0.248886942161D+03 /
      data C(   6)   /   0.719084654648D+03 /
      data C(   7)   /  -0.785231200932D+03 /
      data C(   8)   /   0.268844019872D+04 /
      data C(   9)   /  -0.894783391739D+03 /
      data C(  10)   /  -0.357957764652D+04 /
      data C(  11)   /   0.362675454825D+04 /
      data C(  12)   /   0.102629622676D+04 /
      data C(  13)   /   0.795258994983D+03 /
      data C(  14)   /  -0.946936811683D+04 /
      data C(  15)   /   0.571579383749D+04 /
      data C(  16)   /   0.934372178824D+03 /
      data C(  17)   /  -0.292732443657D+04 /
      data C(  18)   /   0.219656326428D+04 /
      data C(  19)   /  -0.676509263255D+04 /
      data C(  20)   /  -0.282505207332D+04 /
      data C(  21)   /   0.104924792198D+05 /
      data C(  22)   /  -0.658586513846D+04 /
      data C(  23)   /  -0.432459822098D+04 /
      data C(  24)   /   0.316198143890D+04 /
      data C(  25)   /  -0.337163445639D+04 /
      data C(  26)   /  -0.713571319552D+03 /
      data C(  27)   /  -0.134555448022D+04 /
      data C(  28)   /  -0.150137169111D+04 /
      data C(  29)   /  -0.358000122745D+04 /
      data C(  30)   /   0.140590867542D+05 /
      data C(  31)   /  -0.459125633418D+04 /
      data C(  32)   /  -0.550769788422D+04 /
      data C(  33)   /   0.217230838750D+04 /
      data C(  34)   /  -0.632774344127D+03 /
      data C(  35)   /  -0.100999591686D+05 /
      data C(  36)   /   0.169165667764D+04 /
      data C(  37)   /   0.167842440492D+04 /
      data C(  38)   /  -0.257911596712D+05 /
      data C(  39)   /   0.127815645611D+05 /
      data C(  40)   /  -0.537883442310D+03 /
      data C(  41)   /  -0.133764280242D+04 /
      data C(  42)   /  -0.212370345146D+04 /
      data C(  43)   /   0.161754184625D+05 /
      data C(  44)   /   0.355956338618D+04 /
      data C(  45)   /   0.913378256342D+04 /
      data C(  46)   /   0.610950125310D+04 /
      data C(  47)   /   0.515181004953D+04 /
      data C(  48)   /  -0.149547604838D+05 /
      data C(  49)   /   0.447834665993D+04 /
      data C(  50)   /   0.289335443481D+04 /
      data C(  51)   /  -0.430236737245D+04 /
      data C(  52)   /   0.115862975164D+05 /
      data C(  53)   /  -0.223649282872D+04 /
      data C(  54)   /  -0.424281221393D+04 /
      data C(  55)   /   0.607690277921D+04 /
      data C(  56)   /   0.682307953653D+03 /
      data C(  57)   /  -0.146643434299D+05 /
      data C(  58)   /  -0.995963785829D+04 /
      data C(  59)   /   0.151656352036D+05 /
      data C(  60)   /  -0.713080950000D+04 /
      data C(  61)   /   0.992264149235D+04 /
      data C(  62)   /  -0.225192674918D+05 /
      data C(  63)   /  -0.227521700630D+05 /
      data C(  64)   /   0.356736751019D+03 /
      data C(  65)   /  -0.149009596280D+05 /
      data C(  66)   /   0.547202199914D+04 /
      data C(  67)   /  -0.562297176466D+04 /
      data C(  68)   /   0.533482100629D+04 /
      data C(  69)   /   0.173076951226D+05 /
      data C(  70)   /  -0.470318017955D+04 /
      data C(  71)   /   0.252423648786D+04 /
      data C(  72)   /   0.397204506011D+04 /
      data C(  73)   /   0.149604469867D+05 /
      data C(  74)   /   0.261468562194D+04 /
      data C(  75)   /   0.148832634817D+04 /
      data C(  76)   /  -0.525471073088D+04 /
      data C(  77)   /   0.393212581141D+05 /
      data C(  78)   /  -0.306109210725D+05 /
      data C(  79)   /   0.721202884138D+04 /
      data C(  80)   /   0.352945175290D+04 /
      data C(  81)   /  -0.784511969410D+04 /
      data C(  82)   /   0.201473287242D+03 /
      data C(  83)   /   0.124609448656D+05 /
      data C(  84)   /  -0.103099257764D+05 /
      data C(  85)   /   0.123788976352D+04 /
      data C(  86)   /  -0.845440125549D+03 /
      data C(  87)   /  -0.652328346462D+04 /
      data C(  88)   /  -0.397970507273D+04 /
      data C(  89)   /  -0.813327575032D+04 /
      data C(  90)   /  -0.854960333368D+04 /
      data C(  91)   /  -0.146818322655D+05 /
      data C(  92)   /   0.250112509218D+05 /
      data C(  93)   /  -0.120475653865D+04 /
      data C(  94)   /   0.894965121433D+03 /
      data C(  95)   /   0.687552457372D+04 /
      data C(  96)   /  -0.155946821955D+05 /
      data C(  97)   /  -0.556039535294D+03 /
      data C(  98)   /   0.412239459401D+04 /
      data C(  99)   /  -0.103729683198D+05 /
      data C( 100)   /   0.120021137587D+04 /
      data C( 101)   /  -0.106685046856D+05 /
      data C( 102)   /   0.580461159121D+04 /
      data C( 103)   /   0.323628855252D+05 /
      data C( 104)   /   0.801279078904D+03 /
      data C( 105)   /  -0.144743305379D+05 /
      data C( 106)   /   0.245891704416D+05 /
      data C( 107)   /   0.734963482362D+04 /
      data C( 108)   /   0.579008306105D+04 /
      data C( 109)   /  -0.209003437783D+05 /
      data C( 110)   /   0.788390162743D+03 /
      data C( 111)   /   0.795154017688D+03 /
      data C( 112)   /  -0.218632472140D+05 /
      data C( 113)   /   0.172512443045D+05 /
      data C( 114)   /   0.540100143425D+05 /
      data C( 115)   /   0.338674945681D+04 /
      data C( 116)   /  -0.386001569801D+05 /
      data C( 117)   /   0.239712386089D+04 /
      data C( 118)   /  -0.613824367400D+04 /
      data C( 119)   /  -0.488209303900D+04 /
      data C( 120)   /   0.701192258148D+04 /
      data C( 121)   /   0.167984032510D+05 /
      data C( 122)   /  -0.259650106547D+05 /
      data C( 123)   /   0.177986468669D+05 /
      data C( 124)   /   0.150944386869D+05 /
      data C( 125)   /   0.121802007803D+05 /
      data C( 126)   /  -0.978234467399D+04 /
      data C( 127)   /   0.184330692467D+05 /
      data C( 128)   /   0.921531962740D+04 /
      data C( 129)   /   0.183336479381D+05 /
      data C( 130)   /   0.587334293511D+02 /
      data C( 131)   /  -0.167218462260D+05 /
      data C( 132)   /  -0.154061569920D+04 /
      data C( 133)   /  -0.198628866071D+05 /
      data C( 134)   /  -0.595045726927D+04 /
      data C( 135)   /  -0.614631655966D+04 /
      data C( 136)   /  -0.643860187595D+04 /
      data C( 137)   /  -0.328961702219D+05 /
      data C( 138)   /  -0.312220303325D+04 /
      data C( 139)   /  -0.117173908650D+05 /
      data C( 140)   /   0.102108801243D+05 /
      data C( 141)   /  -0.296838514570D+05 /
      data C( 142)   /   0.316743949442D+05 /
      data C( 143)   /  -0.547675942794D+04 /
      data C( 144)   /  -0.428760549879D+03 /
      data C( 145)   /  -0.281670316262D+03 /
      data C( 146)   /  -0.150670294917D+05 /
      data C( 147)   /  -0.251378269785D+05 /
      data C( 148)   /  -0.109642831245D+05 /
      data C( 149)   /   0.452857325822D+04 /
      data C( 150)   /   0.282458960852D+04 /
      data C( 151)   /   0.177576149909D+04 /
      data C( 152)   /   0.461479349791D+05 /
      data C( 153)   /   0.677295737371D+04 /
      data C( 154)   /  -0.544615362371D+04 /
      data C( 155)   /   0.797407925067D+04 /
      data C( 156)   /   0.337166459111D+04 /
      data C( 157)   /   0.531810770795D+04 /
      data C( 158)   /   0.415106816445D+04 /
      data C( 159)   /   0.842950842709D+04 /
      data C( 160)   /   0.131227142369D+05 /
      data C( 161)   /  -0.278620017629D+05 /
      data C( 162)   /  -0.467327134275D+04 /
      data C( 163)   /  -0.274452920819D+04 /
      data C( 164)   /  -0.732936170936D+04 /
      data C( 165)   /   0.128933968259D+05 /
      data C( 166)   /   0.500270166783D+04 /
      data C( 167)   /  -0.832810100467D+04 /
      data C( 168)   /   0.231068303833D+04 /
      data C( 169)   /   0.736066804875D+04 /
      data C( 170)   /  -0.167133346329D+04 /
      data C( 171)   /   0.271382328852D+05 /
      data C( 172)   /  -0.299226497133D+04 /
      data C( 173)   /  -0.300405606047D+05 /
      data C( 174)   /  -0.417234719382D+05 /
      data C( 175)   /  -0.925938352157D+04 /
      data C( 176)   /   0.150738290565D+05 /
      data C( 177)   /   0.621392381200D+04 /
      data C( 178)   /  -0.250466052381D+05 /
      data C( 179)   /   0.189034139443D+05 /
      data C( 180)   /   0.645847033559D+04 /
      data C( 181)   /  -0.122509319704D+04 /
      data C( 182)   /   0.674497164564D+04 /
      data C( 183)   /   0.165909951935D+04 /
      data C( 184)   /   0.314085526715D+05 /
      data C( 185)   /  -0.567978068986D+04 /
      data C( 186)   /   0.804344254374D+04 /
      data C( 187)   /   0.718296953579D+04 /
      data C( 188)   /   0.805081988188D+04 /
      data C( 189)   /   0.959057200505D+04 /
      data C( 190)   /   0.173175419540D+05 /
      data C( 191)   /   0.192224951993D+04 /
      data C( 192)   /  -0.897479623343D+04 /
      data C( 193)   /   0.489185253916D+04 /
      data C( 194)   /  -0.176425900434D+05 /
      data C( 195)   /   0.373806080219D+04 /
      data C( 196)   /  -0.499566648122D+04 /
      data C( 197)   /  -0.453676607117D+04 /
      data C( 198)   /  -0.214019574689D+05 /
      data C( 199)   /  -0.760378970714D+04 /
      data C( 200)   /  -0.291704914069D+04 /
      data C( 201)   /   0.237429310741D+05 /
      data C( 202)   /  -0.111950792650D+05 /
      data C( 203)   /   0.293253496618D+05 /
      data C( 204)   /   0.789023058401D+03 /
      data C( 205)   /  -0.424411435768D+04 /
      data C( 206)   /  -0.117087625656D+05 /
      data C( 207)   /   0.694946236291D+04 /
      data C( 208)   /  -0.149281145977D+05 /
      data C( 209)   /  -0.136623142423D+04 /
      data C( 210)   /  -0.141052998243D+05 /
      data C( 211)   /  -0.102173374749D+05 /
      data C( 212)   /  -0.566906448414D+03 /
      data C( 213)   /   0.789959020997D+04 /
      data C( 214)   /  -0.391936833293D+04 /
      data C( 215)   /   0.999203207828D+04 /
      data C( 216)   /   0.565812874635D+03 /
      data C( 217)   /  -0.882979001006D+04 /
      data C( 218)   /   0.421173930005D+04 /
      data C( 219)   /  -0.625164802951D+04 /
      data C( 220)   /   0.142249936194D+04 /
      data C( 221)   /  -0.122348255800D+05 /
      data C( 222)   /   0.144730122892D+05 /
      data C( 223)   /   0.716888617746D+02 /
      data C( 224)   /  -0.296205282215D+04 /
      data C( 225)   /  -0.135157525654D+05 /
      data C( 226)   /   0.282282739575D+04 /
      data C( 227)   /  -0.427724924831D+04 /
      data C( 228)   /   0.112055827751D+05 /
      data C( 229)   /   0.570635552195D+04 /
      data C( 230)   /  -0.360280649223D+04 /
      data C( 231)   /   0.102208442500D+05 /
      data C( 232)   /   0.102712354122D+05 /
      data C( 233)   /   0.261880233029D+04 /
      data C( 234)   /   0.756346425878D+04 /
      data C( 235)   /   0.164971200139D+05 /
      data C( 236)   /   0.145688784480D+04 /
      data C( 237)   /   0.296086198502D+04 /
      data C( 238)   /  -0.312536908014D+04 /
      data C( 239)   /   0.121650939853D+05 /
      data C( 240)   /  -0.166454887214D+05 /
      data C( 241)   /   0.812412831947D+02 /
      data C( 242)   /  -0.136497078414D+04 /
      data C( 243)   /   0.870474131812D+04 /
      data C( 244)   /  -0.368363201039D+05 /
      data C( 245)   /   0.255189655169D+05 /
      data C( 246)   /   0.206043389737D+02 /
      data C( 247)   /   0.488181964434D+03 /
      data C( 248)   /   0.822653317496D+04 /
      data C( 249)   /   0.115473149939D+04 /
      data C( 250)   /   0.180422606051D+05 /
      data C( 251)   /  -0.187109086247D+04 /
      data C( 252)   /   0.200831958953D+04 /
      data C( 253)   /  -0.144908090443D+05 /
      data C( 254)   /   0.369931716826D+04 /
      data C( 255)   /  -0.251207665591D+04 /
      data C( 256)   /  -0.378840937076D+04 /
      data C( 257)   /  -0.251019035579D+04 /
      data C( 258)   /  -0.352218080889D+05 /
      data C( 259)   /  -0.157646059642D+04 /
      data C( 260)   /   0.135390050658D+04 /
      data C( 261)   /  -0.673592136935D+04 /
      data C( 262)   /   0.513499801117D+04 /
      data C( 263)   /  -0.597893982397D+04 /
      data C( 264)   /  -0.629318425103D+03 /
      data C( 265)   /  -0.327904556063D+04 /
      data C( 266)   /  -0.425819219862D+04 /
      data C( 267)   /   0.138231691435D+05 /
      data C( 268)   /   0.442022386425D+04 /
      data C( 269)   /   0.149938258302D+04 /
      data C( 270)   /   0.348504332442D+04 /
      data C( 271)   /  -0.561845437718D+04 /
      data C( 272)   /  -0.371238889983D+04 /
      data C( 273)   /   0.266879264931D+04 /
      data C( 274)   /  -0.241086213543D+04 /
      data C( 275)   /  -0.256487559688D+04 /
      data C( 276)   /   0.807623594819D+03 /
      data C( 277)   /  -0.146790064277D+05 /
      data C( 278)   /   0.555226660757D+04 /
      data C( 279)   /   0.235811102582D+05 /
      data C( 280)   /  -0.139884844872D+05 /
      data C( 281)   /   0.329907188988D+05 /
      data C( 282)   /  -0.805451894503D+04 /
      data C( 283)   /  -0.103162828961D+05 /
      data C( 284)   /   0.161264355953D+05 /
      data C( 285)   /   0.114990797542D+05 /
      data C( 286)   /  -0.462324369583D+04 /
      data C( 287)   /  -0.368515824108D+04 /
      data C( 288)   /  -0.316340302073D+04 /
      data C( 289)   /   0.515823486141D+03 /
      data C( 290)   /  -0.145836944077D+04 /
      data C( 291)   /  -0.244069295350D+04 /
      data C( 292)   /   0.486854300243D+04 /
      data C( 293)   /   0.840847301500D+04 /
      data C( 294)   /   0.104526133121D+04 /
      data C( 295)   /   0.884893182738D+03 /
      data C( 296)   /  -0.344498293927D+04 /
      data C( 297)   /  -0.208868965702D+05 /
      data C( 298)   /  -0.420935027569D+04 /
      data C( 299)   /  -0.838893564378D+04 /
      data C( 300)   /   0.744963159083D+04 /
      data C( 301)   /  -0.542353889751D+04 /
      data C( 302)   /  -0.195893556763D+04 /
      data C( 303)   /  -0.466724683081D+04 /
      data C( 304)   /   0.275184249484D+04 /
      data C( 305)   /   0.277522069152D+04 /
      data C( 306)   /  -0.662767970556D+04 /
      data C( 307)   /  -0.319233977664D+03 /
      data C( 308)   /  -0.896417468066D+03 /
      data C( 309)   /   0.123175057570D+05 /
      data C( 310)   /  -0.139755586021D+05 /
      data C( 311)   /   0.506568047186D+04 /
      data C( 312)   /  -0.258952868793D+04 /
      data C( 313)   /   0.264676924023D+04 /
      data C( 314)   /  -0.195443921477D+04 /
      data C( 315)   /  -0.127269946878D+05 /
      data C( 316)   /   0.828536232661D+04 /
      data C( 317)   /   0.214896348017D+04 /
      data C( 318)   /  -0.483068390482D+04 /
      data C( 319)   /   0.565494771451D+04 /
      data C( 320)   /   0.446561697944D+03 /
      data C( 321)   /   0.364782524646D+04 /
      data C( 322)   /   0.376605752474D+03 /
      data C( 323)   /   0.579822188327D+03 /
      data C( 324)   /  -0.211839304989D+04 /
      data C( 325)   /  -0.472264011837D+04 /
      data C( 326)   /   0.345025317650D+04 /
      data C( 327)   /  -0.676130963175D+01 /
      data C( 328)   /   0.269359548788D+04 /
      data C( 329)   /   0.276503090233D+04 /
      data C( 330)   /   0.177885809699D+05 /
      data C( 331)   /  -0.586550985360D+03 /
      data C( 332)   /  -0.542475290546D+04 /
      data C( 333)   /   0.509567352435D+04 /
      data C( 334)   /   0.113613052156D+04 /
      data C( 335)   /   0.170984112536D+03 /
      data C( 336)   /   0.150487958530D+04 /
      data C( 337)   /  -0.753940462606D+03 /
      data C( 338)   /   0.169631995550D+04 /
      data C( 339)   /  -0.121674596114D+05 /
      data C( 340)   /   0.285415113429D+04 /
      data C( 341)   /  -0.922400820349D+04 /
      data C( 342)   /   0.204864026100D+04 /
      data C( 343)   /  -0.871423146423D+03 /
      data C( 344)   /   0.248567624983D+04 /
      data C( 345)   /   0.614342185794D+04 /
      data C( 346)   /  -0.187662157284D+04 /
      data C( 347)   /  -0.857978013673D+04 /
      data C( 348)   /  -0.414781715613D+04 /
      data C( 349)   /   0.357346302603D+04 /
      data C( 350)   /   0.214164315004D+04 /
      data C( 351)   /  -0.590074931550D+04 /
      data C( 352)   /   0.120405323962D+04 /
      data C( 353)   /   0.645612000769D+03 /
      data C( 354)   /   0.136516824296D+04 /
      data C( 355)   /  -0.169012829652D+04 /
      data C( 356)   /   0.800622583324D+04 /
      data C( 357)   /   0.992615757371D+03 /
      data C( 358)   /   0.290405325075D+04 /
      data C( 359)   /  -0.598813578191D+03 /
      data C( 360)   /  -0.161195626593D+04 /
      data C( 361)   /   0.545457796504D+03 /
      data C( 362)   /  -0.340800286247D+04 /
      data C( 363)   /  -0.377880997467D+03 /
      data C( 364)   /   0.198426234996D+04 /
      data C( 365)   /   0.705249584289D+04 /
      data C( 366)   /   0.909789266216D+03 /
      data C( 367)   /  -0.119710848734D+04 /
      data C( 368)   /   0.662721663144D+04 /
      data C( 369)   /  -0.388094987001D+04 /
      data C( 370)   /  -0.120538116588D+04 /
      data C( 371)   /   0.121038056363D+04 /
      data C( 372)   /  -0.108996864000D+04 /
      data C( 373)   /  -0.554436730351D+04 /
      data C( 374)   /   0.136234762955D+04 /
      data C( 375)   /  -0.474165011293D+04 /
      data C( 376)   /  -0.106383826047D+04 /
      data C( 377)   /   0.212555444736D+04 /
      data C( 378)   /  -0.236013258835D+04 /
      data C( 379)   /  -0.301482092993D+04 /
      data C( 380)   /   0.264793530210D+02 /
      data C( 381)   /  -0.205207581287D+04 /
      data C( 382)   /  -0.173304927521D+04 /
      data C( 383)   /  -0.581724015357D+03 /
      data C( 384)   /   0.736106916388D+03 /
      data C( 385)   /  -0.862989104295D+03 /
      data C( 386)   /  -0.240480120333D+04 /
      data C( 387)   /   0.344863969249D+04 /
      data C( 388)   /   0.673930809805D+03 /
      data C( 389)   /   0.538758861001D+03 /
      data C( 390)   /   0.147967846147D+03 /
      data C( 391)   /  -0.189845771577D+04 /
      data C( 392)   /  -0.343381266242D+03 /
      data C( 393)   /   0.226151971199D+04 /
      data C( 394)   /  -0.348165990383D+04 /
      data C( 395)   /  -0.137886410095D+04 /
      data C( 396)   /  -0.452426797538D+04 /
      data C( 397)   /  -0.546311644413D+02 /
      data C( 398)   /  -0.480113520257D+02 /
      data C( 399)   /  -0.260536444352D+03 /
      data C( 400)   /  -0.133327190832D+04 /
      data C( 401)   /   0.632765868126D+03 /
      data C( 402)   /  -0.661117799374D+03 /
      data C( 403)   /   0.172506945261D+04 /
      data C( 404)   /  -0.615072273819D+03 /
      data C( 405)   /  -0.383815665429D+03 /
      data C( 406)   /   0.169025584445D+04 /
      data C( 407)   /  -0.808965608264D+03 /
      data C( 408)   /   0.593728899296D+03 /
      data C( 409)   /   0.119569641767D+04 /
      data C( 410)   /   0.430656118447D+03 /
      data C( 411)   /   0.587203407297D+04 /
      data C( 412)   /  -0.773418173876D+03 /
      data C( 413)   /   0.204092911587D+03 /
      data C( 414)   /   0.198285065744D+04 /
      data C( 415)   /  -0.145843638158D+04 /
      data C( 416)   /   0.154180304340D+04 /
      data C( 417)   /  -0.113615015119D+03 /
      data C( 418)   /   0.341159952346D+03 /
      data C( 419)   /   0.128014772782D+03 /
      data C( 420)   /  -0.241151083483D+04 /
      data C( 421)   /  -0.888711069985D+03 /
      data C( 422)   /  -0.255786858400D+03 /
      data C( 423)   /  -0.516753605238D+03 /
      data C( 424)   /   0.952501144399D+03 /
      data C( 425)   /   0.734163721480D+03 /
      data C( 426)   /  -0.103271529318D+04 /
      data C( 427)   /   0.212633590214D+03 /
      data C( 428)   /   0.303411042713D+03 /
      data C( 429)   /   0.409859079006D+03 /
      data C( 430)   /  -0.144512479307D+03 / 
      end 
	  
