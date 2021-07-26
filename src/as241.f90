module quantiles
use iso_c_binding
implicit none

contains

real (C_DOUBLE) function ppnd16(P,Ifault) bind(c)
    implicit none

    real (C_DOUBLE), intent(in) :: P
    Integer (C_INT), intent(out) :: Ifault
    real ZERO , ONE , HALF , SPLIT1 , SPLIT2 , CONST1 ,   &
                     & CONST2 , A0 , A1 , A2 , A3 , A4 , A5 , A6 , A7 , &
                     & B1 , B2 , B3 , B4 , B5 , B6 , B7 , C0 , C1 , C2 ,&
                     & C3 , C4 , C5 , C6 , C7 , D1 , D2 , D3 , D4 , D5 ,&
                     & D6 , D7 , E0 , E1 , E2 , E3 , E4 , E5 , E6 , E7 ,&
                     & F1 , F2 , F3 , F4 , F5 , F6 , F7 , q , r
    PARAMETER (ZERO=0.D0,ONE=1.D0,HALF=0.5D0,SPLIT1=0.425D0,          &
               & SPLIT2=5.D0,CONST1=0.180625D0,CONST2=1.6D0)
!
!	Coefficients for P close to 0.5
!
    PARAMETER (A0=3.3871328727963666080D0,A1=1.3314166789178437745D+2,&
               & A2=1.9715909503065514427D+3,                           &
               & A3=1.3731693765509461125D+4,                           &
               & A4=4.5921953931549871457D+4,                           &
               & A5=6.7265770927008700853D+4,                           &
               & A6=3.3430575583588128105D+4,                           &
               & A7=2.5090809287301226727D+3,                           &
               & B1=4.2313330701600911252D+1,                           &
               & B2=6.8718700749205790830D+2,                           &
               & B3=5.3941960214247511077D+3,                           &
               & B4=2.1213794301586595867D+4,                           &
               & B5=3.9307895800092710610D+4,                           &
               & B6=2.8729085735721942674D+4,                           &
               & B7=5.2264952788528545610D+3)
!	HASH SUM AB    55.88319 28806 14901 4439
!
!	Coefficients for P not close to 0, 0.5 or 1.
!
    PARAMETER (C0=1.42343711074968357734D0,                           &
               & C1=4.63033784615654529590D0,                           &
               & C2=5.76949722146069140550D0,                           &
               & C3=3.64784832476320460504D0,                           &
               & C4=1.27045825245236838258D0,                           &
               & C5=2.41780725177450611770D-1,                          &
               & C6=2.27238449892691845833D-2,                          &
               & C7=7.74545014278341407640D-4,                          &
               & D1=2.05319162663775882187D0,                           &
               & D2=1.67638483018380384940D0,                           &
               & D3=6.89767334985100004550D-1,                          &
               & D4=1.48103976427480074590D-1,                          &
               & D5=1.51986665636164571966D-2,                          &
               & D6=5.47593808499534494600D-4,                          &
               & D7=1.05075007164441684324D-9)
!	HASH SUM CD    49.33206 50330 16102 89036
!
!	Coefficients for P near 0 or 1.
!
    PARAMETER (E0=6.65790464350110377720D0,                           &
               & E1=5.46378491116411436990D0,                           &
               & E2=1.78482653991729133580D0,                           &
               & E3=2.96560571828504891230D-1,                          &
               & E4=2.65321895265761230930D-2,                          &
               & E5=1.24266094738807843860D-3,                          &
               & E6=2.71155556874348757815D-5,                          &
               & E7=2.01033439929228813265D-7,                          &
               & F1=5.99832206555887937690D-1,                          &
               & F2=1.36929880922735805310D-1,                          &
               & F3=1.48753612908506148525D-2,                          &
               & F4=7.86869131145613259100D-4,                          &
               & F5=1.84631831751005468180D-5,                          &
               & F6=1.42151175831644588870D-7,                          &
               & F7=2.04426310338993978564D-15)
!	HASH SUM EF    47.52583 31754 92896 71629
!
    Ifault = 0
    q = P - HALF
    IF ( ABS(q).LE.SPLIT1 ) THEN
        r = CONST1 - q*q
        ppnd16 = q*(((((((A7*r+A6)*r+A5)*r+A4)*r+A3)*r+A2)*r+A1)*r+A0) &
                & /(((((((B7*r+B6)*r+B5)*r+B4)*r+B3)*r+B2)*r+B1)*r+ONE)
        RETURN
    ELSE
        IF ( q.LT.ZERO ) THEN
            r = P
        ELSE
            r = ONE - P
        ENDIF
        IF ( r.LE.ZERO ) THEN
            Ifault = 1
            ppnd16 = ZERO
            RETURN
        ENDIF
        r = SQRT(-LOG(r))
        IF ( r.LE.SPLIT2 ) THEN
            r = r - CONST2
            ppnd16 = (((((((C7*r+C6)*r+C5)*r+C4)*r+C3)*r+C2)*r+C1)*r+C0)&
                   & /(((((((D7*r+D6)*r+D5)*r+D4)*r+D3)*r+D2)*r+D1)     &
                   & *r+ONE)
        ELSE
            r = r - SPLIT2
            ppnd16 = (((((((E7*r+E6)*r+E5)*r+E4)*r+E3)*r+E2)*r+E1)*r+E0)&
                   & /(((((((F7*r+F6)*r+F5)*r+F4)*r+F3)*r+F2)*r+F1)     &
                   & *r+ONE)
        ENDIF
        IF ( q.LT.ZERO ) ppnd16 = -ppnd16
        RETURN
    ENDIF
end function ppnd16
end module quantiles
