module sampling

use mt95

implicit none
private

public algorithm_A, algorithm_D, algorithm_S
contains

subroutine algorithm_A(N_in,m_in,Sample)

integer, intent(in) :: N_in, m_in
double precision :: V, q
integer :: N, m, i, t, S

integer :: P, J
integer, dimension(:), intent(inout) :: Sample

N=N_in
m=m_in
J=0
P=1

do i = m, 2, -1
    t=N-i
!    print*, "i is", i, "t is", t, "N is", N
    call genrand_real2(V)
    S = 0
    q = dble(t)/dble(N)
!    print*, "V", V, "quot", q
    do while (q > V)
        S = S + 1
        t = t - 1
        N = N - 1
        q = q * t / N
!    print*, "V", V, "quot", q
    enddo
    P=P+S
    J = J + 1
    Sample(J)=P
!    print*, "skip", S, "selected", P
!    print*, P
    P=P+1
    N=N-1
enddo
!print*, "now i is", i
if (i==1) then
    call genrand_real2(V)
    S = int(N*V)
    P=P+S
    J = J + 1
    Sample(J)=P
!    print*, "skip", S, "selected", P
!    print*, P
endif

end subroutine algorithm_A

subroutine algorithm_D(N_in,m_in,Sample)
integer, intent(in) :: N_in, m_in
integer, parameter :: negalphainv=-13
integer :: N, m, S
double precision :: Nreal, mreal, minv, mmin1inv
double precision :: Vrand, Vprime, U
integer :: qu1, threshold, t, limit
double precision :: qu1real, X, negSreal, y1, y2
double precision :: bottom, top

integer, dimension(:), intent(inout) :: Sample 
integer :: P, J

N=N_in
m=m_in
J=0
P=1

mreal=dble(m); minv=1.d0/mreal
Nreal=dble(N)
call genrand_real2(Vrand)
Vprime = Vrand**minv
qu1 = -m + 1 + N
qu1real = -mreal + 1.d0 + Nreal
threshold = -negalphainv * m

do while (m > 1 .and. threshold < N)
    mmin1inv = 1.d0/(-1.d0+mreal)
    do
        do 
            X = Nreal*(-Vprime+1.d0)
            S = int(X)
            if (S < qu1) exit
            call genrand_real2(Vrand)
            Vprime = Vrand**minv
        enddo
        call genrand_real2(U)
        negSreal = -dble(S)
        y1 = (U*Nreal/qu1real)**mmin1inv
        Vprime = y1 * (-X/Nreal+1.d0) * (qu1real/(negSreal+qu1real))
        if (Vprime <= 1.d0) exit
        y2=1.d0
        top=-1.d0 + Nreal
        if (-1+n>S) then
            bottom = -mreal + Nreal 
            limit = -S + N
        else
            bottom = -1.d0 + negSreal + Nreal
            limit = qu1
        endif

        do t = -1+N, limit, -1
            y2 = (y2 * top) / bottom
            top = -1.d0 + top 
            bottom = -1.d0 + bottom
        enddo

        if (Nreal/(-X+Nreal) >= y2**mmin1inv) then
            call genrand_real2(Vrand)
            Vprime = Vrand**mmin1inv
            exit
        endif
        call genrand_real2(Vrand)
        Vprime = Vrand**minv
    enddo
    P = P + S
    J = J + 1
    Sample(J)=P
!    print*, "skip", S, "selected", P
!    print*, P
    P=P+1
    N = -S + (-1+N); Nreal = negSreal + (-1.d0 + Nreal)
    m = -1 + m; mreal = -1.d0 + mreal; minv = mmin1inv
    qu1 = -S + qu1; qu1real = negSreal + qu1real
    threshold = threshold + negalphainv
enddo

if (m > 1) then
    call algorithm_A_in_D
else
    S = int(N*Vprime)
    P = P + S
    J = J + 1
    Sample(J)=P
!    print*, "skip", S, "selected", P
!    print*, P
endif

contains

subroutine algorithm_A_in_D

double precision :: V, q
integer :: i, t

do i = m, 2, -1
    t=N-i
!    print*, "i is", i, "t is", t, "N is", N
    call genrand_real2(V)
    S = 0
    q = dble(t)/dble(N)
!    print*, "V", V, "quot", q
    do while (q > V)
        S = S + 1
        t = t - 1
        N = N - 1
        q = q * t / N
!    print*, "V", V, "quot", q
    enddo
    P=P+S
    J = J + 1
    Sample(J)=P
!    print*, "skip", S, "selected", P
!    print*, P
    P=P+1
    N=N-1
enddo
!print*, "now i is", i
if (i==1) then
    call genrand_real2(V)
    S = int(N*V)
    P=P+S
    J = J + 1
    Sample(J)=P
!    print*, "skip", S, "selected", P
!    print*, P
endif

end subroutine algorithm_A_in_D

end subroutine algorithm_D

subroutine algorithm_S(N_in,m_in,Sample)

integer, intent(in) :: N_in, m_in

integer :: n_select, s, m, l
integer :: b_select
double precision :: k

integer, dimension(:), intent(inout) :: Sample

n_select = m_in
b_select = N_in

m=0
do s=1, b_select
    call genrand_real2(k)
    l=int((b_select-s+1)*k)+1
    if (l > (n_select-m)) then
        cycle
    endif
    m=m+1

!    bond_list(m)=s
    Sample(m)=s
!    print*, s
    if (m>=n_select) then
        exit
    endif
enddo

end subroutine algorithm_S

end module
