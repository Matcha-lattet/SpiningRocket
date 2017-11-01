program Simulation

use atmosphere1976 
implicit none

!! ←正しい値未記入、もしくは計算が未実装なこと表す

real(8),parameter :: pi = 3.1415d0

integer i

!----------------------------------
!- 諸元                           -
!----------------------------------
! 機体諸元
real(8) :: m, m_dot                         ! 全重量[kg]
real(8) :: ms = 8.066d0                     ! 機体空虚質量[kg]
 
real(8) :: l = 1.627                        ! 全長[m]
real(8) :: Thrust                           ! 推力[N], ミスアライメントはないものとする
real(8) :: Clp = -0.126d0                   ! 横揺れ減衰モーメント(roll)
real(8) :: Cmq = -2.21d0                    ! 縦揺減衰モーメント(Pitch)
real(8) :: Cnr = -2.21d0                    ! 偏揺れ減衰モーメント(Yaw)
real(8) :: Cd = 0.55d0                      ! 抗力係数 !!- Mach数に対して変動しない -!!
real(8) :: CNa = 7.04d0                     ! 法線力係数[1/rad]
real(8) :: CYb = 7.04d0                     ! 横力係数
real(8) :: Isp = 166.0d0                      ! 比推力[sec]
real(8) :: burn_time = 1.641d0              ! 全燃焼秒時[sec]
real(8) :: Ae = 25.70d0 * 10**(-3)          ! ノズル面積 
real(8) :: lcp = 0.93d0                     ! 空力圧力中心 [m] 
real(8) :: lcgs = 0.897d0                   ! 機体空虚重心位置 [m]
real(8) :: lcg = 0.897d0
real(8) :: S = 0.01863d0                    ! 機体断面積 [m^2]
real(8) :: Ip0 !
real(8) :: lcgox,lcgox0 = 0.441d0           ! 酸化剤重心位置[m](from エンカバ) 
real(8) :: lcgf,lcgf0 = 0.119d0             ! 燃料重心位置[m](from エンカバ)
real(8) :: lcgp                             ! 推進剤重心位置[m]
real(8) :: mp0 
real(8) :: Mfm   !!                         ! フィンミスアライメントによるモーメント
real(8) :: freq = 3000.0d0
real(8) :: mox,mox0 = 0.389d0               ! 酸化剤質量[kg]
real(8) :: mf,mf0 = 0.176 
real(8) :: mf_a = 0.176d0                   ! 全燃料質量[kg]
real(8) :: mf_b = 0.145d0                   ! 燃焼後燃料質量[kg]
real(8) :: mox_dot = 0.238d0                ! 酸化剤質量流量[kg/sec]
real(8) :: mf_dot = (0.176-0.145)/1.64      ! 燃料質量流量[kg/sec] 
real(8) :: Ltank = 0.3d0                    !タンク長さ[m]
real(8) :: Ib(3)
real(8) :: df1=0.044,df2=0.024,ifp,ifr,Is=1.8,lf=0.099,Ir=0.022
real(8) :: CdS = 0.633
! 環境諸元
real(8) :: Ta,Ta0 = 25.0d0
real(8) :: rho, rho0 = 1.225d0 !!
real(8) :: g(3), g0 = 9.80665d0
real(8) :: Gravity(3)
real(8) :: Pa, Pa0 = 101325.0d0             ! 大気圧[Pa]  
real(8) :: Re = 6.371d0*10**3               ! 地球半径
real(8) :: Wh = 4.5d0                       ! 高度分布係数
real(8) :: Hw = 5.0d0                       ! 風速測定点
integer :: Vw_abs                           ! 風速
integer :: Vw_angle                         ! 風向
real(8) :: LL = 5.0d0                       ! ランチャ長[m]

real(8) :: theta,theta0 = -81.0d0 ![deg]
real(8) :: phi,phi0 = 0.0d0  ![deg]
real(8) :: psi,psi0 = 235.0d0  ![deg]

real(8) :: roll = 0
real(8) :: yaw = 0
real(8) :: pitch=0
real(8) :: Mx 
real(8) :: My 
real(8) :: Mz 
real(8) :: Ip !!

logical :: isLiftOff = .FALSE.
logical :: isThrust = .TRUE.
logical :: isLanchClear = .FALSE.

!加速度（機体座標）
real(8) :: accb(3) = 0
real(8) :: acce(3) = 0

! Force
real(8) :: D 
real(8) :: Y
real(8) :: N
real(8) :: Mjp 
real(8) :: Mjq
real(8) :: Mjr
real(8) :: Mp=0
real(8) :: Mq=0
real(8) :: Mr=0
real(8) :: p=0
real(8) :: q=0
real(8) :: r=0
real(8) :: p_dot
real(8) :: q_dot
real(8) :: r_dot

real(8) :: Vwe(3) = 0
real(8) :: Vwb(3) = 0
real(8) :: Vab(3) = 0

! 地上座標
real(8) :: Position(3) = 0

!機体対地速度
real(8) :: Ve(3) = 0

!機体速度（機体座標）
real(8) :: Vb(3) = 0

! 対気速度
real(8) :: Va_abs = 0
real(8) :: Va(3) = 0 

real(8) :: alpha = 0 
real(8) :: beta = 0


real(8) :: Abe(3,3) ! 座標変換行列
real(8) :: Aeb(3,3)
real(8) :: t = 0
real(8) :: dt


real(8) :: topPosition(3)
real(8) :: toptime
real(8) :: topVe(3)


dt = 1/freq

!open(998,file="debug.dat",status="replace")
open(999,file="debug.csv",status="replace")
write(999,*) "t",',', "theta*180/3.1415",',',"phi*180/3.1415",',',"psi*180/3.1415",',',"accb(1)",',',"accb(2)",',',"accb(3)",',',&
"Ve(3)",',',"Position(1)",',',"Position(2)",',',"Position(3)",',',"D",',',"Y",',',"N",',',"Va_abs",',',"alpha",',',"beta",',',&
"m",',',"mf",',',"mox",',',"Ib(1)",',',"Ib(2)",',',"Ib(3)",',',"Vb(1)",',',"Vb(2)",',',"Vb(3)"

!------------------------------------------------------------
!- Main Calculation                                         -
!------------------------------------------------------------
open(300,file="randing.csv")
print *,"風速,風向"
do Vw_abs = 1,1,1
  do Vw_angle = 0,0,45

call Initialize()
open(201,file="thrust.dat",status="old")
do 

  call Environment()
  call CalcMass()
  call CalcLcg()
  if(t < burn_time) then
    read(201,*) thrust
  else 
    Thrust = 0
    isThrust = .FALSE.
  end if
  
  if(sqrt(Position(1)**2+Position(2)**2+Position(3)**2)>LL) then
    isLanchClear = .TRUE.
  end if

  call SetAeb()
  Vwb = matmul(Aeb,Vwe)

  Vab = Vb - Vwb

  Va_abs = sqrt(Vab(1)**2+Vab(2)**2+Vab(3)**2)
  
  if(isLanchClear .eqv. .TRUE.) then
    call CalcMoment()
    Ip = ((m_dot/mp0)*Ip0 - (ms/m)**2*(lcgp-lcgs)**2*m_dot)*dt
    call Rotation()
    alpha = atan(Vab(3)/Vab(1))
    beta = asin(-Vab(2)/Va_abs)
  end if

  call SetAeb()
  Gravity = matmul(Aeb,g)

  D = 0.5*S*rho*Va_abs**2*Cd
  Y = 0.5*S*rho*Va_abs**2*CNa*beta
  N = 0.5*S*rho*Va_abs**2*CNa*alpha
  accb(1) = (Thrust-D)/m + (r*Vb(2) - q*Vb(3)) + Gravity(1)
  accb(2) = -Y/m + (p*Vb(3) - r*Vb(1)) + Gravity(2)
  accb(3) = N/m + (q*Vb(1) - p*Vb(2)) + Gravity(3)

  Vb(1) = Vb(1) + accb(1)*dt
  Vb(2) = Vb(2) + accb(2)*dt
  Vb(3) = Vb(3) + accb(3)*dt

  call SetAbe()
  Ve = matmul(Abe,Vb)

  


  Position(1) = Position(1) + Ve(1)*dt
  Position(2) = Position(2) + Ve(2)*dt
  Position(3) = Position(3) + Ve(3)*dt

  if(isLiftOff .eqv. .FALSE.) then
    if(Position(3) < 0) then
      call Initialize()
    else
      isLiftOff = .TRUE.
    end if
  end if

  if(Position(3)>topPosition(3)) then
    topPosition = Position
    toptime = t
    topVe = Ve
  end if 

  

 write(999,*) t, ',', theta*180/3.1415,',',phi*180/3.1415,',',psi*180/3.1415,',',accb(1),',',accb(2),',',accb(3) &
 ,',',Ve(3),',',Position(1),',',Position(2),',',Position(3),',',D,',',Y,',',N,',',Va_abs,',',&
 alpha*180/pi,',',beta*180/pi,',',&
 m,',',mf,',',mox,',',Ib(1),',',Ib(2),',',Ib(3),',',Vb(1),',',Vb(2),',',Vb(3)

  t = t + dt
  if(Position(3) < 0) exit

end do 


close(201) !Thrust.dat
write(300, '(2(f0.3,a)$)') Position(1),',',Position(2),','

!parachute
Position = topPosition
t = toptime
Ve = topVe
do
  call Environment()
  acce(3) = g(3) + 0.5*rho*CdS*Ve(3)**2/m
  Ve(1) = Vwe(1)
  Ve(2) = Vwe(2)
  Ve(3) = Ve(3) + acce(3)*dt
  Position = Position + Ve*dt

  if(Position(3)<0) exit
end do
write(300,*) Position(1),',',Position(2)


print *, Vw_abs,Vw_angle,topPosition(3)
end do 
end do

close(300)
close(999)
!close(998)
 
contains

subroutine SetAbe()
  Abe(1,1) = cos(theta)*cos(psi)
  Abe(2,1) = cos(theta)*sin(psi)
  Abe(3,1) = -sin(theta)
  Abe(1,2) = sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)
  Abe(2,2) = sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)
  Abe(3,2) = sin(phi)*cos(theta)
  Abe(1,3) = cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)
  Abe(2,3) = cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi)
  Abe(3,3) = cos(phi)*cos(theta)
end subroutine SetAbe


subroutine SetAeb()
  Aeb(1,1) = cos(theta)*cos(psi) 
  Aeb(2,1) = sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)
  Aeb(3,1) = cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)
  Aeb(1,2) = cos(theta)*sin(psi)
  Aeb(2,2) = sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)
  Aeb(3,2) = cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi)
  Aeb(1,3) = -sin(theta)
  Aeb(2,3) = sin(phi)*cos(theta)
  Aeb(3,3) = cos(phi)*cos(theta)
end subroutine SetAeb

subroutine CalcMass()
  if(isThrust .eqv. .TRUE.) then
    mox = mox - mox_dot*dt
    mf = mf - mf_dot*dt
  else
    mox = 0.0d0
    mf = mf_b
  end if
  mp = mox + mf
  m = ms + mp
end subroutine CalcMass


subroutine CalcLcg()
  Lcgox = Lcgox0 - 0.5d0 * (1.0d0 - (mox / mox0)) * Ltank
  ! if (mox <= 0.0d0) lcgox = 0.0d0　この式は間違えてるかも
  Lcgp = (mf * Lcgf + mox * Lcgox) / (mf + mox)
  Lcg = ((mp * Lcgp) + (ms * Lcgs)) / m
end subroutine CalcLcg

subroutine CalcMoment()

real :: Ka(3)

  ! !- ジェットダンピングモーメント -----------------------------
  ! re_dot = ("ノズル半径")**2/2
  ! Mjp = -(m_dot*k**2 + mass*k_dot_bar**2 - mass*re**2)*r
  ! Mjq = -(Ip_dot + m_dot*((lcg - lcgp)**2 - re**2))*r
  ! Mjr = -(Ip_dot + m_dot*((lcg - lcgp)**2 - re**2))*q
  Mjp = 0
  Mjq = 0
  Mjr = 0

  !- モーメント計算 --------------------------------------------------------
  ! 第一項：空気力による加速側のモーメント
  ! 第二項：空気力による減衰モーメント
  ! 第三項：ジェットダンピングモーメント
  ! ※ 推力のミスアライメント、姿勢制御用補助ロケットエンジンなどはないものとする
  !----------------------------------------------------------------------
  Mfm = 0 !フィン取り付け角はとりあえず0deg
  Mp = Mfm + 0.5*rho*S*Va_abs**2*(d**2/0.5/Va_abs)*Clp*p + Mjp
  Mq = -(lcp - lcg)*N + 0.5*rho*S*Va_abs**2*(l**2/0.5/Va_abs)*Cmq*q + Mjq
  Mr = -(lcp - lcg)*Y + 0.5*rho*S*Va_abs**2*(l**2/0.5/Va_abs)*Cnr*r + Mjr



  Ka(1) = 0.5*rho*S*Va_abs**2*(d**2/0.5/Va_abs)*Clp*p
  Ka(2) = 0.5*rho*S*Va_abs**2*(l**2/0.5/Va_abs)*Cmq*q
  Ka(3) = 0.5*rho*S*Va_abs**2*(l**2/0.5/Va_abs)*Cnr*r
  print *,Ka(1),Ka(2),Ka(3)

  Ifp = mf * lf / 12.0d0
  Ifr = mf * (((df1 + df2) / 4.0d0) + ((df1 + df2) / 4.0d0) * (1.0d0 - (mf / mf_b)))
  
  Ib(2) = Ifp + mf * (lcg - lcgf)**2 + Is * ((ms + mox) / ms) + (ms + mox) * (lcg - lcgs)**2
  Ib(3) = Ib(2)
  Ib(1) = Ir + Ifr




end subroutine CalcMoment


subroutine Rotation()
  if(sqrt(Position(1)**2+Position(2)**2+Position(3)**2)<LL) then
    p_dot = 0.0d0
    q_dot = 0.0d0
    r_dot = 0.0d0   
  else
    p_dot = Mp/Ib(1)
    q_dot = ((Ib(3) - Ib(1))*r*p+Mq)/Ib(2)
    r_dot = ((Ib(1) - Ib(2))*p*q+Mr)/Ib(3)
    
  end if

  p = p + p_dot*dt
  q = q + q_dot*dt
  r = r + r_dot*dt

  if(89<theta*180/pi .and. theta*180/pi<91) then
    phi = phi + (0.5*p - (q_dot*sin(phi)+r_dot*cos(phi))/(2*(q*cos(phi) - r*sin(phi))))*dt
    theta = theta + (q*cos(phi) - r*sin(phi))*dt
    psi = psi  + (-0.5*p - (q_dot*sin(phi)+r_dot*cos(phi))/(2*(q*cos(phi) - r*sin(phi))))*dt
  else 
    phi = phi + (p + tan(theta)*(q*sin(phi)+r*cos(phi)))*dt
    theta = theta + (q*cos(phi) - r*sin(phi))*dt
    psi = psi + ((q*sin(phi) + r*cos(phi))/cos(theta))*dt
  end if
end subroutine Rotation 

subroutine Initialize()
  
  
  lcgox = lcgox0                ! 酸化剤重心位置[m](from エンカバ) 
  lcgf = lcgf0                  ! 燃料重心位置[m](from エンカバ)
  mox = mox0                    ! 酸化剤質量[kg]
  mf = mf0 
  
  
  theta = theta0*pi/180 !![deg]
  phi = phi0*pi/180 !! [deg]
  psi = psi0*pi/180 !! [deg]

  roll = 0.0d0!!
  yaw = 0.0d0!!
  pitch = 0.0d0!!

  
  isLiftOff = .FALSE.
  isThrust = .TRUE.
  isLanchClear = .FALSE.

  rho = rho0
  g = g0
  Pa = Pa0
  
  alpha = 0 
  beta = 0
  p=0
  q=0
  r=0

  Position = 0
  Ve = 0
  Vb = 0

  t = 0

  topPosition = 0
  toptime = 0
  topVe = 0

end subroutine Initialize

subroutine Environment()
  call Atmosphere(Position(3)/1000, rho, Pa, Ta)
  Ta = Ta*Ta0
  Pa = Pa*Pa0
  rho = rho*rho0

  g(1) = 0
  g(2) = 0
  g(3) = -g0 * (Re / (Re + Position(3)))**2 ! あるいは
  ! g(3)=G*m/(Re+Position(3))**2

  Vwe(1) = Vw_abs*cos(Vw_angle*pi/180)*(Position(3)/Hw)**(1/Wh)
  Vwe(2) = Vw_abs*sin(Vw_angle*pi/180)*(Position(3)/Hw)**(1/Wh)
  Vwe(3) = 0

end subroutine Environment

end program Simulation

