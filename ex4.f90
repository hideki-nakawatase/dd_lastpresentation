!- 3次元密度一様モデル
!-   Stommel (1948) の設定で実験する場合 (他はデフォルト設定)
!-   im = 100, jm = 62, km = 1
!-   dx = 100km, dy = 100km, dz = 200m
!-   fric_non = 0.02 / dz
program main
  implicit none

  character(10):: &
       expid='ex4-case1' ! 実験名

  integer,parameter:: &
       im = 100, & ! x方向格子数
       jm =  62, & ! y方向格子数
       km =  12   ! z方向格子数

  real(8),parameter:: &
       t0_npm2  = -0.1d0, & ! 東西方向風応力
       gr_mps2  = 9.80d0, & ! 重力加速度
       r0_kgpm3 = 1.d3, & ! 基準海水密度
       f0_psec  = 0.d-4, & ! コリオリ係数
       bt_pmsec = 1.d-11, & ! ベータ
       er_m     = 6370.d3,& ! 地球半径
       om_psec  = 7.29d-5,& ! 地球自転速度
       ah_m2ps  = 1.d5, & ! 水平渦粘性係数
       av_m2ps  = 1.d-3, & ! 鉛直渦粘性係数
       kh_m2ps  = 5.d3,   & ! 水平渦拡散係数
       kv_m2ps  = 1.d-4,  & ! 鉛直渦拡散係数
       gm_psec  = 1.d0/50.d0/86400.d0, & ! 温位緩和係数
       t0_c     = 10.d0,  & ! 基準温位
       s0_psu   = 34.d0,  & ! 基準塩分
       alphat_pk = 1.65d-4,  & ! 熱膨張率
       alphas_ppsu = 7.61d-4,  & ! 塩収縮率
       dx_m     = 200.d3, & ! x方向格子幅
       dy_m     = 100.d3, & ! y方向格子幅
       dz_m     = 50.d0, & ! z方向格子幅
       dt_sec   = 100.d0, & ! 時間刻み幅
       time_to_start_sec   = 0.d0, & ! 実験開始時刻
       time_to_end_sec     = 86400.d0*30.d0*12.d0, &   ! 実験終了時刻（１年）
       !time_to_end_sec     = dt_sec * 10.d0, & ! 実験終了時刻（１年）
       output_interval_sec = 86400.d0*10.d0,       & ! 出力時間間隔（１０日）
       !output_interval_sec = dt_sec , & ! 出力時間間隔（１０日）
       asf      = 0.5d0 ! アセリンフィルター係数
  real(8),parameter :: fric_non = 0.0d0 / dz_m  !- linear friction coefficient (for Stommel 1948)
  !real(8),parameter :: fric_non = 0.02d0 / dz_m  !- linear friction coefficient (for Stommel 1948) !摩擦あり

  real(8),dimension(:,:,:):: &
       ua_mps(0:im+1,0:jm+1,0:km+1), & ! 東西流速
       ub_mps(0:im+1,0:jm+1,0:km+1), &
       uc_mps(0:im+1,0:jm+1,0:km+1), &
       va_mps(0:im+1,0:jm+1,0:km+1), & ! 南北流速
       vb_mps(0:im+1,0:jm+1,0:km+1), &
       vc_mps(0:im+1,0:jm+1,0:km+1), &
       ta_c(0:im+1,0:jm+1,0:km+1),   & ! 温位
       tb_c(0:im+1,0:jm+1,0:km+1),   &
       tc_c(0:im+1,0:jm+1,0:km+1),   &
       sa_psu(0:im+1,0:jm+1,0:km+1), & ! 塩分 （34.0からの偏差）
       sb_psu(0:im+1,0:jm+1,0:km+1), &
       sc_psu(0:im+1,0:jm+1,0:km+1)
  
  real(8),dimension(:,:):: & 
       ea_m(0:im+1,0:jm+1), & ! 水位
       eb_m(0:im+1,0:jm+1), &
       ec_m(0:im+1,0:jm+1)

  real(8),dimension(:,:,:):: &
       ww_mps(0:im+1,0:jm+1,0:km+1), & ! 鉛直流速
       pp_npm2(0:im+1,0:jm+1,0:km+1), &  ! 圧力
       rr_kgpm3(0:im+1,0:jm+1,0:km+1)   ! 密度

  real(8),dimension(:,:):: &
       tx_npm2(0:im+1,0:jm+1), & ! 風応力（東西）
       ty_npm2(0:im+1,0:jm+1), & ! 風応力（南北）
       at_c(0:im+1,0:jm+1),      & ! 緩和温位（気温）
       sf_psumps(0:im+1,0:jm+1)    ! 海面塩分フラックス [psu m/s]

  real(8),dimension(:):: &
       fs_psec(0:jm+1) ! コリオリ係数

  real(8):: &
       time_sec, & ! 時間
       time_to_output_sec, & ! 次期出力時刻
       dpi, & ! π
       gu, gv, ge, gt, gs, &
       pre, adx, ady, adz, cor, dfx, dfy, dfz, frc, sum, kx, ky, tav, sav
  
  integer:: &
       i,j,k,n ! 作業用

  character(80):: &
       buff ! 作業用

  ! コリオリ係数
  do j = 0, jm+1
   	fs_psec(j) = f0_psec + bt_pmsec*dy_m*(dble(j-jm/2)-0.5d0)
     !fs_psec(j) = 2.d0 * om_psec * sin( (dble(j)-0.5d0) * dy_m /er_m )
  end do

  ! 強制力
	! /10忘れずに
  dpi=4.d0*atan(1.d0)
  do j=1,jm
     do i=1,im
        tx_npm2(i,j) = t0_npm2 * cos( dpi*( (dble(j-jm/2)-0.5d0) / dble(jm/2) ) ) !亜熱帯循環
        ty_npm2(i,j) = t0_npm2 * sin( dpi/2 * (dble(j-jm/2)-0.5d0) / dble(jm/2) )
				! 海面温度は真ん中で30℃,端で２５℃になるように設定
        at_c(i,j)      = 5.d0 *cos( dpi*( dble(j-jm/2)-0.5d0 ) / dble(jm/2)) + 25.d0 ! 東西一様水温
        sf_psumps(i,j) = 0.d0
     end do
  end do

  ! 初期化
  uc_mps(:,:,:) = 0.d0
  vc_mps(:,:,:) = 0.d0
  ec_m(:,:)     = 0.d0

  ww_mps(:,:,:)  = 0.d0
  pp_npm2(:,:,:) = 0.d0

  ! 初期設定
  if( time_to_start_sec == 0.d0 ) then ! 初期化
    ua_mps(:,:,:) = 0 
    ub_mps(:,:,:) = 0 
    va_mps(:,:,:) = 0
    vb_mps(:,:,:) = 0
    ea_m(:,:)     = 0
    eb_m(:,:)     = 0
		!　東西方向に温度勾配をつける（エルニーニョ）
		do j=1,jm
			do k=1,km
				ta_c(:,j,k)   = (5.d0 *cos( dpi*( dble(j-jm/2)-0.5d0 ) / dble(jm/2)) + 25.d0)*((k**2)/km/km)
				tb_c(:,j,k)   = (5.d0 *cos( dpi*( dble(j-jm/2)-0.5d0 ) / dble(jm/2)) + 25.d0)*((k**2)/km/km)
			end do
		end do
    sa_psu(:,:,:) = 0.d0 !- 34からの偏差
    sb_psu(:,:,:) = 0.d0
    time_sec      = 0

  else ! 継続

     open(10,file=trim(expid)//'.cnt',form='unformatted',access='stream')
     read(10) time_sec,ua_mps,ub_mps,va_mps,vb_mps,ea_m,eb_m
     close(10)

  endif

  time_to_output_sec = time_sec + output_interval_sec

  do
     if( time_sec > time_to_end_sec ) exit

     ! 診断量の計算
     ! r
     do k = km, 1, -1
        !pre_npm2 = pre_npm2 + gr_mps2 * r0_kgpm3 * dz_m !- 圧力はrrで計算すべきだがr0で近似
        do j = 1, jm
           do i = 1, im
              rr_kgpm3(i,j,k) = r0_kgpm3*(1 - alphat_pk*(tb_c(i,j,k) - t0_c) + alphas_ppsu*(sb_psu(i,j,k) - s0_psu))
           end do
        end do
     end do

     ! w
     do k = 1, km
        do j = 1, jm
           do i = 1, im
              ww_mps(i,j,k) = ww_mps(i,j,k-1) - ((ub_mps(i,j,k) - ub_mps(i-1,j,k))/dx_m + (vb_mps(i,j,k) - vb_mps(i,j-1,k))/dy_m)*dz_m
           end do
        end do
     end do
     ww_mps(:,:,km) = 0.d0

     ! p
     do j = 1, jm
        do i = 1, im
           pp_npm2(i,j,km) = gr_mps2*rr_kgpm3(i,j,km)*(dz_m*0.5d0 + eb_m(i,j))
           do k = km-1, 1, -1
              pp_npm2(i,j,k) = pp_npm2(i,j,k+1) + gr_mps2*(rr_kgpm3(i,j,k) + rr_kgpm3(i,j,k+1))*dz_m*0.5d0
           end do
        end do
     end do

     ! 予報量の計算
     ! u
     do k = 1, km

        do j = 1, jm
           do i = 1, im - 1
              adx = ((ub_mps(i,j,k)+ub_mps(i+1,j,k))**2-(ub_mps(i-1,j,k)+ub_mps(i,j,k))**2)*0.25d0/dx_m
              ady = ((vb_mps(i,j,k)+vb_mps(i+1,j,k))*(ub_mps(i,j,k)+ub_mps(i,j+1,k))-(vb_mps(i,j-1,k)+vb_mps(i+1,j-1,k))*(ub_mps(i,j-1,k)+ub_mps(i,j,k)))*0.25d0/dy_m
              adz = ((ub_mps(i,j,k)+ub_mps(i,j,k+1))*(ww_mps(i,j,k)+ww_mps(i+1,j,k)) - (ub_mps(i,j,k-1)+ub_mps(i,j,k))*(ww_mps(i,j,k-1)+ww_mps(i+1,j,k-1)))*0.25d0/dz_m
              cor = fs_psec(j)*(vb_mps(i,j,k)+vb_mps(i+1,j,k)+vb_mps(i,j-1,k)+vb_mps(i+1,j-1,k))*0.25d0
              pre = -1/rr_kgpm3(i,j,k)*(pp_npm2(i+1,j,k)-pp_npm2(i,j,k))/dx_m
              dfx = ah_m2ps*(ua_mps(i+1,j,k)-2.d0*ua_mps(i,j,k)+ua_mps(i-1,j,k))/(dx_m**2)
              dfy = ah_m2ps*(ua_mps(i,j+1,k)-2.d0*ua_mps(i,j,k)+ua_mps(i,j-1,k))/(dy_m**2)
              dfz = av_m2ps*(ua_mps(i,j,k+1)-2.d0*ua_mps(i,j,k)+ua_mps(i,j,k-1))/(dz_m**2)
              uc_mps(i,j,k) = ua_mps(i,j,k) + 2.d0*dt_sec*(cor+pre-adx-ady-adz+dfx+dfy+dfz)
              !uc_mps(i,j,k) = ua_mps(i,j,k) + 2.d0*dt_sec*(cor+pre) - fric_non*ub_mps(i,j,k) !摩擦あり
           end do
        end do
     end do
     
     ! v
     do k = 1, km
        do j = 1, jm - 1
           do i = 1, im
               adx = ((vb_mps(i,j,k)+vb_mps(i+1,j,k))*(ub_mps(i,j+1,k)+ub_mps(i,j,k))-(vb_mps(i-1,j,k)+vb_mps(i,j,k))*(ub_mps(i-1,j+1,k)+ub_mps(i-1,j,k)))*0.25d0/dx_m
               ady = ((vb_mps(i,j+1,k)+vb_mps(i,j,k))**2-(vb_mps(i,j,k)+vb_mps(i,j-1,k))**2)*0.25d0/dy_m
               adz = ((vb_mps(i,j,k+1)+vb_mps(i,j,k))*(ww_mps(i,j,k)+ww_mps(i,j+1,k))-(vb_mps(i,j,k)+vb_mps(i,j,k-1))*(ww_mps(i,j,k-1)+ww_mps(i,j+1,k-1)))*0.25d0/dz_m
               cor =  fs_psec(j)*(ub_mps(i-1,j,k)+ub_mps(i,j,k)+ub_mps(i-1,j+1,k)+ub_mps(i,j+1,k))*0.25d0
               pre = -1/rr_kgpm3(i,j,k)*(pp_npm2(i,j+1,k)-pp_npm2(i,j,k))/dy_m
               dfx = ah_m2ps*(va_mps(i+1,j,k)-2.d0*va_mps(i,j,k)+va_mps(i-1,j,k))/(dx_m**2)
               dfy = ah_m2ps*(va_mps(i,j+1,k)-2.d0*va_mps(i,j,k)+va_mps(i,j-1,k))/(dy_m**2)
               dfz = av_m2ps*(va_mps(i,j,k+1)-2.d0*va_mps(i,j,k)+va_mps(i,j,k-1))/(dz_m**2)
               vc_mps(i,j,k) = va_mps(i,j,k) + 2.d0*dt_sec*(-cor+pre-adx-ady-adz+dfx+dfy+dfz)
               ! vc_mps(i,j,k) = va_mps(i,j,k) + 2.d0*dt_sec*(-cor+pre) - fric_non*vb_mps(i,j,k) !摩擦あり
      
           end do
        end do
     end do

     ! e
     do j = 1, jm
        do i = 1, im
					sum = 0
					do k = 1, km
						sum = sum + (ub_mps(i,j,k)-ub_mps(i-1,j,k))*dz_m/dx_m + (vb_mps(i,j,k)-vb_mps(i,j-1,k))*dz_m/dy_m
					end do
					kx = kh_m2ps*(ea_m(i+1,j)-2*ea_m(i,j)+ea_m(i-1,j))/dx_m**2
					ky = kh_m2ps*(ea_m(i,j+1)-2*ea_m(i,j)+ea_m(i,j-1))/dy_m**2
					ec_m(i,j) = ea_m(i,j) - 2.d0*dt_sec*(sum + kx + ky)
        end do
     end do

		 ! t
		 do k = 1, km
        do j = 1, jm
           do i = 1, im
              adx=-((tb_c(i+1,j,k) + tb_c(i,j,k))*ub_mps(i,j,k) - (tb_c(i,j,k) + tb_c(i-1,j,k))*ub_mps(i-1,j,k))*0.5d0/dx_m
              ady=-((tb_c(i,j+1,k) + tb_c(i,j,k))*vb_mps(i,j,k) - (tb_c(i,j,k) + tb_c(i,j-1,k))*vb_mps(i,j-1,k))*0.5d0/dy_m
              adz=-((tb_c(i,j,k+1) + tb_c(i,j,k))*ww_mps(i,j,k) - (tb_c(i,j,k) + tb_c(i,j,k-1))*ww_mps(i,j,k-1))*0.5d0/dz_m
              dfx=kh_m2ps*(ta_c(i+1,j,k) - 2*ta_c(i,j,k) + ta_c(i-1,j,k))/dx_m/dx_m
              dfy=kh_m2ps*(ta_c(i,j+1,k) - 2*ta_c(i,j,k) + ta_c(i,j-1,k))/dy_m/dy_m
              dfz=kv_m2ps*(ta_c(i,j,k+1) - 2*ta_c(i,j,k) + ta_c(i,j,k-1))/dz_m/dz_m

              if ( k == km ) then
                frc = gm_psec * ( at_c(i,j) - tb_c(i,j,km) )
              else
                frc = 0.d0
              endif

              gt = adx+ady+adz+dfx+dfy+dfz+frc
              tc_c(i,j,k) = ta_c(i,j,k) + 2.d0 * dt_sec * gt
           end do
        end do
     end do

		 ! s
     do k = 1, km
        do j = 1, jm
           do i = 1, im
              adx=-((sb_psu(i+1,j,k)+sb_psu(i,j,k))*ub_mps(i,j,k) - (sb_psu(i,j,k) + sb_psu(i-1,j,k))*ub_mps(i-1,j,k))*0.5d0/dx_m
              ady=-((sb_psu(i,j+1,k)+sb_psu(i,j,k))*vb_mps(i,j,k) - (sb_psu(i,j,k) + sb_psu(i,j-1,k))*vb_mps(i,j-1,k))*0.5d0/dy_m
              adz=-((sb_psu(i,j,k+1)+sb_psu(i,j,k))*ww_mps(i,j,k) - (sb_psu(i,j,k) + sb_psu(i,j,k-1))*ww_mps(i,j,k-1))*0.5d0/dz_m
              dfx=kh_m2ps*(sa_psu(i+1,j,k) - 2*sa_psu(i,j,k) + sa_psu(i-1,j,k))/dx_m/dx_m
              dfy=kh_m2ps*(sa_psu(i,j+1,k) - 2*sa_psu(i,j,k) + sa_psu(i,j-1,k))/dy_m/dy_m
              dfz=kv_m2ps*(sa_psu(i,j,k+1) - 2*sa_psu(i,j,k) + sa_psu(i,j,k-1))/dz_m/dz_m
              
              gs = adx+ady+adz+dfx+dfy+dfz
              sc_psu(i,j,k) = sa_psu(i,j,k) + 2.d0 * dt_sec * gs
           end do
        end do
     end do

		 ! 対流調節（単純化。2回繰り返す）
     do n = 1, 2
        ! 予報したTSで密度を計算する
        do k = km, 1, -1
           do j = 1, jm
              do i = 1, im
                 rr_kgpm3(i,j,k) = r0_kgpm3 * (1.d0 - alphat_pk *   ( tc_c(i,j,k)   - t0_c ) &
                                &                   + alphas_ppsu * ( sc_psu(i,j,k) - s0_psu ) )
              end do
           end do
        end do

        ! 当該格子の密度と上の格子の密度を比べ、上の方が重ければ混ぜる
        do j = 1, jm
           do i = 1, im
              do k = km-1, 1, -1
                 if( rr_kgpm3(i,j,k+1) <= rr_kgpm3(i,j,k) ) cycle

                 tav = 0.5d0 * ( tc_c(i,j,k) + tc_c(i,j,k+1) )
                 tc_c(i,j,k  ) = tav
                 tc_c(i,j,k+1) = tav

                 sav = 0.5d0 * ( sc_psu(i,j,k) + sc_psu(i,j,k+1) )
                 sc_psu(i,j,k  ) = sav
                 sc_psu(i,j,k+1) = sav
              end do
           end do
        end do
     end do

     ! 東西境界条件
     uc_mps(0   ,:,:) = 0.d0
     uc_mps(im  ,:,:) = 0.d0
     uc_mps(im+1,:,:) = 0.d0  !- dummy

     vc_mps(0   ,:,:) = vc_mps(1,:,:)
     vc_mps(im+1,:,:) = vc_mps(im,:,:)

     ec_m(0   ,:)     = ec_m(1,:)
     ec_m(im+1,:)     = ec_m(im,:)

		 tc_c(0   ,:,:)   = tc_c(1,:,:)
		 tc_c(im+1,:,:)   = tc_c(im,:,:)

		 sc_psu(0   ,:,:) = sc_psu(1,:,:)
		 sc_psu(im+1,:,:) = sc_psu(im,:,:)

     ! 南北境界条件
     uc_mps(:,0   ,:) = uc_mps(:,1,:)
     uc_mps(:,jm+1,:) = uc_mps(:,jm,:)

     vc_mps(:,0   ,:) = 0.d0
     vc_mps(:,jm  ,:) = 0.d0
     vc_mps(:,jm+1,:) = 0.d0

     ec_m(:,0   )     = ec_m(:,1)
     ec_m(:,jm+1)     = ec_m(:,jm)

		 tc_c(:,0   ,:)   = tc_c(:,1,:)
		 tc_c(:,jm+1,:)   = tc_c(:,jm,:)

		 sc_psu(:,0   ,:) = sc_psu(:,1,:)
		 sc_psu(:,jm+1,:) = sc_psu(:,jm,:)

     ! 上下境界条件
     uc_mps(:,:,km+1) = uc_mps(:,:,km)+tx_npm2(:,:)*dz_m/rr_kgpm3(:,:,km)/av_m2ps
     uc_mps(:,:,0   ) = uc_mps(:,:,1)

     vc_mps(:,:,km+1) = vc_mps(:,:,km)+ty_npm2(:,:)*dz_m/rr_kgpm3(:,:,km)/av_m2ps
     vc_mps(:,:,0   ) = vc_mps(:,:,1)

		 tc_c(:,:,0   )   = tc_c(:,:,1)
		 tc_c(:,:,km+1)   = tc_c(:,:,km)

		 sc_psu(:,:,0   ) = sc_psu(:,:,1)
     sc_psu(:,:,km+1) = sc_psu(:,:,km) + dz_m / kv_m2ps * sf_psumps(:,:)

     ! n+1ステップが求まったので配列をシフトさせる
     ! ua_mps=ub_mps, ub_mps=uc_mps 
     ua_mps(:,:,:)    = ub_mps(:,:,:) +0.5d0*asf* ( ua_mps(:,:,:) - 2.d0*ub_mps(:,:,:) + uc_mps(:,:,:) )
     va_mps(:,:,:)    = vb_mps(:,:,:) +0.5d0*asf* ( va_mps(:,:,:) - 2.d0*vb_mps(:,:,:) + vc_mps(:,:,:) )
     ea_m(:,:)        = eb_m(:,:)     +0.5d0*asf* ( ea_m(:,:)     - 2.d0*eb_m(:,:)     + ec_m(:,:)     )
     sa_psu(:,:,:)    = sb_psu(:,:,:) +0.5d0*asf* ( sa_psu(:,:,:) - 2.d0*sb_psu(:,:,:) + sc_psu(:,:,:) )
     ta_c(:,:,:)      = tb_c(:,:,:)   +0.5d0*asf* ( ta_c(:,:,:)   - 2.d0*tb_c(:,:,:)   + tc_c(:,:,:)   )

     ! n+1 step -> ub_mps,vb_mps,eb_m,tb,sb
     ub_mps(:,:,:)    = uc_mps(:,:,:)
     vb_mps(:,:,:)    = vc_mps(:,:,:)
     eb_m(:,:)        = ec_m(:,:)
     sb_psu(:,:,:)    = sc_psu(:,:,:)
     tb_c(:,:,:)      = tc_c(:,:,:)

     time_sec = time_sec + dt_sec

     if( time_sec < time_to_output_sec ) cycle

     n = idnint( time_sec / output_interval_sec )

     ! 診断量の再計算（記録用）
     ! w
     do k = 1, km
        do j = 1, jm
           do i = 1, im
              ww_mps(i,j,k) = ww_mps(i,j,k-1) - ((ub_mps(i,j,k)-ub_mps(i-1,j,k))/dx_m + (vb_mps(i,j,k)-vb_mps(i,j-1,k))/dy_m)*dz_m
           end do
        end do
     end do
     ww_mps(:,:,km) = 0.d0

     ! p
     do j = 1, jm
        do i = 1, im
           pp_npm2(i,j,km) = gr_mps2*rr_kgpm3(i,j,km)*(dz_m*0.5d0 + eb_m(i,j))
           do k = km-1, 1, -1
              pp_npm2(i,j,k) = pp_npm2(i,j,k+1) + gr_mps2*(rr_kgpm3(i,j,k) + rr_kgpm3(i,j,k+1))*dz_m*0.5d0
           end do
        end do
     end do

     ! 記録用（解析・作図用）
     write(buff,'(a,i6.6)') trim(expid)//'.n',n
     open(10,file=trim(buff),form='unformatted',access='stream')
     write(10) time_sec,ub_mps,vb_mps,eb_m,ww_mps,pp_npm2,tb_c,sb_psu,rr_kgpm3
     close(10)


     ! 継続用
     open(10,file=trim(expid)//'.cnt',form='unformatted',access='stream')
     write(10) time_sec,ua_mps,ub_mps,va_mps,vb_mps,ea_m,eb_m,ta_c,tb_c,sa_psu,sb_psu
     close(10)


     ! 画面で監視
     write(6,*) 'time [s] = ',time_sec
     write(6,*) 'umin: ',minval(ub_mps(1:im,1:jm,1:km)),'(',minloc(ub_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'umax: ',maxval(ub_mps(1:im,1:jm,1:km)),'(',maxloc(ub_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'vmin: ',minval(vb_mps(1:im,1:jm,1:km)),'(',minloc(vb_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'vmax: ',maxval(vb_mps(1:im,1:jm,1:km)),'(',maxloc(vb_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'wmin: ',minval(ww_mps(1:im,1:jm,1:km)),'(',minloc(ww_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'wmax: ',maxval(ww_mps(1:im,1:jm,1:km)),'(',maxloc(ww_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'emin: ',minval(eb_m(1:im,1:jm)),'(',minloc(eb_m(1:im,1:jm)),')'
     write(6,*) 'emax: ',maxval(eb_m(1:im,1:jm)),'(',maxloc(eb_m(1:im,1:jm)),')'

     time_to_output_sec = time_to_output_sec + output_interval_sec

  end do

end program main
