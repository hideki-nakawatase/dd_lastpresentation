program main

  !use linear_equation_of_state
  !use nonlinear_equation_of_state

  implicit none

  character(15):: &
       expid='ex4-case1' ! 実験名

  integer,parameter:: &
       im = ???, & ! x方向格子数
       !省略
       
  real(8),parameter:: &
       ! 省略
       kh_m2ps  = 5.d3,   & ! 水平渦拡散係数
       kv_m2ps  = 1.d-4,  & ! 鉛直渦拡散係数
       gm_psec  = 1.d0/50.d0/86400.d0, & ! 温位緩和係数
       ! 省略
       
  real(8),dimension(:,:,:):: &
       ! 省略
       ta_c(0:im+1,0:jm+1,0:km+1),   & ! 温位
       tb_c(0:im+1,0:jm+1,0:km+1),   &
       tc_c(0:im+1,0:jm+1,0:km+1),   &
       sa_psu(0:im+1,0:jm+1,0:km+1), & ! 塩分 （34.0からの偏差）
       sb_psu(0:im+1,0:jm+1,0:km+1), &
       sc_psu(0:im+1,0:jm+1,0:km+1)

  real(8),dimension(:,:,:):: &
       ! 省略
       rr_kgpm3(0:im+1,0:jm+1,0:km+1)   ! 密度

  real(8),dimension(:,:):: &
       ! 省略
       at_c(0:im+1,0:jm+1),      & ! 緩和温位（気温）
       sf_psumps(0:im+1,0:jm+1)    ! 海面塩分フラックス [psu m/s]

  real(8):: &
       ! 省略
       pre, adx, ady, adz, cor, dfx, dfy, dfz, frc, &
       tav, sav ! 作業用

  ! 省略

  ! 強制力
  dpi=4.d0*datan(1.d0)
  do j=1,jm
     do i=1,im
        tx_npm2(i,j) = t0_npm2 * cos( dpi*( (dble(j)-0.5d0) / dble(jm) ) ) !亜熱帯循環
        ty_npm2(i,j) = 0.d0
        at_c(i,j)      = -20.d0 *( dble(j)-0.5d0 ) / dble(jm) + 30.d0 ! 東西一様水温
        sf_psumps(i,j) = 0.d0
     end do
  end do

  ! 初期化
  ! 省略

  ! 初期設定
  if( time_to_start_sec == 0.d0 ) then ! 初期化
     ! 省略
     ta_c(:,:,:)   = 30.d0
     tb_c(:,:,:)   = 30.d0
     sa_psu(:,:,:) = 0.d0 !- 34からの偏差
     sb_psu(:,:,:) = 0.d0
     ! 省略
     
  else ! 継続
     ! 省略
  endif

  ! 省略

  do
     if( time > time_to_end ) exit

     ! 診断量の計算
     
     ! r
     do k = km, 1, -1
        !pre_npm2 = pre_npm2 + gr_mps2 * r0_kgpm3 * dz_m !- 圧力はrrで計算すべきだがr0で近似
        do j = 1, jm
           do i = 1, im
              rr_kgpm3(i,j,k) = ???
           end do
        end do
     end do

     ! p
     do j = 1, jm
        do i = 1, im
           do k = ?
              pp_npm2(i,j,k) = ???
           end do
        end do
     end do

     ! 予報量の計算
     ! u
     ! v
     ! e
     ! 省略

     ! t
     do k = 1, km
        do j = 1, jm
           do i = 1, im
              adx=???
              ady=???
              adz=???
              dfx=???
              dfy=???
              dfz=???

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
              adx=???
              ady=???
              adz=???
              dfx=???
              dfy=???
              dfz=???
              
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
     ! 省略
     
     ! 南北境界条件
     ! 省略

     ! 上下境界条件
     ! 省略

     sc_psu(:,:,km+1) = sc_psu(:,:,km) + dz_m / kv_m2ps * sf_psumps(:,:)


     ! n+1ステップが求まったので配列をシフトさせる
     ! 省略
     
     if( time_sec < time_to_output_sec ) cycle

     n = idnint( time_sec / output_interval_sec )

     ! 診断量の再計算（記録用）
     ! r
     ! 省略
     ! w
     ! 省略
     ! p
     ! 省略

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
     ! 省略

  end do

end program main
