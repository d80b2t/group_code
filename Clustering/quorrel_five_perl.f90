!
! Room 403, Astronomy and Astrophysics Department
! Davey Lab, University Park, PA 16802
! Pennsylvania State University
!
! Okay, so this is the ``new'' and `improved' program that (primarily) should
! take SDSS Quasar data (mainly from DR5) and compute the redshift-space
! correlation function, xi(s), the ``2-D'' correlation function, xi(sigma,pi)
! and the projected correlation function, wp(sigma). wp(sigma) can, of course,
! be used and inverted/interpolated to workout, the real-space correlation, 
! xi(r), via a program such as correl_six.f90
!
! Nicholas P. Ross, Postdoctoral Research Scholar, 18th January 2008   

program correlfive
  
implicit none
  
integer:: number, ra1_rnd, ra2_rnd, dec1_rnd, dec2_rnd, bad_pi_para
integer:: counter, cumulate, counter_two, choice_two, repeat_counter, fld_s01
integer:: count_lt, count_mid, count_gt, q, istat, flag, x_int, y_int, area
integer:: pair_count, jack, z_max_int

integer, allocatable:: ra1(:), ra2(:), dec1(:), dec2(:), template(:), qa(:)
integer, allocatable:: lrgsam(:), run(:), rerun(:), camCol(:), field(:), id(:)
integer, allocatable:: ccd(:), plate(:), twodF_X(:), twodF_Y(:), qm(:), qx(:)
integer, allocatable:: bin(:), z_bin(:), z_rnd_bin(:)
integer, allocatable:: fib(:),qop(:),repeat(:),z_flag(:), f(:)

real*8, allocatable:: ra_J2K(:), dec_J2K(:), objc_rowc(:), objc_colc(:)
real*8, allocatable:: ideV(:), i_deV3(:), reff(:), g_minus_r(:), r_minus_i(:)
real*8, allocatable:: imod(:), ipet(:), rpet(:), fibCts3(:), fibCts2(:)
real*8, allocatable:: sn(:), spmag(:),vsp(:),rsp(:), reddening2(:),seeing2(:)
real*8, allocatable:: ra3(:), dec3(:), z_abs(:), r_mag(:), sky2(:), xcor(:)
real*8, allocatable:: x_coord(:), y_coord(:), z_coord(:), z_fin(:) 
real*8, allocatable:: x_coord_rnd(:), y_coord_rnd(:), z_coord_rnd(:)
real*8, allocatable:: r_c(:), r_fib(:), A_r(:), z_LRG(:)
real*8, allocatable:: bin_RR(:), bin_RR_sigma(:), xi_LS(:), xi_HAM(:)
real*8, allocatable:: xi_lin(:), delta_xi_lin(:), xi_LS_lin(:), xi_HAM_lin(:)
real*8, allocatable:: ra_theta(:), dec_theta(:), one_plus_ws(:), Mpc_sep(:)
real*8, allocatable:: theta_deg(:), theta_ams(:), w_theta_p(:)
real*8, allocatable:: w_theta_z(:), delta_w_theta_p(:), delta_w_theta_z(:) 
real*8, allocatable:: ra_DR5Q(:), dec_DR5Q(:), z_DR5Q(:), rc_DR5Q(:)
real*8, allocatable:: ra_rnd_npr(:), dec_rnd_npr(:), z_rnd_npr(:), rc_rnd(:)
real*8, allocatable:: ra_rnd(:), dec_rnd(:), z_rnd(:)
real*8:: ra3_rnd, dec3_rnd, ra_rad_rnd, dec_rad_rnd, ra_rad, dec_rad
real*8:: x_one,y_one,z_one,x_two,y_two,z_two,zero,x_three,y_three, z_three
real*8:: s,sa,sb,s_sq,s_squared,s_max,s_min,pi_max,pi_min, sigma_max,sigma_min
real*8:: my_sigma, delta_sigma, sigma_rad, sigma_deg, cos_theta, theta
real*8:: pi_para_one, sigma, pi_para_two, pi_para, pi_para_rnd, sigma_rnd
real*8:: z_min, z_max, hundred_z, zbin_width, z_offset, z_LRG_mean, para_one
real*8:: Landy, Szalay, k_max_RR, z_off, z_LRG_sample89_mean, z_sample_mean
real*8:: ra_temp, dec_temp, z_fin_temp
real*8:: ra_min,  ra_max,  ra_min_rad,  ra_max_rad
real*8:: dec_min, dec_max, dec_min_rad, dec_max_rad
real*8:: absdeV_I, num, sdss_id2, model_u, model_g, model_r, model_i, model_z
real*8:: ra_SDSS2, dec_SDSS2,  Landy2d,  Szalay2d
integer:: k, k_max_DD, k_max_DR, k_sigma, k_pi_para
integer:: k_lin 
integer:: N_LRG_sample8, N_LRG_sample9,N_LRG_sample4,N_LRG_sample3
integer:: i,j, i_max, j_max, i_min, j_min, ten_i, ip_counter, l, m, n, o
integer:: N_lns, N_data, N_rnd, N_star, d, x, dummy_i, choice, choice_2, N_LRGs
integer:: N_LRGs_north, N_LRGs_south, N_LRG_sample, N_in_03A, N_LRGs_03A
integer:: N_qop3_sample8, N_qop3_sample9,N_LRGs_sam89north,N_LRGs_sam89south 
integer:: N_LRG_SAMPLE8_QOP3, N_LRG_SAMPLE9_QOP3, N_LRGS_NORTH89,N_LRGS_SOUTH89
integer:: N_LRG_SAMPLE89, N_LRGS_SAMPLE89_03A,N_qop3_north, N_qop3_south
integer:: N_LRGs_sam_north, N_LRGs_sam_south, N_north, N_south, N_weirdos
integer:: N_LRGs_primarysam_north, N_LRGs_primarysam_south, N_star_N, N_star_S
integer:: N_LRGs_rminusi0pnt8, N_LRGs_a01a02s01, N_rnd_total

character:: dec_3char*3, choice_three*1
character, allocatable:: dec_sign(:)*1, dec_sign_rnd(:)*1, dec3char(:)*3
character, allocatable:: flags22(:)*8, flags23(:)*8, flags3(:)*8, flags2(:)*8
character, allocatable:: comments(:)*30,  SDSS_id(:)*18, date(:)*6
character, allocatable:: fld(:)*3, object_id(:)*20, object_sam_id(:)*20

real*8,  allocatable:: bin_DD(:), bin_DR(:), bin_DD_sigma(:), bin_DR_sigma(:)
real*8,  allocatable:: bin_DD_lin(:), bin_DR_lin(:), bin_RR_lin(:)
real*8,  allocatable:: Xi_sigma(:), delta_Xi_sigma(:), bin_sigma_pi(:,:)
real*8,  allocatable:: bin_DD_2d(:,:), bin_DR_2d(:,:), bin_RR_2d(:,:)
real*8,  allocatable:: ra_deg(:), dec_deg(:), xi(:), delta_xi(:) 
real*8,  allocatable:: xi_sigma_pi(:,:),delta_xi_sigma_pi(:,:), Xi_sigma_HAM(:)
real*8,  allocatable:: xi_sigma_pi_HAM(:,:), bin_sigma_pi_HAM(:,:)
real*8,  allocatable:: xi_sigma_pi_LS(:,:), Xi_sigma_LS(:),bin_sigma_pi_LS(:,:)
real*8:: lambda0, omega0, beta, pi, nRD_ratio, ra_off, dec_off, theta_cor
real*8:: rp,wp_norm,wperror,wp_div,wp, rp_lg, wp_lg, wp_div_siglg, wperror_lg
real*8:: rc, zplus1, dlum, dlzcomb, zplus1_rnd, z_mean8, z_mean89
real*8:: z_LRGs_rminusi0pnt8, sum_delta_xi_sigma_pi
real*8:: z_model_1, z_model_2, z_model_2a, z_model_2b, chisq1, chisq2
real*8:: matcher_lo, matcher_hi, matcher2_lo, matcher2_hi

pi = 3.141592654

write(*,*)
lambda0 = 0.7
omega0 = 0.3
beta = 0.0
call tabulate_dist(omega0, lambda0, beta) 
write(*,*)
!write(*,1010) 'omega,lambda,beta,EPS ', omega0, lambda0, beta
1010 format (a23, 2f6.3,2x, f4.2)

!open(108, file='redshift_comov_f90_temp.dat') 
!do i=1,10000
!   zplus1 = (i*0.0005) + 1.0
!   call search_dist(rc, zplus1, dlum, dlzcomb, 2)
!   write(108,*) zplus1, rc, dlum, dlzcomb
!end do
!close(108)

m = 80000   ! Size of data arrays, this could well need to be increased...
n = 1000000 ! Size of random arrays
o = 101     ! Size of xi(s) etc. arrays....
N_lns=0
N_data = 0
N_rnd = 0
N_star=0
N_star_N=0
N_star_S=0
N_north=0
N_south=0
N_LRGs = 0 
N_LRG_sample8 = 0
N_LRG_sample9 = 0
N_qop3_sample8= 0
N_qop3_sample9= 0
N_LRGs_north = 0 
N_LRGs_south = 0 
N_LRG_sample = 0
N_in_03A = 0
N_LRGs_03A = 0
N_LRG_sample8_qop3=0
N_LRG_sample9_qop3=0
N_LRGs_north89=0
N_LRGs_south89=0
N_LRGs_sample89_03A=0
N_LRGs_sam89north=0
N_LRGs_sam89south=0
N_LRG_sample89=0 
N_qop3_north=0
N_qop3_south=0
N_LRGs_sam_north=0
N_LRGs_sam_south=0
N_LRGs_primarysam_north=0
N_LRGs_primarysam_south=0
N_LRGs_a01a02s01=0

allocate(object_id(m),object_sam_id(m), SDSS_id(m), ra_J2K(m), dec_J2K(m))
allocate(ra_DR5Q(m), dec_DR5Q(m), z_DR5Q(m), rc_DR5Q(m))
allocate(x_coord(m), y_coord(m), z_coord(m))
allocate(ra_theta(m), dec_theta(m))
allocate(lrgsam(m), run(m), rerun(m), camCol(m), field(m), repeat(m))
allocate(id(m), objc_rowc(m), objc_colc(m), ideV(m), g_minus_r(m))
allocate(r_minus_i(m), imod(m), ipet(m), rpet(m), fibCts3(m))
allocate(fibCts2(m), i_deV3(m), reff(m), flags3(m), flags2(m))
allocate(flags23(m), flags22(m), reddening2(m), seeing2(m), sky2(m))
allocate(fib(m), ra1(m), ra2(m), ra3(m), dec1(m), z_LRG(m))
allocate(dec2(m), dec3(m), r_mag(m), template(m), xcor(m))
allocate(dec_sign(m), dec_sign_rnd(m), dec3char(m))
allocate(z_abs(m), qa(m), z_fin(m), qm(m), qop(m), qx(m), f(m))
allocate(z_flag(m), fld(m),  date(m), ccd(m), plate(m))
allocate(twodF_X(m), twodF_Y(m), sn(m), spmag(m), vsp(m))
allocate(rsp(m), ra_deg(m), dec_deg(m), comments(m), r_fib(m), A_r(m))

allocate(theta_deg(41), theta_ams(41), w_theta_p(41), w_theta_z(41))
allocate(delta_w_theta_p(41),delta_w_theta_z(41),Mpc_sep(41),one_plus_ws(41))
allocate(ra_rnd_npr(n), dec_rnd_npr(n), z_rnd_npr(n), rc_rnd(n))
allocate(ra_rnd(n), dec_rnd(n), z_rnd(n))
allocate(x_coord_rnd(n), y_coord_rnd(n), z_coord_rnd(n))

allocate(r_c(45000), z_rnd_bin(251))
allocate(bin_DD(101), bin_RR(101), z_bin(251), bin_DR(101), bin(101))
allocate(bin_DD_lin(101), bin_DR_lin(101), bin_RR_lin(101))
allocate(bin_DD_sigma(101),bin_DR_sigma(101),Xi_sigma(101),delta_Xi_sigma(99))
allocate(bin_DD_2d(101,101), bin_DR_2d(101,101), bin_RR_2d(101,101))
allocate(xi(101), delta_xi(101), bin_sigma_pi(101,101))
allocate(xi_LS(101), bin_RR_sigma(101), xi_HAM(101))
allocate(xi_lin(101), delta_xi_lin(101), xi_LS_lin(101), xi_HAM_lin(101))
allocate(xi_sigma_pi(101,101), delta_xi_sigma_pi(101,101))
allocate(xi_sigma_pi_HAM(101,101),bin_sigma_pi_HAM(101,101),Xi_sigma_HAM(101))
allocate(xi_sigma_pi_LS(101,101), Xi_sigma_LS(101), bin_sigma_pi_LS(101,101))

write(*,*)
write(*,*)
open(33, file='w_theta_2SLAQ_Sam8_v3.dat')
i=0
istat=0
do while(istat .eq. 0) 
   i=i+1
   read(33,*, IOSTAT=istat) theta_deg(i), theta_ams(i), w_theta_p(i), &
        & w_theta_z(i), delta_w_theta_p(i), delta_w_theta_z(i), &
        & one_plus_ws(i), Mpc_sep(i) !NB THe Mpc_sep is for LCDM!!
end do
close(33)
write(*,*) 'Read-in w_theta_2SLAQ.dat', i
write(*,*) '******* CURRENTLY NO FIBRE COLLISION IS APPLIED!!!! ***********'

write(*,*) choice_three
write(*,*) 'Do you want to use the Landy-Szalay estimator? y/n'
choice_three = 'y'
write(*,*) choice, '(choice)', choice_two, '(choice_two) ', choice_three, ' (choice_three)'


! I/P FILES
!open(1,file='DR5QSO_uni_data_z_cut_xyz.dat')
!open(2,file='randoms_npr_002_xyz.dat')
!
!
! could do an IF..CHOICE in here...
!
!open(1,file='../data/DR5QSO_uni_data.dat')
!open(2,file='../randoms/randoms_npr_UNIFORM.dat')
! Columns are; ra, dec, redshift, rc, 
! with rc comoing from IDL's lumdist.pro routine
! Omm=0.24, Lambda=0.76, H0=100.


!!  PRIMARY   PRIMARY   PRIMARY   PRIMARY   PRIMARY   PRIMARY 
!!  PRIMARY   PRIMARY   PRIMARY   PRIMARY   PRIMARY   PRIMARY 
!!  PRIMARY   PRIMARY   PRIMARY   PRIMARY   PRIMARY   PRIMARY 
!! For all DRs upto DR5Q:
!open(1,file='../data/DR5QSO_pri_data_20080508_field_filtered.dat')
!open(2,file='../randoms/randoms_npr_PRIMARY_20080519.dat')

!! For EDR/DR1:
!open(1,file='../data/DR5QSO_pri_data_DR1_20080811.dat')
!open(2,file='../randoms/randoms_npr_PRIMARY_DR1_20080811.dat')

!! For DR2345:
!open(1,file='../data/DR5QSO_pri_data_DR2345_20080811.dat')
!open(2,file='../randoms/randoms_npr_PRIMARY_DR2345_20080811.dat')



!!  UNIFORMS   UNIFORMS   UNIFORMS   UNIFORMS   UNIFORMS   UNIFORMS 
!!  UNIFORMS   UNIFORMS   UNIFORMS   UNIFORMS   UNIFORMS   UNIFORMS 
!!  UNIFORMS   UNIFORMS   UNIFORMS   UNIFORMS   UNIFORMS   UNIFORMS 
open(1,file='../data/DR5QSO_uni_data.dat')
!open(1,file='../data/DR5QSO_uni_data_wFibre_v2.dat')
open(2,file='../randoms/randoms_dat/randoms_npr_UNIFORM.dat')


!open(1,file='../data/DR5QSO_uni_data_EBV_le_0pnt0217.dat')
!open(2,file='../randoms/randoms_dat/randoms_npr_UNIFORM_EBV_le_0pnt0217.dat')
!
!open(1,file='../data/DR5QSO_uni_data_EBV_gt_0pnt0217.dat')
!open(2,file='../randoms/randoms_dat/randoms_npr_UNIFORM_EBV_gt_0pnt0217.dat')


!open(1,file='../data/DR5QSO_uni_Xray_data.dat')
!open(2,file='../randoms/randoms_dat/randoms_npr_UNI_Xray.dat')



! O/P FILES
open( 9,file='k_output_UNI22_0pnt92z1pnt13_rv2.dat')
open(10,file='k_output_lin_UNI22_0pnt92z1pnt13_rv2.dat')
open(18,file='k2d_output_UNI22_0pnt92z1pnt13_rv2.dat')
open(27,file='K_wp_output_UNI22_0pnt92z1pnt13_rv2.dat')
!
open(44,file='z_values_UNI22_0pnt92z1pnt13_rv2.dat')
open(55,file='ra_dec_z_DR5Q_used_UNI22_0pnt92z1pnt13_rv2.dat')
open(56,file='ra_dec_z_randoms_used_UNI22_0pnt92z1pnt13_rv2.dat')
!
!!open(2,file='sanity_check.dat')
!!open(2, file='2df_LRG_formatted6970.dat')
!!open(14, file='2df_LRG_formatted.dat')
!!open(89, file='sanity_check.dat')
!!open(91, file='capital_xi.dat')
!!open(77,file='wp_sig_2dfgrs.dat')

count_lt  = 0
count_mid = 0
count_gt = 0
repeat_counter = 0
z_fin = 0.000
z_max = 0.000
z_min = 1.000
z_bin = 0
z_rnd_bin=0
z_LRG_mean = 0
z_mean8 = 0 
z_LRG_sample89_mean=0
z_sample_mean=0
xi = 0.0
bin_DD = 0.0
bin_DR = 0.0
bin_RR = 0.0
bin_DD_lin = 0.0
bin_DR_lin = 0.0
bin_RR_lin = 0.0
bin_DD_2d = 0.0
bin_DR_2d = 0.0
bin_RR_2d = 0.0
bad_pi_para =0

zbin_width = 0.05
z_offset = (zbin_width / 2.0)
zbin_width = (1.0/zbin_width)

ra_min = 0.0000
ra_max = 360.0000
!0.0 !36.47 !151.333 !160.01 !185.386 !201.410! 210.096! 222.335! 324.035 !360
!120.0 !145.0 !175.0 !195.0 !215.0 !240.0 !300.0 !330.0 !360.0 !0.0 !30.0!60.0

dec_min = -20.000
dec_max = 70.000


z_min = 0.920
z_max = 1.130
!! 0.0, 0.300, 0.680, 0.920, 1.130, 1.320, 1.500, 1.660, 1.83, 2.02, 2.20, 2.90 
!!



ra_min_rad = (ra_min/180.0) * pi
ra_max_rad = (ra_max/180.0) * pi
write(*,*)
write(*,*)         'ra_min, ra_max, ra_min_rad, ra_max_rad'
write(*,'(4f11.3)') ra_min, ra_max, ra_min_rad, ra_max_rad

dec_min_rad = (dec_min/180.0) * pi
dec_max_rad = (dec_max/180.0) * pi
write(*,*)         'dec_min, dec_max, dec_min_rad, dec_max_rad'
write(*,'(4f11.3)') dec_min, dec_max, dec_min_rad, dec_max_rad

write(*,*)        'z_min,   z_max'
write(*,'(2f6.3)') z_min, z_max
write(*,*)


i=0
istat=0
N_LRG_sample = 0
do while (istat .eq. 0) 
   i=i+1
   read(1, *, IOSTAT=istat) ra_DR5Q(i), dec_DR5Q(i), z_DR5Q(i), rc_DR5Q(i)

!! Not for jackknifes:
   if(  ((ra_DR5Q(i) .ge. ra_min) .and. (ra_DR5Q(i) .le. ra_max)) .and. & 
        ((dec_DR5Q(i) .ge. dec_min) .and. (dec_DR5Q(i) .le. dec_max))) then
   
!! For jackknifes:
!   if(  ((ra_DR5Q(i) .le. ra_min) .or. (ra_DR5Q(i) .ge. ra_max)) .or. & 
!        ((dec_DR5Q(i) .le. dec_min) .or. (dec_DR5Q(i) .ge. dec_max))) then
   
!! FOR SOUTHERN JACK-KNIFES:
!!
!   if(  ((ra_DR5Q(i) .ge. ra_min) .and. (ra_DR5Q(i) .le. ra_max))) then

      
      if((z_DR5Q(i) .ge. z_min) .and. (z_DR5Q(i) .le. z_max)) then
         
         N_LRG_sample = N_LRG_sample + 1
         ra_J2K(N_LRG_sample)  =  ra_DR5Q(i)
         dec_J2K(N_LRG_sample) = dec_DR5Q(i)
         z_fin(N_LRG_sample)   =   z_DR5Q(i)
         
         d = int(z_fin(N_LRG_sample) * zbin_width) + 1
         if(d .ge. 1 .and. d .le. 170) z_bin(d) = z_bin(d) + 1
         
         ra_rad  = ( ra_J2K(N_LRG_sample)/180.) * pi
         dec_rad = (dec_J2K(N_LRG_sample)/180.) * pi
         x_coord(N_LRG_sample) = rc_DR5Q(i) * cos(dec_rad) * cos(ra_rad)  
         y_coord(N_LRG_sample) = rc_DR5Q(i) * cos(dec_rad) * sin(ra_rad)  
         z_coord(N_LRG_sample) = rc_DR5Q(i) * sin(dec_rad)
         
         write(55,*) ra_J2K(N_LRG_sample), dec_J2K(N_LRG_sample), z_fin(N_LRG_sample) 
         
         !         end if
      end if
   end if
end do
write(*,*) 'Read-in data  points', i, ' of which', &
     & N_LRG_sample, 'are in the sample, i.e.', (N_LRG_sample/38208.)*100.,  &
     & '% of the FULL UNIFORM sample'
close(55)



i=0
istat=0
N_rnd=0
N_rnd_total = N_LRG_Sample * 20.
do while (istat .eq. 0) 
   i=i+1
   
   read(2,*,IOSTAT=istat) ra_rnd_npr(i), dec_rnd_npr(i), z_rnd_npr(i), rc_rnd(i)

   if( ((ra_rnd_npr (i)  .ge. ra_min)  .and. (ra_rnd_npr(i) .le. ra_max)) .and. &
        ((dec_rnd_npr(i) .ge. dec_min) .and. (dec_rnd_npr(i) .le. dec_max))) then
!   if( ((ra_rnd_npr (i) .le. ra_min) .or. (ra_rnd_npr(i) .ge. ra_max)) .or. &
!        ((dec_rnd_npr(i) .le. dec_min) .or. (dec_rnd_npr(i) .ge. dec_max))) then
!   if(  ((ra_DR5Q(i) .ge. ra_min) .and. (ra_DR5Q(i) .le. ra_max))) then

      if((z_rnd_npr(i) .ge. z_min) .and. (z_rnd_npr(i) .le. z_max)) then
         
         N_rnd = N_rnd + 1
         ra_rnd(N_rnd)  =  ra_rnd_npr(i)
         dec_rnd(N_rnd) = dec_rnd_npr(i)
         z_rnd(N_rnd)   =   z_rnd_npr(i)
         
         d = int(z_rnd(N_rnd) * zbin_width) + 1
         if(d .ge. 1 .and. d .le. 170) z_rnd_bin(d) = z_rnd_bin(d) + 1
         
         ra_rad_rnd  = ( ra_rnd(N_rnd)/180.) * pi
         dec_rad_rnd = (dec_rnd(N_rnd)/180.) * pi
         x_coord_rnd(N_rnd) = rc_rnd(i) * cos(dec_rad_rnd) * cos(ra_rad_rnd)  
         y_coord_rnd(N_rnd) = rc_rnd(i) * cos(dec_rad_rnd) * sin(ra_rad_rnd)  
         z_coord_rnd(N_rnd) = rc_rnd(i) * sin(dec_rad_rnd)
         
         write(56,*)  ra_rnd(N_rnd), dec_rnd(N_rnd), z_rnd(N_rnd) 
         
         ! end if
      end if
   end if
   
!   if (N_rnd .ge. N_rnd_total) istat =1
end do
!write(*,*) 'Read-in random points', N_rnd, 'equals', (N_rnd/173120.)*100.,'%'
write(*,*) 'Read-in random points',  i,  'equals', (N_rnd/660220.)*100.,'%'
write(*,*) 'N_rnd', N_rnd, 'N_LRG sample*20=', N_rnd_total
close(56) 









bin = 0
s_max = 0.0
pi_max = 0.0
sigma_max = 0.0
s_min = 10.0
pi_min = 10.0
sigma_min = 10.0
k_max_DD = 0
k_max_DR = 0
pair_count = 0

write(*,*)  'matcher_lo  matcher_hi  matcher2_lo matcher2_hi'
matcher_lo = 0.0
matcher_hi = 0.5
matcher2_lo = 0.0
matcher2_hi = 0.5
!read(*,*)  matcher_lo !, matcher_hi, matcher2_lo, matcher2_hi
write(*,*) matcher_lo, matcher_hi, matcher2_lo, matcher2_hi
write(*,*)
write(*,*)

! 1ST LOOP 1ST LOOP 1ST LOOP ! 1ST LOOP 1ST LOOP 1ST LOOP 
! 1ST LOOP 1ST LOOP 1ST LOOP ! 1ST LOOP 1ST LOOP 1ST LOOP 
! Now we're gonna try and work out the separation between these
! galaxies, these are the DATA-DATA point, DD.
!do i=1, 100
do i=1, N_LRG_sample-1
   if(i .eq. (int(N_LRG_sample/100))) write(*,*) '1 % done on DDs'
   if(i .eq. (int(N_LRG_sample/50))) write(*,*) '2 % done on DDs'
   if(i .eq. (int(N_LRG_sample/20))) write(*,*) '5 % done on DDs'
   if(i .eq. (int(N_LRG_sample/10))) write(*,*) '10 % done on DDs'

   if(i .eq. (int(N_LRG_sample/4))) write(*,*) '1/4 done on DDs'
   if(i .eq. (int((3*N_LRG_sample)/4))) write(*,*) '3/4 done on DDs'
!   pi_para_one = r_c(i)
   x_one = x_coord(i)
   y_one = y_coord(i)
   z_one = z_coord(i)
   zero = 0.0
   pi_para_one = (sqrt(s_squared(x_one,y_one,z_one,zero, zero,zero)))

   do j=i+1,N_LRG_sample
      !write(*,*) i,j
!      pi_para_two = r_c(j)
      x_two = x_coord(j)
      y_two = y_coord(j)
      z_two = z_coord(j)

      pi_para_two = (sqrt(s_squared(x_two,y_two,z_two,zero,zero,zero)))
      s_sq =  (s_squared(x_one,y_one,z_one,x_two,y_two,z_two))
  
      theta = sqrt( ((ra_theta(i) - ra_theta(j))**2) + &
           &   ((dec_theta(i)-dec_theta(j))**2) )

      !theta_cor =  one_plus_ws(int((log10(theta)*5.0)+22-0.5))
      theta_cor =  1.00
      
      s = (sqrt(s_sq))  
      
      pi_para = abs(pi_para_one - pi_para_two)
      if(pi_para .eq. 0.0) pi_para = 0.000001

      my_sigma = (sqrt(s_sq - (pi_para**2)))
      sigma=my_sigma
      
      if( (sigma .gt. matcher_lo .and. sigma .lt. matcher_hi) .and. &
           & (pi_para .gt. matcher2_lo .and. pi_para .lt. matcher2_hi))then
         pair_count = pair_count+1
      end if

      
      if (s .ne. 0 .and. sigma .ne. 0 .and. pi_para .ne. 0) then
         !if (s .ne. 0) then
         if(s.gt.s_max)  s_max = s
         if(s.lt.s_min)  s_min = s
         if(sigma.gt.sigma_max)  sigma_max = sigma
         if(sigma.lt.sigma_min)  sigma_min = sigma
         if(pi_para.gt.pi_max)  pi_max = pi_para
         if(pi_para.lt.pi_min)  pi_min =pi_para
         
         !k         = int((log10(s)*5) + 12)   ! Which "bucket" it goes into
         !!if(s .lt. 10.0) then
         k         = int((log10(s)*10.) + 12)
         k_lin     = int(s/50.) + 1
         !else
         !  k = (int(s/10.) + 21) 
         !end if
         k_sigma   = int((log10(sigma)*5) + 12) ! THIS HAS ALWAYS WORKED!!!!
         k_pi_para = int((log10(pi_para)*5) + 12) ! THIS HAS ALWAYS WORKED!!!!
         !k_sigma   = int((log10(sigma)*10) + 12) 
         !k_pi_para = int((log10(pi_para)*10) + 12)
         
         !k_sigma   = int(sigma)   +1 !USE FOR xi(sigma,pi) "Petal plots" 
         !k_pi_para = int(pi_para) +1 !USE FOR xi(sigma,pi) "Petal plots" 
         
         if (k .gt. k_max_DD)  k_max_DD = k
         !if (k .le. 24 .and.k .gt. 0) then !when k = int((log10(s)*5) + 12
         if (k .le. 56 .and.k .gt. 0) then
            if(theta .gt. 0.1) then
               bin_DD(k) = bin_DD(k) + 2 !WITHOUT w_theta correction!!
            else
               !bin_DD(k) = bin_DD(k) + (2 * theta_cor)
               bin_DD(k) = bin_DD(k) + 2
            end if
            !if (k_sigma .le. 56 .and. k_pi_para .le. 56 .and. &
            !if (k_sigma .le. 40 .and. k_pi_para .le. 40 .and. & !For Petals
            if (k_sigma .le. 22 .and. k_pi_para .le. 22 .and. &
                 & k_sigma .gt. 0 .and. k_pi_para .gt. 0) then                 
               
               if(theta .gt. 0.1) then
                  bin_DD_2d(k_sigma,k_pi_para) = &
                       &     bin_DD_2d(k_sigma,k_pi_para) + 2
               else
                  bin_DD_2d(k_sigma,k_pi_para) = &
                       &     bin_DD_2d(k_sigma,k_pi_para) + (2 * theta_cor)
               end if
               
            end if
         end if !(k .le. 24/56 .and.k .gt. 0)
         
         if (k_lin .le. 56 .and. k_lin .gt. 0) bin_DD_lin(k_lin)=bin_DD_lin(k_lin)+2 
            !WITHOUT w_theta correction!! !k_lin 56 => cut-off at s=275

         counter_two = counter_two +1 
      end if !(s .ne. 0 .and. sigma .ne. 0 .and. pi_para .ne. 0) then
      
      if (s.lt. 3.0) then 
         write(72,*) x_one, y_one, z_one,  pi_para_one
         write(72,*) x_two, y_two, z_two,  pi_para_two
         write(72,*) s,10**((float(k)-12+0.5)/10.), &
              & 10**((float(k)-11+0.5)/10.), k, bin_DD(k),  theta_cor
         write(72,*)
      end if

      if(s .lt. 100) then
         write(73,730)  i, j, ra_J2K(i), dec_J2K(i), z_fin(i), &
              ra_J2K(j), dec_J2K(j), z_fin(j) 
      end if
      
   end do !j = i+1, N_LRG_sample 
   !write(*,*) i
end do !do i=1,NLRGs - 
730 format (i6,1x, i6,1x, &
          &  f15.8,1x, f15.8,1x, f12.8,1x, &
          &  f15.8,1x, f15.8,1x, f12.8)

close(72)
close(73)


write(*,*)
write(*,*) 'Finished working out DDs, k_max_DD', k_max_DD
write(*,*)
write(*,*) 'counter',counter,'counter_two', counter_two
write(*,*) 's_min_DD', s_min, 's_max_DD', s_max
!write(*,*) 'pi_min_DD', pi_min, 'pi_max_DD', pi_max
!write(*,*) 'sigma_min_DD', sigma_min, 'sigma_max_DD', sigma_max
!write(*,*) 'i', i, 'j',j, 'int_s', int_s, 'int_s2', int_s2
write(*,*) 's_max, pi_max, sigma_max'
write(*,*)  s_max, pi_max, sigma_max
write(*,*) 's_min, pi_min, sigma_min'
write(*,*)  s_min, pi_min, sigma_min
write(*,*) 
write(*,*) '************* pair_count **************'
write(*,*) '*************', pair_count, '**************'
write(*,*)
write(*,*) 'N_LRGs', N_LRG_sample, 'N_rnd', N_rnd !, 'int_s', int_s
write(*,*)
!write(*,*) '(float(k)-12+0.5)/5., 10**((float(k)-12+0.5)/5.),  bin_DD(k)'
write(*,*) 'k, (float(k)-12+0.5)/10., 10**((float(k)-12+0.5)/10.),  bin_DD(k)'
do k=1,56
   if(bin_DD(k) .gt. 0) then
      write(*,*) k, (float(k)-12+0.5)/10., 10**((float(k)-12+0.5)/10.),  bin_DD(k)
   end if
end do
write(*,*)




if( choice_two .ne. 8 .and. choice_two .ne. 9) then
write(*,*) 'Now working out DRs, no of randoms', N_rnd  
!! 2ND LOOP 2ND LOOP 2ND LOOP  2ND LOOP 2ND LOOP 2ND LOOP
!! 2ND LOOP 2ND LOOP 2ND LOOP  2ND LOOP 2ND LOOP 2ND LOOP
!! This needs to be modified to include x_rand_coord etc. !!!!
!! We are now going to compare the separations between the real data
!! and a collection of random points, so this is the DR loop.
do i=1,N_LRG_sample
!do i=1,100
   if(i .eq. (int(N_LRG_sample/4))) write(*,*) '1/4 done on DRs'
   if(i .eq. (int((3*N_LRG_sample)/4))) write(*,*) '3/4 done on DRs'
   x_one = x_coord(i)
   y_one = y_coord(i)
   z_one = z_coord(i)
   pi_para_one = (sqrt(s_squared(x_one,y_one,z_one,zero,zero,zero)))
   
   do j=1,N_rnd 
      x_two = x_coord_rnd(j)
      y_two = y_coord_rnd(j)
      z_two = z_coord_rnd(j)
      pi_para_two = (sqrt(s_squared(x_two,y_two,z_two,zero,zero,zero)))
      s_sq =  s_squared(x_one,y_one,z_one,x_two,y_two,z_two)
      s = (sqrt(s_sq))  
      pi_para = abs(pi_para_one - pi_para_two)
      sigma=sqrt(s_sq-(pi_para**2)) 
      
      if(s.gt.s_max) s_max = s
      if(s.lt.s_min) s_min = s 
      
      if (s .ne. 0 .and. sigma .ne. 0 .and. pi_para .ne. 0) then
         k     = int((log10(s)*10) + 12)
         !k     = int((log10(s)*5) + 12)
         k_lin     = int(s/50.) + 1
         !else
         !   k = (int(s/10.) + 21) 
         !end if
         
         k_sigma = int((log10(sigma)*5) + 12)
         k_pi_para = int((log10(pi_para)*5) + 12)
         !k_sigma = int((log10(sigma)*10) + 12)
         !k_pi_para = int((log10(pi_para)*10) + 12)
         !k_sigma   = int(sigma)   +1 !USE FOR xi(sigma,pi) "Petal plots" 
         !k_pi_para = int(pi_para) +1 !USE FOR xi(sigma,pi) "Petal plots" 
         
         !if(k.le. 22) write(*,*) 's_sq', s_sq, 's', s,'DR k values', k
         if (k .gt. k_max_DR)  k_max_DR = k
         !if (k .le. 24 .and. k .gt. 0)   bin_DR(k) = bin_DR(k) + 1 
         if (k .le. 56 .and. k .gt. 0)   bin_DR(k) = bin_DR(k) + 1
         if (k_lin .le. 56 .and. k_lin .gt. 0) bin_DR_lin(k_lin)=bin_DR_lin(k_lin)+1 
         !if (k_sigma .le. 56 .and. k_pi_para .le. 56 .and. &
         !if (k_sigma .le. 40 .and. k_pi_para .le. 40 .and. &
         if (k_sigma .le. 22 .and. k_pi_para .le. 22 .and. &
              & k_sigma .gt. 0 .and. k_pi_para .gt. 0) then
            bin_DR_2d(k_sigma,k_pi_para) =bin_DR_2d(k_sigma,k_pi_para) + 1
         end if
      end if !(s .ne. 0) 
   end do !j=1,N_rnd
end do !i=1,N_LRG_sample 
write(*,*) 'DRs done'
write(*,*) 'i', i, 'j',j, 's_max,',s_max,'s_min', s_min
write(*,*)
write(*,*) 'log_s,   s,   xi(s),    delta_xi,   bin_DD,    bin_DR' 
nRD_ratio = ( real(N_rnd) / real(N_LRG_sample))

do k=1,56
   !if(bin_DR(k) .gt. 0) then
   if (bin_DR(k) .ne. 0 .and. bin_DD(k) .ne. 0) then
      
      xi(k) = (nRD_ratio *(bin_DD(k) / bin_DR(k))) - 1.0
      delta_xi(k) = (1 + xi(k)) * (sqrt(2.0 / bin_DD(k)))
      write(*,*) (float(k)-12+0.5)/10., 10**((float(k)-12+0.5)/10.), & 
!      write(*,*) (float(k)-12+0.5)/5., 10**((float(k)-12+0.5)/5.), & 
           & xi(k),  delta_xi(k), &
           & bin_DD(k), bin_DR(k)
!      write(*,*) (float(k)-12+0.5)/10., 10**((float(k)-12+0.5)/10.),  bin_DR(k), ((N_rnd/N_LRG_sample)*(bin_DD(k)/ bin_DR(k))-1.)
   end if
end do
end if!( choice_two .ne. 8 .or. choice_two .ne. 9) then
write(*,*)



if(choice_three .eq. 'y') then
write(*,*) 'Now working out RRs, no of randoms', N_rnd  
!! 3RD LOOP 3RD LOOP 3RD LOOP  3RD LOOP 3RD LOOP 3RD LOOP  3RD LOOP 3RD LOOP
!! 3RD LOOP 3RD LOOP 3RD LOOP  3RD LOOP 3RD LOOP 3RD LOOP  3RD LOOP 3RD LOOP 
!! 3RD LOOP 3RD LOOP 3RD LOOP  3RD LOOP 3RD LOOP 3RD LOOP  3RD LOOP 3RD LOOP
do i=1,N_rnd-1
   do l=1,10
      if(i .eq. (int((l*N_rnd)/10))) write(*,*) l*10,'% done on RRs'
   end do
   x_one = x_coord_rnd(i)
   y_one = y_coord_rnd(i)
   z_one = z_coord_rnd(i)
   pi_para_one = (sqrt(s_squared(x_one,y_one,z_one,zero,zero,zero)))
   do j=i+1,N_rnd
      x_two = x_coord_rnd(j)
      y_two = y_coord_rnd(j)
      z_two = z_coord_rnd(j)
      pi_para_two = (sqrt(s_squared(x_two,y_two,z_two,zero,zero,zero)))
      s_sq =  s_squared(x_one,y_one,z_one,x_two,y_two,z_two)
      s = (sqrt(s_sq))  
      pi_para = abs(pi_para_one - pi_para_two)
      sigma=(sqrt(s_sq-(pi_para**2)))
      
      if(s.gt.s_max)  s_max = s
      if(s.lt.s_min)  s_min = s 
      
      if (s .ne. 0 .and. sigma .ne. 0 .and. pi_para .ne. 0) then
      !if (s .ne. 0) then
         !k        = int((log10(s)*5) + 12) ! THIS HAS ALWAYS WORKED!!!!
         !if(s .lt. 10.0) then
         k        = int((log10(s)*10) + 12)
         k_lin     = int(s/50.) + 1
         !else
         !   k = (int(s/10.) + 21) 
         !end if

         k_sigma = int((log10(sigma)*5) + 12)
         k_pi_para = int((log10(pi_para)*5) + 12)  
         !k_sigma = int((log10(sigma)*10) + 12)
         !k_pi_para = int((log10(pi_para)*10) + 12)
         !k_sigma   = int(sigma)   +1 !USE FOR xi(sigma,pi) "Petal plots" 
         !k_pi_para = int(pi_para) +1 !USE FOR xi(sigma,pi) "Petal plots" 
         
         !if(k.le. 22) write(*,*) 's_sq', s_sq, 's', s,'RR k values', k
         if (k .gt. k_max_RR)  k_max_RR = k
         !if (k .le. 24 .and. k .gt. 0)   bin_RR(k) = bin_RR(k) + 2 
         if (k     .le. 56 .and. k     .gt. 0) bin_RR(k)     = bin_RR(k)     + 2 
         if (k_lin .le. 56 .and. k_lin .gt. 0) bin_RR_lin(k_lin) = bin_RR_lin(k_lin)+2
         !if (k_sigma .le. 56 .and. k_pi_para .le. 56 .and. &
         !if (k_sigma .le. 40 .and. k_pi_para .le. 40 .and. &
         if (k_sigma .le. 22 .and. k_pi_para .le. 22 .and. &
              & k_sigma .gt. 0 .and. k_pi_para .gt. 0) then
            bin_RR_2d(k_sigma,k_pi_para) =bin_RR_2d(k_sigma,k_pi_para) + 2
         end if
      end if ! if (s .ne. 0 .and. sigma .ne. 0 .and. pi_para .ne. 0) then
      
   end do !j=1,N_rnd
   close(17)
end do !i=1,N_rnd
write(*,*) 'RRs done'
write(*,*) 'i', i, 'j',j
end if !choice_three .eq. 'y'




if( choice_two .ne. 8 .and. choice_two .ne. 9) then
nRD_ratio = ( real(N_rnd) / real(N_LRG_sample))
write(*,*) 'nRD_ratio', nRD_ratio, 'N_rnd', N_rnd, & 
     & 'N_data', N_data, 'N_LRGs', N_LRGs, 'N_LRG_sample', N_LRG_sample
!open(9,file='k_output.dat')
!write(9,*) '# float(k)-0.5)/5, xi, delta_xi, bin_DD, bin_DR'
!write(*,*) '(float(k)-12+0.5)/5., 10**((float(k)-12+0.5)/5.), & 
write(*,*) '(float(k)-12+0.5)/10., 10**((float(k)-12+0.5)/10.), & 
     & xi(k), delta_xi(k), &
     & bin_DD(k), bin_DR(k), bin_RR(k), bin_DR(k)/nRD_ratio'

do k=1,56
   !   if (bin_DR(k) .ne. 0 .and. bin_DD(k) .ne. 0) then
   if (bin_DR(k) .ne. 0) then
      
      xi(k) = (nRD_ratio *(bin_DD(k) / bin_DR(k))) - 1.0
      delta_xi(k) = -1.
      
      if (bin_DD(k) .ne. 0) then
         delta_xi(k) = (1 + xi(k)) * (sqrt(2.0 / bin_DD(k)))
      endif
      
      !write(9,*) (float(k)-12+0.5)/5., 10**((float(k)-12+0.5)/5.), & 
      !    & xi(k),delta_xi(k), bin_DD(k), bin_DR(k), bin_DR(k)/nRD_ratio
      
      if (choice_three .eq. 'y' .and. bin_RR(k) .ne.  0) then 
         Landy = (nRD_ratio**2) * ( bin_DD(k)/bin_RR(k) )  
         Szalay = 2.0 * nRD_ratio * ( bin_DR(k) / bin_RR(k)) 
         xi_LS(k)  = 1 + Landy - Szalay
         xi_HAM(k) = ((bin_DD(k) * bin_RR(k)) / (bin_DR(k)**2)) - 1.0
         
      end if
      
      !if(k .lt. 22) then
!      write(9,1500) (float(k)-12+0.5)/5., 10**((float(k)-12+0.5)/5.), & 
      write(9,1500) (float(k)-12+0.5)/10., 10**((float(k)-12+0.5)/10.), & 
           & xi(k),  delta_xi(k), &
           & bin_DD(k), bin_DR(k), bin_RR(k), &
           & xi_LS(k), xi_HAM(k), bin_DR(k)/nRD_ratio
      !else 
      !   write(9,*) log10((k-21+0.5)*10.), (k-21+0.5)*10., & 
      !        & xi(k),  delta_xi(k), &
      !        & bin_DD(k), bin_DR(k), xi_LS(k), xi_HAM(k), &
      !        & bin_RR(k), bin_DR(k)/nRD_ratio
      !end if
      
      !write(*,*) (float(k)-12+0.5)/5.,xi(k), delta_xi(k), &
      !     & bin_DD(k), bin_DR(k) !, bin_DR(k)/nRD_ratio
      ! and now just write out, k, then the r that gave THAT k, a DD and a DR

   end if ! if (bin_DR(k) .ne. 0 .and. bin_DD(k) .ne. 0) th
end do !k=1,56
close(9)
write(*,*) 'xi(s) HAS BEEN CALCULATED!'
1500 format (f5.2,1x,   f15.6,1x, &
          &  f15.8,1x, f18.13,1x,   &
          &  f18.0,1x,    f18.0,1x, f18.0,1x, &
          &  f18.10,1x, f18.10,1x, f21.5)




do k_lin=1,56
   if (bin_DR_lin(k_lin) .ne. 0 .and. bin_DD_lin(k_lin) .ne. 0) then
      
      xi_lin(k_lin) = (nRD_ratio *(bin_DD_lin(k_lin) / bin_DR_lin(k_lin))) - 1.0
      delta_xi_lin(k_lin) = (1 + xi_lin(k_lin)) * (sqrt(2.0 / bin_DD_lin(k_lin)))

      if (choice_three .eq. 'y' .and. bin_RR_lin(k_lin) .ne.  0) then 
         Landy = (nRD_ratio**2) * ( bin_DD_lin(k_lin)/bin_RR_lin(k_lin) )  
         Szalay = 2.0 * nRD_ratio * ( bin_DR_lin(k_lin) / bin_RR_lin(k_lin)) 
         xi_LS_lin(k_lin)  = 1 + Landy - Szalay
         xi_HAM_lin(k_lin) = ((bin_DD_lin(k_lin) * bin_RR_lin(k_lin)) / (bin_DR_lin(k_lin)**2)) - 1.0
      end if
      
      write(10,1510) k_lin, (k_lin-.5)*5., & 
           & xi_lin(k_lin),  delta_xi_lin(k_lin), &
           & bin_DD_lin(k_lin), bin_DR_lin(k_lin), bin_RR_lin(k_lin), &
           & xi_LS_lin(k_lin), xi_HAM_lin(k_lin), bin_DR_lin(k_lin)/nRD_ratio

   end if ! if (bin_DR_lin(k_lin) .ne. 0 .and. bin_DD_lin(k_lin) .ne. 0) th
end do !k_lin=1,56
write(*,*) 'xi(s) (linear) HAS BEEN CALCULATED!'
1510 format (i4,1x,   f8.4,1x, &
          &  f15.8,1x, f18.13,1x,   &
          &  f18.0,1x,    f18.0,1x, f18.0,1x, &
          &  f18.10,1x, f18.10,1x, f21.5)
close(10)




!open(18,file='k2d_output.dat')
!open(27,file='K_output.dat')
xi_sigma_pi=0.0
xi_sigma_pi_HAM=0.0
xi_sigma_pi_LS=0.0
bin_sigma_pi=0.0
bin_sigma_pi_HAM=0.0
bin_sigma_pi_LS=0.0
Xi_sigma=0.0
Xi_sigma_HAM=0.0
Xi_sigma_LS=0.0

do k_sigma=1,22
!do k_sigma=1,56 
!do k_sigma=1,41
   sum_delta_xi_sigma_pi = 0.00

   do k_pi_para=1,22
   !do k_pi_para =1,56 
   !do k_pi_para =1,41
 
      !! Do we want to calculate xi(si,pi) with the LS estimator
      !! when the DD counts are zero? Sometimes no...
      if( (bin_DR_2d(k_sigma,k_pi_para) .ne. 0) .and. &
           & (bin_DD_2d(k_sigma,k_pi_para) .ne. 0)) then
      !! Sometimes yes...
!      if(bin_DR_2d(k_sigma,k_pi_para) .ne. 0) then
         
         xi_sigma_pi(k_sigma,k_pi_para) = (-1.0) + (nRD_ratio * &
              & (bin_DD_2d(k_sigma,k_pi_para) / bin_DR_2d(k_sigma,k_pi_para)))

         xi_sigma_pi_HAM(k_sigma,k_pi_para) = ((bin_DD_2d(k_sigma,k_pi_para)* &
              & bin_RR_2d(k_sigma,k_pi_para))/(bin_DR_2d(k_sigma,k_pi_para)**2)) - 1.0         
        
         Landy2d = (nRD_ratio**2) * &
              & (bin_DD_2d(k_sigma,k_pi_para)/bin_RR_2d(k_sigma,k_pi_para))
         Szalay2d =  (2.0 * nRD_ratio) * &
              & (bin_DR_2d(k_sigma,k_pi_para)/bin_RR_2d(k_sigma,k_pi_para))
         xi_sigma_pi_LS(k_sigma,k_pi_para) =  1 + Landy2d - Szalay2d

         delta_xi_sigma_pi(k_sigma,k_pi_para) = -1.0
         if(bin_DD_2d(k_sigma,k_pi_para) .ne. 0   .and. &
              k_pi_para .le. 20) then
            delta_xi_sigma_pi(k_sigma,k_pi_para) = (1.0 + &
                 & xi_sigma_pi_LS(k_sigma,k_pi_para)) &
                 & * (sqrt(2.0 / bin_DD_2d(k_sigma,k_pi_para))) 
            
            sum_delta_xi_sigma_pi = sum_delta_xi_sigma_pi + &
                 delta_xi_sigma_pi(k_sigma,k_pi_para)
         endif

         
         if(k_pi_para .le. 20) then 
            !le 19 => Putting a pi_para cut in at 31 h^-1 Mpc,  5 bins per decade
            !le 20 => Putting a pi_para cut in at 63 h^-1 Mpc,  5 bins per decade
            !if(k_pi_para .le. 39) then !for 5 bins/decade=>21; 10/dec => pi=39
            !if(k_pi_para .le. 40) then !Putting a pi_para cut in at 40Mpc 
            
            !if(k_pi_para .eq. 1) then 
            ! bin_sigma_pi(k_sigma,k_pi_para)=xi_sigma_pi(k_sigma,k_pi_para) &
            !      & * (10**((float(k_pi_para)-11)/10.))
            !& * (10**((float(k_pi_para)-12)/10.))
            !else
            bin_sigma_pi(k_sigma,k_pi_para)=xi_sigma_pi(k_sigma,k_pi_para) &
                 & *(10**((float(k_pi_para)-11)/5.) &
                 &  - 10**((float(k_pi_para)-12)/5.)) 
            bin_sigma_pi_HAM(k_sigma,k_pi_para) = &
                 & xi_sigma_pi_HAM(k_sigma,k_pi_para) &
                 & *(10**((float(k_pi_para)-11)/5.) &
                 &  - 10**((float(k_pi_para)-12)/5.)) 
            bin_sigma_pi_LS(k_sigma,k_pi_para) = &
                 & xi_sigma_pi_LS(k_sigma,k_pi_para) &
                 & *(10**((float(k_pi_para)-11)/5.) &
                 &  - 10**((float(k_pi_para)-12)/5.)) 
            
            !bin_sigma_pi(k_sigma,k_pi_para)=xi_sigma_pi(k_sigma,k_pi_para) &
            !     & *(10**((float(k_pi_para)-11)/10.) &
            !     &  - 10**((float(k_pi_para)-12)/10.)) 
            !end if
         end if
         
         Xi_sigma(k_sigma)=Xi_sigma(k_sigma) + &
              & (2*bin_sigma_pi(k_sigma,k_pi_para))
         Xi_sigma_HAM(k_sigma)=Xi_sigma_HAM(k_sigma)+ &
              & (2*bin_sigma_pi_HAM(k_sigma,k_pi_para))
         Xi_sigma_LS(k_sigma)= Xi_sigma_LS(k_sigma) + &
              & (2*bin_sigma_pi_LS(k_sigma,k_pi_para)) 

         bin_DD_sigma(k_sigma)=bin_DD_sigma(k_sigma) + &
              & (bin_DD_2d(k_sigma,k_pi_para))
         bin_DR_sigma(k_sigma)=bin_DR_sigma(k_sigma) + &
              & (bin_DR_2d(k_sigma,k_pi_para))
         
         write(98,*) (10**((float(k_sigma)-12+0.5)/5.)), &
              & k_sigma, k_pi_para,  xi_sigma_pi(k_sigma,k_pi_para), &
              & (10**((float(k_pi_para)-11)/5.)), &
              & (10**((float(k_pi_para)-12)/5.)), &
              & bin_sigma_pi(k_sigma,k_pi_para),  Xi_sigma(k_sigma), &
              & Xi_sigma_HAM(k_sigma),  Xi_sigma_LS(k_sigma)
         
         ! write(*,*) 10**((float(k_sigma)-12+0.5)/5.), &
         !     & 10**((float(k_sigma)-12+0.5)/5.), &
         !    & xi_sigma_pi(k_sigma,k_pi_para), &
         !    & bin_DD_2d(k_sigma,k_pi_para),  bin_DR_2d(k_sigma,k_pi_para)
         
      end if !(bin_DR_2d(k_sigma,k_pi_para).ne.0 .and. bin_DD_2d(k_sigma)...
      
      if (k_sigma .ge. 6 .and. k_pi_para .ge. 6) then    
         write(18,*) 10**((float(k_sigma)-12+0.5)/5.), & !-12+0.5/5.)
              & 10**((float(k_pi_para)-12+0.5)/5.), k_sigma, k_pi_para, &
              !write(18,*) 10**((float(k_sigma)-12+0.5)/10.), & !-12+0.5/5.)
              !    & 10**((float(k_pi_para)-12+0.5)/10.), k_sigma, k_pi_para, &
              !if (k_sigma .ge. 0 .and. k_pi_para .ge. 0) then    
              !  write(18,*) real(k_sigma),real(k_pi_para),k_sigma,k_pi_para, &
              & xi_sigma_pi(k_sigma,k_pi_para), &
              & delta_xi_sigma_pi(k_sigma,k_pi_para), & 
              & bin_DD_2d(k_sigma,k_pi_para),bin_DR_2d(k_sigma,k_pi_para), &   
              & bin_RR_2d(k_sigma,k_pi_para), &
              & xi_sigma_pi_HAM(k_sigma,k_pi_para), &
              & xi_sigma_pi_LS(k_sigma,k_pi_para), &
              & sum_delta_xi_sigma_pi
      end if
   end do !k_pi_para=1,22      
   

   
!   if(bin_DD_sigma(k_sigma) .ne. 0 .and. bin_DR_sigma(k_sigma) .ne. 0) then
!   if(bin_DD_sigma(k_sigma) .ne. 0 ) then
!      delta_Xi_sigma(k_sigma) = (1 + Xi_sigma(k_sigma)) &
!           & * (sqrt(2.0 / bin_DD_sigma(k_sigma)))
!   else
!     delta_Xi_sigma(k_sigma) = 0.
!   end if

     delta_Xi_sigma(k_sigma) =  sum_delta_xi_sigma_pi 

   !   27,file= K_output_2SLAQ.dat
   if(Xi_sigma(k_sigma) .ne. 0.0) then
      write(27,*) 10**((float(k_sigma)-12+0.5)/5.), &
           & (float(k_sigma)-12+0.5)/5., &
           !write(27,*) 10**((float(k_sigma)-12+0.5)/10.), &
           !     & (float(k_sigma)-12+0.5)/10., &
           & Xi_sigma(k_sigma), delta_Xi_sigma(k_sigma), &
           & bin_DD_sigma(k_sigma),bin_DR_sigma(k_sigma), &
           & Xi_sigma_HAM(k_sigma), Xi_sigma_LS(k_sigma)
   end if
   
end do !k_sigma =1,22

end if !if( choice_two .ne. 8 .or. choice_two .ne. 9) then
close(18)
close(27)
close(98)

close(1)
close(2)


write(*,*)
write(*,*)
write(*,*) 'almost there..'
write(*,*) 'z_max', z_max, 'zbin_wdith', zbin_width
N_data = N_LRG_sample
N_star = 1
!if( choice_two .ne. 8 .or. choice_two .ne. 9) z_max = 1.1
z_max_int = int(z_max * zbin_width) + 2
write(*,*) 'z_max_int', z_max_int
do i = z_max_int,1,-1
!do i = 175,1,-1
   hundred_z = real( i / zbin_width)  
   cumulate = cumulate + int(z_bin(i))
   write(*,1600) hundred_z, int(z_bin(i)), int(z_rnd_bin(i))
   write(44,1600) hundred_z, int(z_bin(i)),  int(z_rnd_bin(i))
end do
!write(44,1600) 0.000, N_star,1.000
!write(44,1600) -0.020, 0 ,1.000
write(*,*) ' Closing time.'
1600 format (f6.3,1x, i6,2x, i6,2x)


close(44)
!close(89)


!end do !jack =1,32

write(*,*) 'DONE AND DUSTED'

end program correlfive






function s_squared (a_one,b_one,c_one,a_two,b_two,c_two)
  implicit none
  
  real*8:: s_squared
  real*8:: a_one,b_one,c_one,a_two,b_two,c_two

  !write(*,*) 'in the function', a_one,b_one,c_one,a_two,b_two,c_two
  s_squared = ((a_two-a_one)**2 + (b_two-b_one)**2 + (c_two-c_one)**2) 
  
  return
end function  s_squared








