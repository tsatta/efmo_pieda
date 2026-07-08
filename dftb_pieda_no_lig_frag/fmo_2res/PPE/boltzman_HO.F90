  program main
  implicit none
  integer :: nconf_host,nconf_guest,nconf_complx
  character :: junk0
   double precision, allocatable :: PPE_host(:),   TS_host(:)
   double precision, allocatable :: PPE_guest(:),  TS_guest(:)
   double precision, allocatable :: PPE_complx(:), TS_complx(:)
   double precision :: ave_G_complx, ave_G_host, ave_G_guest
   double precision :: binding_E
   double precision :: E_HO_host, E_HO_guest, E_HO_complx
! read host
  read(5,*) junk0, nconf_host, E_HO_host
  allocate(PPE_host(nconf_host),TS_host(nconf_host))
  call rdfile(nconf_host,PPE_host,TS_host)
! read guest
  read(5,*) junk0, nconf_guest, E_HO_guest
  allocate(PPE_guest(nconf_guest),TS_guest(nconf_guest))
  call rdfile(nconf_guest,PPE_guest,TS_guest)
! read complex
  read(5,*) junk0, nconf_complx, E_HO_complx
  allocate(PPE_complx(nconf_complx),TS_complx(nconf_complx))
  call rdfile(nconf_complx,PPE_complx,TS_complx)

! post-processing for host,guest and complex
      write(*,*) '------------- Host ---------------'
      call post_processing(nconf_host,ave_G_host,TS_host,PPE_host,E_HO_host)
      write(*,*) '------------- Ligands ------------'
      call post_processing(nconf_guest,ave_G_guest,TS_guest,PPE_guest,E_HO_guest)
      write(*,*) '------------- Complex ------------'
      call post_processing(nconf_complx,ave_G_complx,TS_complx,&
                          PPE_complx,E_HO_complx)

! after post-processing, calculate binding E
      call calculate_binding_E(ave_G_host,ave_G_guest,ave_G_complx,&
                              binding_E)

  deallocate(PPE_host,TS_host)
  deallocate(PPE_guest,TS_guest)
  deallocate(PPE_complx,TS_complx)
  stop
  end program main

  subroutine rdfile(nconf,PPE,TS)
  implicit none
  double precision :: G(nconf),E(nconf), TS(nconf), PPE(nconf)
  integer :: i, nconf
  character :: junk0, junk1, junk2
  character :: ch1, ch2, ch3, ch4
  do i = 1, nconf
  read(5,*)  ch1, ch2, ch3, ch4, G(i), junk0, E(i), junk1, TS(i)
  enddo
  read (5,*)
  do i = 1, nconf
  read(5,*) PPE(i)
  enddo
  return
  end subroutine rdfile

!      SUBROUTINE input
!      implicit none
!
!      integer :: unit_inp
!
!      character(len=1000) :: inFile  
!      inFile = input
!      
!      open(unit_inp,file=inFile,status='old')
!
!      read(unit_inp,*) string
! 
!      call to_lowercase(string)
!      length = len_trim(string)
!    
!
!      if (string(1:length) .eq. 'host') then
!
!      else if (string(1:length) .eq. 'guest')
!
!     else if (string(1:length) .eq. 'complex')
!
!     else
!
!     endif
!    
!     RETURN
!      END SUBROUTINE input


      SUBROUTINE post_processing(nconf,ave_G,TS,PP_E,E_HO)
      implicit none

      double precision :: G(nconf), E(nconf),TS(nconf), PP_E(nconf), &
                         PP_G(nconf)

      double precision :: lnzx(nconf),zx_shift(nconf),ave_lnzx,&
                         sum_zx_shift,sum_lnzx,prob(nconf),&
                         sum_prob,ave_G
      double precision, parameter :: Rkcal = 1.98720661358D-03
      double precision :: T,RTkcal
      double precision :: E_HO

      integer :: i,nconf

        T = 300.00
        RTkcal = Rkcal*T
        ave_lnzx = 0.0D00
        lnzx(:)=0.0D00
        zx_shift(:)=0.0D00
        Do i = 1,nconf 
           PP_G(i) = PP_E(i) + TS(i) + E_HO
           lnzx(i) = -PP_G(i)/RTkcal
           ave_lnzx = ave_lnzx + lnzx(i)
        END do
        ave_lnzx = ave_lnzx/dble(nconf) 
       
        sum_zx_shift = 0.0D00 
        Do i = 1,nconf
           zx_shift(i) = exp(lnzx(i)-ave_lnzx)
       !write(*,*) 'ts comment zx_shift', lnzx(i)-ave_lnzx,  lnzx(i), ave_lnzx
           sum_zx_shift = sum_zx_shift + zx_shift(i)
        ENDdo
        sum_lnzx = log(sum_zx_shift)+ave_lnzx
 
        ave_G = - sum_lnzx*RTkcal

        sum_prob = 0.0D00 
        Do i = 1,nconf
           prob(i) = zx_shift(i)/sum_zx_shift*100.0D0
           sum_prob = sum_prob + prob(i)
          write(*,900) 'conf ', i, prob(i)
        Enddo
 900 format(A,I3, F8.2)
      RETURN
      END SUBROUTINE post_processing


      SUBROUTINE calculate_binding_E(E_host,E_guest,E_complx,binding_E)

      implicit none
  
      double precision :: binding_E,E_complx,E_host,E_guest
      write(*,*) 'E_host=',E_host
      write(*,*) 'E_guest=',E_guest
      write(*,*) 'E_complx=',E_complx

      binding_E = E_complx - E_host -E_guest
      
      write(*,*) 'binding energy=', binding_E
      RETURN
      END SUBROUTINE calculate_binding_E

