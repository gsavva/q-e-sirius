!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine stres_har (sigmahar)
  !----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : e2, fpi
  USE cell_base, ONLY: omega, tpiba2
  USE ener,      ONLY: ehart
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft
  USE gvect,     ONLY: ngm, gstart, nl, g, gg, mill
  USE lsda_mod,  ONLY: nspin
  USE scf,       ONLY: rho
  USE control_flags,        ONLY: gamma_only
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum
  USE input_parameters, ONLY : use_sirius
  use sirius

  implicit none
  !
  real(DP) :: sigmahar (3, 3), shart, g2
  real(DP), parameter :: eps = 1.d-8
  integer :: is, ig, l, m, nspin0
  !  
  call sirius_start_timer(c_str("qe|stress_har"))

  sigmahar(:,:) = 0.d0
  psic (:) = (0.d0, 0.d0)

  if (use_sirius) then
    call sirius_get_rho_pw(ngm, mill(1, 1), rho%of_g(1, 1))
    psic(nl(:)) = rho%of_g(:, 1)
  else
    nspin0=nspin
    if (nspin==4) nspin0=1
    do is = 1, nspin0
       call daxpy (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, psic, 2)
    enddo

    CALL fwfft ('Dense', psic, dfftp)
  endif
  ! psic contains now the charge density in G space
  ! the  G=0 component is not computed
  do ig = gstart, ngm
     g2 = gg (ig) * tpiba2
     shart = psic (nl (ig) ) * CONJG(psic (nl (ig) ) ) / g2
     do l = 1, 3
        do m = 1, l
           sigmahar (l, m) = sigmahar (l, m) + shart * tpiba2 * 2 * &
                g (l, ig) * g (m, ig) / g2
        enddo
     enddo
  enddo
  !
  call mp_sum(  sigmahar, intra_bgrp_comm )
  !
  if (gamma_only) then
     sigmahar(:,:) =       fpi * e2 * sigmahar(:,:)
  else
     sigmahar(:,:) = 0.5d0 * fpi * e2 * sigmahar(:,:)
  end if
  do l = 1, 3
     sigmahar (l, l) = sigmahar (l, l) - ehart / omega
  enddo
  do l = 1, 3
     do m = 1, l - 1
        sigmahar (m, l) = sigmahar (l, m)
     enddo
  enddo

  sigmahar(:,:) = -sigmahar(:,:)

  call sirius_stop_timer(c_str("qe|stress_har"))
  return
end subroutine stres_har

