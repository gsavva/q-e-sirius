!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
SUBROUTINE lr_addusddens (drhoscf, dbecsum)
  !---------------------------------------------------------------------------
  !
  ! Calculate the additional charge in reciprocal space due to US PP's
  ! See Eq.(36) in B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007)
  ! Then sum up the normal and ultrasoft charges.
  ! It assumes that the array dbecsum has already been computed.
  ! Inspired by PH/addusddens.f90
  !
  ! Created by Iurii Timrov (2013)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, tau, ntyp => nsp
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : gg, ngm, g, eigts1, eigts2, eigts3, mill
  USE uspp,                 ONLY : okvan
  USE wavefunctions, ONLY : psic
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE qpoint,               ONLY : xq, eigqts
  USE noncollin_module,     ONLY : nspin_mag
  USE mp_world,             ONLY : mpime, nproc
  USE mod_sirius
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout) :: drhoscf(dfftp%nnr, nspin_mag)
  ! input/output : change of the charge density
  COMPLEX(DP), INTENT(in)    :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag)
  ! input : the ultrasoft term
  !
  ! the local variables
  !
  INTEGER :: ig, na, nt, ih, jh, is, ijh, ir, nij, N_nt, na_
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  ! counter on r vectors
  ! counter on spin
  ! counter on combined beta functions
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:,:), aux_new(:,:), tmp0(:,:), tmp01(:,:), tmp1(:,:), tmp11(:,:), tmp2(:,:), tmp2_(:,:)
  COMPLEX(DP) :: z1
  real(8) :: diff
  ! the structure factor
  ! q_lm(G)
  ! auxiliary variable for drho(G)
  !
  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('lr_addusddens')
  !

  ALLOCATE (aux(ngm,nspin_mag))
! ALLOCATE (aux_new(ngm,nspin_mag))
  !
  aux(:,:) = (0.d0, 0.d0)
! aux_new(:,:) = (0.d0, 0.d0)
#if defined(__SIRIUS)
  DO nt = 1, ntyp
     IF (upf(nt)%tvanp) THEN
        CALL sirius_generate_rhoaug_q(gs_handler, nt, nat, ngm, nspin_mag, atom_type(nt)%qpw1, &
            & nh(nt) * (nh(nt) + 1) / 2, eigqts, mill, dbecsum, nhm * (nhm + 1) / 2, aux)
     ENDIF
  ENDDO
#else
  !
  DO nt = 1, ntyp ! loop over atom types
     IF (upf(nt)%tvanp) THEN  ! if US-PP
        ijh = 0
        !======================================================================= native QE implementation below
      !   DO ih = 1, nh (nt)
      !      DO jh = ih, nh (nt)
      !         !
      !         ijh = ijh + 1
      !         DO na = 1, nat
      !            IF (ityp (na) .eq.nt) THEN
      !               !
      !               ! Calculate the second term in Eq.(36) of the ultrasoft paper.
      !               !
      !               !$omp parallel default(shared) private(is, z1)
      !               DO is = 1, nspin_mag
      !                  !$omp do
      !                  DO ig = 1, ngm
      !                     !
      !                     ! Calculate the structure factor
      !                     !
      !                     z1 = eigts1(mill(1,ig),na) * &
      !                          eigts2(mill(2,ig),na) * &
      !                          eigts3(mill(3,ig),na) * &
      !                          eigqts(na)
      !                     !
      !                     aux(ig,is) = aux(ig,is) + 2.0d0 * atom_type(nt)%qpw(ig, ijh) * z1 * dbecsum(ijh,na,is)
      !                     !
      !                  ENDDO
      !                  !$omp end do nowait
      !               ENDDO
      !               !$omp end parallel
      !               !
      !            ENDIF
      !         ENDDO
      !      ENDDO
      !   ENDDO
         !======================================================================= new implementation below
         nij = nh(nt)*(nh(nt)+1)/2 ! max number of (ih,jh) pairs per atom type nt
         N_nt = 0 ! number of atoms of type nt
         DO na = 1, nat
            IF ( ityp(na) == nt ) N_nt = N_nt + 1
         ENDDO
         write(*,*) 'LR ===> nt, N_nt, nij, ngm, mpime, nproc', nt, N_nt, nij, ngm, mpime, nproc
         allocate (tmp2 (nij , N_nt)) ! Density Matrix for all atoms of a give type nt
         allocate (tmp2_(N_nt, nij ))

         allocate (tmp1 (N_nt, ngm)) ! qpw(ig,ijh) * tmp2(ijh,na)
         allocate (tmp11(N_nt, ngm))

         allocate (tmp0(ngm, N_nt))
         allocate (tmp01(ngm, N_nt))
         DO is = 1, nspin_mag ! loop over spins
            na_ = 0 ! count atoms of type nt
            !=============================================================== tmp2
            call start_clock ('lr_tmp2')
            DO na = 1, nat
               IF ( ityp(na) == nt ) then
                  na_ = na_ + 1
                  do ijh = 1, nij
                     tmp2 (ijh, na_) = dbecsum(ijh, na, is)
                     tmp2_(na_, ijh) = dbecsum(ijh, na, is)
                  enddo
               ENDIF
            ENDDO
            call stop_clock ('lr_tmp2')
            !=============================================================== GEMM
            call start_clock ('lr_gemm3')
            call ZGEMM('N', 'N', ngm, N_nt, nij, dcmplx(1.d0, 0.d0), &
                        atom_type(nt)%qpw, ngm, &
                        tmp2, nij, &
                        dcmplx(0.d0, 0.d0), tmp0, ngm)
            call stop_clock ('lr_gemm3')

            call start_clock ('lr_gemm2')
            call ZGEMM('T', 'T', N_nt, ngm, nij, dcmplx(1.d0, 0.d0), &
                          tmp2, nij, &
                          atom_type(nt)%qpw, ngm, &
                          dcmplx(0.d0, 0.d0), tmp11, N_nt)
            call stop_clock ('lr_gemm2')

            call start_clock ('lr_gemm4')
            call ZGEMM('T', 'N', ngm, N_nt, nij, dcmplx(1.d0, 0.d0), &
                        atom_type(nt)%qpw1, nij, &
                        tmp2, nij, &
                        dcmplx(0.d0, 0.d0), tmp01, ngm)
            call stop_clock ('lr_gemm4')

            call start_clock ('lr_gemm1')
            call ZGEMM('T', 'N', N_nt, ngm, nij, dcmplx(1.d0, 0.d0), &
                          tmp2, nij, &
                          atom_type(nt)%qpw1, nij, &
                          dcmplx(0.d0, 0.d0), tmp1, N_nt)
            call stop_clock ('lr_gemm1')

            call start_clock ('lr_gemm7')
            call ZGEMM('N', 'N', N_nt, ngm, nij, dcmplx(1.d0, 0.d0), &
                          tmp2_, N_nt, &
                          atom_type(nt)%qpw1, nij, &
                          dcmplx(0.d0, 0.d0), tmp1, N_nt)
            call stop_clock ('lr_gemm7')

            ! tmp2 is a complex array of dimension (nij, N_nt)
            ! qpw1 is a complex array of dimension (nij, ngm)
            ! ---- (tmp2)^H * qpw1 : (N_nt, nij) x (nij, ngm) = (N_nt, ngm), dimensions of tmp1
            ! ---- (qpw1)^H * tmp2 : (ngm, nij) x (nij, N_nt) = (ngm, N_nt), dimensions of tmp0
! tmp2_ (N_nt, nij)
! qpw ( ngm, nij) <<====== QE uses this one
! qpw1( nij, ngm) <<====== SIRIUS uses this one

!            3. M: rows of matrix op(A) and of matrix C.
!            4. N: cols of matrix op(B) and of matrix C.
!            5. K: cols of matrix op(A) and rows of matrix op(B). "Inner dimensions must agree"
!            6. ALPHA: Scalar multiplier for op(A)*op(B).
!            7. `A`: Array, dimension (LDA, ka), where ka is K when TRANSA = 'N' or 'n', and is M otherwise. The leading M by K part of the array A must contain the matrix A.
!            8. `LDA`: The leading dimension of the array A. When TRANSA = 'N' or 'n' then LDA must be at least max(1, M), otherwise LDA must be at least max(1, K).
!            9. `B`: Array, dimension (LDB, kb), where kb is N when TRANSB = 'N' or 'n', and is K otherwise. The leading K by N part of the array B must contain the matrix B.
!            10. `LDB`: The leading dimension of the array B. When TRANSB = 'N' or 'n' then LDB must be at least max(1, K), otherwise LDB must be at least max(1, N).
!            11. `BETA`: Scalar multiplier for C.
!            12. `C`: Array, dimension (LDC, N). The leading M by N part of the array C must contain the matrix C.
!            13. `LDC`: The leading dimension of the array C. LDC must be at least max(1, M).

            ! qpw is a complex array of dimension (ngm, nij)
            ! tmp is a complex array of dimension (ngm, N_nt)
            ! --- (qpw)^H * tmp : (nij,  ngm) x (ngm, N_nt) = (nij, N_nt), dimensions of res2
            !=============================================================== sum
            call start_clock ('lr_sum')
            na_ = 0
            DO na = 1, nat
               IF ( ityp(na) == nt ) then
                  na_ = na_ + 1
                  z1 = dcmplx(0.d0, 0.d0)
                  diff = 0.0d0
                  do ig = 1, ngm
                     z1 = tmp1(na_, ig) * ( eigts1(mill(1,ig),na) * &
                                            eigts2(mill(2,ig),na) * &
                                            eigts3(mill(3,ig),na) * &
                                            eigqts(na) )
                     aux_new(ig,is) = aux_new(ig,is) + 2.0d0 * z1
                     ! diff = diff + abs( aux_new(ig,is) - aux(ig,is) )
                  enddo
                  ! write(*,*)'nt, na, N_nt, na_, diff = ', nt, na, N_nt, na_, diff
               ENDIF
            ENDDO
            call stop_clock ('lr_sum')
            !===============================================================
            enddo ! loop over spins
            deallocate (tmp2, tmp2_)
            deallocate (tmp1, tmp11)
            deallocate (tmp0, tmp01)
         !======================================================================= new implementation ABOVE
         ENDIF ! if US-PP
  ENDDO ! nt
#endif
  !
  ! Convert aux to real space, and add to the charge density.
  !
!   aux = aux_new
  DO is = 1, nspin_mag
      !
      psic(:) = (0.d0, 0.d0)
      !
      DO ig = 1, ngm
         psic(dfftp%nl(ig)) = aux(ig,is)
      ENDDO
      !
      CALL invfft ('Rho', psic, dfftp)
      !
      DO ir = 1, dfftp%nnr
         drhoscf(ir,is) = drhoscf(ir,is) + psic(ir) 
      ENDDO
      !
  ENDDO
  !
  DEALLOCATE (aux)
  !
  CALL stop_clock ('lr_addusddens')
  !
  RETURN
  !
END SUBROUTINE lr_addusddens
