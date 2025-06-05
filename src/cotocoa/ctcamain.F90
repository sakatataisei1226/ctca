module m_ctcamain
    use ctca
    use oh_type
    use paramt
    use allcom
#define OH_LIB_LEVEL 3
#include "ohhelp_f.h"
    implicit none
    public cotocoa_init, cotocoa_mainstep, cotocoa_finalize

    ! phiデータ
    real(8), allocatable :: phi_data(:)
    integer(kind=8) :: phi_data_size
    integer :: phi_areaid
    

    ! workerから受信するデータ
    real(8), allocatable :: recv_data(:)
    integer(kind=8) :: recv_data_size
    integer :: recv_areaid

contains

    subroutine cotocoa_init

        print *, "ctcainit check 1"
        phi_data_size = 10
        allocate(phi_data(phi_data_size))
        print *, "ctcainit check 2"
        call CTCAR_regarea_real8(phi_data, phi_data_size, phi_areaid)
        print *, "ctcainit check 3"
        recv_data_size = 10
        allocate(recv_data(recv_data_size))
        print *, "ctcainit check 4"
        call CTCAR_regarea_real8(recv_data, recv_data_size, recv_areaid)
        print *, "ctcainit check 5"
    end subroutine cotocoa_init

    subroutine cotocoa_mainstep
        implicit none
        integer :: i, j, from_rank
        from_rank = 0 ! workerのランク

        do i = 1, 10
            ! phiデータを更新してworkerに送る
            phi_data = 10.0d0 * i + [(j, j=1,phi_data_size)]
            call CTCAR_writearea_real8(phi_areaid, from_rank, 0, phi_data_size, phi_data)
            print *, "requester: sent phi_data =", phi_data

            ! workerからデータを受信
            call CTCAR_readarea_real8(recv_areaid, from_rank, 0, recv_data_size, recv_data)
            print *, "requester: received data from worker =", recv_data
        end do
    end subroutine cotocoa_mainstep

    subroutine cotocoa_finalize
    end subroutine cotocoa_finalize

end module m_ctcamain
