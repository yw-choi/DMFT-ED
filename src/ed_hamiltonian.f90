module ed_hamiltonian

    use dmft_params, only: norb, U, Jex, mu 
    use ed_params, only: nbath, kind_basis
    use ed_basis, only: basis_t, get_bitidx, ed_basis_get

    implicit none

    public :: generate_hamiltonian

    double precision, allocatable, public :: &
        ek(:,:),    &   ! ek(nsite,nspin)         impurity/bath onsite energies
        vk(:,:,:)       ! vk(norb,nbath,nspin)    impurity-bath hybridization

    private
contains

    ! generates the Hamiltonian for the sector in basis
    subroutine generate_hamiltonian(basis, H)
        type(basis_t), intent(in) :: basis
        double precision, intent(out) :: H(basis%ntot,basis%ntot)

        ! local variables
        integer(kind=kind_basis) :: basis_i, basis_j
        integer :: i,j,iorb,jorb,ibath,jbath,ispin, bidx

        H = 0.0d0

        ! calculates <i|H|j>. 
        ! 1. find nonvanishing coefficients of H|j> = SUM_k |k> <k|H|j>
        ! 2. set H(i,j) = <i|H|j>
        do j=1,basis%ntot
            basis_j = ed_basis_get(basis,j)
            do ispin=1,2
                do iorb=1,norb
                    bidx = get_bitidx(iorb,ispin)
                    if (BTEST(basis_j,bidx)) then

                    endif
                enddo
            enddo
        enddo

    end subroutine generate_hamiltonian
    
end module ed_hamiltonian
