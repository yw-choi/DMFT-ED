module ed_hamiltonian

    use dmft_params, only: norb, U, Jex, mu, Up, Jp, nspin
    use utils, only: die
    use ed_params, only: nbath, kind_basis, nsite, ek_in, vk_in
    use ed_basis, only: basis_t, ed_basis_get,  ed_basis_idx
    use ed_operator, only: sgn
    implicit none

    public :: generate_hamiltonian
    public :: ed_hamiltonian_init

    double precision, allocatable, public :: &
        ek(:,:),    &   ! ek(nsite,2)         impurity/bath onsite energies
        vk(:,:,:)       ! vk(norb,nbath,2)    impurity-bath hybridization

    private
contains

    subroutine ed_hamiltonian_init
        allocate(ek(nsite,2),vk(norb,nbath,2))
        ! input levels, hybridization
        ek = ek_in
        vk = vk_in
    end subroutine ed_hamiltonian_init

    ! generates the Hamiltonian for the sector in basis
    subroutine generate_hamiltonian(basis, H)
        type(basis_t), intent(in) :: basis
        double precision, intent(out) :: H(basis%nloc,basis%ntot)

        ! local variables
        integer(kind=kind_basis) :: basis_i, basis_j
        integer :: i,j,iorb,jorb,ibath,ispin, isite, jsite
        integer :: ni(2),nj(2),nb,sgntot, itmp

        if (.not.allocated(ek) .or. .not. allocated(vk)) then
            write(*,*) "Cluster Hamiltonian parameters are net initialized."
            call die("generate_hamiltonian", "ek, vk not allocated.")
        endif

        H = 0.0d0

        ! calculates <i|H|j>. 
        ! 1. find nonvanishing coefficients of H|j> = SUM_k |k> <k|H|j>
        ! 2. set H(i,j) = <i|H|j>
        jloop: do j=1,basis%nloc
            basis_j = ed_basis_get(basis,j)
            iorbloop: do iorb=1,norb
                if (BTEST(basis_j,iorb-1)) then
                    ni(1) = 1
                else
                    ni(1) = 0
                endif

                if (BTEST(basis_j,nsite+iorb-1)) then
                    ni(2) = 1 
                else
                    ni(2) = 0
                endif
                
                ! ==================================
                ! 1. onsite energies up/dn, 
                ! 2. intra orbital density-density
                ! ==================================
                H(j,j) = H(j,j) + (ek(iorb,1)-mu)*ni(1) &
                            + (ek(iorb,2)-mu)*ni(2) &
                            + U*ni(1)*ni(2)
                ! ==================================

                ! ==================================
                ! 3. inter orbital density-density
                ! ==================================
                do jorb=iorb+1,norb
                    ! n_j,up
                    if (BTEST(basis_j,jorb-1)) then
                        nj(1) = 1
                    else
                        nj(1) = 0
                    endif
                    ! n_j,dn
                    if (BTEST(basis_j,nsite+jorb-1)) then
                        nj(2) = 1
                    else
                        nj(2) = 0
                    endif

                    ! interorbital density-density
                    H(j,j) = H(j,j) + (Up-Jex)*(ni(1)*nj(1)+ni(2)*nj(2)) &
                                    + Up*(ni(1)*nj(2)+ni(2)*nj(1))
                enddo
                ! ==================================

                ! ==================================
                ! 4. hybridization
                ! ==================================
                do ispin=1,2
                    do ibath=1,nbath
                        isite = (ispin-1)*nsite+iorb
                        jsite = (ispin-1)*nsite+norb+ibath

                        if (ni(ispin).eq.0 .and. BTEST(basis_j,jsite-1)) then
                            ! c^+_{iorb,ispin} b_{ibath,ispin}
                            sgntot = sgn(basis_j,1,jsite-1)
                            basis_i = IBCLR(basis_j,jsite-1)

                            sgntot = sgntot*sgn(basis_i,1,isite-1)
                            basis_i = IBSET(basis_i,isite-1)

                            i = ed_basis_idx(basis, basis_i)

                            H(i,j) = H(i,j) + sgntot*vk(iorb,ibath,ispin) 

                        else if (ni(ispin).ne.0 .and. &
                                 .not.BTEST(basis_j,jsite-1)) then
                            ! b^+_{ibath,ispin} c_{iorb,ispin}
                            sgntot = sgn(basis_j,1,isite-1)
                            basis_i = IBCLR(basis_j,isite-1)

                            sgntot = sgntot*sgn(basis_i,1,jsite-1)
                            basis_i = IBSET(basis_i,jsite-1)

                            i = ed_basis_idx(basis, basis_i)

                            H(i,j) = H(i,j) + sgntot*vk(iorb,ibath,ispin) 

                        endif
                    enddo
                enddo
                ! ==================================

                ! ==================================
                ! 5. spin-flip & pair-hopping
                ! ==================================
                do jorb=1,norb
                    ! n_j,up
                    if (BTEST(basis_j,jorb-1)) then
                        nj(1) = 1
                    else
                        nj(1) = 0
                    endif
                    ! n_j,dn
                    if (BTEST(basis_j,nsite+jorb-1)) then
                        nj(2) = 1
                    else
                        nj(2) = 0
                    endif

                    ! spin-flip
                    ! -Jp c^+_{j,up} c_{j,dn} c^+_{i,dn} c_{i,up}
                    if (ni(1).ne.0.and.ni(2).eq.0.and.&
                        nj(1).eq.0.and.nj(2).ne.0) then

                        sgntot = sgn(basis_j,1,iorb-1)
                        basis_i = IBCLR(basis_j,iorb-1)

                        sgntot = sgntot*sgn(basis_i,1,nsite+iorb-1)
                        basis_i = IBSET(basis_i,nsite+iorb-1)

                        sgntot = sgntot*sgn(basis_i,1,nsite+jorb-1)
                        basis_i = IBCLR(basis_i,nsite+jorb-1)

                        sgntot = sgntot*sgn(basis_i,1,jorb-1)
                        basis_i = IBSET(basis_i,jorb-1)

                        i = ed_basis_idx(basis, basis_i)

                        H(i,j) = H(i,j) - sgntot*Jp

                    endif

                    ! pair-hopping
                    ! -Jp c^+_{j,up} c^+_{j,dn} c_{i,up} c_{i,dn}
                    if (ni(1).eq.0.and.ni(2).eq.0.and.&
                        nj(1).ne.0.and.nj(2).ne.0) then

                        sgntot = sgn(basis_j,1,nsite+iorb-1)
                        basis_i = IBCLR(basis_j,nsite+iorb-1)

                        sgntot = sgntot*sgn(basis_i,1,iorb-1)
                        basis_i = IBCLR(basis_i,iorb-1)

                        sgntot = sgntot*sgn(basis_i,1,nsite+jorb-1)
                        basis_i = IBSET(basis_i,nsite+jorb-1)

                        sgntot = sgntot*sgn(basis_i,1,jorb-1)
                        basis_i = IBSET(basis_i,jorb-1)

                        i = ed_basis_idx(basis, basis_i)

                        H(i,j) = H(i,j) - sgntot*Jp
                    endif
                enddo 
            enddo iorbloop

            ! ==================================
            ! 6. bath onsite
            ! ==================================
            do ibath=1,nbath
                if (BTEST(basis_j,norb+ibath-1)) then
                    H(j,j)=H(j,j)+ek(norb+ibath,1)
                endif
                if (BTEST(basis_j,nsite+norb+ibath-1)) then
                    H(j,j)=H(j,j)+ek(norb+ibath,nspin)
                endif
            enddo
        enddo jloop

    end subroutine generate_hamiltonian
    
end module ed_hamiltonian
