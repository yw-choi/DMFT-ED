program main
    use dmft, only: dmft_init, dmft_loop, dmft_post_processing, dmft_finalize
    implicit none
    
    call dmft_init
    call dmft_loop
    call dmft_post_processing
    call dmft_finalize

end program main
