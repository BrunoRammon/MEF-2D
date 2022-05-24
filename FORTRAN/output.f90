module output_data
    use input_data
    implicit none

contains
    !Subrotina para geraçao de arquivo de saída
    subroutine print_results(filename)
        
        implicit none
        
        character(*),intent(in)::filename
        integer(4)::i,j

        write(*,'(a)',advance='no')'Impressão de resultados...'

        open(unit=21,file=filename)

        write(21,'(a)')'Resutados modelo MEF'
        write(21,'(a)')'nnos nel nlistas'
        write(21,'(a)')'#'
        write(21,'(3i8)')nnos,nel,5

        write(21,*)
        write(21,'(a)')'coord.x coord.y coord.z desl.x desl.y desl.z'
        write(21,'(a)')'#'
        do i = 1,nnos,1
            write(21,'(2f10.5,4f4.1)')coord(i,1),coord(i,2),0.0,0.0,0.0,0.0
        enddo

        write(21,*)
        write(21,'(a)')'tipo gr.aprox no1...no.npe grupo\n'
        write(21,'(a)')'#'
        do i = 1,nel,1
            write(21,'(i2,i2)',advance="no")2,3
            do j = 1,npe,1
                write(21,'(i5)',advance="no")conec(i,j)
            enddo
            write(21,'(i4)')0
        enddo

        write(21,*)
        write(21,'(a)')'listas'
        write(21,'(a)')'desl.x desl.y desl.z valor.cores'

        write(21,'(a)')'deslocamento.x'
        write(21,'(a)')'#'
        write(21,'(a)')'desl.x'
        do i = 1,nnos,1
            write(21,'(2f10.5,1f4.1,1f10.5)')desl(2*i-1),desl(2*i),0.0,desl(2*i-1)
        enddo

        write(21,'(a)')'deslocamento.y'
        write(21,'(a)')'#'
        write(21,'(a)')'desl.y'
        do i = 1,nnos,1
            write(21,'(2f10.5,1f4.1,1f10.5)')desl(2*i-1),desl(2*i),0.0,desl(2*i)
        enddo

        write(21,'(a)')'sigma.x'
        write(21,'(a)')'#'
        write(21,'(a)')'sigma.x'
        do i = 1,nnos,1
            write(21,'(2f10.5,1f4.1,1f15.5)')desl(2*i-1),desl(2*i),0.0,tensao(i,2)
        enddo

        write(21,'(a)')'sigma.y'
        write(21,'(a)')'#'
        write(21,'(a)')'sigma.y'
        do i = 1,nnos,1
            write(21,'(2f10.5,1f4.1,1f15.5)')desl(2*i-1),desl(2*i),0.0,tensao(i,3)
        enddo

        write(21,'(a)')'tau.xy'
        write(21,'(a)')'#'
        write(21,'(a)')'tau.xy'
        do i = 1,nnos,1
            write(21,'(2f10.5,1f4.1,1f15.5)')desl(2*i-1),desl(2*i),0.0,tensao(i,4)
        enddo

        close(21)
        write(*,*)' Ok'
    end subroutine 
end module output_data