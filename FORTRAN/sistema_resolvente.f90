module sistema_resolvente
    use input_data
    use math 
    use elemento_triangular
    implicit none

contains

    !subrotina para montagem do sistema de rigidez global
    subroutine sistema_global()
        implicit none
        
        !variáveis locais
        real(8),dimension(ngln*npe,ngln*npe)::klocal
        real(8),dimension(ngln*npe)::flocal
        integer(4)::iel,il,jl,ic,jc !variáveis para uso no loop
        integer(4)::no,direc
        
        write(*,'(a)',advance='no')'Montagem do sistema global...'
        !Alocação da matriz de rigidez e do vetor de forças globais
        allocate(kglobal(ngln*nnos,ngln*nnos),fglobal(ngln*nnos),desl(ngln*nnos))
        kglobal=0
        fglobal=0
        !loop para calculo da matriz de rigidez global    
        do iel=1,nel,1
            
            !calculo da matriz de rigidez e vetor de força nodal locais do elemento atual do loop
            call sistema_local(iel,klocal,flocal)

            !loop para contribuição do elemento atual na matriz de rigidez e vetor de força nodal globais
            do il=1,npe,1 !varia nas linha
                do jl=1,ngln,1 !varia de acordo com o num de grau de liberdade por nó
                    !contribuição de forças
                    fglobal(ngln*conec(iel,il)-1+jl-1)=&
                    fglobal(ngln*conec(iel,il)-1+jl-1)&
                    +flocal(ngln*il-1+jl-1)
                    do ic=1,npe,1 !varia nas colunas
                        do jc=1,ngln,1
                            !contribuição de rigidez
                            kglobal(ngln*conec(iel,il)-1+jl-1,ngln*conec(iel,ic)-1+jc-1)=&
                            kglobal(ngln*conec(iel,il)-1+jl-1,ngln*conec(iel,ic)-1+jc-1)&
                            +klocal(ngln*il-1+jl-1,ngln*ic-1+jc-1)
                        enddo
                    enddo
                enddo
            enddo

        enddo
        !contribuição das forças nodais aplicadas (dados de entrada)
        do il=1,size(forca_nodal,1),1
            no = int(forca_nodal(il,1))
            direc = int(forca_nodal(il,2))
            fglobal(ngln*no-1+direc)=fglobal(ngln*no-1+direc)+forca_nodal(il,3)
        enddo

        call sist_org()!organiza a matriz global pela técnica de 0 e 1
        write(*,*)' Ok!'
    end subroutine

    !subrotina para criação colocação de 0 e 1
    subroutine sist_org()

        implicit none

        integer(4)::i,no,direc

        
        do i=1,size(vinc,1),1
            no=int(vinc(i,1))
            direc = int(vinc(i,2))
            fglobal=fglobal-vinc(i,3)*kglobal(:,ngln*no-1+direc)
            fglobal(ngln*no-1+direc)=vinc(i,3) !atribui o valor do i-ésimo deslocamento 
            !conhecido à sua posição correspondente no vetor de forças
            kglobal(ngln*no-1+direc,:)=0.0 !zera linha completa associada ao i-ésimo
            !grau de liberdade com deslocamento conhecido
            kglobal(:,ngln*no-1+direc)=0.0 !zera coluna completa associada ao i-ésimo
            !grau de liberdade com deslocamento conhecido
            kglobal(ngln*no-1+direc,ngln*no-1+direc)=1.0 !iguala a 1 o elemento da 
            !diagonal associado ao ivinc-ésimo grau de liberdade com deslocamento conhecido
        enddo

    end subroutine

    !subrotina para resolução do sistema resolvente. Obtenção dos deslocamentos 
    subroutine solve_system_of_equation()

        implicit none
        integer(4)::pivot(nnos*2),ok
        real(4) TpInit1, TpA
        write(*,'(a)',advance='no')'Resolucao do sistema...'
        TpInit1 = SECNDS(0.0)
        !desl = solve_cg(kglobal,fglobal)
        call dgesv(nnos*2,1,kglobal,nnos*2,pivot,fglobal,nnos*2,ok)
        desl=fglobal
        TpA = SECNDS(TpInit1)
        deallocate(kglobal)
        deallocate(fglobal)
        write(*,'(a,f10.5)')' Ok! Tempo decorrido [s]: ',TpA
    end subroutine solve_system_of_equation

end module sistema_resolvente