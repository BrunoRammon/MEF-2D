module calcular_tensoes

    use elemento_triangular
    use input_data

    implicit none

contains
    subroutine calc_tens()
        
        implicit none
        !variáveis locais
        integer(4)::iel,ino
        integer(4)::i,j
        real(8),dimension(npe*ngln)::eldesl
        real(8),dimension(npe)::fi,dfidksi,dfideta
        real(8),dimension(4,ngln*npe)::dfi
        real(8),dimension(2,2)::jac,jacinv
        real(8),dimension(2,4)::dx,dy
        real(8),dimension(4)::vec_dksideta
        real(8),dimension(2)::ddx,ddy
        real(8),dimension(3)::vec_eps, vec_sig
        real(8),dimension(3,3)::md
        real(8)::esp,ksi,eta,det
        real(8)::D11,D12,D22,D33

        write(*,'(a)',advance='no')'Cálculo de tensões... '
        allocate(ksieta(npe,2))
        allocate(tensao(nnos,4))

        !para elemento de 3 nós
        !ksieta = reshape((/real(8)::1.0,0.0,0.0,0.0,1.0,0.0/),(/npe,2/)) !coordenada de elementos de 3 nós

        !para elemento de 10 nós
        ksieta = transpose(reshape((/real(8)::&
                            3.0/3.0,0.0/3.0,2.0/3.0,1.0/3.0,&
                            1.0/3.0,2.0/3.0,0.0/3.0,3.0/3.0,&
                            2.0/3.0,0.0/3.0,1.0/3.0,1.0/3.0,&
                            0.0/3.0,2.0/3.0,1.0/3.0,0.0/3.0,&
                            0.0/3.0,1.0/3.0,0.0/3.0,0.0/3.0 &
                        /),(/2,10/)))!coordenadas de elementos de 10 nós
        
        tensao=0
        !Calculo da matriz de tensões média nos nós. padrão: [[numcontrib,σx,σy,τxy]_iel]
        do iel=1,nel,1 !varia de acordo com o número de elementos iel

            !cálculo de elementos da matriz constitutiva
            call matriz_constitutiva(iel,esp,D11,D12,D22,D33)

            !armazenamento dos deslocamentos do nós do elemento atual
            do i=1,npe,1
                do j=1,ngln,1
                    eldesl(ngln*i-1+j-1) = desl(ngln*conec(iel,i)-1+j-1)
                enddo
            enddo

            !loop para computar a contribuição de tensão de nó do elemento
            do ino=1,npe,1 !Varia de acordo com o número de nós no elemento

                ksi = ksieta(ino,1)
                eta = ksieta(ino,2)

                call fforma(ksi,eta,fi,dfidksi,dfideta) !chama a função de forma
                !salva as funções φ, dφ/dξ, dφ/dη
                !são avaliadas no ino-ésimo nó do elemento

                call matriz_DFI(dfidksi,dfideta,dfi) !calcula DFI do elemento atual 
                !relativo a coordenada do nó atual

                call matriz_jacobiana(iel,dfidksi,dfideta,jac) !calcula jacobiana do elemento atual
                !relativo a coordenada do nó atual
                jacinv = matinv2(jac)
                det = matdet2(jac)

                call matriz_DXDY(jacinv,dx,dy) 

                vec_dksideta = matmul(dfi,eldesl) ![∂u/∂ξ,∂v/∂ξ,∂u/∂η,∂v/∂η]

                ddx = matmul(dx,vec_dksideta) ![∂u/∂x,∂v/∂x]
                ddy = matmul(dy,vec_dksideta) ![∂u/∂y,∂v/∂y]

                vec_eps = (/ddx(1),ddy(2),0.5*(ddy(1)+ddx(2))/) ![εx,εy,εxy]

                md = reshape((/real(8)::D11,D12,0.0,D12,D22,0.0,0.0,0.0,2*D33/),(/3,3/)) !matriz constitutiva

                vec_sig = matmul(md,vec_eps) ![σx,σy,σxy]
                
                !contribuição dos nós de cada elemento
                tensao(conec(iel,ino),1) = tensao(conec(iel,ino),1) + 1 !número de contribuições para o ih_ésimo elemento
                tensao(conec(iel,ino),2) = tensao(conec(iel,ino),2) + vec_sig(1) !contribuição para σx do ih_ésimo nó do iel-ésimo elemento
                tensao(conec(iel,ino),3) = tensao(conec(iel,ino),3) + vec_sig(2) !contribuição para σy do ih_ésimo nó do iel-ésimo elemento
                tensao(conec(iel,ino),4) = tensao(conec(iel,ino),4) + vec_sig(3) !contribuição para σxy do ih_ésimo nó do iel-ésimo elemento
            enddo
        enddo
        !cálculo da tensão média de cada nó.
        tensao(:,2) = tensao(:,2)/tensao(:,1) !σx
        tensao(:,3) = tensao(:,3)/tensao(:,1) !σy
        tensao(:,4) = tensao(:,4)/tensao(:,1) !σxy
        write(*,*)'Ok!'
    end subroutine
end module calcular_tensoes