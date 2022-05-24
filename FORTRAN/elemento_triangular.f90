module elemento_triangular
    use input_data
    use math 
    implicit none   

contains

    !subrotina para os calculo dos elementos da matriz constitutiva do elemento finito
    subroutine matriz_constitutiva(num_el,h,D11,D12,D22,D33)

        implicit none
        !declaração dos parâmetros
        integer(4),intent(in)::num_el
        real(8),intent(out)::h
        real(8),intent(out)::D11,D22,D12,D33
        !variáveis locais
        real(8)::young,ni
    
        young = mat(prop(num_el,1),1) !módulo de young do iel-ésimo elemento
        ni = mat(prop(num_el,1),2) !coeficiente de poisson do iel-ésimo elemento
        h = esp(prop(num_el,2)) !espessura do iel-ésimo elemento
        !decide sobre o uso de expressões para o EPT ou para o EPD
        !e monta as matrizes Mxx,Myy,Mxy,Myx do atual elemento do loop
        if (ep==0) then !EPT
            D11 = young/(1-ni**2)
            D22 = young/(1-ni**2)
            D12 = ni*young/(1-ni**2)
        else !EPD
            D11 = (1-ni)*young/((1+ni)*(1-ni**2))
            D22 = (1-ni)*young/((1+ni)*(1-ni**2))
            D12 = ni*young/((1+ni)*(1-ni**2))
        endif
        D33 = young/(2*(1+ni))
        
    end subroutine

    !subrotina para função de forma e suas derivadas do elemento finito
    subroutine fforma(ksi,eta,fi,dfidksi,dfideta)
        
        implicit none
        !Declaração dos parâmetros dos parãmetros
        real(8),intent(in)::ksi,eta
        real(8),dimension(npe),intent(out)::fi,dfidksi,dfideta
        !Para elementos de 10 nós
        fi=(/ real(8)::&
            1.0*ksi - 4.5*ksi**2 + 4.5*ksi**3,&
            -4.5*eta*ksi + 13.5*eta*ksi**2,&
            -4.5*eta*ksi + 13.5*eta**2*ksi,&
            1.0*eta - 4.5*eta**2 + 4.5*eta**3,&
            -4.5*ksi + 4.5*eta*ksi + 18.0*ksi**2 - 13.5*eta*ksi**2 - 13.5*ksi**3,&
            27.0*eta*ksi - 27.0*eta**2*ksi - 27.0*eta*ksi**2,&
            -4.5*eta + 18.0*eta**2 - 13.5*eta**3 + 4.5*eta*ksi - 13.5*eta**2*ksi,&
            9.0*ksi - 22.5*eta*ksi + 13.5*eta**2*ksi - 22.5*ksi**2 + 27.0*eta*ksi**2 + 13.5*ksi**3,&
            9.0*eta - 22.5*eta**2 + 13.5*eta**3 - 22.5*eta*ksi + 27.0*eta**2*ksi + 13.5*eta*ksi**2,&
            1.0 - 5.5*eta + 9.0*eta**2 - 4.5*eta**3 - 5.5*ksi + 18.0*eta*ksi - 13.5*eta**2*ksi &  
            + 9.0*ksi**2 - 13.5*eta*ksi**2 - 4.5*ksi**3 &
            /)

        dfidksi = (/ real(8)::&
                    13.5*ksi**2 - 9.0*ksi + 1.0, &
                    27.0*eta*ksi - 4.5*eta, &
                    13.5*eta**2 - 4.5*eta, &
                    0,&
                    -27.0*eta*ksi + 4.5*eta - 40.5*ksi**2 + 36.0*ksi - 4.5,&
                    -27.0*eta**2 - 54.0*eta*ksi + 27.0*eta,&
                    -13.5*eta**2 + 4.5*eta,&
                    13.5*eta**2 + 54.0*eta*ksi - 22.5*eta + 40.5*ksi**2 - 45.0*ksi + 9.0,&
                    27.0*eta**2 + 27.0*eta*ksi - 22.5*eta,&
                    -13.5*eta**2 - 27.0*eta*ksi + 18.0*eta - 13.5*ksi**2 + 18.0*ksi - 5.5 &
                  /)
        dfideta = (/real(8):: &
                    0,&
                    13.5*ksi**2 - 4.5*ksi,&
                    27.0*eta*ksi - 4.5*ksi,&
                    13.5*eta**2 - 9.0*eta + 1.0,&
                    -13.5*ksi**2 + 4.5*ksi,&
                    -54.0*eta*ksi - 27.0*ksi**2 + 27.0*ksi,&
                    -40.5*eta**2 - 27.0*eta*ksi + 36.0*eta + 4.5*ksi - 4.5,&
                    27.0*eta*ksi + 27.0*ksi**2 - 22.5*ksi,&
                    40.5*eta**2 + 54.0*eta*ksi - 45.0*eta + 13.5*ksi**2 - 22.5*ksi + 9.0,&
                    -13.5*eta**2 - 27.0*eta*ksi + 18.0*eta - 13.5*ksi**2 + 18.0*ksi - 5.5 &
                  /)

        !Para três nós
        ! fi = (/real(8)::ksi,eta,1-ksi-eta/)
        ! dfidksi = (/real(8)::1,0,-1/)
        ! dfideta = (/real(8)::0,1,-1/)
    end subroutine

    !subrotina para calculo da DFI do elemento finito
    subroutine matriz_DFI(dfidksi,dfideta,dfi)

        implicit none

        !declaração de parâmetros
        real(8),dimension(npe),intent(in)::dfidksi,dfideta
        real(8),dimension(4,ngln*npe),intent(out)::dfi
        !declaração de variáveis locais
        integer(4)::i
        do i=1,npe
            dfi(1,ngln*i-1) = dfidksi(i)
            dfi(2,ngln*i) = dfidksi(i)
            dfi(3,ngln*i-1) = dfideta(i)
            dfi(4,ngln*i) = dfideta(i)
        enddo
    end subroutine

    !subrotina para calculo da jacobiana do elemento finito
    subroutine matriz_jacobiana(num_elem,dfidksi,dfideta,jac)

        implicit none
        !declaração de parâmetros
        integer(4),intent(in)::num_elem
        real(8),dimension(npe),intent(in)::dfidksi,dfideta
        real(8),dimension(2,2),intent(out)::jac
        !declaração de variáveis locais
        integer(4)::i

        jac=0
        do i=1,npe,1
            !Com o número do nó dado pela conectividade
            !busca-se na matriz coord a coordenada x(coluna 1) ou y(coluna 2)
            !desse i-ésimo do elemento atual da iteralção.
            jac(1,1) = jac(1,1) + dfidksi(i)*coord(conec(num_elem,i),1)
            jac(1,2) = jac(1,2) + dfideta(i)*coord(conec(num_elem,i),1)
            jac(2,1) = jac(2,1) + dfidksi(i)*coord(conec(num_elem,i),2)
            jac(2,2) = jac(2,2) + dfideta(i)*coord(conec(num_elem,i),2)
        enddo
    end subroutine

    !subrotina para calculo das matrizes DX e DY
    subroutine matriz_DXDY(jacinv,dx,dy)

        implicit none
        !declaração de parâmetros
        real(8),dimension(2,2),intent(in)::jacinv
        real(8),dimension(2,4),intent(out)::dx,dy

        dx=0
        dy=0
        dx(1,1) = jacinv(1,1)
        dx(1,3) = jacinv(2,1)
        dx(2,2) = jacinv(1,1)
        dx(2,4) = jacinv(2,1)      
        dy(1,1) = jacinv(1,2)
        dy(1,3) = jacinv(2,2)
        dy(2,2) = jacinv(1,2)
        dy(2,4) = jacinv(2,2)
    end subroutine

    !subrotina para calculo de matriz de rigidez global e vetor de força nodal equivalente locais 
    !do elemento finito por integração númerica 
    subroutine sistema_local(iel,klocal,flocal)

        implicit none

        !declaração de parâmetros
        integer(4),intent(in)::iel
        real(8),dimension(ngln*npe,ngln*npe),intent(out)::klocal
        real(8),dimension(ngln*npe),intent(out)::flocal
        !declaração de variáveis locais
        real(8),dimension(ngln*npe)::fiq
        integer(4),parameter::nph=12
        integer(4)::ih,i,j !variáveis para uso em iterações de laços 
        real(8),dimension(npe)::fi,dfidksi,dfideta
        real(8),dimension(4,ngln*npe)::dfi
        real(8),dimension(2,2)::jac,jacinv
        real(8),dimension(2,4)::dx,dy
        real(8),dimension(3,nph)::hammer
        real(8)::ksi,eta,peso,det
        real(8),dimension(2,2)::mxx,myy,mxy,myx
        real(8)::D11,D12,D22,D33
        real(8)::esp

        !Definição da matris de ponto de hammer [[ξ,η,peso]_ihammer]
        hammer=transpose(reshape((/ &
                 0.501426509658179, 0.249286745170910, 0.249286745170910,&
                 0.873821971016996, 0.063089014491502, 0.063089014491502,&
                 0.053145049844816, 0.310352451033785, 0.636502499121399,&
                 0.310352451033785, 0.636502499121399, 0.053145049844816,&

                 0.249286745170910, 0.249286745170910, 0.501426509658179,&
                 0.063089014491502, 0.063089014491502, 0.873821971016996,&
                 0.310352451033785, 0.636502499121399, 0.053145049844816,&
                 0.053145049844816, 0.310352451033785, 0.636502499121399,&

                 0.116786275726379/2.0, 0.116786275726379/2.0, 0.116786275726379/2.0,&
                 0.050844906370207/2.0, 0.050844906370207/2.0, 0.050844906370207/2.0,&
                 0.082851075618374/2.0, 0.082851075618374/2.0, 0.082851075618374/2.0,&
                 0.082851075618374/2.0, 0.082851075618374/2.0, 0.082851075618374/2.0 /),&
               (/12,3/)))
        !calculo de elementos da matriz constitutiva
        call matriz_constitutiva(iel,esp,D11,D12,D22,D33)
        !Definição de matrizes constitutivas: mxx, myy, mxy, myx
        mxx(1,1) = D11
        mxx(1,2) = 0
        mxx(2,1) = 0
        mxx(2,2) = D33
        myy(1,1) = D33
        myy(1,2) = 0
        myy(2,1) = 0
        myy(2,2) = D22
        mxy(1,1) = 0
        mxy(1,2) = D12
        mxy(2,1) = D33
        mxy(2,2) = 0
        myx(1,1) = 0
        myx(1,2) = D33
        myx(2,1) = D12
        myx(2,2) = 0

        klocal=0
        flocal=0
        fiq=0
        do ih=1,nph,1 !varia de acorco com o número de pontos de integração
            ksi = hammer(1,ih) !coordenada ξ do ponto de integração atual do loop
            eta = hammer(2,ih) !coordenada η do ponto de integração atual do loop
            peso = hammer(3,ih)  !peso do ponto de integração atual do loop 


            call fforma(ksi,eta,fi,dfidksi,dfideta)
            !salva as funções φ, dφ/dξ, dφ/dη
            !são avaliadas no ponto de integração atual do loop

            call matriz_DFI(dfidksi,dfideta,dfi)

            call matriz_jacobiana(iel,dfidksi,dfideta,jac)

            
            jacinv = matinv2(jac)!calculo da inversa da jacobiana
            det = matdet2(jac)!calculo do determinante da jacobiana

            call matriz_DXDY(jacinv,dx,dy)

            !contribuição do ponto de integração atual para a matriz de rigidez do elemento considerado
            klocal = klocal + esp * det * peso * matmul(transpose(dfi),&
                                                    matmul(&
                                                        matmul(transpose(dx),matmul(mxx,dx))+&
                                                        matmul(transpose(dy),matmul(myy,dy))+&
                                                        matmul(transpose(dx),matmul(mxy,dy))+&
                                                        matmul(transpose(dy),matmul(myx,dx)),&
                                                        dfi&
                                                    )&
                                                )
                                                
            !criação do vetor fi*q (produto da fi pela carga distribuida 'q')
            fiq=0
            do i=1,npe,1
                do j=1,ngln,1
                    fiq(ngln*i-1+j-1)=fi(i)*dist_elem(iel,j)
                enddo
            enddo
            !contribuição do ponto de integração atual para força distribuída 
            !ao vetor de forças nodal-equivalentes do elemento considerado 
            flocal=flocal + det * peso * fiq

        enddo
    end subroutine

end module elemento_triangular