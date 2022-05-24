module input_data

    implicit none

    !Seção de declaração de variáveis globais de uso no programa
    integer(4)::ngln=2 !número de graus de liberdade por nós
    integer(4)::npe=10 !número de nós por elemento
    integer(4)::ep=0 !Estado Plano de Análise. EPT -> 0 EPD -> 1
    integer(4)::nnos,nmat,nesp,nel,nforca_nodal,nvinc !números de cada entrada
    real(8),dimension(:,:),allocatable::coord,mat,forca_nodal,vinc,forca_dist_elem,dist_elem
    integer(4),dimension(:,:),allocatable::prop,conec
    real(8),dimension(:),allocatable::esp
    real(8),dimension(:,:),allocatable::kglobal
    real(8),dimension(:),allocatable::fglobal,desl
    real(8),dimension(:,:),allocatable::tensao
    real(8),dimension(:,:),allocatable::ksieta

contains

    subroutine read_input(filename)

        implicit none
        
        !Seção de declaração de variáveis locais
        integer(4)::i,j,index,cont_vinc=0
        real(8),dimension(:,:),allocatable::dado_vinc
        character(20),dimension(:),allocatable::vinc_direc
        character(20)::key
        character(*),intent(in)::filename

        write(*,'(a)',advance='no')'Entrada de dados...'
        open(unit=20,file=filename)
        !Entrada das quantidades de cada um dos dados relevantes
        read(20,*)key,nnos !número de nós
        read(20,*)key,nmat !número de materiais distintos
        read(20,*)key,nesp !número de espessuras distintas
        read(20,*) !elemento da lista sem informação útil para entrada de dados
        read(20,*)
        read(20,*)key,nel  !número de elementos
        read(20,*)key,nforca_nodal !número de forças nodais aplicadas
        read(20,*)
        read(20,*)key,nvinc !úmero de vínculo nos nós

        !Alocação do vetores e matrizes 
        allocate(coord(nnos,ngln),conec(nel,npe),mat(nmat,2),esp(nesp),prop(nel,2))
        allocate(forca_nodal(nforca_nodal*2,3),dado_vinc(nvinc,2),vinc_direc(nvinc))
        allocate(forca_dist_elem(nel,3),dist_elem(nel,2))
        !allocate(kglobal(ngln*nnos,ngln*nnos),fglobal(ngln*nnos),desl(ngln*nnos))
        !allocate(tensao(nnos,4))

        read(20,*)
        read(20,*)

        !armazenamento de coordenadas do nós. padrão: [[x,y]_i]
        do i=1,nnos,1
            read(20,*)index,coord(i,1),coord(i,2)
        enddo

        read(20,*)
        read(20,*)

        !criação da matriz de materiais. padrão: [[E,ni]_i]
        do i=1,nmat,1
            read(20,*)index,mat(i,1),mat(i,2)
        enddo
    
        read(20,*)
        read(20,*)

        !criação do vetor de espessuras. padrão: [esp]_i
        do i=1,nesp,1
            read(20,*)index,esp(i)
        enddo

        read(20,*)
        read(20,*)

        !criação das matrizes de conectividade e de associação de propriedades
        !efetivamente a conectividade é a criação dos elementos finitos
        !padrão: conec=[[no 1, ... , npe]_i]
        !Padrão: prop=[[num_mat, num_esp]_i]
        do i=1,nel,1
            read(20,"(i10)",advance="no")index
            do j=1,npe,1
                read(20,"(i10)",advance="no")conec(i,j) !associa espessura ao i-ésimo elemento
                conec(i,j)=conec(i,j)+1
            enddo
            read(20,*)prop(i,1),prop(i,2) !associa material ao i-ésimo elemento
            prop(i,1)=prop(i,1)+1
            prop(i,2)=prop(i,2)+1
        enddo

        read(20,*)
        read(20,*)

        !Criação da matriz de forças nodais com associação aos nós
        !Padrão:[[num_no, dir(0->x ou 1->y), valor]_i]
        do i=1,nforca_nodal,1
            read(20,*)index,forca_nodal(i*2-1,1),forca_nodal(i*2-1,3),forca_nodal(i*2,3) 
            forca_nodal(i*2-1,2)=0
            forca_nodal(i*2,2)=1
            forca_nodal(i*2-1,1)=forca_nodal(i*2-1,1)+1
            forca_nodal(i*2,1)=forca_nodal(i*2-1,1)
        enddo

        read(20,*)
        read(20,*)

        !Guarda dados de vinculação
        do i=1,nvinc,1
            read(20,*) index, dado_vinc(i,1),vinc_direc(i),dado_vinc(i,2)
            dado_vinc(i,1)=dado_vinc(i,1)+1
        enddo

        close(unit=20)

        !Pré-processamento de vinculações
        !Contagem número de graus de liberdade restringidos totais 
        !para alocação da matriz de vinculação
        do i=1,nvinc,1
            if (vinc_direc(i)=='BOTH') then
                cont_vinc = cont_vinc+2
            else if ((vinc_direc(i)=='X').or.((vinc_direc(i)=='Y'))) then
                cont_vinc = cont_vinc+1
            endif
        enddo
        allocate(vinc(cont_vinc,3))!alocação de matriz de vinculações

        !criação da matriz de vinculações.
        !Padrão:[[num_no, dir(0->x ou 1->y), valor]_i]
        i=1
        do j=1,nvinc,1
            if (vinc_direc(j)=='BOTH') then
                vinc(i,1) = dado_vinc(j,1)
                vinc(i,2) = 0
                vinc(i,3) = dado_vinc(j,2)
                i=i+1
                vinc(i,1) = dado_vinc(j,1)
                vinc(i,2) = 1
                vinc(i,3) = dado_vinc(j,2)
                i=i+1
            else if (vinc_direc(j)=='X') then
                vinc(i,1) = dado_vinc(j,1)
                vinc(i,2) = 0
                vinc(i,3) = dado_vinc(j,2)
                i=i+1
            else
                vinc(i,1) = dado_vinc(j,1)
                vinc(i,2) = 1
                vinc(i,3) = dado_vinc(j,2)
                i=i+1
            endif
        enddo

        
        forca_dist_elem=reshape((/1.0d0,0.0d0,0.0d0/),(/1,3/))
        !criação de peso próprio em todos os elementos
        ! do i=1,nel
        !     forca_dist_elem(i,1)=i
        !     forca_dist_elem(i,2)=0 !direção x
        !     forca_dist_elem(i,3)=0.0 !direção y
        ! enddo

        !criação de matriz de cargas distribuídas por elemento padronizada
        !Se houver elementos sem cargas distribuídas em alguma das direções
        !os valores das cargas serão zero nessa direção
        !a matriz criada será: [[qx,qy]_i]
        dist_elem=0
        do i=1,size(forca_dist_elem,1)
            do j=1,ngln
                dist_elem(int(forca_dist_elem(i,1)),j)=forca_dist_elem(i,j+1)
            enddo
        enddo

        write(*,*)' Ok'
    end subroutine read_input    

end module input_data