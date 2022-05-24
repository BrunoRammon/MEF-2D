program main

    use input_data
    use sistema_resolvente
    use calcular_tensoes
    use output_data

    implicit none

    character(50)::nome_entrada
    character(56)::nome_saida
    !RETIRAR
    !integer(4)::i,j

    nome_entrada="validacaopy.txt"
    nome_saida="saida_"//nome_entrada
    
    

    call read_input(nome_entrada)
    
    

    
    call sistema_global()


    call solve_system_of_equation()


    call calc_tens()

   call print_results(nome_saida)

end program main