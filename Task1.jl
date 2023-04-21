function LU_decomposition(Matriz)
   (n,m) = size(Matriz)
    for k = 1: n - 1
        
        let 
            for i in k + 1: n
                Matriz[i, k] = Matriz[i, k]/Matriz[k, k]
            end
        end
        
        let
            for j = k + 1: n
                for i = k + 1: n
                    Matriz[i, j] = Matriz[i, j] - Matriz[i, k]*Matriz[k, j]
                end
            end
        end
    end
    return Matriz
end

function system_solver(Matriz, vetor)
    (n,m) = size(Matriz)
    first_vec = zeros(n)
    results_vec = zeros(n)
    first_vec[1] = vetor[1]
    
    for i = 2: n
        subtraction = vetor[i]
        for j = 1: i-1
            subtraction -= Matriz[i, j]*first_vec[j]
        end
        first_vec[i] = subtraction
    end
    
    results_vec[n] = first_vec[n]/Matriz[n, n]
    
    for k = n - 1:-1:1
        subraction = first_vec[k]
        for l = k + 1:n
            subraction -= Matriz[k, l]*results_vec[l]
        end
        results_vec[k] = subraction/Matriz[k, k]
    end

    return results_vec
end

matriz=[1.0 2.0 2.0 
        4.0 4.0 2.0 
        4.0 6.0 4.0]
b = [3,6,10]
matriz = LU_decomposition(matriz)
resultado = system_solver(matriz, b)
print("$resultado")



