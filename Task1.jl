function matrix_reader(archive)
    f = open(archive, "r")
    row = []
    bigvector = []
    nelements = 0
    while !eof(f)
        linha = chomp(readline(f))
        row = parse.(Float64, split(linha))
        nelements = length(row)
        bigvector = [bigvector;row]
    end

    matrix = zeros(nelements, nelements)
    for i = 1:10
        for j = 1:10
            matrix[i, j] = bigvector[(i-1)*nelements + j]
        end
    end

    close(f)
    return matrix
end

function vector_reader(archive)
    f = open(archive, "r")
    row = []
    bigvector = []
    while !eof(f)
        linha = chomp(readline(f))
        row = parse.(Float64, split(linha))
       push!(bigvector, row[1])
    end
    return bigvector
end

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

function cholesky_decomposition(Matriz)
    (n,m) = size(Matriz)
    for i = 1:n
        soma = 0
        for k = 1:i-1
            soma += (Matriz[i, k])^2
        end
        Matriz[i, i] = sqrt(Matriz[i, i] - soma)
        for j = i+1:n
            soma2 = 0
            for k = 1:i-1
                soma2 += Matriz[i, k]*Matriz[j, k]
            end
            Matriz[j, i] = (1/Matriz[i, i])*(Matriz[i, j] - soma2)
        end
    end
    for i = 1:n
        for l = i:n
            Matriz[i, l] = Matriz[l, i]
        end
    end
    return Matriz
end

function system_solver(Matriz, vetor, method)
    (n,m) = size(Matriz)
    first_vec = zeros(n)
    results_vec = zeros(n)
    if method == 1
        divisor = 1
    else
        divisor = Matriz[1, 1]
    end
    first_vec[1] = vetor[1]/divisor
    
    for i = 2: n
        subtraction = vetor[i]
        for j = 1: i-1
            subtraction -= Matriz[i, j]*first_vec[j]
        end
        if method == 1
            first_vec[i] = subtraction
        else
            first_vec[i] = subtraction/Matriz[i, i]
        end
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
#=
matriz=[4.0 12.0 -16.0 
        12.0 37.0 -43.0 
        -16.0 -43.0 98.0]
b = [3,6,10]
matriz = LU_decomposition(matriz)
resultado = system_solver(matriz, b)
print("$resultado")
=#