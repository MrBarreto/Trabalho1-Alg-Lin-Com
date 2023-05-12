function power_method(matriz, tol)
    (n,m) = size(matriz)
    erro = Inf
    greater_eigenvalue = 1
    greater_eigenvector = ones(1, n)
    while erro > tol
        iter_vec = greater_eigenvector*matriz
        iter_eigenvalue = iter_vec[1]
        for i in 1:n
            if abs(iter_vec[i])> abs(iter_eigenvalue)
                iter_eigenvalue = iter_vec[i]
            end
        end
        for i in 1:n
            iter_vec[i] = iter_vec[i]/iter_eigenvalue
        end
        greater_eigenvector = iter_vec
        erro = (abs(iter_eigenvalue - greater_eigenvalue))/abs(iter_eigenvalue)
        greater_eigenvalue = iter_eigenvalue
    end
    return(greater_eigenvalue, greater_eigenvector)
end

function jacobi_eigenvalues(matriz, tol)