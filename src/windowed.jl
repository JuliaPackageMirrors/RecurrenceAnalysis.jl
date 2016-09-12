matrix_funs2 = [
    :crossdistancematrix,
    :crossrecurrencematrix,
    :jointrecurrencematrix
    ]
    
matrix_funs1 = Dict(
    :recurrencematrix => :crossrecurrencematrix,
    :distancematrix => :crossdistancematrix)

rqa_funs = [
    :recurrencerate,
    :determinism,
    :avgdiag,
    :maxdiag,
    :divergence,
    :entropy,
    :trend,
    :laminarity,
    :trappingtime,
    :maxvert
    ]

macro window(ws, ex)
    # Expression can be of type a = f(x...)
    ret_ex = ex
    if ex.head == :(=)
        left, right = (ex.args...)
        return :($(esc(left)) = @window $ws $right)
    end
    if ex.head == :call
        f = ex.args[1]
        x = ex.args[2]
        # Iteration of RQA functions
        fi = findfirst(rqa_funs, f)
        if fi > 0
            submat = :($x[i+w,i+w])
            ex.args[2] = submat
            ret_ex = quote
                w = 1:$ws
                nw = size($x,1) - $ws
                [$ex for i=0:nw]
            end
            return ret_ex
        end
        if haskey(matrix_funs1, f)
            ex.args[1] = matrix_funs1[f]
            insert!(ex.args, 2, x)
            return :(@window $ws $ex)
        end
        # Iteration of matrix construction (2 inputs)
        fi = findfirst(matrix_funs2, f)
        if fi > 0
            y = ex.args[3]
            ret_ex = Array(Expr,6)
            # Prologue
            ret_ex[1] = quote
                n = size($x,1)
                m = size($y,1)
                mn = min(m,n)
                retmat = map(Bool,spzeros(m,n))
                w = 1:$ws
                max_i = div(mn, $ws)
            end
            # Main diagonal
            ex.args[2] = :($x[i*($ws)+w,:])
            ex.args[3] = :($y[i*($ws)+w,:])
            ret_ex[2] = quote
                for i=0:max_i-1
                    retmat[i*($ws)+w,i*($ws)+w] = $ex
                end
            end
            ex.args[2] = :($x[max_i*($ws)+1:m,:])
            ex.args[3] = :($y[max_i*($ws)+1:n,:])
            ret_ex[3] = quote
                retmat[max_i*($ws)+1:m,max_i*($ws)+1:n] = $ex
            end
            # Upper diagonal
            ex.args[2] = :($x[(i-1)*($ws)+w,:])
            ex.args[3] = :($y[i*($ws)+w,:])
            ret_ex[4] = quote
                for i=1:max_i-1
                    retmat[(i-1)*$ws+w,i*$ws+w] = $ex
                end
            end
            ex.args[2] = :($x[(max_i-1)*($ws)+w,:])
            ex.args[3] = :($y[max_i*($ws)+1:n,:])
            ret_ex[5] = quote
                retmat[(max_i-1)*($ws)+w,max_i*($ws)+1:n] = $ex
            end
            # Lower diagonal
            ex.args[2] = :($x[i*($ws)+w,:])
            ex.args[3] = :($y[(i-1)*($ws)+w,:])
            ret_ex[6] = quote
                for i=1:max_i-1
                    retmat[i*$ws+w,(i-1)*$ws+w] = $ex
                end
            end
            ex.args[2] = :($x[max_i*($ws)+1:m,:])
            ex.args[3] = :($y[(max_i-1)*($ws)+w,:])
            ret_ex[5] = quote
                retmat[max_i*($ws)+1:m,(max_i-1)*($ws)+w] = $ex
            end
        end
        Throw error if it is not a valid function
        error("$(string(ex.args[1])) is not a valid function for windowing")
    end
    # Throw error if it didn't return
    error("Invalid expression for windowing")
end
