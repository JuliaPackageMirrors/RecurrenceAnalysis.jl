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

function blockrecurrencematrix(x, y, bsize, dindex, vargs...; kwargs...)
    n = min(size(x,1), size(y,1))
    rmat = map(Bool, spzeros(n,n))
    brange = 1:bsize
    n_fullblocks = div(n, bsize)
    ix = iy = 0:0
    for i=abs(dindex):n_fullblocks-1
        ix = i*bsize+brange
        iy = i*bsize+brange
        if dindex < 0
            iy -= -dindex*bsize
        elseif dindex > 0
            ix -= dindex*bsize
        end
        rmat[ix,iy] = crossrecurrencematrix(x[ix,:], y[iy,:], vargs...; kwargs...)
    end
    ix1 = ix[end]+1
    iy1 = iy[end]+1
    ix2 = min(ix1+bsize-1, n)
    iy2 = min(iy1+bsize-1, n)
    rmat[ix1:ix2, iy1:iy2] = crossrecurrencematrix(x[ix1:ix2,:], y[iy1:iy2,:], vargs...; kwargs...)
    rmat
end

macro window(ws, ex)
    # Expression can be of type a = f(x...)
    if ex.head == :(=)
        left, right = (ex.args...)
        return :($(esc(left)) = @window $ws $right)
    end
    if ex.head == :call
        f = ex.args[1]
        # Iteration of RQA functions
        fi = findfirst(rqa_funs, f)
        if fi > 0
            x = ex.args[2]
            submat = :($x[i+w,i+w])
            ex.args[2] = submat
            ret_ex = quote
                w = 1:$ws
                nw = size($x,1) - $ws
                [$ex for i=0:nw]
            end
            return ret_ex
        end
        if f == :crossrecurrencematrix
            ex.args[1] = :blockrecurrencematrix
            insert!(ex.args, 4, ws)
            insert!(ex.args, 5, -1)
            exd_lower  = :(bl = $(parse(string(ex))))
            ex.args[5] = 0
            exd_center = :(bd = $(parse(string(ex))))
            ex.args[5] = 1
            exd_upper  = :(bu = $(parse(string(ex))))
            exd_comp = :(bl|bd|bu)
            ret_ex = quote
                $exd_lower
                $exd_center
                $exd_upper
                $exd_comp
            end
            return ret_ex
        elseif f == :recurrencematrix
            ex.args[1] = :blockrecurrencematrix
            insert!(ex.args, 3, ex.args[2])
            insert!(ex.args, 4, ws)
            insert!(ex.args, 5, -1)
            exd_lower  = :(bl = $(parse(string(ex))))
            ex.args[5] = 0
            exd_center = :(bd = $(parse(string(ex))))
            exd_comp = :(bl|bd|bl')
            ret_ex = quote
                $exd_lower
                $exd_center
                $exd_comp
            end
            return ret_ex
        end
        # Throw error if it is not a valid function
        error("$(string(ex.args[1])) is not a valid function for windowing")
    end
    # Throw error if it didn't return
    error("Invalid expression for windowing")
end
