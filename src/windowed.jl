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
    for i=abs(dindex):n_fullblocks-1
        ix = i*bsize+brange
        iy = i*bsize+brange
        if dindex < 0
            ix += -dindex*bsize
        elseif dindex > 0
            iy += dindex*bsize
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
    ret_ex = ex
    if ex.head == :(=)
        left, right = (ex.args...)
        return :($(esc(left)) = @window $ws $right)
    end
    if ex.head == :call
        f = ex.args[1]
        if f == :crossrecurrencematrix
            ex.args[1] = :blockrecurrencematrix
            insert!(ex.args, 4, ws)
            insert!(ex.args, 5, -1)
            exd_lower  = :(bl = $ex)
            ex.args[5] = 0
            exd_center = :(bd = $ex)
            ex.args[5] = 1
            exd_upper  = :(bu = $ex)
            ex_ret = :(rmat = bl|bd|bu)
        elseif f == :recurrencematrix
            ex.args[1] = :blockrecurrencematrix
            insert!(ex.args, 3, ex.args[2])
            insert!(ex.args, 4, ws)
            insert!(ex.args, 5, -1)
            exd_lower  = :(bl = $ex)
            ex.args[5] = 0
            exd_center = :(bd = $ex)
            ex_ret = :(rmat = bl|bd|bl')
        end
        Throw error if it is not a valid function
        error("$(string(ex.args[1])) is not a valid function for windowing")
    end
    # Throw error if it didn't return
    error("Invalid expression for windowing")
end
