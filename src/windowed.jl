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

"""
    @windowed(f(x,...), width)
    @windowed(f(x,...); width, step=1)
    
Calculate windowed RQA parameters with a given window width.

`f(x,...)` may be any call to RQA functions (e.g. `recurrencerate`, `determinism`, etc.),
with `x` being a named variable that designates the recurrence matrix
(do not use 'in-place' calculations of the recurrence matrix).
The results are returned in a vector with one value for each position of the window.
By default the window moves at one-point intervals, but a longer `step` length
may be specified, together with the window `width`, by declaring those options as keyword arguments.

This macro may be also used with recurrence matrix constructors
(`recurrencematrix`, `crossrecurrencematrix`, `jointrecurrencematrix`),
to create 'incomplete' matrices that are suitable for such windowed RQA.
The values of the resulting matrix in the diagonals within the window width will
be equal to those obtained without the `@windowed` macro, if the distances are
not scaled (using the option `scale=1`, see `?recurrencematrix`).
Outside the window width, the values of the recurrence matrix will be undefined (mostly zero).
"""
macro windowed(ex, options...)
    # Expression can be of type a = f(x...)
    if ex.head == :(=)
        left, right = (ex.args...)
        return :($(esc(left)) = @windowed($ex,$(options...)))
    end
    # Parse options
    dict_op = Dict{Symbol,Any}(:step=>1)
    if length(options) == 1
        dict_op[:width] = options[1]
    end
    for op in options
        if typeof(op)<:Expr && in(op.head, [:(=), :kw])
            dict_op[op.args[1]] = op.args[2]
        end
    end
    if ex.head == :call
        f = ex.args[1]
        # Iteration of RQA functions
        # fun(x,...) => [fun(x[i+w,i+w]) for i=0:s:nw]
        if in(f, rqa_funs)
            x = ex.args[2]
            submat = :($x[i+w,i+w])
            ex.args[2] = submat
            ret_ex = quote
                w = 1:$(dict_op[:width])
                s = $(dict_op[:step])
                nw = size($x,1) - $(dict_op[:width])
                [$ex for i=0:s:nw]
            end
            return ret_ex
        end
        if f == :crossrecurrencematrix
            # blockrecurrencematrix(x,y,width,d,...) with d=-1,0,1
            ex.args[1] = :blockrecurrencematrix
            insert!(ex.args, 4, dict_op[:width])
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
            # blockrecurrencematrix(x,x,width,d,...) with d=-1,0
            ex.args[1] = :blockrecurrencematrix
            insert!(ex.args, 3, ex.args[2])
            insert!(ex.args, 4, dict_op[:width])
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
        elseif f == :jointrecurrencematrix
            x_string = string(ex.args[2])
            y_string = string(ex.args[3])
            ex.args[1] = :recurrencematrix
            deleteat!(ex.args, 3)
            ex.args[2] = parse(x_string)
            ex_rmx = :(rm1 = @windowed($(parse(string(ex))),width=$(dict_op[:width])))
            ex.args[2] = parse(y_string)
            ex_rmy = :(rm2 = @windowed($(parse(string(ex))),width=$(dict_op[:width])))
            ex_comp = :(rm1 .* rm2)
            ret_ex = quote
                $ex_rmx
                $ex_rmy
                $ex_comp
            end
            return ret_ex
        end
        # Throw error if it is not a valid function
        error("$(string(ex.args[1])) is not a valid function for windowing")
    end
    # Throw error if it didn't return
    error("Invalid expression for windowing")
end