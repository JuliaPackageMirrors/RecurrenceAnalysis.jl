# Diagonal histograms
function dhist_windowed(x::AbstractMatrix{Bool}, winsize::Integer; kwargs...)
    xw = x[1:winsize,1:winsize]
    n = minimum(size(x))
    nw = n-winsize+1
    dh, fh, lh = dhist_extended(xw; kwargs...)
    # Iterate
    for i=1:nw-1
        # Update the number of diagonals
        # Shift down on the basis of first line histogram
        # Shift up on the basis of the last line histogram
        dh .-= fh+lh
        dh[1:end-1] += fh[2:end]
        dh[2:end] += lh[1:end-1]
        (lh[end] > 0) && push!(dh,lh[end])
        (dh[end] <= 0) && pop!(dh)
    end
    dh
end

function dhist_extended(x::AbstractMatrix{Bool}; theiler::Integer=1, kwargs...)
    theiler < 0 && error("Theiler window length must be greater than 0")
    bins = [0]
    firstbins = [0] #W
    lastbins = [0] #W
    firstsingle = 0 #W
    lastsingle = 0 #W
    nbins = 1
    current_diag = 0
    m, n = size(x)
    # Iterate over diagonals - excluding LOI and Theiler window
    # If the matrix is symmetric, examine only the upper triangle
    diag_collection = collect(theiler:n-2)
    xsym = issym(x) # compat!
    !xsym && prepend!(diag_collection, collect(-(m-2):-max(theiler,1)))
    for d in diag_collection
        increment = (xsym && d > 0) ? 2 : 1
        previous_cell = false
        firstdiag = true #W
        first_c = max(1, d+1)
        last_c = min(n, m+d)
        for c = first_c:last_c
            # Operate on existing diagonal line
            if previous_cell
                # Extend length with current cell
                (extend = x[c-d,c]) && (current_diag += 1)
                # If arrived to the end of a line
                # add the current length to the corresponding bin
                if (!extend || (c == last_c)) && (current_diag > 0)
                    # Append new positions to the bins if needed
                    if current_diag > nbins
                        append!(bins, zeros(current_diag - nbins))
                        append!(firstbins, zeros(current_diag - nbins)) #W
                        append!(lastbins, zeros(current_diag - nbins)) #W
                        nbins = current_diag
                    end
                    bins[current_diag] += increment
                    # Store value of last or first diag if it is the case #W
                    (c == last_c) && extend && (lastbins[current_diag] += increment) #W
                    firstdiag && (firstbins[current_diag] += increment) #W
                    current_diag = 0
                end
            end
            # Store value of previous cell: stop first diagonal search when reaching false
            (previous_cell = x[c-d,c]) || (firstdiag=false) #W
        end
        # Isolated points #W
        (x[first_c-d,first_c] && !x[first_c+1-d,first_c+1]) && (firstsingle += increment) #W
        (!x[last_c-1-d,last_c-1] && x[last_c-d,last_c]) && (lastsingle += increment) #W
    end
    # Add isolated points in first bin
    allpoints = (theiler == 0) ? countnz(x) : countnz(triu(x, theiler)) + countnz(tril(x,-theiler))
    [allpoints - collect(2:nbins+1)'*bins; bins], [firstsingle; firstbins], [lastsingle; lastbins] #W
end
