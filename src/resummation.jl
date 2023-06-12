function resummation(; Nmax, raw_sum_order, nwynn, Ofilepath, Ofilename)

    O = load(Ofilepath * Ofilename * ".jld2")["O"]
    O_size = size(O)
    Nmax_raw = O_size[2]
    NT = O_size[3]
    Nh = O_size[4]
    @assert(Nmax <= Nmax_raw, "maximum order must be smaller than that of the raw data: " * string(Nmax_raw))
    @assert(raw_sum_order <= Nmax, "Euler sum order must be smaller than Nmax: " * string(Nmax))
    # ---------------------------------------------------------
    # O = zeros(Float64, 6, Nmax + 1, NT, Nh)
    # In this part, we add partial sums to obtain the NLCE sums
    # upto different orders. We can do it without any tricks
    # (raw sums), or using numerical resummation tricks (here
    # we have implemented the Wynn algorithm, see Sec. 4 of
    # https://arxiv.org/pdf/1207.3366.pdf for more details).
    # The Euler resummation algorithm is also useful and worth
    # implementing.

    # Doing the raw sums for NLCE
    raw_O = zeros(Float64, 6, Nmax, NT, Nh)
    for N = 1:(Nmax)
        if N != 1
            raw_O[:, N, :, :] = raw_O[:, N-1, :, :] + O[:, N, :, :]
        else
            raw_O[:, N, :, :] = O[:, N, :, :]
        end
    end
    #println("raw sums: ")
    #print(raw_O[2,:,1,1])


    # try to write Wynn algorithm by shifting Array index to start from -1 by using package OffsetArrays
    # nwynn = 5
    epsln_O = Origin(1, 1, 1, 1, -1)(zeros(Float64, 6, Nmax, NT, Nh, 18))
    epsln_O[:, 1:Nmax, :, :, -1] = zeros(Float64, 6, Nmax, NT, Nh)
    epsln_O[:, 1:Nmax, :, :, 0] = raw_O[:, 1:Nmax, :, :]
    for k = 1:(2*nwynn)
        for n = 1:(Nmax-k)
            delta = epsln_O[:, n+1, :, :, k-1] - epsln_O[:, n, :, :, k-1]
            epsln_O[:, n, :, :, k] = (epsln_O[:, n+1, :, :, k-2] + 1 ./ delta)
        end
    end


    # Doing Euler summation
    # raw_sum_order = 10
    m_max = Nmax - raw_sum_order + 1
    Euler_O = zeros(Float64, 6, NT, Nh, m_max)
    # m=1 is the raw sum to the order of raw_sum_order
    Euler_O[:, :, :, 1] = raw_O[:, raw_sum_order, :, :]
    for m = 2:m_max
        Euler_O[:, :, :, m] = Euler_O[:, :, :, m-1]
        for n = 1:(m-1)
            Euler_O[:, :, :, m] +=
                (1 / 2)^(m - 1) * O[:, raw_sum_order+n, :, :] * binomial(m - 2, n - 1)
        end
    end

    # return [epsln_O[:, Nmax-nwynn, :, :, nwynn+1], Euler_O[:, :, :, m_max], raw_O]

    return [OffsetArrays.no_offset_view(epsln_O)[:, Nmax-2*nwynn, :, :, 2*nwynn+2], Euler_O[:, :, :, m_max], raw_O]
end
