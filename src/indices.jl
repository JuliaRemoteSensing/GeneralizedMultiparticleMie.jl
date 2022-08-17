@inline rmax_to_lmax(rmax::Integer) = Int(floor(√(rmax + 1))) - 1

@inline lmax_to_rmax(lmax::Integer) = lmax * (lmax + 2)
