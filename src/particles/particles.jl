abstract type AbstractParticle{T} end

# Methods need to be implemented for each concrete type.
position(::AbstractParticle) = error("Method not implemented")
orientation(::AbstractParticle) = error("Method not implemented")
circumscribed_sphere_radius(::AbstractParticle) = error("Method not implemented")
volume(::AbstractParticle) = error("Method not implemented")
tmatrix(::AbstractParticle) = error("Method not implemented")
Base.:∈(pos, ::AbstractParticle) = error("Method not implemented")

# Methods that can be automatically implemented.
volume_equivalent_radius(p::AbstractParticle) = ∛(3 // 4 * volume(p) / π)

include("sphere.jl")
