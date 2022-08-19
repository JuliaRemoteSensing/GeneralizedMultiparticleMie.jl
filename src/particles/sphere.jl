struct Sphere{T} <: AbstractParticle{T}
    pos::NTuple{3,T}
    r::T
    m::NTuple{2,T}
end

position(sphere::Sphere{T}) = sphere.pos
orientation(::Sphere{T}) = zero(Quaternions.Quaternion{T})
circumscribed_sphere_radius(sphere::Sphere{T}) = sphere.r
volume(sphere::Sphere{T}) = sphere.r^3 * 4 / 3 * π
tmatrix(sphere::Sphere{T}) = error("Method not implemented")
Base.:∈(pos, sphere::Sphere{T}) = LinearAlgebra.norm(pos .- sphere.pos) <= sphere.r
