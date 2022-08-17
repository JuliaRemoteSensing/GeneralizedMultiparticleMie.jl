const Vec3{T} = StaticArrays.SVector{3,T}

Base.complex(::Type{Arblib.Arb}) = Arblib.Acb
Base.complex(::Type{Arblib.ArbRef}) = Arblib.Acb
Base.complex(::Type{Arblib.Acb}) = Arblib.Acb
Base.complex(::Type{Arblib.AcbRef}) = Arblib.Acb
