var documenterSearchIndex = {"docs":
[{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Mishchenko, M.I., Travis, L.D., Lacis, A.A., 2002. Scattering, absorption, and emission of light by small particles. Cambridge university press.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GeneralizedMultiparticleMie","category":"page"},{"location":"#GeneralizedMultiparticleMie","page":"Home","title":"GeneralizedMultiparticleMie","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GeneralizedMultiparticleMie.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GeneralizedMultiparticleMie]","category":"page"},{"location":"#GeneralizedMultiparticleMie.associated_legendre-Tuple{DataType, Integer, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.associated_legendre","text":"Associated Legendre function without the Condon-Shotley phase.\n\nArblib's definition includes the Condon-Shotley phase, so we need to multiply the results by (-1)^m.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.associated_legendre_array-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.associated_legendre_array","text":"Calculate the associated Legendre function for all 0leq nleq n_max and -nleq mleq n, and return the results as a vector.\n\nArblib's Legendre function definition includes the Condon-Shotley phase, so we need to multiply the results by (-1)^m.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.expand_E_cluster-Tuple{DataType, GeneralizedMultiparticleMie.VSWFMode, Number}","page":"Home","title":"GeneralizedMultiparticleMie.expand_E_cluster","text":"Expand the electric field.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.factorial-Tuple{DataType, Any}","page":"Home","title":"GeneralizedMultiparticleMie.factorial","text":"Calculate factorials using the Gamma function.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.gaunt_a-Tuple{DataType, Integer, Integer, Integer, Integer, Integer}","page":"Home","title":"GeneralizedMultiparticleMie.gaunt_a","text":"The Gaunt a-coefficient defined by Gaunt (1929):\n\nbeginaligned\na(m n mu nu p)=(-1)^m+mu(2 p+1)leftfrac(n+m) (nu+mu) (p-m-mu) (n-m) (nu-mu) (p+m+mu) right^1  2 \n timesleft(beginarrayccc\nn  nu  p \n0  0  0\nendarrayright)left(beginarrayccc\nn  nu  p \nm  mu  -m-mu\nendarrayright)\nendaligned\n\nReferences:\n\nGaunt, J.A., 1929. IV. The triplets of helium. Philosophical Transactions of the Royal Society of London. Series A, Containing Papers of a Mathematical or Physical Character 228, 151–196.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.gaunt_b-Tuple{DataType, Integer, Integer, Integer, Integer, Integer}","page":"Home","title":"GeneralizedMultiparticleMie.gaunt_b","text":"beginaligned\nb(m n mu nu p)=(-1)^m+mu(2 p+3)leftfrac(n+m) (nu+mu) (p-m-mu) (n-m) (nu-mu) (p+m+mu) right^1  2 \n timesleft(beginarrayccc\nn  nu  p \n0  0  0\nendarrayright)left(beginarrayccc\nn  nu  p+1 \nm  mu  -m-mu\nendarrayright)\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.init_vswf_cache-Tuple{DataType, Integer}","page":"Home","title":"GeneralizedMultiparticleMie.init_vswf_cache","text":"Precompute required VSWF coefficients.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.pi_func-Tuple{DataType, Integer, Integer, Number}","page":"Home","title":"GeneralizedMultiparticleMie.pi_func","text":"Calculate pi_n^m(theta) defined as\n\npi_n^m(theta)=fracmsinthetaP_n^m(costheta)\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.ricatti_hn1-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.ricatti_hn1","text":"Riccati Hankel function of the first kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.ricatti_hn1_deriv-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.ricatti_hn1_deriv","text":"First-order derivative of Riccati Hankel function of the first kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.ricatti_hn2-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.ricatti_hn2","text":"Riccati Hankel function of the second kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.ricatti_hn2_deriv-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.ricatti_hn2_deriv","text":"First-order derivative of Riccati Hankel function of the second kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.ricatti_jn-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.ricatti_jn","text":"Riccati Bessel function of the first kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.ricatti_jn_deriv-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.ricatti_jn_deriv","text":"First-order derivative of Riccati Bessel function of the first kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.ricatti_yn-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.ricatti_yn","text":"Riccati Bessel function of the second kind.\n\nNote that in miepy, the author used -z⋅y(z) instead of z⋅y(z)\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.ricatti_yn_deriv-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.ricatti_yn_deriv","text":"First-order derivative of Riccati Bessel function of the second kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.spherical_hn1-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.spherical_hn1","text":"Spherical Hankel function of the first kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.spherical_hn1_deriv-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.spherical_hn1_deriv","text":"First-order derivative of spherical Hankel function of the first kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.spherical_hn2-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.spherical_hn2","text":"Spherical Hankel function of the second kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.spherical_hn2_deriv-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.spherical_hn2_deriv","text":"First-order derivative of spherical Hankel function of the second kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.spherical_jn-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.spherical_jn","text":"Spherical Bessel function of the first kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.spherical_jn_deriv-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.spherical_jn_deriv","text":"First-order derivative of spherical Bessel function of the first kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.spherical_yn-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.spherical_yn","text":"Spherical Bessel function of the second kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.spherical_yn_deriv-Tuple{DataType, Integer, Any}","page":"Home","title":"GeneralizedMultiparticleMie.spherical_yn_deriv","text":"First-order derivative of spherical Bessel function of the second kind.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.tau_func-Tuple{DataType, Integer, Integer, Number}","page":"Home","title":"GeneralizedMultiparticleMie.tau_func","text":"Calculate tau_n^m(theta) defined as\n\ntau_n^m(theta)=fracmathrmdmathrmdthetaP_n^m(costheta)\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.vswf_electric-Tuple{DataType, Integer, Integer, GeneralizedMultiparticleMie.VSWFMode, Number, Number, Number, Number}","page":"Home","title":"GeneralizedMultiparticleMie.vswf_electric","text":"Vector spherical wave function, electric (TM) modes.\n\nbeginaligned\nmathbfN_m n^(i)(rho theta phi)=lefthatmathbfe_r n(n+1) P_n^m(cos theta) fracz^(i)_n(rho)rhoright\nleft+lefthatmathbfe_theta tau_m n(cos theta)+hatmathbfe_phi mathrmi pi_m n(cos theta)right frac(rhocdot z_n^(i)(rho))^primerhoright exp (mathrmi m phi)\nendaligned\n\nwhere i=1234 relates to incident, interior, outgoing and ingoing VSWF.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.vswf_magnetic-Tuple{DataType, Integer, Integer, GeneralizedMultiparticleMie.VSWFMode, Number, Number, Number, Number}","page":"Home","title":"GeneralizedMultiparticleMie.vswf_magnetic","text":"Vector spherical wave function, magnetic (TE) modes.\n\nmathbfM_m n^(i)(rho theta phi)=lefthatmathbfe_theta mathrmi pi_m n(cos theta)-hatmathbfe_phi tau_m n(cos theta)right z_n^(i)(rho) exp (mathrmi m phi)\n\nwhere i=1234 relates to incident, interior, outgoing and ingoing VSWF.\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.wigner_3j","page":"Home","title":"GeneralizedMultiparticleMie.wigner_3j","text":"Wigner 3-j symbols.\n\n\n\n\n\n","category":"function"},{"location":"#GeneralizedMultiparticleMie.wigner_D-Tuple{DataType, Integer, Integer, Integer, Number, Number, Number}","page":"Home","title":"GeneralizedMultiparticleMie.wigner_D","text":"Wigner D-function D^j_mn(theta) defined as:\n\nD_m m^prime^n(alpha beta gamma)=mathrme^-mathrmi m alpha d_m m^prime^n(beta) mathrme^-mathrmi m^prime gamma\n\nwhere\n\n0 leq alpha2 pi quad 0 leq beta leq pi quad 0 leq gamma2 pi\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.wigner_d-Tuple{DataType, Integer, Integer, Integer, Number}","page":"Home","title":"GeneralizedMultiparticleMie.wigner_d","text":"Wrapper of WignerD.jl's WignerD.wignerdjmn() function, which is much faster than wigner_d_naive().\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.wigner_d_naive-Tuple{DataType, Integer, Integer, Integer, Number}","page":"Home","title":"GeneralizedMultiparticleMie.wigner_d_naive","text":"Calculate Wigner (small) d-function d_mn^s(theta) for a single (m n s) combination, using Eq. (B.1) of Mishchenko et al. (2002).\n\nbeginaligned\nd_m n^s(vartheta)= sqrt(s+m) (s-m) (s+n) (s-n)  \n times sum_k=max(0m-n)^min(s + m s - n)(-1)^k fracleft(cos frac12 varthetaright)^2 s-2 k+m-nleft(sin frac12 varthetaright)^2 k-m+nk (s+m-k) (s-n-k) (n-m+k) \nendaligned\n\n\n\n\n\n","category":"method"},{"location":"#GeneralizedMultiparticleMie.wigner_d_recursion-Tuple{DataType, Integer, Integer, Integer, Number}","page":"Home","title":"GeneralizedMultiparticleMie.wigner_d_recursion","text":"Calculate Wigner (small) d-function d_mn^s(theta) for sins_min=max(m n)s_max (and also its derivative) via upward recursion, using Eq. (B.22) of Mishchenko et al. (2002).\n\nbeginaligned\nd_m n^s+1(vartheta)= frac1s sqrt(s+1)^2-m^2 sqrt(s+1)^2-n^2left(2 s+1)s(s+1) x-m n d_m n^s(vartheta)right\nleft-(s+1) sqrts^2-m^2 sqrts^2-n^2 d_m n^s-1(vartheta)right quad s geq s_min \nendaligned\n\nThe initial terms are given by Eq. (B.23) and Eq. (B.24).\n\nbeginarrayl\nd_m n^s_min -1(vartheta)=0 \nd_m n^s_min (vartheta)=xi_m n 2^-s_min leftfracleft(2 s_min right) (m-n) (m+n) right^1  2(1-x)^m-n  2(1+x)^m+n  2\nendarray\n\nwhere\n\nxi_m n=leftbeginarrayll\n1  text  for  n geq m \n(-1)^m-n  text  for  nm\nendarrayright\n\n\n\n\n\n","category":"method"}]
}
