import LambertW
import DSP

"""
k is the order of the pole
T is the refractory period
γ is the proliferation rate
"""
pole(k::Int,T::Float64,γ::Float64) = -γ  + LambertW.lambertw(2T*γ* exp(T*γ)+0im,k)/T  # Expression II.11 from Supplementary material


ϕ(s,T,γ) =  s + γ * (1 - 2*exp(-T*s) ) #Expression II.9 from Supplementary material

ϕp(s,T,γ) = 1 + 2*γ*T* exp(-T*s) 



"conv is a convolution function that uses DSP.conv from Julia"

function conv(t,f,g)
    z=zero(t)
    dt=t[2]-t[1]
    c=DSP.conv(vcat(f,z),vcat(g,z))
    c=c[1:length(t)]
    return dt*c
end
               

"""
 nt_var computes population average and variance as function of time using the
 number of poles given by npoles. 
`t` is an array
`T`, `γ`, and `t0` are scalar parameters
"""
function nt_var(t::Array{Float64},npoles::Int,T::Float64,γ::Float64,t0::Float64)
    t=t.+t0
    tmT = t.-T
    tmT=tmT[tmT.>=0]
    f=zero(t)
    ψ=zero(t)

    #
    # compute n(t)
    #
    sn0=pole(0,T,γ)
    f[t.>=T] = 1/ϕ(0,T,γ) .+  exp.(sn0 * tmT) / (sn0 * ϕp(sn0,T,γ)) 
    ψ = real( exp.(sn0 * t) / ϕp(sn0,T,γ) )
    for i in 2:npoles
        sni=pole(i-1,T,γ)
        f[t.>=T] += 2* real( exp.(sni*tmT) / (sni * ϕp(sni,T,γ) ) )
        ψ += 2* real( exp.(sni*t) / ϕp(sni,T,γ) )
    end
    nt = 1 .+ γ*f

    #
    # compute variance
    #
    ntsq=nt.*nt
    h=conv(t,ψ,ntsq)
    h=DSP.shiftsignal(h,length(t[t.<T]))
    var=nt - ntsq +2*γ*h
    return nt,var

end



