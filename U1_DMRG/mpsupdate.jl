#For updating mps by solving locally smallest eigenvalue
function hmap(lt, mpo2, rt, tv)
    @tensor hv[b2l b1r; o1 p11 p21] := lt[b1l; b2l o1] * tv[b1l, p11, p21, b1r]
    @tensor hv[b2l p12 p22; b1r o2] := hv[b2l b1r; o1 p11 p21] * mpo2[o1 p11 p21; p12 p22 o2]
    @tensor hv[b2l, p12, p22, b2r] := hv[b2l p12 p22; b1r o2] * rt[b1r; b2r o2]
    return hv
end
(h::Vector)(v) = hmap(h[1], h[2], h[3], v)