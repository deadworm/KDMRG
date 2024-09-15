#For updating mps by solving locally smallest eigenvalue
function hmap(lt, mpo2, rt, tv)
    hv = 0
    for im in eachindex(mpo2)
        @tensor hvm[b2l b1r; o1 p11 p21] := lt[im][b1l; b2l o1] * tv[b1l, p11, p21, b1r]
        @tensor hvm[b2l p12 p22; b1r o2] := hvm[b2l b1r; o1 p11 p21] * mpo2[im][o1 p11 p21; p12 p22 o2]
        @tensor hvm[b2l, p12, p22, b2r] := hvm[b2l p12 p22; b1r o2] * rt[im][b1r; b2r o2]
        if im == 1
            hv = hvm
        else
            hv = hv + hvm
        end
    end
    return hv
end
(h::Vector)(v) = hmap(h[1], h[2], h[3], v)

#For test hermitian
# for il = 1:l-1
#     for im = 1:length(mpolist)
#         @tensor m1[b1l b1r p11 p21; p12 p22 b2l b2r] := ltlist[im, il][b1l; b2l o1] * mpo2list[im, il][o1 p11 p21; p12 p22 o2] * rtlist[im, il+2][b1r; b2r o2]
#         mm1 = convert(Array, m1)
#         mm1p = convert(Array, permute(m1',(3,4,1,2),(7,8,5,6)))
#         if mm1p ≈ mm1
#         else
#             @show il, im
#             @show maximum(real(mm1p - mm1))
#         end
#     end
# end