#For contracting the fermionic Tensor
function con_left(lt, ten1, ten2)
    nlt=zeros(Float64, space(ten1, 1) ← space(ten2, 3))
    @planar nlt[r1; r2] := lt[l1; l2] * ten1[r1; l1 u] * ten2[l2, u; r2]
    return nlt
end

function con_right(rt, ten1, ten2)
    nrt=zeros(Float64, space(ten2, 1) ← space(ten1, 2))
    @planar nrt[l2; l1] := ten1[r1; l1 u] * ten2[l2, u; r2] * rt[r2; r1]
    return nrt
end

function con_lrt(lt, rt, ten)
    nten=zeros(Float64, space(lt, 1) ⊗ space(ten, 2) ← space(rt, 2))
    @planar nten[l1 u; r1] := lt[l1; l2] * ten[l2 u; r2] * rt[r2; r1]
    return nten
end

#For compute lt and rt
function ltrt!(lt, rt, m, mf, l)
    lt[1] = isometry(Float64, zb ← zb)
    for il = 2:l+1
        lt[il] = con_left(lt[il-1], m[il-1]', mf[il-1])
    end
    rt[l+1] = isometry(Float64, nb ← nb)
    for il = l:-1:1
        rt[il] = con_right(rt[il+1], m[il]', mf[il])
    end
end
