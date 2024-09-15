#For contracting the Tensor
function cmpsleft(lt, ten1, oten, ten2)
    #The index after ' should be changed!!!
    # b0← b0←o←b1←       b1←      p0←           b1←                     b1←
    # o1←    p1←    -->  o0← p0← o0←□←o1←  -->  o1← p0←   p0←      -->  o1←
    # b2→                b2→        p2←         b0→     →b0→o→b2→       b2→
    @tensor ten0[b1; p1 b2 o1] := lt[b0; b2 o1] * ten1[b0 p1; b1]
    @tensor ten0[b1; p2 b2 o1] := ten0[b1; p0 b2 o0] * oten[o0 p0; p2 o1]
    @tensor ten0[b1; b2 o1] := ten0[b1; p0 b0 o1] * ten2[b0 p0; b2]
    return ten0
end

# function cmpoleft(lt, ten1, ten2)
#     #The index after ' should be changed!!!
#     # b0← b0←o←b1←       b1←      p0←           b1←                     b1←
#     # o1←    p1←    -->  o0← p0← o0←□←o1←  -->  o1← p0←   p0←      -->  o1←
#     # b2→                b2→        p2←         b0→     →b0→o→b2→       b2→
#     @tensor ten0[b1; p1 b2 o1] := lt[b0; b2 o1] * ten1[b0 p1; b1]
#     @tensor ten0[b1; p2 b2 o1] := ten0[b1; p0 b2 o0] * oten[o0 p0; p2 o1]
#     @tensor ten0[b1; b2 o1] := ten0[b1; p0 b0 o1] * ten2[b0 p0; b2]
#     return ten0
# end

function cmpsright(rt, ten1, oten, ten2)
    #The index after ' should be changed!!!
    #  ←b1                    ←b1         p0←        ←b1       ←b1←o←b0  ←b0
    #  ←o1  <--    p0←    p0← ←o1  <--  ←o1←□←o0 p0← ←o0  <--      p1←   ←o1
    #  →b2       →b2→o→b0     →b0           p2←      →b2                 →b2
    @tensor ten0[b1; p1 b2 o1] := rt[b0; b2 o1] * ten1[b1 p1; b0]
    @tensor ten0[b1; p2 b2 o1] := ten0[b1; p0 b2 o0] * oten[o1 p0; p2 o0]
    @tensor ten0[b1; b2 o1] := ten0[b1; p0 b0 o1] * ten2[b2 p0; b0]
    return ten0
end

#For compute lt and rt
function ltrt!(lt, rt, mps, mpo, l)
    onelt = TensorMap(ones, space(mps[1],1)', space(mps[1],1)' ⊗ space(mpo[1],1))
    lt[1] = onelt
    for il = 2:l+1
        lt[il] = cmpsleft(lt[il-1], mps[il-1], mpo[il-1], permute(mps[il-1]', (3,), (1, 2)))
    end
    onert = TensorMap(ones, space(mps[l],3)', space(mps[l],3)' ⊗ space(mpo[l],4))
    rt[l+1] = onert
    for il = l:-1:1
        rt[il] = cmpsright(rt[il+1], mps[il], mpo[il], permute(mps[il]', (3,), (1, 2)))
    end
end
