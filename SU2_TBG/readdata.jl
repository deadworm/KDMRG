#read the data from txt
nx = 3
ny = nx
Nk = nx * ny
gx = 11
gy = gx
NG = gx * gy

an = 1.08
theta = an / 180 * pi
a0 = 0.246
d = 40
epsl = 7
Lm = a0 / 2 / sin(theta / 2)
G1 = [-1 / 2, -sqrt(3) / 2]
G2 = [1, 0]

data = []
error = 10^-5;
open("SU2_TBG/lamdaTBG_NC" * string(nx) * "_" * string(an) * ".txt", "r") do file
    for line in eachline(file)
        row = parse.(Float64, split(line, "\t"))
        row = map(x -> abs(x) < error ? 0 : x, row)
        push!(data, row)
    end
end

Lamda = Array{Any}(undef, nx * gx, ny * gy)
Lamcheck = zeros(nx * gx, ny * gy)
for ix = 1:nx*gx
    for iy = 1:ny*gy
        Lamda[ix, iy] = zeros(Complex{Float64}, Nk * 2, Nk * 2)
    end
end
lamda = Array{Complex}(undef, 2, 2)
for i = 1:4:length(data)
    k1 = Int(data[i][3])
    k2 = Int(data[i][4])
    ky1 = mod(k1 - 1, ny) + 1
    kx1 = ceil(k1 / ny)
    ky2 = mod(k2 - 1, ny) + 1
    kx2 = ceil(k2 / ny)
    qy = ky2 - ky1
    qx = kx2 - kx1

    Gn = data[i][5]
    ig1 = mod(Gn - 1, gy) + 1
    ig2 = ceil(Gn / gy)
    igx = ig1 - 1 - (gx - 1) / 2
    igy = ig2 - 1 - (gy - 1) / 2
    qGy = Int(qy + igy * ny)
    qGx = Int(qx + igx * nx)

    lamda[1, 1] = data[i][6] + im * data[i][7]
    lamda[1, 2] = data[i+1][6] + im * data[i+1][7]
    lamda[2, 1] = data[i+2][6] + im * data[i+2][7]
    lamda[2, 2] = data[i+3][6] + im * data[i+3][7]

    Lamda[mod(qGx, nx * gx)+1, mod(qGy, ny * gy)+1][(k1-1)*2+1:(k1-1)*2+2, (k2-1)*2+1:(k2-1)*2+2] = lamda
    Lamda[mod(-qGx, nx * gx)+1, mod(-qGy, ny * gy)+1][(k2-1)*2+1:(k2-1)*2+2, (k1-1)*2+1:(k1-1)*2+2] = lamda'

    Lamcheck[mod(qGx, nx * gx)+1, mod(qGy, ny * gy)+1] += 1
    Lamcheck[mod(-qGx, nx * gx)+1, mod(-qGy, ny * gy)+1] += 1

    #assume hermitian
    # Lamda[mod(qGx, nx * gx)+1, mod(qGy, ny * gy)+1][(k2-1)*2+1:(k2-1)*2+2, (k1-1)*2+1:(k1-1)*2+2] = lamda'
    # Lamda[mod(-qGx, nx * gx)+1, mod(-qGy, ny * gy)+1][(k1-1)*2+1:(k1-1)*2+2, (k2-1)*2+1:(k2-1)*2+2] = lamda
end

hlist = []
aqlist = []
miu = -1
for ix = 1:nx*gx
    for iy = 1:ny*gy
        if Lamda[ix, iy] == zeros(Complex{Float64}, Nk * 2, Nk * 2)
        else
            if iy > ny * gy / 2
                y = iy - ny * gy
            else
                y = iy
            end
            if ix > nx * gx / 2
                x = ix - nx * gx
            else
                x = ix
            end
            qG = norm((x - 1) * G1 + (y - 1) * G2)
            Aq = sqrt((84.93 / epsl / Nk * sin(theta / 2) * sqrt(3) / 4 / pi * sqrt(Nk) / qG * (1 - exp(-4 * pi / sqrt(3) / sqrt(Nk) * d / Lm * qG))) / 2)
            if Lamcheck[ix, iy] == Nk
                push!(hlist, [Aq * Lamda[ix, iy], 0, mod(-(ix - 1), nx), mod(-(iy - 1), ny)])
                push!(aqlist, Aq)
            else
                if ix == 1 && iy == 1
                else
                    push!(hlist, [Aq * Lamda[ix, iy], miu * Aq * tr(Lamda[ix, iy]), mod(-(ix - 1), nx), mod(-(iy - 1), ny)])
                    push!(aqlist, Aq)
                end
            end
        end
    end
end

# aq0=aqlist[2]
# for ih in eachindex(hlist)
#     if hlist[ih][2] != 0
#         println(ih, " ", hlist[ih][3], " ", hlist[ih][4], " ", hlist[ih][2])
#     end

#     # if aqlist[ih] == aq0
#     #     println(ih, " ", hlist[ih][3], " ", hlist[ih][4], " ", hlist[ih][2])
#     # end
# end