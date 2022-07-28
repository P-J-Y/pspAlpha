
using MAT
using Dates
using Plots
using LinearAlgebra
using Statistics
using Optim
using Interpolations
using CurveFit
using JLD2
using XLSX

# 原始磁场数据，每个文件对应的时间
# magFileTime = Dict(
# "2020b"=>[DateTime(2020,7,1),DateTime(2021,1,1)],
# "2021a"=>[DateTime(2021,1,1),DateTime(2021,7,1)],
# "2021b"=>[DateTime(2021,7,1),DateTime(2022,1,1)],
# )

goodevent = [267 284 318 531 562 1224 1231 1247 1253 1290]
notSB = [188,190,198,203,204,205,211,215,216,217,218,
221,222,224,237,245,250,252,258,262,
263,265,272,288,290,296,302,313,318,
320,321,322,325,326,329,330,331,332,
336,339,344,345,346,347,348,457,513,
514,518,519,524,531,533,534,535,536,
537,538,542,543,548,550,555,556,569,
570,582,583,584,588,592,594,603,612,
616,619,620,621,622,623,624,626,630,
631,632,633,785,786,791,794,874,878,
887,889,891,911,913,926,927,928,930,
931,932,933,934,939,940,941,942,943,
948,950,951,952,953,956,958,959,1169,
1170,1171,1172,1173,1174,1181,1182,
1183,1189,1190,1191,1194,1195,1196,
1198,1199,1200,1201,1206,1213,1215,
1222,1224,1225,1226,1227,1232,1241,
1245,1248,1250,1257,1260,1263,1266,
1278,1280,1281,1282,1284,1289,1294,
1295,1303,1309,1316,1317,1321,
] # sb: 1、速度要有跃变 2、br的跳动比较明显
maybenotSB = [191,192,193,196,199,200,
202,206,208,209,210,219,223,232,239,244,
246,259,273,284,299,300,306,307,310,323,
328,334,337,515,532,541,546,547,553,559,
566,568,572,573,574,575,577,578,580,593,
596,598,602,611,617,618,627,628,629,634,
742,790,873,878,879,880,882,883,884,885,
886,888,892,905,906,908,909,914,922,924,
929,944,945,949,1175,1178,1184,1185,1186,
1188,1192,1204,1207,1210,1211,1212,1214,
1216,1217,1218,1219,1220,1228,1229,1233,
1243,1244,1246,1249,1251,1253,1254,1256,
1261,1265,1267,1269,1270,1272,1285,1287,
1291,1292,1293,1298,1300,1308,1319,1320,
]
#790附近有几个事件，磁场有明显的转向，但是速度非常平稳，挺奇怪的
# 884这样的感觉是多个SB

function loadData()
    pVars = matread("data\\psp_spi_sf00.mat")
    αVars = matread("data\\psp_spi_sf0a.mat")
    modifiedVars = matread("data\\psp_modified_alpha.mat")
    modifiedVars_va = matread("data\\psp_modified_va.mat")
    pVars,αVars,modifiedVars,modifiedVars_va
end

function loadList()
    sbList = matread("list//switchback_time_output.mat")
end


"""
epoch(matlab's datenum) to julia datetime
"""
function epoch2datetime(epoch::Real)
    unix2datetime((epoch-719529)*86400)
end

function datetime2epoch(adatetime::DateTime)
    datetime2unix(adatetime)/86400 + 719529
end

function get_circle_model(xs_vec::Matrix)
    circle_model(p) = sum(abs2,sum(abs2,xs_vec.-p[1:2],dims=1).-p[3]^2)
end

function get_circle(xs_vec::Matrix,p0::Vector)
    circle_model = get_circle_model(xs_vec)
    res = optimize(circle_model,p0)
    Optim.minimizer(res)
end

function get_circle_model_rfixed(xs_vec::Matrix,r)
    circle_model(p) = sum(abs2,sum(abs2,xs_vec.-p[1:2],dims=1).-r^2)
end

function get_circle_rfixed(xs_vec::Matrix,p0::Vector,r)
    circle_model = get_circle_model_rfixed(xs_vec,r)
    res = optimize(circle_model,p0)
    Optim.minimizer(res)
end

function circleShape(p::Vector)
    θ = LinRange(0,2π,500)
    p[1] .+ p[3]*sin.(θ), p[2] .+ p[3]*cos.(θ)
end

function spanAngle(v1::Vector,v2::Vector)
    acos(dot(v1,v2)/(norm(v1)*norm(v2)))
end

"""
计算p和α之间的速度差标量（α速度高为正）
"""
function Vαp(vα,vp)
    sign(norm(vα)-norm(vp))*norm(vα-vp)
end

function sbEvent(epoch1,epoch2,pVars,αVars,modifiedVars,magVars;
    figName="test",plotTimeSeries=false,deltaMins=10,figDir="issb")
    # figDir2sbtype = Dict("issb"=>1, "maybesb"=>2, "notsb"=>3)
    output = Dict(
    "αdataFlag"=>1,
    # "ccpb"=>NaN,
    # "ccαb"=>NaN,
    # "Vwp"=>NaN,
    # "rp"=>NaN,
    # "Vwα"=>NaN,
    # "rα"=>NaN,
    # "Cp"=>(NaN,NaN),
    # "Cα"=>(NaN,NaN),
    # "sbtype"=>figDir2sbtype[figDir], #1:issb,2:maybesb,3:notsb
    )
    deltaEpoch = deltaMins/(24*60)
    deltaEpochSmall = 1/(24*60)
    #############
    pPoints = vec((pVars["p_epoch"].>=epoch1) .& (pVars["p_epoch"].<=epoch2))
    pEpoch = pVars["p_epoch"][pPoints]
    pTime = epoch2datetime.(pEpoch)
    pVel = pVars["p_vel_rtn_sun"][pPoints,:]
    pTemp = pVars["p_temp"][pPoints]
    αPoints =  vec((αVars["alpha_epoch"].>=epoch1) .&
                (αVars["alpha_epoch"].<=epoch2))
    αEpoch = αVars["alpha_epoch"][αPoints]
    αTime = epoch2datetime.(αEpoch)
    αVel = αVars["alpha_vel_rtn_sun"][αPoints,:]
    pVel2αEpoch = modifiedVars["p_vel_rtn_sun_alphaEpoch"][αPoints,:]
    αTemp = αVars["alpha_temp"][αPoints]
    va = modifiedVars["va_alphaEpoch"][αPoints]
    va_rtn = modifiedVars["va_rtn_alphaEpoch"][αPoints,:]
    magPoints = vec((magVars["mag_epoch"].>=epoch1) .&
                (magVars["mag_epoch"].<=epoch2))
    magEpoch = magVars["mag_epoch"][magPoints]
    magTime = epoch2datetime.(magEpoch)
    mag_rtn = magVars["mag_rtn"][magPoints,:]
    # SB之前的一段
    pPoints1 = vec((pVars["p_epoch"].>=(epoch1-deltaEpoch)) .& (pVars["p_epoch"].<=(epoch1-deltaEpochSmall)))
    pEpoch1 = pVars["p_epoch"][pPoints1]
    pTime1 = epoch2datetime.(pEpoch1)
    pVel1 = pVars["p_vel_rtn_sun"][pPoints1,:]
    pTemp1 = pVars["p_temp"][pPoints1]
    αPoints1 =  vec((αVars["alpha_epoch"].>=(epoch1-deltaEpoch)) .&
                (αVars["alpha_epoch"].<=(epoch1-deltaEpochSmall)))
    αEpoch1 = αVars["alpha_epoch"][αPoints1]
    αTime1 = epoch2datetime.(αEpoch1)
    αVel1 = αVars["alpha_vel_rtn_sun"][αPoints1,:]
    pVel2αEpoch1 = modifiedVars["p_vel_rtn_sun_alphaEpoch"][αPoints1,:]
    αTemp1 = αVars["alpha_temp"][αPoints1]
    va1 = modifiedVars["va_alphaEpoch"][αPoints1]
    va_rtn1 = modifiedVars["va_rtn_alphaEpoch"][αPoints1,:]
    magPoints1 = vec((magVars["mag_epoch"].>=(epoch1-deltaEpoch)) .&
                (magVars["mag_epoch"].<=(epoch1-deltaEpochSmall)))
    magEpoch1 = magVars["mag_epoch"][magPoints1]
    magTime1 = epoch2datetime.(magEpoch1)
    mag_rtn1 = magVars["mag_rtn"][magPoints1,:]
    # SB之后的一段
    pPoints2 = vec((pVars["p_epoch"].>=(epoch2+deltaEpochSmall)) .& (pVars["p_epoch"].<=(epoch2+deltaEpoch)))
    pEpoch2 = pVars["p_epoch"][pPoints2]
    pTime2 = epoch2datetime.(pEpoch2)
    pVel2 = pVars["p_vel_rtn_sun"][pPoints2,:]
    pTemp2 = pVars["p_temp"][pPoints2]
    αPoints2 =  vec((αVars["alpha_epoch"].>=(epoch2+deltaEpochSmall)) .&
                (αVars["alpha_epoch"].<=(epoch2+deltaEpoch)))
    αEpoch2 = αVars["alpha_epoch"][αPoints2]
    αTime2 = epoch2datetime.(αEpoch2)
    αVel2 = αVars["alpha_vel_rtn_sun"][αPoints2,:]
    pVel2αEpoch2 = modifiedVars["p_vel_rtn_sun_alphaEpoch"][αPoints2,:]
    αTemp2 = αVars["alpha_temp"][αPoints2]
    va2 = modifiedVars["va_alphaEpoch"][αPoints2]
    va_rtn2 = modifiedVars["va_rtn_alphaEpoch"][αPoints2,:]
    magPoints2 = vec((magVars["mag_epoch"].>=(epoch2+deltaEpochSmall)) .&
                (magVars["mag_epoch"].<=(epoch2+deltaEpoch)))
    magEpoch2 = magVars["mag_epoch"][magPoints2]
    magTime2 = epoch2datetime.(magEpoch2)
    mag_rtn2 = magVars["mag_rtn"][magPoints2,:]
    # 判断这段时间是否有足够的阿尔法粒子数据
    if (length(αVel)-sum(isnan.(αVel))) < 10
        output["αdataFlag"] = 0
        return output
    end
    if (length(αVel1)-sum(isnan.(αVel1))) < 10
        output["αdataFlag"] = 0
        return output
    end
    if (length(αVel2)-sum(isnan.(αVel2))) < 10
        output["αdataFlag"] = 0
        return output
    end
    nonNaNpoints = vec(any(isnan.(mag_rtn),dims=2) .!= 1)
    magEpoch = magEpoch[nonNaNpoints]
    magTime = magTime[nonNaNpoints]
    mag_rtn = mag_rtn[nonNaNpoints,:]
    nonNaNpoints = vec(any(isnan.(pVel),dims=2) .!= 1)
    pEpoch = pEpoch[nonNaNpoints]
    pTime = pTime[nonNaNpoints]
    pVel = pVel[nonNaNpoints,:]
    pTemp = pTemp[nonNaNpoints]
    nonNaNpoints = vec(any(isnan.(αVel),dims=2) .!= 1)
    αEpoch = αEpoch[nonNaNpoints]
    αTime = αTime[nonNaNpoints]
    αVel = αVel[nonNaNpoints,:]
    αTemp = αTemp[nonNaNpoints]
    pVel2αEpoch = pVel2αEpoch[nonNaNpoints,:]
    va = va[nonNaNpoints]
    va_rtn = va_rtn[nonNaNpoints,:]

    nonNaNpoints1 = vec(any(isnan.(mag_rtn1),dims=2) .!= 1)
    magEpoch1 = magEpoch1[nonNaNpoints1]
    magTime1 = magTime1[nonNaNpoints1]
    mag_rtn1 = mag_rtn1[nonNaNpoints1,:]
    nonNaNpoints1 = vec(any(isnan.(pVel1),dims=2) .!= 1)
    pEpoch1 = pEpoch1[nonNaNpoints1]
    pTime1 = pTime1[nonNaNpoints1]
    pVel1 = pVel1[nonNaNpoints1,:]
    pTemp1 = pTemp1[nonNaNpoints1]
    nonNaNpoints1 = vec(any(isnan.(αVel1),dims=2) .!= 1)
    αEpoch1 = αEpoch1[nonNaNpoints1]
    αTime1 = αTime1[nonNaNpoints1]
    αVel1 = αVel1[nonNaNpoints1,:]
    αTemp1 = αTemp1[nonNaNpoints1]
    pVel2αEpoch1 = pVel2αEpoch1[nonNaNpoints1,:]
    va1 = va1[nonNaNpoints1]
    va_rtn1 = va_rtn1[nonNaNpoints1,:]

    nonNaNpoints2 = vec(any(isnan.(mag_rtn2),dims=2) .!= 1)
    magEpoch2 = magEpoch2[nonNaNpoints2]
    magTime2 = magTime2[nonNaNpoints2]
    mag_rtn2 = mag_rtn2[nonNaNpoints2,:]
    nonNaNpoints2 = vec(any(isnan.(pVel2),dims=2) .!= 1)
    pEpoch2 = pEpoch2[nonNaNpoints2]
    pTime2 = pTime2[nonNaNpoints2]
    pVel2 = pVel2[nonNaNpoints2,:]
    pTemp2 = pTemp2[nonNaNpoints2]
    nonNaNpoints2 = vec(any(isnan.(αVel2),dims=2) .!= 1)
    αEpoch2 = αEpoch2[nonNaNpoints2]
    αTime2 = αTime2[nonNaNpoints2]
    αVel2 = αVel2[nonNaNpoints2,:]
    αTemp2 = αTemp2[nonNaNpoints2]
    pVel2αEpoch2 = pVel2αEpoch2[nonNaNpoints2,:]
    va2 = va2[nonNaNpoints2]
    va_rtn2 = va_rtn2[nonNaNpoints2,:]

    ### 计算这段时间的磁场最小扰动方向（MVA）
    mag_rtn_scaled = mag_rtn .- mean(mag_rtn,dims=1)
    # Σb = (mag_rtn_scaled'*mag_rtn_scaled)/length(magEpoch)
    # # @show Σb
    # F = svd(Σb)
    # v1 = F.U[:,1]
    # v2 = F.U[:,2]
    # v3 = F.U[:,3]

    mag_rtn_scaled1 = mag_rtn1 .- mean(mag_rtn1,dims=1)
    # Σb1 = (mag_rtn_scaled1'*mag_rtn_scaled1)/length(magEpoch1)
    # # @show Σb
    # F1 = svd(Σb1)
    # v1_1 = F1.U[:,1]
    # v2_1 = F1.U[:,2]
    # v3_1 = F1.U[:,3]

    mag_rtn_scaled2 = mag_rtn2 .- mean(mag_rtn2,dims=1)
    # Σb2 = (mag_rtn_scaled2'*mag_rtn_scaled2)/length(magEpoch2)
    # # @show Σb
    # F2 = svd(Σb2)
    # v1_2 = F2.U[:,1]
    # v2_2 = F2.U[:,2]
    # v3_2 = F2.U[:,3]

    mag_rtn_scaled_tot = [mag_rtn_scaled;mag_rtn_scaled1;mag_rtn_scaled2]
    Σb = (mag_rtn_scaled_tot'*mag_rtn_scaled_tot) /
    (length(magEpoch)+length(magEpoch1)+length(magEpoch2))
    # @show Σb
    F = svd(Σb)
    v1 = F.U[:,1]
    v2 = F.U[:,2]
    v3 = F.U[:,3]
    v1_1 = v1
    v1_2 = v1
    v2_1 = v2
    v2_2 = v2
    v3_1 = v3
    v3_2 = v3
    # @show norm(v1)


    # P = I - v3*v3' # P*x给出x在P的平面上的投影
    mag2P = [(v1'*mag_rtn')' (v2'*mag_rtn')']
    p0_mag = [0.,0.,10.]
    p_mag = get_circle(Matrix(mag2P'),p0_mag)
    pVel2P = [(v1'*pVel')' (v2'*pVel')']
    pVel2Pα = [(v1'*pVel2αEpoch')' (v2'*pVel2αEpoch')']

    va2P = [(v1'*va_rtn')' (v2'*va_rtn')']
    p0_va = [0.,0.,10.]
    p_va = get_circle(Matrix(va2P'),p0_va)
    ~,vaMaxIdx = findmax(vec(
    sum(
    abs2,va2P.-
    reshape(va2P[1,:],(1,2)),
    # mean(αVel2P,dims=1),
    dims=2)))

    # avg_B = mean(sqrt.(sum(abs2,mag_rtn,dims=2)))
    # mag1_interp_linear = LinearInterpolation(magEpoch, mag2P[:,1]/avg_B)
    # mag1_pTime = mag1_interp_linear.(pEpoch)
    # ccpb = cor(mag1_pTime, pVel2P[:,1])
    # coefspb = linear_fit(mag1_pTime,pVel2P[:,1])# pVel1 = c + Vw*mag1_pTime
    p0_p = [0.,0.,30.]
    p_p = get_circle(Matrix(pVel2P'),p0_p)
    # p_p = get_circle_rfixed(Matrix(pVel2P'),p0_p,coefspb[2])
    # push!(p_p,coefspb[2])
    # output["Cp"] = (p_p[1],p_p[2])
    # output["rp"] = p_p[3]
    # ~,pMaxIdx = findmax(vec(
    # sum(
    # abs2,pVel2P.-
    # # mean(pVel2P,dims=1),
    # reshape(pVel2P[1,:],(1,2)),
    # dims=2)))
    # @show pMaxIdx

    αVel2P = [(v1'*αVel')' (v2'*αVel')']
    # p0_α = p_p[1:2]
    p0_α = p_p
    # mag1_αTime = mag1_interp_linear.(αEpoch)
    # ccαb = cor(mag1_αTime, αVel2P[:,1])
    # coefsαb = linear_fit(mag1_αTime,αVel2P[:,1])# pVel1 = c + Vw*mag1_pTime
    p_α = get_circle(Matrix(αVel2P'),p0_α)
    # p_α = get_circle_rfixed(Matrix(αVel2P'),p0_α,coefsαb[2])
    # push!(p_α,coefsαb[2])
    # output["Cα"] = (p_α[1],p_α[2])
    # output["rα"] = p_α[3]
    # ~,αMaxIdx = findmax(vec(
    # sum(
    # abs2,αVel2P.-
    # reshape(αVel2P[1,:],(1,2)),
    # # mean(αVel2P,dims=1),
    # dims=2)))
    # @show αMaxIdx

    mag2P1 = [(v1_1'*mag_rtn1')' (v2_1'*mag_rtn1')']
    p0_mag1 = [0.,0.,10.]
    p_mag1 = get_circle(Matrix(mag2P1'),p0_mag1)
    pVel2P1 = [(v1_1'*pVel1')' (v2_1'*pVel1')']
    pVel2Pα1 = [(v1_1'*pVel2αEpoch1')' (v2_1'*pVel2αEpoch1')']
    va2P1 = [(v1_1'*va_rtn1')' (v2_1'*va_rtn1')']
    p0_va1 = [0.,0.,10.]
    p_va1 = get_circle(Matrix(va2P1'),p0_va1)
    ~,vaMaxIdx1 = findmax(vec(
    sum(
    abs2,va2P1.-
    reshape(va2P1[1,:],(1,2)),
    # mean(αVel2P,dims=1),
    dims=2)))
    p0_p1 = [0.,0.,30.]
    p_p1 = get_circle(Matrix(pVel2P1'),p0_p1)
    αVel2P1 = [(v1_1'*αVel1')' (v2_1'*αVel1')']
    p0_α1 = p_p1
    p_α1 = get_circle(Matrix(αVel2P1'),p0_α1)


    mag2P2 = [(v1_2'*mag_rtn2')' (v2_2'*mag_rtn2')']
    p0_mag2 = [0.,0.,10.]
    p_mag2 = get_circle(Matrix(mag2P2'),p0_mag2)
    pVel2P2 = [(v1_2'*pVel2')' (v2_2'*pVel2')']
    pVel2Pα2 = [(v1_2'*pVel2αEpoch2')' (v2_2'*pVel2αEpoch2')']
    va2P2 = [(v1_2'*va_rtn2')' (v2_2'*va_rtn2')']
    p0_va2 = [0.,0.,10.]
    p_va2 = get_circle(Matrix(va2P2'),p0_va2)
    ~,vaMaxIdx2 = findmax(vec(
    sum(
    abs2,va2P2.-
    reshape(va2P2[1,:],(1,2)),
    # mean(αVel2P,dims=1),
    dims=2)))
    p0_p2 = [0.,0.,30.]
    p_p2 = get_circle(Matrix(pVel2P2'),p0_p2)
    αVel2P2 = [(v1_2'*αVel2')' (v2_2'*αVel2')']
    p0_α2 = p_p2
    p_α2 = get_circle(Matrix(αVel2P2'),p0_α2)

    makerSize = 3

    # 画速度、磁场相关系数
    # sign2sign = Dict(1=>"+",-1=>"-")
    # avg_B = mean(sqrt.(sum(abs2,mag_rtn,dims=2)))
    # mag1_interp_linear = LinearInterpolation(magEpoch, mag2P[:,1]/avg_B)
    # mag1_pTime = mag1_interp_linear.(pEpoch)
    # ccpb = cor(mag1_pTime, pVel2P[:,1])
    # coefspb = linear_fit(mag1_pTime,pVel2P[:,1])# pVel1 = c + Vw*mag1_pTime
    # output["ccpb"] = ccpb
    # output["Vwp"] = coefspb[2]
    # plot_pvelVsB = scatter(
    # mag1_pTime,
    # pVel2P[:,1],
    # label=:none,
    # xlabel = "B1/|B|",
    # ylabel = "Vp1 Km/s",
    # title = "Vp1 = $(round(coefspb[1],digits=3))"*
    # sign2sign[sign(coefspb[2])]*
    # "$(round(abs(coefspb[2]),digits=3))*B1/|B|",
    # )
    # pb_func(x) = coefspb[1]+coefspb[2]*x
    # plot!(
    # plot_pvelVsB,
    # pb_func,
    # ls = :dot,
    # color = :red,
    # label = "cc="*string(round(ccpb,digits=3)),
    # )
    # mag1_αTime = mag1_interp_linear.(αEpoch)
    # ccαb = cor(mag1_αTime, αVel2P[:,1])
    # coefsαb = linear_fit(mag1_αTime,αVel2P[:,1])# pVel1 = c + Vw*mag1_pTime
    # output["ccαb"] = ccαb
    # output["Vwα"] = coefsαb[2]
    # plot_αvelVsB = scatter(
    # mag1_αTime,
    # αVel2P[:,1],
    # label=:none,
    # xlabel = "B1/|B|",
    # ylabel = "Vα1 Km/s",
    # title = "Vα1 = $(round(coefsαb[1],digits=3))"*
    # sign2sign[sign(coefsαb[2])]*
    # "$(round(abs(coefsαb[2]),digits=3))*B1/|B|",
    # )
    # αb_func(x) = coefsαb[1]+coefsαb[2]*x
    # plot!(
    # plot_αvelVsB,
    # αb_func,
    # ls = :dot,
    # color = :red,
    # label = "cc="*string(round(ccαb,digits=3)),
    # )
    # velVsBL = @layout grid(1,2)
    # plot_velVsB = plot(
    # plot_pvelVsB,
    # plot_αvelVsB,
    # layout=velVsBL,
    # size=(800,450),
    # )
    #
    # savefig(plot_velVsB,"figure\\labeledSbEvents\\"*figDir*"\\velVsB\\"*figName*".png")

    # 画投影的速度和磁场
    # plot_vel = scatter(
    # pVel2P[:,1],
    # pVel2P[:,2],
    # xlabel = "V1 Km/s",
    # ylabel = "V2 Km/s",
    # label = "Vp",
    # color = :blue,
    # ms = makerSize,
    # )
    # scatter!(
    # plot_vel,
    # [pVel2P[1,1],],
    # [pVel2P[1,2],],
    # color = :blue,
    # markershape = :square,
    # ms = makerSize+2,
    # label=:none,
    # )
    # scatter!(
    # plot_vel,
    # [pVel2P[pMaxIdx,1],],
    # [pVel2P[pMaxIdx,2],],
    # color = :blue,
    # markershape = :utriangle,
    # ms = makerSize+2,
    # label=:none,
    # )
    # scatter!(
    # plot_vel,
    # [p_p[1],],
    # [p_p[2],],
    # color = :blue,
    # markershape = :x,
    # ms = makerSize+2,
    # label=:none,
    # )
    # plot!(
    # plot_vel,
    # circleShape(p_p),
    # lw= 0.5,
    # linecolor = :blue,
    # linealpha = 0.5,
    # label=:none,
    # )
    # scatter!(
    # plot_vel,
    # αVel2P[:,1],
    # αVel2P[:,2],
    # label = "Vα",
    # color = :red,
    # ms = makerSize,
    # )
    # scatter!(
    # plot_vel,
    # [αVel2P[1,1],],
    # [αVel2P[1,2],],
    # color = :red,
    # markershape = :square,
    # ms = makerSize+2,
    # label=:none,
    # )
    # scatter!(
    # plot_vel,
    # [αVel2P[αMaxIdx,1],],
    # [αVel2P[αMaxIdx,2],],
    # color = :red,
    # markershape = :utriangle,
    # label=:none,
    # ms = makerSize+2,
    # aspect_ratio = :equal,
    # )
    # scatter!(
    # plot_vel,
    # [p_α[1],],
    # [p_α[2],],
    # color = :red,
    # markershape = :x,
    # ms = makerSize+2,
    # label=:none,
    # )
    # plot!(
    # plot_vel,
    # circleShape(p_α),
    # lw= 0.5,
    # linecolor = :red,
    # linealpha = 0.5,
    # label=:none,
    # )
    # # savefig("figure\\sbvel2P\\"*figName*".png")
    # #################
    # #
    # plot_mag = scatter(
    # mag2P[:,1],
    # mag2P[:,2],
    # xlabel = "B1 nT",
    # ylabel = "B2 nT",
    # ms = 1,
    # color = :black,
    # legend = false,
    # aspect_ratio = :equal,
    # )
    # scatter!(
    # plot_mag,
    # [p_mag[1],],
    # [p_mag[2],],
    # color = :blue,
    # markershape = :x,
    # ms = makerSize+2,
    # )
    # plot!(
    # plot_mag,
    # circleShape(p_mag),
    # lw= 0.5,
    # linecolor = :black,
    # linealpha = 0.5,
    # )
    # plot2PL = @layout grid(1,2)
    # plot(
    # plot_vel,
    # plot_mag,
    # layout=plot2PL,
    # size=(800,450),
    # )
    # savefig("figure\\sbevents\\alfvenic\\VelandB2P\\"*figName*".png")
    # savefig("figure\\B2P\\"*figName*".png")
    # F.U 给出svd分解得到的特征向量
    # F.S 给出svd分解得到的特征值，降序排列
    #@show F.U
    #@show F.S
    #@show F.V

    # 画出投影的速度和VA
    function plotVel2P(figsize=(1100,1200))
        fig_velandVA = []
        # 1
        plot_vel = scatter(
        pVel2P1[:,1],
        pVel2P1[:,2],
        xlabel = "V1 Km/s",
        ylabel = "V2 Km/s",
        label = "Vp",
        color = :blue,
        ms = makerSize,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [pVel2P1[1,1],],
        [pVel2P1[1,2],],
        color = :blue,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [pVel2Pα1[vaMaxIdx1,1],],
        [pVel2Pα1[vaMaxIdx1,2],],
        color = :blue,
        markershape = :utriangle,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [p_p1[1],],
        [p_p1[2],],
        color = :blue,
        markershape = :x,
        ms = makerSize+2,
        label=:none,
        )
        plot!(
        plot_vel,
        circleShape(p_p1),
        lw= 0.5,
        linecolor = :blue,
        linealpha = 0.5,
        label=:none,
        )
        scatter!(
        plot_vel,
        αVel2P1[:,1],
        αVel2P1[:,2],
        label = "Vα",
        color = :red,
        ms = makerSize,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [αVel2P1[1,1],],
        [αVel2P1[1,2],],
        color = :red,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [αVel2P1[vaMaxIdx1,1],],
        [αVel2P1[vaMaxIdx1,2],],
        color = :red,
        markershape = :utriangle,
        label=:none,
        ms = makerSize+2,
        aspect_ratio = :equal,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [p_α1[1],],
        [p_α1[2],],
        color = :red,
        markershape = :x,
        ms = makerSize+2,
        label=:none,
        )
        plot!(
        plot_vel,
        circleShape(p_α1),
        lw= 0.5,
        linecolor = :red,
        linealpha = 0.5,
        label=:none,
        )

        plot_va = scatter(
        va2P1[:,1],
        va2P1[:,2],
        xlabel = "V1 Km/s",
        ylabel = "V2 Km/s",
        ms = makerSize,
        color = :black,
        label = "VA",
        aspect_ratio = :equal,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [p_va1[1],],
        [p_va1[2],],
        color = :black,
        markershape = :x,
        ms = makerSize+2,
        label = :none,
        )
        plot!(
        plot_va,
        circleShape(p_va1),
        lw= 0.5,
        linecolor = :black,
        linealpha = 0.5,
        label = :none,
        )
        scatter!(
        plot_va,
        [va2P1[1,1],],
        [va2P1[1,2],],
        color = :black,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [va2P1[vaMaxIdx1,1],],
        [va2P1[vaMaxIdx1,2],],
        color = :black,
        markershape = :utriangle,
        label=:none,
        ms = makerSize+2,
        aspect_ratio = :equal,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        αVel2P1[:,1]-pVel2Pα1[:,1],
        αVel2P1[:,2]-pVel2Pα1[:,2],
        xlabel = "V1 Km/s",
        ylabel = "V2 Km/s",
        label = "Vαp",
        color = :blue,
        ms = makerSize,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [αVel2P1[1,1]-pVel2Pα1[1,1],],
        [αVel2P1[1,2]-pVel2Pα1[1,2],],
        color = :blue,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [αVel2P1[vaMaxIdx1,1]-pVel2Pα1[vaMaxIdx1,1],],
        [αVel2P1[vaMaxIdx1,2]-pVel2Pα1[vaMaxIdx1,2],],
        color = :blue,
        markershape = :utriangle,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        push!(fig_velandVA,plot_vel,plot_va)

        plot_vel = scatter(
        pVel2P[:,1],
        pVel2P[:,2],
        xlabel = "V1 Km/s",
        ylabel = "V2 Km/s",
        label = "Vp",
        color = :blue,
        ms = makerSize,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [pVel2P[1,1],],
        [pVel2P[1,2],],
        color = :blue,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [pVel2Pα[vaMaxIdx,1],],
        [pVel2Pα[vaMaxIdx,2],],
        color = :blue,
        markershape = :utriangle,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [p_p[1],],
        [p_p[2],],
        color = :blue,
        markershape = :x,
        ms = makerSize+2,
        label=:none,
        )
        plot!(
        plot_vel,
        circleShape(p_p),
        lw= 0.5,
        linecolor = :blue,
        linealpha = 0.5,
        label=:none,
        )
        scatter!(
        plot_vel,
        αVel2P[:,1],
        αVel2P[:,2],
        label = "Vα",
        color = :red,
        ms = makerSize,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [αVel2P[1,1],],
        [αVel2P[1,2],],
        color = :red,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [αVel2P[vaMaxIdx,1],],
        [αVel2P[vaMaxIdx,2],],
        color = :red,
        markershape = :utriangle,
        label=:none,
        ms = makerSize+2,
        aspect_ratio = :equal,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [p_α[1],],
        [p_α[2],],
        color = :red,
        markershape = :x,
        ms = makerSize+2,
        label=:none,
        )
        plot!(
        plot_vel,
        circleShape(p_α),
        lw= 0.5,
        linecolor = :red,
        linealpha = 0.5,
        label=:none,
        )

        plot_va = scatter(
        va2P[:,1],
        va2P[:,2],
        xlabel = "V1 Km/s",
        ylabel = "V2 Km/s",
        ms = makerSize,
        color = :black,
        label = "VA",
        aspect_ratio = :equal,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [p_va[1],],
        [p_va[2],],
        color = :black,
        markershape = :x,
        ms = makerSize+2,
        label = :none,
        )
        plot!(
        plot_va,
        circleShape(p_va),
        lw= 0.5,
        linecolor = :black,
        linealpha = 0.5,
        label = :none,
        )
        scatter!(
        plot_va,
        [va2P[1,1],],
        [va2P[1,2],],
        color = :black,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [va2P[vaMaxIdx,1],],
        [va2P[vaMaxIdx,2],],
        color = :black,
        markershape = :utriangle,
        label=:none,
        ms = makerSize+2,
        aspect_ratio = :equal,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        αVel2P[:,1]-pVel2Pα[:,1],
        αVel2P[:,2]-pVel2Pα[:,2],
        xlabel = "V1 Km/s",
        ylabel = "V2 Km/s",
        label = "Vαp",
        color = :blue,
        ms = makerSize,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [αVel2P[1,1]-pVel2Pα[1,1],],
        [αVel2P[1,2]-pVel2Pα[1,2],],
        color = :blue,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [αVel2P[vaMaxIdx,1]-pVel2Pα[vaMaxIdx,1],],
        [αVel2P[vaMaxIdx,2]-pVel2Pα[vaMaxIdx,2],],
        color = :blue,
        markershape = :utriangle,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        push!(fig_velandVA,plot_vel,plot_va)

        plot_vel = scatter(
        pVel2P2[:,1],
        pVel2P2[:,2],
        xlabel = "V1 Km/s",
        ylabel = "V2 Km/s",
        label = "Vp",
        color = :blue,
        ms = makerSize,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [pVel2P2[1,1],],
        [pVel2P2[1,2],],
        color = :blue,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [pVel2Pα2[vaMaxIdx2,1],],
        [pVel2Pα2[vaMaxIdx2,2],],
        color = :blue,
        markershape = :utriangle,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [p_p2[1],],
        [p_p2[2],],
        color = :blue,
        markershape = :x,
        ms = makerSize+2,
        label=:none,
        )
        plot!(
        plot_vel,
        circleShape(p_p2),
        lw= 0.5,
        linecolor = :blue,
        linealpha = 0.5,
        label=:none,
        )
        scatter!(
        plot_vel,
        αVel2P2[:,1],
        αVel2P2[:,2],
        label = "Vα",
        color = :red,
        ms = makerSize,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [αVel2P2[1,1],],
        [αVel2P2[1,2],],
        color = :red,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [αVel2P2[vaMaxIdx2,1],],
        [αVel2P2[vaMaxIdx2,2],],
        color = :red,
        markershape = :utriangle,
        label=:none,
        ms = makerSize+2,
        aspect_ratio = :equal,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_vel,
        [p_α2[1],],
        [p_α2[2],],
        color = :red,
        markershape = :x,
        ms = makerSize+2,
        label=:none,
        )
        plot!(
        plot_vel,
        circleShape(p_α2),
        lw= 0.5,
        linecolor = :red,
        linealpha = 0.5,
        label=:none,
        )

        plot_va = scatter(
        va2P2[:,1],
        va2P2[:,2],
        xlabel = "V1 Km/s",
        ylabel = "V2 Km/s",
        ms = makerSize,
        color = :black,
        label = "VA",
        aspect_ratio = :equal,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [p_va2[1],],
        [p_va2[2],],
        color = :black,
        markershape = :x,
        ms = makerSize+2,
        label = :none,
        )
        plot!(
        plot_va,
        circleShape(p_va2),
        lw= 0.5,
        linecolor = :black,
        linealpha = 0.5,
        label = :none,
        )
        scatter!(
        plot_va,
        [va2P2[1,1],],
        [va2P2[1,2],],
        color = :black,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [va2P2[vaMaxIdx2,1],],
        [va2P2[vaMaxIdx2,2],],
        color = :black,
        markershape = :utriangle,
        label=:none,
        ms = makerSize+2,
        aspect_ratio = :equal,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        αVel2P2[:,1]-pVel2Pα2[:,1],
        αVel2P2[:,2]-pVel2Pα2[:,2],
        xlabel = "V1 Km/s",
        ylabel = "V2 Km/s",
        label = "Vαp",
        color = :blue,
        ms = makerSize,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [αVel2P2[1,1]-pVel2Pα2[1,1],],
        [αVel2P2[1,2]-pVel2Pα2[1,2],],
        color = :blue,
        markershape = :square,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        scatter!(
        plot_va,
        [αVel2P2[vaMaxIdx2,1]-pVel2Pα2[vaMaxIdx2,1],],
        [αVel2P2[vaMaxIdx2,2]-pVel2Pα2[vaMaxIdx2,2],],
        color = :blue,
        markershape = :utriangle,
        ms = makerSize+2,
        label=:none,
        markerstrokewidth = 0,
        )
        push!(fig_velandVA,plot_vel,plot_va)

        plot2PL = @layout grid(3,2)
        plot(
        fig_velandVA...,
        layout=plot2PL,
        size=figsize,
        )
        # savefig("figure\\sbevents\\alfvenic\\VelandVaanddV2P\\"*figName*".png")
        savefig("figure\\labeledSbEvents\\"*figDir*"\\VelandVaanddV2P\\"*figName*".png")
    end
    # plotVel2P()


    # 画出rtn坐标系中的矢量速度数据
    function plotVecRtn(figsize=(1000,1000))
        plot_p = scatter3d(
        pVel1[:,1],
        pVel1[:,2],
        pVel1[:,3],
        label = "1",
        xlabel = "Vr Km/s",
        ylabel = "Vt Km/s",
        zlabel = "Vn Km/s",
        color = :blue,
        ms = makerSize,
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        aspect_ratio = :equal,
        )
        scatter3d!(
        plot_p,
        pVel[:,1],
        pVel[:,2],
        pVel[:,3],
        label = "2",
        color = :red,
        ms = makerSize,
        title = "Vp",
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        )
        scatter3d!(
        plot_p,
        pVel2[:,1],
        pVel2[:,2],
        pVel2[:,3],
        label = "3",
        color = :green,
        ms = makerSize,
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        )
        qVec2 = v3*60
        quiver!(
        plot_p,
        [mean(pVel[:,1]),],
        [mean(pVel[:,2]),],
        [mean(pVel[:,3]),],
        quiver = ([qVec2[1],],[qVec2[2],],[qVec2[3],]),
        color = :red,
        w = 4,
        )
        # qVec1 = v3_1*35
        # quiver!(
        # plot_p,
        # [mean(pVel1[:,1]),],
        # [mean(pVel1[:,2]),],
        # [mean(pVel1[:,3]),],
        # quiver = ([qVec1[1],],[qVec1[2],],[qVec1[3],]),
        # color = :blue,
        # w = 4,
        # )
        # qVec3 = v3_2*35
        # quiver!(
        # plot_p,
        # [mean(pVel2[:,1]),],
        # [mean(pVel2[:,2]),],
        # [mean(pVel2[:,3]),],
        # quiver = ([qVec3[1],],[qVec3[2],],[qVec3[3],]),
        # color = :green,
        # w = 4,
        # )
        plot_α = scatter3d(
        αVel1[:,1],
        αVel1[:,2],
        αVel1[:,3],
        label = "1",
        xlabel = "Vr Km/s",
        ylabel = "Vt Km/s",
        zlabel = "Vn Km/s",
        color = :blue,
        ms = makerSize,
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        aspect_ratio = :equal,
        )
        scatter3d!(
        plot_α,
        αVel[:,1],
        αVel[:,2],
        αVel[:,3],
        label = "2",
        color = :red,
        ms = makerSize,
        title = "Vα",
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        )
        scatter3d!(
        plot_α,
        αVel2[:,1],
        αVel2[:,2],
        αVel2[:,3],
        label = "3",
        color = :green,
        ms = makerSize,
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        )
        quiver!(
        plot_α,
        [mean(αVel[:,1]),],
        [mean(αVel[:,2]),],
        [mean(αVel[:,3]),],
        quiver = ([qVec2[1],],[qVec2[2],],[qVec2[3],]),
        color = :red,
        w = 4,
        )
        # quiver!(
        # plot_α,
        # [mean(αVel1[:,1]),],
        # [mean(αVel1[:,2]),],
        # [mean(αVel1[:,3]),],
        # quiver = ([qVec1[1],],[qVec1[2],],[qVec1[3],]),
        # color = :blue,
        # w = 4,
        # )
        # quiver!(
        # plot_α,
        # [mean(αVel2[:,1]),],
        # [mean(αVel2[:,2]),],
        # [mean(αVel2[:,3]),],
        # quiver = ([qVec3[1],],[qVec3[2],],[qVec3[3],]),
        # color = :green,
        # w = 4,
        # )
        plot_αp = scatter3d(
        αVel1[:,1].-pVel2αEpoch1[:,1],
        αVel1[:,2].-pVel2αEpoch1[:,2],
        αVel1[:,3].-pVel2αEpoch1[:,3],
        label = "1",
        xlabel = "Vr Km/s",
        ylabel = "Vt Km/s",
        zlabel = "Vn Km/s",
        color = :blue,
        ms = makerSize,
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        aspect_ratio = :equal,
        )
        scatter3d!(
        plot_αp,
        αVel[:,1].-pVel2αEpoch[:,1],
        αVel[:,2].-pVel2αEpoch[:,2],
        αVel[:,3].-pVel2αEpoch[:,3],
        label = "2",
        color = :red,
        ms = makerSize,
        title = "Vαp",
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        )
        scatter3d!(
        plot_αp,
        αVel2[:,1].-pVel2αEpoch2[:,1],
        αVel2[:,2].-pVel2αEpoch2[:,2],
        αVel2[:,3].-pVel2αEpoch2[:,3],
        label = "3",
        color = :green,
        ms = makerSize,
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        )
        quiver!(
        plot_αp,
        [mean(αVel[:,1].-pVel2αEpoch[:,1]),],
        [mean(αVel[:,2].-pVel2αEpoch[:,2]),],
        [mean(αVel[:,3].-pVel2αEpoch[:,3]),],
        quiver = ([qVec2[1],],[qVec2[2],],[qVec2[3],]),
        color = :red,
        w = 4,
        )
        plot_va = scatter3d(
        va_rtn1[:,1],
        va_rtn1[:,2],
        va_rtn1[:,3],
        label = "1",
        xlabel = "Vr Km/s",
        ylabel = "Vt Km/s",
        zlabel = "Vn Km/s",
        color = :blue,
        ms = makerSize,
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        aspect_ratio = :equal,
        )
        scatter3d!(
        plot_va,
        va_rtn[:,1],
        va_rtn[:,2],
        va_rtn[:,3],
        label = "2",
        color = :red,
        ms = makerSize,
        title = "VA",
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        )
        scatter3d!(
        plot_va,
        va_rtn2[:,1],
        va_rtn2[:,2],
        va_rtn2[:,3],
        label = "3",
        color = :green,
        ms = makerSize,
        markershape = :circle,
        markerstrokewidth = 0,
        markeralpha = 0.6,
        )
        quiver!(
        plot_va,
        [mean(va_rtn[:,1]),],
        [mean(va_rtn[:,2]),],
        [mean(va_rtn[:,3]),],
        quiver = ([qVec2[1],],[qVec2[2],],[qVec2[3],]),
        color = :red,
        w = 4,
        )
        sbvelL = @layout grid(2,2)
        plot(plot_p,plot_α,plot_αp,plot_va,layout=sbvelL,size=figsize)
        savefig("figure\\labeledSbEvents\\"*figDir*"\\sbvel\\"*figName*".png")
    end
    # plotVecRtn()


    # # 画出时间序列
    # deltaEpoch = deltaMins/(24*60)
    # epoch1_sb = copy(epoch1)
    # epoch2_sb = copy(epoch2)
    # epoch1 -= deltaEpoch
    # epoch2 += deltaEpoch
    # pPoints = vec((pVars["p_epoch"].>=epoch1) .& (pVars["p_epoch"].<=epoch2))
    # pEpoch = pVars["p_epoch"][pPoints]
    # pTime = epoch2datetime.(pEpoch)
    # pVel = pVars["p_vel_rtn_sun"][pPoints,:]
    # pTemp = pVars["p_temp"][pPoints]
    # αPoints =  vec((αVars["alpha_epoch"].>=epoch1) .&
    #             (αVars["alpha_epoch"].<=epoch2))
    # αEpoch = αVars["alpha_epoch"][αPoints]
    # αTime = epoch2datetime.(αEpoch)
    # αVel = αVars["alpha_vel_rtn_sun"][αPoints,:]
    # αTemp = αVars["alpha_temp"][αPoints]
    # va = modifiedVars["va_alphaEpoch"][αPoints]
    # va_rtn = modifiedVars["va_rtn_alphaEpoch"][αPoints,:]
    # magPoints = vec((magVars["mag_epoch"].>=epoch1) .&
    #             (magVars["mag_epoch"].<=epoch2))
    # magEpoch = magVars["mag_epoch"][magPoints]
    # magTime = epoch2datetime.(magEpoch)
    # mag_rtn = magVars["mag_rtn"][magPoints,:]
    # theta = atan.(mag_rtn[:,2]./mag_rtn[:,1])
    # if plotTimeSeries
    #     ts = []
    #     p1 = plot(
    #     pTime,
    #     pVel,
    #     label = ["r" "t" "n"],
    #     ylabel = "Vp Km/s",
    #     xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    #     xticks = nothing,
    #     )
    #     p2 = plot(
    #     pTime,
    #     pTemp,
    #     legend=false,
    #     ylabel = "Tp eV",
    #     xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    #     xticks = nothing,
    #     )
    #     p3 = plot(
    #     αTime,
    #     αVel,
    #     label = ["r" "t" "n"],
    #     ylabel = "Vα Km/s",
    #     xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    #     xticks = nothing,
    #     )
    #     p4 = plot(
    #     αTime,
    #     αTemp,
    #     legend=false,
    #     ylabel = "Tα eV",
    #     xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    #     xticks = nothing,
    #     )
    #     p5 = plot(
    #     αTime,
    #     va_rtn,
    #     label = ["r" "t" "n"],
    #     ylabel = "VA Km/s",
    #     xticks = nothing,
    #     xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    #     )
    #     p6 = plot(
    #     magTime,
    #     mag_rtn,
    #     label = ["r" "t" "n"],
    #     ylabel = "B nT",
    #     xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    #     xticks = nothing,
    #     )
    #     p7 = plot(
    #     magTime,
    #     theta*180/π,
    #     legend = false,
    #     ylabel = "θ °",
    #     ylims = (-90,90),
    #     xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    #     )
    #     plot!(
    #     p7,
    #     [epoch2datetime(epoch1_sb),epoch2datetime(epoch1_sb)],
    #     [-90,90],
    #     ls = :dot,
    #     lw = 2,
    #     color = :black,
    #     )
    #     plot!(
    #     p7,
    #     [epoch2datetime(epoch2_sb),epoch2datetime(epoch2_sb)],
    #     [-90,90],
    #     ls = :dot,
    #     lw = 2,
    #     color = :black,
    #     )
    #     push!(ts,p1,p2,p3,p4,p5,p6,p7)
    #     tsL = @layout grid(7,1)
    #     plot(ts...,layout=tsL,size=(900,1050))
    #     #size=(900,600),
    #     savefig("figure\\labeledSbEvents\\"*figDir*"\\sbtimeplots\\"*figName*".png")
    # end
    output
end

function sbEvent_timeplot(epoch1,epoch2,pVars,αVars,modifiedVars,magVars,vαp,θvαp_va;
    figName="test",deltaMins=180,figDir="issb")
    # 画出时间序列
    deltaEpoch = deltaMins/(24*60)
    epoch1_sb = copy(epoch1)
    epoch2_sb = copy(epoch2)
    epoch1 -= deltaEpoch
    epoch2 += deltaEpoch
    pPoints = vec((pVars["p_epoch"].>=epoch1) .& (pVars["p_epoch"].<=epoch2))
    pEpoch = pVars["p_epoch"][pPoints]
    pTime = epoch2datetime.(pEpoch)
    pVel = pVars["p_vel_rtn_sun"][pPoints,:]
    pTemp = pVars["p_temp"][pPoints]
    αPoints =  vec((αVars["alpha_epoch"].>=epoch1) .&
                (αVars["alpha_epoch"].<=epoch2))
    αEpoch = αVars["alpha_epoch"][αPoints]
    αTime = epoch2datetime.(αEpoch)
    αVel = αVars["alpha_vel_rtn_sun"][αPoints,:]
    αTemp = αVars["alpha_temp"][αPoints]
    va = modifiedVars["va_alphaEpoch"][αPoints]
    va_rtn = modifiedVars["va_rtn_alphaEpoch"][αPoints,:]
    vαp_rtn2α = αVel .- modifiedVars["p_vel_rtn_sun_alphaEpoch"][αPoints,:]
    vαp = vαp[αPoints]
    magPoints = vec((magVars["mag_epoch"].>=epoch1) .&
                (magVars["mag_epoch"].<=epoch2))
    magEpoch = magVars["mag_epoch"][magPoints]
    magTime = epoch2datetime.(magEpoch)
    mag_rtn = magVars["mag_rtn"][magPoints,:]
    theta = atan.(mag_rtn[:,2]./mag_rtn[:,1])
    θvαp_va = θvαp_va[αPoints]
    ts = []
    p1 = plot(
    pTime,
    pVel,
    label = ["r" "t" "n"],
    ylabel = "Vp Km/s",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    xticks = nothing,
    )
    p2 = plot(
    magTime,
    mag_rtn,
    label = ["r" "t" "n"],
    ylabel = "B nT",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    xticks = nothing,
    )
    p3 = plot(
    αTime,
    vαp_rtn2α,
    label = ["r" "t" "n"],
    ylabel = "Vαp Km/s",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    xticks = nothing,
    )
    p4 = plot(
    αTime,
    αVel,
    label = ["r" "t" "n"],
    ylabel = "Vα Km/s",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    xticks = nothing,
    )
    p5 = plot(
    αTime,
    vαp,
    label = :none,
    ylabel = "Vαp Km/s",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    xticks = nothing,
    )
    p6 = plot(
    pTime,
    pTemp,
    legend=false,
    ylabel = "Tp eV",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    xticks = nothing,
    )
    p7 = plot(
    αTime,
    αTemp,
    legend=false,
    ylabel = "Tα eV",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    xticks = nothing,
    )
    # p5 = plot(
    # αTime,
    # va_rtn,
    # label = ["r" "t" "n"],
    # ylabel = "VA Km/s",
    # xticks = nothing,
    # xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    # )
    p8 = plot(
    αTime,
    θvαp_va*180/π,
    legend=false,
    ylabel = "θ(Vαp,B) °",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    ylims = (0,180),
    xticks = nothing,
    )
    p9 = plot(
    magTime,
    theta*180/π,
    legend = false,
    ylabel = "θ °",
    ylims = (-90,90),
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    )
    plot!(
    p9,
    [epoch2datetime(epoch1_sb),epoch2datetime(epoch1_sb)],
    [-90,90],
    ls = :dot,
    lw = 2,
    color = :black,
    )
    plot!(
    p9,
    [epoch2datetime(epoch2_sb),epoch2datetime(epoch2_sb)],
    [-90,90],
    ls = :dot,
    lw = 2,
    color = :black,
    )
    push!(ts,p1,p2,p3,p4,p5,p6,p7,p8,p9)
    tsL = @layout grid(9,1)
    plot(ts...,layout=tsL,size=(800,1050))
    #size=(900,600),
    savefig("figure\\labeledSbEvents\\"*figDir*"\\sbtimeplots_1\\"*figName*".png")
end

function sbEvent_calVectors(epoch1,epoch2,pVars,αVars,modifiedVars,magVars;
    figName="test",plotTimeSeries=false,deltaMins=10,figDir="issb")
    # figDir2sbtype = Dict("issb"=>1, "maybesb"=>2, "notsb"=>3)
    output = Dict(
    "αdataFlag"=>1,
    "e_δb"=>[NaN,NaN,NaN],
    "B0"=>[NaN,NaN,NaN],
    # "ccpb"=>NaN,
    # "ccαb"=>NaN,
    # "Vwp"=>NaN,
    # "rp"=>NaN,
    # "Vwα"=>NaN,
    # "rα"=>NaN,
    # "Cp"=>(NaN,NaN),
    # "Cα"=>(NaN,NaN),
    # "sbtype"=>figDir2sbtype[figDir], #1:issb,2:maybesb,3:notsb
    )
    deltaEpoch = deltaMins/(24*60)
    deltaEpochSmall = 1/(24*60)
    #############
    pPoints = vec((pVars["p_epoch"].>=epoch1) .& (pVars["p_epoch"].<=epoch2))
    pEpoch = pVars["p_epoch"][pPoints]
    pTime = epoch2datetime.(pEpoch)
    pVel = pVars["p_vel_rtn_sun"][pPoints,:]
    pTemp = pVars["p_temp"][pPoints]
    αPoints =  vec((αVars["alpha_epoch"].>=epoch1) .&
                (αVars["alpha_epoch"].<=epoch2))
    αEpoch = αVars["alpha_epoch"][αPoints]
    αTime = epoch2datetime.(αEpoch)
    αVel = αVars["alpha_vel_rtn_sun"][αPoints,:]
    pVel2αEpoch = modifiedVars["p_vel_rtn_sun_alphaEpoch"][αPoints,:]
    αTemp = αVars["alpha_temp"][αPoints]
    va = modifiedVars["va_alphaEpoch"][αPoints]
    va_rtn = modifiedVars["va_rtn_alphaEpoch"][αPoints,:]
    magPoints = vec((magVars["mag_epoch"].>=epoch1) .&
                (magVars["mag_epoch"].<=epoch2))
    magEpoch = magVars["mag_epoch"][magPoints]
    magTime = epoch2datetime.(magEpoch)
    mag_rtn = magVars["mag_rtn"][magPoints,:]
    # 判断这段时间是否有足够的阿尔法粒子数据
    if (length(αVel)-sum(isnan.(αVel))) < 10
        output["αdataFlag"] = 0
        return output
    end
    nonNaNpoints = vec(any(isnan.(mag_rtn),dims=2) .!= 1)
    magEpoch = magEpoch[nonNaNpoints]
    magTime = magTime[nonNaNpoints]
    mag_rtn = mag_rtn[nonNaNpoints,:]
    nonNaNpoints = vec(any(isnan.(pVel),dims=2) .!= 1)
    pEpoch = pEpoch[nonNaNpoints]
    pTime = pTime[nonNaNpoints]
    pVel = pVel[nonNaNpoints,:]
    pTemp = pTemp[nonNaNpoints]
    nonNaNpoints = vec(any(isnan.(αVel),dims=2) .!= 1)
    αEpoch = αEpoch[nonNaNpoints]
    αTime = αTime[nonNaNpoints]
    αVel = αVel[nonNaNpoints,:]
    αTemp = αTemp[nonNaNpoints]
    pVel2αEpoch = pVel2αEpoch[nonNaNpoints,:]
    va = va[nonNaNpoints]
    va_rtn = va_rtn[nonNaNpoints,:]

    ### 计算这段时间的磁场最小扰动方向（MVA）
    mag_rtn_scaled = mag_rtn .- mean(mag_rtn,dims=1)
    Σb = (mag_rtn_scaled'*mag_rtn_scaled) / length(magEpoch)
    # @show Σb
    F = svd(Σb)
    v1 = F.U[:,1]
    v2 = F.U[:,2]
    v3 = F.U[:,3]
    e_δb = v1
    B0 = vec(mean(mag_rtn,dims=1))
    e_r = [1,0,0]
    output["e_δb"] = e_δb
    output["B0"] = B0
    ### 计算这段时间的漂移速度最小扰动方向（MVA）
    vαp_rtn = αVel .- pVel2αEpoch
    vαp_rtn_scaled = vαp_rtn .- mean(vαp_rtn,dims=1)
    Σαp = (vαp_rtn_scaled'*vαp_rtn_scaled) / length(αEpoch)
    if any(isnan.(Σαp))
        output["αdataFlag"] = 0
        return output
    end
    # @show Σαp
    Fαp = svd(Σαp)
    v1αp = Fαp.U[:,1]
    v2αp = Fαp.U[:,2]
    v3αp = Fαp.U[:,3]
    e_δαp = v1αp
    vαp0 = vec(mean(vαp_rtn,dims=1))
    δαp2vαp = sum(abs2,vαp_rtn_scaled) / (length(αEpoch) * sum(abs2,vαp0))
    output["e_δαp"] = e_δαp
    output["vαp0"] = vαp0
    output["δαp2vαp"] = δαp2vαp
    output["θδb_B0"] = spanAngle(e_δb,B0)
    output["θδb_r"] = spanAngle(e_δb,e_r)
    output["θδb_δαp"] = spanAngle(e_δb,e_δαp)
    output["θδb_vαp0"] = spanAngle(e_δb,vαp0)


    # @show norm(v1)
    output
end

function statSbSw(αVars,vαp,vαp2va,θvαp_va,θva,θvαp,sbEpochList;
    notSB=notSB,maybenotSB=maybenotSB,deltaMins=10)
    deltaEpoch = deltaMins/(24*60)
    deltaEpochSmall = 1/(24*60)
    function getpoints(epoch,epoch1,epoch2)
        findall(vec((epoch.>epoch1) .& (epoch.<epoch2)))
    end
    sbpoints = []
    sbpoints1 = []
    sbpoints2 = []
    maybesbpoints = []

    for sbidx in 1:size(sbEpochList)[1]
        if sbidx in notSB
            continue
        elseif sbidx in maybenotSB
            append!(maybesbpoints,
            getpoints(
            αVars["alpha_epoch"],
            sbEpochList[sbidx,1],
            sbEpochList[sbidx,2]))
        else
            append!(sbpoints,
            getpoints(
            αVars["alpha_epoch"],
            sbEpochList[sbidx,1],
            sbEpochList[sbidx,2]))
            append!(sbpoints1,
            getpoints(
            αVars["alpha_epoch"],
            sbEpochList[sbidx,1]-deltaEpoch,
            sbEpochList[sbidx,1]-deltaEpochSmall))
            append!(sbpoints2,
            getpoints(
            αVars["alpha_epoch"],
            sbEpochList[sbidx,1]-deltaEpoch,
            sbEpochList[sbidx,2]-deltaEpochSmall))
        end
    end

    # sbpoints = Int.(sbpoints)
    # maybesbpoints = Int.(maybesbpoints)
    # vαpSb = filter(!isnan,vαp[sbpoints])
    # vαp2vaSb = filter(!isnan,vαp2va[sbpoints])
    #
    # vαpSw = filter(!isnan,vαp)
    # vαp2vaSw = filter(!isnan,vαp2va)




    # histogram(
    # vαpSb,
    # normalize = :pdf,
    # label = "sb",
    # xlabel = "Vαp Km/s",
    # xlims = (-300,300),
    # alpha = 0.4,
    # )
    # histogram!(
    # vαpSw,
    # normalize = :pdf,
    # label = "sw",
    #
    # )
    # savefig("figure\\hist_vαpSb.png")
    #
    # histogram(
    # vαp2vaSb,
    # normalize = :pdf,
    # label = "sb",
    # xlabel = "Vαp/VA",
    # xlims = (-2,2),
    # alpha = 0.4,
    # )
    # histogram!(
    # vαp2vaSw,
    # normalize = :pdf,
    # label = "sw",
    # )
    # savefig("figure\\hist_vαp2vaSb.png")

    # # 统计SB内外，漂移速度与磁场的夹角
    # # θvαp_vaSb1 = filter(!isnan,θvαp_va[sbpoints1])
    # # θvαp_vaSb = filter(!isnan,θvαp_va[sbpoints])
    # # θvαp_vaSb2 = filter(!isnan,θvαp_va[sbpoints2])
    # # histogram(
    # # θvαp_vaSb1*180/π,
    # # normalize = :pdf,
    # # label = "1",
    # # xlabel = "θ(Vαp,VA) °",
    # # alpha = 0.4,
    # # )
    # # histogram!(
    # # θvαp_vaSb*180/π,
    # # normalize = :pdf,
    # # label = "2",
    # # alpha = 0.4,
    # # )
    # # histogram!(
    # # θvαp_vaSb2*180/π,
    # # normalize = :pdf,
    # # label = "3",
    # # alpha = 0.4,
    # # )
    # # savefig("figure\\hist_theta_vap_va.png")
    # 统计SB内外，漂移速度与R方向的夹角
    θvαp_1 = filter(!isnan,θvαp[sbpoints1])
    θvαp_2 = filter(!isnan,θvαp[sbpoints])
    θvαp_3 = filter(!isnan,θvαp[sbpoints2])
    histogram(
    θvαp_1*180/π,
    normalize = :pdf,
    label = "1",
    xlabel = "θ(Vαp,r) °",
    alpha = 0.4,
    )
    histogram!(
    θvαp_2*180/π,
    normalize = :pdf,
    label = "2",
    alpha = 0.4,
    )
    histogram!(
    θvαp_3*180/π,
    normalize = :pdf,
    label = "3",
    alpha = 0.4,
    )
    savefig("figure\\hist_theta_vap.png")
    # 统计SB内外，磁场与R方向的夹角
    θva_1 = filter(!isnan,θva[sbpoints1])
    θva_2 = filter(!isnan,θva[sbpoints])
    θva_3 = filter(!isnan,θva[sbpoints2])
    histogram(
    θva_1*180/π,
    normalize = :pdf,
    label = "1",
    xlabel = "θ(B,r) °",
    alpha = 0.4,
    )
    histogram!(
    θva_2*180/π,
    normalize = :pdf,
    label = "2",
    alpha = 0.4,
    )
    histogram!(
    θva_3*180/π,
    normalize = :pdf,
    label = "3",
    alpha = 0.4,
    )
    savefig("figure\\hist_theta_va.png")

    # SB和全体太阳风，漂移速度与R方向的夹角+磁场与R方向的夹角
    sbpoints = Int.(sbpoints)
    maybesbpoints = Int.(maybesbpoints)
    θvαpSb = filter(!isnan,θvαp[sbpoints])
    θvaSb = filter(!isnan,θva[sbpoints])

    θvαpSw = filter(!isnan,θvαp)
    θvaSw = filter(!isnan,θva)

    histogram(
    θvαpSb*180/π,
    normalize = :pdf,
    label = "sb",
    xlabel = "θ(Vαp,r) °",
    alpha = 0.4,
    )
    histogram!(
    θvαpSw*180/π,
    normalize = :pdf,
    label = "sw",
    alpha = 0.4,
    )
    savefig("figure\\hist_theta_vap_sbsw.png")

    histogram(
    θvaSb*180/π,
    normalize = :pdf,
    label = "sb",
    xlabel = "θ(B,r) °",
    alpha = 0.4,
    )
    histogram!(
    θvaSw*180/π,
    normalize = :pdf,
    label = "sw",
    alpha = 0.4,
    )
    savefig("figure\\hist_theta_va_sbsw.png")
    nothing
end

pVars,αVars,modifiedVars,modifiedVars_va = loadData()
sbList = loadList()
sbEpochList = sbList["switchback_time_output"]
sbList = epoch2datetime.(sbEpochList)
sbNum = size(sbList)[1]
vα = [αVars["alpha_vel_rtn_sun"][i,:] for i in 1:size(αVars["alpha_vel_rtn_sun"])[1]]
vp = [modifiedVars["p_vel_rtn_sun_alphaEpoch"][i,:] for i in 1:size(modifiedVars["p_vel_rtn_sun_alphaEpoch"])[1]]
vαp = Vαp.(vα,vp)
vαp2va = vαp ./ modifiedVars["va_alphaEpoch"]

@load "list\\sbEventInfo.jld2" sbEventInfos
vαp_rtn = vα.-vp
va_rtn = [modifiedVars["va_rtn_alphaEpoch"][i,:] for i in 1:size(modifiedVars["va_rtn_alphaEpoch"])[1]]
θvαp_va = spanAngle.(vαp_rtn,va_rtn)
er_rtn = [[1,0,0],]
θva = spanAngle.(er_rtn,va_rtn)
θvαp = spanAngle.(er_rtn,vαp_rtn)

statSbSw(αVars,vαp,vαp2va,θvαp_va,θva,θvαp,sbEpochList)


#### 看看所有数据的时间序列
# t1 = DateTime(2020,5,1)
# αTimeLst = epoch2datetime.(αVars["alpha_epoch"])
# pTimeLst = epoch2datetime.(pVars["p_epoch"])
# p1 = scatter(
# αTimeLst,
# αVars["alpha_vel_rtn_sun"],
# ms=1,
# xlims=(t1,αTimeLst[end]),
# legend=false,
# ylabel="Vα Km/s",
# )
# p2 = scatter(
# pTimeLst,
# pVars["p_vel_rtn_sun"],
# ms=1,
# xlims=(t1,αTimeLst[end]),
# legend=false,
# ylabel="Vp Km/s",
# )
# p3 = plot(
# ylabel="sb",
# legend=false,
# ylims=(0,1),
# xlims=(t1,αTimeLst[end]),
# )

sbVecInfos = []
θδb_B0s = []
θδb_rs = []
θδb_δαps = []
θδb_vαp0s = []
δαp2vαps = []

sbidx = 187
magVars = matread("data\\psp_fld_mag_rtn_2020b.mat")
tend = DateTime(2021,9,1)
while sbEpochList[sbidx,2]<magVars["mag_epoch"][end]
# while sbidx<200
    println("sbidx=",sbidx)
    # if (sbEventInfos[sbidx]["sbtype"] != 1) |
    #    (sbEventInfos[sbidx]["αdataFlag"] == 0) |
    #    (sbEventInfos[sbidx]["ccpb"] < 0.6) |
    #    (sbEventInfos[sbidx]["ccαb"] < 0.6)
    #    global sbidx += 1
    #     continue
    # end
    if sbidx in notSB
        figDir = "notsb"
    elseif sbidx in maybenotSB
        figDir = "maybesb"
    else
        figDir = "issb"
    end
    # sbEvent_timeplot(
    # sbEpochList[sbidx,1],
    # sbEpochList[sbidx,2],
    # pVars,
    # αVars,
    # modifiedVars,
    # magVars,
    # vαp,
    # θvαp_va;
    # figName="SBevent"*string(sbidx),
    # figDir = figDir,
    # )
    # output = sbEvent(
    # sbEpochList[sbidx,1],
    # sbEpochList[sbidx,2],
    # pVars,
    # αVars,
    # modifiedVars,
    # magVars;
    # figName="SBevent"*string(sbidx),
    # plotTimeSeries=true,
    # figDir = figDir,
    # )
    output = sbEvent_calVectors(
    sbEpochList[sbidx,1],
    sbEpochList[sbidx,2],
    pVars,
    αVars,
    modifiedVars,
    magVars;
    figName="SBevent"*string(sbidx),
    plotTimeSeries=true,
    figDir = figDir,
    )
    output["sbidx"] = sbidx
    push!(sbVecInfos,output)
    if output["αdataFlag"]==1
        push!(θδb_B0s,output["θδb_B0"]*180/π)
        push!(θδb_rs,output["θδb_r"]*180/π)
        push!(θδb_δαps,output["θδb_δαp"]*180/π)
        push!(θδb_vαp0s,output["θδb_vαp0"]*180/π)
        push!(δαp2vαps,output["δαp2vαp"])
    end
    global sbidx += 1
end
magVars = matread("data\\psp_fld_mag_rtn_2021a.mat")
while sbEpochList[sbidx, 2] < magVars["mag_epoch"][end]
    println("sbidx=", sbidx)
    # if (sbEventInfos[sbidx]["sbtype"] != 1) |
    #    (sbEventInfos[sbidx]["αdataFlag"] == 0) |
    #    (sbEventInfos[sbidx]["ccpb"] < 0.6) |
    #    (sbEventInfos[sbidx]["ccαb"] < 0.6)
    #    global sbidx += 1
    #     continue
    # end
    if sbidx in notSB
        figDir = "notsb"
    elseif sbidx in maybenotSB
        figDir = "maybesb"
    else
        figDir = "issb"
    end
    # sbEvent_timeplot(
    # sbEpochList[sbidx,1],
    # sbEpochList[sbidx,2],
    # pVars,
    # αVars,
    # modifiedVars,
    # magVars,
    # vαp,
    # θvαp_va;
    # figName="SBevent"*string(sbidx),
    # figDir = figDir,
    # )

    # output = sbEvent(
    #     sbEpochList[sbidx, 1],
    #     sbEpochList[sbidx, 2],
    #     pVars,
    #     αVars,
    #     modifiedVars,
    #     magVars;
    #     figName = "SBevent" * string(sbidx),
    #     plotTimeSeries = true,
    #     figDir = figDir,
    # )
    # output["sbidx"] = sbidx
    # push!(sbEventInfos, output)
    output = sbEvent_calVectors(
    sbEpochList[sbidx,1],
    sbEpochList[sbidx,2],
    pVars,
    αVars,
    modifiedVars,
    magVars;
    figName="SBevent"*string(sbidx),
    plotTimeSeries=true,
    figDir = figDir,
    )
    output["sbidx"] = sbidx
    push!(sbVecInfos,output)
    if output["αdataFlag"]==1
        push!(θδb_B0s,output["θδb_B0"]*180/π)
        push!(θδb_rs,output["θδb_r"]*180/π)
        push!(θδb_δαps,output["θδb_δαp"]*180/π)
        push!(θδb_vαp0s,output["θδb_vαp0"]*180/π)
        push!(δαp2vαps,output["δαp2vαp"])
    end
    global sbidx += 1
end
magVars = matread("data\\psp_fld_mag_rtn_2021b.mat")
# while sbEpochList[sbidx,2]<magVars["mag_epoch"][end]
while sbEpochList[sbidx,2]<datetime2epoch(tend)
    println("sbidx=",sbidx)
    # if (sbEventInfos[sbidx]["sbtype"] != 1) |
    #    (sbEventInfos[sbidx]["αdataFlag"] == 0) |
    #    (sbEventInfos[sbidx]["ccpb"] < 0.6) |
    #    (sbEventInfos[sbidx]["ccαb"] < 0.6)
    #    global sbidx += 1
    #     continue
    # end
    if sbidx in notSB
        figDir = "notsb"
    elseif sbidx in maybenotSB
        figDir = "maybesb"
    else
        figDir = "issb"
    end
    # sbEvent_timeplot(
    # sbEpochList[sbidx,1],
    # sbEpochList[sbidx,2],
    # pVars,
    # αVars,
    # modifiedVars,
    # magVars,
    # vαp,
    # θvαp_va;
    # figName="SBevent"*string(sbidx),
    # figDir = figDir,
    # )

    # output = sbEvent(
    # sbEpochList[sbidx,1],
    # sbEpochList[sbidx,2],
    # pVars,
    # αVars,
    # modifiedVars,
    # magVars;
    # figName="SBevent"*string(sbidx),
    # plotTimeSeries=true,
    # figDir = figDir,
    # )
    # output["sbidx"] = sbidx
    # push!(sbEventInfos,output)
    output = sbEvent_calVectors(
    sbEpochList[sbidx,1],
    sbEpochList[sbidx,2],
    pVars,
    αVars,
    modifiedVars,
    magVars;
    figName="SBevent"*string(sbidx),
    plotTimeSeries=true,
    figDir = figDir,
    )
    output["sbidx"] = sbidx
    push!(sbVecInfos,output)
    if output["αdataFlag"]==1
        push!(θδb_B0s,output["θδb_B0"]*180/π)
        push!(θδb_rs,output["θδb_r"]*180/π)
        push!(θδb_δαps,output["θδb_δαp"]*180/π)
        push!(θδb_vαp0s,output["θδb_vαp0"]*180/π)
        push!(δαp2vαps,output["δαp2vαp"])
    end
    global sbidx += 1
end

# 统计每个sb里面矢量的特征
θB0_vαp0s = []
θB0_δvαps = []
θvαp0_rs = []
θB0_rs = []
for theinfo in sbVecInfos
    if theinfo["αdataFlag"] == 1

        push!(θB0_δvαps,
        spanAngle(theinfo["B0"],theinfo["e_δαp"])*180/π)

        if theinfo["δαp2vαp"]<1
            push!(θvαp0_rs,
            spanAngle(theinfo["vαp0"],[1,0,0])*180/π)
            push!(θB0_vαp0s,
            spanAngle(theinfo["B0"],theinfo["vαp0"])*180/π)
            push!(θB0_rs,
            spanAngle(theinfo["B0"],[1,0,0])*180/π)
        end
    end
end

histogram(
log10.(δαp2vαps),
xlabel = "log10(δVαp²/Vαp0²)",
normalize = :pdf,
legend = false,
# xlims = (0, 20),
)
savefig("figure\\hist_deltaVap2Vap.png")

histogram(
θδb_δαps,
xlabel = "θ",
label = "δb^δvαp",
normalize = :pdf,
alpha = 0.4,
xlims = (0,180),
# ylims = (0,0.025),
)
histogram!(
θδb_vαp0s,
xlabel = "θ",
label = "δb^vαp0",
normalize = :pdf,
alpha = 0.4,
)
savefig("figure\\hist_theta_deltab_and_others1.png")
histogram(
θδb_B0s,
xlabel = "θ",
label = "δb^B0",
normalize = :pdf,
alpha = 0.4,
xlims = (0,180),
)
histogram!(
θδb_rs,
xlabel = "θ",
label = "δb^r",
normalize = :pdf,
alpha = 0.4,
)
savefig("figure\\hist_theta_deltab_and_others2.png")
histogram(
θB0_vαp0s,
xlabel = "θ",
label = "vαp0^B0",
normalize = :pdf,
alpha = 0.4,
xlims = (0,180),
)
histogram!(
θB0_δvαps,
xlabel = "θ",
label = "δvαp^B0",
normalize = :pdf,
alpha = 0.4,
)
savefig("figure\\hist_theta_deltab_and_others3.png")

histogram(
θvαp0_rs,
xlabel = "θ",
label = "r^Vαp0",
normalize = :pdf,
xlims = (0,180),
alpha = 0.4,
)
histogram!(
θB0_rs,
xlabel = "θ",
label = "r^B0",
normalize = :pdf,
alpha = 0.4,
)
savefig("figure\\hist_theta_vap_and_r.png")

# for sbidx in 1:sbNum
# # for sbidx in 187:200
#
#
#
#     # @show output["αdataFlag"]
#     # plot!(p3,
#     # epoch2datetime.(sbEpochList[sbidx,:]),
#     # [0.5,0.5],
#     # )
# end
# plot(p1,p2,p3,layout=@layout [a;b;c])
# savefig("figure\\timeplot.png")

# @save "list\\sbEventInfo.jld2" sbEventInfos



# sbinfoList = zeros(77,11)
# listi = 1
#
# for (sbidx,sbInfo) in enumerate(sbEventInfos)
#     if (sbInfo["sbtype"] != 1) | (sbInfo["αdataFlag"]==0) | (sbInfo["ccpb"]<0.6) | (sbInfo["ccαb"]<0.6)
#         continue
#     end
#     sbinfoList[listi, 1] = sbidx
#     sbinfoList[listi, 2] = sbInfo["ccpb"]
#     sbinfoList[listi, 3] = sbInfo["ccαb"]
#     sbinfoList[listi, 4] = sbInfo["Vwp"]
#     sbinfoList[listi, 5] = sbInfo["rp"]
#     sbinfoList[listi, 6] = sbInfo["Vwα"]
#     sbinfoList[listi, 7] = sbInfo["rα"]
#     sbinfoList[listi, 8] = sbInfo["Cp"][1]
#     sbinfoList[listi, 9] = sbInfo["Cp"][2]
#     sbinfoList[listi, 10] = sbInfo["Cα"][1]
#     sbinfoList[listi, 11] = sbInfo["Cα"][2]
#
#     global listi += 1
# end
# # 筛选出比较符合SB特征、阿尔芬性较好、阿尔法粒子数据充足的事件
# XLSX.openxlsx("list\\sbevent\\sbevent_alfvenic.xlsx", mode="w") do xf
#     sheet = xf[1]
#     XLSX.rename!(sheet, "data")
#     sheet["A1"] = "idx"
#     sheet["B1"] = "CC_Vp&B"
#     sheet["C1"] = "CC_Vα&B"
#     sheet["D1"] = "Vwp"
#     sheet["E1"] = "rp"
#     sheet["F1"] = "Vwα"
#     sheet["G1"] = "rα"
#     sheet["H1"] = "Cxp"
#     sheet["I1"] = "Cyp"
#     sheet["J1"] = "Cxα"
#     sheet["K1"] = "Cyα"
#
#     sheet["A2"] = sbinfoList
# end
#
# scatter(
# abs.(sbinfoList[:,4]),
# abs.(sbinfoList[:,5]),
# xlabel = "Vwp Km/s",
# ylabel = "rp Km/s",
# legend = false,
# aspect_ratio = :equal,
# )
# plot!(
# x->x,
# ls=:dot,
# )
# savefig("figure\\sbevents\\alfvenic\\vwpVsRp.png")
#
# scatter(
# abs.(sbinfoList[:,6]),
# abs.(sbinfoList[:,7]),
# xlabel = "Vwα Km/s",
# ylabel = "rα Km/s",
# legend = false,
# aspect_ratio = :equal,
# )
# plot!(
# x->x,
# ls=:dot,
# )
# savefig("figure\\sbevents\\alfvenic\\vwαVsRα.png")
#
# markerSize=3
# scatter(
# sbinfoList[:,8],
# sbinfoList[:,10],
# xlabel = "Vp1 Km/s",
# ylabel = "Vα1 Km/s",
# legend = false,
# aspect_ratio = :equal,
# ms = markerSize,
# )
# plot!(
# x->x,
# ls=:dot,
# )
# savefig("figure\\sbevents\\alfvenic\\centre1.png")
# scatter(
# sbinfoList[:,9],
# sbinfoList[:,11],
# xlabel = "Vp2 Km/s",
# ylabel = "Vα2 Km/s",
# legend = false,
# aspect_ratio = :equal,
# ms = markerSize,
# )
# plot!(
# x->x,
# ls=:dot,
# )
# savefig("figure\\sbevents\\alfvenic\\centre2.png")
#
# CpCαDist = sqrt.(
# (sbinfoList[:,8].-sbinfoList[:,10]).^2 .+
# (sbinfoList[:,9].-sbinfoList[:,11]).^2
# )
# 记得改pVel用的哪个，是p的还是阿尔法的
