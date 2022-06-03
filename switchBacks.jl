
using MAT
using Dates
using Plots
using LinearAlgebra
using Statistics
using Optim
using Interpolations
using CurveFit

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

function circleShape(p::Vector)
    θ = LinRange(0,2π,500)
    p[1] .+ p[3]*sin.(θ), p[2] .+ p[3]*cos.(θ)
end

function sbEvent(epoch1,epoch2,pVars,αVars,modifiedVars,magVars;
    figName="test",plotTimeSeries=false,deltaMins=10)
    output = Dict(
    "αdataFlag"=>1,
    )
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


    ### 计算这段时间的磁场最小扰动方向（MVA）
    mag_rtn_scaled = mag_rtn .- mean(mag_rtn,dims=1)
    Σb = (mag_rtn_scaled'*mag_rtn_scaled)/length(magEpoch)
    # @show Σb
    F = svd(Σb)
    v1 = F.U[:,1]
    v2 = F.U[:,2]
    # v3 = F.U[:,3]
    # P = I - v3*v3' # P*x给出x在P的平面上的投影
    mag2P = [(v1'*mag_rtn')' (v2'*mag_rtn')']
    p0_mag = [0.,0.,10.]
    p_mag = get_circle(Matrix(mag2P'),p0_mag)
    pVel2P = [(v1'*pVel')' (v2'*pVel')']
    p0_p = [0.,0.,100.]
    p_p = get_circle(Matrix(pVel2P'),p0_p)
    ~,pMaxIdx = findmax(vec(
    sum(
    abs2,pVel2P.-
    # mean(pVel2P,dims=1),
    reshape(pVel2P[1,:],(1,2)),
    dims=2)))
    @show pMaxIdx
    αVel2P = [(v1'*αVel')' (v2'*αVel')']
    p0_α = [0.,0.,100.]
    p_α = get_circle(Matrix(αVel2P'),p0_α)
    ~,αMaxIdx = findmax(vec(
    sum(
    abs2,αVel2P.-
    reshape(αVel2P[1,:],(1,2)),
    # mean(αVel2P,dims=1),
    dims=2)))
    @show αMaxIdx
    makerSize = 3

    sign2sign = Dict(1=>"+",-1=>"-")
    avg_B = mean(sqrt.(sum(abs2,mag_rtn,dims=2)))
    mag1_interp_linear = LinearInterpolation(magEpoch, mag2P[:,1]/avg_B)
    mag1_pTime = mag1_interp_linear.(pEpoch)
    ccpb = cor(mag1_pTime, pVel2P[:,1])
    coefspb = linear_fit(mag1_pTime,pVel2P[:,1])# pVel1 = c + Vw*mag1_pTime
    plot_velVsB = scatter(
    mag1_pTime,
    pVel2P[:,1],
    label=:none,
    xlabel = "B1/|B|",
    ylabel = "Vp1 Km/s",
    title = "Vp1 = $(round(coefspb[1],digits=3))"*
    sign2sign[sign(coefspb[2])]*
    "$(round(abs(coefspb[2]),digits=3))*B1/|B|",
    )
    pb_func(x) = coefspb[1]+coefspb[2]*x
    plot!(
    plot_velVsB,
    pb_func,
    ls = :dot,
    color = :red,
    label = "cc="*string(round(ccpb,digits=3)),
    )

    savefig(plot_velVsB,"figure\\velVsB\\p_"*figName*".png")




    plot_vel = scatter(
    pVel2P[:,1],
    pVel2P[:,2],
    xlabel = "V1 Km/s",
    ylabel = "V2 Km/s",
    label = "Vp",
    color = :blue,
    ms = makerSize,
    )
    scatter!(
    plot_vel,
    [pVel2P[1,1],],
    [pVel2P[1,2],],
    color = :blue,
    markershape = :square,
    ms = makerSize+2,
    label=:none,
    )
    scatter!(
    plot_vel,
    [pVel2P[pMaxIdx,1],],
    [pVel2P[pMaxIdx,2],],
    color = :blue,
    markershape = :utriangle,
    ms = makerSize+2,
    label=:none,
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
    )
    scatter!(
    plot_vel,
    [αVel2P[1,1],],
    [αVel2P[1,2],],
    color = :red,
    markershape = :square,
    ms = makerSize+2,
    label=:none,
    )
    scatter!(
    plot_vel,
    [αVel2P[αMaxIdx,1],],
    [αVel2P[αMaxIdx,2],],
    color = :red,
    markershape = :utriangle,
    label=:none,
    ms = makerSize+2,
    aspect_ratio = :equal,
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
    # savefig("figure\\sbvel2P\\"*figName*".png")
    #################
    #
    plot_mag = scatter(
    mag2P[:,1],
    mag2P[:,2],
    xlabel = "B1 nT",
    ylabel = "B2 nT",
    ms = 1,
    color = :black,
    legend = false,
    aspect_ratio = :equal,
    )
    scatter!(
    plot_mag,
    [p_mag[1],],
    [p_mag[2],],
    color = :blue,
    markershape = :x,
    ms = makerSize+2,
    )
    plot!(
    plot_mag,
    circleShape(p_mag),
    lw= 0.5,
    linecolor = :black,
    linealpha = 0.5,
    )
    plot2PL = @layout grid(1,2)
    plot(
    plot_vel,
    plot_mag,
    layout=plot2PL,
    size=(800,450),
    )
    savefig("figure\\VelandB2P\\"*figName*".png")
    # savefig("figure\\B2P\\"*figName*".png")
    # F.U 给出svd分解得到的特征向量
    # F.S 给出svd分解得到的特征值，降序排列
    #@show F.U
    #@show F.S
    #@show F.V


    # 画出rtn坐标系中的矢量速度数据
    # makerSize = 5
    # scatter3d(
    # pVel[:,1],
    # pVel[:,2],
    # pVel[:,3],
    # xlabel = "Vr Km/s",
    # ylabel = "Vt Km/s",
    # zlabel = "Vn Km/s",
    # label = "Vp",
    # color = :blue,
    # ms = makerSize,
    # )
    # scatter3d!(
    # αVel[:,1],
    # αVel[:,2],
    # αVel[:,3],
    # label = "Vα",
    # color = :red,
    # ms = makerSize,
    # )
    # savefig("figure\\sbvel\\"*figName*".png")

    # 画出时间序列
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
    #     savefig("figure\\sbs\\"*figName*".png")
    # end
    output
end

pVars,αVars,modifiedVars,modifiedVars_va = loadData()
sbList = loadList()
sbEpochList = sbList["switchback_time_output"]
sbList = epoch2datetime.(sbEpochList)
sbNum = size(sbList)[1]

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
sbidx = 187
magVars = matread("data\\psp_fld_mag_rtn_2020b.mat")
tend = DateTime(2021,9,1)
while sbEpochList[sbidx,2]<magVars["mag_epoch"][end]
# while sbidx<199
    println("sbidx=",sbidx)
    output = sbEvent(
    sbEpochList[sbidx,1],
    sbEpochList[sbidx,2],
    pVars,
    αVars,
    modifiedVars,
    magVars;
    figName="SBevent"*string(sbidx),
    plotTimeSeries=true,
    )
    global sbidx += 1
end
magVars = matread("data\\psp_fld_mag_rtn_2021a.mat")
while sbEpochList[sbidx,2]<magVars["mag_epoch"][end]
    println("sbidx=",sbidx)
    output = sbEvent(
    sbEpochList[sbidx,1],
    sbEpochList[sbidx,2],
    pVars,
    αVars,
    modifiedVars,
    magVars;
    figName="SBevent"*string(sbidx),
    plotTimeSeries=true,
    )
    global sbidx += 1
end
magVars = matread("data\\psp_fld_mag_rtn_2021b.mat")
# while sbEpochList[sbidx,2]<magVars["mag_epoch"][end]
while sbEpochList[sbidx,2]<datetime2epoch(tend)
    println("sbidx=",sbidx)
    output = sbEvent(
    sbEpochList[sbidx,1],
    sbEpochList[sbidx,2],
    pVars,
    αVars,
    modifiedVars,
    magVars;
    figName="SBevent"*string(sbidx),
    plotTimeSeries=true,
    )
    global sbidx += 1
end


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
pPoints = vec((pVars["p_epoch"].>=sbEpochList[sbidx,1]) .& (pVars["p_epoch"].<=sbEpochList[sbidx,2]))
pEpoch = pVars["p_epoch"][pPoints]
pTime = epoch2datetime.(pEpoch)
pVel = pVars["p_vel_rtn_sun"][pPoints,:]
nonNaNpoints = vec(any(isnan.(pVel),dims=2) .!= 1)
pEpoch = pEpoch[nonNaNpoints]
pTime = pTime[nonNaNpoints]
pVel = pVel[nonNaNpoints,:]
