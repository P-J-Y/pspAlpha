
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
using StatsBase
# using CoordinateTransformations

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

function fitHist(x; θedges=-π:0.1π:π)
    result = StatsBase.fit(Histogram, x, θedges)
end

function mat2vec(V::Matrix)
    [V[i,:] for i in 1:size(V)[1]]
end

function cartesian2Polar(x::Vector;dim=2)
    x = x[1:dim]
    r = norm(x)
    θ = acos(x[1]/r)
    if x[2]<0
        θ = -θ
    end
    (r,θ)
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


magVars = matread("data\\psp_fld_mag_rtn_2020b.mat")

sbidx = 234
sbepoch1 = sbEpochList[sbidx,1]
sbepoch2 = sbEpochList[sbidx,2]

pPoints = vec((pVars["p_epoch"].>=sbepoch1) .& (pVars["p_epoch"].<=sbepoch2))
pEpoch = pVars["p_epoch"][pPoints]
pTime = epoch2datetime.(pEpoch)
pVel = pVars["p_vel_rtn_sun"][pPoints,:]
vp_scaled = pVel.-mean(pVel,dims=1)
# pTemp = pVars["p_temp"][pPoints]
αPoints =  vec((αVars["alpha_epoch"].>=sbepoch1) .&
            (αVars["alpha_epoch"].<=sbepoch2))
αEpoch = αVars["alpha_epoch"][αPoints]
αTime = epoch2datetime.(αEpoch)
αVel = αVars["alpha_vel_rtn_sun"][αPoints,:]
# αTemp = αVars["alpha_temp"][αPoints]
va = modifiedVars["va_alphaEpoch"][αPoints]
va0 = mean(va)
va_rtn = modifiedVars["va_rtn_alphaEpoch"][αPoints,:]
vαp_rtn2α = αVel .- modifiedVars["p_vel_rtn_sun_alphaEpoch"][αPoints,:]
vαp_rtn_scaled = vαp_rtn2α .- mean(vαp_rtn2α,dims=1)
vα_scaled = αVel .- mean(αVel,dims=1)
vαp_event = vαp[αPoints]
vαp0 = norm(mean(vαp_rtn2α,dims=1))
drift2va = vαp0/va0
δα2va = sqrt(sum(abs2,vα_scaled)/size(vα_scaled)[1])/va0
δα2vas = sqrt.(sum(abs2,vα_scaled,dims=2))./va0
δαp2va = sqrt(sum(abs2,vαp_rtn_scaled)/size(vαp_rtn_scaled)[1])/va0
δp2va = sqrt(sum(abs2,vp_scaled)/size(vp_scaled)[1])/va0
δαp2drift = δαp2va/drift2va
magPoints = vec((magVars["mag_epoch"].>=sbepoch1) .&
            (magVars["mag_epoch"].<=sbepoch2))
magEpoch = magVars["mag_epoch"][magPoints]
magTime = epoch2datetime.(magEpoch)
mag_rtn = magVars["mag_rtn"][magPoints,:]
theta = atan.(mag_rtn[:,2]./mag_rtn[:,1])
θvαp_va_event = θvαp_va[αPoints]

# 统计几个角度的分布直方图
# 是看磁场方向还是磁场扰动的方向？
vecb = mat2vec(mag_rtn)
pb = cartesian2Polar.(vecb)
θb = [pb[i][2] for i in 1:length(pb)]
# θb = spanAngle.(er_rtn,vecb)
vecVα = mat2vec(vα_scaled)
pVα = cartesian2Polar.(vecVα)
θVα = [pVα[i][2] for i in 1:length(pVα)]
# θvα = spanAngle.(er_rtn,vecVα)
vecVp = mat2vec(vp_scaled)
pVp = cartesian2Polar.(vecVp)
θVp = [pVp[i][2] for i in 1:length(pVp)]
# θvp = spanAngle.(er_rtn,vecVp)
vecVαp_event = mat2vec(vαp_rtn2α)
pVαp = cartesian2Polar.(vecVαp_event)
θVαp_event = [pVαp[i][2] for i in 1:length(pVαp)]

# θvαp_event =  θvαp[αPoints]


δθ = 0.2π
θs = -π+δθ/2:δθ:π-δθ/2
hb = fitHist(θb;θedges=-π:δθ:π)
hVαp_event = fitHist(θVαp_event;θedges=-π:δθ:π)
p1 = plot(θs, hb.weights/sum(hb.weights);
proj = :polar, m = 2,
label="b",
title="Vαp0/VA=$(round(drift2va,digits=3)), "*
"δVαp/VA=$(round(δαp2va,digits=3))")
plot!(p1,θs, hVαp_event.weights/sum(hVαp_event.weights);
proj = :polar, m = 2,
label="Vαp")
hVα = fitHist(θVα;θedges=-π:δθ:π)
hVp = fitHist(θVp;θedges=-π:δθ:π)
p2 = plot(θs, hVα.weights/sum(hVα.weights);
proj = :polar, m = 2,
label="δVα",
title="δVα/VA=$(round(δα2va,digits=3))",
)
plot!(p2,θs, hVp.weights/sum(hVp.weights);
proj = :polar, m = 2,
label="δVp")
plot(p1,p2;
layout=@layout grid(1,2))

# 事件
# 磁场方向、速度方向，看场向漂移or非场向漂移
# VA与Vap0之间相对大小属于哪一类
# Vap0与δVap之间的相对大小属于哪一类
# 分析VA与Vap0接近时，是否就是阿尔法扰动很小的时候
