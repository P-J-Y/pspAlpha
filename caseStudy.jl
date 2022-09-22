
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

include("caseIdxs.jl")
# using CoordinateTransformations

# 原始磁场数据，每个文件对应的时间
# magFileTime = Dict(
# "2020b"=>[DateTime(2020,7,1),DateTime(2021,1,1)],
# "2021a"=>[DateTime(2021,1,1),DateTime(2021,7,1)],
# "2021b"=>[DateTime(2021,7,1),DateTime(2022,1,1)],
# )



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

"""
输入总数据和事件的时间，得到sb事件各种扰动之间的相对大小，用于列表筛选和分类事件
"""
function getCaseInfo(sbepoch1,sbepoch2,pVars,αVars,modifiedVars)
    output = Dict("αdataFlag"=>NaN,
    "δVα2VA"=>NaN,
    "t1"=>epoch2datetime(sbepoch1),
    "t2"=>epoch2datetime(sbepoch2)
    )
    pPoints = vec((pVars["p_epoch"].>=sbepoch1) .& (pVars["p_epoch"].<=sbepoch2))
    pEpoch = pVars["p_epoch"][pPoints]
    pVel = pVars["p_vel_rtn_sun"][pPoints,:]
    # pTemp = pVars["p_temp"][pPoints]
    αPoints =  vec((αVars["alpha_epoch"].>=sbepoch1) .&
                (αVars["alpha_epoch"].<=sbepoch2))
    αEpoch = αVars["alpha_epoch"][αPoints]
    αVel = αVars["alpha_vel_rtn_sun"][αPoints,:]
    if (length(αVel)-sum(isnan.(αVel))) < 10
        output["αdataFlag"] = 0
        return output
    end
    # αTemp = αVars["alpha_temp"][αPoints]
    va = modifiedVars["va_alphaEpoch"][αPoints]
    va_rtn = modifiedVars["va_rtn_alphaEpoch"][αPoints,:]
    vαp_rtn2α = αVel .- modifiedVars["p_vel_rtn_sun_alphaEpoch"][αPoints,:]

    nonNaNpoints = vec(any(isnan.(pVel),dims=2) .!= 1)
    pEpoch = pEpoch[nonNaNpoints]
    pVel = pVel[nonNaNpoints,:]
    nonNaNpoints = vec(any(isnan.(αVel),dims=2) .!= 1)
    αEpoch = αEpoch[nonNaNpoints]
    αVel = αVel[nonNaNpoints,:]
    vαp_rtn2α = vαp_rtn2α[nonNaNpoints,:]
    va = va[nonNaNpoints]
    va_rtn = va_rtn[nonNaNpoints,:]

    pTime = epoch2datetime.(pEpoch)
    vp_scaled = pVel.-mean(pVel,dims=1)
    αTime = epoch2datetime.(αEpoch)
    va0 = mean(va)
    vαp_rtn_scaled = vαp_rtn2α .- mean(vαp_rtn2α,dims=1)
    vα_scaled = αVel .- mean(αVel,dims=1)
    vαp0 = norm(mean(vαp_rtn2α,dims=1))
    drift2va = vαp0/va0
    δα2va = sqrt(sum(abs2,vα_scaled)/size(vα_scaled)[1])/va0
    δα2vas = sqrt.(sum(abs2,vα_scaled,dims=2))./va0
    δαp2va = sqrt(sum(abs2,vαp_rtn_scaled)/size(vαp_rtn_scaled)[1])/va0
    δp2va = sqrt(sum(abs2,vp_scaled)/size(vp_scaled)[1])/va0
    δαp2drift = δαp2va/drift2va

    output["δVα2VA"] = δα2va
    output["δVp2VA"] = δp2va
    output["δVαp2VA"] = δαp2va
    output["Vαp2VA"] = drift2va
    output["δVαp2Vαp"] = δαp2drift
    output
end

function getAllCaseInfo(pVars,αVars,modifiedVars,sbEpochList)
    caseIdx = []
    δVα2VA = []
    δVp2VA = []
    δVαp2VA = []
    Vαp2VA = []
    δVαp2Vαp = []
    t1 = []
    t2 = []
    magFiles = ["data\\psp_fld_mag_rtn_2020b.mat",
    "data\\psp_fld_mag_rtn_2021a.mat",
    "data\\psp_fld_mag_rtn_2021b.mat"]

    sbidx = 187
    tend = DateTime(2021,9,1)
    for magFile in magFiles
        magVars = matread(magFile)
        while sbEpochList[sbidx,2]<min(magVars["mag_epoch"][end],datetime2epoch(tend))
            println("sbidx=",sbidx)
            if sbidx in notSB
                sbidx += 1
                continue
            elseif sbidx in maybenotSB
                sbidx += 1
                continue
            end
            sbepoch1 = sbEpochList[sbidx,1]
            sbepoch2 = sbEpochList[sbidx,2]
            output = getCaseInfo(sbepoch1,sbepoch2,pVars,αVars,modifiedVars)
            if output["αdataFlag"] == 0
                sbidx += 1
                continue
            end
            push!(caseIdx,sbidx)
            push!(δVα2VA,output["δVα2VA"])
            push!(δVp2VA,output["δVp2VA"])
            push!(δVαp2VA,output["δVαp2VA"])
            push!(Vαp2VA,output["Vαp2VA"])
            push!(δVαp2Vαp,output["δVαp2Vαp"])
            push!(t1,output["t1"])
            push!(t2,output["t2"])
            sbidx += 1
        end
    end
    output = Dict(
    "caseIdx"=>caseIdx,
    "t1"=>t1,
    "t2"=>t2,
    "δVα2VA"=>δVα2VA,
    "δVp2VA"=>δVp2VA,
    "δVαp2VA"=>δVαp2VA,
    "Vαp2VA"=>Vαp2VA,
    "δVαp2Vαp"=>δVαp2Vαp,
    )
    return output
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


# magVars = matread("data\\psp_fld_mag_rtn_2020b.mat")
@load "data\\sbCaseMagData.jld2" magVars

# test function zone
# output = getAllCaseInfo(pVars,αVars,modifiedVars,sbEpochList)
# XLSX.openxlsx("data\\sbCaseInfo.xlsx", mode="w") do xf
#     sheet = xf[1]
#     XLSX.rename!(sheet, "sbInfo")
#     sheet["A1"] = "caseIdx"
#     sheet["B1"] = "t1"
#     sheet["C1"] = "t2"
#     sheet["D1"] = "δVα"
#     sheet["E1"] = "δVp"
#     sheet["F1"] = "δVαp"
#     sheet["G1"] = "Vαp"
#     sheet["H1"] = "δVαp/Vαp"
#
#     sheet["A2", dim=1] = output["caseIdx"]
#     sheet["B2", dim=1] = output["t1"]
#     sheet["C2", dim=1] = output["t2"]
#     sheet["D2", dim=1] = Float64.(output["δVα2VA"])
#     sheet["E2", dim=1] = Float64.(output["δVp2VA"])
#     sheet["F2", dim=1] = Float64.(output["δVαp2VA"])
#     sheet["G2", dim=1] = Float64.(output["Vαp2VA"])
#     sheet["H2", dim=1] = Float64.(output["δVαp2Vαp"])
# end
# scatter(output["δVp2VA"],output["δVα2VA"];
# xlims=(0,5),
# ylims=(0,5),
# )
#############################




# case study:
function caseStudy(sbidx,dirName,sbEpochList,pVars,αVars,modifiedVars,magVars)
    # for 1
    # sbepoch1 = sbEpochList[sbidx,1]-1/24
    # sbepoch2 = sbEpochList[sbidx,1]
    # for 2
    sbepoch1 = sbEpochList[sbidx,1]
    sbepoch2 = sbEpochList[sbidx,2]

    # output = getCaseInfo(sbepoch1,sbepoch2,pVars,αVars,modifiedVars)

    pPoints = vec((pVars["p_epoch"].>=sbepoch1) .& (pVars["p_epoch"].<=sbepoch2))
    pEpoch = pVars["p_epoch"][pPoints]
    pVel = pVars["p_vel_rtn_sun"][pPoints,:]
    # pTemp = pVars["p_temp"][pPoints]
    αPoints =  vec((αVars["alpha_epoch"].>=sbepoch1) .&
                (αVars["alpha_epoch"].<=sbepoch2))
    αEpoch = αVars["alpha_epoch"][αPoints]
    αVel = αVars["alpha_vel_rtn_sun"][αPoints,:]
    # αTemp = αVars["alpha_temp"][αPoints]
    # vαp_event = vαp[αPoints]
    va = modifiedVars["va_alphaEpoch"][αPoints]
    va_rtn = modifiedVars["va_rtn_alphaEpoch"][αPoints,:]
    vαp_rtn2α = αVel .- modifiedVars["p_vel_rtn_sun_alphaEpoch"][αPoints,:]
    magPoints = vec((magVars["mag_epoch"].>=sbepoch1) .&
                (magVars["mag_epoch"].<=sbepoch2))
    magEpoch = magVars["mag_epoch"][magPoints]
    mag_rtn = magVars["mag_rtn"][magPoints,:]

    # θvαp_va_event = θvαp_va[αPoints]

    nonNaNpoints = vec(any(isnan.(pVel),dims=2) .!= 1)
    pEpoch = pEpoch[nonNaNpoints]
    pVel = pVel[nonNaNpoints,:]
    pTime = epoch2datetime.(pEpoch)
    vp_scaled = pVel.-mean(pVel,dims=1)
    nonNaNpoints = vec(any(isnan.(αVel),dims=2) .!= 1)
    αEpoch = αEpoch[nonNaNpoints]
    αVel = αVel[nonNaNpoints,:]
    vαp_rtn2α = vαp_rtn2α[nonNaNpoints,:]
    va = va[nonNaNpoints]
    va_rtn = va_rtn[nonNaNpoints,:]
    nonNaNpoints = vec(any(isnan.(mag_rtn),dims=2) .!= 1)
    mag_rtn = mag_rtn[nonNaNpoints,:]
    magEpoch = magEpoch[nonNaNpoints]

    magTime = epoch2datetime.(magEpoch)
    theta = atan.(mag_rtn[:,2]./mag_rtn[:,1])
    αTime = epoch2datetime.(αEpoch)
    va0 = mean(va)
    vαp_rtn_scaled = vαp_rtn2α .- mean(vαp_rtn2α,dims=1)
    vα_scaled = αVel .- mean(αVel,dims=1)
    vαp0 = norm(mean(vαp_rtn2α,dims=1))
    drift2va = vαp0/va0
    δα2va = sqrt(sum(abs2,vα_scaled)/size(vα_scaled)[1])/va0
    δα2vas = sqrt.(sum(abs2,vα_scaled,dims=2))./va0
    δαp2va = sqrt(sum(abs2,vαp_rtn_scaled)/size(vαp_rtn_scaled)[1])/va0
    δp2va = sqrt(sum(abs2,vp_scaled)/size(vp_scaled)[1])/va0
    δαp2drift = δαp2va/drift2va


    # 统计几个角度的分布直方图
    # 是看磁场方向还是磁场扰动的方向？
    vecb = mat2vec(mag_rtn)
    pb = cartesian2Polar.(vecb)
    θb = [pb[i][2] for i in 1:length(pb)]
    vecVa = mat2vec(va_rtn)
    pVa = cartesian2Polar.(vecVa)
    θVa = [pVa[i][2] for i in 1:length(pVa)]
    rVa = [pVa[i][1]/va0 for i in 1:length(pVa)]
    # θb = spanAngle.(er_rtn,vecb)
    # vecVα = mat2vec(vα_scaled)
    vecVα = mat2vec(αVel)
    pVα = cartesian2Polar.(vecVα)
    θVα = [pVα[i][2] for i in 1:length(pVα)]
    rVα = [pVα[i][1]/va0 for i in 1:length(pVα)]
    # θvα = spanAngle.(er_rtn,vecVα)
    # vecVp = mat2vec(vp_scaled)
    vecVp = mat2vec(pVel)
    pVp = cartesian2Polar.(vecVp)
    θVp = [pVp[i][2] for i in 1:length(pVp)]
    rVp = [pVp[i][1]/va0 for i in 1:length(pVp)]
    # θvp = spanAngle.(er_rtn,vecVp)
    vecVαp_event = mat2vec(vαp_rtn2α)
    pVαp = cartesian2Polar.(vecVαp_event)
    θVαp_event = [pVαp[i][2] for i in 1:length(pVαp)]
    rVαp = [pVαp[i][1]/va0 for i in 1:length(pVαp)]

    # θvαp_event =  θvαp[αPoints]


    δθ = 0.2π
    θs = -π+δθ/2:δθ:π-δθ/2
    hb = fitHist(θb;θedges=-π:δθ:π)
    hVαp_event = fitHist(θVαp_event;θedges=-π:δθ:π)
    # p1 = plot(θs, hb.weights/sum(hb.weights);
    # proj = :polar, m = 2,
    # label="b",
    # title="Vαp0/VA0=$(round(drift2va,digits=3)), "*
    # "δVαp/VA=$(round(δαp2va,digits=3))",
    # markerstrokewidth=0,
    # )
    # plot!(p1,θs, hVαp_event.weights/sum(hVαp_event.weights);
    # proj = :polar, m = 2,
    # label="Vαp",
    # markerstrokewidth=0,
    # )
    hVα = fitHist(θVα;θedges=-π:δθ:π)
    hVp = fitHist(θVp;θedges=-π:δθ:π)
    p2 = plot(θs, hVα.weights/sum(hVα.weights);
    proj = :polar, m = 2,
    # label="δVα",
    label="Vα",
    # title="δVα/VA0=$(round(δα2va,digits=3))",
    markerstrokewidth=0,
    )
    plot!(p2,θs, hVp.weights/sum(hVp.weights);
    proj = :polar, m = 2,
    label="Vp",
    # label="δVp",
    markerstrokewidth=0,
    )
    # p3 = scatter(
    # θVa,rVa;
    # proj = :polar, m = 2,
    # label="VA",
    # title="VA0",
    # markerstrokewidth=0,
    # )
    # scatter!(p3,θVαp_event,rVαp;
    # proj = :polar, m = 2,
    # label="Vαp",
    # markerstrokewidth=0,
    # )
    p4 = scatter(
    θVα,rVα;
    proj = :polar, m = 2,
    # label="δVα",
    label="Vα",
    # title="VA0",
    markerstrokewidth=0,
    )
    scatter!(p4,θVp,rVp;
    proj = :polar, m = 2,
    # label="δVp",
    label="Vp",
    markerstrokewidth=0,
    )
    plot(p2,p4;
    layout=@layout grid(1,2))
    savefig("figure\\deltaVapAndVap0\\"*dirName*"\\case"*string(sbidx)*".png")
    nothing
end

# δVap/Vap0<0.4


dirName = "vap0\\VaVp"
for sbidx in vap0Cases
    caseStudy(sbidx,dirName,sbEpochList,pVars,αVars,modifiedVars,magVars)
end
dirName = "deltaVap\\VaVp"
for sbidx in deltaVapCases
    caseStudy(sbidx,dirName,sbEpochList,pVars,αVars,modifiedVars,magVars)
end
dirName = "mid\\VaVp"
for sbidx in midVapCases
    caseStudy(sbidx,dirName,sbEpochList,pVars,αVars,modifiedVars,magVars)
end




# 事件
# 磁场方向、速度方向，看场向漂移or非场向漂移
# VA与Vap0之间相对大小属于哪一类
# Vap0与δVap之间的相对大小属于哪一类
# 分析VA与Vap0接近时，是否就是阿尔法扰动很小的时候
