
# 变量单位：
# 速度：Km/s
# 到太阳距离：Km
# 温度：eV
# 温度tensor：eV
# 数密度：cm-3


using Statistics
using MAT
using Plots
using LinearAlgebra
using Dates

const AU = 149597871 # Km

################ utils ###############
function loadData()
    pVars = matread("data\\psp_spi_sf00.mat")
    αVars = matread("data\\psp_spi_sf0a.mat")
    modifiedVars = matread("data\\psp_modified.mat")
    pVars,αVars,modifiedVars
end

"""
计算p和α之间的速度差标量（α速度高为正）
"""
function Vαp(vα,vp)
    sign(norm(vα)-norm(vp))*norm(vα-vp)
end

"""
计算p和α之间速度差矢量（α-p）
可能还需要转换坐标系（目前是rtn那么似乎也不需要转了）
"""
function vecVαp(vα,vp)
    vα-vp
end

"""
epoch(matlab's datenum) to julia datetime
"""
function epoch2datetime(epoch::Real)
    unix2datetime((epoch-719529)*86400)
end

"""
Vαp VS R
"""
function vαpVsR(αVars::Dict,vp::Vector,vαp::Vector;binNum=5,maxR=0.4)
    R = vec(αVars["alpha_sun_dist"]/AU)
    binEdges = collect(range(minimum(R),maxR;length=(binNum+1)))
    Rs = zeros(binNum)
    vαps = zeros(binNum)
    for i in 1:binNum
        Rs[i] = (binEdges[i]+binEdges[i+1])/2
        vαps[i] = mean(filter(!isnan,vαp[(R.>binEdges[i]) .& (R.<binEdges[i+1])]))
    end
    Rs,vαps
end

pVars,αVars,modifiedVars = loadData()
# icmeList = matread("list\\PSP_ICME_list.mat")
sbList = matread("list\\switchback_time_output.mat")
sbEpochList = sbList["switchback_time_output"]
αTimeLst = epoch2datetime.(αVars["alpha_epoch"])
sbList = epoch2datetime.(sbEpochList)
vα = [αVars["alpha_vel_rtn_sun"][i,:] for i in 1:size(αVars["alpha_vel_rtn_sun"])[1]]
vp = [modifiedVars["p_vel_rtn_sun_alphaEpoch"][i,:] for i in 1:size(modifiedVars["p_vel_rtn_sun_alphaEpoch"])[1]]
vαp = Vαp.(vα,vp)

rs,vαps = vαpVsR(αVars,vp,vαp;binNum=8)


# scatter(
histogram2d(
αVars["alpha_sun_dist"]/AU,
vαp,
# xlabel="R [au]",
# ylabel="Vαp [Km/s]",
# ms=1,
legend=false,
ylims=(0,200),
xlims=(0,0.4),
)
plot!(
rs,
vαps,
xlabel="R [au]",
ylabel="Vαp [Km/s]",
legend=false,
color=:blue,
mark=:cross,
)
savefig("figure\\VαpVsR.png")


# savefig("figure\\VαpVsR.png")
