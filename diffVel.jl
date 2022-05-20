# va为什么大部分都在70附近？很不正常


# 变量单位：
# 速度：Km/s
# 到太阳距离：Km
# 温度：eV
# 温度tensor：eV
# 数密度：cm-3

using Distributions
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
    modifiedVars = matread("data\\psp_modified_alpha.mat")
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
function vαpVsR(αVars::Dict,vp::Vector,vαp::Vector;binNum=5,maxR=0.35)
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

"""
angle between 2 vectors
"""
function spanAngle(v1::Vector,v2::Vector)
    acos(dot(v1,v2)/(norm(v1)*norm(v2)))
end

pVars,αVars,modifiedVars = loadData()
# icmeList = matread("list\\PSP_ICME_list.mat")
αTimeLst = epoch2datetime.(αVars["alpha_epoch"])
vα = [αVars["alpha_vel_rtn_sun"][i,:] for i in 1:size(αVars["alpha_vel_rtn_sun"])[1]]
vp = [modifiedVars["p_vel_rtn_sun_alphaEpoch"][i,:] for i in 1:size(modifiedVars["p_vel_rtn_sun_alphaEpoch"])[1]]
va = vec(modifiedVars["va_alphaEpoch"])
va_rtn = [modifiedVars["va_rtn_alphaEpoch"][i,:] for i in 1:size(modifiedVars["va_rtn_alphaEpoch"])[1]]
vαp_rtn = vecVαp.(vα,vp)
vαp = Vαp.(vα,vp)
vαp2va = vαp./va

rs,absVαps = vαpVsR(αVars,vp,abs.(vαp);binNum=5)
~,vas = vαpVsR(αVars,vp,va;binNum=5)
~,vαp2vas = vαpVsR(αVars,vp,vαp2va;binNum=5)

saVαpVa = spanAngle.(vαp_rtn,va_rtn)

histogram(
saVαpVa*180/π,
legend=false,
xlabel="θ(Vαp,VA) °",
ylabel="counts",
)
savefig("figure\\hist_thetaDiffVA.png")

figRlims = (0.05,0.35)
histogram2d(
αVars["alpha_sun_dist"]/AU,
saVαpVa*180/π,
xlabel="R [au]",
ylabel="θ(Vαp,VA) °",
xlims=figRlims,
legend=false,
colorbar=true,
# colormap=:jet,
colorrange=(2000,20000),
)
savefig("figure\\thetaDiffVAVsR.png")

figRlims = (0.05,0.35)
# scatter(
p1 = histogram2d(
αVars["alpha_sun_dist"]/AU,
abs.(vαp),
# xlabel="R [au]",
# ylabel="Vαp [Km/s]",
# ms=1,
legend=false,
ylims=(0,200),
xlims=figRlims,
)
plot!(
p1,
rs,
absVαps,
# xlabel="R [au]",
ylabel="|Vαp| [Km/s]",
legend=false,
color=:red,
mark=:cross,
)
p2 = histogram2d(
αVars["alpha_sun_dist"]/AU,
va,
xlabel="R [au]",
ylabel="VA [Km/s]",
# ms=1,
legend=false,
ylims=(0,200),
xlims=figRlims,
)
plot!(
p2,
rs,
vas,
# xlabel="R [au]",
# ylabel="Vαp [Km/s]",
legend=false,
color=:red,
mark=:cross,
)
plot(p1,p2,layout=@layout grid(2,1))
savefig("figure\\absVαpandVaVsR.png")

p3 = histogram2d(
αVars["alpha_sun_dist"]/AU,
vαp2va,
xlabel="R [au]",
ylabel="Vαp/VA",
# ms=1,
legend=false,
)
plot!(
p3,
rs,
vαp2vas,
# xlabel="R [au]",
# ylabel="Vαp [Km/s]",
# legend=false,
label="avg Vαp/VA",
color=:red,
mark=:cross,
ylims=(-5.5,5.5),
xlims=figRlims,
)
plot!(
p3,
[0.,0.4],
[1,1],
legend=:topright,
color=:yellow,
label="Vαp=VA",
ls=:dash,
)
savefig("figure\\vαp2vaVsR.png")

# savefig("figure\\VαpVsR.png")
histogram(
vαp2va,
xlims=(-2,2),
xlabel="Vαp/VA",
ylabel="counts",
legend=false,
)
savefig("figure\\hist_vαpva.png")

histogram(
va,
#xlims=(-2,2),
xlabel="VA Km/s",
ylabel="counts",
legend=false,
)
savefig("figure\\hist_va.png")
