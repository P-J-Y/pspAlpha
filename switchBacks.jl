
using MAT
using Dates

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

function sbEvent(epoch1,epoch2,pVars,αVars,modifiedVars;figName="test")
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

    p1 = plot(
    pTime,
    pVel,
    label = ["r" "t" "n"],
    ylabel = "Vp Km/s",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    )
    p2 = plot(
    pTime,
    pTemp,
    legend=false,
    ylabel = "Tp eV",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    )
    p3 = plot(
    αTime,
    αVel,
    label = ["r" "t" "n"],
    ylabel = "Vα Km/s",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    )
    p4 = plot(
    αTime,
    αTemp,
    legend=false,
    ylabel = "Tα eV",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    )
    p5 = plot(
    αTime,
    va_rtn,
    label = ["r" "t" "n"],
    ylabel = "VA Km/s",
    xlims = (epoch2datetime(epoch1),epoch2datetime(epoch2)),
    )
    plot(p1,p2,p3,p4,p5,layout=@layout grid(5,1))
    savefig("figure\\sbs\\"*figName*".png")
    nothing
end

pVars,αVars,modifiedVars,modifiedVars_va = loadData()
sbList = loadList()
sbEpochList = sbList["switchback_time_output"]
sbList = epoch2datetime.(sbEpochList)
sbNum = size(sbList)[1]
for sbidx in 1:sbNum
    sbEvent(
    sbEpochList[sbidx,1],
    sbEpochList[sbidx,2],
    pVars,
    αVars,
    modifiedVars;
    figName="SBevent"*string(sbidx),
    )
end
