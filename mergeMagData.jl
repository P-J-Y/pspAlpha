# 把所有sb case以及附近1天内的磁场数据整合起来，构成一个更小的磁场数据文件
using MAT
using Dates
using JLD2




"""
epoch(matlab's datenum) to julia datetime
"""
function epoch2datetime(epoch::Real)
    unix2datetime((epoch-719529)*86400)
end

function datetime2epoch(adatetime::DateTime)
    datetime2unix(adatetime)/86400 + 719529
end

function loadList()
    sbList = matread("list//switchback_time_output.mat")
end
sbList = loadList()
sbEpochList = sbList["switchback_time_output"]

tend = DateTime(2021,9,1)


global sbidx = 187
magVars = matread("data\\psp_fld_mag_rtn_2020b.mat")
global magpoints = zeros(Bool,size(magVars["mag_epoch"]))
while sbEpochList[sbidx,2]<magVars["mag_epoch"][end]
    if ~(sbidx in isSB)
        global sbidx += 1
        continue
    end
# while sbidx<200

    epoch1 = sbEpochList[sbidx,1]-1/8
    epoch2 = sbEpochList[sbidx,2]+1/8
    thepoints = vec((magVars["mag_epoch"].>=epoch1) .&
                (magVars["mag_epoch"].<=epoch2))
    global magpoints = vec(magpoints .| thepoints)
    println("sbidx=",sbidx," points=",string(sum(magpoints)))
    global sbidx += 1
end
mag_epoch = magVars["mag_epoch"][magpoints]
mag_rtn = magVars["mag_rtn"][magpoints,:]

magVars = matread("data\\psp_fld_mag_rtn_2021a.mat")
global magpoints = zeros(Bool,size(magVars["mag_epoch"]))
while sbEpochList[sbidx,2]<magVars["mag_epoch"][end]
    if ~(sbidx in isSB)
        global sbidx += 1
        continue
    end

    epoch1 = sbEpochList[sbidx,1]-1/8
    epoch2 = sbEpochList[sbidx,2]+1/8
    thepoints = vec((magVars["mag_epoch"].>=epoch1) .&
                (magVars["mag_epoch"].<=epoch2))
    global magpoints = vec(magpoints .| thepoints)
    println("sbidx=",sbidx," points=",string(sum(magpoints)))
    global sbidx += 1
end
mag_epoch = [mag_epoch;magVars["mag_epoch"][magpoints]]
mag_rtn = [mag_rtn;magVars["mag_rtn"][magpoints,:]]

magVars = matread("data\\psp_fld_mag_rtn_2021b.mat")
global magpoints = zeros(Bool,size(magVars["mag_epoch"]))
while sbEpochList[sbidx,2]<datetime2epoch(tend)
    if ~(sbidx in isSB)
        global sbidx += 1
        continue
    end
    epoch1 = sbEpochList[sbidx,1]-1/8
    epoch2 = sbEpochList[sbidx,2]+1/8
    thepoints = vec((magVars["mag_epoch"].>=epoch1) .&
                (magVars["mag_epoch"].<=epoch2))
    global magpoints = vec(magpoints .| thepoints)
    println("sbidx=",sbidx," points=",string(sum(magpoints)))
    global sbidx += 1
end
mag_epoch = [mag_epoch;magVars["mag_epoch"][magpoints]]
mag_rtn = [mag_rtn;magVars["mag_rtn"][magpoints,:]]


file = matopen("data\\sbCaseMagData.mat", "w")
write(file, "mag_epoch", mag_epoch)
write(file, "mag_rtn", mag_rtn)
close(file)

magVars = Dict(
"mag_epoch"=>mag_epoch,
"mag_rtn"=>mag_rtn,
)
@save "data\\sbCaseMagData.jld2" magVars
