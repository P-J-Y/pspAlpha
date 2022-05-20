
using MAT



function loadList()
    sbList = matread("list//switchback_time_output.mat")
end


"""
epoch(matlab's datenum) to julia datetime
"""
function epoch2datetime(epoch::Real)
    unix2datetime((epoch-719529)*86400)
end

sbList = loadList()
sbEpochList = sbList["switchback_time_output"]
sbList = epoch2datetime.(sbEpochList)
