using DataStructures, DataFrames, StatsFuns, GLM , JuMP,NLopt, StatStack, NLsolve, DBAPI, JavaCall, JDBC


# -------------------- Test regf_:mod -----------------------
#dfx[(dfx[:modelType].==:GLM)&(dfx[:model].==:occ) ,[:parameter,:coef,:stderr,:model,:modelType]]
deleterows!(df,find(isna(df[:,Symbol("B")])))

#f=readFixedFile("Occasion/Occ_fixed_effects.csv")
function readFixedFile(dfx::DataFrame,mname::String)
    vout = dfx[(dfx[:modelType].=="GLM")&(dfx[:model].==mname),[:parameter,:coef,:stderr,:zval,:pval]]
    names!(vout,[:x,:estimate, :std_error, :z_value,:pr_z_])
    return vout
end
readFixedFile(dfx,"occ")


#function rkey(key::String, ilevel::String="") ilevel != "" ? key*" ("*ilevel*")" : key  end
#function rkey(row::DataFrameRow) return rkey(row[:class], row[:level]) end

#x =readRandFile("Occasion/Occ_random_effects.csv",cfg)
function readRandFile(dfx::DataFrame,mname::String)
    vout = dfx[(dfx[:modelType].=="GLMM")&(dfx[:model].==mname),[:ranef,:parameter,:coef,:stderr,:zval,:pval]]
    vout[:adj_coef] = 0.0
    vout[:adj_stderr]=  0.0
    vout[:row] = 0
    vout[:levelse] = 0
    vout[:exposed_orig]= 1
    vout[:exposed] = true
    vout[:key] = map( (c,l) ->  c*" ("*l*")"  , vout[:ranef],vout[:parameter] ) 
    
    for row in eachrow(vout)
        n=vout[(vout[:ranef].==row[:ranef])&(vout[:parameter].=="none"),:coef][1]
        s=vout[(vout[:ranef].==row[:ranef])&(vout[:parameter].=="none"),:stderr][1]
        println(row[:ranef]," ~~ ",n," ~~ ",s)
        row[:adj_coef] = row[:coef]-n
        row[:adj_stderr] = sqrt(row[:stderr]^2+s^2)
        row[:coef] = 0
    end

    vout=vout[[  :ranef,:row,:parameter,:coef,:adj_coef,:levelse,:stderr,:adj_stderr,:exposed_orig,:exposed,:key ]]
    names!(vout,[:class,:row,:level,    :B0  ,:B1      ,:levelse,:SE0   ,:SE1       ,:exposed_orig,:exposed,:key ])

    return vout[vout[:level].!="none",:]
end


readRandFile(dfx,"occ")



#----- ok ----
using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib
#root="/mnt/resource/analytics/models/Jenny-o"
root="/mnt/resource/analytics/models/ALL#10"
jmod_fname = root*"/dfd_model.json"
mod_fname = root*"/dfd_model.csv"         
modelsDict = readModels(jmod_fname) 
modelsDict[:occ]=modelsDict[:iocc]; delete!(modelsDict,:iocc)
modelsDict[:dolocc]=modelsDict[:idolocc]; delete!(modelsDict,:idolocc)
modelsDict[:pen]=modelsDict[:ipen]; delete!(modelsDict,:ipen)
campaign_fname = root*"/campaign.csv" 
dfx = readtable(campaign_fname,header=true);
dfd = readtable(mod_fname,header=true);
dfd1 =dfd[(dfd[:buyer_pos_p1].==1),:];

dfx[:adj_coef] = 0.0
dfx[:adj_stderr]=  0.0
#dfx[:Equation] = ""
#dfx[:Equation0] = ""
#dfx[:Equation1] = ""




function regFixedGen(m::Dict, unAdjusted::Bool=false) 
    model=string(m[:modelName])
    modelType="GLM"
    dfname = m[:Buyer_Pos_P1_is1] ? "dfd1" : "dfd"
    adj0 = unAdjusted ?  "($dfname[:group].==0)," : ""
    adj1 = unAdjusted ?  "($dfname[:group].==1)," : ""
    println(dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model) ,:])
    intercept=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="(Intercept)") ,[:coef]][1][1]
    grp=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="group") ,[:coef]][1][1]
    vout0=string(intercept)*""
    vout1=string(intercept)*""
    results=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model) ,[:parameter,:coef]]
    for row in eachrow( results[findin(results[:parameter],map(x->string(x),    setdiff(m[:finalvars],[:group])  )) ,[:parameter,:coef]] )
        vout0=vout0*"+("*dfname*"[$adj0:"*string(row[:parameter])*"].*"*string(row[:coef])*")"
        vout1=vout1*"+("*dfname*"[$adj1:"*string(row[:parameter])*"].*"*string(row[:coef])*")"
    end
    vout1=vout1*"+"*string(grp)
    return vout0, vout1
end
r0,r1=regFixedGen(modelsDict[:occ],true)


viewGLM(mname::Symbol) = dfx[(dfx[:modelType].=="GLM" )& (dfx[:model].==string(mname)),:]
 
"""

hc_dolocc_fx1= dfx[(dfx[:modelType] .== :GLM )& (dfx[:model] .==  :dolocc)  ,:];
equ_dolocc=0
for i in 1: size(hc_dolocc_fx1,1)
    if hc_dolocc_fx1[i,1] == "(Intercept)" 
        equ_dolocc =string(hc_dolocc_fx1[i,2])
    elseif  hc_dolocc_fx1[i,1] == "group"  
        equ_dolocc =equ_dolocc
    else 
        equ_dolocc =string(equ_dolocc,+,"dfd1[:",parse(hc_dolocc_fx1[i,1]),"]",*,hc_dolocc_fx1[i,2])
    end
end
"""



"""
function regGEN()
    model="occ"
    modelType="GLM"
    dfname="dfd"
    println(dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model) ,:])
    intercept=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="(Intercept)") ,[:coef]][1][1]
    grp=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="group") ,[:coef]][1][1]
    vout=string(intercept)*""
    dfselect=""
    for row in eachrow(dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].!="(Intercept)")&(dfx[:parameter].!="group") ,:] )
        vout=vout*"+("*dfname*"[:"*string(row[:parameter])*"].*"*string(row[:coef])*")"
        dfselect=dfselect*",:"*string(row[:parameter])
       #println(row[:parameter]," ~~ ",row[:coef])
    end
    dfselect=dfname*"[["dfselect[2:end]*"]]"
    println(dfselect)
    vout=vout*"+"*string(grp)
    return vout
end
p=regGEN()

eval(parse(p))

p=parse("dfx[:ttt]=(dfx[:stderr].*3)+(dfx[:zval].*4)+(dfx[:pval].*5)")
"""



"""

function regFixedGen()   # Based On campaign.csv variables
    model="occ"
    modelType="GLM"
    dfname="dfd"
    println(dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model) ,:])
    intercept=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="(Intercept)") ,[:coef]][1][1]
    grp=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="group") ,[:coef]][1][1]
    vout=string(intercept)*""
    dfselect=""
    dfreg=""
    for row in eachrow(dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].!="(Intercept)")&(dfx[:parameter].!="group") ,:] )
        vout=vout*"+("*dfname*"[:"*string(row[:parameter])*"].*"*string(row[:coef])*")"
        p=string(row[:parameter])
        c=string(row[:coef])
        vmeta="\"(\"*string(row[:$p])*\"*$c)+\",\n"
        dfreg=dfreg*vmeta
        dfselect=dfselect*",:"*string(row[:parameter])
    end
    dfselect=dfname*"[1,["dfselect[2:end]*"]]"
    println("\n\nfor row in eachrow("*dfselect*")")
    #println("    println("*"\"y="*string(intercept)*"+\"*"*dfreg*"\""*string(grp)*"\")")   # WITH GROUP
    println("    println("*"\"y="*string(intercept)*"+\"*"*dfreg[1:end-4]*"\")")
    println("end\n\n")
    #vout=vout*"+"*string(grp)   # WITH GROUP
    return grp,vout
end
(grp,p)=regFixedGen()

#eval(parse(p))
occReg=DataFrame(fixed=eval(parse(p)))
occReg[:grp]=grp

"""


function regRandGen(dfd::DataFrame, modelsDict::Dict,modelIN::Symbol)   # Based On campaign.csv variables
    m=modelsDict[modelIN]
    model=string(m[:modelName])
    #imodel=Symbol("i"*replace(string(model),"i",""))
    modelType="GLMM"
    dfo=DataFrame()
    for r in m[:raneff]
        sr=string(r)
        dfo[r] = m[:Buyer_Pos_P1_is1] ? dfd[ dfd[:buyer_pos_p1].==1,r] : dfd[r]
        for row in eachrow(dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:ranef].==sr) ,:])
            dfo[ (dfo[r].==row[:parameter]),r] = string(row[:coef])
            println(r," ~~ ",row[:parameter])
        end
        dfo[r] = map(x->parse(Float64,x),dfo[r])
        dfo[r] = convert(Array{Float64},dfo[r])
    end
    return dfo
end
#x = regRandGen(modelsDict,:iocc)



function regFixedGen(dfd::DataFrame, modelsDict::Dict,modelIN::Symbol)   # Based On modelsDict
    m=modelsDict[modelIN]
    model=string(m[:modelName])
    #imodel=Symbol("i"*replace(string(model),"i",""))
    modelType="GLM"
    dfo=DataFrame()
    dfname="dfd"
    println(dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model) ,:])
    intercept=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="(Intercept)") ,[:coef]][1][1]
    grp=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:parameter].=="group") ,[:coef]][1][1]
    vout=string(intercept)*""
    dfselect=""
    dfreg=""
    results=dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model) ,[:parameter,:coef]]
    for row in eachrow( results[findin(results[:parameter],map(x->string(x),    setdiff(m[:finalvars],[:group])  )) ,[:parameter,:coef]] )
        vout=vout*"+("*dfname*"[:"*string(row[:parameter])*"].*"*string(row[:coef])*")"
        p=string(row[:parameter])
        c=string(row[:coef])
        vmeta="\"(\"*string(row[:$p])*\"*$c)+\",\n"
        dfreg=dfreg*vmeta
        dfselect=dfselect*",:"*string(row[:parameter])
    end
    dfselect=dfname*"[1,["dfselect[2:end]*"]]"
    println("\n\nfor row in eachrow("*dfselect*")")
    #println("    println("*"\"y="*string(intercept)*"+\"*"*dfreg*"\""*string(grp)*"\")")   # WITH GROUP
    println("    println("*"\"y="*string(intercept)*"+\"*"*dfreg[1:end-4]*"\")")
    println("end\n\n")
    #vout=vout*"+"*string(grp)   # WITH GROUP
    #dfd = m[:Buyer_Pos_P1_is1] ? dfd[ dfd[:buyer_pos_p1].==1, :] : dfd
    println(  "eval(parse("*vout*"))"   )
    dfo = DataFrame(panid=dfd[:panid], buyer_pos_p1=dfd[:buyer_pos_p1], group=dfd[:group], fixed=eval(parse(vout)))
    dfo = m[:Buyer_Pos_P1_is1] ? dfo[dfo[:buyer_pos_p1].==1,:] : dfo
    dfo[:grp] = grp
    #return grp,vout, dfo   ..... #(grp,p, dftf)=regFixedGen(dfd, modelsDict,:iocc)
    return dfo
end
#dff = regFixedGen(dfd, modelsDict,:iocc)

#modelsDict[:occ][:raneff]
rocc = hcat(regFixedGen(dfd, modelsDict,:occ), regRandGen(dfd, modelsDict,:occ)  )
#t = "rocc[:tot0]="*reduce(*, [".+rocc[:"*string(c)*"]" for c in setdiff(names(rocc),[:buyer_pos_p1,:panid,:grp])])[3:end]
t = "rocc[:tot0]="*reduce(*, [".+rocc[:"*string(c)*"]" for c in vcat([:fixed],modelsDict[:occ][:raneff])])[3:end]
eval(parse(t)) 
rocc[:tot1]=rocc[:tot0].+rocc[:grp]

rdolocc = hcat(regFixedGen(dfd, modelsDict,:dolocc), regRandGen(dfd, modelsDict,:dolocc)  )
t = "rdolocc[:tot0]="*reduce(*, [".+rdolocc[:"*string(c)*"]" for c in setdiff(names(rdolocc),[:buyer_pos_p1,:panid,:grp])])[3:end]
eval(parse(t)) 
rdolocc[:tot1]=rdolocc[:tot0].+rdolocc[:grp]

rpen = hcat(regFixedGen(dfd, modelsDict,:pen), regRandGen(dfd, modelsDict,:pen)  )
t = "rpen[:tot0]="*reduce(*, [".+rpen[:"*string(c)*"]" for c in setdiff(names(rpen),[:buyer_pos_p1,:panid,:grp])])[3:end]
eval(parse(t)) 
rpen[:tot1]=rpen[:tot0].+rpen[:grp]


# *****************  BUSY HERE ***********************
function adjustCoefs(dfx::DataFrame)
    for row in eachrow(dfx)
        if row[:modelType] == "GLMM"
            noneB(model::String,modelType::String,ranef::String) =  dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:ranef].==ranef),:coef][1]
            noneStdErr(model::String,modelType::String,ranef::String) =  dfx[(dfx[:modelType].==modelType)&(dfx[:model].==model)&(dfx[:ranef].==ranef),:stderr][1]
            println(row[:parameter],"... ",row[:coef]," ~~",noneB(row[:model],row[:modelType],row[:ranef]))
            row[:adj_coef] = row[:coef] - noneB(row[:model],row[:modelType],row[:ranef])
            row[:adj_stderr] = sqrt(row[:stderr]^2+ noneStdErr(row[:model],row[:modelType],row[:ranef])^2)
            ##NOTE check if the above actually writes to DF
        end
    end
    #dfx=dfx[  filter()  ,:]  # remove none's
    #dfx[dfx[].!="" ,:] 
end
adjustCoefs(dfx)







# **************************************************************************************************
# ***************************************** END BUSY ***********************************************
# **************************************************************************************************
# **************************************************************************************************
#Totals



# -----------------------------------------------------------
# ---- test -----
using DataStructures, DataArrays , DataFrames, StatsFuns, GLM , Distributions, MixedModels, StatsBase, JSON, StatLib
using DataStructures, DataFrames, StatsFuns, GLM , JuMP,NLopt, JStack, NLsolve
root="/mnt/resource/analytics/models/ALL#10"
jmod_fname = root*"/dfd_model.json"
mod_fname = root*"/dfd_model.csv"
modelsDict = readModels(jmod_fname) 
modelsDict[:occ]=modelsDict[:iocc]; delete!(modelsDict,:iocc)
modelsDict[:dolocc]=modelsDict[:idolocc]; delete!(modelsDict,:idolocc)
modelsDict[:pen]=modelsDict[:ipen]; delete!(modelsDict,:ipen)
campaign_fname = root*"/campaign.csv" 
dfx = readtable(campaign_fname,header=true);
#if isfile("campaign.csv")
#    dfx=readtable("campaign.csv",header=true);
#end

"""
-------- app.cfg --------
excludedBreaks = 
excludedKeys = 
P2_Competitor = true
pvalue_lvl = 0.20
exposed_flag_var = new_exposed_flag
sigLevel = 0.2
random_demos =
random_campaigns = 
dropvars = 
scoring_vars =
TotalModelsOnly=false
"""

cfgDefaults=OrderedDict( :P2_Competitor => true
                        ,:pvalue_lvl => 0.20  #pvalue_lvl = 0.20 
                        ,:excludedBreaks => AbstractString[]    #["estimated_hh_income","hh_age","number_of_children_in_living_un","person_1_gender"]
                        ,:excludedLevels => ["none"]
                        ,:excludedKeys => AbstractString[]
                        ,:exposed_flag_var => :exposed_flag
                        ,:sigLevel => "0.2"
                        ,:random_demos => Symbol[]
                        ,:random_campaigns => Symbol[]
                        ,:dropvars => Symbol[]
                        ,:scoring_vars => Symbol[]
                        ,:TotalModelsOnly=>false
                       )


cfg = JStack.loadCFG(cfgDefaults,pwd()*"/"*"app.cfg")
cfg[:counts] = isfile(pwd()*"/hhcounts.csv") ? true : false 
if !cfg[:counts] println("\n\n WARNING: HH Count file does not exist. No weighting will be applied to the scoring.\n\n\n") end


# ------------------------
"""
- Interested in the role
- I have a clear vision of what needs to be done and what is involved to get us where we need to be
- I don't want to be strictly ETL, and I mentioned that I ddn't want to leave the analytics side - and I have a potential solution that I'll mention in a second
- I don't mind being more focused 
- Having said all this,
"""
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------ SCORING -------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------
"""
function getFixedFormula(idf::DataFrame, logvar::AbstractString="", df_name::AbstractString="df_data")
    v_intercept= idf[idf[:x] .== "(Intercept)" , :estimate][1]
    v_group1 = idf[idf[:x] .== "group1" , :estimate][1]
    v_out = string(v_intercept)
    for fixedeff in setdiff(idf[:x], ["(Intercept)","group1"])
        v_coef=  idf[idf[:x] .== fixedeff, :estimate][1]
        v_out *= " + ("*string(v_coef)*"*"*df_name*"[:"*fixedeff*"])"
    end
    if length(logvar) > 0
        v_out*="+log( "*df_name*"[:"*logvar*"] + 1) "
    end
    v_out2=v_out*"+"*string(v_group1)
    return v_out, v_out2
end
idf = JStack.readFixedFile(dfx,"occ")
getFixedFormula(idf)
"""


#df_data = readtable("final_data.csv",header=true); lowercase(df_data)
df_data = readtable(mod_fname,header=true);  #lowercase(df_data)

mocc = MOcc(df_data,cfg,dfx)
mdolocc = MDolOcc(df_data,cfg,dfx)
mpen = MPen(df_data,cfg,dfx)
mdolhh=MDolHH(mocc,mdolocc,mpen)

cnts = Cnts(df_data, mocc,mdolocc,mpen,cfg[:TotalModelsOnly])

if cfg[:counts]
    hhcounts(cnts)
end

 
function aggregate!(df_data::DataFrame, focc::FOcc, cnts::Cnts)
    sdf = focc.sdf
    #sdf[:M] = cnts.get("Total Campaign")[:M][1]
    #sdf[:Mt] = cnts.get("Total Campaign")[:Mt][1]
    #sdf[:Mc] = cnts.get("Total Campaign")[:Mc][1]
    
    sdf[:M] = getk(cnts.sdf,"Total Campaign")[:M][1]   #cnts.get("Total Campaign")[:M][1]
    sdf[:Mt] = getk(cnts.sdf,"Total Campaign")[:Mt][1]   #cnts.get("Total Campaign")[:Mt][1]
    sdf[:Mc] = getk(cnts.sdf,"Total Campaign")[:Mc][1]   #cnts.get("Total Campaign")[:Mc][1]
    
    sdf[:mean_score0]=mean( df_data[ df_data[:buyer_pos_p1] .== 1 , :occ_score0] )
    sdf[:mean_score1]=mean( df_data[ df_data[:buyer_pos_p1] .== 1 , :occ_score1] ) 
    
    sdf[:unadj_mean_score0]=mean(df_data[(df_data[:buyer_pos_p1] .== 1) & (df_data[:group] .== 0),:occ_score0])
    sdf[:unadj_mean_score1]=mean(df_data[(df_data[:buyer_pos_p1] .== 1) & (df_data[:group] .== 1),:occ_score1])

    sdf[:adj_mean_cntrl_grp ] = mean( df_data[ df_data[:buyer_pos_p1] .== 1 , :occ_score0] )                     
    sdf[:adj_mean_expsd_grp ] = mean( df_data[ df_data[:buyer_pos_p1] .== 1 , :occ_score1] )
    # -------  Occ Score
    sdf[:unadj_avg_cntrl_hh_pre] = mean(df_data[ (df_data[:buyer_pre_p1] .== 1) & (df_data[:group] .== 0), :trps_pre_p1] )  
    sdf[:unadj_avg_expsd_hh_pre] = mean(df_data[ (df_data[:buyer_pre_p1] .== 1) & (df_data[:group] .== 1), :trps_pre_p1] )
    sdf[:unadj_avg_cntrl_hh_pst] = mean(df_data[ (df_data[:buyer_pos_p1] .== 1) & (df_data[:group] .== 0), :trps_pos_p1] )
    sdf[:unadj_avg_expsd_hh_pst] = mean(df_data[ (df_data[:buyer_pos_p1] .== 1) & (df_data[:group] .== 1), :trps_pos_p1] )
    
    #A hack to fix missing data -- 
    sdf[isnan(sdf[:unadj_avg_expsd_hh_pre]),:unadj_avg_expsd_hh_pre] = 0.0
    
end
#aggregate!(df_data, mocc.feff,cnts)



function aggregate!(df_data::DataFrame, fdolocc::FDolOcc, cnts::Cnts)
    sdf = fdolocc.sdf
    #sdf[:M] = cnts.get("Total Campaign")[:M][1]
    #sdf[:Mt] = cnts.get("Total Campaign")[:Mt][1]
    #sdf[:Mc] = cnts.get("Total Campaign")[:Mc][1]
    
    sdf[:M] = getk(cnts.sdf,"Total Campaign")[:M][1]   #cnts.get("Total Campaign")[:M][1]
    sdf[:Mt] = getk(cnts.sdf,"Total Campaign")[:Mt][1]   #cnts.get("Total Campaign")[:Mt][1]
    sdf[:Mc] = getk(cnts.sdf,"Total Campaign")[:Mc][1]   #cnts.get("Total Campaign")[:Mc][1]
    
    sdf[:mean_score0]=mean( df_data[ df_data[:buyer_pos_p1] .== 1 , :dolocc_score0] )
    sdf[:mean_score1]=mean( df_data[ df_data[:buyer_pos_p1] .== 1 ,  :dolocc_score1] )
    
    sdf[:unadj_mean_score0]=mean(df_data[(df_data[:buyer_pos_p1] .== 1) & (df_data[:group] .== 0),:dolocc_score0])
    sdf[:unadj_mean_score1]=mean(df_data[(df_data[:buyer_pos_p1] .== 1) & (df_data[:group] .== 1),:dolocc_score1])
    
    sdf[:adj_mean_cntrl_grp ] =  mean( df_data[ df_data[:buyer_pos_p1] .== 1 , :dolocc_score0] )             
    sdf[:adj_mean_expsd_grp ] =  mean( df_data[ df_data[:buyer_pos_p1] .== 1 ,  :dolocc_score1] )            
    # -------  DolOcc Score
    sdf[:unadj_avg_cntrl_hh_pre] =  mean(df_data[ (df_data[:buyer_pre_p1] .== 1) & (df_data[:group] .== 0), :dol_per_trip_pre_p1] )
    sdf[:unadj_avg_expsd_hh_pre] =  mean(df_data[ (df_data[:buyer_pre_p1] .== 1) & (df_data[:group] .== 1), :dol_per_trip_pre_p1] )
    sdf[:unadj_avg_cntrl_hh_pst] =  mean(df_data[ (df_data[:buyer_pos_p1] .== 1) & (df_data[:group] .== 0), :dol_per_trip_pos_p1] )
    sdf[:unadj_avg_expsd_hh_pst] =  mean(df_data[ (df_data[:buyer_pos_p1] .== 1) & (df_data[:group] .== 1), :dol_per_trip_pos_p1] )
    
    #A hack to fix missing data -- 
    sdf[isnan(sdf[:unadj_avg_expsd_hh_pre]),:unadj_avg_expsd_hh_pre] = 0.0
    
end
#aggregate!(df_data, mdolocc.feff, cnts)


function aggregate!(df_data::DataFrame, fpen::FPen, cnts::Cnts)
    sdf = fpen.sdf
    #sdf[:M] = cnts.get("Total Campaign")[:N][1]
    #sdf[:Mt] = cnts.get("Total Campaign")[:Nt][1]
    #sdf[:Mc] = cnts.get("Total Campaign")[:Nc][1]
    
    sdf[:M] = getk(cnts.sdf,"Total Campaign")[:N][1]   #cnts.get("Total Campaign")[:M][1]
    sdf[:Mt] = getk(cnts.sdf,"Total Campaign")[:Nt][1]   #cnts.get("Total Campaign")[:Mt][1]
    sdf[:Mc] = getk(cnts.sdf,"Total Campaign")[:Nc][1]   #cnts.get("Total Campaign")[:Mc][1]
    
    sdf[:mean_score0]=mean( df_data[ :pen_score0] )
    sdf[:mean_score1]=mean( df_data[ :pen_score1] )
    sdf[:unadj_mean_score0]=mean(df_data[(df_data[:group] .== 0),:pen_score0])
    sdf[:unadj_mean_score1]=mean(df_data[(df_data[:group] .== 1),:pen_score1])
    
    sdf[:adj_mean_cntrl_grp ] = mean( df_data[ :pen_score0] )             
    sdf[:adj_mean_expsd_grp ] = mean( df_data[ :pen_score1] )            
    # -------  Pen Score   - buyer_pre_p1
    sdf[:unadj_avg_cntrl_hh_pre ] = mean(df_data[ (df_data[:group] .== 0)  , :buyer_pre_p1] )
    sdf[:unadj_avg_expsd_hh_pre ] = mean(df_data[ (df_data[:group] .== 1)  , :buyer_pre_p1] )
    sdf[:unadj_avg_cntrl_hh_pst ] = mean(df_data[ (df_data[:group] .== 0)  , :buyer_pos_p1] )
    sdf[:unadj_avg_expsd_hh_pst ] = mean(df_data[ (df_data[:group] .== 1)  , :buyer_pos_p1] )
    
    #A hack to fix missing data -- 
    sdf[isnan(sdf[:unadj_avg_expsd_hh_pre]),:unadj_avg_expsd_hh_pre] = 0.0
    
end
#aggregate!(df_data, mpen.feff, cnts)


function aggregateCommon!(fixedeff::FixedEffect)
    sdf = fixedeff.sdf
    #src = fixedeff.src
    sdf[:adj_dod_effct] = ((sdf[:adj_mean_expsd_grp] - sdf[:adj_mean_cntrl_grp])  ./ sdf[:adj_mean_cntrl_grp] ) * 100   
    sdf[:twotail_pval] = sdf[:P]   #twotail_pval      
    sdf[:onetail_pval] = 1 - (sdf[:twotail_pval] ./ 2)   #onetail_pval
    sdf[:twotail_pval] = 1-sdf[:P]    #redo :twotail_pval
    #For NON-SIG Models, update - adj_mean_expsd_grp (to = adj_mean_cntrl_grp) & adj_dod_effct (= 0)
    #sdf[(sdf[:twotail_pval] .< 0.80) & (sdf[:twotail_pval] .!= 0),:adj_mean_expsd_grp]      =  sdf[(sdf[:twotail_pval] .< 0.80) & (sdf[:twotail_pval] .!= 0),:adj_mean_cntrl_grp]
end
#aggregateCommon!(mocc.feff)
#aggregateCommon!(mdolocc.feff)
#aggregateCommon!(mpen.feff) 

# -------- Confidence Intervals - Total --------
function ConfidenceIntervals(feff::FixedEffect)
    sdf = feff.sdf
    for row in eachrow(sdf)
            k = row[:key]
            md = JStack.df2dict(row)    
            md[:mean_score0] = md[:unadj_mean_score0];
            md[:mean_score1] = md[:unadj_mean_score1];
            cis =  CIs(md, ZDict, 0)
            for (zscore_key,zscore) in  ZDict     
                    ubk = Symbol(zscore_key*"_ub")
                    lbk = Symbol(zscore_key*"_lb")
                    sdf[sdf[:key].==k, lbk] = md[lbk]*100
                    sdf[sdf[:key].==k, ubk] = md[ubk]*100
            end
    end
end


aggregate!(df_data, mocc.feff,cnts)
aggregate!(df_data, mpen.feff, cnts)
aggregate!(df_data, mdolocc.feff, cnts)
aggregateCommon!(mocc.feff)
aggregateCommon!(mdolocc.feff)
aggregateCommon!(mpen.feff) 
ConfidenceIntervals(mocc.feff)
ConfidenceIntervals(mdolocc.feff)
ConfidenceIntervals(mpen.feff)


# --------------------------------------------------------
# ---------- BREAKS --------------------------------------
# --------------------------------------------------------

# === Random Effects ====
if !cfg[:TotalModelsOnly]
    
function init!(cnts::Cnts, m::MModel)  
        src=m.reff.src
        # --- GROUP1 -----
        src[:Z1] = src[:B1] ./ src[:SE1]
        src[:P1] = 2.0 * ccdf(Normal(), abs(src[:Z1]))
        src[:SIG1] = map( x -> ((isnan(x))|(x > 0.2)) ?  false : true, src[:P1] ) 
        # --- GROUP0 -----
        src[:Z0] = src[:B0] ./ src[:SE0]                          
        src[:P0] = 2.0 * ccdf(Normal(), abs(src[:Z0]))                       
        src[:SIG0] = map( x -> (isnan(x)|(x > 0.2)) ?  false : true, src[:P0] )  
        
        m.reff.sdf = join(m.reff.sdf,m.reff.src[[:class,:level,:B0,:B1,:SE0,:SE1,:exposed,:P1,:key]], on = :key)
        sdf=m.reff.sdf        
        sdf[:B_fixed] = m.feff.sdf[1,:B]
        sdf[:P_fixed] = m.feff.sdf[1,:P]
        sdf[:SE_fixed] = m.feff.sdf[1,:SE]
        #sdf[:B1_combo] = sdf[:B1] + sdf[:B_fixed] #- sdf[:B0]
        #sdf[:SE1_combo] = sqrt(sdf[:SE_fixed].^2+sdf[:SE1].^2)
        sdf[:B1_combo] = sdf[:B1] + sdf[:B_fixed] - sdf[:B0]
        sdf[:SE1_combo] = sqrt(sdf[:SE_fixed].^2+sdf[:SE1].^2+sdf[:SE0].^2)
        sdf[:P1_combo] = 2.0 * ccdf(Normal(), abs( (sdf[:B1_combo]) ./ sdf[:SE1_combo])) 
        
        if m.reff.v_model == "pen"
            for k in sdf[:key]     
                sdf[sdf[:key].==k ,:M] = getk(cnts.sdf,k)[:N][1] 
                sdf[sdf[:key].==k ,:Mt] = getk(cnts.sdf,k)[:Nt][1]
                sdf[sdf[:key].==k ,:Mc] = getk(cnts.sdf,k)[:Nc][1]
            end  
        else
            for k in sdf[:key]     
                sdf[sdf[:key].==k ,:M] = getk(cnts.sdf,k)[:M][1]
                sdf[sdf[:key].==k ,:Mt] = getk(cnts.sdf,k)[:Mt][1]
                sdf[sdf[:key].==k ,:Mc] = getk(cnts.sdf,k)[:Mc][1] 
            end   
        end 
end
init!(cnts, mocc)
init!(cnts, mdolocc)
init!(cnts, mpen)


function aggregate!(df_data::DataFrame, mocc::MOcc)
    sdf=mocc.reff.sdf
    for row = eachrow(mocc.reff.sdf)
        k = row[:key]
        ranfx = lowercase(row[:class])
        sranfx = Symbol(ranfx)
        v_level = row[:level]
        exposed = row[:exposed]  
        dest_colname0=Symbol(k*"_occ_0")
        dest_colname1=Symbol(k*"_occ_1")
        mean_score0 = mean(df_data[df_data[:buyer_pos_p1] .== 1, dest_colname0])
        mean_score1 = mean(df_data[ df_data[:buyer_pos_p1] .== 1 , dest_colname1])
        #println(mean_score1)
        sdf[sdf[:key].==k ,:mean_score0] = mean_score0
        sdf[sdf[:key].==k ,:mean_score1] = mean_score1
        sdf[sdf[:key].==k ,:adj_mean_cntrl_grp] = mean_score0
        sdf[sdf[:key].==k ,:adj_mean_expsd_grp] = mean_score1
        
        sdf[sdf[:key].==k ,:unadj_avg_expsd_hh_pre] = mean(df_data[ (df_data[:group] .== 1) & (df_data[:buyer_pre_p1] .== 1) & (df_data[sranfx] .== v_level) , :trps_pre_p1] )  #(df_data[:buyer_pre_p1] .== 1)
        sdf[sdf[:key].==k ,:unadj_avg_expsd_hh_pst] = mean(df_data[ (df_data[:group] .== 1) & (df_data[:buyer_pos_p1] .== 1) & (df_data[sranfx] .== v_level) , :trps_pos_p1] )  #(df_data[:buyer_pos_p1] .== 1)
        if exposed
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pre] = mocc.feff.sdf[1,:unadj_avg_cntrl_hh_pre]
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pst] = mocc.feff.sdf[1,:unadj_avg_cntrl_hh_pst]
        else
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pre] = mean(df_data[ (df_data[:group] .== 0) & (df_data[:buyer_pre_p1] .== 1) & (df_data[sranfx] .== v_level) , :trps_pre_p1] )  #(df_data[:buyer_pre_p1] .== 1)
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pst] = mean(df_data[ (df_data[:group] .== 0) & (df_data[:buyer_pos_p1] .== 1) & (df_data[sranfx] .== v_level) , :trps_pos_p1] )  #(df_data[:buyer_pos_p1] .== 1)
        end       
    end 
end
aggregate!(df_data, mocc)


function aggregate!(df_data::DataFrame, mdolocc::MDolOcc)
    sdf=mdolocc.reff.sdf
    for row = eachrow(mdolocc.reff.sdf)
        k = row[:key]
        ranfx = lowercase(row[:class])
        sranfx = Symbol(ranfx)
        v_level = row[:level]
        exposed = row[:exposed]         
        dest_colname0=Symbol(k*"_dolocc_0")
        dest_colname1=Symbol(k*"_dolocc_1")
        mean_score0 = mean(df_data[df_data[:buyer_pos_p1] .== 1, dest_colname0])
        mean_score1 = mean(df_data[ df_data[:buyer_pos_p1] .== 1 , dest_colname1])
        sdf[sdf[:key].==k ,:mean_score0] = mean_score0
        sdf[sdf[:key].==k ,:mean_score1] = mean_score1   
        sdf[sdf[:key].==k ,:adj_mean_cntrl_grp] = mean_score0
        sdf[sdf[:key].==k ,:adj_mean_expsd_grp] = mean_score1
        sdf[sdf[:key].==k ,:unadj_avg_expsd_hh_pre] = mean(df_data[ (df_data[:group] .== 1) & (df_data[:buyer_pre_p1] .== 1) & (df_data[sranfx] .== v_level) , :dol_per_trip_pre_p1] )
        sdf[sdf[:key].==k ,:unadj_avg_expsd_hh_pst] = mean(df_data[ (df_data[:group] .== 1) & (df_data[:buyer_pos_p1] .== 1) & (df_data[sranfx] .== v_level) , :dol_per_trip_pos_p1] )
        #sdf[sdf[:key].==k ,:unadj_avg_expsd_hh_pre] = mean(df_data[ (df_data[:group] .== 1) & (df_data[sranfx] .== v_level) , :dol_per_trip_pre_p1] )
        #sdf[sdf[:key].==k ,:unadj_avg_expsd_hh_pst] = mean(df_data[ (df_data[:group] .== 1) & (df_data[sranfx] .== v_level) , :dol_per_trip_pos_p1] )
        
        if exposed
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pre] = mdolocc.feff.sdf[1,:unadj_avg_cntrl_hh_pre]
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pst] = mdolocc.feff.sdf[1,:unadj_avg_cntrl_hh_pst]
        else
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pre] = mean(df_data[ (df_data[:group] .== 0) & (df_data[:buyer_pre_p1] .== 1) & (df_data[sranfx] .== v_level) , :dol_per_trip_pre_p1] ) 
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pst] = mean(df_data[ (df_data[:group] .== 0) & (df_data[:buyer_pos_p1] .== 1) & (df_data[sranfx] .== v_level) , :dol_per_trip_pos_p1] )
            #sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pre] = mean(df_data[ (df_data[:group] .== 0) & (df_data[sranfx] .== v_level) , :dol_per_trip_pre_p1] ) 
            #sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pst] = mean(df_data[ (df_data[:group] .== 0) & (df_data[sranfx] .== v_level) , :dol_per_trip_pos_p1] )
        end    
    end   
end
aggregate!(df_data, mdolocc)



function aggregate!(df_data::DataFrame,mpen::MPen)
    sdf=mpen.reff.sdf
    for row = eachrow(mpen.reff.sdf)
        k = row[:key]
        ranfx = lowercase(row[:class])
        sranfx = Symbol(ranfx)
        v_level = row[:level]
        exposed = row[:exposed]         
        dest_colname0=Symbol(k*"_pen_0")
        dest_colname1=Symbol(k*"_pen_1")    
        mean_score0 = mean(df_data[dest_colname0])
        mean_score1 = mean(df_data[dest_colname1])   #println("..... ",mean_score0," ~ ",mean_score1)
        sdf[sdf[:key].==k ,:mean_score0] = mean_score0
        sdf[sdf[:key].==k ,:mean_score1] = mean_score1 
        sdf[sdf[:key].==k ,:adj_mean_cntrl_grp] = mean_score0
        sdf[sdf[:key].==k ,:adj_mean_expsd_grp] = mean_score1
        sdf[sdf[:key].==k ,:unadj_avg_expsd_hh_pre] = mean(df_data[ (df_data[:group] .== 1) & (df_data[Symbol(ranfx)] .== v_level) , :buyer_pre_p1] )
        sdf[sdf[:key].==k ,:unadj_avg_expsd_hh_pst] = mean(df_data[ (df_data[:group] .== 1) & (df_data[Symbol(ranfx)] .== v_level) , :buyer_pos_p1] )
        if exposed
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pre] = mpen.feff.sdf[1,:unadj_avg_cntrl_hh_pre]
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pst] = mpen.feff.sdf[1,:unadj_avg_cntrl_hh_pst]    
        else
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pre] = mean(df_data[ (df_data[:group] .== 0) & (df_data[Symbol(ranfx)] .== v_level) , :buyer_pre_p1] )
            sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pst] = mean(df_data[ (df_data[:group] .== 0) & (df_data[Symbol(ranfx)] .== v_level) , :buyer_pos_p1] )
        end         
    end      
end
aggregate!(df_data, mpen)
     

    #sdf[:B1_combo] = sdf[:B1] + sdf[:B_fixed] - sdf[:B0]
    #sdf[:SE1_combo] = sqrt(sdf[:SE_fixed].^2+sdf[:SE1].^2+sdf[:SE0].^2)
    #sdf[:P1_combo] = 2.0 * ccdf(Normal(), abs( (sdf[:B1_combo]) ./ sdf[:SE1_combo])) 

function aggregateCommon!(raneff::RanEffect,cfg::OrderedDict)
    sdf = raneff.sdf
    sdf[:adj_dod_effct] = ((sdf[:adj_mean_expsd_grp] - sdf[:adj_mean_cntrl_grp])  ./ sdf[:adj_mean_cntrl_grp] ) * 100   
    sdf[:twotail_pval] = (1-sdf[:P1_combo]) * 100
    sdf[:onetail_pval] = (1-(sdf[:P1_combo] ./ 2)) * 100

    sdf[:Braw] = sdf[:B1] - sdf[:B0]
    sdf[:SEraw] = sqrt(sdf[:SE1].^2+sdf[:SE0].^2)
    sdf[:Zraw] = sdf[:Braw] ./ sdf[:SEraw]
    sdf[:Praw] = 2.0 * ccdf(Normal(), abs(sdf[:Zraw]))   
    sdf[:twotail_pval_raw] = (1-sdf[:Praw]) * 100
    sdf[:onetail_pval_raw] = (1-(sdf[:Praw] ./ 2)) * 100
    
    sigLvl=parse(Float64,get(cfg, :sigLevel, "0.2"))
    
    
end
aggregateCommon!(mocc.reff,cfg)
aggregateCommon!(mdolocc.reff,cfg)
aggregateCommon!(mpen.reff,cfg) 

#--add 3 random coefs..... combin STD_ERR
# for non-significant - set to 0... or default to total....


function np_val(iDict::OrderedDict)
    M = get(iDict, :M, NA)
    Mt = get(iDict, :Mt, NA)
    Mc = get(iDict, :Mc, NA)
    mean_score0 = get(iDict, :adj_mean_cntrl_grp, NA)
    mean_score1 = get(iDict, :adj_mean_expsd_grp, NA)
    B = get(iDict, :B1_combo, NA)
    SE = get(iDict, :SE1_combo, NA)
    #zscore = get(iDict, :zscore, NA)
    #println("NLsolve : ",iDict,"\n")
    println(get(iDict,:key, NA),"~",M,"~",Mt,"~",Mc,"~",mean_score0,"~",mean_score1,"~",B,"~",SE)
    
    function f!(x, fvec)
        #fvec[1] = exp(x[1])-45
        Lb_pre = (   (mean_score1*(Mt/M))    -   (mean_score1*exp(-(B-(SE*x[1])))*(Mt/M))   )  +
                 (   (mean_score0*exp((B-(SE*x[1])))*(Mc/M))    -   (mean_score0*(Mc/M))    )
        Lb = Lb_pre/mean_score0
        fvec[1] = Lb
        ### ------------ Lower Bound ---------------
        #Lb_pre = (   (mean_score1*(Mt/M))    -   (mean_score1*exp(-(B-(SE*zscore)))*(Mt/M))   )  +
        #         (   (mean_score0*exp((B-(SE*zscore)))*(Mc/M))    -   (mean_score0*(Mc/M))    )
        #Lb = Lb_pre/mean_score0
        ### ------------ Upper Bound ---------------
        #Ub_pre =  (     ( mean_score1*(Mt/M) )   -   ( mean_score1*exp(-(B+(SE*zscore)))*(Mt/M))   )  + 
        #          (     ( mean_score0*exp((B+(SE*zscore)))*(Mc/M))  - (mean_score0*(Mc/M) )   )
        #Ub = Ub_pre/mean_score0
        #md[Symbol(zscore_key*"_lb")] = Lb
        #md[Symbol(zscore_key*"_ub")] = Ub  
    end
    r=nlsolve(f!,[0.1])
    #println(r)
    zvalue=r.zero[1]
    pvalue=2.0 * ccdf(Normal(), abs(zvalue))
    two_tail = 1-pvalue     
    one_tail = 1-(pvalue/2)
    return one_tail, two_tail
end

function applyWeights(mocc::MOcc,mdolocc::MDolOcc,mpen::MPen)
    for m in [mocc,mdolocc,mpen]
        println("Processing : ",m.modelName)
        m.reff.sdf = join(m.reff.sdf,cnts.sdf[[:key,:weight]], on = :key)
        sdf=m.reff.sdf
        #sdf[:factor]=0.0  #unadj_mean_score0
        fdf=deepcopy(m.reff.sdf)
        for c in [:unadj_avg_expsd_hh_pre,:unadj_avg_cntrl_hh_pre,:unadj_avg_expsd_hh_pst,:unadj_avg_cntrl_hh_pst,:adj_mean_expsd_grp,:adj_mean_cntrl_grp,:mean_score0,:mean_score1]
            fdf[c] = map((x,y) -> (x*y), fdf[c],fdf[:weight])
        end
        for col in [:unadj_avg_expsd_hh_pre,:unadj_avg_cntrl_hh_pre,:unadj_avg_expsd_hh_pst,:unadj_avg_cntrl_hh_pst,:adj_mean_expsd_grp,:adj_mean_cntrl_grp,:mean_score0,:mean_score1]
            tot=m.feff.sdf[1,col]
            scol=string(col)
            sdf[Symbol("f_unweighted_"*scol)] = sdf[col]
            sdf[Symbol("f_"*scol)] = 0.0
            lst=by(fdf[!isnan(fdf[col]),:], :class, df -> sum(df[col]))
            rename!(lst, :x1,:creative_weight)
            lst[:factor] = tot ./ lst[:creative_weight]             
            for row in eachrow(lst)
                f=row[:factor]
                sdf[sdf[:class].==row[:class], Symbol("f_"*scol) ] = f
                sdf[sdf[:class].==row[:class], col ] = sdf[sdf[:class].==row[:class], col ] * f
            end
        end          
        sdf[:f_unweighted_onetail_pval] = sdf[:onetail_pval]
        sdf[:f_unweighted_twotail_pval] = sdf[:twotail_pval]
        sdf[:f_unweighted_adj_dod_effct] = sdf[:adj_dod_effct]
        
        sdf[:adj_dod_effct] = ((sdf[:adj_mean_expsd_grp] - sdf[:adj_mean_cntrl_grp])  ./ sdf[:adj_mean_cntrl_grp] ) * 100   
        for i in 1:length(sdf[1])   
            d=df2dict(sdf[i,:])
            p1, p2 = np_val(d)
            p1=p1*100
            p2=p2*100
            sdf[i,:onetail_pval] = p1
            sdf[i,:twotail_pval] = p2
        end
        
       # #RE-Default exposed Break Ctrl cols to total
        ##exposed
        #for col in [:unadj_avg_cntrl_hh_pre,:unadj_avg_cntrl_hh_pst,:adj_mean_cntrl_grp]
        #    totval=m.feff.sdf[1,col]
        #    sdf[sdf[:exposed].==true,col] = totval
        #end
        ##    sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pre] = mpen.feff.sdf[1,:unadj_avg_cntrl_hh_pre]
        ##    sdf[sdf[:key].==k ,:unadj_avg_cntrl_hh_pst] = mpen.feff.sdf[1,:unadj_avg_cntrl_hh_pst] 
        
        
    end
end
if cfg[:counts]
    applyWeights(mocc,mdolocc,mpen)
end
 
# -------- Confidence Intervals - Random Effects --------
function ConfidenceIntervals(cnts::Cnts,reff::RanEffect)
    sdf = reff.sdf
    for row in eachrow(sdf)
        k = row[:key]
        md = JStack.df2dict(row)    #[:B1_combo,:SE1_combo,:mean_score0,:mean_score1,:M,:Mt,:Mc]
        md[:B] =  md[:B1_combo] # don't really need these...need to clean up
        md[:SE] =  md[:SE1_combo] # ditto
            
#        md[:M_F] = cnts.sdf[1,:M]
#        md[:Mt_F] = cnts.sdf[1,:Mt]
#        md[:Mc_F] = cnts.sdf[1,:Mc]
        #if row[:exposed]
        #    cis =  CIs(md, ZDict, 1)
        #else
        #    cis =  CIs(md, ZDict, 2)
        #end
        cis =  CIs(md, ZDict)
        for (zscore_key,zscore) in  ZDict     
                ubk = Symbol(zscore_key*"_ub")
                lbk=Symbol(zscore_key*"_lb")
                sdf[sdf[:key].==k, lbk] = md[lbk]*100
                sdf[sdf[:key].==k, ubk] = md[ubk]*100
        end
    end
end
ConfidenceIntervals(cnts, mocc.reff)
ConfidenceIntervals(cnts, mdolocc.reff)
ConfidenceIntervals(cnts, mpen.reff)

#mdolhh.sdf=mdolhh.sdf[mdolhh.sdf[:key].!="creative_groups (none)",:]
    
    
end # :TotalModelsOnly    
    


# ==========================================================================================
# =================================== DolHH ================================================
# ==========================================================================================

mdf=MDF(cnts,mocc,mdolocc,mpen)
rdf=RDF(mocc,mdolocc,mpen,mdolhh)

function aggregate!(df_data::DataFrame, mdolhh::MDolHH, smdf::DataFrame)   #mdf::MDF
    sdf=mdolhh.sdf
    xsdf = join(mdolhh.sdf, cnts.sdf[[:key,:class,:level,:exposed]], on = :key)
    for row = eachrow(xsdf)
        k = row[:key]
        #if k != "Total Campaign"
            ranfx = row[:class]
            v_level = row[:level]
            exposed = row[:exposed]
            println("getting : ",k)
            #md = df2dict(mdf.get(k))
            md = JStack.df2dict(getk(smdf,k))
            adjctl = md[:o_mean_score0] * md[:y_mean_score0] * md[:p_mean_score0]
            adjexp = md[:o_mean_score1] * md[:y_mean_score1] * md[:p_mean_score1]
            println(k," ~~ ",adjctl,"~",adjexp) #,md)            
            sdf[sdf[:key].==k,:adj_mean_cntrl_grp] = adjctl
            sdf[sdf[:key].==k,:adj_mean_expsd_grp] = adjexp       
            sdf[sdf[:key].==k,:adj_dod_effct] = ((adjexp - adjctl) ./ adjctl ) * 100
            if k == "Total Campaign"
                sdf[sdf[:key].==k,:unadj_avg_cntrl_hh_pre] =  mean(df_data[ (df_data[:group] .== 0)  , :prd_1_net_pr_pre] )
                sdf[sdf[:key].==k,:unadj_avg_expsd_hh_pre] =  mean(df_data[ (df_data[:group] .== 1)  , :prd_1_net_pr_pre] )
                sdf[sdf[:key].==k,:unadj_avg_cntrl_hh_pst] =  mean(df_data[ (df_data[:group] .== 0)  , :prd_1_net_pr_pos] )
                sdf[sdf[:key].==k,:unadj_avg_expsd_hh_pst] =  mean(df_data[ (df_data[:group] .== 1)  , :prd_1_net_pr_pos] )
            else
                sdf[sdf[:key].== k,:unadj_avg_expsd_hh_pre] = mean(df_data[ (df_data[:group] .== 1) & (df_data[Symbol(lowercase(ranfx))] .== v_level) , :prd_1_net_pr_pre] )
                sdf[sdf[:key].== k,:unadj_avg_expsd_hh_pst] = mean(df_data[ (df_data[:group] .== 1) & (df_data[Symbol(lowercase(ranfx))] .== v_level) , :prd_1_net_pr_pos] )
                if exposed
                    sdf[sdf[:key].== k,:unadj_avg_cntrl_hh_pre] = sdf[1,:unadj_avg_cntrl_hh_pre] 
                    sdf[sdf[:key].== k,:unadj_avg_cntrl_hh_pst] = sdf[1,:unadj_avg_cntrl_hh_pst]     
                else
                    sdf[sdf[:key].== k,:unadj_avg_cntrl_hh_pre] = mean(df_data[ (df_data[:group] .== 0) & (df_data[Symbol(lowercase(ranfx))] .== v_level) , :prd_1_net_pr_pre] )
                    sdf[sdf[:key].== k,:unadj_avg_cntrl_hh_pst] = mean(df_data[ (df_data[:group] .== 0) & (df_data[Symbol(lowercase(ranfx))] .== v_level) , :prd_1_net_pr_pos] ) 
                end            
            end
        #end
    end
end
aggregate!(df_data,mdolhh,genMDF(cnts,mocc,mdolocc,mpen))
#aggregate!(df_data,mdolhh,genMDF(cnts, mocc, mdolocc, mpen))

 


function calc_rawPvals_Opt(iDict::OrderedDict)
    dout = iDict 
    M = get(iDict, :M, NA)
    Mt = get(iDict, :Mt, NA)
    Mc = get(iDict, :Mc, NA)
    N = get(iDict, :N, NA)
    Nt = get(iDict, :Mt, NA)
    Nc = get(iDict, :Nc, NA)
    B1 = get(iDict, :B1, NA)
    B2 = get(iDict, :B2, NA)
    B3 = get(iDict, :B3, NA)
    SE1 = get(iDict, :SE1, NA)
    SE2 = get(iDict, :SE2, NA)
    SE3 = get(iDict, :SE3, NA)
    SEsq=sqrt(SE1^2+SE2^2+SE3^2)
    o_mean_score0 = get(iDict, :o_mean_score0, NA)
    o_mean_score1 = get(iDict, :o_mean_score1, NA)
    y_mean_score0 = get(iDict, :y_mean_score0, NA)
    y_mean_score1 = get(iDict, :y_mean_score1, NA)
    p_mean_score0 = get(iDict, :p_mean_score0, NA)
    p_mean_score1 =get(iDict, :p_mean_score1, NA)
    Bsum=B1+B2+B3
    dout[:Bsum] = Bsum  
    ###### PVALUE - ONE & TWO ########
    m=nothing
    m = Model(solver=NLoptSolver(algorithm=:LD_MMA, maxtime=v_ttl))
    if B1 < 0
        @variable(m, Bocc >= B1)
    else
        @variable(m, Bocc <= B1)
    end
    if B2 < 0
        @variable(m, Bdolocc >= B2)
    else
        @variable(m, Bdolocc <= B2)
    end
    if B3 < 0
        @variable(m, Bpen >= B3)
    else
        @variable(m, Bpen <= B3)
    end
    @objective(m, Max, (((Bocc+Bpen+Bdolocc)-Bsum)/SEsq ))   # z = (x-u)/se
    @NLconstraint(m, 0.00000 <= ((((p_mean_score1*(Nt/N))+(p_mean_score0*exp(Bpen)*(Nc/N)))
	                           * ((o_mean_score1*(Mt/M))+(o_mean_score0*exp(Bocc)*(Mc/M)))
	                           * ((y_mean_score1*(Mt/M))+(y_mean_score0*exp(Bdolocc)*(Mc/M)))
	                             )
	                           -(((p_mean_score1*(Nt/N)*exp(-Bpen))+(p_mean_score0*(Nc/N)))
	                           *((o_mean_score1*(Mt/M)*exp(-Bocc))+(o_mean_score0*(Mc/M)))
	                           *((y_mean_score1*(Mt/M)*exp(-Bdolocc))+(y_mean_score0*(Mc/M)))
	                             )
	                          ) 
	               <= 0.00001
                    )
    #print(m)
    status = solve(m)
    zvalue=getobjectivevalue(m)
    pvalue=2.0 * ccdf(Normal(), abs(zvalue))
    two_tail = 1-pvalue     
    one_tail = 1-(pvalue/2)
    dout[:onetail_pval] = one_tail
    dout[:twotail_pval] = two_tail
    println("z-value: ", string(zvalue)," --> p-value: ",string(two_tail))
    return dout           
end

function rawPvals(mdolhh::MDolHH,x::DataFrame)
    sdf=mdolhh.sdf       
    sdf[:twotail_pval_raw] = 0.0
    sdf[:onetail_pval_raw]= 0.0
    for row = eachrow(x)
        k = row[:key]
        if k != "Total Campaign"
            md = JStack.df2dict(x[x[:key].==k,:])
            md[:metakey] = k
            println("Raw PValue : ",k)
            if :B1_orig in keys(md)  #!!!!!!!!!
                md[:B1] = md[:B1_orig]
                md[:SE1] = md[:SE1_orig]
                md[:B2] = md[:B2_orig]
                md[:SE2] = md[:SE2_orig]
                md[:B3] = md[:B3_orig]
                md[:SE3] = md[:SE3_orig]
            println(md)
                calc_rawPvals_Opt(md)
                sdf[sdf[:key].==k, :onetail_pval_raw]=md[:onetail_pval]*100
                sdf[sdf[:key].==k, :twotail_pval_raw]=md[:twotail_pval]*100
            end
        end
    end
end
#rawPvals(mdolhh,genRndMDF(cnts, mocc, mdolocc, mpen, false))
rawPvals(mdolhh,genMDF(cnts, mocc, mdolocc, mpen))
#mdolhh.sdf[[:key,:twotail_pval,:onetail_pval,:twotail_pval_raw,:onetail_pval_raw]]


function ConfidenceIntervals(mdolhh::MDolHH,x::DataFrame)
    sdf=mdolhh.sdf       
    for row = eachrow(x)
        k = row[:key]
        md = JStack.df2dict(x[x[:key].==k,:])
        md[:metakey] = k
        println("getting : ",k)
        calcPValue_Opt(md)
        CIs_O(md)
        for zk in [ :onetail_pval, :twotail_pval, :onetail_80_pct_intrvl_lb,:onetail_80_pct_intrvl_ub, :onetail_90_pct_intrvl_lb, :onetail_90_pct_intrvl_ub
                    ,:twotail_80_pct_intrvl_lb, :twotail_80_pct_intrvl_ub, :twotail_90_pct_intrvl_lb, :twotail_90_pct_intrvl_ub    
                  ]
            sdf[ sdf[:key].==k, zk]= md[zk]*100
        end
    end
end
#ConfidenceIntervals(mdolhh,genRndMDF(cnts, mocc, mdolocc, mpen, false))
#ConfidenceIntervals(mdolhh,genFixedMDF(cnts, mocc, mdolocc, mpen))
ConfidenceIntervals(mdolhh,genMDF(cnts, mocc, mdolocc, mpen))


#function vcatreff(mocc, mdolocc, mpen)
#    #dfo=DataFrame()
#    for m in [mocc,mdolocc, mpen]
#        if m.hasBreaks
#            
#    end
#end
#    function rmap(m::MModel)
#        return m.hasBreaks
#    end
#        rvs=vcat(map(x->if x.hasBreaks rmap(x) end,[mocc,mdolocc, mpen]))



function extendRDF(rdf::DataFrame)
    dfo=DataFrame(MODEL_DESC=String[], Model=String[], TIME_AGG_PERIOD=Int64[], START_WEEK=Int64[], END_WEEK=Int64[], dependent_variable=String[],
                  CNT_EXPSD_HH=Int64[], UNADJ_AVG_EXPSD_HH_PRE=Float64[], UNADJ_AVG_CNTRL_HH_PRE=Float64[], UNADJ_AVG_EXPSD_HH_PST=Float64[], 
                  UNADJ_AVG_CNTRL_HH_PST=Float64[], UNADJ_DOD_EFFCT=Float64[], UNADJ_DIFF_EFFCT=Float64[], ADJ_MEAN_EXPSD_GRP=Float64[],
                  ADJ_MEAN_CNTRL_GRP=Float64[], ADJ_DOD_EFFCT=Float64[], TWOTAIL_PVAL=Float64[], ONETAIL_PVAL=Float64[], ABS_DIFF=Float64[], 
                  DOL_DIFF=Float64[], ONETAIL_80_PCT_INTRVL_UB=Float64[], ONETAIL_80_PCT_INTRVL_LB=Float64[], ONETAIL_90_PCT_INTRVL_UB=Float64[],
                  ONETAIL_90_PCT_INTRVL_LB=Float64[], TWOTAIL_80_PCT_INTRVL_UB=Float64[], TWOTAIL_80_PCT_INTRVL_LB=Float64[], TWOTAIL_90_PCT_INTRVL_UB=Float64[],
                  TWOTAIL_90_PCT_INTRVL_LB=Float64[], CNT_IMPRESSIONS=Int64[], TWOTAIL_PVAL_to_Campaign=Float64[], ONETAIL_PVAL_to_Campaign=Float64[],
                  CNT_Model_HH=Float64[]
                 )

  ###NBNBNB : this needs to come out - only to hack total only        
  #          rvs = vcat(mocc.reff.sdf[[:key,:model,:twotail_pval_raw,:onetail_pval_raw]],                                    
#               mdolocc.reff.sdf[[:key,:model,:twotail_pval_raw,:onetail_pval_raw]],
#               mpen.reff.sdf[[:key,:model,:twotail_pval_raw,:onetail_pval_raw]],
##               mdolhh.sdf[[:key,:model,:twotail_pval_raw,:onetail_pval_raw]]
#              )
    if (!cfg[:TotalModelsOnly]) & (length(filter(x-> x.hasBreaks ,[mocc,mdolocc, mpen]))  > 0 )
            rvs=vcat( map(x-> x.reff.sdf, filter(x-> x.hasBreaks ,[mocc,mdolocc, mpen])))[[:key,:model,:twotail_pval_raw,:onetail_pval_raw]]
            rvs = vcat(rvs,mdolhh.sdf[[:key,:model,:twotail_pval_raw,:onetail_pval_raw]])
    else
        rvs = mdolhh.sdf[[:key,:model,:twotail_pval_raw,:onetail_pval_raw]]
    end
        
    for i in 1:length(rdf[1])
        push!(dfo, [
                     rdf[i,:key], 
                     NA, #Model=NA, # code
                     NA, #TIME_AGG_PERIOD=NA, 
                     NA, #START_WEEK=NA, 
                     NA, #END_WEEK=NA, 
                     rdf[i,:model],
                     NA, #CNT_EXPSD_HH=Int64[], 
                     rdf[i,:unadj_avg_expsd_hh_pre], 
                     rdf[i,:unadj_avg_cntrl_hh_pre], 
                     rdf[i,:unadj_avg_expsd_hh_pst], 
                     rdf[i,:unadj_avg_cntrl_hh_pst], 
                     NA, #UNADJ_DOD_EFFCT=Float64[], 
                     NA, #UNADJ_DIFF_EFFCT=Float64[], 
                     rdf[i,:adj_mean_expsd_grp],
                     rdf[i,:adj_mean_cntrl_grp], 
                     rdf[i,:adj_dod_effct], 
                     rdf[i,:twotail_pval], 
                     rdf[i,:onetail_pval], 
                     NA, #ABS_DIFF=Float64[], 
                     NA, #DOL_DIFF=Float64[], 
                     rdf[i,:onetail_80_pct_intrvl_ub], 
                     rdf[i,:onetail_80_pct_intrvl_lb], 
                     rdf[i,:onetail_90_pct_intrvl_ub],
                     rdf[i,:onetail_90_pct_intrvl_lb], 
                     rdf[i,:twotail_80_pct_intrvl_ub], 
                     rdf[i,:twotail_80_pct_intrvl_lb], 
                     rdf[i,:twotail_90_pct_intrvl_ub],
                     rdf[i,:twotail_90_pct_intrvl_lb], 
                     NA, #CNT_MPRESSIONS=Float64[], 
                     NA, #TWOTAIL_PVAL_to_Campaign=Float64[], 
                     NA, #ONETAIL_PVAL_to_Campaign=Float64[],
                     NA #CNT_Model_HH=Float64[]      
                   ]
             )
        x=length(dfo[1])
        dfo[x,:UNADJ_DOD_EFFCT] = ( ((dfo[x,:UNADJ_AVG_EXPSD_HH_PST] - dfo[x,:UNADJ_AVG_EXPSD_HH_PRE]) - (dfo[x,:UNADJ_AVG_CNTRL_HH_PST] - dfo[x,:UNADJ_AVG_CNTRL_HH_PRE]))  /  dfo[x,:UNADJ_AVG_CNTRL_HH_PST] ) *100
        
        dfo[x,:UNADJ_DIFF_EFFCT] = ((dfo[x,:UNADJ_AVG_EXPSD_HH_PST] - dfo[x,:UNADJ_AVG_CNTRL_HH_PST]) / dfo[x,:UNADJ_AVG_CNTRL_HH_PST] )* 100 
    
        if dfo[i,:dependent_variable] == "dolhh"
            k=dfo[x,:MODEL_DESC]
            if cfg[:counts]
                dfo[x,:CNT_EXPSD_HH] =  getk(cnts.sdf,k,:hh)   #  cnts.get(k,:hh)
                dfo[x,:CNT_IMPRESSIONS] =  getk(cnts.sdf,k,:impressions)  #cnts.get(k,:impressions)
            end
            dfo[x,:DOL_DIFF] = dfo[x,:ADJ_MEAN_EXPSD_GRP] - dfo[x,:ADJ_MEAN_CNTRL_GRP]
            dfo[x,:CNT_Model_HH] = cnts.sdf[cnts.sdf[:key].==dfo[x,:MODEL_DESC],:M][1]
            
            
        elseif dfo[i,:dependent_variable] == "pen"
            dfo[x,:CNT_Model_HH] = cnts.sdf[cnts.sdf[:key].==dfo[x,:MODEL_DESC],:Nt][1]
        else
            dfo[x,:CNT_Model_HH] = cnts.sdf[cnts.sdf[:key].==dfo[x,:MODEL_DESC],:Mt][1]
        end        
        
        if (dfo[x,:dependent_variable] in ["occ","dolocc","pen","dolhh"]) & (dfo[x,:MODEL_DESC] !== "Total Campaign")
            
            
            pv = rvs[(rvs[:key].== dfo[x,:MODEL_DESC])&(rvs[:model].== dfo[x,:dependent_variable]) ,:]
            if length(pv[1]) >=1
                dfo[i,:TWOTAIL_PVAL_to_Campaign] = pv[:twotail_pval_raw][1]
                dfo[i,:ONETAIL_PVAL_to_Campaign] = pv[:onetail_pval_raw][1]
            end
        end
         
        #format for unify
        dfo[dfo[:dependent_variable].=="occ", :dependent_variable] = "OCC" 
        dfo[dfo[:dependent_variable].=="dolocc", :dependent_variable] = "DOL/OCC" 
        dfo[dfo[:dependent_variable].=="pen", :dependent_variable] = "PEN"
        dfo[dfo[:dependent_variable].=="dolhh", :dependent_variable] = "DOL/HH"
        
    end
    dfo[:ABS_DIFF] = dfo[:ADJ_MEAN_EXPSD_GRP] - dfo[:ADJ_MEAN_CNTRL_GRP]
    
    return dfo
end
xrdf=extendRDF(genRDF(mocc,mdolocc,mpen, mdolhh,cfg) )

#xrdf[[:MODEL_DESC,:dependent_variable,:ADJ_MEAN_EXPSD_GRP,:ADJ_MEAN_CNTRL_GRP]]
#x=genMDF(cnts,mocc,mdolocc,mpen)
#x[[:key,:o_mean_score0,:y_mean_score0,:p_mean_score0,:o_mean_score1,:y_mean_score1,:p_mean_score1]]
#write2disk(xrdf,"/home/rmadmin/g/StatStack/src/xrdf.csv")
write2disk(xrdf,root*"/scored.csv")









