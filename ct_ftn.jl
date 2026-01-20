
"""
2024.10.23
    Contact tracing code with missing information
    IBM Julia language code developed by Min-Kyung Chae (in KIAS).
    Multi-thread code:     julia --threads 2 filename.jl

    contact tracing 의 불확실성. (number_control_infection 보다 감염자수가 많아지면,)
        1. probability of contact tracing : ratio_ct = [1, 0.95, 0.95, 0.8, 0.1] -> [0.9, 0.9, 0.9, 0.4, 0.05]
        2. probability of being reported : contact tracing을 하는 사람의 숫자 자체를 변경. p_ = 1.0 -> 0.5

    10/23 전체적으로 코드 만들기 + 최종본
    10/28 친구 만날 확률 고침 

    고쳐야 할 것.
    - 친구 만날 확률 > 10/28 고침

    고려해야 할 것.
    지하철 인구 이동을 넣어야 하나? < 고민중 : 데이터가 없음.

"""

using Base.Threads;
using Random, Distributions, StatsBase;
using DataFrames, CSV;

function Viral_Shedding()
    vs_ = zeros(10)
    vs_[1] = 0.16098790572752392; vs_[2] = 0.2412521255679755; vs_[3] = 0.2080390352727106; vs_[4] = 0.15005495814415215; vs_[5] = 0.09895714244707202;
    vs_[6] = 0.06182814230417574; vs_[7] = 0.03725825991796377; vs_[8] = 0.021877903740751214; vs_[9] = 0.012598721837235396; vs_[10] = 0.007145805040439527;

    return vs_
end

function Read_CSV_File()
    df1 = DataFrame(CSV.File("./../ct/input/synthetic_population_edueco_v2.csv"))
    df2 = DataFrame(CSV.File("./../ct/input/synthetic_population_eco_v2.csv"))
    df3 = DataFrame(CSV.File("./../ct/input/synthetic_population_frd_v2.csv"))

    #df1 = DataFrame(CSV.File("/home/mkchae/ct/input/synthetic_population_edueco_v2.csv"))
    #df2 = DataFrame(CSV.File("/home/mkchae/ct/input/synthetic_population_eco_v2.csv"))
    #df3 = DataFrame(CSV.File("/home/mkchae/ct/input/synthetic_population_frd_v2.csv"))

    df_ = outerjoin(df1, df2; on=[:id, :hh_id, :rc, :age])
    df_ = transform(df_, :age => ByRow(x -> searchsortedlast(vcat(0:10:80, 1000), x)) => :age_group)
    df_ = outerjoin(df_, df3; on=[:id, :hh_id, :rc, :age])
    sort!(df_)

    df_[!, 1:11] = convert.(Int64, df_[:, 1:11])
    df_[!, :"tmp"] .= [Vector{Int64}()]
    for i in 1:nrow(df_)
        if df_[i, :n_frd] != 0
            list = parse.(Int, split(chop(df_[i, :frd]; head=1, tail=1), ','))
            df_[i, :tmp] = list
            list = nothing
        end
    end
    select!(df_, Not(:frd))
    rename!(df_, :tmp => :frd)

    println(" population # = ", nrow(df1))
    println(" population # = ", nrow(df2))
    println(" population # = ", nrow(df3))
    println(" population # = ", nrow(df_))
    
    return df_
end

# state of agents
function Make_df_State(N_)
    df_ = DataFrame("id" => 1:N_, "state" => 0, 
                      "date_E" => 0, "date_I" => 0, "date_R" => 0,
                      "period_E" => 0, "period_I" => 0,
                      "count_s" => 0,
                      # (self-)quarantine + isolation
                      "state_q" => 0,
                      "date_SQ" => 0, "date_IS" => 0,
                      "period_SQ" => 0, "period_IS" => 0, "period_Q" => 0, # Q는 여러번 될 수 있으니까 날짜를 따로 기록하지 않음.
                      "count_q" => 0,
                      # relative infectiousness
                      "RI" => 0.0,
		      # count infection
		      "count_inf" => 0
		      )

    df_[!, :"I_type"] .= [Vector{Int64}()]
    df_[!, :"I_pop"] .= [Vector{Int64}()]
    sort!(df_, [:id])
    
    return df_
end

# qurantine of agents
function Make_df_Qurantine(N_)
    df_ = DataFrame("id" => 1:N_, "count_q" => 0)

    sort!(df_, [:id])    
    return df_
end


# virus of agents
function Make_df_Viral(N_)
    df_ = DataFrame("id" => 1:N_, "state" => 0, 
                      "A_vs" => 0.0, "period_no" => 0, "period_emit" => 0,
                      "count_vs" => 0)
    sort!(df_, [:id])

    return df_
end

# irregular contact data of agents
function Make_df_Daily(N_)
    df_ = DataFrame("id" => 1:N_, "friends_out" => 0, "region" => 0, "encounters_out" => 0)
    sort!(df_, [:id])
    
    return df_
end

# friends meeting data of agents - record meeting number
function Make_df_Daily_Record(N_, t_)
    df_ = [DataFrame("id" => 1:N_, "num_f" => 0, "num_e" => 0) for i in 1:t_]
    for i in 1:t_
        sort!(df_[i], [:id])
    end

    return df_
end


# no contact list (per day) # 격리되어있는 사람을 적는 것 
function Make_df_List(ntries_)
    df_ = DataFrame("t" => 1:ntries_+1)
    df_[!, :list_SQ] .= [Vector{Int64}()]
    df_[!, :list_Q] .= [Vector{Int64}()]
    df_[!, :list_IS] .= [Vector{Int64}()]
    df_[!, :list_CT] .= [Vector{Int64}()] 
    sort!(df_, [:t])

    return df_
end

# no contact list (per day) # 격리되어있는 사람을 수 적는 것 
function Make_df_ListCount(ntries_)
    df_ = DataFrame("t" => 1:ntries_+1, "list_SQ" => 0, "list_Q" => 0, "list_IS" => 0, "list_CT" => 0)

    sort!(df_, [:t])
    return df_
end

# friends meeting data per meeting (per day)
function Make_df_Contact_List(N_, t_)
    n = ceil(Int64, N_/1.5)

    df_ = [DataFrame("meeting_num" => 1:n) for i in 1:t_]
    for i in 1:t_
        df_[i][!, :list_f] .= [Vector{Int64}()]
        df_[i][!, :time_f] .= 0.0
        df_[i][!, :list_e] .= [Vector{Int64}()]
        df_[i][!, :time_e] .= 0.0
        sort!(df_[i], [:meeting_num])
    end
    
    return df_
end

# dataframe initial setting
function InitSetting_Df(dfs_, dfv_, dft_, dfd_, df_l_, df_cl_, t_)
    dfs_[:, 2:15] .= 0; dfs_[:, 16] .= 0.0; dfs_[:, 17] .= 0; dfs_[:, 18:19] .= [Vector{Int64}()];
    sort!(dfs_, [:id])

    dfv_[:, 2] .= 0; dfv_[:, 3] .= 0.0; dfv_[:, 4:6] .= 0;
    sort!(dfv_, [:id])

    dft_[:, 2:4] .= 0;
    sort!(dft_, [:id])

    for i in 1:t_
        dfd_[i][:, 2:3] .= 0
        sort!(dfd_[i], [:id])
    end

    df_l_[:, 2:5] .= [Vector{Int64}()];
    sort!(df_l_, [:t])

    for i in 1:t_
        df_cl_[i][:, [2, 4]] .= [Vector{Int64}()]
        df_cl_[i][:, [3, 5]] .= 0.0
        sort!(df_cl_[i], [:meeting_num])
    end
end

function Compute_Chance(rngs_, p_)
    rn_ = (rand(rngs_) < p_ ? 1 : 0)
    
    return rn_
end

# td = 0~1 (1 = 12h = 720m)
function Print_Duration(rngs_, n_)
    pf = zeros(5)
    if n_ == 1
        pf = [0.01, 0.027054054054054054, 0.08189189189189189, 0.1599189189189189, 0.7211351351351352]
    elseif n_ == 2
        pf = [0.0155063629558336, 0.059779702705592984, 0.15688161693936478, 0.2802908779809646, 0.487541439418244]
    elseif n_ == 3
        pf = [0.07885250063467886, 0.17065244986037065, 0.2395024117796395, 0.2255902513328256, 0.2854023863924854]
    elseif n_ == 4
        pf = [0.022328548644338118, 0.061765985210961286, 0.22618529795563289, 0.47020443671161377, 0.21951573147745396]
    elseif n_ == 5
        pf = [0.46633499170812603, 0.2620232172470978, 0.16583747927031509, 0.08606965174129354, 0.019734660033167494]
    else
        println("????????? in Print Duration")
    end
    num = sample(rngs_, 1:5, Weights(pf), 1)
    td_ = 0.0

    if num[1] == 1                               #less than 5 min
        td_ = rand(rngs_, Uniform(0, 5))
    elseif num[1] == 2                            # 5 - 15 min
        td_ = rand(rngs_, Uniform(5, 15))
    elseif num[1] == 3                            # 15 - 60 min
        td_ = rand(rngs_, Uniform(15, 60))
    elseif num[1] == 4                            # 1 - 4h 
        td_ = rand(rngs_, Uniform(60, 240))
    elseif num[1] == 5                            # more than 4h # 최대 12시간 까지 가능
        td_ = rand(rngs_, Uniform(240, 720))
    end
    
    return td_/720 
end

# 친구 만날 사람을 정해 주고, 누굴 만날지 정해 준다.
function Friends_Meeting(rngs_, num_, p_, df_, dfs_, dft_, dfd_, df_cl_)
    # 만남을 정해 주기 전에 dataframe 을 reset 한다. 
    dfd_[num_][:, :num_f] .= 0;
    df_cl_[num_][:, :list_f] .= [Vector{Int64}()];
    dft_[:, :friends_out] .= 0;

    # 오늘 친구를 만나러 나갈 사람 결정.
    dft_[:, :friends_out] .= Compute_Chance.(rngs_, p_); 
    # 오늘 외출을 하고, no contact 대상자가 아니야.
    list = shuffle(rngs_, findall((dft_.friends_out .== 1) .& (dfs_.state_q .== 0)))
    #println(" 오늘 친구를 만나러 갈 사람 # = ", length(list))
    #println(" 친구 만나러 나갈 사람 골라 = ", list[1:10])

    outing_status = deepcopy(dft_[:, :friends_out])
    meeting_status = zeros(Int64, length(outing_status))
    
    meeting_num = 1
    for i in 1:length(list)
        id = list[i]; 
        #println(" i = ", i, " id = ", id)
    
        # 내가 아직 약속이 없으면 약속을 정해야지
        if meeting_status[id] == 0

            # 내 친구 목록
            list_frd = df_[id, :frd]
            # 그 중 밖에 나갈 계획이지만, 아직 약속 없음.
            list_frd_O = list_frd[findall((outing_status[list_frd] .== 1) .& (meeting_status[list_frd] .== 0))]

            # 몇 명 만날지 정해
            n = length(list_frd_O)
            n_size = (n > 19 ? 19 : n)

            if n_size != 0
                list_frd_O = shuffle(rngs_, list_frd_O)

                # 만날 사람을 뽑앙.
                list_meeting = Vector{Int64}()
                append!(list_meeting, list_frd_O[1:rand(rngs_, 1:n_size)], id) 
                meeting_status[list_meeting] .= meeting_num

                dfd_[num_][list_meeting, :num_f] .= meeting_num
                df_cl_[num_][meeting_num, :list_f] = list_meeting
                df_cl_[num_][meeting_num, :time_f] = Print_Duration(rngs_, 4)                

                meeting_num += 1
            else
                outing_status[id] = 0
                dft_[id, :friends_out] = 0
            end
        end
    end
end

# 우연히 길에서 모르는 사람을 만남
function Random_Encounter(rngs_, num_, p_, df_, dfs_, dft_, dfd_, df_cl_)

    encounters_avg = 5

    # 만남을 정해 주기 전에 dataframe 을 reset 한다.
    df_cl_[num_][:, :list_e] .= [Vector{Int64}()];
    dfd_[num_][:, :num_e] .= 0;
    dft_[:, :encounters_out] .= 0;
    dft_[:, :region] .= df_[:, :rc]

    # 오늘 외출할 사람 결정.
    dft_[:, :encounters_out] .= Compute_Chance.(rngs_, p_)

    meeting_num = 1
    group = groupby(dft_, [:region])
    for i in 1:length(group)
        tmp = group[i][:, :id]
    
        # 해당 지역에서 우연한 만남을 할 수 있는 사람 목록
        # 오늘 외출을 하고, # no contact 대상자가 아니야.
        list = tmp[findall((dft_[tmp, :encounters_out] .== 1) .& (dfs_[tmp, :state_q] .== 0))]

        # 필요한 모임 개수
        n = round(Int64, length(list)/encounters_avg)
        
        # 각 사람마다 참여할 모임의 수를 할당해
        dfd_[num_][list, :num_e] .= rand(rngs_, 1:n, length(list)) .+ meeting_num
        meeting_num = maximum(dfd_[num_][list, :num_e])    
    end
    
    group = groupby(dfd_[num_], [:num_e])
    for i in 2:length(group) # group[1] 나가지 않는 individuals
        list = group[i][:, :id]
        n = group[i][1, :num_e]
        df_cl_[num_][n, :list_e] = list
        df_cl_[num_][n, :time_e] = Print_Duration(rngs_, 5)
    end
end


function Print_Period_E(rngs_)
    k = 1.926;     theta = 1.775;
    rn = round(Int64, rand(rngs_, truncated(Gamma(k, theta), 1, 10)))
    #rn = round(Int64, rand(rngs_, truncated(Normal(4, 1), 1, 7)))
    
    return rn
end

function Print_Period_I(rngs_)
    rn = 8
    
    return rn
end

function Print_Relative_Infectiousness(rngs_, std_)
    if std_ == 0
        rn = 1
    else
        mean_ = 1;
        k = (mean_*mean_)/(std_*std_);     theta = (std_*std_)/mean_;
        rn = round(rand(rngs_, Gamma(k, theta)), digits=3)
    end

    return rn
end


# New infectors setting.
function Setting_E_State(rngs_, list_, t_, dfs_, dfv_)
    len = length(list_)
    
    dfs_[list_, :state] .= 1
    dfs_[list_, :date_E] .= t_ + 1
    dfs_[list_, :period_E] .= [Print_Period_E(rngs_) for i in 1:len]
    dfs_[list_, :period_I] .= [Print_Period_I(rngs_) for i in 1:len]
    dfs_[list_, :count_s] .= 1
    dfs_[list_, :RI] .= [Print_Relative_Infectiousness(rngs_, 0.5) for i in 1:len]
    
    dfv_[list_, :period_no] .= [(dfs_[list_[i], :period_E] < 3 ? 1 : dfs_[list_[i], :period_E] - 2) for i in 1:len]
    dfv_[list_, :period_emit] .= dfs_[list_, :period_I] .+ dfs_[list_, :period_E] .- dfv_[list_, :period_no]
    dfv_[list_, :A_vs] .= [(dfs_[list_[i], :period_E] == 1 ? 0.6000404822402197 : (dfs_[list_[i], :period_E] == 2 ? 0.8422130085822237 : 1.0)) for i in 1:len]
    dfv_[list_, :count_vs] .= 1
end


function Initial_Infectors_Rand(rngs_, t_, n_, dfs_, dfv_)
    #tmp = shuffle(rngs_, collect(1:nrow(dfs_)))
    #list_init0 = tmp[1:n_]
    list_init0 = [238593, 851606, 970618, 1621341, 2764780, 2868731, 3734217, 3746060, 4429319, 5369163, 5547237, 5649542, 6792047, 7310305, 7741843, 8261747, 8444813, 8465292, 9085879, 9528829]
    #list_init0 = [3860578, 2387352, 8880020, 2349520, 2182159, 6702556, 8430524, 2658550, 3301693, 2827828, 491006, 1973895, 7820447, 6110082, 8581873, 6597353, 598155, 6203774, 1264896, 1649721, 7972214, 6857659, 9034806, 3609916, 4306084, 466003, 6285285, 7351895, 3514875, 6031701, 3129638, 4127783, 9149576, 5107418, 6809405, 556957, 5578488, 5819175, 583229, 9491763]

    dfs_[list_init0, :I_type] .= [[0] for i in 1:n_]
    dfs_[list_init0, :I_pop] .= [[0] for i in 1:n_]
    Setting_E_State(rngs_, list_init0, t_, dfs_, dfv_)
end


# Part 1: Viral Shedding
function Check_Viral_State(id_, dfv_)
    if (dfv_[id_, :state] == 0) .& (dfv_[id_, :period_no] == dfv_[id_, :count_vs])
        dfv_[id_, :state] = 1
        dfv_[id_, :count_vs] = 1
    elseif (dfv_[id_, :state] == 1) .& (dfv_[id_, :period_emit] == dfv_[id_, :count_vs])
        dfv_[id_, :state] = 0
        dfv_[id_, :count_vs] += 1
    else
        dfv_[id_, :count_vs] += 1
    end
end


# Part 2: Disease Progression (E > I / I > R)
function Check_Disease_State(rngs_, id_, t_, Pa_, Psq_, dfs_)
    # state: 0 - susceptible / 1 - exposed / 2 - asymptomatic / 3 - symptomatic / 4 - R

    # E (1) > I (2, 3)
    if (dfs_[id_, :state] == 1) .& (dfs_[id_, :period_E] == dfs_[id_, :count_s])

        dfs_[id_, :date_I] = t_ + 1
        dfs_[id_, :count_s] = 1
        
        rn = rand(rngs_)
        # asymptomatic
        if rn < Pa_
            dfs_[id_, :state] = 2
            
        # symptomatic
        else
            dfs_[id_, :state] = 3

            # Self-Quarantine
            if (rand(rngs_) < Psq_) .& (dfs_[id_, :state_q] == 0)
                dfs_[id_, :period_SQ] = 1
                dfs_[id_, :date_SQ] = t_ + 1
                dfs_[id_, :state_q] = 1
                dfs_[id_, :count_q] = 1
            end
        end
        
    # I (2, 3) > R (4)
    elseif (2 <= dfs_[id_, :state] <= 3) .& (dfs_[id_, :period_I] == dfs_[id_, :count_s])
        dfs_[id_, :state] = 4
        dfs_[id_, :date_R] = t_ + 1

    else 
        dfs_[id_, :count_s] += 1
    end
end

# Part 3: Confinement
function Check_Confinement_State(rngs_, id_, t_, period_is_, dfs_)
    # state_q: 0 - carry on as usual / 1 - self-quarantine / 2 - quarantine / 3 - isolation
    list = Vector{Int64}(); # positive person > The contact tracing policy will be implemented
    
    # 1: self-quarantine > isolation (positive)
    if (dfs_[id_, :state_q] == 1) .& (dfs_[id_, :count_q] == dfs_[id_, :period_SQ])
        dfs_[id_, :state_q] = 3
        dfs_[id_, :date_IS] = t_ + 1
        dfs_[id_, :period_IS] = period_is_
        dfs_[id_, :count_q] = 1
        append!(list, id_)
        
    # 2: quarantine > carry on as usual or isolation
    elseif (dfs_[id_, :state_q] == 2) .& (dfs_[id_, :count_q] == dfs_[id_, :period_Q])
        # > carry on as usual (negative)
        if (dfs_[id_, :state] == 0) .| (dfs_[id_, :state] == 4)
            dfs_[id_, :count_q] = 0
            dfs_[id_, :state_q] = 0
            dfs_[id_, :period_Q] = 0

        # > isolation (positive)
        elseif ( 1 <= dfs_[id_, :state] <= 3 )
            dfs_[id_, :state_q] = 3
            dfs_[id_, :date_IS] = t_ + 1
            dfs_[id_, :period_IS] = period_is_
            dfs_[id_, :count_q] = 1
            append!(list, id_)
        else 
            println("@@@@@ ???????????? @@@@@@")
        end  

    # 3: isolation > carry on as usual
    elseif (dfs_[id_, :state_q] == 3) .& (dfs_[id_, :count_q] == dfs_[id_, :period_IS])
        # 아직 감염이 안 끝났으면, count_q 올리고, period_IS 추가로 늘리기.
        if (dfs_[id_, :state] != 4)
            dfs_[id_, :period_IS] += period_is_
            dfs_[id_, :count_q] += 1
            #println(" id = ", id_, " in check confinement state : ISO Extend the duration = ", dfs_[id_, :state], " q = ", dfs_[id_, :state_q])
        else 
            dfs_[id_, :count_q] = 0
            dfs_[id_, :state_q] = 0
        end
    else
        dfs_[id_, :count_q] += 1
    end

    return list
end


# 감염자 결정 
function Compute_Lambda(rngs_, id_, t_, day_, num_, list_q_, vs_, group_h_, group_c_, group_w_, df_, dfs_, dfv_, dfd_, df_cl_)
    # 새로운 감염자가 누구인지 찾아 볼까용
    list_new = Vector{Int64}()
    
    # 장소별 감염 scale vector : 가정, 학교, 회사, 친구, 우연
    B = [5.0, 5.0, 5.0, 5.0, 5.0];

    # 장소별 접촉자 리스트
    list = [Vector{Int64}() for i in 1:5]
    
    # 평일/주말 구분 (토/일 은 구분 안함)
    wc = (day_ <= 5 ? 0 : 1)
    
    # 가정내 접촉자.
    list_a = group_h_[df_[id_, :hh_id]][:, :id]
    #list[1] = setdiff(list_a[findall(dfs_[list_a, :state] .== 0)], id_, list_q_)
    list[1] = list_a[findall((dfs_[list_a, :state] .== 0) .& (dfs_[list_a, :id] .!= id_))]

    if wc == 0
        # 학교내 접촉자.
        if df_[id_, :class] != 0
            list_a = group_c_[df_[id_, :class]+1][:, :id]
            #list[2] = setdiff(list_a[findall(dfs_[list_a, :state] .== 0)], id_, list_q_)
            list[2] = list_a[findall((dfs_[list_a, :state] .== 0) .& (dfs_[list_a, :id] .!= id_))]
        # 회사내 접촉자.
        elseif df_[id_, :office] != 0
            list_a = group_w_[df_[id_, :office]+1][:, :id]
            #list[3] = setdiff(list_a[findall(dfs_[list_a, :state] .== 0)], id_, list_q_)
            list[3] = list_a[findall((dfs_[list_a, :state] .== 0) .& (dfs_[list_a, :id] .!= id_))]
        end
    end

    # 친구 모임 접촉자.
    if dfd_[num_][id_, :num_f] != 0
        #list_a = setdiff(df_cl_[num_][dfd_[num_][id_, :num_f], :list_f], id_)
        #list[4] = list_a[findall(dfs_[list_a, :state] .== 0)]
        list_a = df_cl_[num_][dfd_[num_][id_, :num_f], :list_f]
        list[4] = list_a[findall((dfs_[list_a, :state] .== 0) .& (dfs_[list_a, :id] .!= id_))]
    end

    # 우연한 모임 접촉자.
    if dfd_[num_][id_, :num_e] != 0
        #list_a = setdiff(df_cl_[num_][dfd_[num_][id_, :num_e], :list_e], id_)
        #list[5] = list_a[findall(dfs_[list_a, :state] .== 0)]
        list_a = df_cl_[num_][dfd_[num_][id_, :num_e], :list_e]
        list[5] = list_a[findall((dfs_[list_a, :state] .== 0) .& (dfs_[list_a, :id] .!= id_))]
    end

    for n in 1:5
        len = length(list[n]); 
        td = zeros(len);
        if len != 0
            if n == 4
                td .= df_cl_[num_][dfd_[num_][id_, :num_f], :time_f]
            else
                td .= Print_Duration.(rngs_, n)
            end

            ifn = vs_[dfv_[id_, :count_vs]-1] / dfv_[id_, :A_vs] * dfs_[id_, :RI]   # infectousness
            llambda = B[n] .* td .* ifn
            x = 1 .- exp.(-1 .* llambda)
            rn = rand(rngs_, len)
            ll = list[n][findall(rn .< x)]

            if length(ll) != 0
                append!(list_new, ll)
                dfs_[id_, :count_inf] += length(ll)
		        for i in 1:length(ll)
                    tmp1 = Vector{Int64}(); 
                    tmp2 = Vector{Int64}();
                    append!(tmp1, dfs_[ll[i], :I_type], n)
                    append!(tmp2, dfs_[ll[i], :I_pop], id_)
                    dfs_[ll[i], :I_type] = tmp1
                    dfs_[ll[i], :I_pop] = tmp2
                end
            end
        end
    end

    if length(list_new) != 0
        # 접촉자 셋팅. E 로 바꿔 줘야행. (그리고 필요한 값들 넣어 줘야 함.)
        Setting_E_State(rngs_, list_new, t_, dfs_, dfv_)
    end
end

# contact tracing : t_ 기준으로 (t_-1-period_dy ~ t_-period_ct-period_dy) 기간 동안 id_가 만난 사람들을 체크해 준다. 
function Contact_Tracing(rngs_, t0_, date_, pi_, pct_, period_ct_, period_dy_, group_h_, group_c_, group_w_, df_, dfs_, dfd_, df_l_, df_cl_, dfq_)
    t_ = t0_ - period_dy_;
    tt = [t_-i for i in 1:period_ct_]; 
    #println(" today = ", t0_, " - delay = ", t_, " : contact tracing 해야 하는 모든 날짜 = ", tt)
    
    # list_inf는 (t_-1-period_d ~ t_-period_cy-period_dy) 기간 동안 contact 한 인구들: list
    # contact tracing 할 사람들은 확률 p_에 의해 결정. (missing information: infector version) 
    list_inf = randsubseq(rngs_, df_l_[t_, :list_CT], pi_); 
    #println(" infector (in ct) = ", list_inf, " vs original ", df_l_[t_, :list_CT])
    list = Vector{Int64}()

    for n in 1:length(list_inf)
        id_ = list_inf[n]
        
        for i in 1:period_ct_
            num_ = findall(date_ .== tt[i])[1]

            # 그 날 격리였던 사람들 확인
            list_q = Vector{Int64}(); append!(list_q, df_l_[tt[i], :list_SQ], df_l_[tt[i], :list_Q], df_l_[tt[i], :list_IS]);
            
            if sum(findall(list_q .== id_)) == 0
                # 가정 내 접촉
                list_a = group_h_[df_[id_,:hh_id]][:, :id]
                tmp = randsubseq(rngs_, setdiff(list_a, list_q, id_), pct_[1])
                append!(list, tmp)
                
                # 평일: 학교/회사 
                if (tt[i]-1) % 7 + 1 <= 5
                    if df_[id_, :class] != 0
                        list_a = group_c_[df_[id_, :class]+1][:, :id]
                        tmp = randsubseq(rngs_, setdiff(list_a, list_q, id_), pct_[2])
                        append!(list, tmp)
                    elseif df_[id_, :office] != 0
                        list_a = group_w_[df_[id_, :office]+1][:, :id]
                        tmp = randsubseq(rngs_, setdiff(list_a, list_q, id_), pct_[3])
                        append!(list, tmp)
                    end
                end
    
                # 친구
                if dfd_[num_][id_, :num_f] != 0
                    tmp = randsubseq(rngs_, setdiff(df_cl_[num_][dfd_[num_][id_, :num_f], :list_f], id_), pct_[4])
                    append!(list, tmp)
                end
                
                # 우연
                if dfd_[num_][id_, :num_e] != 0
                    tmp = randsubseq(rngs_, setdiff(df_cl_[num_][dfd_[num_][id_, :num_e], :list_e], id_), pct_[5])
                    append!(list, tmp)
                end
            end
        end # for period_c
    end # for list_inf

    # id_ 가 period_c 기간 동안 만난 사람들 목록: list > Quarantine 상태로 넣어준다. 단, 이 인구들이 다른 state_q 상태면 안된다.
    unique!(list)
    ll = list[findall(dfs_[list, :state_q] .== 0)]

    dfs_[ll, :state_q] .= 2 # quarantine
    dfq_[ll, :count_q] .+= 1
    dfs_[ll, :period_Q] .= 1
    dfs_[ll, :count_q] .= 1

end
