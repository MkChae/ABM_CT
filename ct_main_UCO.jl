
filepath = pwd() * "/"
include(filepath * "ct_ftn.jl")

function Oneset(thread_size)

    Random.seed!(1)

    tf = 3000   # end of simulation
    dt = 1     # dt = 1 means we check the infection end of days.
    ntries = Int64(tf/dt) 

    println(" read csv file")
    df = Read_CSV_File()
    N = nrow(df)

    group_h = groupby(df, [:hh_id])  
    group_c = groupby(df, [:class])
    group_w = groupby(df, [:office])

    println(" initial parameter seeting ")
    period_ct = 7                         # period of contact tracing policy
    period_is = 7                         # period of isolation
    period_dy = 0                         # period of delay
    period_rd = period_ct + period_dy + 1 # period of record

    Pf = [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7] # 요일별 친구 만날 확률
    Pe = [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7] # 요일별 우연히 누군가를 만날 확률

    n0 = 20 # initial infector

    Pa = 0.2  # probability of asymptomatic 
    Psq = 0.5 # probability of test and self-quarantine after symptoms 

    # missing information
    number_control_infection = 10000;
    Pmi    = 1.0;                           # missing infector                             
    Pmi_or = 0.5;                           # 확진자가 너무 많아지면, missing infector 비율이 더 커짐
    Pmc    = [0.90, 0.90, 0.90, 0.90, 0.50];    # missing information (contact people missing)
    Pmc_or = [0.5, 0.5, 0.5, 0.5, 0.5];     # 확진자가 너무 많아지면, missing contact people 비율이 더 커짐

    vs = Viral_Shedding()

    println(" make dataframe : node num")
    # agents data
    dfs = [Make_df_State(N) for i in 1:thread_size]
    dfq = [Make_df_Qurantine(N) for i in 1:thread_size]
    dfv = [Make_df_Viral(N) for i in 1:thread_size]
    dft = [Make_df_Daily(N) for i in 1:thread_size]
    dfd = [Make_df_Daily_Record(N, period_rd) for i in 1:thread_size]

    # daily data
    df_l = [Make_df_List(ntries) for i in 1:thread_size]
    df_lc = [Make_df_ListCount(ntries) for i in 1:thread_size]
    df_cl = [Make_df_Contact_List(N, period_rd) for i in 1:thread_size]

    ncycles = 20; nt = Int64(thread_size * ncycles);
    rngs = [MersenneTwister(i) for i in 1:thread_size]; 
    println(" code running : total cycle = ", ncycles, " final time = ", ntries+Int64(period_ct/dt))

    Threads.@threads for rank = 1:thread_size
        for n in 1:thread_size:nt

            fn = n + rank - 1
            rs = Int64(round(time() + fn*(fn+1))) # tmp_rs[rank]  #rank #
            Random.seed!(rngs[rank], rs)
            println(" rank = ", rank, " fn = ", fn, " / ", nt ," random seed = ", rs)

            date = zeros(Int64, period_rd); # 과거 기록이 필요해서 몇일 인지 기록해 두는 것. 
            InitSetting_Df(dfs[rank], dfv[rank], dft[rank], dfd[rank], df_l[rank], df_cl[rank], period_rd)

            start_time = time()
            # 감염병 시작하기 전 일상 (기록을 해둬야지 나중에 ct 할 수 있음.)
            for t in 1:Int64(period_rd/dt)
                num = (t-1) % period_rd + 1
                day = (t-1) % 7 + 1
                date[num] = t

                Friends_Meeting(rngs[rank], num, Pf[day], df, dfs[rank], dft[rank], dfd[rank], df_cl[rank])    # 친구 만남
                Random_Encounter(rngs[rank], num, Pe[day], df, dfs[rank], dft[rank], dfd[rank], df_cl[rank])   # 우연한 만남 
            end
            end_time = time()
            println(" rank = ", rank, " : fn = ", fn ," ", period_rd, "을 계산하는데 총 소요 시간 = ", round(end_time-start_time, digits=3), " s")
        
            Initial_Infectors_Rand(rngs[rank], Int64(period_rd/dt), n0, dfs[rank], dfv[rank])

            start_time = time()
            # 감염병 시작
            for t in Int64(period_rd/dt)+1:ntries

                num = (t-1) % period_rd + 1
                day = (t-1) % 7 + 1
                date[num] = t

                df_l[rank][t, :list_SQ] = dfs[rank][dfs[rank].state_q .== 1, :id]; df_lc[rank][t, :list_SQ] = length(df_l[rank][t, :list_SQ])
                df_l[rank][t, :list_Q] = dfs[rank][dfs[rank].state_q .== 2, :id]; df_lc[rank][t, :list_Q] = length(df_l[rank][t, :list_Q])
                df_l[rank][t, :list_IS] = dfs[rank][dfs[rank].state_q .== 3, :id]; df_lc[rank][t, :list_IS] = length(df_l[rank][t, :list_IS])

                # 남을 감염 시킬 수 있는 인구: 바리어스 배출 O + 격리 X
                list_i = findall((dfv[rank].state .== 1) .& (dfs[rank].state_q .== 0))
                # 격리중 인구
                list_q = Vector{Int64}(); append!(list_q, df_l[rank][t, :list_SQ], df_l[rank][t, :list_Q], df_l[rank][t, :list_IS])
                # SEIE : state: 0 - susceptible / 1 - exposed / 2 - asymptomatic / 3 - symptomatic / 4 - R
                list_s = dfs[rank][(1 .<= dfs[rank].state .<= 3), :id]

                Friends_Meeting(rngs[rank], num, Pf[day], df, dfs[rank], dft[rank], dfd[rank], df_cl[rank])    # 친구 만남
                Random_Encounter(rngs[rank], num, Pe[day], df, dfs[rank], dft[rank], dfd[rank], df_cl[rank])   # 우연한 만남 

                [Check_Viral_State(list_s[i], dfv[rank]) for i in 1:length(list_s)]
                [Check_Disease_State(rngs[rank], list_s[i], t, Pa, Psq, dfs[rank]) for i in 1:length(list_s)]
                list_c = Vector{Int64}(); [append!(list_c, Check_Confinement_State(rngs[rank], list_q[i], t, period_is, dfs[rank])) for i in 1:length(list_q)]
                [Compute_Lambda(rngs[rank], list_i[i], t, day, num, list_q, vs, group_h, group_c, group_w, df, dfs[rank], dfv[rank], dfd[rank], df_cl[rank]) for i in 1:length(list_i)]
                
                #=
                if length(list_s) > number_control_infection
                    pif = Pmi_or;   pct = Pmc_or;
                else
                    pif = Pmi;        pct = Pmc;
                end
                =#
                pif = Pmi;        
                pct = Pmc;

                # 오늘 positive 가 나와서, contact tracing 을 할 사람
                df_l[rank][t+1, :list_CT] = unique!(list_c); df_lc[rank][t+1, :list_CT] = length(df_l[rank][t+1, :list_CT])
                Contact_Tracing(rngs[rank], t, date, pif, pct, period_ct, period_dy, group_h, group_c, group_w, df, dfs[rank], dfd[rank], df_l[rank], df_cl[rank], dfq[rank])

                if (length(list_s) == 0) .& (length(list_q) == 0) .& (length(list_i) == 0)
                    break;
                end
            end
            end_time = time()
            println(" rank = ", rank, " : fn = ", fn , " ", ntries, "일 계산하는데 총 소요 시간 = ", round(end_time-start_time, digits=3), " s")    
            
            CSV.write("./output/IBM_CT_rs$(rs)_state.csv", dfs[rank][dfs[rank].state .!= 0, :]; bufsize=60_000_000)
            CSV.write("./output/IBM_CT_rs$(rs)_qlist.csv", df_lc[rank]; bufsize=60_000_000)
            CSV.write("./output/IBM_CT_rs$(rs)_qurantine.csv", dfq[rank][dfq[rank].count_q .!= 0, :]; bufsize=60_000_000)
        end
    end

    df = nothing;
    group_h = nothing; group_c = nothing; group_w = nothing;

    dfs = nothing; dfv = nothing; dfd = nothing; dft = nothing; dfq = nothing;
    df_cl = nothing; df_l = nothing; rngs = nothing;
    GC.gc()

end



function main()
    Random.seed!(1)
    thread_size = Threads.nthreads()
    println(" contact tracing : threads # = ", thread_size)
    Oneset(thread_size)     
    println(" All Complete!")
end

main()