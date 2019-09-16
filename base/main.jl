# SemiDefinite Programmin Optimal Power Flow (SDP-OPF)
# Tested in Julia 0.6.4.
#
# Alvaro González, Moscow, August 2018.
# alvaro.gonzalez@skolkovotech.ru

#-------------------------------- 0. SETUP ------------------------------------#
# 0.1. Importing packages
using Plots, LaTeXStrings, CSV, DataFrames, JuMP, Mosek, Gurobi, Missings
  pyplot()
  # Setting font for plots
  font = Plots.font("Times New Roman", 17)
    pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

# Including functions
include("base\\readMP.jl");
include("base\\runsdpopf.jl"); # Function with the SDP OPF formulation
include("base\\readElData.jl");

#Loading System Data based on MATPOWER-like files
BusData, BranchData, GenData, GenCostData = readMP(["BusData9" "BranchData9" "GenData9" "GenCostData9"], sep = ["\t" " " " " ","]);

# Processing netwrok data
c, Mk, Mlm, Nodes, Pd, Pmin, Pmax, Qd, Qmin, Qmax, Smax, Vmin, Vmax, Yk, Yk_over, Ylm, Ylm_over, Ωg, Ωlm =
   readElData(BusData, BranchData, GenData, gencostdata = GenCostData, Approx = "SDP");

# Decreasing the STR of the line 6-5
STR = Smax[2,8]
Bus_RESULTS, Objective_RESULTS, EigenRatio, SNew = Dict(), Array{Float64}(1:9), Array{Float64}(1:9), Array{Float64}(1:9);
  W = Dict()
  for i in 1:9
    SNew[i] = STR-(i-1)*25/100
    Smax[2,8] = SNew[i]
    Smax[8,2] = SNew[i]
    # running the SDP-OPF in the modified network
    Bus_RESULTS[i], Objective_RESULTS[i], EigenRatio[i] =
    runsdpopf(c, Mk, Mlm, Nodes, Pd, Pmin, Pmax, Qd, Qmin, Qmax, Smax, Vmin, Vmax, Yk, Yk_over, Ylm, Ylm_over, Ωg, Ωlm)
  end
scatter(collect((STR-(i-1)*25/100)*100 for i in 1:9), EigenRatio,label = "", ylabel = L"Eigen ratio: $\frac{ρ_2}{ρ_3}$", xlabel = L"$|S_{2,8}^{max}| [MVA]$")
  savefig("figs\\Ratio9")
  df9 = DataFrame(Smax_28 = SNew, SDP_Objective_Value = Objective_RESULTS, Eigen_Ratio = EigenRatio);
  CSV.write("results_case9_sdpopf.csv",df9)
  for i in 1:9
    CSV.write("results\\Bus_results_case9_Smax28@"*string(100*SNew[i])*"MVA.csv",Bus_RESULTS[i])
  end

# 3-BUS SYSTEM
BusData, BranchData, GenData, GenCostData = readMP(["BusData3T" "BranchData3T" "GenData3T" "GenCostData3T"], sep = ["\t" " " " " ","]);

  # Processing network data
  c, Mk, Mlm, Nodes, Pd, Pmin, Pmax, Qd, Qmin, Qmax, Smax, Vmin, Vmax, Yk, Yk_over, Ylm, Ylm_over, Ωg, Ωlm =
   readElData(BusData, BranchData, GenData, gencostdata = GenCostData, Approx = "SDP");

# Decreasing the STR of the line 6-5
STR3 = Smax[3,2]
Bus_RESULTS3, Objective_RESULTS3, EigenRatio3= Dict(), Array{Float64}(1:5), Array{Float64}(1:5);
  for i in 1:5
    Smax[3,2] = STR3-(i-1)*2.5/100
    Smax[2,3] = STR3-(i-1)*2.5/100
    # running the SDP-OPF in the modified network
    Bus_RESULTS3[i], Objective_RESULTS3[i], EigenRatio3[i] =
    runsdpopf(c, Mk, Mlm, Nodes, Pd, Pmin, Pmax, Qd, Qmin, Qmax, Smax, Vmin, Vmax, Yk, Yk_over, Ylm, Ylm_over, Ωg, Ωlm)
  end
  scatter(collect((STR3-(i-1)*2.5/100)*100 for i in 1:5), EigenRatio3, label = "", ylabel = L"Eigen ratio: $\frac{ρ_2}{ρ_3}$", xlabel = L"$|S^{max}_{2,3}|[MVA]$", yaxis = :log10)
    savefig("figs\\Ratio3")
    df3 = DataFrame(Smax_23 = collect((STR3-(i-1)*2.5/100)*100 for i in 1:5), SDP_Objective_Value = Objective_RESULTS3, Eigen_Ratio = EigenRatio3)
    CSV.write("results_case3_sdpopf.csv",df3)

    for i in 1:5
      CSV.write("results\\Bus_results_case3_Smax23@"*string((STR3-(i-1)*2.5/100)*100)*"MVA.csv",Bus_RESULTS3[i])
    end
