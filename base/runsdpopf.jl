function runsdpopf(c, Mk, Mlm, Nodes, Pd, Pmin, Pmax, Qd, Qmin, Qmax, Smax, Vmin, Vmax, Yk, Yk_over, Ylm, Ylm_over, Ωg, Ωlm)

# Optimization model
sdpOPF = Model(solver = MosekSolver(MSK_IPAR_LOG=0))

  #PARAMETERS
  a=zeros(Nodes[end],1)
    b=zeros(Nodes[end],1)
    for k in collect(k for k in Nodes if Ωg[k]==1)
      a[k]=c[k,1] + c[k,2]*Pd[k]
      b[k]=sqrt(c[k,3])*Pd[k]
    end

  #VARIABLES
  @variable(sdpOPF, W[1:2*Nodes[end],1:2*Nodes[end]], SDP)
  @variable(sdpOPF, α[k in Nodes; Ωg[k]==1])

  #CONSTRAINTS
  @constraint(sdpOPF, PConstMin[k in Nodes],
      Pmin[k] - Pd[k] ≤ trace(Yk[k]*W))
  @constraint(sdpOPF, PConstMax[k in Nodes],
      Pmax[k] - Pd[k] ≥ trace(Yk[k]*W))
  @constraint(sdpOPF, QConstMin[k in Nodes],
      Qmin[k] - Qd[k] ≤ trace(Yk_over[k]*W))
  @constraint(sdpOPF, QConstMax[k in Nodes],
      Qmax[k] - Qd[k] ≥ trace(Yk_over[k]*W))
  @constraint(sdpOPF, VConstMin[k in Nodes],
      (Vmin[k])^2 ≤ trace(Mk[k]*W))
  @constraint(sdpOPF, VConstMax[k in Nodes],
      (Vmax[k])^2 ≥ trace(Mk[k]*W))
  # Limits on Slm
  for L in collect([l,m] for l in Nodes, m in Nodes if Ωlm[l,m]==1)
    @SDconstraint(sdpOPF,
      [-(Smax[L[1],L[2]])^2           trace(Ylm[L[1],L[2]]*W) trace(Ylm_over[L[1],L[2]]*W);
      trace(Ylm[L[1],L[2]]*W)       -1                0;
      trace(Ylm_over[L[1],L[2]]*W)  0                 -1                    ]   ≤ 0)
  end
  # Cost function limits
  for k in collect(k for k in Nodes if Ωg[k]==1)
    @SDconstraint(sdpOPF,
      [c[k,2]*trace(Yk[k]*W)-α[k]+a[k]           sqrt(c[k,3])*trace(Yk[k]*W)+b[k];
      sqrt(c[k,3])*trace(Yk[k]*W)+b[k]          -1                             ]   ≤ 0)
  end

  #OBJECTIVE
  @objective(sdpOPF, Min, sum(α[k] for k in Nodes if Ωg[k]==1))

  status = solve(sdpOPF);
  Wopt =getvalue(W)
  eigenvals = sort(eigvals(Wopt), rev=true)
  eigenvects = eigvecs(Wopt)
  EigenRatio = eigenvals[2]/eigenvals[3]

    # RESULTS AT THE BUSES
    Xopt = sqrt(eigenvals[1])*eigenvects[:,end] + sqrt(eigenvals[2])*eigenvects[:,end-1]
    Vopt = Xopt[1:trunc(Int,length(Xopt)/2)] + im* Xopt[trunc(Int,length(Xopt)/2)+1:end]
    PGopt = Array{Float64}(Nodes[end])
    QGopt = Array{Float64}(Nodes[end])
    for k in collect(k for k in Nodes if Ωg[k]==1)
      PGopt[k] = 100*(trace(Yk[k]*Wopt)+Pd[k])
      QGopt[k] = 100*(trace(Yk_over[k]*Wopt)+Qd[k])
    end

    Bus_RESULTS =DataFrame(Bus = Nodes,
      VMag_pu = round.(abs.(Vopt),3),
      VAng_deg = round.(rad2deg.(angle.(Vopt)-angle(Vopt[1])),3),
      Pg_MW = round.(PGopt,2),
      Qg_MVAr = round.(QGopt,2),
      Pd_MW = 100*Pd[:],
      Qd_MVAr = 100*Qd[:])

    # RESULTS FOR THE BRANCHES
    Objective_RESULTS = sum(c[k,2+1]*(trace(Yk[k]*Wopt)+Pd[k])^2
                            + c[k,1+1]*(trace(Yk[k]*Wopt)+Pd[k])
                            + c[k,0+1] for k in Nodes if Ωg[k]==1)

  return Bus_RESULTS, Objective_RESULTS, EigenRatio
end
