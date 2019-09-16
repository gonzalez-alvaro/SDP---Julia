function readElData(BusData, BranchData, GenData; Approx = "DC", gencostdata = [], Sbase = 100)
    SystemData = Dict()

    N = size(BusData,1) # Number of nodes
    Nodes = 1:N
    L = size(BranchData,1) # Number of lines

    # Processing Generation Data
    NGen = size(GenData,1)

    # Location of Gen
    Ωg = zeros(N,1);
    for g in 1:NGen
        Ωg[GenData[g,:bus]] = 1;
    end

    # Generation Limits
    Pmin, Pmax, Qmin, Qmax = zeros(N,1), zeros(N,1), zeros(N,1), zeros(N,1);
    for g in 1:NGen
        Pmin[GenData[g,:bus]], Pmax[GenData[g,:bus]], Qmin[GenData[g,:bus]], Qmax[GenData[g,:bus]] =
        GenData[g, :Pmin], GenData[g, :Pmax],
        GenData[g, :Qmin], GenData[g, :Qmax]
    end

    # Cost coefficients
    c = zeros(N,3);
    for g in 1:NGen
        c[GenData[g,:bus],:] = [float(gencostdata[g,:c0]), Sbase*float(gencostdata[g,:c1]), (Sbase^2)*float(gencostdata[g,:c2])]
    end
    # Nodal information: Voltage limits and loads
    Vmin, Vmax, Pd, Qd = zeros(N,1), zeros(N,1), zeros(N,1), zeros(N,1);
    for k in 1:N
        Vmin[k], Vmax[k], Pd[k], Qd[k] = float(BusData[k,:Vmin]), float(BusData[k,:Vmax]), float(BusData[k,:Pd]), float(BusData[k,:Qd]);
    end


    # Maximum apparent power through line
    Smax = zeros(N,N);
    # Existance of lines
    Ωlm = zeros(N,N);
    for l in 1:L
        Smax[BranchData[l,:fbus],BranchData[l,:tbus]], Smax[BranchData[l,:tbus],BranchData[l,:fbus]]  =
        BranchData[l,:rateA],BranchData[l,:rateA]
        Ωlm[BranchData[l,:fbus],BranchData[l,:tbus]], Ωlm[BranchData[l,:tbus],BranchData[l,:fbus]]  =
        1,1
    end

    # Branch Impedance and Admittance Parameters
    BranchData[:Z] = BranchData[:r] + im*BranchData[:x];
    BranchData[:Y] = 1 ./ BranchData[:Z];

    if Approx == "SDP"
        # Creating Admittance Matrix - Y
        Y = zeros(Complex,N,N);
        # Line information
        for l in 1:L
            Y[BranchData[l,:fbus],BranchData[l,:tbus]] = -BranchData[l,:Y]
            Y[BranchData[l,:tbus],BranchData[l,:fbus]] = -BranchData[l,:Y]
        end
        # Diagonal elements
        for k = 1:N, l =1:L
            if (k==BranchData[l,:fbus])||(k==BranchData[l,:tbus])
                # Admittance Matrix
                Y[k,k] += BranchData[l, :Y] + 0.5*im*BranchData[l,:b]
            end
        end

        # Shunt element
        ylm_over = zeros(Complex,N,N);
        for l in 1:L
            ylm_over[BranchData[l,:fbus],BranchData[l,:tbus]] = 0.5*im*BranchData[l,:b]
            ylm_over[BranchData[l,:tbus],BranchData[l,:fbus]] = 0.5*im*BranchData[l,:b]
        end

        ek = eye(N);
        Yk_pre = Dict() #Dict(N)[1:N,1:N]
        Yk = Dict() #Dict(N)[1:2N,1:2N]
        Yk_over = Dict() #Dict(N)[1:2N,1:2N]
        for k in 1:N
            # Calculating Yk
            Yk_pre[k] = ek[:,k]*ek[:,k]'*Y

            # Calculating *bold* Yk
            Yk[k] = 0.5*[real.(Yk_pre[k]+Yk_pre[k].') imag.(Yk_pre[k].'-Yk_pre[k]);
            imag.(Yk_pre[k]-Yk_pre[k].') real.(Yk_pre[k]+Yk_pre[k].')]

            # Calulcating ̅Yk
            Yk_over[k] = -0.5*[imag.(Yk_pre[k]+Yk_pre[k].') real.(Yk_pre[k]-Yk_pre[k].');
            real.(Yk_pre[k].'-Yk_pre[k]) imag.(Yk_pre[k]+Yk_pre[k].')]
        end

        # Shunt admittance matrices
        Ylm_pre = Dict() #Dict(N,N)[1:N,1:N]
        Ylm = Dict() #Dict(N,N)[1:2N,1:2N]
        Ylm_over = Dict() #Dict(N,N)[1:2N,1:2N]
        for l in 1:N, m in 1:N
            if Ωlm[l,m] == 1
                # Calculating Ylm
                Ylm_pre[l,m] = (ylm_over[l,m]-Y[l,m])*ek[:,l]*ek[:,l]'+Y[l,m]*ek[:,l]*ek[:,m]'

                # Calculating *bold* Ylm
                Ylm[l,m] =
                0.5*[real.(Ylm_pre[l,m]+Ylm_pre[l,m].') imag.(Ylm_pre[l,m].'-Ylm_pre[l,m]);
                imag.(Ylm_pre[l,m]-Ylm_pre[l,m].') real.(Ylm_pre[l,m]+Ylm_pre[l,m].')]

                # Calulcating ̅Ylm
                Ylm_over[l,m] =
                -0.5*[imag.(Ylm_pre[l,m]+Ylm_pre[l,m].') real.(Ylm_pre[l,m]-Ylm_pre[l,m].');
                real.(Ylm_pre[l,m].'-Ylm_pre[l,m]) imag.(Ylm_pre[l,m]+Ylm_pre[l,m].')]
            end
        end

        Mk = Dict() #Dict(N)(2*N,2*N)
        for k = 1:N
            Mk[k] = [ek[:,k]*ek[:,k]' zeros(N,N)
            zeros(N,N) ek[:,k]*ek[:,k]']
        end

        Mlm = Dict() #Dict(N,N)(2*N,2*N)
        for l = 1:N, m in 1:N
            if Ωlm[l,m] == 1
                Mlm[l,m] = [(ek[:,l]-ek[:,m])*(ek[:,l]-ek[:,m])' zeros(N,N)
                zeros(N,N) (ek[:,l]-ek[:,m])*(ek[:,l]-ek[:,m])']
            end
        end

        return c, Mk, Mlm, Nodes, Pd./Sbase, Pmin./Sbase, Pmax./Sbase, Qd./Sbase, Qmin./Sbase, Qmax./Sbase, Smax./Sbase, Vmin, Vmax, Yk, Yk_over, Ylm, Ylm_over, Ωg, Ωlm
        end

end
