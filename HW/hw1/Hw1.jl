
using PyPlot;

# Problem 6

# Solve the IVP
# 
# u′(t) = u(t)² - sin(t) + cos²(t),
# u(0) = 1
# 
# numerically using (a) Forward Euler, and (b) 2-stage RK method (5.32).


# u′ = f(u,t)
function f(u,t)
    return u^2 - sin(t) - cos(t)^2; 
end

# Exact solution (for plotting)
uₑ = t -> cos(t);

# Parameters:
T = 8;         # Final time
u₀ = 1;         # Initial guess (in this case the initial condition)
NVals = 25 * (2.^linspace(0,7,8));      # Number of grid points to be used in each approximation


# (a) Forward Euler

FE_ErrorVector = zeros(8);    # Vector in which to store error at t=T for each choice of Δt


for k = 1:8
    N = convert(Int64,NVals[k]);     # Number of grid points
    Δt = T / N;       # Timestep size
    u = zeros(N);     # Vector in which to store values of solution as time is advanced

    # Iterate from t=0 to t=T in increments of Δt
    u[1] = u₀ + Δt * f(u₀,0);
    for j=2:N
        u[j] = u[j-1] + Δt * f(u[j-1],(j-1) * Δt);
    end
    
    # Add entry to vector keeping track of error at t=8
    FE_ErrorVector[k] = uₑ(T) - u[end];

    # Plot one of the approximations against the true solution
    if N==50
        PyPlot.figure()
        t = linspace(0,T,10000)
        PyPlot.plot(t,uₑ(t),"k")
        PyPlot.plot(linspace(0,T,N+1),[u₀; u],"bo")
        PyPlot.legend(["Exact Solution","Forward Euler Approximation"]);
        PyPlot.xlabel("time")
        PyPlot.ylabel("u(x)")
    end
end

# Plot Error for FE method
PyPlot.figure()
PyPlot.loglog(T ./ NVals,abs(FE_ErrorVector),"bo-")
PyPlot.xlabel("Time Step")
PyPlot.ylabel("Error at t = 8")



# (b) 2-stage Runge-Kutta

RK_ErrorVector = zeros(8);    # Vector in which to store error at t=T for each choice of Δt

for k = 1:8
    N = convert(Int64,NVals[k]);     # Number of grid points
    Δt = T / N;       # Timestep size
    u = zeros(N);     # Vector in which to store values of solution as time is advanced

    # Iterate from t=0 to t=T in increments of Δt
    ū = u₀ + (1/2) * Δt * f(u₀,0);
    u[1] = u₀ + Δt * f(ū, Δt / 2);
    for j=2:N
        ū = u[j-1] + (1/2) * Δt * f(u[j-1], (j-1) * Δt);
        u[j] = u[j-1] + Δt * f(ū,(j-1) * Δt + Δt/2);
    end
    
    # Add entry to vector keeping track of error at t=10
    RK_ErrorVector[k] = uₑ(T) - u[end];
    
    # Plot one of the approximations against the true solution
    if N==50
        t = linspace(0,T,10000)
        PyPlot.figure()
        PyPlot.plot(t,uₑ(t),"k")
        PyPlot.plot(linspace(0,T,N+1),[u₀; u],"bo")
        PyPlot.legend(["Exact Solution","Runge-Kutta Approximation"]);
        PyPlot.xlabel("time")
        PyPlot.ylabel("u(x)")
    end
end

# Plot Error for RK method
PyPlot.figure()
PyPlot.loglog(T ./ NVals,abs(RK_ErrorVector),"bo-")
PyPlot.xlabel("Time Step")
PyPlot.ylabel("Error at t = 8");


