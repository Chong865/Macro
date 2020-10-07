using Printf
using Statistics
using Plots
using Roots
using QuadGK
using TypedTables
#function to consider
u1(x) = log(x)
u2(x) = sqrt(x)
u3(x,y) =  (x^(1-y))/(1-y)
global xaxis = range(0.05,2;length=1000) ;
logf = map(u1, xaxis)

function newton(x::Array,y::Array,z)
    m=length(x) #here m is the number of data points, not the degree # of the polynomial
    a=diff(x,y)
    sum=a[1]
    pr=1.0
    for j in 1:(m-1)
        pr=pr*(z-x[j])
        sum=sum+a[j+1]*pr #newton's formula
    end
    return sum
end



function diff(x::Array,y::Array) #coefficient
    m = length(x) #here m is the number of data points. #the degree of the polynomial is m-1 a=Array{Float64}(undef,m)
    a = zeros(m)
    for i in 1:m
        a[i]=y[i]
    end
    for j in 2:m
        for i in reverse(collect(j:m))
            #println(i)
            a[i]=(a[i]-a[i-1])/(x[i]-x[i-(j-1)])
        end
    end
    return(a)
end

function RMSE(logf, interp)
    sum = 0
    for i in 2:length(xaxis)
        sum = sum + (logf[i] .-interp[i])^2
    end
    sum = sqrt(sum ./length(xaxis))
    return sum
end

Error = zeros(4)
global i = 1
for n = (4,8,12,20)
    global i# Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)) ; # Collect makes it an array instead of a collection
    yi = map(u1,xi) # the corresponding y-coordinates
    # Interpolation
    interp=map(z->newton(xi,yi,z),xaxis)
    println(n)
    Error[i]= RMSE(logf, interp)
    global i += 1
    println(i)
    # Plot
    gr()
    plot(title="Interpolation n=$n - Newton Polynomial")
    plot!(xaxis,logf,linewidth=3,label = "Log Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(xaxis,interp,linewidth=3,label="Interpolation")
    plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
    #savefig("HW3_Log_Newton_n_$n")
    #println(n)


end

#1b

t = Table(name = ["Newton Polynomial"], n4 = [Error[1]], n8 = [Error[2]],n12 = [Error[3]],n20 = [Error[4]])

x_axis = 1:4
plot(x_axis, Error,
label = "RMSE", legend = :bottomright, widths = [10], title = "Accuracy vs grid size (Newton)")
savefig("HW3_Log_Newton_Acc")
