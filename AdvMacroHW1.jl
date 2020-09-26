using Printf
using Statistics
using Plots
run = 19
a = 1/3
b = 0.8
z = 1
function parameter(a,b,v)
    global k_star = (a*b*v)^(1/(1-a))
    global A = a/(1-a*b)
    global k_prime = x -> (b*A*v*x^a)/(1+b*A)
    global c = (x,y) -> v*x^a-y
    global r = x -> a*v*x^(a-1)
    global w = x -> (1-a)*v*x^a
    global y = x -> v*x^a
end
parameter(a,b,z)
ss_k = k_star
k = k_star*0.8
data = zeros(run+1, 5)
ss_data_a = [ones(run+1,1)*k_star ones(run+1,1)*r(k_star) (
ones(run+1,1)*w(k_star)) ones(run+1,1)*y(k_star) (
ones(run+1,1)*c(k_star, k_star))]
data[1,1] = k
#println(k)
for  i in 1:run
    #r
    global data[i,2] = r(k)
    #w
    global data[i,3] = w(k)
    #y
    global data[i,4] = y(k)
    global k = k_prime(k)
    #capital
    global data[i+1,1] = k
    #c
    global data[i,5] = c(data[i,1],data[i+1,1])
    if i== run
        global data[i+1,2] =r(k)
        global data[i+1,3] =w(k)
        global data[i+1,4] =y(k)
    else
    end
    println(k)
end
x_axis = 1:10
println("--------------------------")
plot(x_axis, data[1:10,1:5], layout = (3, 2),
label = ["k" "r" "w" "y" "c"], legend = :bottomright)

plot!(x_axis, ss_data_a[1:10,1:5], layout = (3, 2),
label = ["ss k" "ss r" "ss w" "ss y" "ss c"], linestyle = :dot)
xlabel!("Iteration")
#xt = 1:20
#yt = rand(20, 5)
#plot(xt, yt, layout = (5, 1))
#plot(xt, data, layout = (5, 1))
savefig("HW17a")

parameter(a,b,1.05)

data_b = zeros(run+1, 5)
ss_data_b = [ones(run+1,1)*k_star ones(run+1,1)*r(k_star) (
ones(run+1,1)*w(k_star)) ones(run+1,1)*y(k_star) (
ones(run+1,1)*c(k_star, k_star))]
data_b[1,1] = ss_k
println(ss_k)
for  i in 1:run
    #r
    global data_b[i,2] = r(ss_k)
    #w
    global data_b[i,3] = w(ss_k)
    #y
    global data_b[i,4] = y(ss_k)
    global ss_k = k_prime(ss_k)
    #capital
    global data_b[i+1,1] = ss_k
    #c
    global data_b[i,5] = c(data_b[i,1],data_b[i+1,1])
    if i== run
        global data_b[i+1,2] =r(ss_k)
        global data_b[i+1,3] =w(ss_k)
        global data_b[i+1,4] =y(ss_k)
    else
    end
    println(ss_k)
end
x_axis = 1:10
println("--------------------------")
plot(x_axis, data_b[1:10,1:5], layout = (3, 2),
label = ["k" "r" "w" "y" "c"], legend = :bottomright)

plot!(x_axis, ss_data_b[1:10,1:5], layout = (3, 2),
label = ["ss k" "ss r" "ss w" "ss y" "ss c"], linestyle = :dot)
xlabel!("Iteration")
savefig("HW17b")
