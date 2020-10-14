### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ 63202204-0e51-11eb-1a8f-f9dc0e3cad29
begin
	using PyCall # Call PyCall wrapper
	using PyPlot # Call Matplotlib for plotting purposes
	using LaTeXStrings # Call LaTeX formatting library
	using LinearAlgebra # Call algebra library for vector manipulation
end

# ╔═╡ 6492430e-0e50-11eb-3419-ff3b319119b7
md"# Amplitude evolution on a circular array of harmonic elements coupled by nearest neighbour interaction"

# ╔═╡ 708fa56e-0e55-11eb-07b6-690d171d4d3b
md"The evolution of the initial amplitude $\sum_{l} q_{l}(0)\ket{l}$ is given by the equation

$$\begin{equation}\boxed{q_{n}(t)=\frac{1}{N+1}\sum\limits_{m,l=0}^{N}\cos\left(2t
        \sqrt{k}\sin\left(\frac{\pi m}{N+1}\right)\right)\cos\left(\frac{2\pi m}{N+1}(n-l)\right)q_{l}(0)}
\end{equation}$$

where $N$ is the dimension of the array, and $k$ is the strength of the coupling, in this case $k_{l}\equiv k,\forall l$.
"

# ╔═╡ a65d2a08-0e51-11eb-0085-49db35b81bcb
rc("text", usetex=true) # TeX formatting for the labels and titles

# ╔═╡ ad1523a0-0e51-11eb-18cf-ddb9466c94af
begin
	N = 32 # Coupling's dimension
	dim = N+1 # Coupled element's dimension
	k = 1 # Strenght of the coupling
end;

# ╔═╡ c0ad476c-0e51-11eb-2f41-b5a6a61b49be
function F1(t,m,n,l) # Define the analytic function of the evolution, thesis equation (1.13) 
    return cos(2*t*sqrt(k)*sin(pi*m/(N+1)))*cos((2*pi*m/(N+1))*(n-l))
end;

# ╔═╡ cbaa2d5e-0e51-11eb-3fde-71c180be815c
function pint(a) # Initial condition as a vector with an amplitude at position 0<=a<=N+1
    pv = zeros(Float64,dim)
    for j in 0:N
        if a == j
            pv[j+1] = 1
        else
            pv[j+1] = 0
        end
    end
    return pv
end;

# ╔═╡ 39e9cff6-0e52-11eb-3bd8-499b3cc3274e
begin
	icond1 = normalize(pint(0)) # Initial condition one

	function q1(t,n)
    	mevol = zeros(Float64,(dim,dim))
    	for m in 0:N
        	for l in 0:N
            	mevol[m+1,l+1] = F1(t,m,n,l)*icond1[l+1]
        	end
    	end
    	return (1/(N+1))*sum(mevol)
	end
end;

# ╔═╡ 666e31fc-0e52-11eb-2cdc-afde9f0b8a63
begin
	icond2 = normalize(pint(8)+pint(24)) # Initial condition two

	function q2(t,n)
    	mevol = zeros(Float64,(dim,dim))
    	for m in 0:N
        	for l in 0:N
            	mevol[m+1,l+1] = F1(t,m,n,l)*icond2[l+1]
        	end
    	end
    	return (1/(N+1))*sum(mevol)
	end
end;

# ╔═╡ 7fcb5e98-0e52-11eb-108c-156668a650c7
begin
	icond3 = normalize(pint(16)) # Initial condition three

	function q3(t,n)
    	mevol = zeros(Float64,(dim,dim))
    	for m in 0:N
        	for l in 0:N
            	mevol[m+1,l+1] = F1(t,m,n,l)*icond3[l+1]
        	end
    	end
    	return (1/(N+1))*sum(mevol)
	end
end;

# ╔═╡ 9292a5d8-0e52-11eb-3281-bfa3f42a0fc1
t_list = range(0,stop=61,length=61); # Set the vector list of the time sampling

# ╔═╡ 9bc82306-0e52-11eb-286b-ddeff0f34b8b
begin
	evolt1 = zeros(Float64,(dim,length(t_list))) # Time evolution of the initial condition one
	for t in 1:length(t_list)
    	for n in 0:N    
       		evolt1[n+1,t] = q1(t_list[t],n)
    	end
	end
end

# ╔═╡ d1ad92c8-0e52-11eb-26ca-0df20f1729c6
begin
	evolt2 = zeros(Float64,(dim,length(t_list))) # Time evolution of the initial condition two
	for t in 1:length(t_list)
    	for n in 0:N    
       		evolt2[n+1,t] = q2(t_list[t],n)
    	end
	end
end

# ╔═╡ e051dfd2-0e52-11eb-2b8a-67a9754c30fd
begin
	evolt3 = zeros(Float64,(N+1,length(t_list))) # Time evolution of the initial condition three
	for t in 1:length(t_list)
    	for n in 0:N    
       		evolt3[n+1,t] = q3(t_list[t],n)
    	end
	end
end

# ╔═╡ f9bb3afe-0e52-11eb-21de-5d468a99b5c1
begin
	vminh = -1 #minimum((minimum(evolt1),minimum(evolt2),minimum(evolt3))) # Normalized minima of the three evolutions
	vmaxh = maximum((maximum(evolt1),maximum(evolt2),maximum(evolt3))) # Normalized maxima of the three evolutions
end;

# ╔═╡ 103c4e8a-0e53-11eb-1b3e-11ef86f9e9b4
begin
# FIGURE PLOT

cm = "seismic" # Diverging colormap for plotting negative and positive values

ww = 6.20 # Width size of the figure
hh = ww # Height size of the figure

fig,(ax1,ax2,ax3)=plt.subplots(3,1,figsize=(ww,hh),sharex=true)
plt.subplots_adjust(hspace = 0.1)

ax1.imshow(evolt1,vmin=vminh,vmax=vmaxh,cmap=cm)
ax2.imshow(evolt2,vmin=vminh,vmax=vmaxh,cmap=cm)
ax3.imshow(evolt3,vmin=vminh,vmax=vmaxh,cmap=cm)

ax1.set_aspect("auto")
ax2.set_aspect("auto")
ax3.set_aspect("auto")

ax1.tick_params(direction="out",length=5,width=1,labelsize=10)
ax2.tick_params(direction="out",length=5,width=1,labelsize=10)
ax3.tick_params(direction="out",length=5,width=1,labelsize=10)

ax1.set_yticks(0:8:N, minor = false)
ax2.set_yticks(0:8:N, minor = false)
ax3.set_yticks(0:8:N, minor = false)

ax1.set_yticks(1:1:N, minor = true)
ax1.grid(which = "minor", color = "gray", linestyle = ":", linewidth = 0.5, alpha = 0.25)
ax2.set_yticks(1:1:N, minor = true)
ax2.grid(which = "minor", color = "gray", linestyle = ":", linewidth = 0.5, alpha = 0.25)
ax3.set_yticks(1:1:N, minor = true)
ax3.grid(which = "minor", color = "gray", linestyle = ":", linewidth = 0.5, alpha = 0.25)

ax1.set_ylabel(L"q_{n}(t)",fontsize=10)
ax2.set_ylabel(L"q_{n}(t)",fontsize=10)
ax3.set_ylabel(L"q_{n}(t)",fontsize=10)

ax3.set_xlabel(L"t\ [s]",fontsize=10)

pcm = ax1.get_children()[10]
cb = colorbar(pcm,ax=(ax1,ax2,ax3),extend="both",ticks=[-1,-0.5,0,0.5,1],orientation="vertical",shrink=0.7,aspect=35,fraction=0.015)
cb.ax.tick_params(labelsize=10,length=5,width=1,direction="inout")
#cb.ax.set_ylabel("Amplitud",fontsize=12,labelpad=0)

tight_layout(rect=(0, 0, 0.9, 1))

savefig("circ_vd.pdf", transparent = "true", dpi=300, bbox_inches="tight", pad_inches=0)
end

# ╔═╡ 76160e14-0e53-11eb-346b-eb6cddfa034e
print("Julia:"," ", VERSION, "\nPyCall:"," ", PyCall.pyversion, "\nPyPlot:"," ", PyPlot.version, "\nOS Kernel:"," ", Sys.KERNEL, "\nArchitecture:"," ", Sys.ARCH)

# ╔═╡ Cell order:
# ╟─6492430e-0e50-11eb-3419-ff3b319119b7
# ╟─708fa56e-0e55-11eb-07b6-690d171d4d3b
# ╠═63202204-0e51-11eb-1a8f-f9dc0e3cad29
# ╠═a65d2a08-0e51-11eb-0085-49db35b81bcb
# ╠═ad1523a0-0e51-11eb-18cf-ddb9466c94af
# ╠═c0ad476c-0e51-11eb-2f41-b5a6a61b49be
# ╠═cbaa2d5e-0e51-11eb-3fde-71c180be815c
# ╠═39e9cff6-0e52-11eb-3bd8-499b3cc3274e
# ╠═666e31fc-0e52-11eb-2cdc-afde9f0b8a63
# ╠═7fcb5e98-0e52-11eb-108c-156668a650c7
# ╠═9292a5d8-0e52-11eb-3281-bfa3f42a0fc1
# ╠═9bc82306-0e52-11eb-286b-ddeff0f34b8b
# ╠═d1ad92c8-0e52-11eb-26ca-0df20f1729c6
# ╠═e051dfd2-0e52-11eb-2b8a-67a9754c30fd
# ╠═f9bb3afe-0e52-11eb-21de-5d468a99b5c1
# ╠═103c4e8a-0e53-11eb-1b3e-11ef86f9e9b4
# ╠═76160e14-0e53-11eb-346b-eb6cddfa034e
