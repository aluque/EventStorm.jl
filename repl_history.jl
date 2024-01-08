# time: 2023-12-20 15:10:39 CET
# mode: pkg
	activate .
# time: 2023-12-20 15:11:10 CET
# mode: pkg
	add DifferentialEquations StaticArrays CSV DataFrames Interpolations LinearAlgebra Polyester Printf
# time: 2023-12-20 15:17:47 CET
# mode: julia
	using StaticArrays
# time: 2023-12-20 15:17:47 CET
# mode: julia
	using CSV
# time: 2023-12-20 15:17:48 CET
# mode: julia
	using DataFrames
# time: 2023-12-20 15:17:49 CET
# mode: julia
	using Interpolations
# time: 2023-12-20 15:17:49 CET
# mode: julia
	using LinearAlgebra
# time: 2023-12-20 15:17:50 CET
# mode: julia
	using Polyester
# time: 2023-12-20 15:17:50 CET
# mode: julia
	using DifferentialEquations
# time: 2023-12-20 15:18:01 CET
# mode: julia
	using Printf
# time: 2023-12-20 15:18:17 CET
# mode: julia
	using StaticArrays
# time: 2023-12-20 15:18:17 CET
# mode: julia
	using CSV
# time: 2023-12-20 15:18:17 CET
# mode: julia
	using DataFrames
# time: 2023-12-20 15:18:18 CET
# mode: julia
	using Interpolations
# time: 2023-12-20 15:18:18 CET
# mode: julia
	using LinearAlgebra
# time: 2023-12-20 15:18:18 CET
# mode: julia
	using Polyester
# time: 2023-12-20 15:18:19 CET
# mode: julia
	using DifferentialEquations
# time: 2023-12-20 15:18:29 CET
# mode: julia
	using Printf
# time: 2023-12-20 15:18:41 CET
# mode: pkg
	dev /Users/luque/projects/julia/Constants/
# time: 2023-12-20 15:19:02 CET
# mode: pkg
	dev /Users/luque/projects/dongshuai/dipoles/DipoleRadiators/
# time: 2023-12-20 15:19:11 CET
# mode: pkg
	dev /Users/luque/projects/Chemise.jl
# time: 2023-12-20 15:43:32 CET
# mode: julia
	includet("event_storm.jl")
# time: 2023-12-20 15:45:48 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 15:46:00 CET
# mode: julia
	includet("event_storm.jl")
# time: 2023-12-20 15:46:16 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 15:51:01 CET
# mode: julia
	includet("event_storm.jl")
# time: 2023-12-20 15:51:16 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 15:54:45 CET
# mode: shell
	cp ../pde/pde_storm_plot.jl ./plot_event_storm.jl
# time: 2023-12-20 15:54:50 CET
# mode: julia
	includet("plot_event_storm.jl")
# time: 2023-12-20 15:54:57 CET
# mode: pkg
	add PyPlot
# time: 2023-12-20 15:56:01 CET
# mode: julia
	includet("plot_event_storm.jl")
# time: 2023-12-20 15:56:37 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/")
# time: 2023-12-20 16:07:21 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:09:11 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/")
# time: 2023-12-20 16:10:16 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:11:48 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/")
# time: 2023-12-20 16:13:21 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:13:39 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/")
# time: 2023-12-20 16:13:51 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 16:13:56 CET
# mode: julia
	using PyPlot
# time: 2023-12-20 16:13:57 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 16:14:00 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/");
# time: 2023-12-20 16:14:38 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:14:47 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/");
# time: 2023-12-20 16:14:57 CET
# mode: julia
	out.q
# time: 2023-12-20 16:15:34 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:15:39 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/");
# time: 2023-12-20 16:16:46 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:16:57 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/");
# time: 2023-12-20 16:17:14 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:17:23 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/");
# time: 2023-12-20 16:34:55 CET
# mode: julia
	EventStorm.load_waccm("data/waccm_fg_l38.dat")
# time: 2023-12-20 16:35:05 CET
# mode: julia
	#EventStorm.load_waccm("data/waccm_fg_l38.dat")
# time: 2023-12-20 16:35:16 CET
# mode: julia
	head data/waccm_fg_l38.dat
# time: 2023-12-20 16:35:18 CET
# mode: shell
	head data/waccm_fg_l38.dat
# time: 2023-12-20 16:35:52 CET
# mode: julia
	EventStorm.load_waccm("data/waccm_fg_l38.dat")
# time: 2023-12-20 16:36:34 CET
# mode: shell
	head data/waccm_fg_l38.dat
# time: 2023-12-20 16:37:57 CET
# mode: julia
	EventStorm.load_waccm("data/waccm_fg_l38.dat")
# time: 2023-12-20 16:39:31 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 16:39:40 CET
# mode: julia
	waccm = EventStorm.load_waccm("data/waccm_fg_l38.dat");
# time: 2023-12-20 16:40:00 CET
# mode: julia
	plt.plot(waccm.O, waccm.z ./ co.kilo)
# time: 2023-12-20 16:40:05 CET
# mode: julia
	using Constants: co
# time: 2023-12-20 16:40:06 CET
# mode: julia
	plt.plot(waccm.O, waccm.z ./ co.kilo)
# time: 2023-12-20 16:40:26 CET
# mode: julia
	plt.clf(); plt.plot(waccm.O ./ co.centi^-3, waccm.z ./ co.kilo)
# time: 2023-12-20 16:41:56 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:43:24 CET
# mode: julia
	waccm.z[1]
# time: 2023-12-20 16:45:35 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:46:08 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/");
# time: 2023-12-20 16:46:25 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:46:41 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/");
# time: 2023-12-20 16:48:14 CET
# mode: julia
	maximum(out.no)
# time: 2023-12-20 16:48:32 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 16:48:43 CET
# mode: julia
	plt.plot(out.no3, out.z)
# time: 2023-12-20 16:49:02 CET
# mode: julia
	plt.plot(out.no3 ./ co.centi^-3, out.z)
# time: 2023-12-20 16:49:07 CET
# mode: julia
	plt.clf(); plt.plot(out.no3 ./ co.centi^-3, out.z)
# time: 2023-12-20 16:50:50 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/");
# time: 2023-12-20 16:55:05 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 16:57:01 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/");
# time: 2023-12-20 16:57:34 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", c="k");
# time: 2023-12-20 16:57:59 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "O-", c="k");
# time: 2023-12-20 16:58:33 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 16:58:41 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g");
# time: 2023-12-20 16:58:59 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 16:59:01 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g");
# time: 2023-12-20 16:59:22 CET
# mode: julia
	#PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "O-", c="");
# time: 2023-12-20 16:59:25 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="b");
# time: 2023-12-20 16:59:34 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "O-", c="g");
# time: 2023-12-20 16:59:56 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "O2-", c="k");
# time: 2023-12-20 17:01:27 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="b");
# time: 2023-12-20 17:17:01 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 17:17:05 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="b", lw=0.5);
# time: 2023-12-20 17:17:12 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="b", lw=0.75);
# time: 2023-12-20 17:17:16 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 17:17:17 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="b", lw=0.75);
# time: 2023-12-20 17:17:46 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "O2-", c="r", lw=0.75);
# time: 2023-12-20 17:17:57 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "2-", c="b", lw=0.75);
# time: 2023-12-20 17:18:01 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "O-", c="b", lw=0.75);
# time: 2023-12-20 17:19:47 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 17:19:58 CET
# mode: julia
	plt.plot(out.no, out.z)
# time: 2023-12-20 17:22:19 CET
# mode: julia
	plt.plot(out.no3, out.z)
# time: 2023-12-20 17:28:08 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 17:28:53 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 17:29:02 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 17:29:04 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 17:34:33 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 17:35:09 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 17:35:22 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 17:35:24 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 17:38:15 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "O2-", c="r", lw=0.75);
# time: 2023-12-20 17:38:39 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "O-", c="b", lw=0.75);
# time: 2023-12-20 17:40:00 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 17:40:23 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 17:40:27 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 17:44:54 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 17:45:41 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 17:45:44 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 18:17:01 CET
# mode: julia
	Vector
# time: 2023-12-20 18:20:01 CET
# mode: julia
	using DiffEqsCallbacks
# time: 2023-12-20 18:20:15 CET
# mode: julia
	using DiffEqCallbacks
# time: 2023-12-20 18:26:25 CET
# mode: pkg
	add Distributions
# time: 2023-12-20 18:31:52 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:32:14 CET
# mode: julia
	includet("plot_event_storm.jl")
# time: 2023-12-20 18:32:25 CET
# mode: julia
	includet("event_storm.jl")
# time: 2023-12-20 18:32:41 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:33:39 CET
# mode: julia
	includet("event_storm.jl")
# time: 2023-12-20 18:33:55 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:34:15 CET
# mode: julia
	includet("plot_event_storm.jl")
# time: 2023-12-20 18:34:35 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 18:39:41 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:40:17 CET
# mode: julia
	out.storm_rate
# time: 2023-12-20 18:40:37 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:40:40 CET
# mode: julia
	out.ν
# time: 2023-12-20 18:41:40 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:44:07 CET
# mode: julia
	out.sol.t
# time: 2023-12-20 18:44:19 CET
# mode: julia
	out.sol.u[2]
# time: 2023-12-20 18:44:22 CET
# mode: julia
	out.sol.u[1]
# time: 2023-12-20 18:45:19 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:45:24 CET
# mode: julia
	out.sol.t
# time: 2023-12-20 18:45:29 CET
# mode: julia
	out.sol.u
# time: 2023-12-20 18:46:37 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:48:57 CET
# mode: julia
	out.conf.ne1
# time: 2023-12-20 18:49:04 CET
# mode: julia
	out.conf.latt
# time: 2023-12-20 18:49:15 CET
# mode: julia
	out.conf.latt[20:end]
# time: 2023-12-20 18:49:29 CET
# mode: julia
	out.conf.ne1
# time: 2023-12-20 18:49:44 CET
# mode: julia
	out.conf.nz
# time: 2023-12-20 18:49:46 CET
# mode: julia
	out.conf.z
# time: 2023-12-20 18:49:58 CET
# mode: julia
	out.conf.ngas
# time: 2023-12-20 18:50:51 CET
# mode: julia
	EventStorm.logattenuation!(out.conf.latt, out.conf.z, out.conf.ne1, out.conf.ngas)
# time: 2023-12-20 18:50:55 CET
# mode: julia
	out.conf.latt[20:end]
# time: 2023-12-20 18:52:14 CET
# mode: julia
	EventStorm.logattenuation!(out.conf.latt, out.conf.z, out.conf.ne1, out.conf.ngas)
# time: 2023-12-20 18:52:16 CET
# mode: julia
	out.conf.latt[20:end]
# time: 2023-12-20 18:52:44 CET
# mode: julia
	log(-0.2)
# time: 2023-12-20 18:53:27 CET
# mode: julia
	EventStorm.logattenuation!(out.conf.latt, out.conf.z, out.conf.ne1, out.conf.ngas)
# time: 2023-12-20 18:56:23 CET
# mode: julia
	300197.49629063014 / 1.7e-11
# time: 2023-12-20 18:57:43 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:57:46 CET
# mode: julia
	out.ngas
# time: 2023-12-20 18:58:07 CET
# mode: julia
	out.gas_profiles
# time: 2023-12-20 18:58:34 CET
# mode: julia
	out.waccm
# time: 2023-12-20 18:58:37 CET
# mode: julia
	out.waccm_profiles
# time: 2023-12-20 18:58:54 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 18:58:57 CET
# mode: julia
	out.ngas
# time: 2023-12-20 18:59:07 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 19:04:28 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 19:04:44 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 19:04:51 CET
# mode: julia
	using PyPlot
# time: 2023-12-20 19:04:52 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 19:04:58 CET
# mode: julia
	PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 19:08:13 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 19:10:07 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 19:45:08 CET
# mode: julia
	includet("event_storm.jl")
# time: 2023-12-20 19:46:36 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 19:47:52 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 19:47:56 CET
# mode: julia
	using PyPlot
# time: 2023-12-20 19:48:04 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 19:48:13 CET
# mode: julia
	includet("plot_event_storm.jl")
# time: 2023-12-20 19:48:17 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="g", lw=0.75);
# time: 2023-12-20 20:07:56 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 20:09:51 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="k", lw=0.75);
# time: 2023-12-20 20:11:52 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 20:13:40 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="k", lw=0.75);
# time: 2023-12-20 20:16:38 CET
# mode: julia
	@edit PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="k", lw=0.75);
# time: 2023-12-20 20:19:48 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/data/storm/events00/", "e", c="k", lw=0.75);
# time: 2023-12-20 20:20:34 CET
# mode: julia
	plt.colorbar()
# time: 2023-12-20 20:21:27 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 20:23:18 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/data/storm/events00/", "e", c="k", lw=0.75);
# time: 2023-12-20 20:24:00 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 20:25:13 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/data/storm/events00/", "e", c="k", lw=0.75);
# time: 2023-12-20 20:26:16 CET
# mode: julia
	plt.axvline(3600)
# time: 2023-12-20 20:26:20 CET
# mode: julia
	plt.axvline(2 * 3600)
# time: 2023-12-20 20:28:39 CET
# mode: julia
	@edit PlotEventStorm.plot_profiles("/Users/luque/data/storm/events00/", "e", c="k", lw=0.75);
# time: 2023-12-20 20:30:34 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_slice("/Users/luque/data/storm/events00/", 80e3, "e", c="k", lw=0.75);
# time: 2023-12-20 20:30:53 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_slice("/Users/luque/data/storm/events00/", 82e3, "e", c="k", lw=0.75);
# time: 2023-12-20 21:45:36 CET
# mode: julia
	out.rs
# time: 2023-12-20 21:45:53 CET
# mode: julia
	using Chemise
# time: 2023-12-20 21:46:01 CET
# mode: julia
	init(out.rs, Dict(), 0.0)
# time: 2023-12-20 21:46:07 CET
# mode: julia
	Chemise.init(out.rs, Dict(), 0.0)
# time: 2023-12-20 21:46:18 CET
# mode: julia
	n0 = Chemise.init(out.rs, Dict(), 0.0)
# time: 2023-12-20 21:46:34 CET
# mode: julia
	Chemise.derivs(n0, outrs)
# time: 2023-12-20 21:46:37 CET
# mode: julia
	Chemise.derivs(n0, out.rs)
# time: 2023-12-20 21:46:52 CET
# mode: help
	Chemise.derivs(n0, out.rs)
# time: 2023-12-20 21:47:21 CET
# mode: help
	Chemise.derivs(n0, out.rs, 50)
# time: 2023-12-20 21:47:25 CET
# mode: help
	Chemise.derivs(n0, out.rs; x=50)
# time: 2023-12-20 21:47:27 CET
# mode: julia
	Chemise.derivs(n0, out.rs; x=50)
# time: 2023-12-20 21:47:36 CET
# mode: julia
	using Benchmarktools
# time: 2023-12-20 21:47:43 CET
# mode: julia
	using BenchmarkTools
# time: 2023-12-20 21:47:53 CET
# mode: julia
	@benchmark Chemise.derivs($n0, $out.rs; x=50)
# time: 2023-12-20 21:52:20 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 21:53:58 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_slice("/Users/luque/data/storm/nag00/", 82e3, "e", c="k", lw=0.75);
# time: 2023-12-20 21:54:06 CET
# mode: julia
	plt.clf()
# time: 2023-12-20 21:54:07 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_slice("/Users/luque/data/storm/nag00/", 82e3, "e", c="k", lw=0.75);
# time: 2023-12-20 21:54:21 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/data/storm/nag00/", "e", c="k", lw=0.75);
# time: 2023-12-20 21:54:54 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 21:56:03 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/data/storm/berger00/", "e", c="k", lw=0.75);
# time: 2023-12-20 21:56:22 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_slice("/Users/luque/data/storm/nag00/", 85e3, "e", c="k", lw=0.75);
# time: 2023-12-20 21:57:19 CET
# mode: julia
	co.hour
# time: 2023-12-20 21:57:23 CET
# mode: julia
	using Constants: co
# time: 2023-12-20 21:57:25 CET
# mode: julia
	co.hour
# time: 2023-12-20 21:59:05 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 22:04:15 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/data/storm/berger00/", "e", c="k", lw=0.75);
# time: 2023-12-20 22:04:54 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_slice("/Users/luque/data/storm/berger00/", 90e3, "e", c="k", lw=0.75);
# time: 2023-12-20 22:05:39 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-20 22:07:19 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_slice("/Users/luque/data/storm/berger00/", 90e3, "e", c="k", lw=0.75);
# time: 2023-12-21 18:38:23 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/data/storm/berger00/", "e", c="k", lw=0.75);
# time: 2023-12-21 22:59:50 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_slice("/Users/luque/data/storm/berger00/", 85e3, "e", c="k", lw=0.75);
# time: 2023-12-21 23:00:05 CET
# mode: julia
	PlotEventStorm.plot_slice("/Users/luque/data/storm/berger00/", 87e3, "e", c="k", lw=0.75);
# time: 2023-12-21 23:00:15 CET
# mode: julia
	PlotEventStorm.plot_slice("/Users/luque/data/storm/berger00/", 88e3, "e", c="k", lw=0.75);
# time: 2023-12-21 23:00:20 CET
# mode: julia
	PlotEventStorm.plot_slice("/Users/luque/data/storm/berger00/", 89e3, "e", c="k", lw=0.75);
# time: 2023-12-21 23:00:29 CET
# mode: julia
	PlotEventStorm.plot_slice("/Users/luque/data/storm/berger00/", 90e3, "e", c="k", lw=0.75);
# time: 2023-12-23 20:20:18 CET
# mode: julia
	@edit EventStorm.main();
# time: 2023-12-23 20:21:36 CET
# mode: julia
	out = EventStorm.main();
# time: 2023-12-23 20:22:40 CET
# mode: julia
	PlotEventStorm.plot_slice("/Users/luque/data/storm/nag00/", 90e3, "e", c="k", lw=0.75);
# time: 2023-12-23 20:22:47 CET
# mode: julia
	plt.clf()
# time: 2023-12-23 20:22:48 CET
# mode: julia
	PlotEventStorm.plot_slice("/Users/luque/data/storm/nag00/", 90e3, "e", c="k", lw=0.75);
# time: 2023-12-23 20:24:33 CET
# mode: julia
	out.tspan
# time: 2023-12-23 20:25:13 CET
# mode: julia
	length(out.sol.t)
# time: 2023-12-23 20:25:23 CET
# mode: julia
	out.tspan / 300
# time: 2023-12-23 20:25:27 CET
# mode: julia
	out.tspan[2] / 300
# time: 2023-12-23 20:25:54 CET
# mode: julia
	out.sol.t
# time: 2023-12-23 20:26:18 CET
# mode: julia
	#PlotEventStorm.plot_slice("/Users/luque/data/storm/nag00/", 90e3, "e", c="k", lw=0.75);
# time: 2023-12-23 20:26:52 CET
# mode: shell
	cat /Users/luque/data/storm/nag00/times.csv
# time: 2023-12-24 00:24:33 CET
# mode: julia
	PlotEventStorm.plot_slice("/Users/luque/data/storm/nag01/", 90e3, "e", c="k", lw=0.75);
# time: 2023-12-24 00:24:41 CET
# mode: julia
	plt.clf
# time: 2023-12-24 00:24:42 CET
# mode: julia
	plt.clf()
# time: 2023-12-24 00:24:43 CET
# mode: julia
	plt.clf
# time: 2023-12-24 00:24:45 CET
# mode: julia
	PlotEventStorm.plot_slice("/Users/luque/data/storm/nag01/", 90e3, "e", c="k", lw=0.75);
# time: 2023-12-24 00:25:23 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/data/storm/nag01/", "e", c="k", lw=0.75);
# time: 2023-12-24 00:26:49 CET
# mode: julia
	plt.clf(); PlotEventStorm.plot_slice("/Users/luque/data/storm/nag01/", 90e3, "e", c="k", lw=0.75);
# time: 2024-01-01 13:20:11 CET
# mode: shell
	cp ~/projects/amr/Steamroll/src/lxcat.jl .
# time: 2024-01-01 13:21:06 CET
# mode: shell
	ls
# time: 2024-01-01 13:21:09 CET
# mode: shell
	ls data
# time: 2024-01-01 13:21:14 CET
# mode: shell
	ls data/swarm/
# time: 2024-01-01 13:22:03 CET
# mode: shell
	cp ~/projects/amr/Steamroll/data/LxCat_Phelps_20230914.txt data/swarm/
# time: 2024-01-01 13:22:19 CET
# mode: julia
	edit("lxcat.jl")
# time: 2024-01-01 13:22:33 CET
# mode: julia
	includet("lxcat.jl")
# time: 2024-01-01 13:23:10 CET
# mode: julia
	LxCatSwarmData.load("data/swarm/LxCat_Phelps_20230914.txt")
# time: 2024-01-01 13:29:08 CET
# mode: julia
	using Infilitrator
# time: 2024-01-01 13:29:15 CET
# mode: julia
	using Infilitrate
# time: 2024-01-01 13:29:20 CET
# mode: julia
	using Infiltrator
# time: 2024-01-01 13:29:46 CET
# mode: julia
	LxCatSwarmData.load("data/swarm/LxCat_Phelps_20230914.txt")
# time: 2024-01-01 13:35:55 CET
# mode: julia
	lx = LxCatSwarmData.load("data/swarm/LxCat_Phelps_20230914.txt")
# time: 2024-01-01 13:37:49 CET
# mode: julia
	lx.explain
# time: 2024-01-01 13:38:03 CET
# mode: julia
	lx.data.C15
# time: 2024-01-01 13:38:12 CET
# mode: julia
	lx.data.en
# time: 2024-01-01 13:38:28 CET
# mode: julia
	using Chemise
# time: 2024-01-01 13:50:28 CET
# mode: julia
	eachcol(lx.data)[:I20]
# time: 2024-01-01 13:51:11 CET
# mode: julia
	loadtable
# time: 2024-01-01 13:53:11 CET
# mode: julia
	loadtable(eachcol(lx.data), xcol=:en)
# time: 2024-01-01 13:55:53 CET
# mode: julia
	typeof(loadtable(eachcol(lx.data), xcol=:en))
# time: 2024-01-01 13:56:23 CET
# mode: julia
	tlb = loadtable(eachcol(lx.data), xcol=:en)
# time: 2024-01-01 13:56:53 CET
# mode: julia
	RateLookup
# time: 2024-01-01 13:57:05 CET
# mode: julia
	tbl = loadtable(eachcol(lx.data), xcol=:en)
# time: 2024-01-01 13:57:22 CET
# mode: julia
	RateLookup(tbl, :A11)
# time: 2024-01-01 13:57:37 CET
# mode: julia
	RateLookup(tbl, :A12)
# time: 2024-01-01 13:57:57 CET
# mode: julia
	lx.explain
# time: 2024-01-01 14:07:37 CET
# mode: julia
	sort(pairs(lx.explain))
# time: 2024-01-01 14:09:41 CET
# mode: julia
	lx
# time: 2024-01-01 14:10:50 CET
# mode: help
	pad
# time: 2024-01-01 14:10:56 CET
# mode: help
	lpad
# time: 2024-01-01 14:12:07 CET
# mode: julia
	lx
# time: 2024-01-01 14:43:19 CET
# mode: julia
	@ref A => b "ref"
# time: 2024-01-01 16:39:44 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-01 16:47:52 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-01 16:49:06 CET
# mode: julia
	out.frs
# time: 2024-01-01 16:49:15 CET
# mode: julia
	using Chemise
# time: 2024-01-01 16:50:01 CET
# mode: julia
	Chemise.init(out.frs, Dict(:e => 1e20))
# time: 2024-01-01 16:50:16 CET
# mode: julia
	n0 = Chemise.init(out.frs, Dict(:e => 1e20))
# time: 2024-01-01 16:50:20 CET
# mode: julia
	derivs
# time: 2024-01-01 16:50:22 CET
# mode: help
	derivs
# time: 2024-01-01 16:52:53 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-01 16:53:13 CET
# mode: pkg
	resolve
# time: 2024-01-01 16:53:24 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-01 16:53:47 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-01 16:54:03 CET
# mode: julia
	using Chemise
# time: 2024-01-01 16:54:06 CET
# mode: help
	derivs
# time: 2024-01-01 16:54:45 CET
# mode: help
	Chemise.init(out.frs, Dict(:e => 1e20))
# time: 2024-01-01 16:54:51 CET
# mode: julia
	n0 = Chemise.init(out.frs, Dict(:e => 1e20))
# time: 2024-01-01 16:55:13 CET
# mode: julia
	derivs(n0, out.frs, 100.0)
# time: 2024-01-01 16:56:01 CET
# mode: julia
	derivs(n0, out.frs, 100.0, x=100)
# time: 2024-01-01 16:56:05 CET
# mode: julia
	derivs(n0, out.frs, 100.0, x=10)
# time: 2024-01-01 16:56:18 CET
# mode: julia
	n0
# time: 2024-01-01 16:56:39 CET
# mode: julia
	Chemise.species(out.frs)
# time: 2024-01-01 16:58:24 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-01 16:58:30 CET
# mode: julia
	derivs(n0, out.frs, 100.0, x=10)
# time: 2024-01-01 16:58:44 CET
# mode: julia
	derivs(n0, out.frs, 150.0, x=10)
# time: 2024-01-01 16:58:50 CET
# mode: julia
	derivs(n0, out.frs, 180.0, x=10)
# time: 2024-01-01 16:58:58 CET
# mode: julia
	derivs(n0, out.frs, 200.0, x=10)
# time: 2024-01-01 16:59:12 CET
# mode: julia
	derivs(n0, out.frs, 500.0, x=10)
# time: 2024-01-01 17:01:28 CET
# mode: julia
	Chemise.species(out.frs)
# time: 2024-01-01 17:01:47 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-01 17:01:50 CET
# mode: julia
	derivs(n0, out.frs, 500.0, x=10)
# time: 2024-01-01 17:01:53 CET
# mode: julia
	derivs(n0, out.frs, 100.0, x=10)
# time: 2024-01-01 17:01:58 CET
# mode: julia
	derivs(n0, out.frs, 130.0, x=10)
# time: 2024-01-01 17:02:03 CET
# mode: julia
	derivs(n0, out.frs, 120.0, x=10)
# time: 2024-01-01 17:02:15 CET
# mode: julia
	derivs(n0, out.frs, 115.0, x=10)
# time: 2024-01-01 17:02:23 CET
# mode: julia
	@code_native derivs(n0, out.frs, 115.0, x=10)
# time: 2024-01-01 20:27:15 CET
# mode: julia
	species
# time: 2024-01-01 20:27:17 CET
# mode: help
	species
# time: 2024-01-01 20:27:23 CET
# mode: julia
	species(out.rs)
# time: 2024-01-01 20:27:26 CET
# mode: julia
	species(out.frs)
# time: 2024-01-01 20:42:18 CET
# mode: julia
	idx(out.rs, :e)
# time: 2024-01-01 20:42:21 CET
# mode: julia
	idx(out.rs, :elec)
# time: 2024-01-01 20:42:26 CET
# mode: julia
	typeof(idx(out.rs, :elec))
# time: 2024-01-01 20:48:29 CET
# mode: julia
	fn0 = Chemise.init(out.frs, Dict(:e => 1e20))
# time: 2024-01-01 20:48:35 CET
# mode: julia
	n0 = Chemise.init(out.rs, Dict(:e => 1e20))
# time: 2024-01-01 20:48:58 CET
# mode: julia
	EventStorm.mapspecs(out.frs, out.rs, n0)
# time: 2024-01-01 20:49:19 CET
# mode: julia
	@benchmark EventStorm.mapspecs($out.frs, $out.rs, $n0)
# time: 2024-01-01 20:49:27 CET
# mode: julia
	using BenchmarkTools
# time: 2024-01-01 20:49:30 CET
# mode: julia
	@benchmark EventStorm.mapspecs($out.frs, $out.rs, $n0)
# time: 2024-01-02 13:36:44 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-02 13:48:23 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-02 13:49:28 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-02 13:49:45 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-02 14:05:18 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-02 14:06:02 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-02 14:23:58 CET
# mode: julia
	show(err)
# time: 2024-01-02 14:27:14 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-02 14:29:41 CET
# mode: julia
	includet("plot_event_storm.jl")
# time: 2024-01-02 14:30:21 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/localdata/storm/nag01/", "e", c="k", lw=0.75);
# time: 2024-01-02 14:30:27 CET
# mode: julia
	using PyPlot
# time: 2024-01-02 14:30:28 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/localdata/storm/nag01/", "e", c="k", lw=0.75);
# time: 2024-01-02 14:34:21 CET
# mode: pkg
	add ODEInterfaceDiffEq
# time: 2024-01-02 14:43:14 CET
# mode: julia
	out.rs
# time: 2024-01-02 20:45:25 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-03 21:33:08 CET
# mode: julia
	using Cthulhu
# time: 2024-01-03 21:35:54 CET
# mode: julia
	includet("plot_event_storm.jl")
# time: 2024-01-03 21:36:06 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-03 21:36:21 CET
# mode: julia
	using Cthulhu
# time: 2024-01-03 21:36:53 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-03 21:37:50 CET
# mode: julia
	EventStorm.flash!(out.integrator)
# time: 2024-01-03 21:43:52 CET
# mode: julia
	using Infiltrator
# time: 2024-01-03 21:44:20 CET
# mode: julia
	EventStorm.flash!(out.integrator)
# time: 2024-01-03 23:26:44 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-03 23:27:27 CET
# mode: julia
	EventStorm.flash!(out.integrator)
# time: 2024-01-03 23:27:40 CET
# mode: julia
	@descend EventStorm.flash!(out.integrator)
# time: 2024-01-04 15:17:32 CET
# mode: help
	using DifferentialEquations
# time: 2024-01-04 15:17:36 CET
# mode: julia
	using DifferentialEquations
# time: 2024-01-04 15:17:41 CET
# mode: help
	init
# time: 2024-01-04 15:25:14 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-04 15:32:38 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/localdata/storm/nag01/", "e", c="k", lw=0.75);
# time: 2024-01-04 15:32:42 CET
# mode: julia
	using PyPlot
# time: 2024-01-04 15:32:43 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/localdata/storm/nag01/", "e", c="k", lw=0.75);
# time: 2024-01-04 15:33:38 CET
# mode: julia
	out = EventStorm.main(pre_relax=0);
# time: 2024-01-04 15:34:08 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/localdata/storm/nag01/", "e", c="k", lw=0.75);
# time: 2024-01-04 19:37:56 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-04 19:38:24 CET
# mode: julia
	plt.clf(); PlotEventStorm.pcolor_profiles("/Users/luque/localdata/storm/nag01/", "e", c="k", lw=0.75);
# time: 2024-01-04 19:38:43 CET
# mode: julia
	plt.colorbar()
# time: 2024-01-04 20:00:23 CET
# mode: julia
	nrmsis = EventStorm.load_nrlmsis("data/nlrmsis2.txt"")
# time: 2024-01-04 20:00:28 CET
# mode: julia
	nrmsis = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:02:31 CET
# mode: julia
	nrmsis = EventStorm.load_pnrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:03:47 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:04:05 CET
# mode: julia
	nrmsis = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:06:43 CET
# mode: shell
	cat startup.jl
# time: 2024-01-04 20:07:03 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:07:19 CET
# mode: julia
	nrmsis = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:08:45 CET
# mode: pkg
	rm Infiltrator
# time: 2024-01-04 20:08:58 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:09:28 CET
# mode: pkg
	rm Infiltrator
# time: 2024-01-04 20:09:32 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:09:41 CET
# mode: julia
	nrmsis = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:12:36 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:13:00 CET
# mode: julia
	nrmsis = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:15:32 CET
# mode: julia
	rlogger = Revise.debug_logger()
# time: 2024-01-04 20:15:46 CET
# mode: julia
	nrmsis = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:16:06 CET
# mode: julia
	logs = Revise.actions(rlogger)
# time: 2024-01-04 20:17:42 CET
# mode: julia
	Revise.revise()
# time: 2024-01-04 20:17:46 CET
# mode: julia
	nrmsis = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:18:26 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:19:10 CET
# mode: julia
	EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:21:02 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:21:18 CET
# mode: julia
	EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:23:41 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:23:59 CET
# mode: julia
	EventStorm.main()
# time: 2024-01-04 20:25:01 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:26:04 CET
# mode: julia
	EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:27:12 CET
# mode: julia
	df = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:27:27 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:27:46 CET
# mode: julia
	df = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:28:53 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:29:59 CET
# mode: julia
	df = EventStorm.load_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:30:04 CET
# mode: julia
	df = EventStorm.pload_nrlmsis("data/nlrmsis2.txt")
# time: 2024-01-04 20:31:05 CET
# mode: julia
	df = EventStorm.testrevise()
# time: 2024-01-04 20:33:04 CET
# mode: julia
	includet("load_data.jl")
# time: 2024-01-04 20:33:11 CET
# mode: julia
	testrevise()
# time: 2024-01-04 20:33:57 CET
# mode: pkg
	st Revise
# time: 2024-01-04 20:34:47 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:35:22 CET
# mode: julia
	EventStorm.testrevise()
# time: 2024-01-04 20:39:33 CET
# mode: shell
	mv data/nlrmsis2.txt data/nrlmsis2.dat
# time: 2024-01-04 20:39:45 CET
# mode: julia
	includet("load_data.jl")
# time: 2024-01-04 20:40:00 CET
# mode: julia
	load_nrlmsis("data/nrlmsis2.dat"")
# time: 2024-01-04 20:40:02 CET
# mode: julia
	load_nrlmsis("data/nrlmsis2.dat")
# time: 2024-01-04 20:40:15 CET
# mode: julia
	using DataFrames, CSV
# time: 2024-01-04 20:40:17 CET
# mode: julia
	load_nrlmsis("data/nrlmsis2.dat")
# time: 2024-01-04 20:41:52 CET
# mode: julia
	df = load_nrlmsis("data/nrlmsis2.dat")
# time: 2024-01-04 20:41:58 CET
# mode: julia
	names(df)
# time: 2024-01-04 20:42:28 CET
# mode: julia
	df = load_nrlmsis("data/nrlmsis2.dat")
# time: 2024-01-04 20:42:54 CET
# mode: help
	CSV.read
# time: 2024-01-04 20:43:24 CET
# mode: julia
	df = load_nrlmsis("data/nrlmsis2.dat")
# time: 2024-01-04 20:43:43 CET
# mode: julia
	using Constants:co
# time: 2024-01-04 20:43:45 CET
# mode: julia
	df = load_nrlmsis("data/nrlmsis2.dat")
# time: 2024-01-04 20:43:53 CET
# mode: julia
	using Interpolations
# time: 2024-01-04 20:43:55 CET
# mode: julia
	df = load_nrlmsis("data/nrlmsis2.dat")
# time: 2024-01-04 20:45:01 CET
# mode: julia
	includet("load_data.jl")
# time: 2024-01-04 20:48:44 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:48:56 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-04 20:49:30 CET
# mode: julia
	includet("event_storm.jl")
# time: 2024-01-04 20:51:51 CET
# mode: julia
	out = EventStorm.main();
# time: 2024-01-04 20:52:48 CET
# mode: julia
	out.nrlmsis.o2
# time: 2024-01-08 11:12:38 CET
# mode: julia
	out.nrlmsis.o2(100e3)
