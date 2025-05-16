# import Manifolds: AbstractManifold
@kwdef struct ManifoldProblem
	name::String = ""
	M::AbstractManifold
	f::Function
	∇f::Function
	n::Int
	λL::Function
	x0::Union{Array{Float64}}
end

@kwdef struct Stats
	iterations::Number = 0.0
	backtracking_iterations::Number = 0.0
	function_evaluations::Number = 0.0
	grad_evaluations::Number = 0.0
	time::Number = 0.0
end
Stats(s::Stats, t::Float64) = Stats(
	iterations = s.iterations,
	backtracking_iterations = s.backtracking_iterations,
	function_evaluations = s.function_evaluations,
	grad_evaluations = s.grad_evaluations,
	time = t,
)
abstract type Message end
struct SuccessMessage <: Message
	message::String
end
struct FailureMessage <: Message
	message::String
end
@kwdef struct Solution
	x::Any
	problem::ManifoldProblem
	search_direction::String
	step_size::String
	nonmonotone_parameter::String
	stats::Stats = Stats()
	message::Message
	norm_x_d::Float64
end

function Base.show(io::IO, sol::Solution)
	println(io, "Solution:")
	println(io, " Problem Name: ", sol.problem.name)
	println(io, "  Best value: ", isnothing(sol.x) ? "" : sol.problem.f(sol.x))
	println(io, "  ||x-d|| : ", sol.norm_x_d)
	println(io, "  Search Direction: ", sol.search_direction)
	println(io, "  Step Size: ", sol.step_size)
	println(io, "  Nonmonotone Parameter: ", sol.nonmonotone_parameter)
	println(io, "  Stats: ")
	println(io, "       Time (seconds): ", sol.stats.time)
	println(io, "       Iterations: ", sol.stats.iterations)
	println(io, "       Backtracking Iterations: ", sol.stats.backtracking_iterations)
	println(io, "       Function Evaluations: ", sol.stats.function_evaluations)
	println(io, "       Gradient Evaluations: ", sol.stats.grad_evaluations)
	println(io, "  Message: ", sol.message.message)
end
