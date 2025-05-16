# Numerical Experiments for  
**"Iteration Complexity and Asymptotic Analysis of Gradient-Type Methods with Non-Monotone Line Search on Riemannian Manifolds"**

This repository contains the implementation of the numerical experiments reported in our paper:

**Authors:**  
Qamrul Hasan Ansari, Orizon Pereira Ferreira, Moin Uddin, Mohammed Alshahrani

---

## 📄 About

This codebase accompanies the numerical section of the paper and evaluates various **gradient-type optimization algorithms** using **non-monotone line search techniques** on **Riemannian manifolds**.  

Experiments test multiple configurations of:
- **Search directions**
- **Step size strategies**
- **Non-monotone parameters**
- **Problem sizes**
- **Optimization problems** 

All results are automatically saved in `.csv`, `.txt`, and `.xlsx` formats, and performance profiles are generated.

---

## 📁 Folder Structure

```

├── src/
│   ├── main.jl             # Main script to run all experiments
│   ├── depts.jl            # Required dependencies
│   ├── types.jl            # Custom types and structs
│   ├── utils.jl            # Utility functions (printing, saving, etc.)
│   ├── functions.jl        # Algorithmic components (GMP1, line search, etc.)
│   ├── problems.jl         # Optimization problem definitions
├── Project.toml            # Julia environment configuration
├── results/                # Output directory (generated automatically)

````

---

## 🚀 Getting Started

### Prerequisites

- Julia 1.1 or later
- Recommended: Use a virtual environment via `Pkg`

### Setup

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate("")
```


## ▶️ Running the Experiments

Simply run the main experiment script:

```julia
include("src/main.jl")
```

This will:

* Solve all defined problems with all method combinations
* Print solution logs to the terminal
* Save the results to:

  * `results/solutions.csv`
  * `results/solutions.txt`
  * `results/solutions.xlsx`
* Generate performance profiles into:

  * `results/profiles/`

---

## 📊 Output Files

* `solutions.csv` / `.txt` / `.xlsx`: Tabulated solution data (iteration counts, timing, status, etc.)
* `profiles/`: Performance profiles grouped by:

  * Step size strategy
  * Non-monotone parameter scheme
  * Search direction

---


## 🛠️ Tuning Parameters

The main script `src/main.jl` is designed to make it easy to **experiment with different algorithmic settings**. You can adjust the following parameters at the top of the script inside the `begin` block:

### 🔧 Key Parameters

| Parameter            | Description                                                          |
| -------------------- | -------------------------------------------------------------------- |
| `szs`                | List of problem sizes (dimensions `n`) to test                       |
| `SearchDirections`   | Array of search direction strategies (e.g., `SD1`, `SD2`)            |
| `StepSizes`          | Array of step size rules (e.g., `St1`, `St2`, `St3`, `St4`)          |
| `MPs`                | Array of non-monotone parameter strategies (e.g., `MP0`, `MP1`, ...) |
| `problems`           | Array of test problems (e.g., `spd_problem1`, `po_problem2`, ...)    |
| `ϵ`                  | Stopping tolerance                                                   |
| `ν0`                 | Initial non-monotonicity control parameter                           |
| `max_iters`          | Maximum number of iterations per run                                 |
| `backtracking_iters` | Max number of backtracking steps                                     |
| `trycatch`           | Enable/disable error catching per experiment (`true` recommended)    |

### 🧪 Example: Focus on One Step Size

To only test with a specific step size strategy (e.g., `St2`), you can modify this line:

```julia
StepSizes = [St2]
```

### 🧪 Example: Single Problem, Small Size

To restrict experiments to one problem and a small dimension:

```julia
problems = [spd_problem2]
szs = [20]
```

### 🧪 Example: Use All Defaults

To test all combinations across all defined search directions, step sizes, and non-monotone parameters with the full problem set:

```julia
# Leave the original lists unchanged:
SearchDirections = [SD1, SD2]
StepSizes = [St1, St2, St3, St4]
MPs = [MP0, MP1, MP2, MP3, MP4, MP5]
```

## 📜 Citation

If you use this code or its results in your own research, please cite our paper:

Ansari, Q.H., Ferreira, O.P., Uddin, M., & Alshahrani, M.  
Iteration Complexity and Asymptotic Analysis of Gradient-Type Methods with Non-Monotone Line Search on Riemannian Manifolds.  [submitted]

---

## 🔒 License

This project is released under the MIT License.

---

## 🙏 Acknowledgments

This work is part of the collaborative research efforts between:

* **King Fahd University of Petroleum and Minerals (KFUPM)**
* **Universidade Federal de Goiás (UFG)**
* **Aligarh Muslim University (AMU)**



---

