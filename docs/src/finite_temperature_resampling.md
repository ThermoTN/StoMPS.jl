# 从归一化热纯化做 Born sampling

本节推导 StoMPS 的热纯化 Born 工作流：用 FiniteMPS 的 SETTN 和
tanTRG 把带对称性的热 MPO 演化到 `beta0`，对纯化的辅助物理腿做带
对称性的顺序 Born 测量，然后把每个归一化样本继续演化到更大的
逆温。

## 1. 热 MPO 的归一化

定义

```text
Rβ = exp(-βH/2)
Z(β) = Tr exp(-βH)
R̄β = Rβ / sqrt(Z(β))
```

因此 `norm(R̄β) = 1`，并且

```text
R̄β R̄β† = exp(-βH) / Z(β) = ρβ.
```

SETTN/tanTRG 每到达一个新的 beta 点，都必须在归一化前累加

```text
logZ ← logZ + 2 log(norm(R)),
```

再对 MPO 调用 `normalize!`。传给 Born API 的 `thermal_mpo` 必须是这个
归一化的 `R̄β₀`。

## 2. 纯化 Born 概率

将局域 thermal MPO 张量记为

```text
R[left, physical_out; physical_in, right].
```

在 `physical_in` 上选择一组与对称性相容的完备正交投影
`{P_a}`。对于完整局域配置 `a`，有

```text
p(a) = norm(R̄β₀ P_a)²
ρ_a = R̄β₀ P_a R̄β₀† / p(a).
```

完备性给出

```text
Σa p(a) ρ_a = R̄β₀ R̄β₀† = ρβ₀.
```

下面的 Z₂ 和 SU2SU2 sampler 共用从左到右移动 canonical center 的顺序条件
概率引擎。在第 `i` 个点，每个局域投影候选张量的范数平方是条件
概率；抽中后先归一化当前张量，再移动 center。返回的 rank-4 MPS
保留 thermal MPO 的虚拟键，第三腿保留被选子空间。

### 2.1 Z₂ rank-one 投影

`z2_purification_born_sample` 投影到 charge 0 和 charge 1 的一维
sector marker。此时 `P_a = |a⟩⟨a|`，所以条件对象可以写成归一化
纯态 `R̄|a⟩/sqrt(p(a))`。用 `z2_born_space(0, 1)` 构造物理空间，
用 `z2_born_space(1)` 构造 Ising `sigma_z` 算符的奇辅助空间，以保证
TensorKit storage 一致。

### 2.2 SO(4) 退化投影

半满 Hubbard 的局域空间是

```text
(η=1/2, S=0) ⊕ (η=0, S=1/2).
```

`su2su2_purification_born_sample` 使用 `{P_eta,P_spin}`。每个 sector 的
多重度为一，但 SU(2) carrier 维数为二；因此它是退化投影，条件
对象仍是带二维 marker 的纯化。不再选择 `m_eta` 或 `m_S`，否则会
破坏对应的 SU(2)。这些完备的局域对称性投影重构原热密度矩阵。

## 3. 样本的继续演化

Born 概率已经由样本出现频率表示，不能再乘入热力学权重。因此
开始 continuation 时

```text
relative_log_weight(β₀) = 0.
```

对于相邻的逆温点，样本按 `exp(-delta_beta H / 2)` 演化，所以 TDVP
时间步是 `-delta_beta / 2`。每步结束后，在归一化之前更新

```text
relative_log_weight ← relative_log_weight + 2 log(norm(sample)).
```

`continue_born_sample!` 统一实现 beta grid、时间步、权重与归一化约定。
调用者通过

```julia
step!(sample, -delta_beta / 2, beta_index, rng)
observe(sample, beta, beta_index)
```

配置 FiniteMPS 演化算法和标量观测量。`continue_born_sample!` 把显式
RNG 传给 `step!`；回调负责 FiniteMPS 所需的默认 RNG 设种。
共享驱动在操作边界设置种子，具体约定见 [Born API](design/born_sampling.md)。

## 4. 带 log weight 的估计

对于样本观测量 `Oa`，热平均是

```text
⟨O⟩β = Σa exp(ℓa) Oa / Σa exp(ℓa).
```

`logweighted_mean` 使用 `m = maximum(ℓa)` 做 log shift，避免指数上溢。
绝对配分函数估计为

```text
logZ(β) = logZ(β₀) + log(mean(exp(ℓa))),
```

对应 `estimate_logpartition(logZ0, relative_log_weights)`。

一般 Monte Carlo 中，配置已经按 `p(a)` 抽出，不得再把
`log(p(a))` 加到样本权重。仅在小尺寸确定性完整枚举中，因为每个
配置只被等频访问一次，oracle 才使用
`log(p(a)) + relative_log_weight` 做求和。

## 5. 验证要求

公共的 Born 回归要求同时验证：

1. 所有完整输入构型的直接列概率之和为一；
2. 顺序条件概率之积等于直接列概率；
3. 顺序塌缩态与独立直接收缩列的归一化 overlap 为一；
4. Z₂ 的 `Σa p(a)|ψa⟩⟨ψa|` 重构 `R̄R̄†` 并给出相同的 TFIM
   热能；退化投影则使用 `Σa p(a)ρ_a`；
5. continuation 使用严格递增 beta grid、`-delta_beta/2`、归一化前的
   norm 以及初始零相对权重；
6. 两站点 TFIM 从 `beta0=0` 完整枚举 Born 构型并做稠密精确虚时
   演化后，每站点能量与 ED 的绝对误差不超过 `1e-9`，或对于有限
   抽样回归不超过其统计 error bar。完整枚举的 error bar 为零。
7. SO(4) 后端额外检查 `P_eta+P_spin=I`、两投影正交、二维
   marker 被保留，以及两站点 Hubbard 完整四配置的每站点能量
   与 ED 绝对误差不超过 `1e-9`。
