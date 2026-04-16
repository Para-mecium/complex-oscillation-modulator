# `extract_periodic_orbit` I/O Contract

更新时间：2026-04-16

本文件定义 `PO_extract/extract_periodic_orbit.m` 的正式输入输出 contract。

该 contract 采用两层接口策略：

- 稳定层：面向公共调用方，后续版本应保持兼容
- 实验层：当前实现会返回，但不承诺字段稳定性或长期兼容

`PROJECT_STATUS.md` 用于记录工程背景和阶段决策；正式 I/O 定义以本文件为准。

## 1. 稳定输入层

稳定调用签名：

```matlab
result = extract_periodic_orbit(odefun, y0, parameter, opts)
```

稳定输入要求：

- `odefun`：函数句柄，支持 `odefun(t, y)` 或 `odefun(t, y, parameter)`
- `y0`：有限实数初值向量，行向量会被归一化为列向量
- `parameter`：透传到 `odefun`、`observableFcn`、`sectionFcn` 的参数对象
- `opts`：结构体；传空时使用默认配置

正式支持的稳定选项：

- `opts.orbit_solver`
- `opts.solver` 或 `opts.solver_name`
- `opts.solver_options`
- `opts.odeOptions`
- `opts.solver_tol`
- `opts.tspan`
- `opts.single_timespan`
- `opts.max_windows`
- `opts.observableFcn`
- `opts.autoSection`
- `opts.sectionFcn`
- `opts.sectionLevel`
- `opts.sectionDirection`
- `opts.transientFraction`
- `opts.transientTime`
- `opts.minCrossings`
- `opts.consecutiveCycles`
- `opts.samplesPerCycle`
- `opts.extractNumPoints`
- `opts.computeQualityDiagnostics`

优先级规则：

- 若提供 `opts.orbit_solver`，则统一入口按该模式路由；缺省值为 `'direct'`
- 当前正式支持的 solver 仅有 `'direct'` 与 `'matcont'`
- 若 `opts.tspan` 已提供，则优先使用 `opts.tspan`；否则才由 `single_timespan/max_windows/total_timespan` 派生总时间窗
- 若同时提供 `transientTime` 与 `transientFraction`，则优先使用 `transientTime`
- 若同时提供 `solver` 与 `solver_name`，则按归一化后的 solver handle 执行；`solver_name` 仅作为名称记录
- `solver_options` 与 `odeOptions` 会合并，内部事件函数会覆盖外部传入的 `Events`
- `computeQualityDiagnostics` 缺省为 `false`；仅在显式开启时才计算并返回 `diagnostics.quality`

## 2. 弃用兼容输入

以下选项当前仍受支持，但仅作为兼容入口保留，新调用方不应继续使用：

- `opts.event`

兼容语义：

- `event` 可映射为状态索引截面或自定义 section 句柄

已移除路由：

- `orbit_solver='direct_then_matcont'`
- `orbit_solver='matcont_internal'`

除 `opts.event` 外，任何不在正式支持范围内的选项字段都会按通用 unsupported-options 语义拒绝。
这包括历史兼容字段和实验字段，例如 `opts.backend`、`opts.min_events`、`opts.match_tol`、`opts.period_tol`。

后续如果移除这些兼容入口，应先更新本文件并同步迁移调用方。

## 3. 稳定输出层

返回值 `result` 中以下字段属于稳定 contract，并要求始终存在：

- `success`
- `has_orbit`
- `status`
- `code`
- `message`
- `period`
- `orbit_t`
- `orbit_y`
- `amplitude`
- `max_variable`
- `min_variable`
- `solver_used`

字段语义：

- `success`：仅当检测到严格收敛周期轨时为 `true`
- `has_orbit`：当返回了可用提取周期轨时为 `true`
- `status`：稳定状态字符串，调用方可以据此分支
- `code`：稳定状态码，便于旧调用方或脚本做数值判断
- `message`：面向人的说明文本；可读但不保证逐字稳定
- `period`：提取周期轨对应的周期
- `orbit_t`：从 `0` 开始、终点接近 `period` 的相对时间网格
- `orbit_y`：提取周期轨上的状态轨迹
- `amplitude`：按 `orbit_y` 每个状态分量计算得到的半峰峰值，等于 `(max_variable - min_variable) / 2`
- `max_variable`：按 `orbit_y` 每个状态分量计算得到的最大值
- `min_variable`：按 `orbit_y` 每个状态分量计算得到的最小值
- `solver_used`：顶层入口最终返回结果来自哪个 solver；固定为 `"direct"` 或 `"matcont"`

说明：

- `amplitude/max_variable/min_variable` 针对状态轨迹 `orbit_y`，不是 `observableFcn` 的统计量
- 稳定层不直接保证输出观测量轨迹；若调用方需要观测量，请基于 `orbit_y` 自行计算，或使用实验层字段

## 4. 稳定状态码与保证

### `converged_periodic_orbit`

- `code == 2`
- `success == true`
- `has_orbit == true`
- `period/orbit_t/orbit_y/amplitude/max_variable/min_variable` 有意义

### `candidate_periodic_orbit_not_converged`

- `code == 1`
- `success == false`
- `has_orbit == true`
- `period/orbit_t/orbit_y/amplitude/max_variable/min_variable` 有意义

### `no_periodic_orbit_detected_on_tspan`

- `code == 0`
- `success == false`
- `has_orbit == false`

### `decaying_to_equilibrium_or_nonoscillatory`

- `code == -1`
- `success == false`
- `has_orbit == false`

### `invalid_options`

- `code` 保持未定义数值状态，当前实现表现为 `NaN`
- `success == false`
- `has_orbit == false`

### `invalid_input`

- `code` 保持未定义数值状态，当前实现表现为 `NaN`
- `success == false`
- `has_orbit == false`

### `solver_failed`

- `code` 保持未定义数值状态，当前实现表现为 `NaN`
- `success == false`
- `has_orbit == false`

### `detection_failed`

- `code` 保持未定义数值状态，当前实现表现为 `NaN`
- `success == false`
- `has_orbit == false`

对所有 `has_orbit == false` 的状态，稳定层只保证这些字段存在；不保证其内容可继续解释。调用方应先检查 `has_orbit`，再消费轨道与特征字段。

## 5. 实验层

以下字段当前会返回，但不纳入正式稳定 contract：

- `event_times`
- `event_states`
- `backend_used`
- `orbit`
- `diagnostics`
- `raw`

实验层约束：

- 顶层字段是否存在可以保留，但其内部 schema 不冻结
- `diagnostics` 可用于调试、日志和实验分析，但字段名、嵌套结构、数值定义都可能变化
- `diagnostics.quality` 仅在 `opts.computeQualityDiagnostics=true` 时返回
- `raw` 属于内部调试输出，不应作为公共 API 使用
- `diagnostics` 当前还承载 direct 检测细节、MATCONT parameter-return、轨道质量指标等实验信息

## 6. 当前不纳入稳定 contract 的求解器专用输入

以下输入项当前存在于实现中，但属于求解器专用扩展，不纳入统一稳定 contract：

- `matcont_root`
- `matcont_odefile`
- `matcont_active_parameter`
- `matcont_ntst`
- `matcont_ncol`
- `matcont_tolerance`
- `matcont_parameter_values`
- `matcont_window_timespan`
- `matcont_cycle_window_factor`
- `matcont_options`
- `matcont_parameter_tolerance`
- `matcont_return_max_points`
- `matcont_return_scan_both_directions`

这些输入会影响 MATCONT 后端的真实求解流程，但不保证被 direct 后端消费。

未纳入支持范围的 MATCONT 扩展字段同样按通用 unsupported-options 语义拒绝。
当前实现不再支持 external seed 初始化；`seed_orbit` 和 `seed_quality_reference` 都属于不支持字段，MATCONT solver 固定使用 internal preprocessing。
