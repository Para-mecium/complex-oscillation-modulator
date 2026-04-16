# 振荡特征精确计算工程记录

更新时间：2026-04-16

正式 I/O contract 见 [IO_CONTRACT.md](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/IO_CONTRACT.md)。本文件保留工程背景、决策和开发路径，不再承担公共接口定义。

## 1. 工程目标

目标是建立一个可维护的代码库，用于对“各类振荡调控算法”得到的参数，稳定且尽可能精确地计算其对应周期解的目标振荡特征。这里的“振荡特征”至少应包括：

- 周期 `period`
- 指定观测量的振幅 `amplitude`
- 最大值 / 最小值 `max / min`
- 必要时扩展到相位差、均值、占空比、峰宽、峰时序等

从工程角度看，这个目标至少包含两层能力：

1. 给定参数，可靠求出对应的周期轨 `orbit solver`
2. 在周期轨上精确评估目标特征 `feature evaluator`

当前代码已经有周期轨提取和 MATCONT 基础，但尚未形成一条“生产级”的统一计算链路。

## 2. 当前代码基础盘点

### 2.1 `PO_extract`

`PO_extract/extract_periodic_orbit.m` 当前已经是一个可复用的周期轨提取入口，但它的核心方法不是“精确周期轨求解”，而是：

1. 先对 ODE 做正向积分
2. 在尾部时间窗上构造观测量
3. 用 Poincare 截面事件检测过零
4. 用相邻截面穿越之间的时间差和返回点差异判断是否收敛到稳定极限环
5. 将最后一个周期从数值轨迹中插值提取出来

代码证据：

- 主入口和主流程见 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:1)
- `backend_used` 被固定写成 `"direct"`，见 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:42)
- 事件检测与收敛判据见 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:115) 和 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:167)
- 轨道提取是对最后一个周期做 `interp1(...,'pchip')` 插值，见 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:282) 和 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:295)
- 振幅 / 极值是对离散采样后的轨道直接取 `max/min`，见 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:269) 和 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:762)

结论：

- 它适合“检测是否存在可吸引的周期轨，并返回一个近似末周期”
- 它不等价于“高精度求解周期边值问题”
- 它当前更像一个 `detector + extractor`，不是严格意义上的 `precision orbit solver`

### 2.2 `matcont`

仓库中已经存在 `MatCont7p6/` 目录，说明 MATCONT 依赖已经放在仓库内。

同时，上层模块会把 `MatCont7p6` 加入路径：

- [flexible_modulators/+flexmod/ensure_paths.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/flexible_modulators/+flexmod/ensure_paths.m:10)

但是，`PO_extract` 当前并没有真正调用 MATCONT 做周期轨精修：

- 选项层面已经预留 `backend`, `matcont_root`, `matcont_odefile`, `matcont_active_parameter`, `matcont_ntst`, `matcont_ncol`, `matcont_tolerance`，见 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:321)
- 诊断里明确写死 `MATCONT refinement is disabled in this implementation.`，见 [extract_periodic_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/extract_periodic_orbit.m:775)

结论：

- `matcont` 目前是“依赖已存在、接口预留、但未真正接入求解流程”
- 这意味着当前仓库已经具备“接入 MATCONT 后端”的良好起点，但没有完成最后一段工程闭环

## 3. 当前上层模块如何使用周期轨

目前多个模块已经把 `PO_extract` 当成统一的周期轨来源：

- `flexible_modulators` 的轨道搜索接口见 [find_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/flexible_modulators/+flexmod/find_orbit.m:1)
- `flexible_modulators` 的状态再精修接口见 [refine_orbit_from_state.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/flexible_modulators/+flexmod/refine_orbit_from_state.m:1)
- `Circadian` 的轨道搜索接口见 [find_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/Circadian/+circadian/find_orbit.m:1)

这些上层模块消费的主要是：

- 一段离散采样的周期轨 `orbit_t, orbit_y`
- 周期 `period`
- 某个选定观测量或状态分量的振幅
- 极值 `yMax / yMin`

这说明一个很关键的事实：

- 现有工程已经隐含形成了“周期轨求解 -> 特征评估 -> 调控算法使用”的数据流
- 但“特征评估”还散落在上层模块里，没有独立成统一的、可验证的科学核心

## 4. 当前精度与适用性评估

### 4.1 已有优点

- 接口已经统一，`extract_periodic_orbit` 可以被多个系统复用
- 支持自定义 observable 和 section
- 有严格 / 宽松两层判据，可区分 `candidate` 与 `converged`
- 输出了较完整的 `diagnostics`，后续适合做误差分析和调参

### 4.2 主要局限

#### A. 当前方法只能抓“可吸引周期轨”

由于它依赖正向积分收敛，所以它本质上只适合稳定极限环。对以下情形能力不足：

- 不稳定周期轨
- 多稳态下很依赖初值的周期轨
- 非常慢收敛的周期轨
- 临界分岔附近、截面回归很慢的周期轨

#### B. 当前“精确特征”并没有误差控制

当前周期和振幅来自：

- 事件时刻差分
- 对离散轨道插值得到的末周期
- 对插值轨道的离散样本做 `max/min`

这可以得到“可用近似值”，但还不能称为“精确特征计算”，因为目前没有明确给出：

- 周期误差估计
- 极值误差估计
- observable 采样误差估计
- 截面选择对结果的敏感性

#### C. 当前默认容差偏工程可用，不偏高精度

例如当前两个主要上层模块默认使用：

- `RelTol = 1e-6`
- `AbsTol = 1e-9`

代码位置：

- [flexible_modulators/default_config.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/flexible_modulators/default_config.m:54)
- [Circadian/default_config.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/Circadian/default_config.m:14)

这对一般绘图和轨道粗提取通常够用，但如果目标是“精确比较不同调控算法得到的参数对应特征”，后续通常需要更系统的精度分层：

- 快速评估档
- 生产计算档
- 基准验证档

#### D. 当前特征定义没有统一抽象

现在不同模块对 amplitude 的定义实际上是各自写在上层里的：

- `flexible_modulators` 直接取第二个状态分量的半峰峰值，见 [find_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/flexible_modulators/+flexmod/find_orbit.m:49)
- `Circadian` 则先定义 `observable = y(2)+y(3)` 再取振幅，见 [find_orbit.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/Circadian/+circadian/find_orbit.m:36)

这意味着“目标振荡特征”还没有上升为统一的数据模型。

## 5. 当前代码健康度评估

### 5.1 测试的正面信号

`PO_extract/test_unified_periodic_orbit_entrypoint.m` 已经覆盖了这些场景：

- 收敛极限环
- 非振荡平衡态
- 时间窗过短
- 带漂移系统
- 慢收敛候选周期轨
- 输入校验
- 请求 MATCONT 但回退 direct 的行为

见 [test_unified_periodic_orbit_entrypoint.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/test_unified_periodic_orbit_entrypoint.m:1)

这说明 `PO_extract` 已经不是一次性脚本，而是正在向稳定模块演化。

### 5.2 已确认的问题

我在 2026-04-16 直接运行了：

```matlab
/Applications/MATLAB_R2024b.app/bin/matlab -batch "cd('/Users/caiyutong/Desktop/CYT/Codes/FMAM_code'); addpath('PO_extract'); test_unified_periodic_orbit_entrypoint"
```

结果测试失败，报错为：

```text
函数或变量 'POinfo' 无法识别。
出错 test_unified_periodic_orbit_entrypoint (第 87 行)
```

对应代码位置：

- [test_unified_periodic_orbit_entrypoint.m](/Users/caiyutong/Desktop/CYT/Codes/FMAM_code/PO_extract/test_unified_periodic_orbit_entrypoint.m:87)

当前 `PO_extract/` 目录下只有两个文件：

- `extract_periodic_orbit.m`
- `test_unified_periodic_orbit_entrypoint.m`

说明现状至少存在一个明确断裂：

- 测试仍然依赖历史兼容接口 `POinfo`
- 但仓库中已经没有 `POinfo.m`

这会影响后续所有“基于现有基础继续开发”的可信度，因为测试链没有完全打通。

后续决策：

- 第一阶段正式废弃 `POinfo`
- `PO_extract` 的唯一受支持公共入口统一为 `extract_periodic_orbit`
- 旧的 warning 语义不再作为受支持 contract；状态表达统一迁移到结构化返回值

### 5.3 已验证的运行信号

我额外用 Van der Pol 例子直接调用了 `extract_periodic_orbit`。当前实现可以返回：

- `status = "converged_periodic_orbit"`
- `code = 2`
- `period ≈ 6.6633`
- `backend_used = "direct"`
- `diagnostics.matcont.status = "disabled"`

这说明：

- direct 后端当前是能工作的
- MATCONT 后端当前确实没有接入

## 6. 当前阶段的工程判断

截至 2026-04-16，我对现有基础的判断是：

1. `PO_extract` 已经能承担“稳定周期轨检测与末周期提取”的职责
2. 它还不能单独承担“精确振荡特征计算核心”的职责
3. `matcont` 是实现高精度周期轨能力的最有价值现成基础
4. 当前最缺的不是更多脚本，而是统一的科学接口、精度分层和验证基准

换句话说，现有代码不是从零开始，但也还没有到可以直接作为“精确计算平台”的状态。

## 7. 建议的完整开发路径（第一版）

下面给出当前最合理的一版开发路径。后续如果目标定义变化，需要修订。

### 阶段 0. 先把问题定义收紧

先明确以下工程输入输出：

- 输入：系统 ODE、参数向量、初值 / 种子轨道、需要计算的特征列表
- 输出：周期轨、特征值、误差估计、求解状态、诊断信息
- 作用范围：先只覆盖稳定周期轨，还是一开始就要求支持不稳定周期轨

建议先把第一版目标限定为：

- 面向稳定周期轨
- 支持多个观测量
- 输出周期、振幅、最大值、最小值、均值、峰时刻
- 每个特征都附带误差或可信度指标

### 阶段 1. 把科学核心拆成两个模块

建议把未来核心明确拆成：

1. `periodic_orbit_solver`
2. `oscillation_feature_evaluator`

理想接口形态：

```matlab
orbit = solve_periodic_orbit(system, params, seed, options)
features = evaluate_oscillation_features(orbit, observables, featureSpec, options)
```

这样做的好处是：

- 周期轨求解与特征评估解耦
- direct / MATCONT 可以共享同一特征评估层
- 后续 FMAM、Circadian、network_modulatability 都能复用同一科学核心

### 阶段 2. 先把 `PO_extract` 稳定成可靠 direct 后端

这一步不求最高精度，目标是打好统一接口基础。

优先事项：

- 正式废弃 `POinfo`，让测试重新全绿，并统一迁移到 `extract_periodic_orbit`
- 明确 `extract_periodic_orbit` 的 contract
- 增加误差相关 diagnostics
- 统一返回结构，让上层不再各自重复计算 amplitude / extrema

这一步完成后，`PO_extract` 应该被重命名理解为：

- `direct periodic orbit backend`

而不是最终的高精度后端。

### 阶段 3. 真正接入 MATCONT 周期轨后端

这是本工程能否达到“精确计算”目标的关键阶段。

目标不是简单“调用 MATCONT”，而是把 MATCONT 纳入统一 contract：

- direct 后端负责找初始近似轨
- MATCONT 后端负责做周期轨精修
- 两者输出统一的 orbit struct

至少需要解决：

- ODE 文件 / active parameter 映射
- 初值轨道到 MATCONT 周期解初始化的转换
- 网格参数 `ntst`, `ncol` 的配置
- MATCONT 失败时的回退策略
- 求解后轨道重采样与统一格式化

### 阶段 4. 建立统一的振荡特征库

把现在散落在各模块里的特征计算统一抽成单独模块，最少先做：

- period
- observable amplitude
- observable max
- observable min
- state-wise amplitude / max / min
- mean value
- peak time / trough time

并支持：

- 自定义 observable
- 多个 observable 并行评估
- 统一命名和输出 schema

### 阶段 5. 建立“精度声明”和基准验证

这是研究代码走向可发表、可对比、可复现实验平台的关键。

建议建立三类基准：

- 有解析或半解析参考值的标准振子，例如 Van der Pol、Hopf normal form
- direct 与 MATCONT 的交叉比较基准
- 同一参数下容差、采样密度、截面选择变化的敏感性基准

最终输出应能回答：

- 当前 feature 值是多少
- 这个值的数值可信度如何
- 它对求解容差 / 轨道采样 / observable 定义有多敏感

### 阶段 6. 再把上层算法逐步迁移到统一核心

等核心稳定后，再把：

- `flexible_modulators`
- `Circadian`
- `network_modulatability`
- 以及后续调控算法

逐步迁移到统一的 `orbit + feature` 计算接口上。

这一步不要过早做，否则会把未稳定的数值核心扩散到整个仓库。

## 8. 我建议的近期任务顺序

如果按最稳妥的路线推进，下一轮建议按下面顺序做：

1. 修复 `PO_extract` 测试断裂，正式移除 `POinfo` 旧接口依赖
2. 为 `extract_periodic_orbit` 写一份正式 I/O contract
3. 设计统一的 `orbit struct` 和 `feature struct`
4. 在 `PO_extract` 内先实现一个独立的 feature evaluator
5. 评估 MATCONT 周期轨初始化接口，做第一版最小接入
6. 建立 direct vs MATCONT 的基准案例

## 9. 当前结论

当前仓库最有价值的现成资产有两个：

- `PO_extract` 已经提供了统一的 direct 检测接口
- `MatCont7p6` 已经在仓库内，可作为高精度周期轨后端

但要把它们变成“精确计算目标振荡特征”的生产级代码库，还需要补上三件事：

- 统一科学接口
- 真正接入 MATCONT
- 建立精度验证体系

在这三件事完成前，当前代码更适合“提取近似周期轨并做工程可用的特征估计”，还不适合直接宣称“精确计算平台”。

## 10. 后续维护约定

后续如果继续推进本工程，建议持续在本文件追加：

- 新增模块
- 已完成里程碑
- 基准测试结果
- 数值误差评估
- 关键设计决策

这样可以把 `PO_extract` 逐步变成这个工程的真实中枢，而不只是一个脚本目录。
