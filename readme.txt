%%%%%%%%%%%%%%%%%%%%%%%%非线性最优化问题主要算法Matlab程序设计%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%线搜索技术%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1.golds.m 【0.618法程序】 用0.618法求单变量函数在单峰区间[a,b]上的近似极小点。
2.qmin.m 【抛物线算法程序】 求函数在区间[a,b]上的局部最小值，从初始点s开始，然后在[a,s],[s,b]上进行搜索。
3.armijo.m 【Armijo准则程序】 Armijo搜索规则是许多非线性优化算法都必须执行的步骤，把它编制成可重复利用的程序模块是很有意义的。

%%%%%%%%%%%%%%%%%%%%%%%%最速下降法及牛顿法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
4.grad.m 【最速下降法程序】 基于Armijo非精确线搜索的最速下降法Matlab程序。
5.dampnm.m 【阻尼牛顿法程序】 基于 Armijo非精确线搜索的阻尼牛顿法 Matlab 程序。
6.revisenm.m  【修正牛顿法程序】修正牛顿法克服了牛顿法要求Hesse 阵正定的缺陷。

%%%%%%%%%%%%%%%%%%%%%%%%共轭梯度法法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
7.frcg.m 【FR共轭梯度法程序】 基于Armijo非精确线搜索的再开始FR共轭梯度法的Matlab程序。

%%%%%%%%%%%%%%%%%%%%%%%%拟牛顿法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
8.sr1.m 【对称秩 1 算法程序】 基于Armijo搜索的对称秩1算法的Matlab程序。
9.bfgs.m 【BFGS 算法程序】 基于Armijo搜索的BFGS算法的Matlab程序。 
10.dfp.m 【DFP算法程序】 基于Armijo搜索的DFP算法的Matlab程序。
11.broyden.m  【Broyden族算法程序】 基于 Armijo 搜索的 Broyden 族算法的 Matlab 程序。

%%%%%%%%%%%%%%%%%%%%%%%%信赖域方法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
12.trustq.m 【求解子问题的光滑牛顿法程序1】 利用光滑牛顿法求解信赖域子问题, 一般适用于 (近似) Hesse 阵正定的情形。
13.trustm.m 【牛顿型信赖域方法程序】

%%%%%%%%%%%%%%%%%%%%%%%%非线性最小二乘问题%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
14.lmm.m 【L-M算法程序】 利用 LM 方法求解非线性方程组F(x)= 0, 可适用于未知数的个数与方程的个数不相等的情形。

%%%%%%%%%%%%%%%%%%%%%%%%罚函数法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
15.multphr.m 【PHR乘子法程序】 可利用该程序求解约束优化问题。

%%%%%%%%%%%%%%%%%%%%%%%%二次规划法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
16.qlag.m 【Lagrange算法程序】 用拉格朗日方法求解等式约束条件的二次规划问题。
17.qpact.m 【有效集方法程序】 主要适用于求解一般约束条件下的凸二次规划问题。

%%%%%%%%%%%%%%%%%%%%%%%%序列二次规划法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
18.qpsubp.m 【求解子问题的光滑牛顿法程序2 】利用光滑牛顿法求解二次规划子问题。
19.sqpm.m 【SQP 方法程序】求解一般约束优化问题，该程序在每一迭代步调用了程序求解二次规划子问题.
