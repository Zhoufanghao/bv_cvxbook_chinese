# 第 4 章 Convex optimization problems

## 4.1 优化问题

### 4.1.1 基本术语

我们用记号

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & f_i(x)\le 0,\quad i=1,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p
\end{array}
\tag{4.1}
$$

来描述这样的问题：在所有满足条件 $f_i(x)\le 0,\ i=1,\ldots,m$ 和 $h_i(x)=0,\ i=1,\ldots,p$ 的 $x$ 中，寻找使 $f_0(x)$ 最小的 $x$。我们称 $x\in\mathbf{R}^n$ 为优化变量，称函数 $f_0:\mathbf{R}^n\to\mathbf{R}$ 为目标函数或成本函数。不等式 $f_i(x)\le0$ 称为不等式约束，相应函数 $f_i:\mathbf{R}^n\to\mathbf{R}$ 称为不等式约束函数。等式 $h_i(x)=0$ 称为等式约束，函数 $h_i:\mathbf{R}^n\to\mathbf{R}$ 称为等式约束函数。如果没有约束（即 $m=p=0$），则称问题 (4.1) 为无约束问题。

目标函数和所有约束函数都有定义的点集

$$
D=\bigcap_{i=0}^m \operatorname{dom} f_i
\cap
\bigcap_{i=1}^p \operatorname{dom} h_i
$$

称为优化问题 (4.1) 的定义域。若点 $x\in D$ 满足约束 $f_i(x)\le0,\ i=1,\ldots,m$ 和 $h_i(x)=0,\ i=1,\ldots,p$，则称 $x$ 是可行的。如果至少存在一个可行点，则称问题 (4.1) 可行；否则称其不可行。所有可行点组成的集合称为可行集或约束集。

问题 (4.1) 的最优值 $p^\star$ 定义为

$$
p^\star=
\inf\{f_0(x)\mid f_i(x)\le0,\ i=1,\ldots,m,\ h_i(x)=0,\ i=1,\ldots,p\}.
$$

我们允许 $p^\star$ 取扩展值 $\pm\infty$。若问题不可行，则 $p^\star=\infty$（遵循空集下确界为 $\infty$ 的标准约定）。若存在可行点 $x_k$ 使得 $f_0(x_k)\to-\infty$（当 $k\to\infty$），则 $p^\star=-\infty$，并称问题 (4.1) 下无界。

#### 最优点与局部最优点

若 $x^\star$ 可行且 $f_0(x^\star)=p^\star$，则称 $x^\star$ 是最优点，或称它求解了问题 (4.1)。所有最优点的集合称为最优集，记为

$$
X_\mathrm{opt}=
\{x\mid f_i(x)\le0,\ i=1,\ldots,m,\ h_i(x)=0,\ i=1,\ldots,p,\ f_0(x)=p^\star\}.
$$

若问题 (4.1) 存在最优点，则称最优值被达到，或称问题可解。若 $X_\mathrm{opt}$ 为空，则称最优值未被达到。（当问题下无界时总会出现这种情形。）若可行点 $x$ 满足 $f_0(x)\le p^\star+\epsilon$，其中 $\epsilon>0$，则称 $x$ 是 $\epsilon$-次优的；所有 $\epsilon$-次优点组成的集合称为问题 (4.1) 的 $\epsilon$-次优集。

如果可行点 $x$ 满足存在 $R>0$ 使得

$$
f_0(x)=
\inf\{f_0(z)\mid f_i(z)\le0,\ i=1,\ldots,m,\ h_i(z)=0,\ i=1,\ldots,p,\ \|z-x\|_2\le R\},
$$

则称 $x$ 是局部最优的。换言之，$x$ 在可行集中邻近它的点上最小化 $f_0$。有时使用“全局最优”来区别“局部最优”和“最优”。不过在本书中，最优总是指全局最优。

如果 $x$ 可行且 $f_i(x)=0$，则称第 $i$ 个不等式约束 $f_i(x)\le0$ 在 $x$ 处活跃。如果 $f_i(x)<0$，则称该约束在 $x$ 处不活跃。（等式约束在所有可行点处都是活跃的。）如果删除某个约束不改变可行集，则称该约束是冗余的。

**例 4.1** 下面用几个简单的无约束优化问题说明这些定义，其中变量 $x\in\mathbf{R}$，且 $\operatorname{dom}f_0=\mathbf{R}_{++}$。

- $f_0(x)=1/x$：$p^\star=0$，但最优值未被达到。
- $f_0(x)=-\log x$：$p^\star=-\infty$，所以该问题下无界。
- $f_0(x)=x\log x$：$p^\star=-1/e$，在唯一最优点 $x^\star=1/e$ 处达到。

#### 可行性问题

如果目标函数恒为零，则最优值要么为零（若可行集非空），要么为 $\infty$（若可行集为空）。我们称这类问题为可行性问题，有时写作

$$
\begin{array}{ll}
\text{find} & x\\
\text{subject to} & f_i(x)\le0,\quad i=1,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p.
\end{array}
$$

因此，可行性问题就是判断这些约束是否相容；若相容，则找到一个满足它们的点。

### 4.1.2 用标准形式表示问题

我们把 (4.1) 称为标准形式的优化问题。在标准形式中，我们约定不等式和等式约束的右端为零。这总可以通过减去非零右端来实现：例如，等式约束 $g_i(x)=\tilde g_i(x)$ 可表示为 $h_i(x)=0$，其中 $h_i(x)=g_i(x)-\tilde g_i(x)$。类似地，形如 $f_i(x)\ge0$ 的不等式可表示为 $-f_i(x)\le0$。

**例 4.2 盒约束。** 考虑优化问题

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & l_i\le x_i\le u_i,\quad i=1,\ldots,n,
\end{array}
$$

其中 $x\in\mathbf{R}^n$ 是变量。这些约束称为变量界（因为它们给出每个 $x_i$ 的下界和上界）或盒约束（因为可行集是一个盒子）。

它可以写成标准形式：

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & l_i-x_i\le0,\quad i=1,\ldots,n\\
& x_i-u_i\le0,\quad i=1,\ldots,n.
\end{array}
$$

因此有 $2n$ 个不等式约束函数：

$$
f_i(x)=l_i-x_i,\quad i=1,\ldots,n,
$$

以及

$$
f_i(x)=x_{i-n}-u_{i-n},\quad i=n+1,\ldots,2n.
$$

#### 最大化问题

按惯例，我们集中讨论最小化问题。最大化问题

$$
\begin{array}{ll}
\text{maximize} & f_0(x)\\
\text{subject to} & f_i(x)\le0,\quad i=1,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p
\end{array}
\tag{4.2}
$$

可以通过在同一约束下最小化 $-f_0$ 来求解。借助这种对应关系，上述所有术语都可以为最大化问题 (4.2) 定义。例如，(4.2) 的最优值定义为

$$
p^\star=
\sup\{f_0(x)\mid f_i(x)\le0,\ i=1,\ldots,m,\ h_i(x)=0,\ i=1,\ldots,p\},
$$

而可行点 $x$ 是 $\epsilon$-次优的，当且仅当 $f_0(x)\ge p^\star-\epsilon$。在最大化问题中，目标有时称为效用或满意程度，而不是成本。

### 4.1.3 等价问题

本书将以非正式方式使用优化问题的等价概念。如果可以很容易地从一个问题的解得到另一个问题的解，反之亦然，则称两个问题等价。（给出等价的形式化定义是可能的，但会比较复杂。）

一个简单例子是

$$
\begin{array}{ll}
\text{minimize} & \tilde f_0(x)=\alpha_0 f_0(x)\\
\text{subject to} & \tilde f_i(x)=\alpha_i f_i(x)\le0,\quad i=1,\ldots,m\\
& \tilde h_i(x)=\beta_i h_i(x)=0,\quad i=1,\ldots,p,
\end{array}
\tag{4.3}
$$

其中 $\alpha_i>0,\ i=0,\ldots,m$，且 $\beta_i\ne0,\ i=1,\ldots,p$。这个问题由标准形式问题 (4.1) 通过用正数缩放目标函数和不等式约束函数、用非零常数缩放等式约束函数得到。因此，问题 (4.3) 与原问题 (4.1) 的可行集相同。点 $x$ 是原问题 (4.1) 的最优点，当且仅当它是缩放后问题 (4.3) 的最优点，所以我们称这两个问题等价。不过，除非所有 $\alpha_i$ 和 $\beta_i$ 都等于 1，否则两个问题并不相同，因为它们的目标函数和约束函数不同。

#### 变量变换

设 $\phi:\mathbf{R}^n\to\mathbf{R}^n$ 是一一映射，并且其像覆盖问题定义域 $D$，即 $\phi(\operatorname{dom}\phi)\supseteq D$。定义函数

$$
\tilde f_i(z)=f_i(\phi(z)),\quad i=0,\ldots,m,\qquad
\tilde h_i(z)=h_i(\phi(z)),\quad i=1,\ldots,p.
$$

考虑变量为 $z$ 的问题

$$
\begin{array}{ll}
\text{minimize} & \tilde f_0(z)\\
\text{subject to} & \tilde f_i(z)\le0,\quad i=1,\ldots,m\\
& \tilde h_i(z)=0,\quad i=1,\ldots,p.
\end{array}
\tag{4.4}
$$

我们称标准形式问题 (4.1) 与问题 (4.4) 由变量变换或变量代换 $x=\phi(z)$ 联系起来。两个问题显然等价：若 $x$ 求解 (4.1)，则 $z=\phi^{-1}(x)$ 求解 (4.4)；若 $z$ 求解 (4.4)，则 $x=\phi(z)$ 求解 (4.1)。

#### 目标函数和约束函数变换

设 $\psi_0:\mathbf{R}\to\mathbf{R}$ 单调递增，$\psi_1,\ldots,\psi_m:\mathbf{R}\to\mathbf{R}$ 满足 $\psi_i(u)\le0$ 当且仅当 $u\le0$，并且 $\psi_{m+1},\ldots,\psi_{m+p}:\mathbf{R}\to\mathbf{R}$ 满足 $\psi_i(u)=0$ 当且仅当 $u=0$。定义复合函数

$$
\tilde f_i(x)=\psi_i(f_i(x)),\quad i=0,\ldots,m,\qquad
\tilde h_i(x)=\psi_{m+i}(h_i(x)),\quad i=1,\ldots,p.
$$

相关问题与标准形式问题 (4.1) 等价；事实上，它们的可行集相同，最优点也相同。前面的缩放例子 (4.3) 是所有 $\psi_i$ 都为线性函数的特例。

**例 4.3 最小范数与最小范数平方问题。** 考虑无约束欧几里得范数最小化问题

$$
\text{minimize}\quad \|Ax-b\|_2,
\tag{4.5}
$$

变量为 $x\in\mathbf{R}^n$。由于范数总是非负的，也可以求解

$$
\text{minimize}\quad \|Ax-b\|_2^2=(Ax-b)^T(Ax-b).
\tag{4.6}
$$

问题 (4.5) 和 (4.6) 显然等价；它们的最优点相同。但二者并不相同。例如，(4.5) 中的目标在任何满足 $Ax-b=0$ 的 $x$ 处不可微，而 (4.6) 中的目标对所有 $x$ 都可微（事实上是二次函数）。

#### 松弛变量

一个简单变换基于如下观察：$f_i(x)\le0$ 当且仅当存在 $s_i\ge0$ 使得 $f_i(x)+s_i=0$。由此得到问题

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & s_i\ge0,\quad i=1,\ldots,m\\
& f_i(x)+s_i=0,\quad i=1,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p,
\end{array}
\tag{4.7}
$$

变量为 $x\in\mathbf{R}^n$ 和 $s\in\mathbf{R}^m$。该问题有 $n+m$ 个变量、$m$ 个不等式约束（$s_i$ 的非负性约束）以及 $m+p$ 个等式约束。新变量 $s_i$ 称为原不等式约束 $f_i(x)\le0$ 对应的松弛变量。引入松弛变量会把每个不等式约束替换为一个等式约束和一个非负性约束。

问题 (4.7) 与原标准形式问题 (4.1) 等价。若 $(x,s)$ 对 (4.7) 可行，则 $x$ 对原问题可行，因为 $s_i=-f_i(x)\ge0$。反过来，若 $x$ 对原问题可行，则取 $s_i=-f_i(x)$ 时 $(x,s)$ 对 (4.7) 可行。同样，$x$ 是原问题 (4.1) 的最优点，当且仅当 $(x,s)$ 是问题 (4.7) 的最优点，其中 $s_i=-f_i(x)$。

#### 消去等式约束

如果能用参数 $z\in\mathbf{R}^k$ 显式参数化所有等式约束

$$
h_i(x)=0,\quad i=1,\ldots,p
\tag{4.8}
$$

的解，则可以从问题中消去这些等式约束。设函数 $\phi:\mathbf{R}^k\to\mathbf{R}^n$ 满足：$x$ 满足 (4.8) 当且仅当存在某个 $z\in\mathbf{R}^k$ 使得 $x=\phi(z)$。则优化问题

$$
\begin{array}{ll}
\text{minimize} & \tilde f_0(z)=f_0(\phi(z))\\
\text{subject to} & \tilde f_i(z)=f_i(\phi(z))\le0,\quad i=1,\ldots,m
\end{array}
$$

与原问题 (4.1) 等价。变换后的问题变量为 $z\in\mathbf{R}^k$，有 $m$ 个不等式约束，没有等式约束。若 $z$ 对变换后问题最优，则 $x=\phi(z)$ 对原问题最优。反过来，若 $x$ 对原问题最优，则由于 $x$ 可行，至少存在一个 $z$ 使得 $x=\phi(z)$；任何这样的 $z$ 都对变换后问题最优。

#### 消去线性等式约束

当等式约束全是线性的，即具有形式 $Ax=b$ 时，消去变量的过程可以更显式地描述，并且容易数值执行。如果 $Ax=b$ 不相容，即 $b\notin\mathcal{R}(A)$，则原问题不可行。假设不是这种情况，令 $x_0$ 为等式约束的任意一个解。令 $F\in\mathbf{R}^{n\times k}$ 是任意满足 $\mathcal{R}(F)=\mathcal{N}(A)$ 的矩阵，则线性方程 $Ax=b$ 的通解为 $Fz+x_0$，其中 $z\in\mathbf{R}^k$。（可以选择 $F$ 满秩，此时 $k=n-\operatorname{rank}A$。）

将 $x=Fz+x_0$ 代入原问题得到

$$
\begin{array}{ll}
\text{minimize} & f_0(Fz+x_0)\\
\text{subject to} & f_i(Fz+x_0)\le0,\quad i=1,\ldots,m,
\end{array}
$$

变量为 $z$。该问题与原问题等价，没有等式约束，并且变量数减少了 $\operatorname{rank}A$。

#### 引入等式约束

也可以向问题中引入等式约束和新变量。一般情形描述起来复杂且启发性不强，我们只给出一个后面会用到的典型例子。考虑问题

$$
\begin{array}{ll}
\text{minimize} & f_0(A_0x+b_0)\\
\text{subject to} & f_i(A_ix+b_i)\le0,\quad i=1,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p,
\end{array}
$$

其中 $x\in\mathbf{R}^n$，$A_i\in\mathbf{R}^{k_i\times n}$，并且 $f_i:\mathbf{R}^{k_i}\to\mathbf{R}$。在这个问题中，目标和约束函数被给成 $f_i$ 与由 $A_ix+b_i$ 定义的仿射变换的复合。

引入新变量 $y_i\in\mathbf{R}^{k_i}$，以及新等式约束 $y_i=A_ix+b_i$，$i=0,\ldots,m$，可得到等价问题

$$
\begin{array}{ll}
\text{minimize} & f_0(y_0)\\
\text{subject to} & f_i(y_i)\le0,\quad i=1,\ldots,m\\
& y_i=A_ix+b_i,\quad i=0,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p.
\end{array}
$$

这个问题引入 $k_0+\cdots+k_m$ 个新变量 $y_0\in\mathbf{R}^{k_0},\ldots,y_m\in\mathbf{R}^{k_m}$，以及 $k_0+\cdots+k_m$ 个新等式约束 $y_0=A_0x+b_0,\ldots,y_m=A_mx+b_m$。该问题中的目标和不等式约束是相互独立的，即它们涉及不同的优化变量。

#### 对部分变量优化

总有

$$
\inf_{x,y} f(x,y)=\inf_x \tilde f(x),
$$

其中 $\tilde f(x)=\inf_y f(x,y)$。换言之，总可以先对部分变量最小化一个函数，然后再对剩余变量最小化。这个简单而一般的原则可用于把问题变换为等价形式。

设变量 $x\in\mathbf{R}^n$ 被分块为 $x=(x_1,x_2)$，其中 $x_1\in\mathbf{R}^{n_1}$、$x_2\in\mathbf{R}^{n_2}$，且 $n_1+n_2=n$。考虑问题

$$
\begin{array}{ll}
\text{minimize} & f_0(x_1,x_2)\\
\text{subject to} & f_i(x_1)\le0,\quad i=1,\ldots,m_1\\
& \tilde f_i(x_2)\le0,\quad i=1,\ldots,m_2,
\end{array}
\tag{4.9}
$$

其中约束是独立的，也就是说每个约束函数只依赖于 $x_1$ 或 $x_2$。先对 $x_2$ 最小化。定义 $x_1$ 的函数

$$
\tilde f_0(x_1)=
\inf\{f_0(x_1,z)\mid \tilde f_i(z)\le0,\ i=1,\ldots,m_2\}.
$$

于是问题 (4.9) 等价于

$$
\begin{array}{ll}
\text{minimize} & \tilde f_0(x_1)\\
\text{subject to} & f_i(x_1)\le0,\quad i=1,\ldots,m_1.
\end{array}
\tag{4.10}
$$

**例 4.4 对部分变量有约束的二次函数最小化。** 考虑一个具有严格凸二次目标的问题，其中部分变量无约束：

$$
\begin{array}{ll}
\text{minimize} & x_1^TP_{11}x_1+2x_1^TP_{12}x_2+x_2^TP_{22}x_2\\
\text{subject to} & f_i(x_1)\le0,\quad i=1,\ldots,m,
\end{array}
$$

其中 $P_{11}$ 和 $P_{22}$ 对称。这里可以解析地对 $x_2$ 最小化：

$$
\inf_{x_2}
\left(x_1^TP_{11}x_1+2x_1^TP_{12}x_2+x_2^TP_{22}x_2\right)
=x_1^T(P_{11}-P_{12}P_{22}^{-1}P_{12}^T)x_1.
$$

因此原问题等价于

$$
\begin{array}{ll}
\text{minimize} & x_1^T(P_{11}-P_{12}P_{22}^{-1}P_{12}^T)x_1\\
\text{subject to} & f_i(x_1)\le0,\quad i=1,\ldots,m.
\end{array}
$$

#### 上图形式问题

标准问题 (4.1) 的上图形式是

$$
\begin{array}{ll}
\text{minimize} & t\\
\text{subject to} & f_0(x)-t\le0\\
& f_i(x)\le0,\quad i=1,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p,
\end{array}
\tag{4.11}
$$

变量为 $x\in\mathbf{R}^n$ 和 $t\in\mathbf{R}$。它与原问题等价：$(x,t)$ 对 (4.11) 最优，当且仅当 $x$ 对 (4.1) 最优且 $t=f_0(x)$。注意，上图形式问题的目标函数是变量 $x,t$ 的线性函数。

上图形式问题 (4.11) 可以几何地解释为“图空间” $(x,t)$ 中的优化问题：在满足 $x$ 上约束的同时，在 $f_0$ 的上图上最小化 $t$。这如图 4.1 所示。

**图 4.1** 对无约束问题，上图形式问题的几何解释。问题是在上图（阴影区域）中找到使 $t$ 最小的点，即上图中的“最低”点。最优点为 $(x^\star,t^\star)$。

#### 隐式约束和显式约束

通过 §3.1.2 中已经提到的简单技巧，可以把任何约束隐式地包含在目标函数中，即重新定义其定义域。作为极端例子，标准形式问题可表示为无约束问题

$$
\text{minimize}\quad F(x),
\tag{4.12}
$$

其中 $F$ 定义为 $f_0$，但定义域限制为可行集：

$$
\operatorname{dom}F=
\{x\in\operatorname{dom}f_0\mid f_i(x)\le0,\ i=1,\ldots,m,\ h_i(x)=0,\ i=1,\ldots,p\},
$$

并且对 $x\in\operatorname{dom}F$ 有 $F(x)=f_0(x)$。（等价地，也可以定义 $F(x)$ 在不可行点取值 $\infty$。）问题 (4.1) 与 (4.12) 显然等价：它们有相同的可行集、最优点和最优值。

当然，这种变换只是一个记号技巧。尽管问题 (4.12) 至少名义上是无约束的，把约束隐式化并没有让问题更容易分析或求解。在某些方面，这种变换反而让问题更困难。例如，假设原问题中的目标 $f_0$ 可微，因此其定义域是开集；受限目标函数 $F$ 很可能不可微，因为其定义域很可能不是开集。

反过来，我们也会遇到带隐式约束的问题，并可以把这些约束显式化。作为简单例子，考虑无约束问题

$$
\text{minimize}\quad f(x),
\tag{4.13}
$$

其中

$$
f(x)=
\begin{cases}
x^Tx, & Ax=b,\\
\infty, & \text{otherwise}.
\end{cases}
$$

因此，目标函数在由 $Ax=b$ 定义的仿射集上等于二次型 $x^Tx$，在该仿射集之外为 $\infty$。由于显然只需关注满足 $Ax=b$ 的点，我们说问题 (4.13) 的目标中隐藏了一个隐式等式约束 $Ax=b$。可以把它显式化，形成等价问题

$$
\begin{array}{ll}
\text{minimize} & x^Tx\\
\text{subject to} & Ax=b.
\end{array}
\tag{4.14}
$$

问题 (4.13) 和 (4.14) 显然等价，但并不相同。问题 (4.13) 是无约束的，但其目标函数不可微。问题 (4.14) 有等式约束，但目标函数和约束函数都是可微的。

### 4.1.4 参数描述和 oracle 描述

对于标准形式 (4.1) 的问题，还需要说明目标函数和约束函数如何给出。在许多情况下，这些函数具有某种解析形式或闭式形式，即由包含变量 $x$ 和若干参数的公式或表达式给出。例如，假设目标是二次函数，形式为 $f_0(x)=(1/2)x^TPx+q^Tx+r$。为了指定该目标函数，需要给出系数（也称为问题参数或问题数据）$P\in\mathbf{S}^n$、$q\in\mathbf{R}^n$ 和 $r\in\mathbf{R}$。我们称这为参数化问题描述，因为要指定待求解的具体问题（即问题实例），只需给出目标函数和约束函数表达式中出现的参数值。

在其他情况下，目标和约束函数由 oracle 模型（也称黑箱模型或子程序模型）描述。在 oracle 模型中，我们并不显式知道 $f$，但可以在任意 $x\in\operatorname{dom}f$ 处计算 $f(x)$（通常还可以计算某些导数）。这称为查询 oracle，通常会有某种成本，例如时间成本。我们还会得到一些关于函数的先验信息，例如凸性以及函数值的界。作为 oracle 模型的具体例子，考虑一个无约束问题，目标是最小化函数 $f$。函数值 $f(x)$ 和梯度 $\nabla f(x)$ 由一个子程序计算。可以在任意 $x\in\operatorname{dom}f$ 处调用该子程序，但不能访问其源代码。用参数 $x$ 调用子程序后，会返回 $f(x)$ 和 $\nabla f(x)$。注意，在 oracle 模型中，我们从未真正知道函数；我们只知道已经查询 oracle 的点处的函数值（和某些导数）。（同时也知道一些给定的函数先验信息，例如可微性和凸性。）

在实践中，参数化问题描述和 oracle 问题描述之间的区别并不十分尖锐。如果给定了参数化问题描述，就可以为它构造一个 oracle，在查询时简单计算所需函数和导数。第三部分中研究的大多数算法使用 oracle 模型，但当它们被限制为求解某个具体参数化问题族时，可以做得更高效。

## 4.2 凸优化

### 4.2.1 标准形式的凸优化问题

凸优化问题具有如下形式：

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & f_i(x)\le0,\quad i=1,\ldots,m\\
& a_i^Tx=b_i,\quad i=1,\ldots,p,
\end{array}
\tag{4.15}
$$

其中 $f_0,\ldots,f_m$ 是凸函数。与一般标准形式问题 (4.1) 相比，凸问题有三个额外要求：

- 目标函数必须是凸的；
- 不等式约束函数必须是凸的；
- 等式约束函数 $h_i(x)=a_i^Tx-b_i$ 必须是仿射的。

我们立刻得到一个重要性质：凸优化问题的可行集是凸的，因为它是问题定义域

$$
D=\bigcap_{i=0}^m\operatorname{dom}f_i
$$

与 $m$ 个凸下水平集 $\{x\mid f_i(x)\le0\}$ 和 $p$ 个超平面 $\{x\mid a_i^Tx=b_i\}$ 的交集。（可以不失一般性地假设 $a_i\ne0$：若某个 $i$ 有 $a_i=0$ 且 $b_i=0$，则第 $i$ 个等式约束可删除；若 $a_i=0$ 且 $b_i\ne0$，则第 $i$ 个等式约束不相容，问题不可行。）因此，在凸优化问题中，我们是在凸集上最小化凸目标函数。

如果 $f_0$ 是拟凸而不是凸的，则称问题 (4.15) 为（标准形式的）拟凸优化问题。由于凸函数或拟凸函数的下水平集是凸的，可以得出：对于凸优化或拟凸优化问题，$\epsilon$-次优集是凸的。特别地，最优集是凸的。如果目标函数严格凸，则最优集至多包含一个点。

#### 凹最大化问题

稍微滥用记号，当目标函数 $f_0$ 为凹函数、而不等式约束函数 $f_1,\ldots,f_m$ 为凸函数时，我们也把

$$
\begin{array}{ll}
\text{maximize} & f_0(x)\\
\text{subject to} & f_i(x)\le0,\quad i=1,\ldots,m\\
& a_i^Tx=b_i,\quad i=1,\ldots,p
\end{array}
\tag{4.16}
$$

称为凸优化问题。这个凹最大化问题可以通过最小化凸目标函数 $-f_0$ 直接求解。我们为最小化问题描述的所有结果、结论和算法，都很容易转置到最大化情形。类似地，如果 $f_0$ 是拟凹的，则称最大化问题 (4.16) 为拟凸的。

#### 抽象形式的凸优化问题

需要注意我们对凸优化问题定义中的一个细微之处。考虑 $x\in\mathbf{R}^2$ 的例子

$$
\begin{array}{ll}
\text{minimize} & f_0(x)=x_1^2+x_2^2\\
\text{subject to} & f_1(x)=x_1/(1+x_2^2)\le0\\
& h_1(x)=(x_1+x_2)^2=0.
\end{array}
\tag{4.17}
$$

这个问题是标准形式 (4.1)，但不是标准形式的凸优化问题，因为等式约束函数 $h_1$ 不是仿射的，不等式约束函数 $f_1$ 也不是凸的。然而，可行集 $\{x\mid x_1\le0,\ x_1+x_2=0\}$ 是凸的。因此，虽然这个问题是在凸集上最小化凸函数 $f_0$，但按我们的定义它不是凸优化问题。

当然，该问题很容易改写为

$$
\begin{array}{ll}
\text{minimize} & f_0(x)=x_1^2+x_2^2\\
\text{subject to} & \tilde f_1(x)=x_1\le0\\
& \tilde h_1(x)=x_1+x_2=0,
\end{array}
\tag{4.18}
$$

这就是标准凸优化形式，因为 $f_0$ 和 $\tilde f_1$ 是凸的，$\tilde h_1$ 是仿射的。

有些作者使用抽象凸优化问题来描述“在凸集上最小化凸函数”这个抽象问题。按这种术语，问题 (4.17) 是抽象凸优化问题。本书不会使用这个术语。对我们而言，凸优化问题不仅仅是在凸集上最小化凸函数；还要求可行集具体由一组涉及凸函数的不等式和一组线性等式约束描述。问题 (4.17) 不是凸优化问题，但问题 (4.18) 是凸优化问题。（不过两者等价。）

我们采用较严格的凸优化问题定义，在实践中影响不大。要解决在凸集上最小化凸函数的抽象问题，仍需找到用凸不等式和线性等式约束来描述该集合的方法。正如上例所示，这通常是直接的。

### 4.2.2 局部最优与全局最优

凸优化问题的一个基本性质是：任何局部最优点也是（全局）最优点。为说明这一点，设 $x$ 是某凸优化问题的局部最优点，即 $x$ 可行且

$$
f_0(x)=\inf\{f_0(z)\mid z\ \text{可行},\ \|z-x\|_2\le R\}
\tag{4.19}
$$

对某个 $R>0$ 成立。现在假设 $x$ 不是全局最优，即存在可行点 $y$ 使得 $f_0(y)<f_0(x)$。显然 $\|y-x\|_2>R$，否则会有 $f_0(x)\le f_0(y)$。考虑点

$$
z=(1-\theta)x+\theta y,\qquad
\theta=\frac{R}{2\|y-x\|_2}.
$$

则 $\|z-x\|_2=R/2<R$，并且由于可行集凸，$z$ 可行。由 $f_0$ 的凸性，有

$$
f_0(z)\le(1-\theta)f_0(x)+\theta f_0(y)<f_0(x),
$$

这与 (4.19) 矛盾。因此不存在满足 $f_0(y)<f_0(x)$ 的可行点 $y$，即 $x$ 是全局最优的。

拟凸优化问题的局部最优点不一定是全局最优点；见 §4.2.5。

### 4.2.3 可微 $f_0$ 的最优性准则

假设凸优化问题中的目标 $f_0$ 可微，因此对所有 $x,y\in\operatorname{dom}f_0$ 有

$$
f_0(y)\ge f_0(x)+\nabla f_0(x)^T(y-x)
\tag{4.20}
$$

（见 §3.1.3）。令 $X$ 表示可行集，即

$$
X=\{x\mid f_i(x)\le0,\ i=1,\ldots,m,\ h_i(x)=0,\ i=1,\ldots,p\}.
$$

则 $x$ 最优当且仅当 $x\in X$ 且

$$
\nabla f_0(x)^T(y-x)\ge0\quad \text{对所有 }y\in X.
\tag{4.21}
$$

这个最优性准则可以几何理解：若 $\nabla f_0(x)\ne0$，则它表示 $-\nabla f_0(x)$ 在 $x$ 处定义了可行集的一个支撑超平面（见图 4.2）。

**图 4.2** 最优性条件 (4.21) 的几何解释。可行集 $X$ 显示为阴影区域。$f_0$ 的若干水平曲线显示为虚线。点 $x$ 是最优的：$-\nabla f_0(x)$ 在 $x$ 处定义了 $X$ 的一个支撑超平面（实线）。

#### 最优性条件证明

首先假设 $x\in X$ 并满足 (4.21)。若 $y\in X$，由 (4.20) 有 $f_0(y)\ge f_0(x)$。因此 $x$ 是问题 (4.1) 的最优点。

反过来，设 $x$ 最优，但条件 (4.21) 不成立，即存在某个 $y\in X$ 使得

$$
\nabla f_0(x)^T(y-x)<0.
$$

考虑 $z(t)=ty+(1-t)x$，其中 $t\in[0,1]$。由于 $z(t)$ 在线段 $x$ 与 $y$ 上，且可行集凸，$z(t)$ 可行。我们断言，当 $t$ 为足够小的正数时，有 $f_0(z(t))<f_0(x)$，从而证明 $x$ 不是最优的。事实上，

$$
\left.\frac{d}{dt}f_0(z(t))\right|_{t=0}
=\nabla f_0(x)^T(y-x)<0,
$$

所以对足够小的正 $t$，有 $f_0(z(t))<f_0(x)$。

第 5 章将更深入地研究最优性条件。这里先考察几个简单例子。

#### 无约束问题

对无约束问题（即 $m=p=0$），条件 (4.21) 化为著名的充要条件

$$
\nabla f_0(x)=0
\tag{4.22}
$$

用于判断 $x$ 是否最优。虽然我们已经见过这个最优性条件，但看看它如何由 (4.21) 推出仍有帮助。设 $x$ 最优，也就是说 $x\in\operatorname{dom}f_0$，并且对所有可行 $y$ 有 $\nabla f_0(x)^T(y-x)\ge0$。由于 $f_0$ 可微，其定义域（按定义）是开集，因此所有充分接近 $x$ 的 $y$ 都可行。取 $y=x-t\nabla f_0(x)$，其中 $t\in\mathbf{R}$。当 $t$ 很小且为正时，$y$ 可行，于是

$$
\nabla f_0(x)^T(y-x)=-t\|\nabla f_0(x)\|_2^2\ge0,
$$

从而推出 $\nabla f_0(x)=0$。

根据方程 (4.22) 的解的数量，可能出现几种情形。如果 (4.22) 无解，则没有最优点；问题的最优值未被达到。这里可以区分两种情况：问题下无界，或者最优值有限但未达到。另一方面，方程 (4.22) 可以有多个解，此时每个解都是 $f_0$ 的极小点。

**例 4.5 无约束二次优化。** 考虑最小化二次函数

$$
f_0(x)=(1/2)x^TPx+q^Tx+r,
$$

其中 $P\in\mathbf{S}_+^n$（因此 $f_0$ 凸）。$x$ 是 $f_0$ 极小点的充要条件是

$$
\nabla f_0(x)=Px+q=0.
$$

根据这个线性方程无解、唯一解或多解，会出现几种情况。

- 如果 $q\notin\mathcal{R}(P)$，则无解。此时 $f_0$ 下无界。
- 如果 $P\succ0$（这是 $f_0$ 严格凸的条件），则存在唯一极小点 $x^\star=-P^{-1}q$。
- 如果 $P$ 奇异但 $q\in\mathcal{R}(P)$，则最优点集合是仿射集 $X_\mathrm{opt}=-P^\dagger q+\mathcal{N}(P)$，其中 $P^\dagger$ 表示 $P$ 的伪逆（见 §A.5.4）。

**例 4.6 解析中心。** 考虑最小化凸函数 $f_0:\mathbf{R}^n\to\mathbf{R}$ 的无约束问题，其中

$$
f_0(x)=-\sum_{i=1}^m\log(b_i-a_i^Tx),
\qquad
\operatorname{dom}f_0=\{x\mid Ax\prec b\},
$$

$a_1^T,\ldots,a_m^T$ 是 $A$ 的各行。函数 $f_0$ 可微，所以 $x$ 最优的充要条件为

$$
Ax\prec b,\qquad
\nabla f_0(x)=\sum_{i=1}^m\frac{1}{b_i-a_i^Tx}a_i=0.
\tag{4.23}
$$

条件 $Ax\prec b$ 只是 $x\in\operatorname{dom}f_0$。如果 $Ax\prec b$ 不可行，则 $f_0$ 的定义域为空。假设 $Ax\prec b$ 可行，仍有几种可能（见习题 4.2）：

- (4.23) 无解，因此该问题没有最优点。这当且仅当 $f_0$ 下无界时发生。
- (4.23) 有多个解。此时可以证明这些解构成一个仿射集。
- (4.23) 有唯一解，即 $f_0$ 有唯一极小点。这当且仅当开多面体 $\{x\mid Ax\prec b\}$ 非空且有界时发生。

#### 只有等式约束的问题

考虑有等式约束但没有不等式约束的情形：

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & Ax=b.
\end{array}
$$

这里可行集是仿射的。假设它非空；否则问题不可行。可行点 $x$ 的最优性条件是

$$
\nabla f_0(x)^T(y-x)\ge0
$$

对所有满足 $Ay=b$ 的 $y$ 成立。由于 $x$ 可行，每个可行 $y$ 都可写为 $y=x+v$，其中 $v\in\mathcal{N}(A)$。因此最优性条件可表示为

$$
\nabla f_0(x)^Tv\ge0\quad \text{对所有 }v\in\mathcal{N}(A).
$$

若一个线性函数在子空间上非负，则它必须在该子空间上为零，所以 $\nabla f_0(x)^Tv=0$ 对所有 $v\in\mathcal{N}(A)$ 成立。换言之，

$$
\nabla f_0(x)\perp\mathcal{N}(A).
$$

利用 $\mathcal{N}(A)^\perp=\mathcal{R}(A^T)$，该最优性条件可表示为 $\nabla f_0(x)\in\mathcal{R}(A^T)$，即存在 $\nu\in\mathbf{R}^p$ 使得

$$
\nabla f_0(x)+A^T\nu=0.
$$

再加上要求 $Ax=b$（即 $x$ 可行），这就是经典的 Lagrange 乘子最优性条件，我们将在第 5 章更详细地研究它。

#### 在非负正交锥上最小化

再考虑问题

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & x\succeq0,
\end{array}
$$

其中唯一的不等式约束是变量的非负性约束。

最优性条件 (4.21) 为

$$
x\succeq0,\qquad
\nabla f_0(x)^T(y-x)\ge0\quad \text{对所有 }y\succeq0.
$$

项 $\nabla f_0(x)^Ty$ 是 $y$ 的线性函数，除非 $\nabla f_0(x)\succeq0$，否则它在 $y\succeq0$ 上下无界。于是条件化为 $-\nabla f_0(x)^Tx\ge0$。但 $x\succeq0$ 且 $\nabla f_0(x)\succeq0$，所以必须有 $\nabla f_0(x)^Tx=0$，即

$$
\sum_{i=1}^n(\nabla f_0(x))_i x_i=0.
$$

该求和中每一项都是两个非负数的乘积，所以每一项都必须为零，即 $(\nabla f_0(x))_i x_i=0$，$i=1,\ldots,n$。因此最优性条件可表示为

$$
x\succeq0,\qquad
\nabla f_0(x)\succeq0,\qquad
x_i(\nabla f_0(x))_i=0,\quad i=1,\ldots,n.
$$

最后一个条件称为互补性，因为它意味着向量 $x$ 和 $\nabla f_0(x)$ 的稀疏模式（即非零分量对应的指标集合）是互补的（交集为空）。我们将在第 5 章再次遇到互补性条件。

### 4.2.4 等价凸问题

下面考察 §4.1.3 中描述的哪些变换会保持凸性。

#### 消去等式约束

对凸问题，等式约束必须是线性的，即形如 $Ax=b$。这时可以通过寻找一个特解 $x_0$ 以及一个值域为 $A$ 的零空间的矩阵 $F$ 来消去等式约束，得到问题

$$
\begin{array}{ll}
\text{minimize} & f_0(Fz+x_0)\\
\text{subject to} & f_i(Fz+x_0)\le0,\quad i=1,\ldots,m,
\end{array}
$$

变量为 $z$。由于凸函数与仿射函数的复合仍然凸，消去等式约束保持问题的凸性。此外，消去等式约束以及从变换后问题的解重构原问题解的过程，都只涉及标准线性代数运算。

至少在原则上，这意味着我们可以把注意力限制在没有等式约束的凸优化问题上。不过在许多情况下，保留等式约束更好，因为消去它们可能使问题更难理解和分析，或破坏求解算法的效率。例如，当变量 $x$ 维度很大时，消去等式约束可能破坏问题的稀疏性或其他有用结构。

#### 引入等式约束

可以在凸优化问题中引入新变量和等式约束，只要等式约束是线性的，则所得问题仍是凸的。例如，如果目标或约束函数具有形式 $f_i(A_ix+b_i)$，其中 $A_i\in\mathbf{R}^{k_i\times n}$，则可以引入新变量 $y_i\in\mathbf{R}^{k_i}$，用 $f_i(y_i)$ 替换 $f_i(A_ix+b_i)$，并加入线性等式约束 $y_i=A_ix+b_i$。

#### 松弛变量

引入松弛变量会带来新约束 $f_i(x)+s_i=0$。由于凸问题中的等式约束函数必须是仿射的，所以必须要求 $f_i$ 仿射。换言之：为线性不等式引入松弛变量会保持问题的凸性。

#### 上图形式问题

凸优化问题 (4.15) 的上图形式为

$$
\begin{array}{ll}
\text{minimize} & t\\
\text{subject to} & f_0(x)-t\le0\\
& f_i(x)\le0,\quad i=1,\ldots,m\\
& a_i^Tx=b_i,\quad i=1,\ldots,p.
\end{array}
$$

目标是线性的（因而是凸的），新约束函数 $f_0(x)-t$ 在 $(x,t)$ 中也是凸的，所以该上图形式问题也是凸的。

有时会说线性目标对凸优化是“通用的”，因为任何凸优化问题都可以很容易地变换为具有线性目标的问题。凸问题的上图形式有若干实际用途。假设凸优化问题的目标是线性的，可以简化理论分析。它也可以简化算法开发，因为一个能够求解线性目标凸优化问题的算法，通过上述变换就能求解任何凸优化问题（只要它能处理约束 $f_0(x)-t\le0$）。

#### 对部分变量最小化

对部分变量最小化凸函数保持凸性。因此，如果 (4.9) 中的 $f_0$ 关于 $x_1$ 和 $x_2$ 联合凸，并且 $f_i,\ i=1,\ldots,m_1$ 与 $\tilde f_i,\ i=1,\ldots,m_2$ 都是凸的，则等价问题 (4.10) 是凸的。

### 4.2.5 拟凸优化

回忆拟凸优化问题具有标准形式

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & f_i(x)\le0,\quad i=1,\ldots,m\\
& Ax=b,
\end{array}
\tag{4.24}
$$

其中不等式约束函数 $f_1,\ldots,f_m$ 是凸的，而目标 $f_0$ 是拟凸的（不像凸优化问题中那样要求凸）。（拟凸约束函数可以替换为等价的凸约束函数，即具有相同 0-下水平集的凸约束函数，见 §3.4.5。）

本节指出凸优化和拟凸优化之间的一些基本差别，并说明如何把求解拟凸优化问题化为求解一系列凸优化问题。

#### 局部最优解和最优性条件

凸优化与拟凸优化之间最重要的差别是：拟凸优化问题可以有非全局最优的局部最优解。即使在 $\mathbf{R}$ 上无约束最小化拟凸函数的简单情形中，也会出现这种现象，如图 4.3 所示。

不过，对于目标函数可微的拟凸优化问题，§4.2.3 中最优性条件 (4.21) 的一个变体仍然成立。令 $X$ 表示拟凸优化问题 (4.24) 的可行集。由拟凸性的一阶条件 (3.20) 可知，如果

$$
x\in X,\qquad
\nabla f_0(x)^T(y-x)>0\quad \text{对所有 }y\in X\setminus\{x\},
\tag{4.25}
$$

则 $x$ 是最优的。这个准则与凸优化中类似准则 (4.21) 有两个重要差别：

- 条件 (4.25) 只是最优性的充分条件；简单例子表明，最优点不一定满足它。相比之下，条件 (4.21) 是 $x$ 求解凸问题的充要条件。
- 条件 (4.25) 要求 $f_0$ 的梯度非零，而条件 (4.21) 不要求。事实上，在凸情形中，当 $\nabla f_0(x)=0$ 时，条件 (4.21) 成立，且 $x$ 最优。

**图 4.3** $\mathbf{R}$ 上的一个拟凸函数 $f$，其局部最优点 $x$ 不是全局最优点。这个例子表明，对凸函数有效的简单最优性条件 $f'(x)=0$ 对拟凸函数并不成立。

#### 通过凸可行性问题求解拟凸优化

拟凸优化的一种一般方法依赖于用一族凸不等式表示拟凸函数下水平集，如 §3.4.5 所述。令 $\phi_t:\mathbf{R}^n\to\mathbf{R},\ t\in\mathbf{R}$，是一族凸函数，满足

$$
f_0(x)\le t\quad \Longleftrightarrow\quad \phi_t(x)\le0,
$$

并且对每个 $x$，$\phi_t(x)$ 是 $t$ 的非增函数，即当 $s\ge t$ 时 $\phi_s(x)\le\phi_t(x)$。

令 $p^\star$ 表示拟凸优化问题 (4.24) 的最优值。如果可行性问题

$$
\begin{array}{ll}
\text{find} & x\\
\text{subject to} & \phi_t(x)\le0\\
& f_i(x)\le0,\quad i=1,\ldots,m\\
& Ax=b
\end{array}
\tag{4.26}
$$

可行，则有 $p^\star\le t$。反过来，若问题 (4.26) 不可行，则可得 $p^\star\ge t$。问题 (4.26) 是凸可行性问题，因为所有不等式约束函数都是凸的，等式约束是线性的。因此，可以通过求解凸可行性问题 (4.26)，检查拟凸优化问题的最优值 $p^\star$ 是小于还是大于给定值 $t$。若该凸可行性问题可行，则 $p^\star\le t$，且任何可行点 $x$ 对拟凸问题可行并满足 $f_0(x)\le t$。若该凸可行性问题不可行，则知道 $p^\star\ge t$。

这一观察可以作为一个简单算法的基础：用二分法求解拟凸优化问题 (4.24)，每一步求解一个凸可行性问题。假设问题可行，并从已知包含最优值 $p^\star$ 的区间 $[l,u]$ 开始。然后在中点 $t=(l+u)/2$ 处求解凸可行性问题，判断最优值在区间的下半部分还是上半部分，并相应更新区间。这样得到的新区间仍包含最优值，但长度为原区间的一半。重复此过程直到区间长度足够小。

**算法 4.1 拟凸优化的二分法。**

给定 $l\le p^\star$、$u\ge p^\star$、容差 $\epsilon>0$。

重复：

1. $t:=(l+u)/2$。
2. 求解凸可行性问题 (4.26)。
3. 如果 (4.26) 可行，令 $u:=t$；否则令 $l:=t$。

直到 $u-l\le\epsilon$。

区间 $[l,u]$ 保证在每一步都包含 $p^\star$，即 $l\le p^\star\le u$。每次迭代都把区间二等分，所以 $k$ 次迭代后区间长度为 $2^{-k}(u-l)$，其中 $u-l$ 是初始区间长度。因此，算法终止前恰好需要 $\lceil\log_2((u-l)/\epsilon)\rceil$ 次迭代。每一步都涉及求解凸可行性问题 (4.26)。

## 4.3 线性优化问题

当目标函数和约束函数全是仿射函数时，问题称为线性规划（linear program, LP）。一般线性规划具有形式

$$
\begin{array}{ll}
\text{minimize} & c^Tx+d\\
\text{subject to} & Gx\preceq h\\
& Ax=b,
\end{array}
\tag{4.27}
$$

其中 $G\in\mathbf{R}^{m\times n}$ 且 $A\in\mathbf{R}^{p\times n}$。线性规划当然是凸优化问题。

目标函数中的常数 $d$ 通常会省略，因为它不影响最优集（或可行集）。由于最大化仿射目标 $c^Tx+d$ 等价于最小化 $-c^Tx-d$（仍是凸的），我们也把具有仿射目标和仿射约束函数的最大化问题称为 LP。

LP 的几何解释如图 4.4 所示。LP (4.27) 的可行集是多面体 $P$；问题是在 $P$ 上最小化仿射函数 $c^Tx+d$（或等价地，最小化线性函数 $c^Tx$）。

**图 4.4** LP 的几何解释。可行集 $P$ 是一个多面体，显示为阴影区域。目标 $c^Tx$ 是线性的，所以它的水平集是与 $c$ 正交的超平面（虚线）。点 $x^\star$ 是最优的；它是 $P$ 中沿 $-c$ 方向尽可能远的点。

#### 标准形式和不等式形式的线性规划

LP (4.27) 的两个特例非常常见，因此有专门名称。标准形式 LP 中唯一的不等式是逐分量非负性约束 $x\succeq0$：

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & Ax=b\\
& x\succeq0.
\end{array}
\tag{4.28}
$$

如果 LP 没有等式约束，则称为不等式形式 LP，通常写作

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & Ax\preceq b.
\end{array}
\tag{4.29}
$$

#### 将 LP 转换为标准形式

有时将一般 LP (4.27) 转换为标准形式 (4.28) 很有用（例如为了使用标准形式 LP 的算法）。第一步是为不等式引入松弛变量 $s_i$，得到

$$
\begin{array}{ll}
\text{minimize} & c^Tx+d\\
\text{subject to} & Gx+s=h\\
& Ax=b\\
& s\succeq0.
\end{array}
$$

第二步是把变量 $x$ 表示为两个非负变量 $x^+$ 和 $x^-$ 的差，即 $x=x^+-x^-$，$x^+,x^-\succeq0$。这得到

$$
\begin{array}{ll}
\text{minimize} & c^Tx^+-c^Tx^-+d\\
\text{subject to} & Gx^+-Gx^-+s=h\\
& Ax^+-Ax^-=b\\
& x^+\succeq0,\quad x^-\succeq0,\quad s\succeq0,
\end{array}
$$

这是变量为 $x^+,x^-,s$ 的标准形式 LP。（关于该问题与原问题 (4.27) 的等价性，见习题 4.10。）

这些问题变换技术（以及后续例子和习题中的许多其他技术）可用于把许多问题表述为线性规划。稍微滥用术语时，即使一个问题不具有形式 (4.27)，只要它能表述为 LP，也常称它为 LP。

### 4.3.1 例子

LP 出现在大量领域和应用中；这里给出几个典型例子。

#### 饮食问题

健康饮食包含 $m$ 种不同营养素，含量至少分别为 $b_1,\ldots,b_m$。可以通过选择 $n$ 种不同食物的非负数量 $x_1,\ldots,x_n$ 来组成这样的饮食。单位数量的食物 $j$ 含有营养素 $i$ 的数量为 $a_{ij}$，成本为 $c_j$。我们希望确定满足营养需求的最便宜饮食。该问题可表述为 LP：

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & Ax\succeq b\\
& x\succeq0.
\end{array}
$$

该问题的若干变体也可以表述为 LP。例如，可以要求某种营养素的含量精确等于某个值（得到线性等式约束），或者除上述下界外再加上营养素含量的上界。

#### 多面体的 Chebyshev 中心

考虑寻找位于由线性不等式描述的多面体

$$
P=\{x\in\mathbf{R}^n\mid a_i^Tx\le b_i,\ i=1,\ldots,m\}
$$

内部的最大欧几里得球。（最优球的中心称为该多面体的 Chebyshev 中心；它是多面体内部最深的点，即离边界最远的点；见 §8.5.1。）把球表示为

$$
B=\{x_c+u\mid \|u\|_2\le r\}.
$$

问题变量是中心 $x_c\in\mathbf{R}^n$ 和半径 $r$；目标是在约束 $B\subseteq P$ 下最大化 $r$。

先考虑更简单的约束：$B$ 位于一个半空间 $a_i^Tx\le b_i$ 中，即

$$
\|u\|_2\le r\quad\Longrightarrow\quad a_i^T(x_c+u)\le b_i.
\tag{4.30}
$$

由于

$$
\sup\{a_i^Tu\mid \|u\|_2\le r\}=r\|a_i\|_2,
$$

可以把 (4.30) 写成

$$
a_i^Tx_c+r\|a_i\|_2\le b_i,
\tag{4.31}
$$

这是关于 $x_c$ 和 $r$ 的线性不等式。换言之，要求球位于由不等式 $a_i^Tx\le b_i$ 确定的半空间内，可以写成一个线性不等式。

因此，$B\subseteq P$ 当且仅当 (4.31) 对所有 $i=1,\ldots,m$ 成立。所以 Chebyshev 中心可通过求解 LP

$$
\begin{array}{ll}
\text{maximize} & r\\
\text{subject to} & a_i^Tx_c+r\|a_i\|_2\le b_i,\quad i=1,\ldots,m
\end{array}
$$

来确定，变量为 $r$ 和 $x_c$。（更多关于 Chebyshev 中心的内容见 §8.5.1。）

#### 动态活动规划

考虑在 $N$ 个时期内选择或规划经济中 $n$ 个活动或部门的活动水平。令 $x_j(t)\ge0,\ t=1,\ldots,N$，表示部门 $j$ 在时期 $t$ 的活动水平。活动按其活动水平比例消费并生产产品或商品。单位活动 $j$ 生产商品 $i$ 的数量为 $a_{ij}$，消费商品 $i$ 的数量为 $b_{ij}$。时期 $t$ 的商品总产出为 $Ax(t)\in\mathbf{R}^m$，消费量为 $Bx(t)\in\mathbf{R}^m$。（虽然称为“商品”，它们也可以包括污染物等不希望出现的产物。）

某时期消费的商品不能超过上一时期生产的商品：必须有 $Bx(t+1)\preceq Ax(t)$，$t=1,\ldots,N-1$。给定初始商品向量 $g_0\in\mathbf{R}^m$，它约束第一时期活动水平：$Bx(1)\preceq g_0$。活动未消费的剩余商品向量为

$$
\begin{aligned}
s(0)&=g_0-Bx(1),\\
s(t)&=Ax(t)-Bx(t+1),\quad t=1,\ldots,N-1,\\
s(N)&=Ax(N).
\end{aligned}
$$

目标是最大化剩余商品的折现总价值：

$$
c^Ts(0)+\gamma c^Ts(1)+\cdots+\gamma^N c^Ts(N),
$$

其中 $c\in\mathbf{R}^m$ 给出商品价值，$\gamma>0$ 是折现因子。（若第 $i$ 个产品是不希望出现的，例如污染物，则 $c_i$ 为负；$|c_i|$ 是单位处置成本。）

综合起来得到 LP：

$$
\begin{array}{ll}
\text{maximize} & c^Ts(0)+\gamma c^Ts(1)+\cdots+\gamma^N c^Ts(N)\\
\text{subject to} & x(t)\succeq0,\quad t=1,\ldots,N\\
& s(t)\succeq0,\quad t=0,\ldots,N\\
& s(0)=g_0-Bx(1)\\
& s(t)=Ax(t)-Bx(t+1),\quad t=1,\ldots,N-1\\
& s(N)=Ax(N),
\end{array}
$$

变量为 $x(1),\ldots,x(N),s(0),\ldots,s(N)$。这是一个标准形式 LP；变量 $s(t)$ 是约束 $Bx(t+1)\preceq Ax(t)$ 相关的松弛变量。

#### Chebyshev 不等式

考虑集合 $\{u_1,\ldots,u_n\}\subseteq\mathbf{R}$ 上离散随机变量 $x$ 的概率分布。用向量 $p\in\mathbf{R}^n$ 描述 $x$ 的分布，其中

$$
p_i=\operatorname{prob}(x=u_i),
$$

所以 $p$ 满足 $p\succeq0$ 且 $\mathbf{1}^Tp=1$。反过来，若 $p$ 满足这些条件，则它定义了 $x$ 的一个概率分布。假设 $u_i$ 已知且固定，但分布 $p$ 未知。

如果 $f$ 是 $x$ 的任意函数，则

$$
\mathbf{E}f=\sum_{i=1}^n p_i f(u_i)
$$

是 $p$ 的线性函数。如果 $S$ 是 $\mathbf{R}$ 的任意子集，则

$$
\operatorname{prob}(x\in S)=\sum_{u_i\in S}p_i
$$

也是 $p$ 的线性函数。

虽然不知道 $p$，但给定如下形式的先验知识：某些 $x$ 的函数期望以及某些 $\mathbf{R}$ 子集概率的上下界。这些先验知识可表示为关于 $p$ 的线性不等式约束

$$
\alpha_i\le a_i^Tp\le\beta_i,\quad i=1,\ldots,m.
$$

问题是给出 $\mathbf{E}f_0(x)=a_0^Tp$ 的下界和上界，其中 $f_0$ 是 $x$ 的某个函数。

为寻找下界，求解 LP

$$
\begin{array}{ll}
\text{minimize} & a_0^Tp\\
\text{subject to} & p\succeq0,\quad \mathbf{1}^Tp=1\\
& \alpha_i\le a_i^Tp\le\beta_i,\quad i=1,\ldots,m,
\end{array}
$$

变量为 $p$。该 LP 的最优值给出所有与先验信息一致的分布中 $\mathbf{E}f_0(X)$ 的最低可能值。而且该界是紧的：最优解给出一个与先验信息一致并达到该下界的分布。类似地，可在同样约束下最大化 $a_0^Tp$ 来寻找最佳上界。（§7.4.1 将更详细地讨论 Chebyshev 不等式。）

#### 分段线性最小化

考虑最小化分段线性凸函数的无约束问题

$$
f(x)=\max_{i=1,\ldots,m}(a_i^Tx+b_i).
$$

先形成上图问题，再把一个最大值不等式写成 $m$ 个单独不等式，可将该问题变换为等价 LP：

$$
\begin{array}{ll}
\text{minimize} & t\\
\text{subject to} & a_i^Tx+b_i\le t,\quad i=1,\ldots,m.
\end{array}
$$

这是变量为 $x,t$ 的不等式形式 LP。

### 4.3.2 线性分式规划

在多面体上最小化两个仿射函数之比的问题称为线性分式规划：

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & Gx\preceq h\\
& Ax=b,
\end{array}
\tag{4.32}
$$

其中目标函数为

$$
f_0(x)=\frac{c^Tx+d}{e^Tx+f},
\qquad
\operatorname{dom}f_0=\{x\mid e^Tx+f>0\}.
$$

该目标函数是拟凸的（事实上是拟线性的），所以线性分式规划是拟凸优化问题。

#### 变换为线性规划

如果可行集

$$
\{x\mid Gx\preceq h,\ Ax=b,\ e^Tx+f>0\}
$$

非空，则线性分式规划 (4.32) 可变换为等价 LP

$$
\begin{array}{ll}
\text{minimize} & c^Ty+dz\\
\text{subject to} & Gy-hz\preceq0\\
& Ay-bz=0\\
& e^Ty+fz=1\\
& z\ge0,
\end{array}
\tag{4.33}
$$

变量为 $y,z$。

为说明等价性，先注意若 $x$ 对 (4.32) 可行，则

$$
y=\frac{x}{e^Tx+f},\qquad
z=\frac{1}{e^Tx+f}
$$

对 (4.33) 可行，且具有相同目标值 $c^Ty+dz=f_0(x)$。因此 (4.32) 的最优值不小于 (4.33) 的最优值。

反过来，若 $(y,z)$ 对 (4.33) 可行且 $z\ne0$，则 $x=y/z$ 对 (4.32) 可行，且具有相同目标值 $f_0(x)=c^Ty+dz$。若 $(y,z)$ 对 (4.33) 可行且 $z=0$，并且 $x_0$ 对 (4.32) 可行，则 $x=x_0+ty$ 对所有 $t\ge0$ 都对 (4.32) 可行。此外，$\lim_{t\to\infty}f_0(x_0+ty)=c^Ty+dz$，所以可以在 (4.32) 中找到目标值任意接近 $(y,z)$ 目标值的可行点。因此 (4.32) 的最优值不大于 (4.33) 的最优值。

#### 广义线性分式规划

线性分式规划 (4.32) 的一个推广是广义线性分式规划，其中

$$
f_0(x)=
\max_{i=1,\ldots,r}
\frac{c_i^Tx+d_i}{e_i^Tx+f_i},
\qquad
\operatorname{dom}f_0=
\{x\mid e_i^Tx+f_i>0,\ i=1,\ldots,r\}.
$$

目标函数是 $r$ 个拟凸函数的逐点最大值，因此是拟凸的，所以这个问题是拟凸的。当 $r=1$ 时，它退化为标准线性分式规划。

**例 4.7 Von Neumann 增长问题。** 考虑一个有 $n$ 个部门的经济，当前时期活动水平为 $x_i>0$，下一时期活动水平为 $x_i^+>0$。（这里仅考虑一个时期。）有 $m$ 种商品由这些活动消费并生产：活动水平 $x$ 消费商品 $Bx\in\mathbf{R}^m$，生产商品 $Ax$。下一时期消费的商品不能超过当前时期生产的商品，即 $Bx^+\preceq Ax$。部门 $i$ 在该时期的增长率为 $x_i^+/x_i$。

Von Neumann 增长问题是寻找活动水平向量 $x$，使经济所有部门中的最小增长率最大。该问题可表示为广义线性分式问题

$$
\begin{array}{ll}
\text{maximize} & \min_{i=1,\ldots,n} x_i^+/x_i\\
\text{subject to} & x^+\succeq0\\
& Bx^+\preceq Ax
\end{array}
$$

定义域为 $\{(x,x^+)\mid x\succ0\}$。注意该问题关于 $x$ 和 $x^+$ 是齐次的，所以可以用显式约束 $x\succeq\mathbf{1}$ 替换隐式约束 $x\succ0$。

## 4.4 二次优化问题

如果目标函数是（凸）二次函数，约束函数是仿射函数，则凸优化问题 (4.15) 称为二次规划（quadratic program, QP）。二次规划可表示为

$$
\begin{array}{ll}
\text{minimize} & (1/2)x^TPx+q^Tx+r\\
\text{subject to} & Gx\preceq h\\
& Ax=b,
\end{array}
\tag{4.34}
$$

其中 $P\in\mathbf{S}_+^n$，$G\in\mathbf{R}^{m\times n}$，$A\in\mathbf{R}^{p\times n}$。在二次规划中，我们在多面体上最小化凸二次函数，如图 4.5 所示。

**图 4.5** QP 的几何示意。可行集 $P$ 是一个多面体，显示为阴影区域。目标函数为凸二次函数，其等值线显示为虚线曲线。点 $x^\star$ 是最优的。

如果 (4.15) 中的目标和不等式约束函数都是（凸）二次函数，如

$$
\begin{array}{ll}
\text{minimize} & (1/2)x^TP_0x+q_0^Tx+r_0\\
\text{subject to} & (1/2)x^TP_ix+q_i^Tx+r_i\le0,\quad i=1,\ldots,m\\
& Ax=b,
\end{array}
\tag{4.35}
$$

其中 $P_i\in\mathbf{S}_+^n,\ i=0,1,\ldots,m$，则称该问题为二次约束二次规划（QCQP）。在 QCQP 中，我们在由若干椭球（当 $P_i\succ0$ 时）交成的可行区域上最小化凸二次函数。

二次规划包含线性规划作为特例，只需在 (4.34) 中取 $P=0$。二次约束二次规划包含二次规划（因而也包含线性规划）作为特例，只需在 (4.35) 中对 $i=1,\ldots,m$ 取 $P_i=0$。

### 4.4.1 例子

#### 最小二乘和回归

最小化凸二次函数

$$
\|Ax-b\|_2^2=x^TA^TAx-2b^TAx+b^Tb
$$

是一个（无约束）QP。它出现在许多领域，并有许多名称，例如回归分析或最小二乘逼近。这个问题足够简单，具有著名解析解 $x=A^\dagger b$，其中 $A^\dagger$ 是 $A$ 的伪逆（见 §A.5.4）。

加入线性不等式约束后，该问题称为带约束回归或带约束最小二乘，此时不再有简单解析解。例如，可以考虑带变量上下界的回归：

$$
\begin{array}{ll}
\text{minimize} & \|Ax-b\|_2^2\\
\text{subject to} & l_i\le x_i\le u_i,\quad i=1,\ldots,n,
\end{array}
$$

这是一个 QP。（第 6 章和第 7 章将更深入研究最小二乘和回归问题。）

#### 多面体之间的距离

$\mathbf{R}^n$ 中多面体 $P_1=\{x\mid A_1x\preceq b_1\}$ 和 $P_2=\{x\mid A_2x\preceq b_2\}$ 之间的欧几里得距离定义为

$$
\operatorname{dist}(P_1,P_2)=
\inf\{\|x_1-x_2\|_2\mid x_1\in P_1,\ x_2\in P_2\}.
$$

如果两个多面体相交，则距离为零。

为寻找 $P_1$ 和 $P_2$ 之间的距离，可以求解 QP

$$
\begin{array}{ll}
\text{minimize} & \|x_1-x_2\|_2^2\\
\text{subject to} & A_1x_1\preceq b_1,\quad A_2x_2\preceq b_2,
\end{array}
$$

变量为 $x_1,x_2\in\mathbf{R}^n$。该问题不可行当且仅当其中一个多面体为空。最优值为零当且仅当两个多面体相交，此时最优的 $x_1$ 和 $x_2$ 相等（并且是交集 $P_1\cap P_2$ 中的点）。否则，最优的 $x_1$ 和 $x_2$ 分别是 $P_1$ 和 $P_2$ 中彼此最近的点。（第 8 章将更详细研究涉及距离的几何问题。）

#### 方差界

再次考虑 Chebyshev 不等式例子（第 150 页），其中变量是未知概率分布 $p\in\mathbf{R}^n$，并且有一些先验信息。随机变量 $f(x)$ 的方差为

$$
\mathbf{E}f^2-(\mathbf{E}f)^2
=\sum_{i=1}^n f_i^2p_i-\left(\sum_{i=1}^n f_ip_i\right)^2,
$$

其中 $f_i=f(u_i)$。这是 $p$ 的凹二次函数。

因此，在给定先验信息下，可以通过求解 QP 来最大化 $f(x)$ 的方差：

$$
\begin{array}{ll}
\text{maximize} & \sum_{i=1}^n f_i^2p_i-\left(\sum_{i=1}^n f_ip_i\right)^2\\
\text{subject to} & p\succeq0,\quad \mathbf{1}^Tp=1\\
& \alpha_i\le a_i^Tp\le\beta_i,\quad i=1,\ldots,m.
\end{array}
$$

最优值给出所有与先验信息一致的分布中 $f(x)$ 的最大可能方差；最优 $p$ 给出达到该最大方差的分布。

#### 随机成本的线性规划

考虑 LP

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & Gx\preceq h\\
& Ax=b,
\end{array}
$$

变量为 $x\in\mathbf{R}^n$。假设成本向量 $c\in\mathbf{R}^n$ 是随机的，均值为 $\bar c$，协方差为 $\mathbf{E}(c-\bar c)(c-\bar c)^T=\Sigma$。（为简单起见，假设其他问题参数是确定性的。）对给定 $x\in\mathbf{R}^n$，成本 $c^Tx$ 是一个标量随机变量，均值为 $\mathbf{E}c^Tx=\bar c^Tx$，方差为

$$
\operatorname{var}(c^Tx)=\mathbf{E}(c^Tx-\mathbf{E}c^Tx)^2=x^T\Sigma x.
$$

一般而言，较小期望成本和较小成本方差之间存在折中。考虑方差的一种方法是最小化成本期望和方差的线性组合：

$$
\mathbf{E}c^Tx+\gamma\operatorname{var}(c^Tx),
$$

称为风险敏感成本。参数 $\gamma\ge0$ 称为风险厌恶参数，因为它设定了成本方差和期望值的相对重要性。（当 $\gamma>0$ 时，我们愿意用期望成本的增加换取足够大的成本方差下降。）

为最小化风险敏感成本，求解 QP

$$
\begin{array}{ll}
\text{minimize} & \bar c^Tx+\gamma x^T\Sigma x\\
\text{subject to} & Gx\preceq h\\
& Ax=b.
\end{array}
$$

#### Markowitz 投资组合优化

考虑一个经典投资组合问题：在一个时期内持有 $n$ 个资产或股票。令 $x_i$ 表示整个时期中持有的资产 $i$ 数量，以时期开始时的价格按美元计。资产 $i$ 的普通多头头寸对应 $x_i>0$；空头头寸（即期末买入该资产的义务）对应 $x_i<0$。令 $p_i$ 表示资产 $i$ 在该时期内的相对价格变化，即价格变化除以期初价格。投资组合的总收益为 $r=p^Tx$（以美元计）。优化变量是投资组合向量 $x\in\mathbf{R}^n$。

可以考虑大量关于投资组合的约束。最简单的约束是 $x_i\ge0$（即不允许做空）和 $\mathbf{1}^Tx=B$（即要投资的总预算为 $B$，通常取为 1）。

对价格变化采用随机模型：$p\in\mathbf{R}^n$ 是随机向量，均值 $\bar p$ 和协方差 $\Sigma$ 已知。因此对投资组合 $x\in\mathbf{R}^n$，收益 $r$ 是标量随机变量，均值为 $\bar p^Tx$，方差为 $x^T\Sigma x$。投资组合 $x$ 的选择涉及收益均值和方差之间的折中。

Markowitz 引入的经典投资组合优化问题是 QP

$$
\begin{array}{ll}
\text{minimize} & x^T\Sigma x\\
\text{subject to} & \bar p^Tx\ge r_\mathrm{min}\\
& \mathbf{1}^Tx=1,\quad x\succeq0,
\end{array}
$$

变量为投资组合 $x$。这里寻找一个组合，在满足最低可接受均值收益 $r_\mathrm{min}$、预算约束和不做空约束的同时，使收益方差（与投资组合风险相关）最小。

许多扩展都是可能的。例如，一个标准扩展是允许空头头寸，即 $x_i<0$。为此引入变量 $x^\mathrm{long}$ 和 $x^\mathrm{short}$，满足

$$
x^\mathrm{long}\succeq0,\quad
x^\mathrm{short}\succeq0,\quad
x=x^\mathrm{long}-x^\mathrm{short},\quad
\mathbf{1}^Tx^\mathrm{short}\le\eta\,\mathbf{1}^Tx^\mathrm{long}.
$$

最后一个约束把期初总空头头寸限制为期初总多头头寸的某个比例 $\eta$。

另一个扩展是在投资组合优化问题中包含线性交易成本。从给定初始组合 $x^\mathrm{init}$ 出发，买入和卖出资产以达到组合 $x$，然后在整个时期中持有该组合。买卖资产会产生交易费用，费用与买入或卖出金额成正比。为处理这一点，引入变量 $u^\mathrm{buy}$ 和 $u^\mathrm{sell}$，确定持有期前买入和卖出的每种资产数量。约束为

$$
x=x^\mathrm{init}+u^\mathrm{buy}-u^\mathrm{sell},\qquad
u^\mathrm{buy}\succeq0,\quad u^\mathrm{sell}\succeq0.
$$

用“初始买卖及交易费用造成净现金为零”的条件替换简单预算约束 $\mathbf{1}^Tx=1$：

$$
(1-f^\mathrm{sell})\mathbf{1}^Tu^\mathrm{sell}
=(1+f^\mathrm{buy})\mathbf{1}^Tu^\mathrm{buy}.
$$

左侧是卖出资产所得总额减去卖出交易费，右侧是买入资产的总成本（包括交易费）。常数 $f^\mathrm{buy}\ge0$ 和 $f^\mathrm{sell}\ge0$ 是买入和卖出的交易费率（为简单起见，假设所有资产相同）。在最小化收益方差、满足最低均值收益、预算和交易约束的问题中，变量为 $x,u^\mathrm{buy},u^\mathrm{sell}$，该问题是 QP。

### 4.4.2 二阶锥规划

与二次规划密切相关的一类问题是二阶锥规划（second-order cone program, SOCP）：

$$
\begin{array}{ll}
\text{minimize} & f^Tx\\
\text{subject to} & \|A_ix+b_i\|_2\le c_i^Tx+d_i,\quad i=1,\ldots,m\\
& Fx=g,
\end{array}
\tag{4.36}
$$

其中 $x\in\mathbf{R}^n$ 是优化变量，$A_i\in\mathbf{R}^{n_i\times n}$，$F\in\mathbf{R}^{p\times n}$。形如

$$
\|Ax+b\|_2\le c^Tx+d
$$

的约束称为二阶锥约束，因为它等价于要求仿射函数 $(Ax+b,c^Tx+d)$ 落在 $\mathbf{R}^{k+1}$ 中的二阶锥内。

当 $c_i=0,\ i=1,\ldots,m$ 时，SOCP (4.36) 等价于一个 QCQP（通过对每个约束平方得到）。类似地，如果 $A_i=0,\ i=1,\ldots,m$，则二阶锥约束退化为线性不等式，SOCP 退化为 LP。因此，SOCP 包含 LP 和许多 QCQP 作为特例。

#### 鲁棒线性规划

考虑不等式形式的线性规划

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & a_i^Tx\le b_i,\quad i=1,\ldots,m,
\end{array}
$$

其中参数 $c,a_i,b_i$ 有一些不确定性或变化。为简化说明，假设 $c$ 和 $b_i$ 固定，而 $a_i$ 已知位于给定椭球中：

$$
a_i\in\mathcal{E}_i=\{\bar a_i+P_iu\mid \|u\|_2\le1\},
$$

其中 $P_i\in\mathbf{R}^{n\times n}$。（若 $P_i$ 奇异，则得到维数为 $\operatorname{rank}P_i$ 的“扁”椭球；$P_i=0$ 表示 $a_i$ 完全已知。）

要求约束对参数 $a_i$ 的所有可能值都成立，得到鲁棒线性规划

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & a_i^Tx\le b_i\quad \text{对所有 }a_i\in\mathcal{E}_i,\quad i=1,\ldots,m.
\end{array}
\tag{4.37}
$$

鲁棒线性约束“对所有 $a_i\in\mathcal{E}_i$ 有 $a_i^Tx\le b_i$”可表示为

$$
\sup\{a_i^Tx\mid a_i\in\mathcal{E}_i\}\le b_i,
$$

其左侧为

$$
\sup\{a_i^Tx\mid a_i\in\mathcal{E}_i\}
=\bar a_i^Tx+\sup\{u^TP_i^Tx\mid \|u\|_2\le1\}
=\bar a_i^Tx+\|P_i^Tx\|_2.
$$

因此鲁棒线性约束可写为

$$
\bar a_i^Tx+\|P_i^Tx\|_2\le b_i,
$$

这显然是二阶锥约束。因此鲁棒 LP (4.37) 可表示为 SOCP

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & \bar a_i^Tx+\|P_i^Tx\|_2\le b_i,\quad i=1,\ldots,m.
\end{array}
$$

注意，额外的范数项起正则化作用；它们防止 $x$ 在参数 $a_i$ 有显著不确定性的方向上过大。

#### 随机约束的线性规划

上述鲁棒 LP 也可以放在统计框架下理解。这里假设参数 $a_i$ 是独立高斯随机向量，均值为 $\bar a_i$，协方差为 $\Sigma_i$。要求每个约束 $a_i^Tx\le b_i$ 以超过 $\eta$ 的概率（或置信度）成立，其中 $\eta\ge0.5$，即

$$
\operatorname{prob}(a_i^Tx\le b_i)\ge\eta.
\tag{4.38}
$$

可以证明该概率约束可表示为二阶锥约束。令 $u=a_i^Tx$，其方差记为 $\sigma^2$，则约束可写为

$$
\operatorname{prob}\left(\frac{u-\bar u}{\sigma}\le\frac{b_i-\bar u}{\sigma}\right)\ge\eta.
$$

由于 $(u-\bar u)/\sigma$ 是零均值、单位方差的高斯变量，上述概率为 $\Phi((b_i-\bar u)/\sigma)$，其中

$$
\Phi(z)=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^z e^{-t^2/2}\,dt
$$

是标准高斯随机变量的累积分布函数。因此概率约束 (4.38) 可表示为

$$
\frac{b_i-\bar u}{\sigma}\ge\Phi^{-1}(\eta),
$$

或等价地，

$$
\bar u+\Phi^{-1}(\eta)\sigma\le b_i.
$$

由 $\bar u=\bar a_i^Tx$ 和 $\sigma=(x^T\Sigma_i x)^{1/2}$ 得到

$$
\bar a_i^Tx+\Phi^{-1}(\eta)\|\Sigma_i^{1/2}x\|_2\le b_i.
$$

由于假设 $\eta\ge1/2$，有 $\Phi^{-1}(\eta)\ge0$，所以这是二阶锥约束。

综上，问题

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & \operatorname{prob}(a_i^Tx\le b_i)\ge\eta,\quad i=1,\ldots,m
\end{array}
$$

可表示为 SOCP

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & \bar a_i^Tx+\Phi^{-1}(\eta)\|\Sigma_i^{1/2}x\|_2\le b_i,\quad i=1,\ldots,m.
\end{array}
$$

（第 6 章将更深入讨论鲁棒凸优化问题。另见习题 4.13、4.28 和 4.59。）

**例 4.8 带损失风险约束的投资组合优化。** 再考虑前述经典 Markowitz 投资组合问题（第 155 页）。这里假设价格变化向量 $p\in\mathbf{R}^n$ 是高斯随机变量，均值为 $\bar p$，协方差为 $\Sigma$。因此收益 $r$ 是高斯随机变量，均值为 $\bar r=\bar p^Tx$，方差为 $\sigma_r^2=x^T\Sigma x$。

考虑损失风险约束

$$
\operatorname{prob}(r\le\alpha)\le\beta,
\tag{4.39}
$$

其中 $\alpha$ 是给定的不希望出现的收益水平（例如大额亏损），$\beta$ 是给定的最大概率。

和上述鲁棒 LP 的随机解释一样，可以用单位高斯随机变量的累积分布函数 $\Phi$ 表示该约束。不等式 (4.39) 等价于

$$
\bar p^Tx+\Phi^{-1}(\beta)\|\Sigma^{1/2}x\|_2\ge\alpha.
$$

当 $\beta\le1/2$（即 $\Phi^{-1}(\beta)\le0$）时，该损失风险约束是二阶锥约束。（若 $\beta>1/2$，损失风险约束关于 $x$ 变为非凸。）

因此，在 $\beta\le1/2$ 时，最大化期望收益并约束损失风险的问题可以表述为带一个二阶锥约束的 SOCP：

$$
\begin{array}{ll}
\text{maximize} & \bar p^Tx\\
\text{subject to} & \bar p^Tx+\Phi^{-1}(\beta)\|\Sigma^{1/2}x\|_2\ge\alpha\\
& x\succeq0,\quad \mathbf{1}^Tx=1.
\end{array}
$$

该问题有许多扩展。例如，可以施加多个损失风险约束

$$
\operatorname{prob}(r\le\alpha_i)\le\beta_i,\quad i=1,\ldots,k,
$$

其中 $\beta_i\le1/2$，表达我们愿意为不同损失水平 $\alpha_i$ 接受的风险 $\beta_i$。

#### 极小曲面

考虑可微函数 $f:\mathbf{R}^2\to\mathbf{R}$，其定义域为 $C$。其图像的曲面面积为

$$
A=\int_C \sqrt{1+\|\nabla f(x)\|_2^2}\,dx
=\int_C \|(\nabla f(x),1)\|_2\,dx,
$$

这是 $f$ 的凸泛函。极小曲面问题是寻找使 $A$ 最小的函数 $f$，并满足一些约束，例如边界上 $f$ 的某些给定值。

通过离散化函数 $f$ 来近似该问题。令 $C=[0,1]\times[0,1]$，令 $f_{ij}$ 表示 $f$ 在点 $(i/K,j/K)$ 处的值，其中 $i,j=0,\ldots,K$。使用前向差分可得到 $f$ 在点 $x=(i/K,j/K)$ 处梯度的近似：

$$
\nabla f(x)\approx
K\begin{bmatrix}
f_{i+1,j}-f_{ij}\\
f_{i,j+1}-f_{ij}
\end{bmatrix}.
$$

将其代入图像面积表达式，并把积分近似为求和，得到图像面积的近似：

$$
A\approx A_\mathrm{disc}
=\frac{1}{K^2}\sum_{i,j=0}^{K-1}
\left\|
\begin{bmatrix}
K(f_{i+1,j}-f_{ij})\\
K(f_{i,j+1}-f_{ij})\\
1
\end{bmatrix}
\right\|_2.
$$

离散面积近似 $A_\mathrm{disc}$ 是 $f_{ij}$ 的凸函数。

可以对 $f_{ij}$ 考虑多种约束，例如对任意条目（例如边界值）的等式或不等式约束，或对其矩的约束。作为例子，考虑寻找正方形左右边界具有固定边界值的极小面积曲面：

$$
\begin{array}{ll}
\text{minimize} & A_\mathrm{disc}\\
\text{subject to} & f_{0j}=l_j,\quad j=0,\ldots,K\\
& f_{Kj}=r_j,\quad j=0,\ldots,K,
\end{array}
\tag{4.40}
$$

其中 $f_{ij},\ i,j=0,\ldots,K$ 是变量，$l_j,r_j$ 是左右两侧给定边界值。

通过引入新变量 $t_{ij},\ i,j=0,\ldots,K-1$，可把问题 (4.40) 转换为 SOCP：

$$
\begin{array}{ll}
\text{minimize} & \frac{1}{K^2}\sum_{i,j=0}^{K-1} t_{ij}\\
\text{subject to} &
\left\|
\begin{bmatrix}
K(f_{i+1,j}-f_{ij})\\
K(f_{i,j+1}-f_{ij})\\
1
\end{bmatrix}
\right\|_2\le t_{ij},\quad i,j=0,\ldots,K-1\\
& f_{0j}=l_j,\quad j=0,\ldots,K\\
& f_{Kj}=r_j,\quad j=0,\ldots,K.
\end{array}
$$

## 4.5 几何规划

本节描述一族在自然形式下并非凸的问题。不过，通过变量变换以及目标函数和约束函数变换，这些问题可以转化为凸优化问题。

### 4.5.1 单项式和正项式

定义域为 $\mathbf{R}_{++}^n$ 的函数 $f:\mathbf{R}^n\to\mathbf{R}$

$$
f(x)=cx_1^{a_1}x_2^{a_2}\cdots x_n^{a_n},
\tag{4.41}
$$

其中 $c>0$ 且 $a_i\in\mathbf{R}$，称为单项式函数，或简称单项式。单项式的指数 $a_i$ 可以是任意实数，包括分数或负数，但系数 $c$ 只能为正。（这里“单项式”与代数中的标准定义冲突，后者要求指数为非负整数，但这不应造成混淆。）单项式之和，即形如

$$
f(x)=\sum_{k=1}^K c_k x_1^{a_{1k}}x_2^{a_{2k}}\cdots x_n^{a_{nk}},
\tag{4.42}
$$

其中 $c_k>0$，称为正项式函数（有 $K$ 项），或简称正项式。

正项式对加法、乘法和非负缩放封闭。单项式对乘法和除法封闭。若正项式乘以单项式，结果仍是正项式；类似地，正项式除以单项式，结果仍是正项式。

### 4.5.2 几何规划

形如

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & f_i(x)\le1,\quad i=1,\ldots,m\\
& h_i(x)=1,\quad i=1,\ldots,p
\end{array}
\tag{4.43}
$$

的优化问题称为几何规划（geometric program, GP），其中 $f_0,\ldots,f_m$ 是正项式，$h_1,\ldots,h_p$ 是单项式。该问题的定义域为 $D=\mathbf{R}_{++}^n$；约束 $x\succ0$ 是隐式的。

#### 几何规划的扩展

若 $f$ 是正项式且 $h$ 是单项式，则约束 $f(x)\le h(x)$ 可通过写成 $f(x)/h(x)\le1$ 来处理（因为 $f/h$ 是正项式）。这包括形如 $f(x)\le a$ 的约束作为特例，其中 $f$ 是正项式且 $a>0$。类似地，若 $h_1$ 和 $h_2$ 都是非零单项式函数，则等式约束 $h_1(x)=h_2(x)$ 可写为 $h_1(x)/h_2(x)=1$（因为 $h_1/h_2$ 是单项式）。可以通过最小化非零单项式目标函数的倒数（仍是单项式）来最大化该目标。

例如，考虑问题

$$
\begin{array}{ll}
\text{maximize} & x/y\\
\text{subject to} & 2\le x\le3\\
& x^2+3y/z\le\sqrt{y}\\
& x/y=z^2,
\end{array}
$$

变量为 $x,y,z\in\mathbf{R}$（并有隐式约束 $x,y,z>0$）。使用上述简单变换，得到等价标准形式 GP：

$$
\begin{array}{ll}
\text{minimize} & x^{-1}y\\
\text{subject to} & 2x^{-1}\le1,\quad (1/3)x\le1\\
& x^2y^{-1/2}+3y^{1/2}z^{-1}\le1\\
& xy^{-1}z^{-2}=1.
\end{array}
$$

像这样很容易转换为标准形式 (4.43) 的等价 GP 的问题，我们也称为 GP。（这类似于把容易转换为 LP 的问题称为 LP。）

### 4.5.3 凸形式的几何规划

几何规划一般不是凸优化问题，但可以通过变量变换以及目标函数和约束函数变换转化为凸问题。

使用变量 $y_i=\log x_i$，于是 $x_i=e^{y_i}$。若 $f$ 是 (4.41) 中的单项式

$$
f(x)=cx_1^{a_1}x_2^{a_2}\cdots x_n^{a_n},
$$

则

$$
f(x)=f(e^{y_1},\ldots,e^{y_n})
=c(e^{y_1})^{a_1}\cdots(e^{y_n})^{a_n}
=e^{a^Ty+b},
$$

其中 $b=\log c$。变量变换 $y_i=\log x_i$ 把单项式函数变成仿射函数的指数。

类似地，若 $f$ 是 (4.42) 中的正项式

$$
f(x)=\sum_{k=1}^K c_k x_1^{a_{1k}}x_2^{a_{2k}}\cdots x_n^{a_{nk}},
$$

则

$$
f(x)=\sum_{k=1}^K e^{a_k^Ty+b_k},
$$

其中 $a_k=(a_{1k},\ldots,a_{nk})$，$b_k=\log c_k$。变量变换后，正项式变成若干仿射函数指数的和。

几何规划 (4.43) 可用新变量 $y$ 表示为

$$
\begin{array}{ll}
\text{minimize} & \sum_{k=1}^{K_0} e^{a_{0k}^Ty+b_{0k}}\\
\text{subject to} & \sum_{k=1}^{K_i} e^{a_{ik}^Ty+b_{ik}}\le1,\quad i=1,\ldots,m\\
& e^{g_i^Ty+h_i}=1,\quad i=1,\ldots,p,
\end{array}
$$

其中 $a_{ik}\in\mathbf{R}^n,\ i=0,\ldots,m$，包含原几何规划中正项式不等式约束的指数，$g_i\in\mathbf{R}^n,\ i=1,\ldots,p$，包含单项式等式约束的指数。

现在通过取对数来变换目标和约束函数，得到

$$
\begin{array}{ll}
\text{minimize} & \tilde f_0(y)=\log\left(\sum_{k=1}^{K_0}e^{a_{0k}^Ty+b_{0k}}\right)\\
\text{subject to} &
\tilde f_i(y)=\log\left(\sum_{k=1}^{K_i}e^{a_{ik}^Ty+b_{ik}}\right)\le0,\quad i=1,\ldots,m\\
& \tilde h_i(y)=g_i^Ty+h_i=0,\quad i=1,\ldots,p.
\end{array}
\tag{4.44}
$$

由于函数 $\tilde f_i$ 是凸的，$\tilde h_i$ 是仿射的，该问题是凸优化问题。我们称它为凸形式的几何规划。为区别原问题，称 (4.43) 为正项式形式的几何规划。

注意，从正项式形式 GP (4.43) 到凸形式 GP (4.44) 的变换不涉及任何计算；两个问题的问题数据相同。它只是改变目标函数和约束函数的形式。

如果所有正项式目标和约束函数都只有一项，即都是单项式，则凸形式 GP (4.44) 退化为一个（一般）线性规划。因此，可以把几何规划看作线性规划的推广或扩展。

### 4.5.4 例子

#### Frobenius 范数对角缩放

考虑矩阵 $M\in\mathbf{R}^{n\times n}$ 以及相关线性函数 $u\mapsto y=Mu$。假设缩放坐标，即变换变量为 $\tilde u=Du$、$\tilde y=Dy$，其中 $D$ 是对角矩阵，$D_{ii}>0$。在新坐标中，线性函数为 $\tilde y=DMD^{-1}\tilde u$。

现在希望选择缩放，使所得矩阵 $DMD^{-1}$ 较小。使用 Frobenius 范数平方来度量矩阵大小：

$$
\begin{aligned}
\|DMD^{-1}\|_F^2
&=\operatorname{tr}\left((DMD^{-1})^T(DMD^{-1})\right)\\
&=\sum_{i,j=1}^n (DMD^{-1})_{ij}^2\\
&=\sum_{i,j=1}^n M_{ij}^2 d_i^2/d_j^2,
\end{aligned}
$$

其中 $D=\operatorname{diag}(d)$。由于这是 $d$ 的正项式，选择缩放 $d$ 以最小化 Frobenius 范数的问题是无约束 GP：

$$
\text{minimize}\quad
\sum_{i,j=1}^n M_{ij}^2d_i^2/d_j^2,
$$

变量为 $d$。该几何规划中的指数只有 $0,2,-2$。

#### 悬臂梁设计

考虑一根悬臂梁的设计，它由 $N$ 段组成，从右到左编号为 $1,\ldots,N$，如图 4.6 所示。每段长度为 1，且具有均匀矩形截面，宽度为 $w_i$，高度为 $h_i$。在梁的右端施加垂直载荷（力）$F$。该载荷使梁向下挠曲，并在梁的每段中产生应力。假设挠度很小，材料为线弹性，Young 模量为 $E$。

**图 4.6** 由 4 段组成的分段悬臂梁。每段长度为 1，且具有矩形截面。在梁的右端施加垂直力 $F$。

问题的设计变量是 $N$ 段的宽度 $w_i$ 和高度 $h_i$。目标是最小化梁的总体积（与重量成正比）

$$
w_1h_1+\cdots+w_Nh_N,
$$

并满足若干设计约束。对各段宽度和高度施加上下界：

$$
w_\mathrm{min}\le w_i\le w_\mathrm{max},\qquad
h_\mathrm{min}\le h_i\le h_\mathrm{max},\quad i=1,\ldots,N,
$$

以及长宽比约束

$$
S_\mathrm{min}\le h_i/w_i\le S_\mathrm{max}.
$$

此外，还限制材料中的最大允许应力以及梁端的垂直挠度。

先考虑最大应力约束。第 $i$ 段中的最大应力记为 $\sigma_i$，由

$$
\sigma_i=\frac{6iF}{w_ih_i^2}
$$

给出。施加约束

$$
\frac{6iF}{w_ih_i^2}\le\sigma_\mathrm{max},\quad i=1,\ldots,N,
$$

以确保梁中任何位置的应力不超过最大允许值 $\sigma_\mathrm{max}$。

最后一个约束是梁端垂直挠度的限制，记为 $y_1$：

$$
y_1\le y_\mathrm{max}.
$$

挠度 $y_1$ 可由递推得到，该递推涉及梁各段的挠度和斜率：

$$
v_i=\frac{12(i-1/2)F}{Ew_ih_i^3}+v_{i+1},\qquad
y_i=\frac{6(i-1/3)F}{Ew_ih_i^3}+v_{i+1}+y_{i+1},
\tag{4.45}
$$

其中 $i=N,N-1,\ldots,1$，初值为 $v_{N+1}=y_{N+1}=0$。在这个递推中，$y_i$ 是第 $i$ 段右端的挠度，$v_i$ 是该点的斜率。

可以用递推 (4.45) 证明这些挠度和斜率量实际上是变量 $w,h$ 的正项式。首先 $v_{N+1}$ 和 $y_{N+1}$ 为零，因此是正项式。假设 $v_{i+1}$ 和 $y_{i+1}$ 是 $w,h$ 的正项式。式 (4.45) 左式表明 $v_i$ 是一个单项式和一个正项式（即 $v_{i+1}$）之和，因此是正项式。右式表明挠度 $y_i$ 是一个单项式和两个正项式（$v_{i+1}$ 与 $y_{i+1}$）之和，因此也是正项式。特别地，梁端挠度 $y_1$ 是正项式。

问题为

$$
\begin{array}{ll}
\text{minimize} & \sum_{i=1}^N w_ih_i\\
\text{subject to} & w_\mathrm{min}\le w_i\le w_\mathrm{max},\quad i=1,\ldots,N\\
& h_\mathrm{min}\le h_i\le h_\mathrm{max},\quad i=1,\ldots,N\\
& S_\mathrm{min}\le h_i/w_i\le S_\mathrm{max},\quad i=1,\ldots,N\\
& 6iF/(w_ih_i^2)\le\sigma_\mathrm{max},\quad i=1,\ldots,N\\
& y_1\le y_\mathrm{max},
\end{array}
\tag{4.46}
$$

变量为 $w,h$。这是 GP，因为目标是正项式，约束都可以表示为正项式不等式。（事实上，除挠度限制这个复杂正项式不等式外，所有约束都可以表示为单项式不等式。）

当段数 $N$ 很大时，正项式 $y_1$ 中出现的单项式项数近似按 $N^2$ 增长。习题 4.31 探索了该问题的另一种形式：引入 $v_1,\ldots,v_N$ 和 $y_1,\ldots,y_N$ 作为变量，并把递推的修改版本作为一组约束。该形式避免了单项式项数的增长。

#### 通过 Perron-Frobenius 理论最小化谱半径

设矩阵 $A\in\mathbf{R}^{n\times n}$ 逐元素非负，即 $A_{ij}\ge0$，并且不可约，即矩阵 $(I+A)^{n-1}$ 逐元素为正。Perron-Frobenius 定理说明，$A$ 有一个正实特征值 $\lambda_\mathrm{pf}$，等于其谱半径，即特征值模的最大值。Perron-Frobenius 特征值 $\lambda_\mathrm{pf}$ 决定 $A^k$ 当 $k\to\infty$ 时的渐近增长或衰减速率；事实上，矩阵 $((1/\lambda_\mathrm{pf})A)^k$ 收敛。粗略地说，这意味着当 $k\to\infty$ 时，如果 $\lambda_\mathrm{pf}>1$，则 $A^k$ 像 $\lambda_\mathrm{pf}^k$ 一样增长；如果 $\lambda_\mathrm{pf}<1$，则像 $\lambda_\mathrm{pf}^k$ 一样衰减。

非负矩阵理论中的一个基本结果是 Perron-Frobenius 特征值由

$$
\lambda_\mathrm{pf}=\inf\{\lambda\mid Av\preceq \lambda v\ \text{对某个 }v\succ0\}
$$

给出（并且下确界会达到）。不等式 $Av\preceq\lambda v$ 可表示为

$$
\sum_{j=1}^n \frac{A_{ij}v_j}{\lambda v_i}\le1,\quad i=1,\ldots,n,
\tag{4.47}
$$

这是关于变量 $A_{ij},v_i,\lambda$ 的一组正项式不等式。因此，条件 $\lambda_\mathrm{pf}\le\lambda$ 可用关于 $A,v,\lambda$ 的正项式不等式表示。这使我们可以用几何规划求解一些涉及 Perron-Frobenius 特征值的优化问题。

假设矩阵 $A$ 的元素是某些底层变量 $x\in\mathbf{R}^k$ 的正项式函数。在这种情况下，不等式 (4.47) 是变量 $x\in\mathbf{R}^k$、$v\in\mathbf{R}^n$ 和 $\lambda\in\mathbf{R}$ 的正项式不等式。考虑选择 $x$ 以最小化 $A$ 的 Perron-Frobenius 特征值（或谱半径），并可能满足关于 $x$ 的正项式不等式：

$$
\begin{array}{ll}
\text{minimize} & \lambda_\mathrm{pf}(A(x))\\
\text{subject to} & f_i(x)\le1,\quad i=1,\ldots,p,
\end{array}
$$

其中 $f_i$ 是正项式。利用上述刻画，可把该问题表示为 GP

$$
\begin{array}{ll}
\text{minimize} & \lambda\\
\text{subject to} &
\sum_{j=1}^n A_{ij}v_j/(\lambda v_i)\le1,\quad i=1,\ldots,n\\
& f_i(x)\le1,\quad i=1,\ldots,p,
\end{array}
$$

变量为 $x,v,\lambda$。

作为具体例子，考虑一个细菌种群动力学简单模型，时间或时期记为 $t=0,1,2,\ldots$（单位为小时）。向量 $p(t)\in\mathbf{R}_+^4$ 表征时期 $t$ 的种群年龄分布：$p_1(t)$ 是 0 到 1 小时龄细菌总数，$p_2(t)$ 是 1 到 2 小时龄细菌总数，依此类推。任意假设细菌寿命不超过 4 小时。种群按 $p(t+1)=Ap(t)$ 演化，其中

$$
A=
\begin{bmatrix}
b_1 & b_2 & b_3 & b_4\\
s_1 & 0 & 0 & 0\\
0 & s_2 & 0 & 0\\
0 & 0 & s_3 & 0
\end{bmatrix}.
$$

这里 $b_i$ 是年龄组 $i$ 中细菌的出生率，$s_i$ 是从年龄组 $i$ 存活到年龄组 $i+1$ 的存活率。假设 $b_i>0$ 且 $0<s_i<1$，这意味着矩阵 $A$ 不可约。

$A$ 的 Perron-Frobenius 特征值决定种群的渐近增长或衰减速率。如果 $\lambda_\mathrm{pf}<1$，种群以 $\lambda_\mathrm{pf}^t$ 的速率收敛到零，半衰期为 $-1/\log_2\lambda_\mathrm{pf}$ 小时。如果 $\lambda_\mathrm{pf}>1$，种群以 $\lambda_\mathrm{pf}^t$ 几何增长，倍增时间为 $1/\log_2\lambda_\mathrm{pf}$ 小时。最小化 $A$ 的谱半径对应于寻找种群最快衰减率或最慢增长率。

作为矩阵 $A$ 所依赖的底层变量，取环境中影响细菌出生率和存活率的两种化学物质浓度 $c_1,c_2$。把出生率和存活率建模为两种浓度的单项式函数：

$$
\begin{aligned}
b_i&=b_i^\mathrm{nom}(c_1/c_1^\mathrm{nom})^{\alpha_i}
\, (c_2/c_2^\mathrm{nom})^{\beta_i},\quad i=1,\ldots,4,\\
s_i&=s_i^\mathrm{nom}(c_1/c_1^\mathrm{nom})^{\gamma_i}
\, (c_2/c_2^\mathrm{nom})^{\delta_i},\quad i=1,\ldots,3.
\end{aligned}
$$

这里 $b_i^\mathrm{nom}$ 是标称出生率，$s_i^\mathrm{nom}$ 是标称存活率，$c_i^\mathrm{nom}$ 是化学物质 $i$ 的标称浓度。常数 $\alpha_i,\beta_i,\gamma_i,\delta_i$ 给出化学物质浓度偏离标称值时对出生率和存活率的影响。例如，$\alpha_2=-0.3$ 且 $\gamma_1=0.5$ 表示化学物质 1 的浓度高于标称浓度时，会降低 1 到 2 小时龄细菌的出生率，并提高 0 到 1 小时龄细菌的存活率。

假设浓度 $c_1,c_2$ 可以通过给药独立增加或降低（例如在 2 倍范围内），提出寻找药物混合以最大化种群衰减率（即最小化 $\lambda_\mathrm{pf}(A)$）的问题。使用上述方法，该问题可表述为 GP：

$$
\begin{array}{ll}
\text{minimize} & \lambda\\
\text{subject to} & b_1v_1+b_2v_2+b_3v_3+b_4v_4\le\lambda v_1\\
& s_1v_1\le\lambda v_2\\
& s_2v_2\le\lambda v_3\\
& s_3v_3\le\lambda v_4\\
& 1/2\le c_i/c_i^\mathrm{nom}\le2,\quad i=1,2\\
& b_i=b_i^\mathrm{nom}(c_1/c_1^\mathrm{nom})^{\alpha_i}
\, (c_2/c_2^\mathrm{nom})^{\beta_i},\quad i=1,\ldots,4\\
& s_i=s_i^\mathrm{nom}(c_1/c_1^\mathrm{nom})^{\gamma_i}
\, (c_2/c_2^\mathrm{nom})^{\delta_i},\quad i=1,\ldots,3,
\end{array}
$$

变量为 $b_i,s_i,c_i,v_i,\lambda$。

## 4.6 广义不等式约束

标准形式凸优化问题 (4.15) 的一个非常有用的推广，是允许不等式约束函数取向量值，并在约束中使用广义不等式：

$$
\begin{array}{ll}
\text{minimize} & f_0(x)\\
\text{subject to} & f_i(x)\preceq_{K_i}0,\quad i=1,\ldots,m\\
& Ax=b,
\end{array}
\tag{4.48}
$$

其中 $f_0:\mathbf{R}^n\to\mathbf{R}$，$K_i\subseteq\mathbf{R}^{k_i}$ 是正常锥，$f_i:\mathbf{R}^n\to\mathbf{R}^{k_i}$ 是 $K_i$-凸的。我们称该问题为带广义不等式约束的（标准形式）凸优化问题。问题 (4.15) 是特例，其中 $K_i=\mathbf{R}_+,\ i=1,\ldots,m$。

普通凸优化问题的许多结果同样适用于带广义不等式的问题。例如：

- 可行集、任意下水平集和最优集都是凸的。
- 对问题 (4.48) 局部最优的任意点也是全局最优的。
- §4.2.3 中给出的可微 $f_0$ 的最优性条件不需任何改变仍然成立。

第 11 章还将看到，带广义不等式约束的凸优化问题通常可以像普通凸优化问题一样容易求解。

### 4.6.1 锥形式问题

带广义不等式的最简单凸优化问题之一是锥形式问题（或锥规划），它有线性目标和一个仿射（因而 $K$-凸）不等式约束函数：

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & Fx+g\preceq_K0\\
& Ax=b.
\end{array}
\tag{4.49}
$$

当 $K$ 是非负正交锥时，锥形式问题退化为线性规划。可以把锥形式问题看作线性规划的推广，其中逐分量不等式被广义线性不等式取代。

延续与线性规划的类比，我们把

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & x\succeq_K0\\
& Ax=b
\end{array}
$$

称为标准形式的锥形式问题。类似地，问题

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & Fx+g\preceq_K0
\end{array}
$$

称为不等式形式的锥形式问题。

### 4.6.2 半定规划

当 $K$ 是 $\mathbf{S}_+^k$，即 $k\times k$ 正半定矩阵锥时，相关锥形式问题称为半定规划（semidefinite program, SDP），具有形式

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & x_1F_1+\cdots+x_nF_n+G\preceq0\\
& Ax=b,
\end{array}
\tag{4.50}
$$

其中 $G,F_1,\ldots,F_n\in\mathbf{S}^k$，$A\in\mathbf{R}^{p\times n}$。这里的不等式是线性矩阵不等式（见例 2.10）。

如果矩阵 $G,F_1,\ldots,F_n$ 全部为对角矩阵，则 (4.50) 中的 LMI 等价于一组线性不等式，SDP (4.50) 退化为线性规划。

#### 标准形式和不等式形式的半定规划

按与 LP 类似的方式，标准形式 SDP 有线性等式约束，以及变量 $X\in\mathbf{S}^n$ 的（矩阵）非负性约束：

$$
\begin{array}{ll}
\text{minimize} & \operatorname{tr}(CX)\\
\text{subject to} & \operatorname{tr}(A_iX)=b_i,\quad i=1,\ldots,p\\
& X\succeq0,
\end{array}
\tag{4.51}
$$

其中 $C,A_1,\ldots,A_p\in\mathbf{S}^n$。（回忆 $\operatorname{tr}(CX)=\sum_{i,j=1}^n C_{ij}X_{ij}$ 是 $\mathbf{S}^n$ 上一般实值线性函数的形式。）这个形式应与标准形式线性规划 (4.28) 比较。在 LP 和 SDP 标准形式中，我们都在线性等式约束和变量非负性约束下最小化变量的线性函数。

与不等式形式 LP (4.29) 类似，不等式形式 SDP 没有等式约束，并有一个 LMI：

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & x_1A_1+\cdots+x_nA_n\preceq B,
\end{array}
$$

变量为 $x\in\mathbf{R}^n$，参数为 $B,A_1,\ldots,A_n\in\mathbf{S}^k$、$c\in\mathbf{R}^n$。

#### 多个 LMI 和线性不等式

具有线性目标、线性等式和不等式约束以及若干 LMI 约束的问题，也常被称为 SDP：

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} &
F^{(i)}(x)=x_1F^{(i)}_1+\cdots+x_nF^{(i)}_n+G^{(i)}\preceq0,\quad i=1,\ldots,K\\
& Gx\preceq h,\quad Ax=b.
\end{array}
$$

这类问题很容易转换为 SDP，只需把各个 LMI 和线性不等式组成一个大的块对角 LMI：

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & \operatorname{diag}(Gx-h,F^{(1)}(x),\ldots,F^{(K)}(x))\preceq0\\
& Ax=b.
\end{array}
$$

### 4.6.3 例子

#### 二阶锥规划

SOCP (4.36) 可表示为锥形式问题

$$
\begin{array}{ll}
\text{minimize} & c^Tx\\
\text{subject to} & -(A_ix+b_i,c_i^Tx+d_i)\preceq_{K_i}0,\quad i=1,\ldots,m\\
& Fx=g,
\end{array}
$$

其中

$$
K_i=\{(y,t)\in\mathbf{R}^{n_i+1}\mid \|y\|_2\le t\},
$$

即 $\mathbf{R}^{n_i+1}$ 中的二阶锥。这解释了优化问题 (4.36) 被称为二阶锥规划的原因。

#### 矩阵范数最小化

令 $A(x)=A_0+x_1A_1+\cdots+x_nA_n$，其中 $A_i\in\mathbf{R}^{p\times q}$。考虑无约束问题

$$
\text{minimize}\quad \|A(x)\|_2,
$$

其中 $\|\cdot\|_2$ 表示谱范数（最大奇异值），$x\in\mathbf{R}^n$ 是变量。由于 $\|A(x)\|_2$ 是 $x$ 的凸函数，这是一个凸问题。

利用 $\|A\|_2\le s$ 当且仅当 $A^TA\preceq s^2I$（且 $s\ge0$），该问题可表示为

$$
\begin{array}{ll}
\text{minimize} & s\\
\text{subject to} & A(x)^TA(x)\preceq sI,
\end{array}
$$

变量为 $x,s$。由于函数 $A(x)^TA(x)-sI$ 关于 $(x,s)$ 矩阵凸，这是一个带单个 $q\times q$ 矩阵不等式约束的凸优化问题。

也可以使用一个大小为 $(p+q)\times(p+q)$ 的线性矩阵不等式来表述该问题，利用事实

$$
A^TA\preceq t^2I\ \text{且 }t\ge0
\quad\Longleftrightarrow\quad
\begin{bmatrix}
tI & A\\
A^T & tI
\end{bmatrix}\succeq0
$$

（见 §A.5.5）。于是得到 SDP

$$
\begin{array}{ll}
\text{minimize} & t\\
\text{subject to} &
\begin{bmatrix}
tI & A(x)\\
A(x)^T & tI
\end{bmatrix}\succeq0,
\end{array}
$$

变量为 $x,t$。

#### 矩问题

令 $t$ 是 $\mathbf{R}$ 上的随机变量。期望值 $\mathbf{E}t^k$（假设存在）称为 $t$ 分布的（幂）矩。以下经典结果给出了矩序列的刻画。

如果存在 $\mathbf{R}$ 上的概率分布，使得 $x_k=\mathbf{E}t^k,\ k=0,\ldots,2n$，则 $x_0=1$ 且

$$
H(x_0,\ldots,x_{2n})=
\begin{bmatrix}
x_0 & x_1 & x_2 & \cdots & x_{n-1} & x_n\\
x_1 & x_2 & x_3 & \cdots & x_n & x_{n+1}\\
x_2 & x_3 & x_4 & \cdots & x_{n+1} & x_{n+2}\\
\vdots & \vdots & \vdots & & \vdots & \vdots\\
x_{n-1} & x_n & x_{n+1} & \cdots & x_{2n-2} & x_{2n-1}\\
x_n & x_{n+1} & x_{n+2} & \cdots & x_{2n-1} & x_{2n}
\end{bmatrix}\succeq0.
\tag{4.52}
$$

矩阵 $H$ 称为与 $x_0,\ldots,x_{2n}$ 关联的 Hankel 矩阵。这一点很容易看出：令 $x_i=\mathbf{E}t^i,\ i=0,\ldots,2n$ 是某分布的矩，且 $y=(y_0,y_1,\ldots,y_n)\in\mathbf{R}^{n+1}$。则

$$
y^TH(x_0,\ldots,x_{2n})y
=\sum_{i,j=0}^n y_iy_j\mathbf{E}t^{i+j}
=\mathbf{E}(y_0+y_1t+\cdots+y_nt^n)^2\ge0.
$$

下面的部分逆命题不那么显然：若 $x_0=1$ 且 $H(x)\succ0$，则存在 $\mathbf{R}$ 上的概率分布，使得 $x_i=\mathbf{E}t^i,\ i=0,\ldots,2n$。（证明见习题 2.37。）现在假设 $x_0=1$ 且 $H(x)\succeq0$（但可能 $H(x)\not\succ0$），即线性矩阵不等式 (4.52) 成立但可能不严格。此时存在 $\mathbf{R}$ 上的一列分布，其矩收敛到 $x$。总结起来：$x_0,\ldots,x_{2n}$ 是某个 $\mathbf{R}$ 上分布的矩（或一列分布矩的极限）这一条件，可以表示为变量 $x$ 的线性矩阵不等式 (4.52)，以及线性等式 $x_0=1$。利用这一事实，可以把一些涉及矩的有趣问题表述为 SDP。

假设 $t$ 是 $\mathbf{R}$ 上的随机变量。我们不知道其分布，但知道一些矩的界：

$$
\underline\mu_k\le \mathbf{E}t^k\le \overline\mu_k,\quad k=1,\ldots,2n
$$

（这包括某些矩的精确值已知作为特例）。令 $p(t)=c_0+c_1t+\cdots+c_{2n}t^{2n}$ 是给定多项式。$p(t)$ 的期望是矩 $\mathbf{E}t^i$ 的线性函数：

$$
\mathbf{E}p(t)=\sum_{i=0}^{2n}c_i\mathbf{E}t^i
=\sum_{i=0}^{2n}c_ix_i.
$$

在所有满足给定矩界的概率分布上，计算 $\mathbf{E}p(t)$ 的上下界，可以通过求解 SDP

$$
\begin{array}{ll}
\text{minimize (maximize)} & c_1x_1+\cdots+c_{2n}x_{2n}\\
\text{subject to} & \underline\mu_k\le x_k\le\overline\mu_k,\quad k=1,\ldots,2n\\
& H(1,x_1,\ldots,x_{2n})\succeq0
\end{array}
$$

来完成，变量为 $x_1,\ldots,x_{2n}$。这给出所有满足已知矩约束的概率分布上 $\mathbf{E}p(t)$ 的界。该界是紧的：存在一列分布，其矩满足给定矩界，并且 $\mathbf{E}p(t)$ 收敛到这些 SDP 给出的上下界。

#### 用不完整协方差信息界定投资组合风险

再次考虑经典 Markowitz 投资组合问题的设置（见第 155 页）。有 $n$ 个资产或股票的投资组合，$x_i$ 表示在某投资期内持有的资产 $i$ 的数量，$p_i$ 表示资产 $i$ 在该时期的相对价格变化。投资组合总价值变化为 $p^Tx$。价格变化向量 $p$ 建模为随机向量，均值和协方差为

$$
\bar p=\mathbf{E}p,\qquad
\Sigma=\mathbf{E}(p-\bar p)(p-\bar p)^T.
$$

因此，投资组合价值变化是一个均值为 $\bar p^Tx$、标准差为 $\sigma=(x^T\Sigma x)^{1/2}$ 的随机变量。大额损失风险，即投资组合价值变化显著低于期望值的风险，与标准差 $\sigma$ 直接相关，并随其增大而增大。因此，标准差 $\sigma$（或方差 $\sigma^2$）被用作投资组合风险的度量。

在经典投资组合优化问题中，投资组合 $x$ 是优化变量，目标是在最低均值收益和其他约束下最小化风险。价格变化统计量 $\bar p$ 和 $\Sigma$ 是已知问题参数。这里考虑的风险界定问题把问题反过来：假设投资组合 $x$ 已知，但关于协方差矩阵 $\Sigma$ 只知道部分信息。例如，可能有每个元素的上下界：

$$
L_{ij}\le\Sigma_{ij}\le U_{ij},\quad i,j=1,\ldots,n,
$$

其中 $L,U$ 给定。现在提出问题：在所有与给定界一致的协方差矩阵中，我们的投资组合最大风险是多少？定义投资组合最坏情况方差为

$$
\sigma_\mathrm{wc}^2=
\sup\{x^T\Sigma x\mid L_{ij}\le\Sigma_{ij}\le U_{ij},\ i,j=1,\ldots,n,\ \Sigma\succeq0\}.
$$

加入条件 $\Sigma\succeq0$，这是协方差矩阵当然必须满足的。

可以通过求解 SDP

$$
\begin{array}{ll}
\text{maximize} & x^T\Sigma x\\
\text{subject to} & L_{ij}\le\Sigma_{ij}\le U_{ij},\quad i,j=1,\ldots,n\\
& \Sigma\succeq0
\end{array}
$$

来找到 $\sigma_\mathrm{wc}$，变量为 $\Sigma\in\mathbf{S}^n$（问题参数为 $x,L,U$）。最优 $\Sigma$ 是与给定元素界一致的最坏协方差矩阵，其中“最坏”指对给定组合 $x$ 风险最大。由 SDP 的最优 $\Sigma$，可以很容易构造一个与给定界一致并达到最坏情况方差的 $p$ 分布。例如，可以取 $p=\bar p+\Sigma^{1/2}v$，其中 $v$ 是任意满足 $\mathbf{E}v=0$ 和 $\mathbf{E}vv^T=I$ 的随机向量。

显然，对于任何关于 $\Sigma$ 的凸先验信息，都可用同样方法确定 $\sigma_\mathrm{wc}$。下面列出一些例子。

- 已知某些投资组合的方差。可能有等式约束

$$
u_k^T\Sigma u_k=\sigma_k^2,
$$

其中 $u_k,\sigma_k$ 给定。这对应于先验知道某些已知投资组合（由 $u_k$ 给出）具有已知（或估计非常准确）的方差。

- 包含估计误差影响。如果协方差 $\Sigma$ 由经验数据估计，估计方法会给出估计值 $\hat\Sigma$ 以及关于估计可靠性的一些信息，例如置信椭球。这可表示为

$$
C(\Sigma-\hat\Sigma)\le\alpha,
$$

其中 $C$ 是 $\mathbf{S}^n$ 上的正定二次型，常数 $\alpha$ 确定置信水平。

- 因子模型。协方差可能具有形式

$$
\Sigma=F\Sigma_\mathrm{factor}F^T+D,
$$

其中 $F\in\mathbf{R}^{n\times k}$，$\Sigma_\mathrm{factor}\in\mathbf{S}^k$，$D$ 为对角矩阵。这对应于价格变化模型

$$
p=Fz+d,
$$

其中 $z$ 是随机变量（影响价格变化的底层因子），$d_i$ 相互独立（每个资产价格的额外波动）。假设因子已知。由于 $\Sigma$ 与 $\Sigma_\mathrm{factor}$ 和 $D$ 线性相关，可以对它们施加任何凸约束（表示先验信息），并仍然用凸优化计算 $\sigma_\mathrm{wc}$。

- 相关系数信息。最简单情况下，$\Sigma$ 的对角元素（即每个资产价格的波动率）已知，并且价格变化之间相关系数的界已知：

$$
l_{ij}\le \rho_{ij}=
\frac{\Sigma_{ij}}{\Sigma_{ii}^{1/2}\Sigma_{jj}^{1/2}}
\le u_{ij},\quad i,j=1,\ldots,n.
$$

由于 $\Sigma_{ii}$ 已知，而 $i\ne j$ 时 $\Sigma_{ij}$ 未知，这些是线性不等式。

#### 图上的最快混合 Markov 链

考虑一个无向图，节点为 $1,\ldots,n$，边集为

$$
E\subseteq\{1,\ldots,n\}\times\{1,\ldots,n\}.
$$

这里 $(i,j)\in E$ 表示节点 $i$ 和 $j$ 之间有边。由于图是无向的，$E$ 是对称的：$(i,j)\in E$ 当且仅当 $(j,i)\in E$。允许自环，即可能有 $(i,i)\in E$。

定义一个 Markov 链，状态为 $X(t)\in\{1,\ldots,n\}$，其中 $t\in\mathbf{Z}_+$（非负整数集合）。对每条边 $(i,j)\in E$，关联一个概率 $P_{ij}$，表示 $X$ 在节点 $i$ 和 $j$ 之间转移的概率。状态转移只能沿边发生；对 $(i,j)\notin E$ 有 $P_{ij}=0$。与边相关的概率必须非负，并且对每个节点，与该节点相连的边（包括自环，如果有）概率之和必须为 1。

该 Markov 链的转移概率矩阵为

$$
P_{ij}=\operatorname{prob}(X(t+1)=i\mid X(t)=j),\quad i,j=1,\ldots,n.
$$

该矩阵必须满足

$$
P_{ij}\ge0,\quad i,j=1,\ldots,n,\qquad
\mathbf{1}^TP=\mathbf{1}^T,\quad P=P^T,
\tag{4.53}
$$

并且

$$
P_{ij}=0\quad \text{for }(i,j)\notin E.
\tag{4.54}
$$

由于 $P$ 对称且 $\mathbf{1}^TP=\mathbf{1}^T$，可得 $P\mathbf{1}=\mathbf{1}$，所以均匀分布 $(1/n)\mathbf{1}$ 是该 Markov 链的平衡分布。$X(t)$ 的分布收敛到 $(1/n)\mathbf{1}$ 的速度由 $P$ 的第二大（按模）特征值决定，即由 $r=\max\{\lambda_2,-\lambda_n\}$ 决定，其中

$$
1=\lambda_1\ge\lambda_2\ge\cdots\ge\lambda_n
$$

是 $P$ 的特征值。称 $r$ 为 Markov 链的混合率。如果 $r=1$，则 $X(t)$ 的分布不一定收敛到 $(1/n)\mathbf{1}$（这意味着 Markov 链不混合）。当 $r<1$ 时，$X(t)$ 的分布以 $r^t$ 的渐近速率趋近 $(1/n)\mathbf{1}$（当 $t\to\infty$）。因此，$r$ 越小，Markov 链混合越快。

最快混合 Markov 链问题是在约束 (4.53) 和 (4.54) 下寻找使 $r$ 最小的 $P$。（问题数据是图，即 $E$。）下面说明该问题可表述为 SDP。

由于特征值 $\lambda_1=1$ 对应特征向量 $\mathbf{1}$，混合率可表示为矩阵 $P$ 限制在子空间 $\mathbf{1}^\perp$ 上的范数：$r=\|QPQ\|_2$，其中 $Q=I-(1/n)\mathbf{1}\mathbf{1}^T$ 是到 $\mathbf{1}^\perp$ 的正交投影矩阵。利用性质 $P\mathbf{1}=\mathbf{1}$，有

$$
\begin{aligned}
r&=\|QPQ\|_2\\
&=\|(I-(1/n)\mathbf{1}\mathbf{1}^T)P(I-(1/n)\mathbf{1}\mathbf{1}^T)\|_2\\
&=\|P-(1/n)\mathbf{1}\mathbf{1}^T\|_2.
\end{aligned}
$$

这表明混合率 $r$ 是 $P$ 的凸函数，所以最快混合 Markov 链问题可写成凸优化问题

$$
\begin{array}{ll}
\text{minimize} & \|P-(1/n)\mathbf{1}\mathbf{1}^T\|_2\\
\text{subject to} & P\mathbf{1}=\mathbf{1}\\
& P_{ij}\ge0,\quad i,j=1,\ldots,n\\
& P_{ij}=0\quad \text{for }(i,j)\notin E,
\end{array}
$$

变量为 $P\in\mathbf{S}^n$。引入标量变量 $t$ 来界定 $\|P-(1/n)\mathbf{1}\mathbf{1}^T\|_2$，可把问题表示为 SDP：

$$
\begin{array}{ll}
\text{minimize} & t\\
\text{subject to} & -tI\preceq P-(1/n)\mathbf{1}\mathbf{1}^T\preceq tI\\
& P\mathbf{1}=\mathbf{1}\\
& P_{ij}\ge0,\quad i,j=1,\ldots,n\\
& P_{ij}=0\quad \text{for }(i,j)\notin E.
\end{array}
\tag{4.55}
$$

## 4.7 向量优化

### 4.7.1 一般和凸向量优化问题

在 §4.6 中，我们把标准形式问题 (4.1) 扩展为包含向量值约束函数。本节研究向量值目标函数的含义。把一般向量优化问题记为

$$
\begin{array}{ll}
\text{minimize (with respect to }K\text{)} & f_0(x)\\
\text{subject to} & f_i(x)\le0,\quad i=1,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p.
\end{array}
\tag{4.56}
$$

这里 $x\in\mathbf{R}^n$ 是优化变量，$K\subseteq\mathbf{R}^q$ 是正常锥，$f_0:\mathbf{R}^n\to\mathbf{R}^q$ 是目标函数，$f_i:\mathbf{R}^n\to\mathbf{R}$ 是不等式约束函数，$h_i:\mathbf{R}^n\to\mathbf{R}$ 是等式约束函数。该问题与标准优化问题 (4.1) 的唯一区别在于：这里目标函数取值于 $\mathbf{R}^q$，且问题描述中包含一个正常锥 $K$，用于比较目标值。在向量优化语境中，标准优化问题 (4.1) 有时称为标量优化问题。

如果目标函数 $f_0$ 是 $K$-凸的，不等式约束函数 $f_1,\ldots,f_m$ 是凸的，等式约束函数 $h_1,\ldots,h_p$ 是仿射的，则称向量优化问题 (4.56) 为凸向量优化问题。（和标量情形一样，通常把等式约束写成 $Ax=b$，其中 $A\in\mathbf{R}^{p\times n}$。）

如何理解向量优化问题 (4.56)？设 $x,y$ 是两个可行点。它们对应的目标值 $f_0(x)$ 和 $f_0(y)$ 要用广义不等式 $\preceq_K$ 比较。我们把 $f_0(x)\preceq_K f_0(y)$ 解释为：$x$ 的值“好于或等于” $y$（由目标 $f_0$ 相对于 $K$ 判断）。向量优化令人困惑之处在于，两个目标值 $f_0(x)$ 和 $f_0(y)$ 不一定可比较；可能既没有 $f_0(x)\preceq_K f_0(y)$，也没有 $f_0(y)\preceq_K f_0(x)$，即二者谁也不优于另一个。这在标量目标优化问题中不会发生。

### 4.7.2 最优点和最优值

先考虑一个特例，其中向量优化问题的含义很清楚。考虑可行点的目标值集合

$$
O=\{f_0(x)\mid \exists x\in D,\ f_i(x)\le0,\ i=1,\ldots,m,\ h_i(x)=0,\ i=1,\ldots,p\}\subseteq\mathbf{R}^q,
$$

称为可达到目标值集合。如果该集合有最小元素（见 §2.4.2），即存在可行 $x$ 使得对所有可行 $y$ 有 $f_0(x)\preceq_K f_0(y)$，则称 $x$ 是问题 (4.56) 的最优点，并称 $f_0(x)$ 为该问题的最优值。（当向量优化问题有最优值时，最优值唯一。）如果 $x^\star$ 是最优点，则 $x^\star$ 处的目标 $f_0(x^\star)$ 可以与每个其他可行点的目标比较，并且好于或等于它。粗略地说，在可行点中，$x^\star$ 明确是 $x$ 的最佳选择。

点 $x^\star$ 最优当且仅当它可行且

$$
O\subseteq f_0(x^\star)+K.
\tag{4.57}
$$

集合 $f_0(x^\star)+K$ 可解释为差于或等于 $f_0(x^\star)$ 的值集合，所以条件 (4.57) 表示每个可达到值都落在该集合中。这如图 4.7 所示。大多数向量优化问题没有最优点和最优值，但某些特殊情况下会有。

**图 4.7** 目标值在 $\mathbf{R}^2$ 中、锥 $K=\mathbf{R}_+^2$ 的向量优化问题的可达到值集合 $O$（阴影区域）。这里标记为 $f_0(x^\star)$ 的点是问题的最优值，$x^\star$ 是最优点。目标值 $f_0(x^\star)$ 可与每个其他可达到值 $f_0(y)$ 比较，并好于或等于它。（这里“好于或等于”表示“位于左下方”。）浅色区域是 $f_0(x^\star)+K$，即所有对应于比 $f_0(x^\star)$ 更差（或相等）目标值的 $z\in\mathbf{R}^2$。

**例 4.9 最佳线性无偏估计。** 假设 $y=Ax+v$，其中 $v\in\mathbf{R}^m$ 是测量噪声，$y\in\mathbf{R}^m$ 是测量向量，$x\in\mathbf{R}^n$ 是要由测量 $y$ 估计的向量。假设 $A$ 的秩为 $n$，且测量噪声满足 $\mathbf{E}v=0$、$\mathbf{E}vv^T=I$，即其分量零均值且不相关。

$x$ 的线性估计器具有形式 $\hat x=Fy$。若对所有 $x$ 都有 $\mathbf{E}\hat x=x$，即 $FA=I$，则称估计器无偏。无偏估计器的误差协方差为

$$
\mathbf{E}(\hat x-x)(\hat x-x)^T
=\mathbf{E}Fvv^TF^T=FF^T.
$$

目标是寻找误差协方差矩阵“小”的无偏估计器。可以使用矩阵不等式，即相对于 $\mathbf{S}_+^n$，来比较误差协方差。这有如下解释：设 $\hat x_1=F_1y$、$\hat x_2=F_2y$ 是两个无偏估计器。那么第一个估计器至少与第二个一样好，即 $F_1F_1^T\preceq F_2F_2^T$，当且仅当对所有 $c$ 有

$$
\mathbf{E}(c^T\hat x_1-c^Tx)^2
\le
\mathbf{E}(c^T\hat x_2-c^Tx)^2.
$$

换言之，对 $x$ 的任意线性函数，估计器 $F_1$ 给出的估计至少与 $F_2$ 一样好。

寻找 $x$ 的无偏估计器的问题可表示为向量优化问题

$$
\begin{array}{ll}
\text{minimize (w.r.t. }\mathbf{S}_+^n\text{)} & FF^T\\
\text{subject to} & FA=I,
\end{array}
\tag{4.58}
$$

变量为 $F\in\mathbf{R}^{n\times m}$。目标 $FF^T$ 关于 $\mathbf{S}_+^n$ 是凸的，所以问题 (4.58) 是凸向量优化问题。一个简单看法是，对任意固定 $v$，$v^TFF^Tv=\|F^Tv\|_2^2$ 是 $F$ 的凸函数。

一个著名结果是，问题 (4.58) 有最优解，即最小二乘估计器或伪逆

$$
F^\star=A^\dagger=(A^TA)^{-1}A^T.
$$

对任何满足 $FA=I$ 的 $F$，都有 $FF^T\succeq F^\star F^{\star T}$。矩阵

$$
F^\star F^{\star T}=A^\dagger A^{\dagger T}=(A^TA)^{-1}
$$

是问题 (4.58) 的最优值。

### 4.7.3 Pareto 最优点和值

现在考虑更常见的情形：可达到目标值集合没有最小元素，因此问题没有最优点或最优值。在这些情形中，可达到值集合的极小元素起重要作用。如果可行点 $x$ 满足 $f_0(x)$ 是可达到值集合 $O$ 的极小元素，则称 $x$ 是 Pareto 最优的（或有效的）。此时称 $f_0(x)$ 为向量优化问题 (4.56) 的 Pareto 最优值。因此，点 $x$ 是 Pareto 最优的，如果它可行，并且对任意可行 $y$，$f_0(y)\preceq_K f_0(x)$ 都推出 $f_0(y)=f_0(x)$。换言之：任何好于或等于 $x$ 的可行点 $y$（即 $f_0(y)\preceq_K f_0(x)$）都具有与 $x$ 完全相同的目标值。

点 $x$ 是 Pareto 最优的，当且仅当它可行且

$$
(f_0(x)-K)\cap O=\{f_0(x)\}.
\tag{4.59}
$$

集合 $f_0(x)-K$ 可解释为好于或等于 $f_0(x)$ 的值集合，所以条件 (4.59) 表示唯一好于或等于 $f_0(x)$ 的可达到值就是 $f_0(x)$ 本身。这如图 4.8 所示。

一个向量优化问题可以有许多 Pareto 最优值（和点）。Pareto 最优值集合记为 $\mathcal{P}$，满足

$$
\mathcal{P}\subseteq O\cap\operatorname{bd}O,
$$

即每个 Pareto 最优值都是位于可达到目标值集合边界上的可达到目标值（见习题 4.52）。

**图 4.8** 目标值在 $\mathbf{R}^2$ 中、锥 $K=\mathbf{R}_+^2$ 的向量优化问题的可达到值集合 $O$（阴影区域）。该问题没有最优点或最优值，但有一组 Pareto 最优点，对应值显示为 $O$ 左下边界上的深色曲线。标记为 $f_0(x^\mathrm{po})$ 的点是 Pareto 最优值，$x^\mathrm{po}$ 是 Pareto 最优点。浅色区域是 $f_0(x^\mathrm{po})-K$，即所有对应于比 $f_0(x^\mathrm{po})$ 更好（或相等）目标值的 $z\in\mathbf{R}^2$。

### 4.7.4 标量化

标量化是寻找向量优化问题 Pareto 最优点（或最优点）的标准技术，基于 §2.6.3 中通过对偶广义不等式刻画最小点和极小点的结果。选择任意 $\lambda\succ_{K^*}0$，即在对偶广义不等式中为正的任意向量。考虑标量优化问题

$$
\begin{array}{ll}
\text{minimize} & \lambda^Tf_0(x)\\
\text{subject to} & f_i(x)\le0,\quad i=1,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p,
\end{array}
\tag{4.60}
$$

并令 $x$ 是其最优点。则 $x$ 是向量优化问题 (4.56) 的 Pareto 最优点。这可由 §2.6.3 中给出的极小点对偶不等式刻画推出，也很容易直接证明。若 $x$ 不是 Pareto 最优，则存在可行 $y$，满足 $f_0(y)\preceq_K f_0(x)$ 且 $f_0(x)\ne f_0(y)$。由于 $f_0(x)-f_0(y)\succeq_K0$ 且非零，有 $\lambda^T(f_0(x)-f_0(y))>0$，即 $\lambda^Tf_0(x)>\lambda^Tf_0(y)$。这与 $x$ 是标量问题 (4.60) 的最优点矛盾。

使用标量化，可以通过求解普通标量优化问题 (4.60) 来寻找任意向量优化问题的 Pareto 最优点。向量 $\lambda$ 有时称为权重向量，必须满足 $\lambda\succ_{K^*}0$。权重向量是自由参数；改变它会得到向量优化问题 (4.56) 的（可能）不同 Pareto 最优解。如图 4.9 所示。该图也显示了一个无法通过任何权重向量 $\lambda\succ_{K^*}0$ 的标量化得到的 Pareto 最优点。

**图 4.9** 标量化。锥 $K=\mathbf{R}_+^2$ 的向量优化问题的可达到值集合 $O$。图中显示三个 Pareto 最优值 $f_0(x_1),f_0(x_2),f_0(x_3)$。前两个值可由标量化得到：$f_0(x_1)$ 在所有 $u\in O$ 上最小化 $\lambda_1^Tu$，$f_0(x_2)$ 最小化 $\lambda_2^Tu$，其中 $\lambda_1,\lambda_2\succ0$。值 $f_0(x_3)$ 是 Pareto 最优的，但不能由标量化得到。

标量化方法可以几何解释。点 $x$ 对标量化问题最优，即在可行集上最小化 $\lambda^Tf_0$，当且仅当对所有可行 $y$ 有 $\lambda^T(f_0(y)-f_0(x))\ge0$。这等价于说 $\{u\mid -\lambda^T(u-f_0(x))=0\}$ 是可达到目标值集合 $O$ 在点 $f_0(x)$ 处的支撑超平面；特别地，

$$
\{u\mid \lambda^T(u-f_0(x))<0\}\cap O=\emptyset.
\tag{4.61}
$$

因此，当找到标量化问题的最优点时，不仅找到了原向量优化问题的 Pareto 最优点，还找到了 $\mathbf{R}^q$ 中由 (4.61) 给出的整个不可达到目标值半空间。

#### 凸向量优化问题的标量化

现在假设向量优化问题 (4.56) 是凸的。则标量化问题 (4.60) 也是凸的，因为 $\lambda^Tf_0$ 是一个标量值凸函数（由 §3.6 的结果）。这意味着可以通过求解凸标量优化问题来寻找凸向量优化问题的 Pareto 最优点。对权重向量 $\lambda\succ_{K^*}0$ 的每个选择，通常得到不同的 Pareto 最优点。

对凸向量优化问题，有一个部分逆命题：对每个 Pareto 最优点 $x^\mathrm{po}$，存在某个非零 $\lambda\succeq_{K^*}0$，使得 $x^\mathrm{po}$ 是标量化问题 (4.60) 的解。因此，粗略地说，对凸问题，当权重向量 $\lambda$ 在 $K^*$-非负非零值上变化时，标量化方法产生所有 Pareto 最优点。这里需要小心，因为并非标量化问题在 $\lambda\succeq_{K^*}0$ 且 $\lambda\ne0$ 下的每个解都是向量问题的 Pareto 最优点。（相比之下，当 $\lambda\succ_{K^*}0$ 时，标量化问题的每个解都是 Pareto 最优点。）

在某些情况下，可以用这个部分逆命题找到凸向量优化问题的所有 Pareto 最优点。用 $\lambda\succ_{K^*}0$ 标量化得到一组 Pareto 最优点（对非凸向量优化问题也一样）。为了找到剩余的 Pareto 最优解，必须考虑满足 $\lambda\succeq_{K^*}0$ 的非零权重向量。对每个这样的权重向量，先确定标量化问题的所有解；然后在这些解中检查哪些实际上是向量优化问题的 Pareto 最优点。这些“极端”Pareto 最优点也可以作为由正权重向量得到的 Pareto 最优点的极限找到。

为证明这个部分逆命题，考虑集合

$$
A=O+K=\{t\in\mathbf{R}^q\mid f_0(x)\preceq_K t\ \text{对某个可行 }x\},
\tag{4.62}
$$

它由所有差于或等于（相对于 $\preceq_K$）某个可达到目标值的值组成。当问题是凸的时，虽然可达到目标值集合 $O$ 不一定凸，但集合 $A$ 是凸的。此外，$A$ 的极小元素与可达到值集合 $O$ 的极小元素完全相同，即与 Pareto 最优值相同。（见习题 4.53。）现在使用 §2.6.3 的结果可知，$A$ 的任意极小元素都会对某个非零 $\lambda\succeq_{K^*}0$ 最小化 $A$ 上的 $\lambda^Tz$。这意味着向量优化问题的每个 Pareto 最优点，都是某个非零权重 $\lambda\succeq_{K^*}0$ 下标量化问题的最优点。

**例 4.10 一组矩阵的极小上界。** 考虑相对于正半定锥的凸向量优化问题

$$
\begin{array}{ll}
\text{minimize (w.r.t. }\mathbf{S}_+^n\text{)} & X\\
\text{subject to} & X\succeq A_i,\quad i=1,\ldots,m,
\end{array}
\tag{4.63}
$$

其中 $A_i\in\mathbf{S}^n,\ i=1,\ldots,m$ 给定。约束意味着 $X$ 是给定矩阵 $A_1,\ldots,A_m$ 的上界；(4.63) 的 Pareto 最优解是这些矩阵的极小上界。

为寻找 Pareto 最优点，应用标量化：选择任意 $W\in\mathbf{S}_{++}^n$，并形成问题

$$
\begin{array}{ll}
\text{minimize} & \operatorname{tr}(WX)\\
\text{subject to} & X\succeq A_i,\quad i=1,\ldots,m,
\end{array}
\tag{4.64}
$$

这是一个 SDP。不同的 $W$ 一般会给出不同的极小解。

部分逆命题告诉我们，若 $X$ 是向量问题 (4.63) 的 Pareto 最优点，则它是某个非零权重矩阵 $W\succeq0$ 下 SDP (4.64) 的最优点。（不过在这个例子中，并非 (4.64) 的每个解都是向量优化问题的 Pareto 最优点。）

该问题有一个简单几何解释。对每个 $A\in\mathbf{S}_{++}^n$，关联一个以原点为中心的椭球

$$
E_A=\{u\mid u^TA^{-1}u\le1\},
$$

于是 $A\preceq B$ 当且仅当 $E_A\subseteq E_B$。问题 (4.63) 的 Pareto 最优点 $X$ 对应于包含与 $A_1,\ldots,A_m$ 关联的椭球的极小椭球。图 4.10 给出了一个例子。

**图 4.10** 问题 (4.63) 的几何解释。三个阴影椭球对应数据 $A_1,A_2,A_3\in\mathbf{S}_{++}^2$；Pareto 最优点对应包含它们的极小椭球。边界标记为 $X_1$ 和 $X_2$ 的两个椭球，显示由两个不同权重矩阵 $W_1,W_2$ 求解 SDP (4.64) 得到的两个极小椭球。

### 4.7.5 多准则优化

当向量优化问题涉及锥 $K=\mathbf{R}_+^q$ 时，称为多准则或多目标优化问题。$f_0$ 的分量记为 $F_1,\ldots,F_q$，可解释为 $q$ 个不同的标量目标，每个都希望最小化。称 $F_i$ 为该问题的第 $i$ 个目标。若 $f_1,\ldots,f_m$ 凸、$h_1,\ldots,h_p$ 仿射，并且目标 $F_1,\ldots,F_q$ 都凸，则多准则优化问题是凸的。

由于多准则问题是向量优化问题，§4.7.1-§4.7.4 的所有内容都适用。不过，对多准则问题，可以给出更具体的解释。如果 $x$ 可行，则可把 $F_i(x)$ 看作按第 $i$ 个目标衡量的得分或值。如果 $x,y$ 都可行，$F_i(x)\le F_i(y)$ 表示按第 $i$ 个目标，$x$ 至少与 $y$ 一样好；$F_i(x)<F_i(y)$ 表示按第 $i$ 个目标，$x$ 优于 $y$，或 $x$ 击败 $y$。如果 $x,y$ 都可行，并且对 $i=1,\ldots,q$ 有 $F_i(x)\le F_i(y)$，且至少对一个 $j$ 有 $F_j(x)<F_j(y)$，则称 $x$ 优于 $y$，或称 $x$ 支配 $y$。粗略地说，$x$ 在所有目标上不差于 $y$，并且在至少一个目标上优于 $y$。

在多准则问题中，最优点 $x^\star$ 满足对每个可行 $y$ 都有

$$
F_i(x^\star)\le F_i(y),\quad i=1,\ldots,q.
$$

换言之，$x^\star$ 同时对每个标量问题

$$
\begin{array}{ll}
\text{minimize} & F_j(x)\\
\text{subject to} & f_i(x)\le0,\quad i=1,\ldots,m\\
& h_i(x)=0,\quad i=1,\ldots,p
\end{array}
$$

最优，其中 $j=1,\ldots,q$。当存在最优点时，称这些目标非竞争，因为目标之间不需折中；即使忽略其他目标，每个目标也都已经尽可能小。

Pareto 最优点 $x^\mathrm{po}$ 满足：若 $y$ 可行且 $F_i(y)\le F_i(x^\mathrm{po})$，$i=1,\ldots,q$，则 $F_i(x^\mathrm{po})=F_i(y)$，$i=1,\ldots,q$。可重述为：一个点是 Pareto 最优的，当且仅当它可行并且不存在更好的可行点。特别地，若一个可行点不是 Pareto 最优的，则至少存在另一个更好的可行点。因此，在寻找好点时，显然可以把搜索限制在 Pareto 最优点上。

#### 折中分析

现在假设 $x$ 和 $y$ 是 Pareto 最优点，并且例如

$$
\begin{aligned}
F_i(x)&<F_i(y),\quad i\in A,\\
F_i(x)&=F_i(y),\quad i\in B,\\
F_i(x)&>F_i(y),\quad i\in C,
\end{aligned}
$$

其中 $A\cup B\cup C=\{1,\ldots,q\}$。换言之，$A$ 是 $x$ 击败 $y$ 的目标指标集合，$B$ 是 $x$ 和 $y$ 打平的目标指标集合，$C$ 是 $y$ 击败 $x$ 的目标指标集合。如果 $A$ 和 $C$ 都为空，则两个点 $x,y$ 具有完全相同的目标值。若不是这种情况，则 $A$ 和 $C$ 都必须非空。换言之，比较两个 Pareto 最优点时，它们要么性能完全相同（所有目标相等），要么彼此都在至少一个目标上击败对方。

比较点 $x$ 和 $y$ 时，我们说用 $i\in A$ 的更好目标值交换了 $i\in C$ 的更差目标值。最优折中分析（或简称折中分析）研究为了在某些目标上做得更好，必须在一个或多个其他目标上做得多差；更一般地，它研究哪些目标值集合是可达到的。

作为例子，考虑双准则问题。设 $x$ 是 Pareto 最优点，目标为 $F_1(x)$ 和 $F_2(x)$。可以问：为了得到一个可行点 $z$ 满足 $F_1(z)\le F_1(x)-a$（其中 $a>0$），$F_2(z)$ 必须增大多少？粗略地说，这是在问：为在第一个目标上改进 $a$，必须在第二个目标上付出多大代价。如果为了实现 $F_1$ 的小幅下降必须接受 $F_2$ 的大幅增加，则称在 Pareto 最优值 $(F_1(x),F_2(x))$ 附近两个目标存在强折中。反之，如果 $F_1$ 的大幅下降只需 $F_2$ 的小幅增加，则称两个目标之间的折中较弱。

也可以考虑牺牲第一个目标的性能来改进第二个目标的情况。这里寻找 $F_2(z)$ 可以降低多少，同时得到一个满足 $F_1(z)\le F_1(x)+a$ 的可行点 $z$，其中 $a>0$。此时我们在第二个目标上获得收益，即相对于 $F_2(x)$ 的下降。如果该收益很大（即只需小幅增加 $F_1$ 就能大幅降低 $F_2$），则称目标表现出强折中；如果收益很小，则称在 Pareto 最优值 $(F_1(x),F_2(x))$ 附近目标折中较弱。

#### 最优折中曲面

多准则问题的 Pareto 最优值集合称为最优折中曲面（一般当 $q>2$ 时）或最优折中曲线（当 $q=2$ 时）。（由于接受任何非 Pareto 最优点都是不明智的，可以把折中分析限制在 Pareto 最优点上。）折中分析有时也称为探索最优折中曲面。（最优折中曲面通常但并不总是通常意义上的曲面。例如，如果问题有最优点，则最优折中曲面只包含一个点，即最优值。）

最优折中曲线很容易解释。图 4.11 给出了一个（凸）双准则问题的例子。从该曲线可以直观理解两个目标之间的折中。

- 右端点显示不考虑 $F_1$ 时 $F_2$ 的最小可能值。
- 左端点显示不考虑 $F_2$ 时 $F_1$ 的最小可能值。
- 通过寻找曲线与竖线 $F_1=\alpha$ 的交点，可以看出为了达到 $F_1\le\alpha$，$F_2$ 必须多大。
- 通过寻找曲线与横线 $F_2=\beta$ 的交点，可以看出为了达到 $F_2\le\beta$，$F_1$ 必须多大。
- 曲线上某点（即 Pareto 最优值）处的斜率显示两个目标之间的局部最优折中。斜率陡峭处，$F_1$ 的小变化伴随 $F_2$ 的大变化。
- 曲率大的点是指一个目标的小幅降低只能通过另一个目标的大幅增加实现的点。这就是折中曲线中常说的“膝点”，在许多应用中代表一个良好的折中解。

这些解释都可以简单扩展到折中曲面，不过当目标超过三个时，曲面可视化会很困难。

#### 标量化多准则问题

通过形成加权和目标来标量化多准则问题：

$$
\lambda^Tf_0(x)=\sum_{i=1}^q \lambda_iF_i(x),
$$

其中 $\lambda\succ0$。可以把 $\lambda_i$ 解释为赋予第 $i$ 个目标的权重。权重 $\lambda_i$ 可以被看作量化我们让 $F_i$ 变小的愿望（或对 $F_i$ 很大的反感）。特别地，如果希望 $F_i$ 小，就应取较大的 $\lambda_i$；如果对 $F_i$ 不太在意，则可取较小的 $\lambda_i$。比值 $\lambda_i/\lambda_j$ 可解释为第 $i$ 个目标相对于第 $j$ 个目标的相对权重或相对重要性。也可以把 $\lambda_i/\lambda_j$ 看作两个目标之间的汇率，因为在加权和目标中，$F_i$ 减少 $\alpha$ 与 $F_j$ 增加 $(\lambda_i/\lambda_j)\alpha$ 被视为相同。

这些解释给出在探索最优折中曲面时如何设置或改变权重的直觉。例如，假设权重向量 $\lambda\succ0$ 得到 Pareto 最优点 $x^\mathrm{po}$，其目标值为 $F_1(x^\mathrm{po}),\ldots,F_q(x^\mathrm{po})$。为寻找一个（可能）新的 Pareto 最优点，它用其他目标值（可能）变差来换取第 $k$ 个目标值更好，可构造新权重向量 $\tilde\lambda$，满足

$$
\tilde\lambda_k>\lambda_k,\qquad
\tilde\lambda_j=\lambda_j,\quad j\ne k,\ j=1,\ldots,q,
$$

即增大第 $k$ 个目标的权重。这会得到新的 Pareto 最优点 $\tilde x^\mathrm{po}$，满足 $F_k(\tilde x^\mathrm{po})\le F_k(x^\mathrm{po})$（通常严格小于），即第 $k$ 个目标改进的新 Pareto 最优点。

还可以看出，在最优折中曲面光滑的任意点处，$\lambda$ 给出该 Pareto 最优点处曲面的内法向量。特别地，当选择权重向量 $\lambda$ 并应用标量化时，得到的 Pareto 最优点处，$\lambda$ 给出目标之间的局部折中。

实践中，通常基于上述直觉临时调整权重来探索最优折中曲面。后面（第 5 章）将看到，标量化的基本思想，即最小化目标的加权和，然后调整权重以获得合适解，正是对偶性的核心。

### 4.7.6 例子

#### 正则化最小二乘

给定 $A\in\mathbf{R}^{m\times n}$ 和 $b\in\mathbf{R}^m$，希望选择 $x\in\mathbf{R}^n$，同时考虑两个二次目标：

- $F_1(x)=\|Ax-b\|_2^2=x^TA^TAx-2b^TAx+b^Tb$，度量 $Ax$ 与 $b$ 的失配；
- $F_2(x)=\|x\|_2^2=x^Tx$，度量 $x$ 的大小。

目标是寻找既拟合良好（即 $F_1$ 小）又不太大（即 $F_2$ 小）的 $x$。可以把该问题表述为关于锥 $\mathbf{R}_+^2$ 的向量优化问题，即无约束双准则问题。

通过取 $\lambda_1>0$ 和 $\lambda_2>0$，并最小化标量加权和目标，可对该问题标量化：

$$
\begin{aligned}
\lambda^Tf_0(x)
&=\lambda_1F_1(x)+\lambda_2F_2(x)\\
&=x^T(\lambda_1A^TA+\lambda_2I)x-2\lambda_1b^TAx+\lambda_1b^Tb,
\end{aligned}
$$

得到

$$
x(\mu)=(\lambda_1A^TA+\lambda_2I)^{-1}\lambda_1A^Tb
=(A^TA+\mu I)^{-1}A^Tb,
$$

其中 $\mu=\lambda_2/\lambda_1$。对任意 $\mu>0$，该点都是双准则问题的 Pareto 最优点。可以把 $\mu=\lambda_2/\lambda_1$ 解释为相对于 $F_1$ 赋予 $F_2$ 的相对权重。

该方法产生除两个极限点之外的所有 Pareto 最优点，这两个极限点对应 $\mu\to\infty$ 和 $\mu\to0$。第一种情形得到 Pareto 最优解 $x=0$，它可由权重 $\lambda=(0,1)$ 的标量化得到。另一个极端得到 Pareto 最优解 $A^\dagger b$，其中 $A^\dagger$ 是 $A$ 的伪逆。该 Pareto 最优解作为标量化问题最优解在 $\mu\to0$（即 $\lambda\to(1,0)$）时的极限得到。（§6.3.2 将再次遇到正则化最小二乘问题。）

图 4.11 显示了一个正则化最小二乘问题的最优折中曲线和可达到值集合，其中问题数据 $A\in\mathbf{R}^{100\times10}$、$b\in\mathbf{R}^{100}$。（更多讨论见习题 4.50。）

**图 4.11** 正则化最小二乘问题的最优折中曲线。阴影集合是可达到值集合 $(\|Ax-b\|_2^2,\|x\|_2^2)$。最优折中曲线显示为边界左下方的深色部分。

#### 投资组合优化中的风险-收益折中

第 155 页描述的经典 Markowitz 投资组合优化问题自然可表示为双准则问题，其中目标是负均值收益（因为希望最大化均值收益）和收益方差：

$$
\begin{array}{ll}
\text{minimize (w.r.t. }\mathbf{R}_+^2\text{)} &
(F_1(x),F_2(x))=(-\bar p^Tx,x^T\Sigma x)\\
\text{subject to} & \mathbf{1}^Tx=1,\quad x\succeq0.
\end{array}
$$

形成相应标量化问题时，可以不失一般性地取 $\lambda_1=1$、$\lambda_2=\mu>0$：

$$
\begin{array}{ll}
\text{minimize} & -\bar p^Tx+\mu x^T\Sigma x\\
\text{subject to} & \mathbf{1}^Tx=1,\quad x\succeq0,
\end{array}
$$

这是一个 QP。在这个例子中，也会得到除 $\mu\to0$ 和 $\mu\to\infty$ 两个极限情形外的所有 Pareto 最优投资组合。粗略地说，第一种情形得到不考虑收益方差的最大均值收益；第二种情形得到不考虑均值收益的最小方差收益。假设 $\bar p_k>\bar p_i$ 对 $i\ne k$ 成立，即资产 $k$ 是唯一均值收益最大的资产，则组合 $x=e_k$ 是唯一对应 $\mu\to0$ 的组合。（换言之，把投资组合完全集中在均值收益最大的资产中。）在许多投资组合问题中，资产 $n$ 对应无风险投资，其确定性收益为 $r_\mathrm{rf}$。假设去掉最后一行和最后一列（它们为零）后的 $\Sigma$ 满秩，则另一个极端 Pareto 最优组合为 $x=e_n$，即投资组合完全集中在无风险资产中。

作为具体例子，考虑一个有 4 个资产的简单投资组合优化问题，其价格变化均值和标准差如下表所示。

| 资产 | $\bar p_i$ | $\Sigma_{ii}^{1/2}$ |
| --- | ---: | ---: |
| 1 | 12% | 20% |
| 2 | 10% | 10% |
| 3 | 7% | 5% |
| 4 | 3% | 0% |

资产 4 是无风险资产，具有确定的 3% 收益。资产 3、2、1 的均值收益依次增大，从 7% 到 12%，标准差也依次增大，从 5% 到 20%。资产之间的相关系数为 $\rho_{12}=30\%$、$\rho_{13}=-40\%$ 和 $\rho_{23}=0\%$。

图 4.12 显示该投资组合优化问题的最优折中曲线。图按惯例绘制，横轴为标准差（即方差平方根），纵轴为期望收益。下图显示每个 Pareto 最优点对应的最优资产配置向量 $x$。

**图 4.12** 上图：简单投资组合优化问题的最优风险-收益折中曲线。左端点对应把全部资源投入无风险资产，因此标准差为零。右端点对应把全部资源投入资产 1，它具有最高均值收益。下图：对应的最优配置。

## Bibliography

线性规划自 20 世纪 40 年代以来得到广泛研究，并有许多优秀书籍讨论，包括 Dantzig [Dan63]、Luenberger [Lue84]、Schrijver [Sch86]、Papadimitriou 和 Steiglitz [PS98]、Bertsimas 和 Tsitsiklis [BT97]、Vanderbei [Van96]，以及 Roos、Terlaky 和 Vial [RTV97]。Dantzig 和 Schrijver 还详细叙述了线性规划的历史。近期综述见 Todd [Tod02]。

Schaible [Sch82, Sch83] 综述了分式规划，其中包括线性分式问题以及凸-凹分式问题等扩展（见习题 4.7）。例 4.7 中增长经济模型见 von Neumann [vN46]。

二次规划研究始于 20 世纪 50 年代（例如见 Frank 和 Wolfe [FW56]、Markowitz [Mar56]、Hildreth [Hil57]），部分动机来自第 155 页讨论的投资组合优化问题（Markowitz [Mar52]），以及第 154 页讨论的随机成本 LP（见 Freund [Fre56]）。

对二阶锥规划的兴趣较新，始于 Nesterov 和 Nemirovski [NN94, §6.2.3]。SOCP 的理论和应用综述见 Alizadeh 和 Goldfarb [AG03]，Ben-Tal 和 Nemirovski [BTN01, lecture 3]（其中该问题称为锥二次规划），以及 Lobo、Vandenberghe、Boyd 和 Lebret [LVBL98]。

鲁棒线性规划以及一般鲁棒凸优化起源于 Ben-Tal 和 Nemirovski [BTN98, BTN99] 以及 El Ghaoui 和 Lebret [EL97]。Goldfarb 和 Iyengar [GI03a, GI03b] 讨论鲁棒 QCQP 及其在投资组合优化中的应用。El Ghaoui、Oustry 和 Lebret [EOL98] 关注鲁棒半定规划。

几何规划自 20 世纪 60 年代以来已为人所知。Duffin、Peterson 和 Zener [DPZ67] 以及 Zener [Zen71] 最早倡导其在工程设计中的使用。Peterson [Pet76] 和 Ecker [Eck80] 描述了 20 世纪 70 年代取得的进展。这些文章和书籍也包含工程应用例子，特别是化学工程和土木工程中的应用。Fishburn 和 Dunlop [FD85]，Sapatnekar、Rao、Vaidya 和 Kang [SRVK93]，以及 Hershenson、Boyd 和 Lee [HBL01] 将几何规划应用于集成电路设计问题。悬臂梁设计例子（第 163 页）来自 Vanderplaats [Van84, page 147]。Perron-Frobenius 特征值的变分刻画（第 165 页）在 Berman 和 Plemmons [BP94, page 31] 中证明。

Nesterov 和 Nemirovski [NN94, chapter 4] 将锥形式问题 (4.49) 引入为非线性凸优化中的标准问题格式。Ben-Tal 和 Nemirovski [BTN01] 进一步发展了锥规划方法，并描述了大量应用。

Alizadeh [Ali91] 以及 Nesterov 和 Nemirovski [NN94, §6.4] 最早系统研究半定规划，并指出其在凸优化中的广泛应用。20 世纪 90 年代半定规划的后续研究由组合优化（Goemans 和 Williamson [GW95]）、控制（Boyd、El Ghaoui、Feron 和 Balakrishnan [BEFB94]，Scherer、Gahinet 和 Chilali [SGC97]，Dullerud 和 Paganini [DP00]）、通信和信号处理（Luo [Luo03]，Davidson、Luo、Wong 和 Ma [DLW00, MDW+02]）以及其他工程领域的应用推动。Wolkowicz、Saigal 和 Vandenberghe 编辑的书 [WSV00]，以及 Todd [Tod01]、Lewis 和 Overton [LO96]、Vandenberghe 和 Boyd [VB95] 的文章提供了综述和大量参考文献。SDP 与矩问题之间的联系，本章第 170 页给出了一个简单例子，Bertsimas 和 Sethuraman [BS00]、Nesterov [Nes00] 以及 Lasserre [Las02] 对其作了详细研究。最快混合 Markov 链问题来自 Boyd、Diaconis 和 Xiao [BDX04]。

多准则优化和 Pareto 最优性是经济学中的基本工具；见 Pareto [Par71]、Debreu [Deb59] 和 Luenberger [Lue95]。例 4.9 中的结果称为 Gauss-Markov 定理（Kailath、Sayed 和 Hassibi [KSH00, page 97]）。
