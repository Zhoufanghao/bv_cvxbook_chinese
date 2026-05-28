# 第 3 章 Convex functions

## 3.1 基本性质和例子

### 3.1.1 定义

如果 $\operatorname{dom} f$ 是凸集，并且对所有 $x,y\in\operatorname{dom} f$ 以及 $0\le\theta\le 1$，都有

$$
f(\theta x+(1-\theta)y)\le \theta f(x)+(1-\theta)f(y), \tag{3.1}
$$

则函数 $f:\mathbf{R}^n\to\mathbf{R}$ 称为凸函数。从几何上看，这个不等式表示 $(x,f(x))$ 与 $(y,f(y))$ 之间的线段，即从 $x$ 到 $y$ 的弦，位于 $f$ 的图像上方（图 3.1）。如果只要 $x\ne y$ 且 $0<\theta<1$，(3.1) 中的不等式都严格成立，则称 $f$ 为严格凸函数。如果 $-f$ 是凸函数，则称 $f$ 为凹函数；如果 $-f$ 是严格凸函数，则称 $f$ 为严格凹函数。

对于仿射函数，(3.1) 中总是取等号，因此所有仿射函数（因而所有线性函数）既是凸的也是凹的。反过来，任何既凸又凹的函数都是仿射函数。

一个函数是凸函数，当且仅当它在任何与其定义域相交的直线上的限制都是凸函数。换句话说，$f$ 是凸函数，当且仅当对所有 $x\in\operatorname{dom} f$ 和所有 $v$，函数 $g(t)=f(x+tv)$ 在其定义域 $\{t\mid x+tv\in\operatorname{dom} f\}$ 上是凸函数。这个性质很有用，因为它允许我们通过把函数限制到直线来检查凸性。

**图 3.1** 凸函数的图像。图像上任意两点之间的弦（即线段）位于图像上方。

凸函数分析是一个成熟领域，本书不会深入展开。例如，一个简单结果是：凸函数在其定义域的相对内部连续；它只能在相对边界上有不连续点。

### 3.1.2 扩展值延拓

通常可以方便地把凸函数扩展到整个 $\mathbf{R}^n$，方法是在定义域外把函数值定义为 $\infty$。如果 $f$ 是凸函数，定义其扩展值延拓 $\tilde f:\mathbf{R}^n\to\mathbf{R}\cup\{\infty\}$ 为

$$
\tilde f(x)=
\begin{cases}
f(x), & x\in\operatorname{dom} f,\\
\infty, & x\notin\operatorname{dom} f.
\end{cases}
$$

扩展函数 $\tilde f$ 定义在整个 $\mathbf{R}^n$ 上，并在 $\mathbf{R}\cup\{\infty\}$ 中取值。我们可以由 $\tilde f$ 恢复原函数的定义域：

$$
\operatorname{dom} f=\{x\mid \tilde f(x)<\infty\}.
$$

这种扩展可以简化记号，因为不必每次都显式写出定义域，或者在每次使用 $f(x)$ 时补充“对所有 $x\in\operatorname{dom} f$”。例如，基本定义不等式 (3.1) 可写成：对任意 $x,y$ 和 $0<\theta<1$，

$$
\tilde f(\theta x+(1-\theta)y)\le
\theta\tilde f(x)+(1-\theta)\tilde f(y).
$$

当 $x,y\in\operatorname{dom} f$ 时，这与 (3.1) 相同；如果其中一个点不在 $\operatorname{dom} f$ 中，则右侧为 $\infty$，不等式自动成立。这里的不等式按扩展实数的运算和排序解释。对 $\theta=0$ 或 $\theta=1$，该不等式恒成立。

另一个例子是两个凸函数 $f_1$ 和 $f_2$ 的和。通常，$f=f_1+f_2$ 的定义域为

$$
\operatorname{dom} f=\operatorname{dom} f_1\cap\operatorname{dom} f_2,
$$

且对 $x\in\operatorname{dom} f$ 有 $f(x)=f_1(x)+f_2(x)$。使用扩展值延拓后，可以简单写作

$$
\tilde f(x)=\tilde f_1(x)+\tilde f_2(x),
$$

因为只要 $x\notin\operatorname{dom} f_1$ 或 $x\notin\operatorname{dom} f_2$，右侧就为 $\infty$，定义域自动变成两者的交集。本书中只要没有歧义，将用同一符号表示凸函数及其扩展值延拓；这等价于假设所有凸函数都已在定义域外取值 $\infty$。类似地，凹函数可通过在定义域外取值 $-\infty$ 来扩展。

**例 3.1 凸集的指示函数** 设 $C\subseteq\mathbf{R}^n$ 为凸集，并考虑定义域为 $C$ 的函数 $I_C$，其中对所有 $x\in C$，$I_C(x)=0$。它的扩展值延拓为

$$
I_C(x)=
\begin{cases}
0, & x\in C,\\
\infty, & x\notin C.
\end{cases}
$$

凸函数 $I_C$ 称为集合 $C$ 的指示函数。利用指示函数可以做一些记号上的简化。例如，在集合 $C$ 上最小化函数 $f$，等价于在整个 $\mathbf{R}^n$ 上最小化 $f+I_C$；后者自动把可行点限制在 $C$ 内。

### 3.1.3 一阶条件

假设 $f$ 可微，即梯度 $\nabla f$ 在 $\operatorname{dom} f$ 的每一点都存在，并且 $\operatorname{dom} f$ 是开集。那么 $f$ 是凸函数，当且仅当 $\operatorname{dom} f$ 是凸集，并且对所有 $x,y\in\operatorname{dom} f$ 有

$$
f(y) \ge f(x)+\nabla f(x)^T(y-x) \tag{3.2}
$$

成立。这个不等式如**图 3.2** 所示。关于 $y$ 的仿射函数 $f(x)+\nabla f(x)^T(y-x)$ 是 $f$ 在 $x$ 附近的一阶 Taylor 近似。不等式 (3.2) 表明，对凸函数来说，一阶 Taylor 近似实际上是函数的全局下估计。反过来，如果函数的一阶 Taylor 近似总是其全局下估计，那么该函数是凸的。

不等式 (3.2) 表明，从凸函数的局部信息，即某一点处的函数值和导数，可以推出全局信息，即一个全局下估计。这也许是凸函数最重要的性质，并解释了凸函数和凸优化问题的一些显著性质。例如，如果 $\nabla f(x)=0$，则对所有 $y\in\operatorname{dom} f$ 有 $f(y)\ge f(x)$，也就是说 $x$ 是 $f$ 的全局极小点。

严格凸性也可以用一阶条件刻画：$f$ 是严格凸函数，当且仅当 $\operatorname{dom} f$ 是凸集，并且对 $x,y\in\operatorname{dom} f$、$x\ne y$ 有

$$
f(y)>f(x)+\nabla f(x)^T(y-x). \tag{3.3}
$$

对凹函数有相应刻画：$f$ 是凹函数，当且仅当 $\operatorname{dom} f$ 是凸集，并且对所有 $x,y\in\operatorname{dom} f$ 有

$$
f(y)\le f(x)+\nabla f(x)^T(y-x).
$$

**一阶凸性条件的证明** 为证明 (3.2)，先考虑 $n=1$ 的情形。我们证明，可微函数 $f:\mathbf{R}\to\mathbf{R}$ 是凸的，当且仅当对 $\operatorname{dom} f$ 中所有 $x$ 和 $y$ 有

$$
f(y) \ge f(x)+f'(x)(y-x). \tag{3.4}
$$

首先假设 $f$ 是凸的，且 $x,y\in\operatorname{dom} f$。由于 $\operatorname{dom} f$ 是凸集（即区间），对所有 $0<t\le1$，有 $x+t(y-x)\in\operatorname{dom} f$，并且由凸性得

$$
f(x+t(y-x))\le (1-t)f(x)+tf(y).
$$

两边整理并令 $t\to0$，得到 (3.4)。反过来，假设 (3.4) 对定义域中所有 $x,y$ 成立。取任意 $x\ne y$ 和 $0\le\theta\le1$，令 $z=\theta x+(1-\theta)y$。两次应用 (3.4) 得

$$
f(x)\ge f(z)+f'(z)(x-z),\qquad
f(y)\ge f(z)+f'(z)(y-z).
$$

将第一个不等式乘以 $\theta$，第二个乘以 $1-\theta$，再相加，可得 $\theta f(x)+(1-\theta)f(y)\ge f(z)$，这证明 $f$ 是凸的。

一般情形可通过限制到直线来证明。给定 $x,y\in\mathbf{R}^n$，定义 $g(t)=f(ty+(1-t)x)$，则

$$
g'(t)=\nabla f(ty+(1-t)x)^T(y-x).
$$

如果 $f$ 是凸的，则 $g$ 是凸的，由一维情形得 $g(1)\ge g(0)+g'(0)$，即 (3.2)。反过来，如果 (3.2) 对任意 $x,y$ 成立，则任意直线限制 $g$ 都满足一维条件 (3.4)，因此 $g$ 是凸的。由于凸性可由任意直线限制刻画，$f$ 是凸的。

### 3.1.4 二阶条件

现在假设 $f$ 二次可微，即 Hessian 矩阵 $\nabla^2 f$ 在 $\operatorname{dom} f$ 的每一点都存在，并且 $\operatorname{dom} f$ 是开集。那么 $f$ 是凸函数，当且仅当 $\operatorname{dom} f$ 是凸集，并且对所有 $x\in\operatorname{dom} f$，

$$
\nabla^2 f(x)\succeq0.
$$

对 $\mathbf{R}$ 上的函数，这化为简单条件 $f''(x)\ge0$（并且定义域为区间），也就是说导数非减。条件 $\nabla^2 f(x)\succeq0$ 可从几何上理解为函数图像在 $x$ 处具有非负曲率。二阶条件的证明留作练习（练习 3.8）。类似地，$f$ 是凹函数，当且仅当 $\operatorname{dom} f$ 是凸集，并且对所有 $x\in\operatorname{dom} f$ 有 $\nabla^2 f(x)\preceq0$。

二阶条件也给出严格凸性的一个充分条件：如果对所有 $x\in\operatorname{dom} f$ 有 $\nabla^2 f(x)\succ0$，则 $f$ 严格凸。但反过来不成立，例如 $f(x)=x^4$ 在 $\mathbf{R}$ 上严格凸，而在 $x=0$ 处二阶导数为零。

**例 3.2 二次函数** 考虑定义在 $\mathbf{R}^n$ 上的二次函数

$$
f(x)=\frac{1}{2}x^TPx+q^Tx+r,
$$

其中 $P\in\mathbf{S}^n$、$q\in\mathbf{R}^n$、$r\in\mathbf{R}$。由于对所有 $x$ 都有 $\nabla^2 f(x)=P$，所以 $f$ 凸当且仅当 $P\succeq0$，凹当且仅当 $P\preceq0$。对二次函数，严格凸性也容易刻画：$f$ 严格凸当且仅当 $P\succ0$；严格凹当且仅当 $P\prec0$。

**备注 3.1** 在一阶或二阶凸性刻画中，不能去掉“$\operatorname{dom} f$ 为凸集”这一要求。例如，函数 $f(x)=1/x^2$，定义域为 $\{x\in\mathbf{R}\mid x\ne0\}$，满足 $f''(x)>0$，但它不是凸函数。

### 3.1.5 例子

我们已经提到，所有线性函数和仿射函数都是凸函数（也是凹函数），并且已经刻画了凸二次函数和凹二次函数。本节给出更多凸函数和凹函数的例子。

先看定义在 $\mathbf{R}$ 或其子集上的一元函数。

- **指数函数。** 对任意 $a\in\mathbf{R}$，$e^{ax}$ 在 $\mathbf{R}$ 上是凸函数。
- **幂函数。** 当 $a\ge1$ 或 $a\le0$ 时，$x^a$ 在 $\mathbf{R}_{++}$ 上是凸函数；当 $0\le a\le1$ 时，$x^a$ 在 $\mathbf{R}_{++}$ 上是凹函数。
- **绝对值的幂。** 当 $p\ge1$ 时，$|x|^p$ 在 $\mathbf{R}$ 上是凸函数。
- **对数函数。** $\log x$ 在 $\mathbf{R}_{++}$ 上是凹函数。
- **负熵。** $x\log x$ 在 $\mathbf{R}_{++}$ 上是凸函数；若定义 $0\log0=0$，它在 $\mathbf{R}_+$ 上也是凸函数。

这些例子的凸性或凹性可直接验证定义不等式 (3.1)，也可检查二阶导数。例如，对 $f(x)=x\log x$，有

$$
f'(x)=\log x+1,\qquad f''(x)=1/x,
$$

所以当 $x>0$ 时 $f''(x)>0$，负熵函数严格凸。

下面是一些定义在 $\mathbf{R}^n$ 或矩阵空间上的例子。

- **范数。** $\mathbf{R}^n$ 上的任意范数都是凸函数。
- **最大函数。** $f(x)=\max\{x_1,\ldots,x_n\}$ 在 $\mathbf{R}^n$ 上是凸函数。
- **二次线性函数。** 函数 $f(x,y)=x^2/y$ 在定义域 $\mathbf{R}\times\mathbf{R}_{++}$ 上是凸函数（图 3.3）。
- **Log-sum-exp。** 函数

$$
f(x)=\log\left(e^{x_1}+\cdots+e^{x_n}\right)
$$

在 $\mathbf{R}^n$ 上是凸函数。它可看作 $\max$ 函数的光滑近似，因为

$$
\max\{x_1,\ldots,x_n\}\le f(x)
\le \max\{x_1,\ldots,x_n\}+\log n.
$$

当 $x$ 的所有分量都相等时，右侧不等式取等号。图 3.4 显示了 $n=2$ 时的情形。

- **几何平均。** 几何平均函数

$$
f(x)=\left(\prod_{i=1}^n x_i\right)^{1/n}
$$

在 $\mathbf{R}_{++}^n$ 上是凹函数。

- **对数行列式。** 函数 $f(X)=\log\det X$ 在 $\mathbf{S}_{++}^n$ 上是凹函数。

**图 3.3** 函数 $f(x,y)=x^2/y$ 的图像。

**图 3.4** 函数 $f(x,y)=\log(e^x+e^y)$ 的图像。

这些例子的凸性或凹性可以用多种方法验证，例如直接检查 (3.1)，验证 Hessian 半正定或半负定，或者把函数限制到任意直线后验证所得一元函数的凸性。

**范数。** 若 $f$ 是 $\mathbf{R}^n$ 上的范数，且 $0\le\theta\le1$，则

$$
f(\theta x+(1-\theta)y)
\le f(\theta x)+f((1-\theta)y)
=\theta f(x)+(1-\theta)f(y),
$$

其中不等式来自三角不等式，等式来自范数的齐次性。

**最大函数。** 对 $f(x)=\max_i x_i$，有

$$
\begin{aligned}
f(\theta x+(1-\theta)y)
&=\max_i\bigl(\theta x_i+(1-\theta)y_i\bigr)\\
&\le \theta\max_i x_i+(1-\theta)\max_i y_i\\
&=\theta f(x)+(1-\theta)f(y).
\end{aligned}
$$

**二次线性函数。** 对 $f(x,y)=x^2/y$，当 $y>0$ 时，

$$
\nabla^2 f(x,y)
=\frac{2}{y^3}
\begin{bmatrix}
y^2 & -xy\\
-xy & x^2
\end{bmatrix}
=\frac{2}{y^3}
\begin{bmatrix}y\\-x\end{bmatrix}
\begin{bmatrix}y\\-x\end{bmatrix}^T
\succeq0.
$$

所以该函数是凸的。

**Log-sum-exp。** 设 $z=(e^{x_1},\ldots,e^{x_n})$。Log-sum-exp 函数的 Hessian 为

$$
\nabla^2 f(x)=\frac{1}{(\mathbf{1}^Tz)^2}
\left((\mathbf{1}^Tz)\operatorname{diag}(z)-zz^T\right).
$$

要证明 $\nabla^2 f(x)\succeq0$，只需证明对任意 $v$ 有

$$
\begin{aligned}
v^T\nabla^2 f(x)v
=\frac{1}{(\mathbf{1}^Tz)^2}
\left(\left(\sum_{i=1}^n z_i\right)
\left(\sum_{i=1}^n v_i^2z_i\right)
-\left(\sum_{i=1}^n v_iz_i\right)^2\right)
\ge0.
\end{aligned}
$$

这正是 Cauchy-Schwarz 不等式应用于 $a_i=v_i\sqrt{z_i}$、$b_i=\sqrt{z_i}$ 的结果。

**几何平均。** 类似地，可证明几何平均 $f(x)=(\prod_{i=1}^n x_i)^{1/n}$ 在 $\mathbf{R}_{++}^n$ 上是凹的。其 Hessian 可写为

$$
\nabla^2 f(x)
=-\frac{\left(\prod_{i=1}^n x_i\right)^{1/n}}{n^2}
\left(n\operatorname{diag}(1/x_1^2,\ldots,1/x_n^2)-qq^T\right),
\qquad q_i=1/x_i.
$$

因此对任意 $v$，

$$
v^T\nabla^2 f(x)v
=-\frac{\left(\prod_{i=1}^n x_i\right)^{1/n}}{n^2}
\left(n\sum_{i=1}^n\frac{v_i^2}{x_i^2}
-\left(\sum_{i=1}^n\frac{v_i}{x_i}\right)^2\right)
\le0,
$$

其中最后一步同样来自 Cauchy-Schwarz 不等式。

**对数行列式。** 对 $f(X)=\log\det X$，可通过限制到任意直线来验证凹性。令 $X=Z+tV$，其中 $Z,V\in\mathbf{S}^n$，并考虑 $Z+tV\succ0$ 的 $t$ 区间。若 $Z\succ0$，则

$$
\begin{aligned}
g(t)&=\log\det(Z+tV)\\
&=\log\det\left(Z^{1/2}(I+tZ^{-1/2}VZ^{-1/2})Z^{1/2}\right)\\
&=\sum_{i=1}^n\log(1+t\lambda_i)+\log\det Z,
\end{aligned}
$$

其中 $\lambda_i$ 是 $Z^{-1/2}VZ^{-1/2}$ 的特征值。于是

$$
g'(t)=\sum_{i=1}^n\frac{\lambda_i}{1+t\lambda_i},
\qquad
 g''(t)=-\sum_{i=1}^n\frac{\lambda_i^2}{(1+t\lambda_i)^2}\le0.
$$

所以 $g$ 是凹的，因而 $\log\det X$ 是凹函数。

### 3.1.6 下水平集

函数 $f:\mathbf{R}^n\to\mathbf{R}$ 的 $\alpha$-下水平集定义为

$$
C_\alpha=\{x\in\operatorname{dom} f\mid f(x)\le\alpha\}.
$$

凸函数的任意下水平集都是凸集。证明很直接：若 $x,y\in C_\alpha$，则 $f(x)\le\alpha$ 且 $f(y)\le\alpha$，所以对 $0\le\theta\le1$，

$$
f(\theta x+(1-\theta)y)\le \theta f(x)+(1-\theta)f(y)\le\alpha,
$$

从而 $\theta x+(1-\theta)y\in C_\alpha$。

反过来并不成立：一个函数可以有所有下水平集都是凸集，但自身不是凸函数。例如 $f(x)=-e^x$ 在 $\mathbf{R}$ 上不是凸函数（它严格凹），但它的所有下水平集都是凸的。若 $f$ 是凹函数，则其上水平集 $\{x\in\operatorname{dom} f\mid f(x)\ge\alpha\}$ 是凸集。

下水平集性质常用于证明集合凸性：把集合表示为凸函数的下水平集，或凹函数的上水平集。

**例 3.3** 对 $x\in\mathbf{R}_+^n$，几何平均和算术平均分别为

$$
G(x)=\left(\prod_{i=1}^n x_i\right)^{1/n},
\qquad
A(x)=\frac{1}{n}\sum_{i=1}^n x_i,
$$

其中在 $G$ 中约定 $0^{1/n}=0$。算术-几何平均不等式说明 $G(x)\le A(x)$。给定 $0\le\alpha\le1$，集合

$$
\{x\in\mathbf{R}_+^n\mid G(x)\ge\alpha A(x)\}
$$

是凸集，因为它是凹函数 $G(x)-\alpha A(x)$ 的 $0$-上水平集。事实上该集合还是正齐次的，因此是凸锥。

### 3.1.7 上图

函数 $f:\mathbf{R}^n\to\mathbf{R}$ 的图像是集合

$$
\{(x,f(x))\mid x\in\operatorname{dom} f\}\subseteq\mathbf{R}^{n+1}.
$$

它的上图（epigraph）定义为

$$
\operatorname{epi} f=
\{(x,t)\mid x\in\operatorname{dom} f,\ f(x)\le t\}
\subseteq\mathbf{R}^{n+1}.
$$

“epi”意为“在上方”，所以上图就是图像上方的集合，如图 3.5 所示。凸集和凸函数之间的联系可由上图精确刻画：函数是凸函数，当且仅当其上图是凸集。类似地，函数是凹函数，当且仅当其下图

$$
\operatorname{hypo} f=\{(x,t)\mid t\le f(x)\}
$$

是凸集。

**图 3.5** 函数 $f$ 的上图。阴影区域表示 $\operatorname{epi} f$，其下边界是 $f$ 的图像。

**例 3.4 矩阵分式函数** 函数

$$
f(x,Y)=x^T Y^{-1}x,
\qquad
\operatorname{dom} f=\mathbf{R}^n\times\mathbf{S}_{++}^n,
$$

是凸函数。它推广了二次线性函数 $x^2/y$。一种证明方法是考察上图：

$$
\begin{aligned}
\operatorname{epi} f
&=\{(x,Y,t)\mid Y\succ0,
\ x^TY^{-1}x\le t\}\\
&=\left\{(x,Y,t)\ \middle|\
\begin{bmatrix}Y & x\\ x^T & t\end{bmatrix}\succeq0,
\ Y\succ0\right\},
\end{aligned}
$$

其中第二个等价式来自 Schur 补条件。最后一个条件是 $(x,Y,t)$ 中的线性矩阵不等式，因此上图是凸集。

当 $n=1$ 时，它化为二次线性函数 $x^2/y$，其 LMI 表示为

$$
\begin{bmatrix}y & x\\ x & t\end{bmatrix}\succeq0,
\qquad y>0.
$$

凸函数的许多性质都可以用上图作几何解释。例如一阶条件 (3.2) 可解释为：对任意 $(y,t)\in\operatorname{epi} f$，

$$
t\ge f(y)\ge f(x)+\nabla f(x)^T(y-x),
$$

等价于

$$
\begin{bmatrix}\nabla f(x)\\ -1\end{bmatrix}^T
\left(
\begin{bmatrix}y\\ t\end{bmatrix}
-
\begin{bmatrix}x\\ f(x)\end{bmatrix}
\right)\le0.
$$

因此向量 $(\nabla f(x),-1)$ 在边界点 $(x,f(x))$ 处定义了 $\operatorname{epi} f$ 的支持超平面（图 3.6）。

**图 3.6** 对可微凸函数 $f$，向量 $(\nabla f(x),-1)$ 定义了 $x$ 处上图的支持超平面。

### 3.1.8 Jensen 不等式

基本凸性不等式 (3.1)

$$
f(\theta x+(1-\theta)y)\le \theta f(x)+(1-\theta)f(y)
$$

也称为 Jensen 不等式。它可推广到有限多个点的凸组合：若 $f$ 是凸函数，$x_1,\ldots,x_k\in\operatorname{dom} f$，$\theta_i\ge0$ 且 $\sum_{i=1}^k\theta_i=1$，则

$$
f\left(\sum_{i=1}^k\theta_i x_i\right)
\le \sum_{i=1}^k\theta_i f(x_i).
$$

同样，它也可推广到积分和期望。例如，若 $p(x)\ge0$ 且 $\int_S p(x)dx=1$，则在积分存在时

$$
f\left(\int_S p(x)x\,dx\right)
\le \int_S f(x)p(x)\,dx.
$$

更一般地，如果随机变量 $x$ 以概率 $1$ 落在 $\operatorname{dom} f$ 中，则

$$
f(\mathbf{E}x) \le \mathbf{E}f(x). \tag{3.5}
$$

基本不等式 (3.1) 可由 (3.5) 取两点分布得到。反过来，(3.5) 也刻画凸性：若 $f$ 不是凸函数，则存在随机变量 $x$ 使 $f(\mathbf{E}x)>\mathbf{E}f(x)$。

**备注 3.2** 若 $z$ 是零均值随机向量，且 $x+z\in\operatorname{dom} f$ 几乎处处成立，则

$$
\mathbf{E}f(x+z)\ge f(x).
$$

因此，对凸函数而言，给参数加入零均值随机扰动不会降低函数值的平均值。

### 3.1.9 不等式

许多著名不等式都可由 Jensen 不等式应用于合适的凸函数导出。作为简单例子，考虑算术-几何平均不等式：对 $a,b\ge0$，

$$
\sqrt{ab}\le (a+b)/2. \tag{3.6}
$$

函数 $-\log x$ 是凸函数。用 $\theta=1/2$ 应用 Jensen 不等式得

$$
-\log\left(\frac{a+b}{2}\right)
\le -\frac{\log a+\log b}{2}.
$$

两边取指数即可得 (3.6)。

作为稍复杂的例子，可证明 Hölder 不等式：若 $p>1$，$1/p+1/q=1$，则对 $x,y\in\mathbf{R}^n$，

$$
\sum_{i=1}^n |x_i y_i|
\le
\left(\sum_{i=1}^n |x_i|^p\right)^{1/p}
\left(\sum_{i=1}^n |y_i|^q\right)^{1/q}.
$$

由 $-\log x$ 的凸性可得加权算术-几何平均不等式

$$
a^\theta b^{1-\theta}\le \theta a+(1-\theta)b,
\qquad a,b\ge0,
\quad 0\le\theta\le1.
$$

令

$$
a=\frac{|x_i|^p}{\sum_j |x_j|^p},
\qquad
b=\frac{|y_i|^q}{\sum_j |y_j|^q},
\qquad
\theta=1/p,
$$

并对 $i$ 求和，即可推出 Hölder 不等式。

## 3.2 保持凸性的运算

在本节中，我们将描述一些保留函数凸性或凹性的操作，或者允许我们构造新的凸函数和凹函数。我们从一些简单的操作开始，例如加法、缩放和逐点求上，然后描述一些更复杂的操作（其中一些包括作为特殊情况的简单操作）。

### 3.2.1 非负加权和

显然，如果 $f$ 是凸函数且 $\alpha\ge0$，则函数 $\alpha f$ 是凸函数。如果 $f_1$ 和 $f_2$ 都是凸函数，那么它们的和 $f_1+f_2$ 也是凸函数。结合非负缩放和加法可知，凸函数集本身构成一个凸锥：凸函数的非负加权和

$$
f=w_1f_1+\cdots+w_mf_m
$$

是凸的。类似地，凹函数的非负加权和是凹的。严格凸（或严格凹）函数的非负、非零加权和是严格凸（或严格凹）的。

这些性质可扩展到无限和与积分。例如，如果对每个 $y\in A$，$f(x,y)$ 关于 $x$ 是凸的，并且对每个 $y\in A$，$w(y)\ge0$，则由

$$
g(x)=\int_A w(y)f(x,y)\,dy
$$

定义的函数 $g$ 关于 $x$ 是凸的（假设积分存在）。非负缩放和加法保持凸性这一事实很容易直接验证，也可以从相应的上图看出。例如，如果 $w\ge0$ 且 $f$ 是凸函数，则

$$
\operatorname{epi}(wf)=
\begin{bmatrix}
I & 0\\
0 & w
\end{bmatrix}
\operatorname{epi}f,
$$

它是凸的，因为凸集在线性映射下的像仍是凸集。

### 3.2.2 与仿射映射复合

假设 $f:\mathbf{R}^n\to\mathbf{R}$，$A\in\mathbf{R}^{n\times m}$，$b\in\mathbf{R}^n$。定义 $g:\mathbf{R}^m\to\mathbf{R}$ 为

$$
g(x)=f(Ax+b),
\qquad
\operatorname{dom}g=\{x\mid Ax+b\in\operatorname{dom}f\}.
$$

那么如果 $f$ 是凸的，$g$ 也是凸的；如果 $f$ 是凹的，$g$ 也是凹的。

### 3.2.3 逐点最大值和上确界

如果 $f_1$ 和 $f_2$ 是凸函数，则逐点最大值

$$
f(x)=\max\{f_1(x),f_2(x)\},
\qquad
\operatorname{dom} f=\operatorname{dom} f_1\cap\operatorname{dom} f_2
$$

也是凸函数。事实上，若 $0\le\theta\le1$ 且 $x,y\in\operatorname{dom} f$，则

$$
\begin{aligned}
f(\theta x+(1-\theta)y)
&=\max\{f_1(\theta x+(1-\theta)y),f_2(\theta x+(1-\theta)y)\}\\
&\le\max\{\theta f_1(x)+(1-\theta)f_1(y),
\theta f_2(x)+(1-\theta)f_2(y)\}\\
&\le \theta\max\{f_1(x),f_2(x)\}
+(1-\theta)\max\{f_1(y),f_2(y)\}\\
&=\theta f(x)+(1-\theta)f(y).
\end{aligned}
$$

同理，有限个凸函数的逐点最大值也是凸函数。

**例 3.5 分段线性函数** 函数

$$
f(x)=\max\{a_1^Tx+b_1,\ldots,a_L^Tx+b_L\}
$$

定义了一个分段线性（或分段仿射）凸函数。它是凸的，因为它是仿射函数的逐点最大值。反过来，任何具有 $L$ 个或更少区域的分段线性凸函数都可写成这种形式（见练习 3.29）。

**例 3.6 最大的 $r$ 个分量之和** 对 $x\in\mathbf{R}^n$，记 $x_{[i]}$ 为 $x$ 的第 $i$ 大分量，即

$$
x_{[1]}\ge x_{[2]}\ge\cdots\ge x_{[n]}.
$$

函数

$$
f(x)=\sum_{i=1}^r x_{[i]}
$$

是凸函数，因为它可写成

$$
f(x)=\max\{x_{i_1}+\cdots+x_{i_r}
\mid 1\le i_1<i_2<\cdots<i_r\le n\},
$$

即若干线性函数的逐点最大值。进一步地，若 $w_1\ge w_2\ge\cdots\ge w_r\ge0$，则 $\sum_{i=1}^r w_i x_{[i]}$ 也是凸函数（见练习 3.19）。

逐点最大值的性质可推广到任意凸函数族的逐点上确界。若对每个 $y\in A$，函数 $f(x,y)$ 关于 $x$ 是凸的，则

$$
g(x)=\sup_{y\in A} f(x,y) \tag{3.7}
$$

关于 $x$ 是凸函数。这里

$$
\operatorname{dom} g=
\left\{x\mid (x,y)\in\operatorname{dom} f\ \text{对所有 }y\in A,
\ \sup_{y\in A} f(x,y)<\infty\right\}.
$$

类似地，凹函数族的逐点下确界是凹函数。从上图角度看，逐点上确界对应上图的交集：若 $g$ 如 (3.7) 定义，则

$$
\operatorname{epi} g=\bigcap_{y\in A}\operatorname{epi} f(\cdot,y),
$$

因此 $g$ 的上图是凸集族的交集。

**例 3.7 集合的支持函数** 令 $C\subseteq\mathbf{R}^n$ 且 $C\ne\emptyset$。集合 $C$ 的支持函数定义为

$$
S_C(x)=\sup\{x^Ty\mid y\in C\},
$$

其定义域为 $\operatorname{dom} S_C=\{x\mid \sup_{y\in C}x^Ty<\infty\}$。对每个 $y\in C$，$x^Ty$ 是 $x$ 的线性函数，因此 $S_C$ 是线性函数族的逐点上确界，因而是凸函数。

**例 3.8 到集合最远点的距离** 设 $C\subseteq\mathbf{R}^n$。在任意范数下，到集合 $C$ 最远点的距离

$$
f(x)=\sup_{y\in C}\|x-y\|
$$

是凸函数，因为对任意固定 $y$，$\|x-y\|$ 关于 $x$ 是凸函数，而 $f$ 是这些凸函数的逐点上确界。

**例 3.9 最小二乘成本作为权重的函数** 令 $a_1,\ldots,a_n\in\mathbf{R}^m$。在加权最小二乘问题中，考虑最小化

$$
\sum_{i=1}^n w_i(a_i^Tx-b_i)^2
$$

其中变量为 $x\in\mathbf{R}^m$，权重 $w_i$ 可为负。定义最优加权最小二乘成本

$$
g(w)=\inf_x\sum_{i=1}^n w_i(a_i^Tx-b_i)^2.
$$

由于 $g$ 是关于 $w$ 的一族线性函数的下确界（由 $x$ 参数化），所以 $g$ 是凹函数。若 $W=\operatorname{diag}(w)$，矩阵 $A$ 的第 $i$ 行为 $a_i^T$，则

$$
g(w)=\inf_x (Ax-b)^TW(Ax-b).
$$

当 $A^TWA\not\succeq0$ 时，该二次函数对 $x$ 下无界，所以 $g(w)=-\infty$。当 $A^TWA\succ0$ 时，解析最小化得

$$
g(w)=b^TWb-b^TWA(A^TWA)^{-1}A^TWb.
$$

这个表达式中凹性并不明显，但从“下确界保持凹性”的结论可立即得到。

**例 3.10 对称矩阵的最大特征值** 函数 $f(X)=\lambda_{\max}(X)$（定义域为 $\mathbf{S}^m$）是凸函数，因为

$$
\lambda_{\max}(X)=\sup\{y^TXy\mid \|y\|_2=1\},
$$

它是关于 $X$ 的一族线性函数的逐点上确界。

**例 3.11 矩阵范数** 谱范数 $f(X)=\|X\|_2$ 是凸函数，因为

$$
\|X\|_2=\sup\{u^TXv\mid \|u\|_2=1,
\|v\|_2=1\}.
$$

更一般地，若 $\|\cdot\|_a$ 和 $\|\cdot\|_b$ 分别是 $\mathbf{R}^p$ 和 $\mathbf{R}^q$ 上的范数，则诱导范数

$$
\|X\|_{a,b}=\sup_{v\ne0}\frac{\|Xv\|_a}{\|v\|_b}
$$

可写为

$$
\|X\|_{a,b}=\sup\{u^TXv\mid \|u\|_{a,*}=1,
\|v\|_b=1\},
$$

其中 $\|\cdot\|_{a,*}$ 是 $\|\cdot\|_a$ 的对偶范数。由于它是关于 $X$ 的线性函数族的逐点上确界，所以是凸函数。

**表示为仿射函数的逐点上确界** 上述例子说明了一种建立凸性的常用方法：把函数表示为仿射函数族的逐点上确界。除技术条件外，逆命题也成立：几乎所有凸函数都可表示为仿射函数族的逐点上确界。例如，若 $f:\mathbf{R}^n\to\mathbf{R}$ 是凸函数且 $\operatorname{dom} f=\mathbf{R}^n$，则

$$
f(x)=\sup\{g(x)\mid g\ \text{仿射，且对所有 }z,
\ g(z)\le f(z)\}.
$$

显然，右侧不超过 $f(x)$。为证明等号，只需证明对每个 $x$，存在一个仿射全局下估计 $g$ 满足 $g(x)=f(x)$。由于 $\operatorname{epi} f$ 是凸集，在边界点 $(x,f(x))$ 处存在支持超平面，即存在 $(a,b)\ne0$，使得对所有 $(z,t)\in\operatorname{epi} f$，

$$
\begin{bmatrix}a\\ b\end{bmatrix}^T
\begin{bmatrix}x-z\\ f(x)-t\end{bmatrix}\le0.
$$

这意味着对所有 $z\in\operatorname{dom} f=\mathbf{R}^n$ 和所有 $s\ge0$，

$$
a^T(x-z)+b(f(x)-f(z)-s)\le0. \tag{3.8}
$$

为了使 (3.8) 对所有 $s\ge0$ 成立，必须有 $b\ge0$。如果 $b=0$，则对所有 $z$ 有 $a^T(x-z)\le0$，这推出 $a=0$，与 $(a,b)\ne0$ 矛盾。因此 $b>0$。令 $s=0$ 并整理可得

$$
g(z)=f(x)+(a/b)^T(x-z)\le f(z),
$$

且 $g(x)=f(x)$。这说明 $g$ 是所需的仿射全局下估计。

### 3.2.4 复合

本节讨论函数 $h:\mathbf{R}^k\to\mathbf{R}$ 和 $g:\mathbf{R}^n\to\mathbf{R}^k$ 在什么条件下能保证复合函数 $f=h\circ g$ 的凸性或凹性，其中

$$
f(x)=h(g(x)),\qquad
\operatorname{dom} f=\{x\in\operatorname{dom} g\mid g(x)\in\operatorname{dom} h\}.
$$

**标量复合** 首先考虑 $k=1$ 的情形，即 $h:\mathbf{R}\to\mathbf{R}$、$g:\mathbf{R}^n\to\mathbf{R}$。由于凸性可由函数在任意直线上的限制刻画，我们可以先考虑 $n=1$。为发现复合规则，假设 $h$ 和 $g$ 二次可微，且 $\operatorname{dom} g=\operatorname{dom} h=\mathbf{R}$。此时 $f$ 的凸性化为 $f''\ge0$。复合函数 $f=h\circ g$ 的二阶导数为

$$
f''(x)=h''(g(x))g'(x)^2+h'(g(x))g''(x). \tag{3.9}
$$

例如，如果 $g$ 是凸的（所以 $g''\ge0$），$h$ 是凸且非减的（所以 $h''\ge0$ 且 $h'\ge0$），则由 (3.9) 得 $f''\ge0$，即 $f$ 是凸的。类似地，由 (3.9) 可得以下规则：

$$
\begin{array}{ll}
h\ \text{凸且非减},\ g\ \text{凸} &\Longrightarrow h\circ g\ \text{凸},\\
h\ \text{凸且非增},\ g\ \text{凹} &\Longrightarrow h\circ g\ \text{凸},\\
h\ \text{凹且非减},\ g\ \text{凹} &\Longrightarrow h\circ g\ \text{凹},\\
h\ \text{凹且非增},\ g\ \text{凸} &\Longrightarrow h\circ g\ \text{凹}.
\end{array} \tag{3.10}
$$

当 $g$ 和 $h$ 二次可微且定义域都是整个实线时，这些陈述成立。在一般情形中，不要求 $h$ 和 $g$ 可微，也不要求 $\operatorname{dom} g=\mathbf{R}^n$ 或 $\operatorname{dom} h=\mathbf{R}$，类似规则仍然成立：

$$
\begin{array}{ll}
h\ \text{凸},\ \tilde h\ \text{非减},\ g\ \text{凸} &\Longrightarrow h\circ g\ \text{凸},\\
h\ \text{凸},\ \tilde h\ \text{非增},\ g\ \text{凹} &\Longrightarrow h\circ g\ \text{凸},\\
h\ \text{凹},\ \tilde h\ \text{非减},\ g\ \text{凹} &\Longrightarrow h\circ g\ \text{凹},\\
h\ \text{凹},\ \tilde h\ \text{非增},\ g\ \text{凸} &\Longrightarrow h\circ g\ \text{凹}.
\end{array} \tag{3.11}
$$

这里 $\tilde h$ 表示 $h$ 的扩展值延拓：当 $h$ 凸（凹）时，它在 $\operatorname{dom} h$ 外取值 $\infty$（$-\infty$）。(3.11) 与 (3.10) 的差别在于，单调性要求施加在扩展值函数 $\tilde h$ 上。

为了理解这一点，设 $h$ 是凸函数，所以 $\tilde h$ 在 $\operatorname{dom} h$ 外取值 $\infty$。说 $\tilde h$ 非减，意味着对任意 $x<y$ 都有 $\tilde h(x)\le\tilde h(y)$。特别地，如果 $y\in\operatorname{dom} h$，则必须有 $x\in\operatorname{dom} h$。换句话说，$h$ 的定义域必须向负方向无限延伸，可以是 $\mathbf{R}$、$(-\infty,a)$ 或 $(-\infty,a]$。类似地，如果 $h$ 是凸的且 $\tilde h$ 非增，则 $h$ 非增，并且 $\operatorname{dom} h$ 向正方向无限延伸。见**图 3.7**。

例 3.12 的几个简单函数说明了复合定理中对 $h$ 的条件。函数 $h(x)=\log x$，$\operatorname{dom} h=\mathbf{R}_{++}$，是凹函数，并且 $\tilde h$ 非减。函数 $h(x)=x^{1/2}$，$\operatorname{dom} h=\mathbf{R}_{+}$，是凹函数，并且 $\tilde h$ 非减。函数 $h(x)=x^{3/2}$，$\operatorname{dom} h=\mathbf{R}_{+}$，是凸函数，但不满足 $\tilde h$ 非减；若改为 $h(x)=\max\{x,0\}^{3/2}$ 且 $\operatorname{dom} h=\mathbf{R}$，则该函数是凸的，并且满足 $\tilde h$ 非减。

**图 3.7** 左图函数 $x^2$ 在定义域 $\mathbf{R}_{+}$ 上是凸的，但其扩展值延拓不是非递减的。右图函数 $\max\{x,0\}^2$ 定义在 $\mathbf{R}$ 上，是凸且非递减的。

可以不借助可微性直接证明复合结果 (3.11)。作为例子，我们证明：若 $g$ 是凸的，$h$ 是凸的，并且 $\tilde h$ 非减，则 $f=h\circ g$ 是凸的。设 $x,y\in\operatorname{dom} f$，且 $0\le\theta\le1$。由于 $x,y\in\operatorname{dom} f$，有 $x,y\in\operatorname{dom} g$ 且 $g(x),g(y)\in\operatorname{dom} h$。


由于 $\operatorname{dom} g$ 是凸的，可得 $\theta x+(1-\theta)y\in\operatorname{dom} g$，并且

$$
g(\theta x+(1-\theta)y) \le \theta g(x)+(1-\theta)g(y) \tag{3.12}
$$


由于 $g(x),g(y)\in\operatorname{dom} h$，可得 $\theta g(x)+(1-\theta)g(y)\in\operatorname{dom} h$，即 (3.12) 的右侧在 $\operatorname{dom} h$ 中。又因为 $\tilde h$ 非减，其定义域向负方向无限延伸，所以 (3.12) 的左侧 $g(\theta x+(1-\theta)y)$ 也在 $\operatorname{dom} h$ 中。至此已经证明 $\operatorname{dom} f$ 是凸的。


由 $\tilde h$ 的非减性和不等式 (3.12)，得到

$$
h(g(\theta x+(1-\theta)y)) \le h(\theta g(x)+(1-\theta)g(y)) \tag{3.13}
$$


根据 $h$ 的凸性，我们有

$$
h(\theta g(x)+(1-\theta)g(y))
\le \theta h(g(x))+(1-\theta)h(g(y)). \tag{3.14}
$$

由 (3.13) 和 (3.14) 得

$$
h(g(\theta x+(1-\theta)y))
\le \theta h(g(x))+(1-\theta)h(g(y)),
$$

因此 $h\circ g$ 是凸的。

例 3.13 给出几个简单复合结果。如果 $g$ 是凸的，则 $\exp g(x)$ 是凸的。如果 $g$ 是凹且为正的，则 $1/g(x)$ 是凸的。如果 $g$ 是凸且非负的，并且 $p\ge1$，则 $g(x)^p$ 是凸的。如果 $g$ 是凸的，则 $-\log(-g(x))$ 在 $\{x\mid g(x)<0\}$ 上是凸的。

备注 3.3 不能只要求函数 $h$ 单调，而忽略其扩展值延拓 $\tilde h$ 的单调性。例如令 $g(x)=x^2$，$\operatorname{dom} g=\mathbf{R}$，并令 $h(x)=0$，$\operatorname{dom} h=[1,2]$。这里 $g$ 是凸的，$h$ 是凸且非减的，但 $f=h\circ g$ 的定义域为 $[-\sqrt2,-1]\cup[1,\sqrt2]$，不是凸集，因此 $f$ 不是凸函数。

**向量复合** 现在转向 $k\ge1$ 的情形。设

$$
f(x)=h(g(x))=h(g_1(x),\ldots,g_k(x)),
$$

其中 $h:\mathbf{R}^k\to\mathbf{R}$，$g_i:\mathbf{R}^n\to\mathbf{R}$。同样，不失一般性，可先考虑 $n=1$。若函数二次可微，$\operatorname{dom} g=\mathbf{R}$ 且 $\operatorname{dom} h=\mathbf{R}^k$，则

$$
f''(x)=g'(x)^T\nabla^2 h(g(x))g'(x)+\nabla h(g(x))^Tg''(x), \tag{3.15}
$$

这是 (3.9) 的向量形式。一般情形中，仍需要把 $h$ 的单调性条件理解为扩展值延拓 $\tilde h$ 的单调性。例如，当 $h:\mathbf{R}^k\to\mathbf{R}$ 是凸函数且 $\tilde h$ 非递减时，只要 $u\preceq v$ 就有 $\tilde h(u)\le\tilde h(v)$。这意味着如果 $v\in\operatorname{dom} h$，则 $u\in\operatorname{dom} h$；也就是说，$h$ 的定义域必须沿 $-\mathbf{R}^k_+$ 方向无限延伸，可简洁写成 $\operatorname{dom} h-\mathbf{R}^k_+=\operatorname{dom} h$。

例 3.14 给出一些向量复合示例。若 $g_1,\ldots,g_k$ 是 $\mathbf{R}^n$ 上的凸函数，则 $r$ 个最大 $g_i$ 的逐点和是凸函数。函数 $h(z)=\log(\sum_{i=1}^k e^{z_i})$ 是凸函数且在每个参数中非减，因此当 $g_i$ 凸时，$\log(\sum_{i=1}^k e^{g_i})$ 是凸函数。当 $0<p\le1$ 时，$\mathbf{R}^k_+$ 上的函数 $h(z)=(\sum_{i=1}^k z_i^p)^{1/p}$ 是凹函数，并且其扩展（在 $z\not\succeq0$ 时取值 $-\infty$）在每个分量中非减。若 $p\ge1$，$g_1,\ldots,g_k$ 为非负凸函数，则 $(\sum_{i=1}^k g_i(x)^p)^{1/p}$ 是凸函数。$\mathbf{R}^k_+$ 上的几何平均值 $h(z)=(\prod_{i=1}^k z_i)^{1/k}$ 是凹函数，并且其扩展在每个参数中非减；因此若 $g_1,\ldots,g_k$ 是非负凹函数，则它们的几何平均值 $(\prod_{i=1}^k g_i)^{1/k}$ 也是非负凹函数。

### 3.2.5 最小化

我们已经看到，任意凸函数族的最大值或上确界仍是凸函数。某些特殊形式的最小化也会产生凸函数。

若 $f(x,y)$ 关于 $(x,y)$ 联合凸，且 $C$ 是凸集，则函数

$$
g(x)=\inf_{y\in C} f(x,y) \tag{3.16}
$$

关于 $x$ 是凸函数，前提是对所有 $x$ 有 $g(x)>-\infty$。$g$ 的定义域是 $\operatorname{dom} f$ 在 $x$ 坐标上的投影：

$$
\operatorname{dom} g=\{x\mid \text{存在 }y\in C\text{ 使 }(x,y)\in\operatorname{dom} f\}.
$$

证明可直接验证 Jensen 不等式。取 $x_1,x_2\in\operatorname{dom} g$ 和 $\epsilon>0$，存在 $y_1,y_2\in C$ 使得

$$
f(x_i,y_i)\le g(x_i)+\epsilon,
\qquad i=1,2.
$$

对 $0\le\theta\le1$，由于 $C$ 凸，$\theta y_1+(1-\theta)y_2\in C$，于是

$$
\begin{aligned}
g(\theta x_1+(1-\theta)x_2)
&\le f(\theta x_1+(1-\theta)x_2,
\theta y_1+(1-\theta)y_2)\\
&\le \theta f(x_1,y_1)+(1-\theta)f(x_2,y_2)\\
&\le \theta g(x_1)+(1-\theta)g(x_2)+\epsilon.
\end{aligned}
$$

令 $\epsilon\to0$，得到 $g$ 的凸性。

这个结果也可从上图解释。如果下确界总能在某个 $y\in C$ 处达到，则

$$
\operatorname{epi} g
=\{(x,t)\mid \text{存在 }y\in C\text{ 使 }(x,y,t)\in\operatorname{epi} f\},
$$

即 $\operatorname{epi} f$ 在部分坐标上的投影，所以是凸集。

**例 3.15 Schur 补** 设二次函数

$$
f(x,y)=x^TAx+2x^TBy+y^TCy
$$

关于 $(x,y)$ 凸，即

$$
\begin{bmatrix}A & B\\ B^T & C\end{bmatrix}\succeq0.
$$

则

$$
g(x)=\inf_y f(x,y)=x^T(A-BC^\dagger B^T)x,
$$

其中 $C^\dagger$ 是 $C$ 的伪逆。由最小化规则，$g$ 是凸函数，因此 $A-BC^\dagger B^T\succeq0$。若 $C\succ0$，矩阵 $A-BC^{-1}B^T$ 称为块矩阵中 $C$ 的 Schur 补。

**例 3.16 到集合的距离** 点 $x$ 到集合 $S\subseteq\mathbf{R}^n$ 的距离定义为

$$
\operatorname{dist}(x,S)=\inf_{y\in S}\|x-y\|.
$$

函数 $\|x-y\|$ 关于 $(x,y)$ 是凸的。因此若 $S$ 是凸集，距离函数 $\operatorname{dist}(x,S)$ 关于 $x$ 是凸函数。

**例 3.17** 若 $h$ 是凸函数，则

$$
g(x)=\inf\{h(y)\mid Ay=x\}
$$

是凸函数。可令

$$
f(x,y)=
\begin{cases}
h(y), & Ay=x,\\
\infty, & \text{否则},
\end{cases}
$$

则 $f$ 关于 $(x,y)$ 凸，而 $g$ 是对 $y$ 的最小化。

### 3.2.6 函数的透视

函数 $f:\mathbf{R}^n\to\mathbf{R}$ 的透视函数定义为

$$
g(x,t)=t f(x/t),
$$

其定义域为

$$
\operatorname{dom} g=
\{(x,t)\mid x/t\in\operatorname{dom} f,
\ t>0\}.
$$

透视运算保持凸性：若 $f$ 是凸函数，则 $g$ 是凸函数；若 $f$ 是凹函数，则 $g$ 是凹函数。可直接验证定义不等式，也可用上图和透视映射证明。对 $t>0$，

$$
\begin{aligned}
(x,t,s)\in\operatorname{epi} g
&\Longleftrightarrow t f(x/t)\le s\\
&\Longleftrightarrow f(x/t)\le s/t\\
&\Longleftrightarrow (x/t,s/t)\in\operatorname{epi} f.
\end{aligned}
$$

因此 $\operatorname{epi} g$ 是 $\operatorname{epi} f$ 在透视映射下的逆像，所以是凸集。

**例 3.18 欧几里得范数平方** 凸函数 $f(x)=x^Tx$ 的透视为

$$
g(x,t)=t(x/t)^T(x/t)=\frac{x^Tx}{t},
\qquad t>0.
$$

因此 $x^Tx/t$ 关于 $(x,t)$ 是凸函数。它也可看作二次线性函数 $x_i^2/t$ 的和，或矩阵分式函数的特例。

**例 3.19 负对数** 对 $\mathbf{R}_{++}$ 上的凸函数 $f(x)=-\log x$，透视函数为

$$
g(x,t)=-t\log(x/t)=t\log(t/x)=t\log t-t\log x,
$$

在 $\mathbf{R}_{++}^2$ 上凸。该函数称为 $t$ 和 $x$ 的相对熵；当 $x=1$ 时化为负熵。

两个正向量 $u,v\in\mathbf{R}_{++}^n$ 的相对熵

$$
\sum_{i=1}^n u_i\log(u_i/v_i)
$$

关于 $(u,v)$ 是凸函数。相关的 Kullback-Leibler 散度定义为

$$
D_{\mathrm{kl}}(u,v)=\sum_{i=1}^n
\left(u_i\log(u_i/v_i)-u_i+v_i\right). \tag{3.17}
$$

它也是凸函数，因为它等于相对熵加上线性项。Kullback-Leibler 散度满足 $D_{\mathrm{kl}}(u,v)\ge0$，且当且仅当 $u=v$ 时取零，可作为两个正向量差异的度量。当 $u$ 和 $v$ 都是概率向量时，相对熵和 Kullback-Leibler 散度相同。

若在相对熵中取 $v_i=\mathbf{1}^Tu$，可得到 $u\in\mathbf{R}_{++}^n$ 的归一化熵函数

$$
\sum_{i=1}^n u_i\log\frac{\mathbf{1}^Tu}{u_i}
=(\mathbf{1}^Tu)\sum_{i=1}^n z_i\log(1/z_i),
\qquad z=\frac{u}{\mathbf{1}^Tu}.
$$

它是凹且齐次的。

**例 3.20** 设 $f:\mathbf{R}^m\to\mathbf{R}$ 是凸函数，$A\in\mathbf{R}^{m\times n}$、$b\in\mathbf{R}^m$、$c\in\mathbf{R}^n$、$d\in\mathbf{R}$。定义

$$
g(x)=(c^Tx+d)f\left(\frac{Ax+b}{c^Tx+d}\right),
$$

其中

$$
\operatorname{dom} g=\left\{x\mid c^Tx+d>0,
\ \frac{Ax+b}{c^Tx+d}\in\operatorname{dom} f\right\}.
$$

则 $g$ 是凸函数。

## 3.3 共轭函数

本节介绍一个后面会经常用到的操作。给定函数 $f:\mathbf{R}^n\to\mathbf{R}$，其共轭函数 $f^*:\mathbf{R}^n\to\mathbf{R}$ 定义为

$$
f^*(y)=\sup_{x\in\operatorname{dom} f}(y^Tx-f(x)). \tag{3.18}
$$

共轭函数的定义域由使上确界有限的 $y$ 组成。无论 $f$ 是否凸，$f^*$ 总是凸函数，因为它是关于 $y$ 的一族仿射函数的逐点上确界。若 $f$ 已按扩展值方式定义，则上式中对 $x\in\operatorname{dom} f$ 的限制可省略。

图 3.8 给出了一维情形的几何解释：$f^*(y)$ 是直线函数 $yx$ 与 $f(x)$ 之间的最大差距；若 $f$ 可微，最大值出现在满足 $f'(x)=y$ 的点。

**例 3.21 一元函数的共轭**

- **仿射函数。** 若 $f(x)=ax+b$，则 $yx-ax-b$ 关于 $x$ 有界当且仅当 $y=a$。因此 $\operatorname{dom} f^*=\{a\}$，且 $f^*(a)=-b$。
- **负对数。** 若 $f(x)=-\log x$，$\operatorname{dom} f=\mathbf{R}_{++}$，则

$$
\operatorname{dom} f^*=-\mathbf{R}_{++},
\qquad
f^*(y)=-\log(-y)-1.
$$

- **指数函数。** 若 $f(x)=e^x$，则

$$
\operatorname{dom} f^*=\mathbf{R}_+,
\qquad
f^*(y)=y\log y-y,
$$

其中约定 $0\log0=0$。

- **负熵。** 若 $f(x)=x\log x$，$\operatorname{dom} f=\mathbf{R}_+$，且 $f(0)=0$，则

$$
\operatorname{dom} f^*=\mathbf{R},
\qquad
f^*(y)=e^{y-1}.
$$

- **倒数。** 若 $f(x)=1/x$，$\operatorname{dom} f=\mathbf{R}_{++}$，则

$$
\operatorname{dom} f^*=-\mathbf{R}_+,
\qquad
f^*(y)=-2(-y)^{1/2}.
$$

**例 3.22 严格凸二次函数** 对

$$
f(x)=\frac{1}{2}x^TQx,
\qquad Q\in\mathbf{S}_{++}^n,
$$

函数 $y^Tx-(1/2)x^TQx$ 在 $x=Q^{-1}y$ 处达到最大值，因此

$$
f^*(y)=\frac{1}{2}y^TQ^{-1}y.
$$

**例 3.23 对数行列式** 对 $\mathbf{S}_{++}^n$ 上的函数

$$
f(X)=-\log\det X,
$$

其共轭为

$$
f^*(Y)=\sup_{X\succ0}\left(\operatorname{tr}(YX)+\log\det X\right).
$$

若 $Y\not\prec0$，该上确界无界；若 $Y\prec0$，一阶条件 $Y+X^{-1}=0$ 给出 $X=-Y^{-1}$，从而

$$
\operatorname{dom} f^*=-\mathbf{S}_{++}^n,
\qquad
f^*(Y)=-\log\det(-Y)-n.
$$

**例 3.24 指示函数** 设 $I_S$ 是集合 $S\subseteq\mathbf{R}^n$ 的指示函数，即在 $S$ 上取 $0$，在 $S$ 外取 $\infty$。则

$$
I_S^*(y)=\sup_{x\in S} y^Tx,
$$

即集合 $S$ 的支持函数。

**例 3.25 Log-sum-exp 函数** 对

$$
f(x)=\log\left(\sum_{i=1}^n e^{x_i}\right),
$$

其共轭为概率单纯形上的负熵：

$$
f^*(y)=
\begin{cases}
\sum_{i=1}^n y_i\log y_i,
& y\succeq0,
\ \mathbf{1}^Ty=1,\\
\infty, & \text{否则}.
\end{cases}
$$

这里约定 $0\log0=0$。

**例 3.26 范数** 设 $\|\cdot\|$ 是 $\mathbf{R}^n$ 上的范数，对偶范数为 $\|\cdot\|_*$。函数 $f(x)=\|x\|$ 的共轭为

$$
f^*(y)=
\begin{cases}
0, & \|y\|_*\le1,\\
\infty, & \text{否则}.
\end{cases}
$$

也就是说，范数的共轭是对偶范数单位球的指示函数。

**例 3.27 范数平方** 对 $f(x)=(1/2)\|x\|^2$，其共轭为

$$
f^*(y)=\frac{1}{2}\|y\|_*^2.
$$

这个结论来自不等式 $y^Tx\le\|y\|_*\|x\|$ 以及一维二次函数的最大值。

**例 3.28 收入和利润函数** 设 $r$ 表示资源消耗向量，$S(r)$ 表示销售收入，$p$ 表示资源价格。最大利润为

$$
M(p)=\sup_r(S(r)-p^Tr).
$$

因此 $M(p)=(-S)^*(-p)$。也就是说，最大利润函数与收入函数的共轭密切相关。

### 3.3.2 基本性质

**Fenchel 不等式** 由共轭函数定义立即得到，对所有 $x,y$，

$$
f(x)+f^*(y)\ge x^Ty.
$$

这称为 Fenchel 不等式；当 $f$ 可微时也称为 Young 不等式。例如，若 $f(x)=(1/2)x^TQx$ 且 $Q\succ0$，则

$$
x^Ty\le \frac{1}{2}x^TQx+\frac{1}{2}y^TQ^{-1}y.
$$

**共轭的共轭** 如果 $f$ 是闭凸函数，则

$$
f^{**}=f.
$$

例如，当 $\operatorname{dom} f=\mathbf{R}^n$ 且 $f$ 是闭凸函数时，$f$ 是其所有仿射全局下估计的逐点上确界。

**可微函数** 若 $f$ 是凸且可微的，并且 $\operatorname{dom} f=\mathbf{R}^n$，则 $y^Tx-f(x)$ 的最大化点 $x^*$ 满足

$$
y=\nabla f(x^*).
$$

反过来，若 $x^*$ 满足该方程，则它最大化 $y^Tx-f(x)$，并且

$$
f^*(y)=(x^*)^T\nabla f(x^*)-f(x^*).
$$

**缩放和平移** 若 $a>0$、$b\in\mathbf{R}$，且

$$
g(x)=af(x)+b,
$$

则

$$
g^*(y)=a f^*(y/a)-b.
$$

若 $A$ 非奇异、$b\in\mathbf{R}^n$，且 $g(x)=f(Ax+b)$，则

$$
g^*(y)=f^*(A^{-T}y)-b^TA^{-T}y,
\qquad
\operatorname{dom} g^*=A^T\operatorname{dom} f^*.
$$

**独立函数之和** 若

$$
f(u,v)=f_1(u)+f_2(v),
$$

则

$$
f^*(w,z)=f_1^*(w)+f_2^*(z).
$$

## 3.4 拟凸函数

### 3.4.1 定义和例子

函数 $f:\mathbf{R}^n\to\mathbf{R}$ 称为拟凸函数（quasiconvex），如果其定义域以及所有下水平集

$$
S_\alpha=\{x\in\operatorname{dom} f\mid f(x)\le\alpha\}
$$

都是凸集。若 $-f$ 是拟凸函数，则称 $f$ 为拟凹函数；这等价于每个上水平集 $\{x\mid f(x)\ge\alpha\}$ 都是凸集。既拟凸又拟凹的函数称为拟线性函数。拟线性函数的定义域和每个水平集 $\{x\mid f(x)=\alpha\}$ 都是凸集。

对 $\mathbf{R}$ 上的函数，拟凸性要求所有下水平集都是区间（可以是无限区间）。凸函数一定是拟凸函数，但反过来不成立。

**图 3.9** $\mathbf{R}$ 上的拟凸函数。每个 $\alpha$-下水平集 $S_\alpha$ 都是区间。

**例 3.29** $\mathbf{R}_{++}$ 上的对数函数是拟凸函数，也是拟凹函数，因此是拟线性函数。

**例 3.30 向量长度** 定义 $x\in\mathbf{R}^n$ 的长度为非零分量的最大下标：

$$
f(x)=\max\{i\mid x_i\ne0\},
$$

并定义零向量长度为 $0$。该函数是拟凸的，因为其下水平集是子空间：

$$
f(x)\le\alpha
\Longleftrightarrow
x_i=0,
\qquad i=\lfloor\alpha\rfloor+1,\ldots,n.
$$

**例 3.31 乘积函数** 函数 $f(x_1,x_2)=x_1x_2$ 在 $\mathbf{R}_+^2$ 上是拟凹函数，因为集合

$$
\{x\in\mathbf{R}_+^2\mid x_1x_2\ge\alpha\}
$$

对所有 $\alpha$ 都是凸集。该函数既不是凸函数也不是凹函数，因为 Hessian

$$
\nabla^2 f(x)=\begin{bmatrix}0&1\\1&0\end{bmatrix}
$$

是不定的。注意它在整个 $\mathbf{R}^2$ 上不是拟凹的。

**例 3.32 线性分式函数** 函数

$$
f(x)=\frac{a^Tx+b}{c^Tx+d},
\qquad
\operatorname{dom} f=\{x\mid c^Tx+d>0\},
$$

是拟凸且拟凹的，即拟线性的。其 $\alpha$-下水平集为

$$
\begin{aligned}
S_\alpha
&=\left\{x\mid c^Tx+d>0,
\frac{a^Tx+b}{c^Tx+d}\le\alpha\right\}\\
&=\{x\mid c^Tx+d>0,
\ a^Tx+b\le\alpha(c^Tx+d)\},
\end{aligned}
$$

是开半空间和闭半空间的交集，因此是凸集。上水平集同理。

**例 3.33 距离比函数** 给定 $a,b\in\mathbf{R}^n$，定义

$$
f(x)=\frac{\|x-a\|_2}{\|x-b\|_2}.
$$

在半空间 $\{x\mid \|x-a\|_2\le\|x-b\|_2\}$ 上，$f$ 是拟凸函数。对 $\alpha\le1$，下水平集满足

$$
\|x-a\|_2\le\alpha\|x-b\|_2.
$$

平方并整理得

$$
(1-\alpha^2)x^Tx-2(a-\alpha^2b)^Tx+a^Ta-\alpha^2b^Tb\le0,
$$

这是一个欧几里得球（或半空间），因此是凸集。

**例 3.34 内部收益率** 令 $x=(x_0,x_1,\ldots,x_n)$ 表示现金流，$x_i>0$ 表示第 $i$ 期收到现金，$x_i<0$ 表示支付现金。利率 $r\ge0$ 下的现值为

$$
\operatorname{PV}(x,r)=\sum_{i=0}^n (1+r)^{-i}x_i.
$$

考虑满足 $x_0<0$ 且 $x_0+x_1+\cdots+x_n>0$ 的现金流。内部收益率定义为使现值为零的最小非负利率：

$$
\operatorname{IRR}(x)=\inf\{r\ge0\mid \operatorname{PV}(x,r)=0\}.
$$

它是 $x$ 的拟凹函数。因为

$$
\operatorname{IRR}(x)\ge R
\Longleftrightarrow
\operatorname{PV}(x,r)>0,
\qquad 0\le r<R.
$$

右侧是由 $r$ 索引的一族开半空间的交集，因此是凸集。

### 3.4.2 基本性质

拟凸性是凸性的显著推广，但许多凸函数的性质仍有类似物。拟凸性的 Jensen 型刻画为：$f$ 是拟凸函数，当且仅当 $\operatorname{dom} f$ 是凸集，并且

$$
f(\theta x+(1-\theta)y) \le \max\{f(x),f(y)\} \tag{3.19}
$$

对所有 $x,y\in\operatorname{dom} f$ 和 $0\le\theta\le1$ 成立。换句话说，函数在线段上的值不超过两端点处较大的值。

**例 3.35 非负向量的基数** 向量 $x\in\mathbf{R}^n$ 的基数 $\operatorname{card}(x)$ 是非零分量的个数。函数 $\operatorname{card}$ 在 $\mathbf{R}_+^n$ 上是拟凹的，因为对 $x,y\succeq0$，

$$
\operatorname{card}(x+y)\ge
\min\{\operatorname{card}(x),\operatorname{card}(y)\}.
$$

**例 3.36 半正定矩阵的秩** 函数 $\operatorname{rank}X$ 在 $\mathbf{S}_+^n$ 上是拟凹的，因为对 $X,Y\in\mathbf{S}_+^n$，

$$
\operatorname{rank}(X+Y)\ge
\min\{\operatorname{rank}X,\operatorname{rank}Y\}.
$$

与凸性一样，拟凸性也可由任意直线限制来刻画。特别地，可通过将函数限制到任意直线并检查所得一元函数的拟凸性来验证拟凸性。

对连续的一元函数 $f:\mathbf{R}\to\mathbf{R}$，拟凸性等价于以下三种情形之一：$f$ 非减；$f$ 非增；或者存在点 $c$，使得 $f$ 在 $t\le c$ 时非增，在 $t\ge c$ 时非减。点 $c$ 可取为任一全局极小点。

### 3.4.3 可微拟凸函数

**一阶条件** 假设 $f:\mathbf{R}^n\to\mathbf{R}$ 可微。那么 $f$ 是拟凸的，当且仅当 $\operatorname{dom} f$ 是凸集，并且对所有 $x,y\in\operatorname{dom} f$ 有

$$
f(y)\le f(x)\Longrightarrow \nabla f(x)^T(y-x)\le 0. \tag{3.20}
$$

几何上，当 $\nabla f(x)\ne0$ 时，$\nabla f(x)$ 在点 $x$ 处定义了子水平集 $\{y\mid f(y)\le f(x)\}$ 的支持超平面。

**二阶条件** 假设 $f$ 二次可微。如果 $f$ 是拟凸函数，则对所有 $x\in\operatorname{dom} f$ 和 $y\in\mathbf{R}^n$，有

$$
y^T\nabla f(x)=0 \Longrightarrow y^T\nabla^2 f(x)y\ge 0. \tag{3.21}
$$

对一元函数，这化为 $f'(x)=0\Longrightarrow f''(x)\ge0$。当 $\nabla f(x)\ne0$ 时，(3.21) 表示 $\nabla^2 f(x)$ 在子空间 $\nabla f(x)^\perp$ 上半正定，因此 Hessian 最多有一个负特征值。

一个部分逆命题是：若对所有 $x\in\operatorname{dom} f$ 和所有非零 $y\in\mathbf{R}^n$，

$$
y^T\nabla f(x)=0 \Longrightarrow y^T\nabla^2 f(x)y>0, \tag{3.22}
$$

则 $f$ 是拟凸函数。证明可限制到任意直线，从一元情形推出。

### 3.4.4 保持拟凸性的运算

**非负加权最大值** 若 $f_i$ 拟凸且 $w_i\ge0$，则

$$
f(x)=\max\{w_1f_1(x),\ldots,w_mf_m(x)\}
$$

是拟凸函数。更一般地，拟凸函数族的逐点上确界仍是拟凸函数，因为其下水平集是一族凸下水平集的交集。

**例 3.37 广义特征值** 对 $Y\succ0$，最大广义特征值

$$
\lambda_{\max}(X,Y)=\sup_{u\ne0}\frac{u^TXu}{u^TYu}
$$

在 $\mathbf{S}^n\times\mathbf{S}_{++}^n$ 上是拟凸函数。对固定 $u\ne0$，商 $u^TXu/u^TYu$ 关于 $(X,Y)$ 是线性分式函数，因而拟凸；上确界保持拟凸性。

**复合** 若 $g$ 拟凸，$h:\mathbf{R}\to\mathbf{R}$ 非减，则 $h\circ g$ 拟凸。拟凸函数与仿射映射或线性分式变换复合仍得到拟凸函数。例如，若 $f$ 拟凸，则 $g(x)=f(Ax+b)$ 拟凸；在 $c^Tx+d>0$ 的定义域上，

$$
\tilde g(x)=f\left(\frac{Ax+b}{c^Tx+d}\right)
$$

也是拟凸函数。

**最小化** 若 $f(x,y)$ 关于 $(x,y)$ 联合拟凸，且 $C$ 是凸集，则

$$
g(x)=\inf_{y\in C} f(x,y)
$$

是拟凸函数。证明思路是：$g$ 的 $\alpha$-下水平集可由 $f$ 的近似 $\alpha$-下水平集投影得到，并保持凸性。

### 3.4.5 通过凸函数族表示

把拟凸函数的下水平集表示为凸函数的不等式通常很方便。我们寻找一族由 $t\in\mathbf{R}$ 索引的凸函数 $\phi_t:\mathbf{R}^n\to\mathbf{R}$，使得

$$
f(x)\le t \Longleftrightarrow \phi_t(x)\le0. \tag{3.23}
$$

也就是说，$f$ 的 $t$-下水平集是凸函数 $\phi_t$ 的 $0$-下水平集。若对每个 $x$，$\phi_t(x)$ 关于 $t$ 非增，则该表示与下水平集的嵌套性质一致。

这种表示总是存在：可令 $\phi_t$ 为集合 $\{x\mid f(x)\le t\}$ 的指示函数。若下水平集闭，也可取

$$
\phi_t(x)=\operatorname{dist}(x,\{z\mid f(z)\le t\}).
$$

**例 3.38 凹函数上的凸函数** 设 $p$ 是凸函数，$q$ 是凹函数，并且在凸集 $C$ 上 $p(x)\ge0$、$q(x)>0$。则

$$
f(x)=\frac{p(x)}{q(x)}
$$

在 $C$ 上拟凸。因为

$$
f(x)\le t \Longleftrightarrow p(x)-tq(x)\le0,
$$

当 $t\ge0$ 时，$\phi_t(x)=p(x)-tq(x)$ 是凸函数，且对每个 $x$ 关于 $t$ 非增。

## 3.5 对数凹函数和对数凸函数

### 3.5.1 定义

函数 $f:\mathbf{R}^n\to\mathbf{R}$ 称为对数凹函数，如果 $f(x)>0$ 对所有 $x\in\operatorname{dom} f$ 成立，并且 $\log f$ 是凹函数。若 $\log f$ 是凸函数，则称 $f$ 为对数凸函数。显然，$f$ 对数凸当且仅当 $1/f$ 对数凹。

允许 $f$ 取零值也很方便；此时取 $\log f(x)=-\infty$。若扩展值函数 $\log f$ 是凹函数，则称 $f$ 为对数凹函数。

对数凹性也可不通过对数来表述：若 $\operatorname{dom} f$ 凸且 $f(x)>0$，则 $f$ 对数凹当且仅当对所有 $x,y\in\operatorname{dom} f$ 和 $0\le\theta\le1$，

$$
f(\theta x+(1-\theta)y)
\ge f(x)^\theta f(y)^{1-\theta}.
$$

因此，对数凹函数在两点凸组合处的值至少为两端点函数值的几何平均。由复合规则可知，对数凸函数是凸函数；非负凹函数是对数凹函数。对数凸函数是拟凸的，对数凹函数是拟凹的。

**例 3.39 简单例子**

- 仿射函数 $f(x)=a^Tx+b$ 在 $\{x\mid a^Tx+b>0\}$ 上是对数凹函数。
- 幂函数 $f(x)=x^a$ 在 $\mathbf{R}_{++}$ 上，当 $a\le0$ 时对数凸，当 $a\ge0$ 时对数凹。
- 指数函数 $e^{ax}$ 既对数凸又对数凹。
- 高斯分布函数

$$
\Phi(x)=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^x e^{-u^2/2}\,du
$$

是对数凹函数（见练习 3.54）。

- Gamma 函数

$$
\Gamma(x)=\int_0^\infty u^{x-1}e^{-u}\,du
$$

在 $x\ge1$ 时对数凸（见练习 3.52）。

- $\det X$ 在 $\mathbf{S}_{++}^n$ 上对数凹。
- $\det X/\operatorname{tr}X$ 在 $\mathbf{S}_{++}^n$ 上对数凹（见练习 3.49）。

**例 3.40 对数凹密度函数** 很多常见概率密度都是对数凹的。例如，多元正态密度

$$
f(x)=\frac{1}{\sqrt{(2\pi)^n\det\Sigma}}
\exp\left(-\frac{1}{2}(x-\bar x)^T\Sigma^{-1}(x-\bar x)\right)
$$

是对数凹的，其中 $\Sigma\in\mathbf{S}_{++}^n$。指数分布

$$
f(x)=\left(\prod_{i=1}^n\lambda_i\right)e^{-\lambda^Tx},
\qquad \lambda\succ0,
$$

在 $\mathbf{R}_+^n$ 上也是对数凹的。

凸集 $C$ 上的均匀分布也是对数凹的：

$$
f(x)=
\begin{cases}
1/\alpha, & x\in C,\\
0, & x\notin C,
\end{cases}
$$

其中 $\alpha=\operatorname{vol}(C)$。这是因为 $\log f$ 在 $C$ 上为常数，在 $C$ 外为 $-\infty$。

Wishart 密度

$$
f(X)=a(\det X)^{(p-n-1)/2}
\exp\left(-\frac{1}{2}\operatorname{tr}(\Sigma^{-1}X)\right),
\qquad X\in\mathbf{S}_{++}^n,
$$

在 $p>n$ 且 $p-n-1\ge0$ 时也是对数凹的，因为 $\log\det X$ 凹，而迹项是仿射函数。

### 3.5.2 性质

**二阶条件** 假设 $f$ 二次可微且 $\operatorname{dom} f$ 凸。由于

$$
\nabla^2\log f(x)
=\frac{1}{f(x)}\nabla^2 f(x)
-\frac{1}{f(x)^2}\nabla f(x)\nabla f(x)^T,
$$

可得：$f$ 对数凸当且仅当对所有 $x\in\operatorname{dom} f$，

$$
f(x)\nabla^2 f(x)\succeq
\nabla f(x)\nabla f(x)^T;
$$

$f$ 对数凹当且仅当对所有 $x\in\operatorname{dom} f$，

$$
f(x)\nabla^2 f(x)\preceq
\nabla f(x)\nabla f(x)^T.
$$

**乘法、加法和积分** 对数凸性和对数凹性在乘法和正比例缩放下封闭。若 $f$ 和 $g$ 对数凹，则 $fg$ 对数凹，因为

$$
\log(fg)=\log f+\log g.
$$

对数凸函数的和仍对数凸。若 $f=e^F$、$g=e^G$，其中 $F,G$ 凸，则

$$
\log(f+g)=\log(e^F+e^G)
$$

由 log-sum-exp 的凸性可知是凸函数。更一般地，若对每个 $y\in C$，$f(x,y)$ 关于 $x$ 对数凸，则

$$
g(x)=\int_C f(x,y)\,dy
$$

也是对数凸函数。

**例 3.41 Laplace 变换和矩母函数** 若 $p(x)\ge0$，其 Laplace 变换

$$
P(z)=\int p(x)e^{-z^Tx}\,dx
$$

是对数凸函数。若 $p$ 是概率密度，则 $M(z)=P(-z)$ 称为矩母函数。函数 $\log M(z)$ 是凸函数，称为累积量生成函数。

**对数凹函数的积分** 一个重要但证明较深的事实是：若 $f:\mathbf{R}^n\times\mathbf{R}^m\to\mathbf{R}$ 对数凹，则

$$
g(x)=\int f(x,y)\,dy
$$

关于 $x$ 对数凹。其直接后果包括：对数凹密度的边缘分布仍对数凹；对数凹函数的卷积仍对数凹。

若 $C\subseteq\mathbf{R}^n$ 是凸集，随机向量 $w$ 具有对数凹密度 $p$，则

$$
f(x)=\operatorname{prob}(x+w\in C)
$$

关于 $x$ 对数凹。

**例 3.42 累积分布函数** 若随机向量 $w$ 的密度 $f$ 对数凹，则其累积分布函数

$$
F(x)=\operatorname{prob}(w\preceq x)
=\int_{-\infty}^{x_n}\cdots\int_{-\infty}^{x_1}
 f(z)\,dz_1\cdots dz_n
$$

也是对数凹函数。

**例 3.43 良率函数** 设 $x$ 是产品参数的目标值，制造误差为随机向量 $w$，可接受参数集合为 $S$。良率函数为

$$
Y(x)=\operatorname{prob}(x+w\in S).
$$

若 $w$ 的密度对数凹且 $S$ 是凸集，则 $Y$ 对数凹。因此任意 $\alpha$-良率区域 $\{x\mid Y(x)\ge\alpha\}$ 是凸集。

**例 3.44 多面体体积** 设 $A\in\mathbf{R}^{m\times n}$，定义

$$
P_u=\{x\in\mathbf{R}^n\mid Ax\preceq u\}.
$$

则 $\operatorname{vol}(P_u)$ 是 $u$ 的对数凹函数。证明可令

$$
\Psi(x,u)=
\begin{cases}
1, & Ax\preceq u,\\
0, & \text{否则},
\end{cases}
$$

并应用对数凹函数积分保持对数凹的结论。

## 3.6 关于广义不等式的凸性

本节用广义不等式替代通常的实数序关系，推广单调性和凸性的概念。

### 3.6.1 关于广义不等式的单调性

设 $K\subseteq\mathbf{R}^n$ 是与广义不等式 $\preceq_K$ 相关的真锥。如果

$$
x\preceq_K y \Longrightarrow f(x)\le f(y),
$$

则 $f:\mathbf{R}^n\to\mathbf{R}$ 称为 $K$-非减函数。如果 $x\preceq_K y$ 且 $x\ne y$ 蕴含 $f(x)<f(y)$，则称 $f$ 为 $K$-增函数。$K$-非增和 $K$-减函数类似定义。

**例 3.45 单调向量函数** 对 $K=\mathbf{R}_+^n$，$f$ 是 $K$-非减的，当且仅当只要 $x_i\le y_i$ 对所有 $i$ 成立，就有 $f(x)\le f(y)$。这等价于 $f$ 对每个分量变量单调非减。

**例 3.46 矩阵单调函数** 若 $f:\mathbf{S}^n\to\mathbf{R}$ 关于半正定锥单调，则称为矩阵单调函数。例如：

- $\operatorname{tr}(WX)$ 在 $W\succeq0$ 时矩阵非减，在 $W\succ0$ 时矩阵增；
- $\operatorname{tr}(X^{-1})$ 在 $\mathbf{S}_{++}^n$ 上矩阵递减；
- $\det X$ 在 $\mathbf{S}_{++}^n$ 上矩阵递增，在 $\mathbf{S}_+^n$ 上矩阵非减。

**梯度条件** 对具有凸定义域的可微函数，$f$ 是 $K$-非减函数，当且仅当对所有 $x\in\operatorname{dom} f$，

$$
\nabla f(x) \succeq_{K^*} 0. \tag{3.24}
$$

这里 $K^*$ 是对偶锥。若对所有 $x\in\operatorname{dom} f$，

$$
\nabla f(x) \succ_{K^*} 0, \tag{3.25}
$$

则 $f$ 是 $K$-增函数。证明与一元单调性的导数判别类似，只是把非负性换成对偶锥中的非负性。

### 3.6.2 关于广义不等式的凸性

设 $K\subseteq\mathbf{R}^m$ 是真锥，并定义广义不等式 $\preceq_K$。若对所有 $x,y$ 和 $0\le\theta\le1$，

$$
f(\theta x+(1-\theta)y)
\preceq_K \theta f(x)+(1-\theta)f(y),
$$

则称 $f:\mathbf{R}^n\to\mathbf{R}^m$ 是 $K$-凸函数。若对 $x\ne y$ 和 $0<\theta<1$ 上式严格成立，则称 $f$ 严格 $K$-凸。当 $m=1$ 且 $K=\mathbf{R}_+$ 时，这就是通常的凸性定义。

**例 3.47 分量不等式下的凸性** 当 $K=\mathbf{R}_+^m$ 时，$f$ 是 $K$-凸，当且仅当每个分量函数 $f_i$ 都是凸函数。

**例 3.48 矩阵凸性** 若 $f:\mathbf{R}^n\to\mathbf{S}^m$，并且对任意 $x,y$、$0\le\theta\le1$，

$$
f(\theta x+(1-\theta)y)
\preceq \theta f(x)+(1-\theta)f(y),
$$

则称 $f$ 为矩阵凸函数。等价地，对所有向量 $z$，标量函数 $z^Tf(x)z$ 都是凸函数。

典型例子包括：

- $f(X)=XX^T$ 是矩阵凸函数，因为 $z^TXX^Tz=\|X^Tz\|_2^2$ 关于 $X$ 是凸二次函数；
- $X^p$ 在 $\mathbf{S}_{++}^n$ 上，当 $1\le p\le2$ 或 $-1\le p\le0$ 时矩阵凸；当 $0\le p\le1$ 时矩阵凹；
- 当 $n\ge2$ 时，$e^X$ 不是 $\mathbf{S}^n$ 上的矩阵凸函数。

**对偶刻画** 函数 $f$ 是 $K$-凸的，当且仅当对每个 $w\succeq_{K^*}0$，实值函数 $w^Tf$ 是凸函数。$f$ 严格 $K$-凸，当且仅当对每个非零 $w\succeq_{K^*}0$，$w^Tf$ 严格凸。

**可微 $K$-凸函数** 可微函数 $f$ 是 $K$-凸的，当且仅当其定义域凸，并且对所有 $x,y\in\operatorname{dom} f$，

$$
f(y)\succeq_K f(x)+Df(x)(y-x),
$$

其中 $Df(x)\in\mathbf{R}^{m\times n}$ 是 $f$ 在 $x$ 处的导数（Jacobian）。严格 $K$-凸有对应的严格不等式版本。

**复合定理** 若 $g:\mathbf{R}^n\to\mathbf{R}^p$ 是 $K$-凸函数，$h:\mathbf{R}^p\to\mathbf{R}$ 是凸函数，并且 $h$ 的扩展值延拓 $\tilde h$ 是 $K$-非减的，则 $h\circ g$ 是凸函数。条件 $\tilde h$ 是 $K$-非减的可写为

$$
\operatorname{dom} h-K=\operatorname{dom} h.
$$

**例 3.49** 定义二次矩阵函数

$$
g(X)=X^TAX+B^TX+X^TB+C,
$$

其中 $A\in\mathbf{S}^m$、$B\in\mathbf{R}^{m\times n}$、$C\in\mathbf{S}^n$。若 $A\succeq0$，则 $g$ 是矩阵凸函数。函数

$$
h(Y)=-\log\det(-Y)
$$

在 $\operatorname{dom} h=-\mathbf{S}_{++}^n$ 上凸且递增。因此

$$
f(X)=-\log\det\left(-(X^TAX+B^TX+X^TB+C)\right)
$$

在定义域

$$
\{X\mid X^TAX+B^TX+X^TB+C\prec0\}
$$

上是凸函数。

## Bibliography

凸分析的标准参考是 Rockafellar [Roc70]。其他关于凸函数的书籍有 Stoer 和 Witzgall [ SW70]、Roberts 和 $Varberg [ R V73]$、$Van Tiel [ vT84]$、$Hiriart-Urruty$ 和 Lemar´ echal [HUL93]、Ekeland 和 T´ emam [ET99]、Borwein 和 Lewis [BL00]、Florenzano 和 Le Van [ FL01]、Barvinok [ Bar02] 和Bertsekas、Nedi´ c 和 Ozdaglar [ Ber03]。大多数非线性编程教科书还包括有关凸函数的章节（例如，请参见 Mangasarian [ Man94]、Bazaraa、Sherali 和 Shetty [ BSS93]、Bertsekas [Ber99]、Polyak [ Pol87] 以及 Peressini、Sullivan 和 Uhl [ PSU88]）。 Jensen 不等式出现在 [Jen06] 中。 Hardy、Little wood 和 P´olya [HLP52] 以及 Beckenbach 和 Bellman [BB65] 提出了对不平等的一般研究，其中 Jensen 不平等发挥了核心作用。术语透视函数来自 $Hiriart-Urruty$ 和 Lemar´ echal [ HUL93，第 1 卷，第 100 页]。对于**例 3.19（相对熵和 $Kullback-Leibler$ 散度）中的定义以及相关练习 3.13，请参阅 Cover 和 Thomas [CT91]。关于拟凸函数（以及凸性的其他扩展）的一些重要的早期参考文献包括 $Nikaid^ o [Nik54]$、$Mangasarian [ Man94$，第 9 章]、Arrow 和 Enthoven [AE61]、Ponstein [ Pon67] 和 Luenberger [ Lue68]。如需更全面的参考列表，请参阅 Bazaraa、Sherali 和 Shetty [BSS93，第 126 页]。 Pr´ ekopa [Pr´ e80] 给出了对数凹函数的调查。$Barndorff-Nielsen [BN78, §7]$ 中提到了拉普拉斯变换的对数凸性。有关对数凹函数积分结果的证明，请参阅 Pr´ ekopa [Pr´ e71, Pr´ e73]。从 Nesterov 和 Nemirovski 开始，广义不等式在最近关于锥规划的文献中被广泛使用 [NN94，第 156 页]；另请参见 $Ben-Tal$ 和 Nemirovski [BTN01] 以及第 4 章末尾的参考文献。关于广义不等式的凸性也出现在 Luenberger [ Lue69, §8.2] 和 Isii [ Isi64] 的著作中。矩阵单调性和矩阵凸性归因于 L´owner [L´ow34]，并且由 Davis [Dav63]、Roberts 和 Varberg [R V73，第 216 页] 以及 Marshall 和 Olkin [MO79，§16E] 详细讨论。例 3.48 中函数 X p 的凸凹性结果，参见 Bondar [ Bon94，定理 16.1]。有关证明 eX 不是凸矩阵的简单示例，请参阅 Marshall 和 Olkin [MO79，第 474 页]。

## Exercises


### 凸性的定义

3.1 假设 $f:\mathbf{R}\to\mathbf{R}$ 是凸的，并且 $a,b\in\operatorname{dom}f$ 且 $a<b$。 (a) 证明对所有 $x\in[a,b]$，

$$
f(x)\le \frac{b-x}{b-a}f(a)+\frac{x-a}{b-a}f(b).
$$

(b) 证明对所有 $x\in(a,b)$，

$$
\frac{f(x)-f(a)}{x-a}\le
\frac{f(b)-f(a)}{b-a}\le
\frac{f(b)-f(x)}{b-x}.
$$

画一张草图说明这个不等式。 (c) 假设 $f$ 可微。使用 (b) 中的结果证明 $f'(a)\le (f(b)-f(a))/(b-a)\le f'(b)$。注意这些不等式也可由 (3.2) 得出：$f(b)\ge f(a)+f'(a)(b-a)$，$f(a)\ge f(b)+f'(b)(a-b)$。 (d) 假设 $f$ 二次可微。使用 (c) 中的结果证明 $f''(a)\ge0$ 且 $f''(b)\ge0$。

3.2 凸函数、凹函数、拟凸函数和拟凹函数的水平集。函数 f 的一些水平集如下所示。标记为 1 的曲线显示 ${x |f (x) = 1 }$ 等。 1 2 3 f 可以是凸的（凹的、拟凸的、拟凹的）吗？解释一下你的答案。对如下所示的水平曲线重复此操作。 1 2 3 4 5 6

3.3 递增凸函数的反函数。假设 $f : R \to R$ 在其域 ( a,b ) 上递增且凸。令 g 表示其反函数，即具有定义域 ( f (a),f (b)) 且 g(f (x)) = x 的函数，其中 $a<x<b$。关于 g 的凸性或凹性你能说什么？

3.4 [RV73，第 15 页] 证明连续函数 $f:\mathbf{R}^n\to\mathbf{R}$ 是凸函数，当且仅当对于每个线段，其在线段上的平均值小于或等于其在线段端点处的平均值：对每个 $x,y\in\mathbf{R}^n$，

$$
\int_0^1 f(x+\lambda(y-x))\,d\lambda\le \frac{f(x)+f(y)}{2}.
$$

3.5 [RV73，第 22 页] 凸函数的运行平均值。假设 $f:\mathbf{R}\to\mathbf{R}$ 是凸的，且 $\mathbf{R}_+\subseteq\operatorname{dom}f$。证明其运行平均值 $F$ 是凸的，其中

$$
F(x)=\frac{1}{x}\int_0^x f(t)\,dt,\qquad \operatorname{dom}F=\mathbf{R}_{++}.
$$

提示：对每个 $s$，$f(sx)$ 关于 $x$ 是凸的，因此 $\int_0^1 f(sx)\,ds$ 是凸的。

3.6 功能和上图。函数的上图什么时候是半空间？函数的上图什么时候是凸锥体？函数的上图什么时候是多面体？

3.7 假设$f : \mathbf{R}^n \to R$ 是凸的，且 $\operatorname{dom} f = \mathbf{R}^n$，并且在 $\mathbf{R}^n$ 上有界。证明 f 是常数。

3.8 凸性的二阶条件。证明二次可微函数 f 是凸函数当且仅当其域是凸函数并且对于所有 $x \in \operatorname{dom} f\nabla^2 f (x) \succeq 0$。暗示。首先考虑情况：$R \to R$。您可以使用凸性的一阶条件（在第 70 页已证明）。

3.9 仿射集凸性的二阶条件。令 $F\in\mathbf{R}^{n\times m}$，$\hat x\in\mathbf{R}^n$。函数 $f:\mathbf{R}^n\to\mathbf{R}$ 在仿射集 $\{Fz+\hat x\mid z\in\mathbf{R}^m\}$ 上的限制定义为函数 $\tilde f:\mathbf{R}^m\to\mathbf{R}$，其中 $\tilde f(z)=f(Fz+\hat x)$，$\operatorname{dom}\tilde f=\{z\mid Fz+\hat x\in\operatorname{dom}f\}$。假设 $f$ 在凸域上二次可微。 (a) 证明 $\tilde f$ 是凸的，当且仅当对所有 $z\in\operatorname{dom}\tilde f$ 有 $F^T\nabla^2 f(Fz+\hat x)F\succeq0$。 (b) 假设 $A\in\mathbf{R}^{p\times n}$ 是一个矩阵，其零空间等于 $F$ 的值域，即 $AF=0$ 且 $\operatorname{rank}A=n-\operatorname{rank}F$。证明如果对所有 $z\in\operatorname{dom}\tilde f$ 都存在 $\lambda\in\mathbf{R}$ 使得 $\nabla^2 f(Fz+\hat x)+\lambda A^TA\succeq0$，则 $\tilde f$ 是凸的。提示：使用以下结果：如果 $B\in\mathbf{S}^n$ 且 $A\in\mathbf{R}^{p\times n}$，则 $x^TBx\ge0$ 对所有 $x\in\mathcal{N}(A)$ 成立，当且仅当存在 $\lambda$ 使得 $B+\lambda A^TA\succeq0$。

3.10 Jensen 不等式的扩展。Jensen 不等式的一种解释是随机化或抖动会造成损害，即提高凸函数的平均值：对于凸函数 $f$ 和零均值随机变量 $v$，有 $\mathbf{E}f(x_0+v)\ge f(x_0)$。这引出下面的猜想：如果 $f$ 是凸函数，那么 $v$ 的方差越大，$\mathbf{E}f(x_0+v)$ 就越大。 (a) 举出一个反例证明这个猜想是错误的。找到零均值随机变量 $v$ 和 $w$，其中 $\operatorname{var}(v)>\operatorname{var}(w)$，以及凸函数 $f$ 和点 $x_0$，使得 $\mathbf{E}f(x_0+v)<\mathbf{E}f(x_0+w)$。


(b) 当 $v$ 和 $w$ 是彼此的缩放版本时，猜想成立。证明当 $f$ 为凸且 $v$ 为零均值时，$\mathbf{E}f(x_0+tv)$ 关于 $t\ge0$ 单调递增。

3.11 单调映射。如果对所有 $x,y\in\operatorname{dom}\psi$ 有 $(\psi(x)-\psi(y))^T(x-y)\ge0$，则函数 $\psi:\mathbf{R}^n\to\mathbf{R}^n$ 称为单调的。（注意，此处定义的“单调”与第 3.6.1 节中给出的定义不同；这两个定义都被广泛使用。）假设 $f:\mathbf{R}^n\to\mathbf{R}$ 是可微凸函数。证明它的梯度 $\nabla f$ 是单调映射。反之是否成立，即每个单调映射都是某个凸函数的梯度吗？

3.12 假设 $f:\mathbf{R}^n\to\mathbf{R}$ 是凸的，$g:\mathbf{R}^n\to\mathbf{R}$ 是凹的，$\operatorname{dom}f=\operatorname{dom}g=\mathbf{R}^n$，并且对所有 $x$ 有 $g(x)\le f(x)$。证明存在一个仿射函数 $h$，使得对所有 $x$ 有 $g(x)\le h(x)\le f(x)$。换句话说，如果凹函数 $g$ 是凸函数 $f$ 的低估量，那么可以在 $f$ 和 $g$ 之间拟合一个仿射函数。

3.13 Kullback-Leibler 散度和信息不等式。令 $D_{\mathrm{kl}}$ 为 (3.17) 中定义的 Kullback-Leibler 散度。证明信息不等式：对所有 $u,v\in\mathbf{R}^n_{++}$，$D_{\mathrm{kl}}(u,v)\ge0$。还要证明 $D_{\mathrm{kl}}(u,v)=0$ 当且仅当 $u=v$。提示：Kullback-Leibler 散度可以表示为 $D_{\mathrm{kl}}(u,v)=f(u)-f(v)-\nabla f(v)^T(u-v)$，其中 $f(v)=\sum_{i=1}^n v_i\log v_i$ 是 $v$ 的负熵。

3.14 凸凹函数和鞍点。如果对每个固定的 $x$，$f(x,z)$ 是 $z$ 的凹函数，并且对每个固定的 $z$，$f(x,z)$ 是 $x$ 的凸函数，则称函数 $f:\mathbf{R}^n\times\mathbf{R}^m\to\mathbf{R}$ 是凸凹函数。还要求其定义域具有乘积形式 $\operatorname{dom}f=A\times B$，其中 $A\subseteq\mathbf{R}^n$ 和 $B\subseteq\mathbf{R}^m$ 是凸集。 (a) 用 Hessian $\nabla^2 f(x,z)$ 给出二次可微函数 $f:\mathbf{R}^n\times\mathbf{R}^m\to\mathbf{R}$ 为凸凹函数的二阶条件。 (b) 假设 $f:\mathbf{R}^n\times\mathbf{R}^m\to\mathbf{R}$ 是凸凹且可微的，且 $\nabla f(\tilde x,\tilde z)=0$。证明鞍点性质成立：对所有 $x,z$，有 $f(\tilde x,z)\le f(\tilde x,\tilde z)\le f(x,\tilde z)$。说明这意味着 $f$ 满足强 max-min 性质：$\sup_z\inf_x f(x,z)=\inf_x\sup_z f(x,z)$，其共同值为 $f(\tilde x,\tilde z)$。 (c) 现在假设 $f:\mathbf{R}^n\times\mathbf{R}^m\to\mathbf{R}$ 可微，但不一定是凸凹的，并且对所有 $x,z$，鞍点性质在 $(\tilde x,\tilde z)$ 处成立：$f(\tilde x,z)\le f(\tilde x,\tilde z)\le f(x,\tilde z)$。证明 $\nabla f(\tilde x,\tilde z)=0$。

3.15 一族凹效用函数。对于 $0<\alpha\le1$，设 $u_\alpha(x)=(x^\alpha-1)/\alpha$，其中 $\operatorname{dom}u_\alpha=\mathbf{R}_+$。还定义 $u_0(x)=\log x$，其中 $\operatorname{dom}u_0=\mathbf{R}_{++}$。 (a) 证明当 $x>0$ 时，$u_0(x)=\lim_{\alpha\to0}u_\alpha(x)$。 (b) 证明 $u_\alpha$ 是凹的、单调递增的，并且都满足 $u_\alpha(1)=0$。这些函数在经济学中常用于建模一定数量商品或货币的收益或效用。$u_\alpha$ 的凹性意味着边际效用（即固定增加商品所获得的效用增加）随着商品数量的增加而减少。换句话说，凹性模拟了饱和效应。

3.16 对于下列每个函数，确定它是否是凸函数、凹函数、拟凸函数或拟凹函数。 (a) $\mathbf{R}$ 上的 $f(x)=e^x-1$。 (b) $\mathbf{R}^2_{++}$ 上的 $f(x_1,x_2)=x_1x_2$。 (c) $\mathbf{R}^2_{++}$ 上的 $f(x_1,x_2)=1/(x_1x_2)$。 (d) $\mathbf{R}^2_{++}$ 上的 $f(x_1,x_2)=x_1/x_2$。 (e) $\mathbf{R}\times\mathbf{R}_{++}$ 上的 $f(x_1,x_2)=x_1^2/x_2$。 (f) $\mathbf{R}^2_{++}$ 上的 $f(x_1,x_2)=x_1^\alpha x_2^{1-\alpha}$，其中 $0\le\alpha\le1$。

3.17 假设 $p<1$，$p\ne0$。证明函数 $f(x)=(\sum_{i=1}^n x_i^p)^{1/p}$ 且 $\operatorname{dom}f=\mathbf{R}^n_{++}$ 是凹函数。这包括特殊情况 $f(x)=(\sum_{i=1}^n x_i^{1/2})^2$ 和调和平均值 $f(x)=(\sum_{i=1}^n 1/x_i)^{-1}$。提示：修改第 3.1.5 节中 log-sum-exp 函数和几何平均值的证明。

3.18 修改第 3.1.5 节中对数行列式函数凹性的证明，以说明以下内容。 (a) $f(X)=\operatorname{tr}(X^{-1})$ 在 $\operatorname{dom}f=\mathbf{S}^n_{++}$ 上是凸的。 (b) $f(X)=(\det X)^{1/n}$ 在 $\operatorname{dom}f=\mathbf{S}^n_{++}$ 上是凹的。

3.19 非负加权和和积分。 (a) 证明 $f(x)=\sum_{i=1}^r\alpha_i x_{[i]}$ 是 $x$ 的凸函数，其中 $\alpha_1\ge\alpha_2\ge\cdots\ge\alpha_r\ge0$，$x_{[i]}$ 表示 $x$ 的第 $i$ 个最大分量。（可以使用 $f(x)=\sum_{i=1}^k x_{[i]}$ 在 $\mathbf{R}^n$ 上凸的事实。） (b) 令 $T(x,\omega)$ 表示三角多项式 $T(x,\omega)=x_1+x_2\cos\omega+x_3\cos2\omega+\cdots+x_n\cos(n-1)\omega$。证明函数 $f(x)=-\int_0^{2\pi}\log T(x,\omega)\,d\omega$ 在 $\{x\in\mathbf{R}^n\mid T(x,\omega)>0,\ 0\le\omega\le2\pi\}$ 上是凸的。

3.20 与仿射函数复合。证明以下函数 $f:\mathbf{R}^n\to\mathbf{R}$ 是凸函数。 (a) $f(x)=\|Ax-b\|$，其中 $A\in\mathbf{R}^{m\times n}$，$b\in\mathbf{R}^m$，$\|\cdot\|$ 是 $\mathbf{R}^m$ 上的范数。 (b) $f(x)=-(\det(A_0+x_1A_1+\cdots+x_nA_n))^{1/m}$，定义在 $\{x\mid A_0+x_1A_1+\cdots+x_nA_n\succ0\}$ 上，其中 $A_i\in\mathbf{S}^m$。 (c) $f(x)=\operatorname{tr}(A_0+x_1A_1+\cdots+x_nA_n)^{-1}$，定义在 $\{x\mid A_0+x_1A_1+\cdots+x_nA_n\succ0\}$ 上，其中 $A_i\in\mathbf{S}^m$。（利用 $\operatorname{tr}(X^{-1})$ 在 $\mathbf{S}^m_{++}$ 上凸的事实；参见练习 3.18。）


3.21 逐点最大值和上确界。证明以下函数 $f:\mathbf{R}^n\to\mathbf{R}$ 是凸函数。 (a) $f(x)=\max_{i=1,\ldots,k}\|A^{(i)}x-b^{(i)}\|$，其中 $A^{(i)}\in\mathbf{R}^{m\times n}$、$b^{(i)}\in\mathbf{R}^m$，且 $\|\cdot\|$ 是 $\mathbf{R}^m$ 上的范数。 (b) $\mathbf{R}^n$ 上的 $f(x)=\sum_{i=1}^r |x|_{[i]}$，其中 $|x|$ 表示向量 $|x|_i=|x_i|$，并且 $|x|_{[i]}$ 是 $|x|$ 的第 $i$ 个最大分量。换句话说，$|x|_{[1]}, |x|_{[2]}, \ldots, |x|_{[n]}$ 是 $x$ 分量的绝对值按非递增顺序排序后的结果。

3.22 组合规则。证明下列函数是凸函数。 (a) $f(x)=-\log\left(-\log\left(\sum_{i=1}^m e^{a_i^Tx+b_i}\right)\right)$，定义域为 $\operatorname{dom}f=\{x\mid \sum_{i=1}^m e^{a_i^Tx+b_i}<1\}$。可以使用 $\log(\sum_{i=1}^n e^{y_i})$ 是凸函数这一事实。 (b) $f(x,u,v)=-\sqrt{uv-x^Tx}$，定义域为 $\operatorname{dom}f=\{(x,u,v)\mid uv>x^Tx,\ u,v>0\}$。当 $u>0$ 时，$x^Tx/u$ 关于 $(x,u)$ 是凸的，并且 $-\sqrt{x_1x_2}$ 在 $\mathbf{R}^2_{++}$ 上是凸的。 (c) $f(x,u,v)=-\log(uv-x^Tx)$，定义域为 $\operatorname{dom}f=\{(x,u,v)\mid uv>x^Tx,\ u,v>0\}$。 (d) $f(x,t)=-(t^p-\|x\|_p^p)^{1/p}$，其中 $p>1$，$\operatorname{dom}f=\{(x,t)\mid t\ge\|x\|_p\}$。可以使用如下事实：当 $u>0$ 时，$\|x\|_p^p/u^{p-1}$ 关于 $(x,u)$ 是凸的（参见练习 3.23），并且 $-x^{1/p}y^{1-1/p}$ 在 $\mathbf{R}^2_+$ 上是凸的（参见练习 3.16）。 (e) $f(x,t)=-\log(t^p-\|x\|_p^p)$，其中 $p>1$，$\operatorname{dom}f=\{(x,t)\mid t>\|x\|_p\}$。可以使用如下事实：当 $u>0$ 时，$\|x\|_p^p/u^{p-1}$ 关于 $(x,u)$ 是凸的（参见练习 3.23）。

3.23 函数的透视。 (a) 证明对于 $p>1$，$f(x,t)=(|x_1|^p+\cdots+|x_n|^p)/t^{p-1}=\|x\|_p^p/t^{p-1}$ 在 $\{(x,t)\mid t>0\}$ 上是凸的。 (b) 证明 $f(x)=\|Ax+b\|_2^2/(c^Tx+d)$ 在 $\{x\mid c^Tx+d>0\}$ 上是凸的，其中 $A\in\mathbf{R}^{m\times n}$、$b\in\mathbf{R}^m$、$c\in\mathbf{R}^n$ 和 $d\in\mathbf{R}$。

3.24 概率单纯形上的一些函数。令 $x$ 为实值随机变量，其取值于 $\{a_1,\ldots,a_n\}$，其中 $a_1<a_2<\cdots<a_n$，且 $\operatorname{prob}(x=a_i)=p_i$，$i=1,\ldots,n$。对于 $p$ 的以下每个函数（在概率单纯形 $\{p\in\mathbf{R}^n_+\mid \mathbf{1}^Tp=1\}$ 上），确定该函数是凸函数、凹函数、拟凸函数还是拟凹函数。 (a) $\mathbf{E}x$。 (b) $\operatorname{prob}(x\ge\alpha)$。 (c) $\operatorname{prob}(\alpha\le x\le\beta)$。 (d) $\sum_{i=1}^n p_i\log p_i$，即分布的负熵。 (e) $\operatorname{var}x=\mathbf{E}(x-\mathbf{E}x)^2$。 (f) $\operatorname{quartile}(x)=\inf\{\beta\mid \operatorname{prob}(x\le\beta)\ge0.25\}$。 (g) 最小集合 $A\subseteq\{a_1,\ldots,a_n\}$ 的基数，使其概率至少为 $90\%$。（这里基数指 $A$ 中元素的数量。） (h) 包含 $90\%$ 概率的最小宽度区间，即 $\inf\{\beta-\alpha\mid \operatorname{prob}(\alpha\le x\le\beta)\ge0.9\}$。

3.25 分布之间的最大概率距离。令 $p,q\in\mathbf{R}^n$ 表示 $\{1,\ldots,n\}$ 上的两个概率分布（因此 $p,q\succeq0$，$\mathbf{1}^Tp=\mathbf{1}^Tq=1$）。我们将 $p$ 和 $q$ 之间的最大概率距离 $d_{\mathrm{mp}}(p,q)$ 定义为 $p$ 和 $q$ 在所有事件上分配概率的最大差异：

$$
d_{\mathrm{mp}}(p,q)=\max\{|\operatorname{prob}(p,C)-\operatorname{prob}(q,C)|\mid C\subseteq\{1,\ldots,n\}\}.
$$

这里 $\operatorname{prob}(p,C)$ 是分布 $p$ 下事件 $C$ 的概率，即 $\operatorname{prob}(p,C)=\sum_{i\in C}p_i$。找到 $d_{\mathrm{mp}}$ 的一个简单表达式，其中包含 $\|p-q\|_1=\sum_{i=1}^n |p_i-q_i|$，并证明 $d_{\mathrm{mp}}$ 是 $\mathbf{R}^n\times\mathbf{R}^n$ 上的凸函数。（它的定义域是 $\{(p,q)\mid p,q\succeq0,\ \mathbf{1}^Tp=\mathbf{1}^Tq=1\}$，但它有一个到所有 $\mathbf{R}^n\times\mathbf{R}^n$ 的自然扩展。）

3.26 更多特征值函数。令 $\lambda_1(X)\ge\lambda_2(X)\ge\cdots\ge\lambda_n(X)$ 表示矩阵 $X\in\mathbf{S}^n$ 的特征值。我们已经看到几个关于特征值的函数是 $X$ 的凸函数或凹函数。最大特征值 $\lambda_1(X)$ 是凸函数（例 3.10），最小特征值 $\lambda_n(X)$ 是凹函数；特征值（或迹）之和 $\operatorname{tr}X=\lambda_1(X)+\cdots+\lambda_n(X)$ 是线性的；特征值倒数之和（或逆矩阵的迹）$\operatorname{tr}(X^{-1})=\sum_{i=1}^n 1/\lambda_i(X)$ 在 $\mathbf{S}^n_{++}$ 上是凸的（练习 3.18）；特征值的几何平均 $(\det X)^{1/n}=(\prod_{i=1}^n\lambda_i(X))^{1/n}$ 和特征值乘积的对数 $\log\det X=\sum_{i=1}^n\log\lambda_i(X)$ 在 $X\in\mathbf{S}^n_{++}$ 上是凹的（练习

3.18 和第 74 页）。在这个问题中，我们利用变分表征来探索更多特征值函数。 (a) $k$ 个最大特征值的和。证明 $\sum_{i=1}^k\lambda_i(X)$ 在 $\mathbf{S}^n$ 上是凸的。提示：[HJ85，第 191 页] 使用变分表征 $\sum_{i=1}^k\lambda_i(X)=\sup\{\operatorname{tr}(V^TXV)\mid V\in\mathbf{R}^{n\times k},\ V^TV=I\}$。 (b) $k$ 个最小特征值的几何平均。证明 $(\prod_{i=n-k+1}^n\lambda_i(X))^{1/k}$ 在 $\mathbf{S}^n_{++}$ 上是凹的。提示：[MO79，第 513 页] 对于 $X\succ0$，有 $(\prod_{i=n-k+1}^n\lambda_i(X))^{1/k}=(1/k)\inf\{\operatorname{tr}(V^TXV)\mid V\in\mathbf{R}^{n\times k},\ \det V^TV=1\}$。 (c) $k$ 个最小特征值乘积的对数。证明 $\sum_{i=n-k+1}^n\log\lambda_i(X)$ 在 $\mathbf{S}^n_{++}$ 上是凹的。提示：[MO79，第 513 页] 对于 $X\succ0$，$\prod_{i=n-k+1}^n\lambda_i(X)=\inf\{\prod_{i=1}^k(V^TXV)_{ii}\mid V\in\mathbf{R}^{n\times k},\ V^TV=I\}$。

3.27 Cholesky 因子的对角线元素。每个 $X\in\mathbf{S}^n_{++}$ 具有唯一的 Cholesky 分解 $X=LL^T$，其中 $L$ 是下三角矩阵，且 $L_{ii}>0$。证明 $L_{ii}$ 是 $X$ 的凹函数（定义域为 $\mathbf{S}^n_{++}$）。提示：$L_{ii}$ 可表示为 $L_{ii}=(w-z^TY^{-1}z)^{1/2}$，其中 $\begin{bmatrix}Y&z\\ z^T&w\end{bmatrix}$ 是 $X$ 的前导 $i\times i$ 子矩阵。


### 保持凸性的运算

3.28 将凸函数表示为仿射函数族的逐点上确界。在这个问题中，我们将第 83 页证明的结果扩展到 $\operatorname{dom}f\ne\mathbf{R}^n$ 的情况。设 $f:\mathbf{R}^n\to\mathbf{R}$ 是凸函数。将 $\tilde f:\mathbf{R}^n\to\mathbf{R}$ 定义为所有仿射函数的逐点上确界，这些仿射函数是 $f$ 的全局低估量：

$$
\tilde f(x)=\sup\{g(x)\mid g \text{ 仿射，且对所有 }z,\ g(z)\le f(z)\}.
$$

(a) 证明对 $x\in\operatorname{int}\operatorname{dom}f$ 有 $f(x)=\tilde f(x)$。 (b) 如果 $f$ 是闭的（即 $\operatorname{epi}f$ 是闭集；参见 §A.3.3），则证明 $f=\tilde f$。

3.29 分段线性凸函数的表示。如果存在 $\mathbf{R}^n$ 的划分 $\mathbf{R}^n=X_1\cup X_2\cup\cdots\cup X_L$，其中 $\operatorname{int}X_i\ne\emptyset$，且当 $i\ne j$ 时 $\operatorname{int}X_i\cap\operatorname{int}X_j=\emptyset$，并且存在一组仿射函数 $a_1^Tx+b_1,\ldots,a_L^Tx+b_L$，使得对 $x\in X_i$ 有 $f(x)=a_i^Tx+b_i$，则凸函数 $f:\mathbf{R}^n\to\mathbf{R}$（$\operatorname{dom}f=\mathbf{R}^n$）称为分段线性函数。证明这样的函数可表示为

$$
f(x)=\max\{a_1^Tx+b_1,\ldots,a_L^Tx+b_L\}.
$$

3.30 函数的凸包或凸包络。函数 $f:\mathbf{R}^n\to\mathbf{R}$ 的凸包或凸包络定义为 $g(x)=\inf\{t\mid (x,t)\in\operatorname{conv}\operatorname{epi}f\}$。在几何上，$g$ 的上图是 $f$ 的上图的凸包。证明 $g$ 是 $f$ 的最大凸低估量。换句话说，证明如果 $h$ 是凸的，并且对所有 $x$ 满足 $h(x)\le f(x)$，则对所有 $x$ 满足 $h(x)\le g(x)$。

3.31 [Roc70，第 35 页] 最大齐次低估量。令 $f$ 为凸函数。将函数 $g$ 定义为

$$
g(x)=\inf_{\alpha>0}\frac{f(\alpha x)}{\alpha}.
$$

(a) 证明 $g$ 是齐次的，即对所有 $t\ge0$，$g(tx)=tg(x)$。 (b) 证明 $g$ 是 $f$ 的最大齐次低估量：如果 $h$ 是齐次的，并且对所有 $x$ 有 $h(x)\le f(x)$，那么对所有 $x$ 有 $h(x)\le g(x)$。 (c) 证明 $g$ 是凸的。

3.32 凸函数的乘积和比率。一般来说，两个凸函数的乘积或比率不是凸函数。然而，有一些结果适用于 R 上的函数。证明以下内容。 (a) 如果 f 和 g 是凸函数，且均为非递减（或非递增）且区间上的正函数，则 fg 是凸函数。 (b) 如果 f 、 g 是凹的、正的，其中一个不减，另一个不增，则 fg 是凹的。 (c) 如果 f 是凸的、非减的、正的，并且 g 是凹的、非增的、正的，则 $f/g$ 是凸的。

3.33 透视定理的直接证明。直接证明凸函数 $f$ 的透视函数 $g$（如第 3.2.6 节所定义）是凸的：证明 $\operatorname{dom}g$ 是凸集，并且对于 $(x,t),(y,s)\in\operatorname{dom}g$ 与 $0\le\theta\le1$，有

$$
g(\theta x+(1-\theta)y,\theta t+(1-\theta)s)
\le \theta g(x,t)+(1-\theta)g(y,s).
$$

3.34 Minkowski 函数。凸集 $C$ 的 Minkowski 函数定义为 $M_C(x)=\inf\{t>0\mid t^{-1}x\in C\}$。 (a) 画图给出如何求 $M_C(x)$ 的几何解释。 (b) 证明 $M_C$ 是齐次的，即当 $\alpha\ge0$ 时 $M_C(\alpha x)=\alpha M_C(x)$。 (c) $\operatorname{dom}M_C$ 是什么？ (d) 证明 $M_C$ 是凸函数。 (e) 假设 $C$ 也是闭的、有界的、对称的（如果 $x\in C$ 则 $-x\in C$），并且具有非空内部。证明 $M_C$ 是一个范数。对应的单位球是什么？

3.35 支持函数演算。回想一下，集合 $C\subseteq\mathbf{R}^n$ 的支持函数定义为 $S_C(y)=\sup\{y^Tx\mid x\in C\}$。在第 81 页，我们证明了 $S_C$ 是凸函数。 (a) 证明 $S_B=S_{\operatorname{conv}B}$。 (b) 证明 $S_{A+B}=S_A+S_B$。 (c) 证明 $S_{A\cup B}=\max\{S_A,S_B\}$。 (d) 设 $B$ 是闭凸集。证明 $A\subseteq B$ 当且仅当对所有 $y$ 有 $S_A(y)\le S_B(y)$。

### 共轭函数

3.36 推导下列函数的共轭。 (a) 最大函数：$\mathbf{R}^n$ 上的 $f(x)=\max_{i=1,\ldots,n}x_i$。 (b) 最大元素之和：$\mathbf{R}^n$ 上的 $f(x)=\sum_{i=1}^r x_{[i]}$。 (c) $\mathbf{R}$ 上的分段线性函数：$f(x)=\max_{i=1,\ldots,m}(a_ix+b_i)$。可以假设 $a_i$ 按升序排序，即 $a_1\le\cdots\le a_m$，并且函数 $a_ix+b_i$ 中没有一个是冗余的，即对每个 $k$ 至少有一个 $x$ 使 $f(x)=a_kx+b_k$。 (d) 幂函数：$\mathbf{R}_{++}$ 上的 $f(x)=x^p$，其中 $p>1$；并对 $p<0$ 重复。 (e) 负几何平均值：$\mathbf{R}^n_{++}$ 上的 $f(x)=-(\prod_i x_i)^{1/n}$。 (f) 二阶锥的负广义对数：$f(x,t)=-\log(t^2-x^Tx)$，定义域为 $\{(x,t)\in\mathbf{R}^n\times\mathbf{R}\mid \|x\|_2<t\}$。

3.37 证明 $f(X)=\operatorname{tr}(X^{-1})$、$\operatorname{dom}f=\mathbf{S}^n_{++}$ 的共轭为 $f^*(Y)=-2\operatorname{tr}((-Y)^{1/2})$，$\operatorname{dom}f^*=-\mathbf{S}^n_+$。提示：$f$ 的梯度为 $\nabla f(X)=-X^{-2}$。

3.38 Young 不等式。令 $f:\mathbf{R}\to\mathbf{R}$ 为增函数，其中 $f(0)=0$，并令 $g$ 为其反函数。定义

$$
F(x)=\int_0^x f(a)\,da,\qquad
G(y)=\int_0^y g(a)\,da.
$$

证明 $F$ 和 $G$ 是共轭函数。给出 Young 不等式 $xy\le F(x)+G(y)$ 的简单图形解释。

3.39 共轭函数的性质。 (a) 凸函数加仿射函数的共轭。定义 $g(x)=f(x)+c^Tx+d$，其中 $f$ 是凸的。用 $f^*$（以及 $c,d$）表示 $g^*$。 (b) 透视函数的共轭。用 $f^*$ 表示凸函数 $f$ 的透视函数的共轭。


(c) 共轭和最小化。设 $f(x,z)$ 关于 $(x,z)$ 是凸函数，并定义 $g(x)=\inf_z f(x,z)$。用 $f^*$ 表示 $g^*$。作为一个应用，用 $h^*$、$A$ 和 $b$ 表示 $g(x)=\inf_z\{h(z)\mid Az+b=x\}$ 的共轭，其中 $h$ 是凸的。 (d) 共轭的共轭。证明闭凸函数的共轭的共轭是它本身：如果 $f$ 是闭凸函数，则 $f=f^{\ast\ast}$。（如果函数的上图是闭的，则该函数是闭的；参见 §A.3.3。）提示：证明 $f^{\ast\ast}$ 是 $f$ 的所有仿射全局低估量的逐点上确界。然后应用练习 3.28 的结果。

3.40 共轭函数的梯度和 Hessian。假设 $f:\mathbf{R}^n\to\mathbf{R}$ 是凸且两次连续可导的。假设 $\bar y$ 和 $\bar x$ 通过 $\bar y=\nabla f(\bar x)$ 相关，并且 $\nabla^2 f(\bar x)\succ0$。 (a) 证明 $\nabla f^*(\bar y)=\bar x$。 (b) 证明 $\nabla^2 f^*(\bar y)=\nabla^2 f(\bar x)^{-1}$。

3.41 负归一化熵的共轭。证明负归一化熵

$$
f(x)=\sum_{i=1}^n x_i\log\frac{x_i}{\mathbf{1}^Tx},
\qquad \operatorname{dom}f=\mathbf{R}^n_{++},
$$

的共轭为

$$
f^*(y)=
\begin{cases}
0, & \sum_{i=1}^n e^{y_i}\le1,\\
\infty, & \text{否则}.
\end{cases}
$$

### 拟凸函数

3.42 近似宽度。设 $f_0,\ldots,f_n:\mathbf{R}\to\mathbf{R}$ 为连续函数。我们考虑用 $f_1,\ldots,f_n$ 的线性组合近似 $f_0$ 的问题。对于 $x\in\mathbf{R}^n$，称 $f=x_1f_1+\cdots+x_nf_n$ 在区间 $[0,T]$ 上以容差 $\epsilon>0$ 近似 $f_0$，如果对 $0\le t\le T$ 有 $|f(t)-f_0(t)|\le\epsilon$。现在固定一个容差 $\epsilon>0$，并将近似宽度定义为使 $f$ 在区间 $[0,T]$ 上近似 $f_0$ 的最大 $T$：

$$
W(x)=\sup\{T\mid |x_1f_1(t)+\cdots+x_nf_n(t)-f_0(t)|\le\epsilon,\ 0\le t\le T\}.
$$

证明 $W$ 是拟凹的。

3.43 拟凸性的一阶条件。证明第 3.4.3 节中给出的拟凸性一阶条件：具有凸定义域的可微函数 $f:\mathbf{R}^n\to\mathbf{R}$ 是拟凸的，当且仅当对所有 $x,y\in\operatorname{dom}f$，

$$
f(y)\le f(x)\Longrightarrow \nabla f(x)^T(y-x)\le0.
$$

提示：足以证明 $\mathbf{R}$ 上函数的结果；一般结果可通过限制到任意直线得到。

3.44 拟凸性的二阶条件。在这个问题中，我们推导第 3.4.3 节中给出的拟凸性二阶条件的替代表示。

(a) 证明：如果存在 $\sigma$ 使得

$$
\nabla^2 f(x)+\sigma\nabla f(x)\nabla f(x)^T\succeq 0, \tag{3.26}
$$

则点 $x\in\operatorname{dom} f$ 满足 (3.21)。对所有 $y\ne0$，它满足 (3.22) 当且仅当存在 $\sigma$ 使得

$$
\nabla^2 f(x)+\sigma\nabla f(x)\nabla f(x)^T\succ 0. \tag{3.27}
$$

提示：可以不失一般性地假设 $\nabla^2 f(x)$ 是对角矩阵。

(b) 点 $x\in\operatorname{dom} f$ 满足 (3.21)，当且仅当 $\nabla f(x)=0$ 且 $\nabla^2 f(x)\succeq0$，或者 $\nabla f(x)\ne0$ 且矩阵

$$
H(x)=\begin{bmatrix}
\nabla^2 f(x) & \nabla f(x)\\
\nabla f(x)^T & 0
\end{bmatrix}
$$

恰好有一个负特征值。对所有 $y\ne0$，它满足 (3.22) 当且仅当 $H(x)$ 恰好有一个非正特征值。提示：可以使用 (a) 部分的结果。线性代数中的特征值交错定理也可能有用：若 $B\in\mathbf{S}^n$ 且 $a\in\mathbf{R}^n$，则

$$
\lambda_n\left(\begin{bmatrix}B & a\\ a^T & 0\end{bmatrix}\right)\ge \lambda_n(B).
$$

3.45 使用第 3.4.3 节中给出的拟凸性一阶和二阶条件来验证函数 $f(x)=-x_1x_2$ 的拟凸性，其中 $\operatorname{dom}f=\mathbf{R}^2_{++}$。

3.46 定义域为 $\mathbf{R}^n$ 的拟线性函数。$\mathbf{R}$ 上的拟线性函数（即拟凸且拟凹）是单调的，即非减或非增。在这个问题中，我们考虑将此结果推广到 $\mathbf{R}^n$ 上的函数。假设函数 $f:\mathbf{R}^n\to\mathbf{R}$ 是拟线性且连续的，且 $\operatorname{dom}f=\mathbf{R}^n$。证明它可以表示为 $f(x)=g(a^Tx)$，其中 $g:\mathbf{R}\to\mathbf{R}$ 是单调函数，$a\in\mathbf{R}^n$。换句话说，定义域为 $\mathbf{R}^n$ 的拟线性函数必须是线性函数的单调函数。（反之亦然。）

### 对数凹函数和对数凸函数

3.47 假设 $f:\mathbf{R}^n\to\mathbf{R}$ 是可微的，$\operatorname{dom}f$ 是凸的，并且对所有 $x\in\operatorname{dom}f$ 有 $f(x)>0$。证明 $f$ 是对数凹的，当且仅当对所有 $x,y\in\operatorname{dom}f$，

$$
\frac{f(y)}{f(x)}
\le \exp\left(\frac{\nabla f(x)^T(y-x)}{f(x)}\right).
$$

3.48 证明如果 $f:\mathbf{R}^n\to\mathbf{R}$ 是对数凹且 $a\ge0$，则函数 $g=f-a$ 是对数凹的，其中 $\operatorname{dom}g=\{x\in\operatorname{dom}f\mid f(x)>a\}$。

3.49 证明下列函数是对数凹函数。 (a) Logistic 函数：$f(x)=e^x/(1+e^x)$，$\operatorname{dom}f=\mathbf{R}$。 (b) 调和平均值：$f(x)=1/(1/x_1+\cdots+1/x_n)$，$\operatorname{dom}f=\mathbf{R}^n_{++}$。 (c) 乘积：$f(x)=\prod_{i=1}^n x_i/(\sum_{i=1}^n x_i)^n$，$\operatorname{dom}f=\mathbf{R}^n_{++}$。 (d) 迹上的行列式：$f(X)=\det X/\operatorname{tr}X$，$\operatorname{dom}f=\mathbf{S}^n_{++}$。


3.50 多项式的系数作为根的函数。证明具有实数负根的多项式的系数是根的对数凹函数。换句话说，由恒等式

$$
s^n+a_1(\lambda)s^{n-1}+\cdots+a_{n-1}(\lambda)s+a_n(\lambda)
=(s-\lambda_1)(s-\lambda_2)\cdots(s-\lambda_n)
$$

定义的函数 $a_i:\mathbf{R}^n\to\mathbf{R}$，在 $-\mathbf{R}^n_{++}$ 上是对数凹函数。提示：函数

$$
S_k(x)=\sum_{1\le i_1<i_2<\cdots<i_k\le n}x_{i_1}x_{i_2}\cdots x_{i_k},
$$

其中 $\operatorname{dom}S_k=\mathbf{R}^n_+$ 且 $1\le k\le n$，称为 $\mathbf{R}^n$ 上的第 $k$ 个初等对称函数。可以证明 $S_k^{1/k}$ 是凹的（参见 [ML57]）。

3.51 [BL00，第 41 页] 设 $p$ 为 $\mathbf{R}$ 上的多项式，其所有根均为实数。证明它在任何为正的区间上都是对数凹的。

3.52 [MO79，§3.E.2] 矩函数的对数凸性。假设 $f:\mathbf{R}\to\mathbf{R}$ 是非负函数，且 $\mathbf{R}_+\subseteq\operatorname{dom}f$。对于 $x\ge0$，定义 $\phi(x)=\int_0^\infty u^x f(u)\,du$。证明 $\phi$ 是对数凸函数。（如果 $x$ 是正整数，$f$ 是概率密度函数，则 $\phi(x)$ 是分布的 $x$ 阶矩。）用它来证明 Gamma 函数 $\Gamma(x)=\int_0^\infty u^{x-1}e^{-u}\,du$ 在 $x\ge1$ 时是对数凸函数。

3.53 假设 x 和 y 是 $\mathbf{R}^n$ 中的独立随机向量，分别具有对数凹概率密度函数 f 和 g。证明 $z =x +y$ 之和的概率密度函数是对数凹的。

3.54 高斯累积分布函数的对数凹性。高斯随机变量的累积分布函数 $f(x)=(1/\sqrt{2\pi})\int_{-\infty}^x e^{-t^2/2}\,dt$ 是对数凹函数。这是从两个对数凹函数的卷积是对数凹函数的一般结果得出的。在这个问题中，我们将给出一个简单的独立证明，证明 $f$ 是对数凹的。回想一下，$f$ 是对数凹的，当且仅当对所有 $x$ 有 $f''(x)f(x)\le f'(x)^2$。 (a) 验证当 $x\ge0$ 时 $f''(x)f(x)\le f'(x)^2$。这留下了困难的部分，即证明 $x<0$ 时的不等式。 (b) 验证对任意 $t$ 和 $x$，有 $t^2/2\ge -x^2/2+xt$。 (c) 利用 (b) 表明 $e^{-t^2/2}\le e^{x^2/2-xt}$。得出结论：当 $x<0$ 时，$\int_{-\infty}^x e^{-t^2/2}\,dt\le e^{x^2/2}\int_{-\infty}^x e^{-xt}\,dt$。 (d) 利用 (c) 验证当 $x\le0$ 时 $f''(x)f(x)\le f'(x)^2$。

3.55 对数凹概率密度的累积分布函数的对数凹性。在本题中，我们扩展练习 3.54 的结果。令 $g(t)=\exp(-h(t))$ 为可微对数凹概率密度函数，并令 $f(x)=\int_{-\infty}^x g(t)\,dt=\int_{-\infty}^x e^{-h(t)}\,dt$ 为其累积分布函数。我们将证明 $f$ 是对数凹的，即对所有 $x$，它满足 $f''(x)f(x)\le (f'(x))^2$。 (a) 用函数 $h$ 表示 $f$ 的导数。验证如果 $h'(x)\ge0$，则 $f''(x)f(x)\le(f'(x))^2$。 (b) 假设 $h'(x)<0$。使用由 $h$ 的凸性得出的不等式 $h(t)\ge h(x)+h'(x)(t-x)$，证明 $\int_{-\infty}^x e^{-h(t)}\,dt\le e^{-h(x)}/(-h'(x))$。使用此不等式验证如果 $h'(x)<0$，则 $f''(x)f(x)\le(f'(x))^2$。

3.56 更多的对数凹密度。证明下列密度是对数凹的。 (a) [MO79，第 493 页] 伽玛密度，定义为 $f(x)=\alpha^\lambda\Gamma(\lambda)^{-1}x^{\lambda-1}e^{-\alpha x}$，其中 $\operatorname{dom}f=\mathbf{R}_+$。参数 $\lambda$ 和 $\alpha$ 满足 $\lambda\ge1$、$\alpha>0$。 (b) [MO79，第 306 页] Dirichlet 密度

$$
f(x)=\frac{\Gamma(\mathbf{1}^T\lambda)}{\Gamma(\lambda_1)\cdots\Gamma(\lambda_{n+1})}
x_1^{\lambda_1-1}\cdots x_n^{\lambda_n-1}
\left(1-\sum_{i=1}^n x_i\right)^{\lambda_{n+1}-1},
$$

其中 $\operatorname{dom}f=\{x\in\mathbf{R}^n_{++}\mid \mathbf{1}^Tx<1\}$。参数 $\lambda$ 满足 $\lambda\succeq\mathbf{1}$。

3.57 证明函数 $f(X)=X^{-1}$ 在 $\mathbf{S}^n_{++}$ 上是矩阵凸的。

3.58 Schur 补。假设 $X\in\mathbf{S}^n$ 划分为 $X=\begin{bmatrix}A&B\\B^T&C\end{bmatrix}$，其中 $A\in\mathbf{S}^k$。$X$ 的 Schur 补（相对于 $A$）是 $S=C-B^TA^{-1}B$（参见第 A.5.5 节）。证明 Schur 补（被视为从 $\mathbf{S}^n$ 到 $\mathbf{S}^{n-k}$ 的函数）在 $\mathbf{S}^n_{++}$ 上是矩阵凹的。

3.59 K 凸性的二阶条件。令 $K\subseteq\mathbf{R}^m$ 为真凸锥，并具有相关的广义不等式 $\preceq_K$。证明具有凸域的二次可微函数 $f:\mathbf{R}^n\to\mathbf{R}^m$ 是 $K$-凸的，当且仅当对所有 $x\in\operatorname{dom}f$ 和所有 $y\in\mathbf{R}^n$，

$$
\sum_{i,j=1}^n \frac{\partial^2 f(x)}{\partial x_i\partial x_j}y_i y_j \succeq_K 0,
$$

也就是说，二阶导数是 $K$-非负双线性形式。（这里 $\partial^2 f/\partial x_i\partial x_j\in\mathbf{R}^m$，其分量为 $\partial^2 f_k/\partial x_i\partial x_j$，$k=1,\ldots,m$；参见 §A.4.1。）


3.60 K-凸函数的下水平集和上图。令 $K\subseteq\mathbf{R}^m$ 为真凸锥，并具有广义不等式 $\preceq_K$，并令 $f:\mathbf{R}^n\to\mathbf{R}^m$。对于 $\alpha\in\mathbf{R}^m$，$f$ 的 $\alpha$ 下水平集（相对于 $\preceq_K$）定义为 $C_\alpha=\{x\in\mathbf{R}^n\mid f(x)\preceq_K\alpha\}$。$f$ 相对于 $\preceq_K$ 的上图定义为集合 $\operatorname{epi}_K f=\{(x,t)\in\mathbf{R}^n\times\mathbf{R}^m\mid f(x)\preceq_K t\}$。证明以下内容： (a) 如果 $f$ 是 $K$-凸的，则它的下水平集 $C_\alpha$ 对所有 $\alpha$ 都是凸的。 (b) $f$ 是 $K$-凸的，当且仅当 $\operatorname{epi}_K f$ 是凸集。
