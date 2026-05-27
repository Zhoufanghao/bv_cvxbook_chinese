# 第 2 章 Convex sets

## 2.1 仿射集与凸集

### 2.1.1 直线与线段

设 $x_1\neq x_2$ 是 $\mathbf{R}^n$ 中的两个点。形如

$$
y=\theta x_1+(1-\theta)x_2
$$

的点，其中 $\theta\in\mathbf{R}$，构成经过 $x_1$ 和 $x_2$ 的直线。参数值 $\theta=0$ 对应 $y=x_2$，$\theta=1$ 对应 $y=x_1$。介于 0 和 1 之间的 $\theta$ 对应 $x_1$ 与 $x_2$ 之间的（闭）线段。

把 $y$ 写成

$$
y=x_2+\theta(x_1-x_2)
$$

给出了另一个解释：$y$ 是基点 $x_2$（对应 $\theta=0$）加上方向 $x_1-x_2$（从 $x_2$ 指向 $x_1$）按参数 $\theta$ 缩放后的结果。因此，$\theta$ 给出了 $y$ 位于从 $x_2$ 到 $x_1$ 路径上的比例。随着 $\theta$ 从 0 增加到 1，点 $y$ 从 $x_2$ 移动到 $x_1$；当 $\theta>1$ 时，点 $y$ 位于越过 $x_1$ 的直线上。图 2.1 展示了这一点。

**图 2.1** 经过 $x_1$ 和 $x_2$ 的直线可由 $\theta x_1+(1-\theta)x_2$ 参数化描述，其中 $\theta$ 在 $\mathbf{R}$ 中变化。$x_1$ 与 $x_2$ 之间的线段对应 $0\leq\theta\leq1$。

### 2.1.2 仿射集

若对集合 $C\subseteq\mathbf{R}^n$ 中任意两个不同点，经过它们的整条直线都包含在 $C$ 中，则称 $C$ 是仿射集；也就是说，对任意 $x_1,x_2\in C$ 和 $\theta\in\mathbf{R}$，都有 $\theta x_1+(1-\theta)x_2\in C$。换言之，只要线性组合中的系数和为 1，$C$ 就包含其任意两个点的这种线性组合。

这个思想可以推广到两个以上的点。形如 $\theta_1x_1+\cdots+\theta_kx_k$ 且 $\theta_1+\cdots+\theta_k=1$ 的点，称为点 $x_1,\ldots,x_k$ 的仿射组合。由仿射集的定义作归纳可知，仿射集包含其点的所有仿射组合：若 $C$ 是仿射集，$x_1,\ldots,x_k\in C$，并且 $\theta_1+\cdots+\theta_k=1$，则 $\theta_1x_1+\cdots+\theta_kx_k\in C$。

若 $C$ 是仿射集且 $x_0\in C$，则集合

$$
V=C-x_0=\{x-x_0\mid x\in C\}
$$

是一个子空间，即对加法和标量乘法封闭。为说明这一点，设 $v_1,v_2\in V$ 且 $\alpha,\beta\in\mathbf{R}$。于是 $v_1+x_0\in C$ 且 $v_2+x_0\in C$，并且

$$
\alpha v_1+\beta v_2+x_0
=\alpha(v_1+x_0)+\beta(v_2+x_0)+(1-\alpha-\beta)x_0\in C,
$$

因为 $C$ 是仿射的，且 $\alpha+\beta+(1-\alpha-\beta)=1$。由 $\alpha v_1+\beta v_2+x_0\in C$ 可得 $\alpha v_1+\beta v_2\in V$。

因此，仿射集 $C$ 可以表示为

$$
C=V+x_0=\{v+x_0\mid v\in V\},
$$

即一个子空间加上一个偏移。与仿射集 $C$ 关联的子空间 $V$ 不依赖于 $x_0$ 的选择，所以 $x_0$ 可以取为 $C$ 中任意一点。我们把仿射集 $C$ 的维数定义为子空间 $V=C-x_0$ 的维数，其中 $x_0$ 是 $C$ 的任意元素。

**例 2.1 线性方程组的解集。** 线性方程组的解集 $C=\{x\mid Ax=b\}$，其中 $A\in\mathbf{R}^{m\times n}$、$b\in\mathbf{R}^m$，是仿射集。事实上，若 $x_1,x_2\in C$，即 $Ax_1=b$、$Ax_2=b$，则对任意 $\theta$，

$$
A(\theta x_1+(1-\theta)x_2)
=\theta Ax_1+(1-\theta)Ax_2
=\theta b+(1-\theta)b=b.
$$

因此仿射组合 $\theta x_1+(1-\theta)x_2$ 也在 $C$ 中。与仿射集 $C$ 关联的子空间是 $A$ 的零空间。反过来也成立：每个仿射集都可以表示为某个线性方程组的解集。

集合 $C\subseteq\mathbf{R}^n$ 中所有点的仿射组合组成的集合称为 $C$ 的仿射包，记为 $\operatorname{aff} C$：

$$
\operatorname{aff} C
=\{\theta_1x_1+\cdots+\theta_kx_k\mid
x_1,\ldots,x_k\in C,\ \theta_1+\cdots+\theta_k=1\}.
$$

仿射包是包含 $C$ 的最小仿射集：若 $S$ 是任意满足 $C\subseteq S$ 的仿射集，则 $\operatorname{aff} C\subseteq S$。

### 2.1.3 仿射维数与相对内部

集合 $C$ 的仿射维数定义为其仿射包的维数。仿射维数在凸分析和优化中很有用，但并不总是与其他维数定义一致。例如，考虑 $\mathbf{R}^2$ 中的单位圆 $\{x\in\mathbf{R}^2\mid x_1^2+x_2^2=1\}$。它的仿射包是整个 $\mathbf{R}^2$，所以仿射维数为 2；但按大多数维数定义，$\mathbf{R}^2$ 中的单位圆维数为 1。

若集合 $C\subseteq\mathbf{R}^n$ 的仿射维数小于 $n$，则该集合位于仿射集 $\operatorname{aff}C\neq\mathbf{R}^n$ 中。我们定义集合 $C$ 的相对内部，记为 $\operatorname{relint} C$，为其相对于 $\operatorname{aff} C$ 的内部：

$$
\operatorname{relint} C
=\{x\in C\mid B(x,r)\cap\operatorname{aff}C\subseteq C\ \text{for some } r>0\},
$$

其中 $B(x,r)=\{y\mid \|y-x\|\leq r\}$ 是以 $x$ 为中心、半径为 $r$ 的球。（这里 $\|\cdot\|$ 可以是任意范数；所有范数给出相同的相对内部。）集合 $C$ 的相对边界定义为 $\operatorname{cl}C\setminus\operatorname{relint}C$，其中 $\operatorname{cl}C$ 是 $C$ 的闭包。

**例 2.2** 考虑 $\mathbf{R}^3$ 中 $(x_1,x_2)$ 平面上的一个正方形：

$$
C=\{x\in\mathbf{R}^3\mid -1\leq x_1\leq1,\ -1\leq x_2\leq1,\ x_3=0\}.
$$

它的仿射包是 $(x_1,x_2)$ 平面，即 $\operatorname{aff}C=\{x\in\mathbf{R}^3\mid x_3=0\}$。$C$ 的内部为空，但相对内部为

$$
\operatorname{relint}C
=\{x\in\mathbf{R}^3\mid -1<x_1<1,\ -1<x_2<1,\ x_3=0\}.
$$

它在 $\mathbf{R}^3$ 中的边界是它自身；它的相对边界是线框轮廓

$$
\{x\in\mathbf{R}^3\mid \max\{|x_1|,|x_2|\}=1,\ x_3=0\}.
$$

### 2.1.4 凸集

若集合 $C$ 中任意两点之间的线段都包含在 $C$ 中，则称 $C$ 是凸集；也就是说，对任意 $x_1,x_2\in C$ 和任意 $0\leq\theta\leq1$，都有

$$
\theta x_1+(1-\theta)x_2\in C.
$$

粗略地说，如果集合中任意一点都能沿着一条无阻挡的直线路径“看到”集合中的任意另一点，并且“无阻挡”表示路径位于集合中，则该集合是凸的。每个仿射集也是凸的，因为它包含任意两点之间的整条直线，因此也包含这两点之间的线段。图 2.2 给出了 $\mathbf{R}^2$ 中一些简单的凸集和非凸集。

**图 2.2** 一些简单的凸集和非凸集。左：包含边界的六边形是凸的。中：肾形集合不是凸的，因为图中两个点之间的线段不全在集合内。右：正方形包含部分边界点但不包含另一些边界点，因此不是凸的。

若 $\theta_1x_1+\cdots+\theta_kx_k$ 满足 $\theta_1+\cdots+\theta_k=1$ 且 $\theta_i\geq0$，$i=1,\ldots,k$，则称其为点 $x_1,\ldots,x_k$ 的凸组合。和仿射集一样，可以证明一个集合是凸的，当且仅当它包含其点的所有凸组合。点的凸组合可看作这些点的混合或加权平均，其中 $\theta_i$ 是混合中 $x_i$ 的比例。

集合 $C$ 的凸包记为 $\operatorname{conv}C$，是 $C$ 中点的所有凸组合的集合：

$$
\operatorname{conv}C
=\{\theta_1x_1+\cdots+\theta_kx_k\mid
x_i\in C,\ \theta_i\geq0,\ i=1,\ldots,k,\ \theta_1+\cdots+\theta_k=1\}.
$$

顾名思义，凸包 $\operatorname{conv}C$ 总是凸的。它是包含 $C$ 的最小凸集：若 $B$ 是任意包含 $C$ 的凸集，则 $\operatorname{conv}C\subseteq B$。图 2.3 展示了凸包的定义。

**图 2.3** $\mathbf{R}^2$ 中两个集合的凸包。左：15 个点的凸包是阴影五边形。右：图 2.2 中肾形集合的凸包是阴影集合。

凸组合的思想可以推广到无穷求和、积分，以及最一般形式的概率分布。假设 $\theta_1,\theta_2,\ldots$ 满足

$$
\theta_i\geq0,\quad i=1,2,\ldots,\qquad
\sum_{i=1}^{\infty}\theta_i=1,
$$

并且 $x_1,x_2,\ldots\in C$，其中 $C\subseteq\mathbf{R}^n$ 是凸集。若级数收敛，则

$$
\sum_{i=1}^{\infty}\theta_i x_i\in C.
$$

更一般地，设 $p:\mathbf{R}^n\to\mathbf{R}$ 满足对所有 $x\in C$ 都有 $p(x)\geq0$，并且 $\int_Cp(x)\,dx=1$，其中 $C\subseteq\mathbf{R}^n$ 是凸集。若积分存在，则

$$
\int_C p(x)x\,dx\in C.
$$

最一般地，设 $C\subseteq\mathbf{R}^n$ 是凸集，随机向量 $x$ 以概率 1 属于 $C$。则 $\mathbf{E}x\in C$。事实上，这一形式包含前述所有情形。例如，若随机变量 $x$ 只取两个值 $x_1$ 和 $x_2$，且 $\operatorname{prob}(x=x_1)=\theta$、$\operatorname{prob}(x=x_2)=1-\theta$，其中 $0\leq\theta\leq1$，则 $\mathbf{E}x=\theta x_1+(1-\theta)x_2$，这正是两个点的简单凸组合。

### 2.1.5 锥

若对每个 $x\in C$ 和 $\theta\geq0$ 都有 $\theta x\in C$，则称集合 $C$ 是锥，或称它非负齐次。若集合 $C$ 既是凸集又是锥，则称其为凸锥；这意味着对任意 $x_1,x_2\in C$ 和 $\theta_1,\theta_2\geq0$，都有

$$
\theta_1x_1+\theta_2x_2\in C.
$$

从几何上看，这类点形成以 0 为顶点、边通过 $x_1$ 和 $x_2$ 的二维扇形区域。（见图 2.4。）

**图 2.4** 扇形区域表示所有形如 $\theta_1x_1+\theta_2x_2$ 的点，其中 $\theta_1,\theta_2\geq0$。顶点对应 $\theta_1=\theta_2=0$，位于 0；两条边对应 $\theta_1=0$ 或 $\theta_2=0$，分别通过 $x_1$ 和 $x_2$。

形如 $\theta_1x_1+\cdots+\theta_kx_k$ 且 $\theta_1,\ldots,\theta_k\geq0$ 的点，称为 $x_1,\ldots,x_k$ 的锥组合（或非负线性组合）。若 $x_i$ 在凸锥 $C$ 中，则 $x_i$ 的任意锥组合也在 $C$ 中。反过来，集合 $C$ 是凸锥，当且仅当它包含其元素的所有锥组合。与凸组合（或仿射组合）类似，锥组合的思想也可以推广到无穷求和和积分。

集合 $C$ 的锥包是 $C$ 中点的所有锥组合的集合，即

$$
\{\theta_1x_1+\cdots+\theta_kx_k\mid
x_i\in C,\ \theta_i\geq0,\ i=1,\ldots,k\},
$$

它也是包含 $C$ 的最小凸锥（见图 2.5）。

**图 2.5** 图 2.3 中两个集合的锥包。

## 2.2 一些重要例子

本节描述一些重要的凸集例子，这些集合将在本书其余部分反复出现。先从若干简单例子开始。

- 空集 $\emptyset$、任意单点集 $\{x_0\}$、以及整个空间 $\mathbf{R}^n$，都是 $\mathbf{R}^n$ 的仿射子集（因而也是凸集）。
- 任意直线都是仿射的。若它经过原点，则它是子空间，因此也是凸锥。
- 线段是凸的，但不是仿射的（除非它退化为一个点）。
- 射线形如 $\{x_0+\theta v\mid \theta\geq0\}$，其中 $v\neq0$。射线是凸的，但不是仿射的；若其基点 $x_0$ 为 0，则它是凸锥。
- 任意子空间都是仿射的，并且是凸锥（因此也是凸集）。

### 2.2.1 超平面与半空间

超平面是形如

$$
\{x\mid a^Tx=b\}
$$

的集合，其中 $a\in\mathbf{R}^n$、$a\neq0$，$b\in\mathbf{R}$。从解析角度看，它是 $x$ 的分量之间一个非平凡线性方程的解集（因此是仿射集）。从几何角度看，超平面 $\{x\mid a^Tx=b\}$ 可以解释为与给定向量 $a$ 的内积为常数的点集，或解释为法向量为 $a$ 的超平面；常数 $b\in\mathbf{R}$ 决定该超平面相对于原点的偏移。这个几何解释可通过把超平面写成

$$
\{x\mid a^T(x-x_0)=0\}
$$

来理解，其中 $x_0$ 是超平面上的任意点，即满足 $a^Tx_0=b$ 的任意点。该表示又可写为

$$
\{x\mid a^T(x-x_0)=0\}=x_0+a^\perp,
$$

其中 $a^\perp$ 表示 $a$ 的正交补，即所有与 $a$ 正交的向量集合：

$$
a^\perp=\{v\mid a^Tv=0\}.
$$

这说明超平面由偏移 $x_0$ 加上所有与法向量 $a$ 正交的向量组成。图 2.6 展示了这些几何解释。

**图 2.6** $\mathbf{R}^2$ 中的超平面，具有法向量 $a$ 和超平面上一点 $x_0$。对超平面上的任意点 $x$，$x-x_0$ 与 $a$ 正交。

一个超平面把 $\mathbf{R}^n$ 分成两个半空间。（闭）半空间是形如

$$
\{x\mid a^Tx\leq b\},
\tag{2.1}
$$

的集合，其中 $a\neq0$，即一个非平凡线性不等式的解集。半空间是凸的，但不是仿射的。图 2.7 展示了这一点。

**图 2.7** $\mathbf{R}^2$ 中由 $a^Tx=b$ 定义的超平面决定两个半空间。$a^Tx\geq b$ 决定的半空间沿 $a$ 的方向延伸；$a^Tx\leq b$ 决定的半空间沿 $-a$ 的方向延伸。向量 $a$ 是后者的外法向量。

半空间 (2.1) 也可以表示为

$$
\{x\mid a^T(x-x_0)\leq0\},
\tag{2.2}
$$

其中 $x_0$ 是关联超平面上的任意点，即满足 $a^Tx_0=b$。表示 (2.2) 给出一个简单几何解释：半空间由 $x_0$ 加上所有与（外法向量）$a$ 成钝角或直角的向量组成。图 2.8 展示了这一点。

**图 2.8** 阴影集合是由 $a^T(x-x_0)\leq0$ 决定的半空间。向量 $x_1-x_0$ 与 $a$ 成锐角，所以 $x_1$ 不在半空间中；向量 $x_2-x_0$ 与 $a$ 成钝角，所以 $x_2$ 在半空间中。

半空间 (2.1) 的边界是超平面 $\{x\mid a^Tx=b\}$。集合 $\{x\mid a^Tx<b\}$ 是半空间 $\{x\mid a^Tx\leq b\}$ 的内部，称为开半空间。

### 2.2.2 欧几里得球与椭球

$\mathbf{R}^n$ 中的（欧几里得）球，或简称球，具有形式

$$
B(x_c,r)=\{x\mid \|x-x_c\|_2\leq r\}
=\{x\mid (x-x_c)^T(x-x_c)\leq r^2\},
$$

其中 $r>0$，$\|\cdot\|_2$ 表示欧几里得范数，即 $\|u\|_2=(u^Tu)^{1/2}$。向量 $x_c$ 是球心，标量 $r$ 是半径；$B(x_c,r)$ 由所有到球心 $x_c$ 距离不超过 $r$ 的点组成。欧几里得球的另一种常用表示为

$$
B(x_c,r)=\{x_c+ru\mid \|u\|_2\leq1\}.
$$

欧几里得球是凸集：若 $\|x_1-x_c\|_2\leq r$、$\|x_2-x_c\|_2\leq r$ 且 $0\leq\theta\leq1$，则

$$
\begin{aligned}
\|\theta x_1+(1-\theta)x_2-x_c\|_2
&=\|\theta(x_1-x_c)+(1-\theta)(x_2-x_c)\|_2\\
&\leq \theta\|x_1-x_c\|_2+(1-\theta)\|x_2-x_c\|_2\\
&\leq r.
\end{aligned}
$$

这里使用了 $\|\cdot\|_2$ 的齐次性和三角不等式；见 §A.1.2。

一类相关的凸集是椭球，其形式为

$$
\mathcal{E}=\{x\mid (x-x_c)^TP^{-1}(x-x_c)\leq1\},
\tag{2.3}
$$

其中 $P=P^T\succ0$，即 $P$ 是对称正定矩阵。向量 $x_c\in\mathbf{R}^n$ 是椭球中心。矩阵 $P$ 决定椭球从 $x_c$ 向各方向延伸的距离；$\mathcal{E}$ 的半轴长度为 $\sqrt{\lambda_i}$，其中 $\lambda_i$ 是 $P$ 的特征值。球是满足 $P=r^2I$ 的椭球。图 2.9 展示了 $\mathbf{R}^2$ 中的一个椭球。

**图 2.9** $\mathbf{R}^2$ 中的椭球。点表示中心 $x_c$，线段表示两个半轴。

椭球的另一种常用表示为

$$
\mathcal{E}=\{x_c+Au\mid \|u\|_2\leq1\},
\tag{2.4}
$$

其中 $A$ 是方阵且非奇异。在这种表示中，可以不失一般性地假设 $A$ 对称正定。取 $A=P^{1/2}$ 就得到 (2.3) 中定义的椭球。当 (2.4) 中的矩阵 $A$ 对称正半定但奇异时，集合 (2.4) 称为退化椭球；它的仿射维数等于 $A$ 的秩。退化椭球也是凸的。

### 2.2.3 范数球与范数锥

设 $\|\cdot\|$ 是 $\mathbf{R}^n$ 上的任意范数（见 §A.1.2）。由范数的一般性质可知，半径为 $r$、中心为 $x_c$ 的范数球 $\{x\mid \|x-x_c\|\leq r\}$ 是凸的。与范数 $\|\cdot\|$ 关联的范数锥是集合

$$
C=\{(x,t)\mid \|x\|\leq t\}\subseteq\mathbf{R}^{n+1}.
$$

它（顾名思义）是一个凸锥。

**例 2.3** 二阶锥是欧几里得范数的范数锥，即

$$
\begin{aligned}
C
&=\{(x,t)\in\mathbf{R}^{n+1}\mid \|x\|_2\leq t\}\\
&=\left\{(x,t)\ \middle|\ 
\begin{bmatrix}x\\ t\end{bmatrix}^T
\begin{bmatrix}I&0\\0&-1\end{bmatrix}
\begin{bmatrix}x\\ t\end{bmatrix}\leq0,\ t\geq0
\right\}.
\end{aligned}
$$

二阶锥还有若干其他名称。由于它由二次不等式定义，也称为二次锥；它也称为 Lorentz 锥或冰激凌锥。图 2.10 展示了 $\mathbf{R}^3$ 中的二阶锥。

**图 2.10** $\mathbf{R}^3$ 中二阶锥 $\{(x_1,x_2,t)\mid (x_1^2+x_2^2)^{1/2}\leq t\}$ 的边界。

### 2.2.4 多面体

多面体定义为有限个线性等式和不等式的解集：

$$
P=\{x\mid a_j^Tx\leq b_j,\ j=1,\ldots,m,\quad c_j^Tx=d_j,\ j=1,\ldots,p\}.
\tag{2.5}
$$

因此，多面体是有限个半空间和超平面的交。仿射集（例如子空间、超平面、直线）、射线、线段和半空间都是多面体。很容易证明多面体是凸集。有界多面体有时称为多胞形，但有些作者采用相反约定（即把形如 (2.5) 的任意集合称为多胞形，而有界时称为多面体）。图 2.11 展示了由五个半空间交成的多面体。

**图 2.11** 多面体 $P$ 是五个半空间的交，外法向量为 $a_1,\ldots,a_5$。

为方便起见，我们对 (2.5) 使用紧凑记号

$$
P=\{x\mid Ax\preceq b,\ Cx=d\},
\tag{2.6}
$$

其中

$$
A=\begin{bmatrix}a_1^T\\ \vdots\\ a_m^T\end{bmatrix},\qquad
C=\begin{bmatrix}c_1^T\\ \vdots\\ c_p^T\end{bmatrix},
$$

符号 $\preceq$ 表示 $\mathbf{R}^m$ 中的向量不等式或逐分量不等式：$u\preceq v$ 表示 $u_i\leq v_i$，$i=1,\ldots,m$。

**例 2.4** 非负正交锥是所有分量非负的点集合，即

$$
\mathbf{R}^n_+=\{x\in\mathbf{R}^n\mid x_i\geq0,\ i=1,\ldots,n\}
=\{x\in\mathbf{R}^n\mid x\succeq0\}.
$$

这里 $\mathbf{R}_+$ 表示非负实数集合：$\mathbf{R}_+=\{x\in\mathbf{R}\mid x\geq0\}$。非负正交锥是多面体且是锥，因此称为多面锥。

#### 单纯形

单纯形是另一类重要的多面体。设 $k+1$ 个点 $v_0,\ldots,v_k\in\mathbf{R}^n$ 仿射无关，也就是说 $v_1-v_0,\ldots,v_k-v_0$ 线性无关。由它们决定的单纯形为

$$
C=\operatorname{conv}\{v_0,\ldots,v_k\}
=\{\theta_0v_0+\cdots+\theta_kv_k\mid \theta\succeq0,\ \mathbf{1}^T\theta=1\}.
\tag{2.7}
$$

其中 $\mathbf{1}$ 表示所有分量均为 1 的向量。该单纯形的仿射维数为 $k$，所以有时称为 $\mathbf{R}^n$ 中的 $k$ 维单纯形。

**例 2.5 一些常见单纯形。** 1 维单纯形是线段；2 维单纯形是三角形（包括内部）；3 维单纯形是四面体。

单位单纯形是由零向量和单位向量 $0,e_1,\ldots,e_n\in\mathbf{R}^n$ 决定的 $n$ 维单纯形。它可以表示为满足

$$
x\succeq0,\qquad \mathbf{1}^Tx\leq1
$$

的向量集合。

概率单纯形是由单位向量 $e_1,\ldots,e_n\in\mathbf{R}^n$ 决定的 $(n-1)$ 维单纯形。它是满足

$$
x\succeq0,\qquad \mathbf{1}^Tx=1
$$

的向量集合。概率单纯形中的向量对应于具有 $n$ 个元素的集合上的概率分布，其中 $x_i$ 解释为第 $i$ 个元素的概率。

为了把单纯形 (2.7) 描述为形如 (2.6) 的多面体，可作如下处理。按定义，$x\in C$ 当且仅当存在 $\theta\succeq0$ 且 $\mathbf{1}^T\theta=1$，使得 $x=\theta_0v_0+\theta_1v_1+\cdots+\theta_kv_k$。等价地，若定义 $y=(\theta_1,\ldots,\theta_k)$ 以及

$$
B=[\,v_1-v_0\ \cdots\ v_k-v_0\,]\in\mathbf{R}^{n\times k},
$$

则 $x\in C$ 当且仅当

$$
x=v_0+By
\tag{2.8}
$$

对某个满足 $y\succeq0$ 且 $\mathbf{1}^Ty\leq1$ 的 $y$ 成立。点 $v_0,\ldots,v_k$ 的仿射无关性意味着矩阵 $B$ 的秩为 $k$。因此存在非奇异矩阵 $A=(A_1,A_2)\in\mathbf{R}^{n\times n}$，使得

$$
\begin{bmatrix}A_1\\A_2\end{bmatrix}B
=\begin{bmatrix}I\\0\end{bmatrix}.
$$

将 (2.8) 左乘 $A$，得到

$$
A_1x=A_1v_0+y,\qquad A_2x=A_2v_0.
$$

因此，$x\in C$ 当且仅当 $A_2x=A_2v_0$，并且向量 $y=A_1x-A_1v_0$ 满足 $y\succeq0$ 和 $\mathbf{1}^Ty\leq1$。换言之，$x\in C$ 当且仅当

$$
A_2x=A_2v_0,\qquad
A_1x\succeq A_1v_0,\qquad
\mathbf{1}^TA_1x\leq1+\mathbf{1}^TA_1v_0,
$$

这是关于 $x$ 的一组线性等式和不等式，因此描述了一个多面体。

#### 多面体的凸包描述

有限集 $\{v_1,\ldots,v_k\}$ 的凸包是

$$
\operatorname{conv}\{v_1,\ldots,v_k\}
=\{\theta_1v_1+\cdots+\theta_kv_k\mid
\theta\succeq0,\ \mathbf{1}^T\theta=1\}.
$$

这个集合是多面体且有界，但除特殊情形（例如单纯形）外，要把它表示为 (2.5) 的形式，即用一组线性等式和不等式描述，并不简单。

这一凸包描述的一个推广是

$$
\{\theta_1v_1+\cdots+\theta_kv_k\mid
\theta_1+\cdots+\theta_m=1,\ \theta_i\geq0,\ i=1,\ldots,k\},
\tag{2.9}
$$

其中 $m\leq k$。这里考虑 $v_i$ 的非负线性组合，但只要求前 $m$ 个系数之和为 1。也可以把 (2.9) 解释为点 $v_1,\ldots,v_m$ 的凸包，加上点 $v_{m+1},\ldots,v_k$ 的锥包。集合 (2.9) 定义了一个多面体；反过来，每个多面体都可以用这种形式表示（虽然这里不证明）。

多面体如何表示是一个微妙问题，并具有非常重要的实际后果。作为简单例子，考虑 $\mathbf{R}^n$ 中 $\ell_\infty$ 范数的单位球

$$
C=\{x\mid |x_i|\leq1,\ i=1,\ldots,n\}.
$$

集合 $C$ 可以用 $2n$ 个线性不等式 $\pm e_i^Tx\leq1$ 表示为 (2.5) 的形式，其中 $e_i$ 是第 $i$ 个单位向量。若要用凸包形式 (2.9) 描述它，则至少需要 $2^n$ 个点：

$$
C=\operatorname{conv}\{v_1,\ldots,v_{2^n}\},
$$

其中 $v_1,\ldots,v_{2^n}$ 是所有分量均为 1 或 $-1$ 的 $2^n$ 个向量。因此当 $n$ 很大时，两种描述的规模差别巨大。

### 2.2.5 正半定锥

我们用 $\mathbf{S}^n$ 表示对称 $n\times n$ 矩阵集合：

$$
\mathbf{S}^n=\{X\in\mathbf{R}^{n\times n}\mid X=X^T\},
$$

它是维数为 $n(n+1)/2$ 的向量空间。用 $\mathbf{S}^n_+$ 表示对称正半定矩阵集合：

$$
\mathbf{S}^n_+=\{X\in\mathbf{S}^n\mid X\succeq0\},
$$

用 $\mathbf{S}^n_{++}$ 表示对称正定矩阵集合：

$$
\mathbf{S}^n_{++}=\{X\in\mathbf{S}^n\mid X\succ0\}.
$$

这个记号意在类比 $\mathbf{R}_+$ 和 $\mathbf{R}_{++}$，前者表示非负实数，后者表示正实数。

集合 $\mathbf{S}^n_+$ 是凸锥：若 $\theta_1,\theta_2\geq0$ 且 $A,B\in\mathbf{S}^n_+$，则 $\theta_1A+\theta_2B\in\mathbf{S}^n_+$。这可以直接由正半定性的定义看出：对任意 $x\in\mathbf{R}^n$，

$$
x^T(\theta_1A+\theta_2B)x
=\theta_1x^TAx+\theta_2x^TBx\geq0,
$$

只要 $A\succeq0$、$B\succeq0$ 且 $\theta_1,\theta_2\geq0$。

**例 2.6 $\mathbf{S}^2$ 中的正半定锥。** 有

$$
X=\begin{bmatrix}x&y\\y&z\end{bmatrix}\in\mathbf{S}^2_+
\quad\Longleftrightarrow\quad
x\geq0,\ z\geq0,\ xz\geq y^2.
$$

图 2.12 以 $(x,y,z)$ 作为 $\mathbf{R}^3$ 中坐标，展示了这个锥的边界。

**图 2.12** $\mathbf{S}^2$ 中正半定锥的边界。

## 2.3 保持凸性的运算

本节描述一些保持集合凸性的运算，或可用于从已有集合构造凸集的运算。这些运算与 §2.2 中的简单例子一起，形成一种凸集演算，可用于判断或证明集合的凸性。

### 2.3.1 交

凸性在取交下保持：若 $S_1$ 和 $S_2$ 是凸的，则 $S_1\cap S_2$ 是凸的。这个性质可推广到无穷多个集合的交：若对每个 $\alpha\in A$，$S_\alpha$ 都是凸的，则 $\bigcap_{\alpha\in A}S_\alpha$ 是凸的。（子空间、仿射集和凸锥也对任意交封闭。）作为简单例子，多面体是半空间和超平面（它们都是凸的）的交，因此是凸的。

**例 2.7** 正半定锥 $\mathbf{S}^n_+$ 可以表示为

$$
\bigcap_{z\neq0}\{X\in\mathbf{S}^n\mid z^TXz\geq0\}.
$$

对每个 $z\neq0$，$z^TXz$ 是 $X$ 的一个非零线性函数，所以集合

$$
\{X\in\mathbf{S}^n\mid z^TXz\geq0\}
$$

实际上是 $\mathbf{S}^n$ 中的半空间。因此，正半定锥是无穷多个半空间的交，所以是凸的。

**例 2.8** 考虑集合

$$
S=\{x\in\mathbf{R}^m\mid |p(t)|\leq1\ \text{for } |t|\leq\pi/3\},
\tag{2.10}
$$

其中 $p(t)=\sum_{k=1}^m x_k\cos kt$。集合 $S$ 可以表示为无穷多个板的交：$S=\bigcap_{|t|\leq\pi/3}S_t$，其中

$$
S_t=\{x\mid -1\leq(\cos t,\ldots,\cos mt)^Tx\leq1\},
$$

因此是凸的。图 2.13 和图 2.14 以 $m=2$ 为例展示了定义和集合。

**图 2.13** 与集合 (2.10) 中点对应的三个三角多项式，$m=2$。虚线绘制的三角多项式是另外两个的平均。

**图 2.14** 集合 (2.10) 在 $m=2$ 时的图像。该集合是无穷多个板的交（图中显示了其中 20 个），因此是凸的。

在以上例子中，我们通过把集合表示为半空间的（可能无穷）交来证明其凸性。§2.5.1 将看到一个反向结论：每个闭凸集 $S$ 都是半空间的（通常无穷）交。事实上，一个闭凸集 $S$ 是所有包含它的半空间的交：

$$
S=\bigcap\{H\mid H\text{ is a halfspace},\ S\subseteq H\}.
$$

### 2.3.2 仿射函数

回忆若函数 $f:\mathbf{R}^n\to\mathbf{R}^m$ 可写为线性函数与常数之和，即 $f(x)=Ax+b$，其中 $A\in\mathbf{R}^{m\times n}$、$b\in\mathbf{R}^m$，则称 $f$ 是仿射函数。设 $S\subseteq\mathbf{R}^n$ 是凸集，$f:\mathbf{R}^n\to\mathbf{R}^m$ 是仿射函数。则 $S$ 在 $f$ 下的像

$$
f(S)=\{f(x)\mid x\in S\}
$$

是凸的。类似地，若 $f:\mathbf{R}^k\to\mathbf{R}^n$ 是仿射函数，则 $S$ 在 $f$ 下的逆像

$$
f^{-1}(S)=\{x\mid f(x)\in S\}
$$

也是凸的。

两个简单例子是缩放和平移。若 $S\subseteq\mathbf{R}^n$ 是凸的，$\alpha\in\mathbf{R}$，$a\in\mathbf{R}^n$，则集合 $\alpha S$ 和 $S+a$ 是凸的，其中

$$
\alpha S=\{\alpha x\mid x\in S\},\qquad
S+a=\{x+a\mid x\in S\}.
$$

凸集投影到某些坐标上仍然是凸的：若 $S\subseteq\mathbf{R}^m\times\mathbf{R}^n$ 是凸的，则

$$
T=\{x_1\in\mathbf{R}^m\mid (x_1,x_2)\in S\ \text{for some }x_2\in\mathbf{R}^n\}
$$

是凸的。

两个集合的和定义为

$$
S_1+S_2=\{x+y\mid x\in S_1,\ y\in S_2\}.
$$

若 $S_1$ 和 $S_2$ 是凸的，则 $S_1+S_2$ 是凸的。理由是：若 $S_1$ 和 $S_2$ 是凸的，则其直积

$$
S_1\times S_2=\{(x_1,x_2)\mid x_1\in S_1,\ x_2\in S_2\}
$$

也是凸的；该集合在线性函数 $f(x_1,x_2)=x_1+x_2$ 下的像就是 $S_1+S_2$。

还可以考虑 $S_1,S_2\subseteq\mathbf{R}^n\times\mathbf{R}^m$ 的部分和，定义为

$$
S=\{(x,y_1+y_2)\mid (x,y_1)\in S_1,\ (x,y_2)\in S_2\},
$$

其中 $x\in\mathbf{R}^n$、$y_i\in\mathbf{R}^m$。当 $m=0$ 时，部分和给出 $S_1$ 与 $S_2$ 的交；当 $n=0$ 时，它就是集合加法。凸集的部分和是凸的（见习题 2.16）。

**例 2.9 多面体。** 多面体 $\{x\mid Ax\preceq b,\ Cx=d\}$ 可以表示为非负正交锥与原点的笛卡尔积在仿射函数 $f(x)=(b-Ax,d-Cx)$ 下的逆像：

$$
\{x\mid Ax\preceq b,\ Cx=d\}=\{x\mid f(x)\in\mathbf{R}_+^m\times\{0\}\}.
$$

**例 2.10 线性矩阵不等式的解集。** 条件

$$
A(x)=x_1A_1+\cdots+x_nA_n\preceq B,
\tag{2.11}
$$

其中 $B,A_i\in\mathbf{S}^m$，称为关于 $x$ 的线性矩阵不等式（LMI）。注意它与普通线性不等式

$$
a^Tx=x_1a_1+\cdots+x_na_n\leq b
$$

的相似性，其中 $b,a_i\in\mathbf{R}$。

线性矩阵不等式的解集 $\{x\mid A(x)\preceq B\}$ 是凸的。事实上，它是正半定锥在仿射函数 $f:\mathbf{R}^n\to\mathbf{S}^m$，$f(x)=B-A(x)$ 下的逆像。

**例 2.11 双曲锥。** 集合

$$
\{x\mid x^TPx\leq(c^Tx)^2,\ c^Tx\geq0\},
$$

其中 $P\in\mathbf{S}^n_+$、$c\in\mathbf{R}^n$，是凸的，因为它是二阶锥

$$
\{(z,t)\mid z^Tz\leq t^2,\ t\geq0\}
$$

在仿射函数 $f(x)=(P^{1/2}x,c^Tx)$ 下的逆像。

**例 2.12 椭球。** 椭球

$$
\mathcal{E}=\{x\mid (x-x_c)^TP^{-1}(x-x_c)\leq1\},
$$

其中 $P\in\mathbf{S}^n_{++}$，是单位欧几里得球 $\{u\mid \|u\|_2\leq1\}$ 在仿射映射 $f(u)=P^{1/2}u+x_c$ 下的像。（它也是单位球在仿射映射 $g(x)=P^{-1/2}(x-x_c)$ 下的逆像。）

### 2.3.3 线性分式函数与透视函数

本节考察一类称为线性分式函数的函数。它比仿射函数更一般，但仍保持凸性。

#### 透视函数

定义透视函数 $P:\mathbf{R}^{n+1}\to\mathbf{R}^n$，其定义域为 $\operatorname{dom}P=\mathbf{R}^n\times\mathbf{R}_{++}$，且 $P(z,t)=z/t$。（这里 $\mathbf{R}_{++}$ 表示正数集合：$\mathbf{R}_{++}=\{x\in\mathbf{R}\mid x>0\}$。）透视函数对向量进行缩放或归一化，使最后一个分量为 1，然后去掉最后一个分量。

**注 2.1** 可以把透视函数解释为针孔相机的作用。$\mathbf{R}^3$ 中的针孔相机由一个不透明水平平面 $x_3=0$ 构成，该平面在原点有唯一针孔，光可以通过；另有水平像平面 $x_3=-1$。位于相机上方（即 $x_3>0$）的物体在像平面上形成点 $-(x_1/x_3,x_2/x_3,1)$。去掉图像点的最后一个分量（因为它总是 $-1$），位于 $x$ 的点在像平面上的图像为 $y=-(x_1/x_3,x_2/x_3)=-P(x)$。图 2.15 展示了这一点。

**图 2.15** 透视函数的针孔相机解释。深色水平线表示 $\mathbf{R}^3$ 中的平面 $x_3=0$，除原点针孔外不透光。平面上方的物体或光源显示在像平面 $x_3=-1$ 上。光源位置到图像位置的映射与透视函数相关。

若 $C\subseteq\operatorname{dom}P$ 是凸的，则其像

$$
P(C)=\{P(x)\mid x\in C\}
$$

也是凸的。这一结果很直观：通过针孔相机观察一个凸物体，会得到一个凸图像。为证明这一点，我们说明透视函数把线段映射为线段。设 $x=(\tilde{x},x_{n+1})$、$y=(\tilde{y},y_{n+1})\in\mathbf{R}^{n+1}$，且 $x_{n+1}>0$、$y_{n+1}>0$。则对 $0\leq\theta\leq1$，

$$
P(\theta x+(1-\theta)y)
=\frac{\theta\tilde{x}+(1-\theta)\tilde{y}}{\theta x_{n+1}+(1-\theta)y_{n+1}}
=\mu P(x)+(1-\mu)P(y),
$$

其中

$$
\mu=\frac{\theta x_{n+1}}{\theta x_{n+1}+(1-\theta)y_{n+1}}\in[0,1].
$$

$\theta$ 与 $\mu$ 的这种对应关系是单调的：当 $\theta$ 在 0 和 1 之间变化（扫过线段 $[x,y]$）时，$\mu$ 也在 0 和 1 之间变化（扫过线段 $[P(x),P(y)]$）。因此 $P([x,y])=[P(x),P(y)]$。

现在设 $C$ 是凸集且 $C\subseteq\operatorname{dom}P$，即对所有 $x\in C$ 都有 $x_{n+1}>0$，并取 $x,y\in C$。要证明 $P(C)$ 的凸性，需要证明线段 $[P(x),P(y)]$ 包含在 $P(C)$ 中。但该线段正是线段 $[x,y]$ 在 $P$ 下的像，而 $[x,y]\subseteq C$，所以它位于 $P(C)$ 中。

凸集在透视函数下的逆像也是凸的：若 $C\subseteq\mathbf{R}^n$ 是凸的，则

$$
P^{-1}(C)=\{(x,t)\in\mathbf{R}^{n+1}\mid x/t\in C,\ t>0\}
$$

是凸的。为证明这一点，设 $(x,t)\in P^{-1}(C)$、$(y,s)\in P^{-1}(C)$，且 $0\leq\theta\leq1$。需要证明

$$
\theta(x,t)+(1-\theta)(y,s)\in P^{-1}(C),
$$

也就是

$$
\frac{\theta x+(1-\theta)y}{\theta t+(1-\theta)s}\in C
$$

（$\theta t+(1-\theta)s>0$ 显然成立）。这由

$$
\frac{\theta x+(1-\theta)y}{\theta t+(1-\theta)s}
=\mu(x/t)+(1-\mu)(y/s),
$$

其中

$$
\mu=\frac{\theta t}{\theta t+(1-\theta)s}\in[0,1],
$$

以及 $C$ 的凸性得到。

#### 线性分式函数

线性分式函数由透视函数与仿射函数复合而成。设 $g:\mathbf{R}^n\to\mathbf{R}^{m+1}$ 是仿射函数，即

$$
g(x)=
\begin{bmatrix}A\\c^T\end{bmatrix}x+
\begin{bmatrix}b\\d\end{bmatrix},
\tag{2.12}
$$

其中 $A\in\mathbf{R}^{m\times n}$、$b\in\mathbf{R}^m$、$c\in\mathbf{R}^n$、$d\in\mathbf{R}$。函数 $f:\mathbf{R}^n\to\mathbf{R}^m$ 定义为 $f=P\circ g$，即

$$
f(x)=\frac{Ax+b}{c^Tx+d},\qquad
\operatorname{dom}f=\{x\mid c^Tx+d>0\},
\tag{2.13}
$$

称为线性分式（或射影）函数。若 $c=0$ 且 $d>0$，则 $f$ 的定义域是 $\mathbf{R}^n$，并且 $f$ 是仿射函数。因此，可以把仿射函数和线性函数看作线性分式函数的特例。

**注 2.2 射影解释。** 常常方便地把线性分式函数表示为矩阵

$$
Q=\begin{bmatrix}A&b\\c^T&d\end{bmatrix}\in\mathbf{R}^{(m+1)\times(n+1)}
\tag{2.14}
$$

作用于形如 $(x,1)$ 的点，得到 $(Ax+b,c^Tx+d)$。然后把结果缩放或归一化，使最后一个分量为 1，从而得到 $(f(x),1)$。

这种表示可以通过如下方式作几何解释：把 $\mathbf{R}^n$ 与 $\mathbf{R}^{n+1}$ 中的一组射线关联。对 $\mathbf{R}^n$ 中每个点 $z$，关联 $\mathbf{R}^{n+1}$ 中的开射线 $\mathcal{P}(z)=\{t(z,1)\mid t>0\}$。该射线的最后一个分量取正值。反过来，$\mathbf{R}^{n+1}$ 中以原点为基点且最后一个分量取正值的任意射线，都可写为 $\mathcal{P}(v)=\{t(v,1)\mid t\geq0\}$，其中 $v\in\mathbf{R}^n$。这种 $\mathbf{R}^n$ 与最后分量为正的射线半空间之间的射影对应是一一且满的。

线性分式函数 (2.13) 可表示为

$$
f(x)=\mathcal{P}^{-1}(Q\mathcal{P}(x)).
$$

也就是说，从 $x\in\operatorname{dom}f$ 开始，即 $c^Tx+d>0$。先形成 $\mathbf{R}^{n+1}$ 中的射线 $\mathcal{P}(x)$。矩阵 $Q$ 对该射线作线性变换，产生另一条射线 $Q\mathcal{P}(x)$。由于 $x\in\operatorname{dom}f$，该射线的最后一个分量取正值。最后取逆射影变换得到 $f(x)$。

与透视函数一样，线性分式函数保持凸性。若 $C$ 是凸集且位于 $f$ 的定义域内（即对 $x\in C$ 有 $c^Tx+d>0$），则其像 $f(C)$ 是凸的。这直接来自前面的结果：$C$ 在仿射映射 (2.12) 下的像是凸的，而该集合再经透视函数 $P$ 的像，即 $f(C)$，也是凸的。类似地，若 $C\subseteq\mathbf{R}^m$ 是凸的，则逆像 $f^{-1}(C)$ 是凸的。

**例 2.13 条件概率。** 设随机变量 $u$ 和 $v$ 分别取值于 $\{1,\ldots,n\}$ 和 $\{1,\ldots,m\}$，并令 $p_{ij}$ 表示 $\operatorname{prob}(u=i,v=j)$。则条件概率 $f_{ij}=\operatorname{prob}(u=i\mid v=j)$ 为

$$
f_{ij}=\frac{p_{ij}}{\sum_{k=1}^n p_{kj}}.
$$

因此 $f$ 由 $p$ 通过线性分式映射得到。于是，若 $C$ 是 $(u,v)$ 的联合概率的一个凸集合，则由其得到的 $u$ 在给定 $v$ 条件下的条件概率集合也是凸的。

图 2.16 展示了集合 $C\subseteq\mathbf{R}^2$ 及其在线性分式函数

$$
f(x)=\frac{x}{x_1+x_2+1},\qquad
\operatorname{dom}f=\{(x_1,x_2)\mid x_1+x_2+1>0\}
$$

下的像。

**图 2.16** 左：集合 $C\subseteq\mathbf{R}^2$。虚线显示线性分式函数 $f(x)=x/(x_1+x_2+1)$ 的定义域边界。右：$C$ 在 $f$ 下的像。虚线显示 $f^{-1}$ 的定义域边界。

## 2.4 广义不等式

### 2.4.1 正则锥与广义不等式

若锥 $K\subseteq\mathbf{R}^n$ 满足以下条件，则称为正则锥：

- $K$ 是凸的。
- $K$ 是闭的。
- $K$ 是实心的，即内部非空。
- $K$ 是尖的，即不包含直线（等价地，$x\in K$、$-x\in K\Rightarrow x=0$）。

正则锥 $K$ 可用于定义广义不等式，即 $\mathbf{R}^n$ 上的一个偏序，它具有实数上标准序的许多性质。与正则锥 $K$ 关联的 $\mathbf{R}^n$ 上的偏序定义为

$$
x\preceq_K y\quad\Longleftrightarrow\quad y-x\in K.
$$

也写作 $x\succeq_K y$ 表示 $y\preceq_K x$。类似地，定义关联的严格偏序为

$$
x\prec_K y\quad\Longleftrightarrow\quad y-x\in\operatorname{int}K,
$$

并写作 $x\succ_K y$ 表示 $y\prec_K x$。（为区分广义不等式 $\preceq_K$ 和严格广义不等式，有时称 $\preceq_K$ 为非严格广义不等式。）

当 $K=\mathbf{R}_+$ 时，偏序 $\preceq_K$ 就是 $\mathbf{R}$ 上通常的序 $\leq$，严格偏序 $\prec_K$ 就是通常的严格序 $<$。因此，广义不等式包含了 $\mathbf{R}$ 上普通（非严格和严格）不等式作为特例。

**例 2.14 非负正交锥与逐分量不等式。** 非负正交锥 $K=\mathbf{R}^n_+$ 是正则锥。关联的广义不等式 $\preceq_K$ 对应向量之间的逐分量不等式：$x\preceq_K y$ 表示 $x_i\leq y_i$，$i=1,\ldots,n$。关联的严格不等式对应逐分量严格不等式：$x\prec_K y$ 表示 $x_i<y_i$，$i=1,\ldots,n$。

非负正交锥关联的非严格和严格偏序出现得非常频繁，因此我们省略下标 $\mathbf{R}^n_+$；当符号 $\preceq$ 或 $\prec$ 出现在向量之间时，默认其含义如此。

**例 2.15 正半定锥与矩阵不等式。** 正半定锥 $\mathbf{S}^n_+$ 是 $\mathbf{S}^n$ 中的正则锥。关联的广义不等式 $\preceq_K$ 就是通常的矩阵不等式：$X\preceq_KY$ 表示 $Y-X$ 是正半定矩阵。$\mathbf{S}^n_+$ 在 $\mathbf{S}^n$ 中的内部由正定矩阵组成，所以严格广义不等式也与对称矩阵之间通常的严格不等式一致：$X\prec_KY$ 表示 $Y-X$ 是正定矩阵。

这里偏序也出现得非常频繁，因此同样省略下标：对称矩阵之间简单写作 $X\preceq Y$ 或 $X\prec Y$。默认这些广义不等式是相对于正半定锥而言的。

**例 2.16 在 $[0,1]$ 上非负的多项式锥。** 定义

$$
K=\{c\in\mathbf{R}^n\mid c_1+c_2t+\cdots+c_nt^{n-1}\geq0\ \text{for }t\in[0,1]\},
\tag{2.15}
$$

即 $K$ 是在区间 $[0,1]$ 上非负的 $n-1$ 次多项式（系数）锥。可以证明 $K$ 是正则锥；其内部是在区间 $[0,1]$ 上为正的多项式的系数集合。

两个向量 $c,d\in\mathbf{R}^n$ 满足 $c\preceq_Kd$ 当且仅当

$$
c_1+c_2t+\cdots+c_nt^{n-1}
\leq d_1+d_2t+\cdots+d_nt^{n-1}
$$

对所有 $t\in[0,1]$ 成立。

#### 广义不等式的性质

广义不等式 $\preceq_K$ 满足许多性质，例如：

- 对加法保持：若 $x\preceq_K y$ 且 $u\preceq_K v$，则 $x+u\preceq_K y+v$。
- 传递性：若 $x\preceq_K y$ 且 $y\preceq_K z$，则 $x\preceq_K z$。
- 对非负缩放保持：若 $x\preceq_K y$ 且 $\alpha\geq0$，则 $\alpha x\preceq_K \alpha y$。
- 自反性：$x\preceq_K x$。
- 反对称性：若 $x\preceq_K y$ 且 $y\preceq_K x$，则 $x=y$。
- 对极限保持：若 $x_i\preceq_K y_i$，$i=1,2,\ldots$，且 $x_i\to x$、$y_i\to y$，则 $x\preceq_K y$。

相应的严格广义不等式 $\prec_K$ 例如满足：

- 若 $x\prec_K y$，则 $x\preceq_K y$。
- 若 $x\prec_K y$ 且 $u\preceq_K v$，则 $x+u\prec_K y+v$。
- 若 $x\prec_K y$ 且 $\alpha>0$，则 $\alpha x\prec_K\alpha y$。
- $x\not\prec_K x$。
- 若 $x\prec_K y$，则对足够小的 $u$ 和 $v$，有 $x+u\prec_K y+v$。

这些性质来自 $\preceq_K$ 和 $\prec_K$ 的定义以及正则锥的性质；见习题 2.30。

### 2.4.2 最小元素与极小元素

广义不等式记号（即 $\preceq_K$、$\prec_K$）意在提示它与 $\mathbf{R}$ 上普通不等式（即 $\leq$、$<$）的类比。尽管普通不等式的许多性质对广义不等式也成立，但有些重要性质并不成立。最明显的区别是 $\mathbf{R}$ 上的 $\leq$ 是线性序：任意两点都可比较，也就是说要么 $x\leq y$，要么 $y\leq x$。这个性质对其他广义不等式并不成立。一个后果是，在广义不等式背景下，最小值和最大值这样的概念更复杂。本节简要讨论这一点。

若对每个 $y\in S$ 都有 $x\preceq_K y$，则称 $x\in S$ 是 $S$ 的最小元素（相对于广义不等式 $\preceq_K$）。最大元素类似定义。若集合有最小（最大）元素，则它是唯一的。相关概念是极小元素。若 $y\in S$ 且 $y\preceq_K x$ 只在 $y=x$ 时成立，则称 $x\in S$ 是 $S$ 的极小元素（相对于 $\preceq_K$）。极大元素类似定义。一个集合可以有许多不同的极小（极大）元素。

可以用简单集合记号描述最小元素和极小元素。点 $x\in S$ 是 $S$ 的最小元素，当且仅当

$$
S\subseteq x+K.
$$

这里 $x+K$ 表示所有与 $x$ 可比较且大于等于 $x$ 的点（按 $\preceq_K$）。点 $x\in S$ 是极小元素，当且仅当

$$
(x-K)\cap S=\{x\}.
$$

这里 $x-K$ 表示所有与 $x$ 可比较且小于等于 $x$ 的点；它与 $S$ 的公共点只有 $x$。

当 $K=\mathbf{R}_+$ 时，它诱导 $\mathbf{R}$ 上通常的序，极小与最小概念相同，并与集合最小元素的通常定义一致。

**例 2.17** 考虑锥 $\mathbf{R}^2_+$，它在 $\mathbf{R}^2$ 中诱导逐分量不等式。这里可以给出极小元素和最小元素的简单几何描述。不等式 $x\preceq y$ 表示 $y$ 位于 $x$ 的右上方。说 $x\in S$ 是集合 $S$ 的最小元素，意味着 $S$ 中所有其他点都位于其右上方。说 $x$ 是集合 $S$ 的极小元素，意味着 $S$ 中没有其他点位于 $x$ 的左下方。图 2.17 展示了这一点。

**图 2.17** 左：集合 $S_1$ 相对于 $\mathbf{R}^2$ 中逐分量不等式具有最小元素 $x_1$。右：点 $x_2$ 是 $S_2$ 的极小点。

**例 2.18 对称矩阵集合的最小元素与极小元素。** 对每个 $A\in\mathbf{S}^n_{++}$，关联以原点为中心的椭球

$$
E_A=\{x\mid x^TA^{-1}x\leq1\}.
$$

有 $A\preceq B$ 当且仅当 $E_A\subseteq E_B$。

给定 $v_1,\ldots,v_k\in\mathbf{R}^n$，定义

$$
S=\{P\in\mathbf{S}^n_{++}\mid v_i^TP^{-1}v_i\leq1,\ i=1,\ldots,k\},
$$

它对应于包含点 $v_1,\ldots,v_k$ 的椭球集合。集合 $S$ 没有最小元素：对任意包含这些点的椭球，都可以找到另一个也包含这些点但与其不可比较的椭球。若一个椭球包含这些点，但没有更小的椭球包含它们，则该椭球是极小的。图 2.18 给出了 $\mathbf{R}^2$ 中 $k=2$ 的例子。

**图 2.18** $\mathbf{R}^2$ 中以原点为中心并包含给定点的三个椭球。$E_1$ 和 $E_3$ 不是极小的；$E_2$ 是极小的。

## 2.5 分离超平面与支撑超平面

### 2.5.1 分离超平面定理

本节描述一个后面很重要的思想：用超平面或仿射函数分离不相交的凸集。基本结果是分离超平面定理：设 $C$ 和 $D$ 是非空不相交凸集，即 $C\cap D=\emptyset$。则存在 $a\neq0$ 和 $b$，使得对所有 $x\in C$ 有 $a^Tx\leq b$，对所有 $x\in D$ 有 $a^Tx\geq b$。换言之，仿射函数 $a^Tx-b$ 在 $C$ 上非正，在 $D$ 上非负。超平面 $\{x\mid a^Tx=b\}$ 称为集合 $C$ 和 $D$ 的分离超平面，或称它分离 $C$ 和 $D$。图 2.19 展示了这一点。

**图 2.19** 超平面 $\{x\mid a^Tx=b\}$ 分离不相交凸集 $C$ 和 $D$。仿射函数 $a^Tx-b$ 在 $C$ 上非正，在 $D$ 上非负。

#### 分离超平面定理的证明

这里考虑一个特殊情形，并把推广到一般情形留作习题（习题 2.22）。我们假设 $C$ 和 $D$ 之间的欧几里得距离

$$
\operatorname{dist}(C,D)
=\inf\{\|u-v\|_2\mid u\in C,\ v\in D\}
$$

为正，并且存在点 $c\in C$ 和 $d\in D$ 达到最小距离，即 $\|c-d\|_2=\operatorname{dist}(C,D)$。（例如，当 $C$ 和 $D$ 都闭且其中一个有界时，这些条件成立。）

定义

$$
a=d-c,\qquad b=\frac{\|d\|_2^2-\|c\|_2^2}{2}.
$$

我们将证明仿射函数

$$
f(x)=a^Tx-b=(d-c)^T\left(x-\frac{1}{2}(d+c)\right)
$$

在 $C$ 上非正，在 $D$ 上非负，即超平面 $\{x\mid a^Tx=b\}$ 分离 $C$ 和 $D$。这个超平面垂直于 $c$ 与 $d$ 之间的线段，并经过其中点，如图 2.20 所示。

**图 2.20** 两个凸集之间分离超平面的构造。点 $c\in C$ 和 $d\in D$ 是两个集合中彼此最近的一对点。分离超平面垂直平分 $c$ 与 $d$ 之间的线段。

先证明 $f$ 在 $D$ 上非负。证明 $f$ 在 $C$ 上非正是类似的（或交换 $C$ 和 $D$ 并考虑 $-f$）。假设存在 $u\in D$ 使得

$$
f(u)=(d-c)^T\left(u-\frac{1}{2}(d+c)\right)<0.
\tag{2.16}
$$

可将 $f(u)$ 写成

$$
f(u)=(d-c)^T(u-d+(1/2)(d-c))
=(d-c)^T(u-d)+(1/2)\|d-c\|_2^2.
$$

由 (2.16) 可知 $(d-c)^T(u-d)<0$。又有

$$
\left.\frac{d}{dt}\|d+t(u-d)-c\|_2^2\right|_{t=0}
=2(d-c)^T(u-d)<0,
$$

所以对某个足够小且 $t\leq1$ 的 $t>0$，

$$
\|d+t(u-d)-c\|_2<\|d-c\|_2,
$$

也就是说点 $d+t(u-d)$ 比 $d$ 更接近 $c$。由于 $D$ 是凸的并且包含 $d$ 和 $u$，有 $d+t(u-d)\in D$。这与 $d$ 是 $D$ 中离 $C$ 最近的点相矛盾。

**例 2.19 仿射集与凸集的分离。** 设 $C$ 是凸集，$D$ 是仿射集，即 $D=\{Fu+g\mid u\in\mathbf{R}^m\}$，其中 $F\in\mathbf{R}^{n\times m}$。假设 $C$ 和 $D$ 不相交，则由分离超平面定理，存在 $a\neq0$ 和 $b$，使得对所有 $x\in C$ 有 $a^Tx\leq b$，对所有 $x\in D$ 有 $a^Tx\geq b$。

现在 $a^Tx\geq b$ 对所有 $x\in D$ 成立，意味着 $a^TFu\geq b-a^Tg$ 对所有 $u\in\mathbf{R}^m$ 成立。但一个线性函数只有在为零时才可能在 $\mathbf{R}^m$ 上有下界，所以 $a^TF=0$（并且 $b\leq a^Tg$）。

因此，存在 $a\neq0$ 使得 $F^Ta=0$，并且对所有 $x\in C$ 有 $a^Tx\leq a^Tg$。

#### 严格分离

上面构造的分离超平面满足更强条件：对所有 $x\in C$ 有 $a^Tx<b$，对所有 $x\in D$ 有 $a^Tx>b$。这称为集合 $C$ 和 $D$ 的严格分离。简单例子表明，一般而言，不相交凸集未必能被某个超平面严格分离（即使这些集合是闭的；见习题 2.23）。不过，在许多特殊情况下可以建立严格分离。

**例 2.20 点与闭凸集的严格分离。** 设 $C$ 是闭凸集，且 $x_0\notin C$。则存在一个超平面严格分离 $x_0$ 与 $C$。

为说明这一点，注意对某个 $\epsilon>0$，集合 $C$ 和 $B(x_0,\epsilon)$ 不相交。由分离超平面定理，存在 $a\neq0$ 和 $b$，使得对 $x\in C$ 有 $a^Tx\leq b$，对 $x\in B(x_0,\epsilon)$ 有 $a^Tx\geq b$。

利用 $B(x_0,\epsilon)=\{x_0+u\mid \|u\|_2\leq\epsilon\}$，第二个条件可写为

$$
a^T(x_0+u)\geq b\quad\text{for all }\|u\|_2\leq\epsilon.
$$

使左端最小的 $u$ 是 $u=-\epsilon a/\|a\|_2$；代入得到

$$
a^Tx_0-\epsilon\|a\|_2\geq b.
$$

因此仿射函数

$$
f(x)=a^Tx-b-\epsilon\|a\|_2/2
$$

在 $C$ 上为负，在 $x_0$ 处为正。

作为直接结果，可以证明前面提到的事实：闭凸集是所有包含它的半空间的交。事实上，设 $C$ 闭且凸，$S$ 是所有包含 $C$ 的半空间的交。显然 $x\in C\Rightarrow x\in S$。反过来，若存在 $x\in S$ 且 $x\notin C$，则严格分离结果给出一个严格分离 $x$ 与 $C$ 的超平面，也就是说存在一个包含 $C$ 但不包含 $x$ 的半空间。这意味着 $x\notin S$，矛盾。

#### 反向分离超平面定理

分离超平面定理的反命题（即存在分离超平面推出 $C$ 和 $D$ 不相交）并不成立，除非对 $C$ 或 $D$ 施加额外条件，甚至仅有凸性也不够。一个简单反例是 $C=D=\{0\}\subseteq\mathbf{R}$。这里超平面 $x=0$ 分离 $C$ 和 $D$。

通过对 $C$ 和 $D$ 加条件，可以得到各种反向分离定理。一个很简单的例子是：设 $C$ 和 $D$ 是凸集，且 $C$ 是开集，并且存在一个仿射函数 $f$，它在 $C$ 上非正，在 $D$ 上非负。则 $C$ 和 $D$ 不相交。（证明：先注意 $f$ 在 $C$ 上必须为负；若 $f$ 在 $C$ 中某点为零，则它会在该点附近取正值，矛盾。于是 $f$ 在 $C$ 上为负、在 $D$ 上非负，所以 $C$ 和 $D$ 必不相交。）将这个反命题与分离超平面定理合并，可得：任意两个凸集 $C$ 和 $D$，只要至少一个是开集，则它们不相交当且仅当存在分离超平面。

**例 2.21 严格线性不等式的择一定理。** 推导严格线性不等式组

$$
Ax\prec b
\tag{2.17}
$$

可解的充要条件。该不等式组不可行，当且仅当两个凸集

$$
C=\{b-Ax\mid x\in\mathbf{R}^n\},\qquad
D=\mathbf{R}^m_{++}=\{y\in\mathbf{R}^m\mid y\succ0\}
$$

不相交。集合 $D$ 是开的，$C$ 是仿射集。因此由上述结果，$C$ 和 $D$ 不相交当且仅当存在分离超平面，即存在非零 $\lambda\in\mathbf{R}^m$ 和 $\mu\in\mathbf{R}$，使得 $\lambda^Ty\leq\mu$ 在 $C$ 上成立，$\lambda^Ty\geq\mu$ 在 $D$ 上成立。

这些条件都可简化。第一个条件表示 $\lambda^T(b-Ax)\leq\mu$ 对所有 $x$ 成立。这意味着（如例 2.19）$A^T\lambda=0$ 且 $\lambda^Tb\leq\mu$。第二个不等式表示 $\lambda^Ty\geq\mu$ 对所有 $y\succ0$ 成立。这意味着 $\mu\leq0$ 且 $\lambda\succeq0$、$\lambda\neq0$。

综合起来，严格不等式组 (2.17) 不可行，当且仅当存在 $\lambda\in\mathbf{R}^m$ 使得

$$
\lambda\neq0,\qquad
\lambda\succeq0,\qquad
A^T\lambda=0,\qquad
\lambda^Tb\leq0.
\tag{2.18}
$$

这同样是关于变量 $\lambda\in\mathbf{R}^m$ 的线性不等式和线性方程组。我们说 (2.17) 和 (2.18) 构成一对择一系统：对任意数据 $A$ 和 $b$，二者恰有一个可解。

### 2.5.2 支撑超平面

设 $C\subseteq\mathbf{R}^n$，且 $x_0$ 是其边界 $\operatorname{bd}C$ 上的点，即

$$
x_0\in\operatorname{bd}C=\operatorname{cl}C\setminus\operatorname{int}C.
$$

若 $a\neq0$ 满足对所有 $x\in C$ 有 $a^Tx\leq a^Tx_0$，则超平面 $\{x\mid a^Tx=a^Tx_0\}$ 称为 $C$ 在点 $x_0$ 处的支撑超平面。这等价于说点 $x_0$ 与集合 $C$ 被超平面 $\{x\mid a^Tx=a^Tx_0\}$ 分离。几何解释是，超平面 $\{x\mid a^Tx=a^Tx_0\}$ 在 $x_0$ 处与 $C$ 相切，并且半空间 $\{x\mid a^Tx\leq a^Tx_0\}$ 包含 $C$。图 2.21 展示了这一点。

**图 2.21** 超平面 $\{x\mid a^Tx=a^Tx_0\}$ 在 $x_0$ 处支撑 $C$。

一个基本结果称为支撑超平面定理：对任意非空凸集 $C$ 和任意 $x_0\in\operatorname{bd}C$，都存在 $C$ 在 $x_0$ 处的支撑超平面。支撑超平面定理可由分离超平面定理直接证明。分两种情况。若 $C$ 的内部非空，则把分离超平面定理应用于集合 $\{x_0\}$ 和 $\operatorname{int}C$，即可得到结果。若 $C$ 的内部为空，则 $C$ 必位于维数小于 $n$ 的仿射集中，任何包含该仿射集的超平面都包含 $C$ 和 $x_0$，并且是一个（平凡的）支撑超平面。

支撑超平面定理也有一个部分反命题：若集合闭、内部非空，并且在其边界每一点都有支撑超平面，则该集合是凸的。（见习题 2.27。）

## 2.6 对偶锥与广义不等式

### 2.6.1 对偶锥

设 $K$ 是一个锥。集合

$$
K^*=\{y\mid x^Ty\geq0\ \text{for all }x\in K\}
\tag{2.19}
$$

称为 $K$ 的对偶锥。顾名思义，$K^*$ 是一个锥，并且总是凸的，即使原锥 $K$ 不是凸的也如此（见习题 2.31）。

从几何上看，$y\in K^*$ 当且仅当 $-y$ 是一个在原点支撑 $K$ 的超平面的法向量。图 2.22 展示了这一点。

**图 2.22** 左：以内法向量 $y$ 定义的半空间包含锥 $K$，所以 $y\in K^*$。右：以内法向量 $z$ 定义的半空间不包含 $K$，所以 $z\notin K^*$。

**例 2.22 子空间。** 子空间 $V\subseteq\mathbf{R}^n$（它也是锥）的对偶锥是其正交补 $V^\perp=\{y\mid v^Ty=0\ \text{for all }v\in V\}$。

**例 2.23 非负正交锥。** 锥 $\mathbf{R}^n_+$ 是自对偶的：

$$
x^Ty\geq0\ \text{for all }x\succeq0
\quad\Longleftrightarrow\quad y\succeq0.
$$

这样的锥称为自对偶锥。

**例 2.24 正半定锥。** 在对称 $n\times n$ 矩阵集合 $\mathbf{S}^n$ 上，我们使用标准内积 $\operatorname{tr}(XY)=\sum_{i,j=1}^n X_{ij}Y_{ij}$（见 §A.1.1）。正半定锥 $\mathbf{S}^n_+$ 是自对偶的，即对 $X,Y\in\mathbf{S}^n$，

$$
\operatorname{tr}(XY)\geq0\ \text{for all }X\succeq0
\quad\Longleftrightarrow\quad Y\succeq0.
$$

证明如下。若 $Y\notin\mathbf{S}^n_+$，则存在 $q\in\mathbf{R}^n$ 使得

$$
q^TYq=\operatorname{tr}(qq^TY)<0.
$$

于是正半定矩阵 $X=qq^T$ 满足 $\operatorname{tr}(XY)<0$，所以 $Y\notin(\mathbf{S}^n_+)^*$。

反过来，设 $X,Y\in\mathbf{S}^n_+$。将 $X$ 写成特征值分解 $X=\sum_{i=1}^n\lambda_iq_iq_i^T$，其中特征值 $\lambda_i\geq0$，$i=1,\ldots,n$。则

$$
\operatorname{tr}(YX)
=\operatorname{tr}\left(Y\sum_{i=1}^n\lambda_iq_iq_i^T\right)
=\sum_{i=1}^n\lambda_i q_i^TYq_i\geq0.
$$

这说明 $Y\in(\mathbf{S}^n_+)^*$。

**例 2.25 范数锥的对偶。** 设 $\|\cdot\|$ 是 $\mathbf{R}^n$ 上的范数。关联锥 $K=\{(x,t)\in\mathbf{R}^{n+1}\mid \|x\|\leq t\}$ 的对偶是由对偶范数定义的锥：

$$
K^*=\{(u,v)\in\mathbf{R}^{n+1}\mid \|u\|_*\leq v\},
$$

其中对偶范数为 $\|u\|_*=\sup\{u^Tx\mid \|x\|\leq1\}$（见 (A.1.6)）。

要证明这个结果，需要说明

$$
x^Tu+tv\geq0\ \text{whenever }\|x\|\leq t
\quad\Longleftrightarrow\quad
\|u\|_*\leq v.
\tag{2.20}
$$

先证明右侧条件推出左侧条件。设 $\|u\|_*\leq v$，且对某个 $t>0$ 有 $\|x\|\leq t$。（若 $t=0$，则 $x$ 必为零，显然 $u^Tx+vt\geq0$。）由对偶范数定义及 $\|-x/t\|\leq1$，

$$
u^T(-x/t)\leq\|u\|_*\leq v,
$$

从而 $u^Tx+vt\geq0$。

再证明 (2.20) 左侧条件推出右侧条件。若 $\|u\|_*>v$，即右侧条件不成立，则由对偶范数定义，存在 $\|x\|\leq1$ 且 $x^Tu>v$。取 $t=1$，有

$$
u^T(-x)+v<0,
$$

这与 (2.20) 左侧条件矛盾。

对偶锥满足若干性质，例如：

- $K^*$ 是闭凸锥。
- $K_1\subseteq K_2$ 蕴含 $K_2^*\subseteq K_1^*$。
- 若 $K$ 内部非空，则 $K^*$ 是尖的。
- 若 $K$ 的闭包是尖的，则 $K^*$ 内部非空。
- $K^{**}$ 是 $K$ 的凸包的闭包。（因此若 $K$ 凸且闭，则 $K^{**}=K$。）

见习题 2.31。这些性质说明，若 $K$ 是正则锥，则其对偶 $K^*$ 也是正则锥，并且 $K^{**}=K$。

### 2.6.2 对偶广义不等式

现在设凸锥 $K$ 是正则锥，因此它诱导广义不等式 $\preceq_K$。则其对偶锥 $K^*$ 也是正则的，因此也诱导一个广义不等式。我们把广义不等式 $\preceq_{K^*}$ 称为广义不等式 $\preceq_K$ 的对偶。

联系一个广义不等式及其对偶的一些重要性质是：

- $x\preceq_K y$ 当且仅当对所有 $\lambda\succeq_{K^*}0$ 都有 $\lambda^Tx\leq\lambda^Ty$。
- $x\prec_K y$ 当且仅当对所有 $\lambda\succeq_{K^*}0$、$\lambda\neq0$ 都有 $\lambda^Tx<\lambda^Ty$。

由于 $K=K^{**}$，与 $\preceq_{K^*}$ 关联的对偶广义不等式就是 $\preceq_K$，所以交换广义不等式与其对偶后，这些性质仍成立。作为具体例子，$\lambda\preceq_{K^*}\mu$ 当且仅当对所有 $x\succeq_K0$ 都有 $\lambda^Tx\leq\mu^Tx$。

**例 2.26 严格线性广义不等式的择一定理。** 设 $K\subseteq\mathbf{R}^m$ 是正则锥。考虑严格广义不等式

$$
Ax\prec_K b,
\tag{2.21}
$$

其中 $x\in\mathbf{R}^n$。

我们推导该不等式的择一定理。假设它不可行，即仿射集 $\{b-Ax\mid x\in\mathbf{R}^n\}$ 与开凸集 $\operatorname{int}K$ 不相交。则存在分离超平面，即存在非零 $\lambda\in\mathbf{R}^m$ 和 $\mu\in\mathbf{R}$，使得对所有 $x$ 有 $\lambda^T(b-Ax)\leq\mu$，对所有 $y\in\operatorname{int}K$ 有 $\lambda^Ty\geq\mu$。第一个条件推出 $A^T\lambda=0$ 且 $\lambda^Tb\leq\mu$。第二个条件推出对所有 $y\in K$ 有 $\lambda^Ty\geq\mu$，这只有在 $\lambda\in K^*$ 且 $\mu\leq0$ 时才可能。

综合起来，若 (2.21) 不可行，则存在 $\lambda$ 使得

$$
\lambda\neq0,\qquad
\lambda\succeq_{K^*}0,\qquad
A^T\lambda=0,\qquad
\lambda^Tb\leq0.
\tag{2.22}
$$

现在证明反向：若 (2.22) 成立，则不等式组 (2.21) 不可能可行。若两组不等式同时成立，则由于 $\lambda\neq0$、$\lambda\succeq_{K^*}0$ 且 $b-Ax\succ_K0$，有 $\lambda^T(b-Ax)>0$。但由 $A^T\lambda=0$ 可得 $\lambda^T(b-Ax)=\lambda^Tb\leq0$，矛盾。

因此，不等式组 (2.21) 和 (2.22) 是择一系统：对任意数据 $A,b$，恰有一个可行。（这推广了特殊情形 $K=\mathbf{R}^m_+$ 下的择一系统 (2.17)、(2.18)。）

### 2.6.3 通过对偶不等式刻画最小元素与极小元素

可以使用对偶广义不等式，刻画集合 $S\subseteq\mathbf{R}^m$（可以非凸）中相对于正则锥 $K$ 诱导的广义不等式的最小元素和极小元素。

#### 最小元素的对偶刻画

先考虑最小元素的刻画：$x$ 是 $S$ 关于广义不等式 $\preceq_K$ 的最小元素，当且仅当对所有 $\lambda\succ_{K^*}0$，$x$ 是在 $z\in S$ 上最小化 $\lambda^Tz$ 的唯一解。从几何上看，这意味着对任意 $\lambda\succ_{K^*}0$，超平面

$$
\{z\mid \lambda^T(z-x)=0\}
$$

都是 $S$ 在 $x$ 处的严格支撑超平面。（严格支撑超平面指该超平面与 $S$ 只在点 $x$ 相交。）注意这里不要求集合 $S$ 凸。图 2.23 展示了这一点。

**图 2.23** 最小元素的对偶刻画。点 $x$ 是集合 $S$ 关于 $\mathbf{R}^2_+$ 的最小元素。这等价于：对每个 $\lambda\succ0$，超平面 $\{z\mid \lambda^T(z-x)=0\}$ 在 $x$ 处严格支撑 $S$。

证明如下。若 $x$ 是 $S$ 的最小元素，即对所有 $z\in S$ 都有 $x\preceq_K z$，并取 $\lambda\succ_{K^*}0$。若 $z\in S$ 且 $z\neq x$，则 $z-x\succeq_K0$ 且 $z-x\neq0$。由 $\lambda\succ_{K^*}0$ 得 $\lambda^T(z-x)>0$。由于 $z$ 是任意不等于 $x$ 的 $S$ 中元素，$x$ 是在 $S$ 上最小化 $\lambda^Tz$ 的唯一解。

反过来，若对所有 $\lambda\succ_{K^*}0$，$x$ 都是 $S$ 上 $\lambda^Tz$ 的唯一最小化点，但 $x$ 不是最小元素，则存在 $z\in S$ 使得 $z\not\succeq_K x$。由于 $z-x\not\succeq_K0$，存在 $\tilde{\lambda}\succeq_{K^*}0$ 使得 $\tilde{\lambda}^T(z-x)<0$。于是对 $\tilde{\lambda}$ 附近的某个 $\lambda\succ_{K^*}0$，有 $\lambda^T(z-x)<0$，这与 $x$ 是 $S$ 上 $\lambda^Tz$ 的唯一最小化点矛盾。

#### 极小元素的对偶刻画

再来看极小元素的类似刻画。这里必要条件和充分条件之间存在差距。若 $\lambda\succ_{K^*}0$ 且 $x$ 在 $z\in S$ 上最小化 $\lambda^Tz$，则 $x$ 是极小的。图 2.24 展示了这一点。

**图 2.24** 集合 $S\subseteq\mathbf{R}^2$ 的极小点集合（相对于 $\mathbf{R}^2_+$）是其左下边界的深色部分。若 $\lambda_1\succ0$，则 $\lambda_1^Tz$ 的最小化点 $x_1$ 是极小的；若 $\lambda_2\succ0$，则 $\lambda_2^Tz$ 的最小化点 $x_2$ 也是极小的。

证明如下。假设 $\lambda\succ_{K^*}0$，且 $x$ 在 $S$ 上最小化 $\lambda^Tz$，但 $x$ 不是极小的，即存在 $z\in S$、$z\neq x$，且 $z\preceq_Kx$。则 $\lambda^T(x-z)>0$，这与 $x$ 是 $\lambda^Tz$ 在 $S$ 上的最小化点相矛盾。

反命题一般不成立：点 $x$ 可以是 $S$ 中的极小元素，但不对任何 $\lambda$ 成为 $\lambda^Tz$ 在 $z\in S$ 上的最小化点，如图 2.25 所示。该图提示凸性在反命题中起重要作用，事实确实如此。若集合 $S$ 是凸的，则对任意极小元素 $x$，存在非零 $\lambda\succeq_{K^*}0$，使得 $x$ 在 $z\in S$ 上最小化 $\lambda^Tz$。

**图 2.25** 点 $x$ 是 $S\subseteq\mathbf{R}^2$ 关于 $\mathbf{R}^2_+$ 的极小元素，但不存在任何 $\lambda$ 使得 $x$ 在 $z\in S$ 上最小化 $\lambda^Tz$。

为证明这一点，设 $x$ 是极小的，这意味着 $((x-K)\setminus\{x\})\cap S=\emptyset$。对凸集 $(x-K)\setminus\{x\}$ 和 $S$ 应用分离超平面定理，可得存在 $\lambda\neq0$ 和 $\mu$，使得对所有 $y\in K$ 有 $\lambda^T(x-y)\leq\mu$，对所有 $z\in S$ 有 $\lambda^Tz\geq\mu$。由第一个不等式可得 $\lambda\succeq_{K^*}0$。由于 $x\in S$ 且 $x\in x-K$，有 $\lambda^Tx=\mu$，所以第二个不等式说明 $\mu$ 是 $\lambda^Tz$ 在 $S$ 上的最小值。因此 $x$ 是 $\lambda^Tz$ 在 $S$ 上的最小化点，其中 $\lambda\neq0$、$\lambda\succeq_{K^*}0$。

这个反向定理不能加强为 $\lambda\succ_{K^*}0$。例子表明，点 $x$ 可以是凸集 $S$ 的极小点，但不对任何 $\lambda\succ_{K^*}0$ 成为 $\lambda^Tz$ 在 $S$ 上的最小化点（见图 2.26 左）。同样，并非任意满足 $\lambda\succeq_{K^*}0$ 的 $\lambda^Tz$ 最小化点都是极小点（见图 2.26 右）。

**图 2.26** 左：点 $x_1\in S_1$ 是极小的，但不是任何 $\lambda\succ0$ 的 $\lambda^Tz$ 最小化点（它确实对 $\lambda=(1,0)$ 最小化）。右：点 $x_2\in S_2$ 不是极小的，但它对 $\lambda=(0,1)\succeq0$ 最小化 $\lambda^Tz$。

**例 2.27 Pareto 最优生产前沿。** 考虑一种产品，其制造需要 $n$ 种资源（例如劳动力、电力、天然气、水）。产品可以用多种方式制造或生产。对每种生产方法，关联一个资源向量 $x\in\mathbf{R}^n$，其中 $x_i$ 表示该方法制造产品所消耗的第 $i$ 种资源量。假设 $x_i\geq0$（即生产方法消耗资源），并且资源是有价值的（所以少用任意一种资源都更好）。

生产集 $P\subseteq\mathbf{R}^n$ 定义为所有对应于某种生产方法的资源向量 $x$ 的集合。

资源向量为 $P$ 中关于逐分量不等式的极小元素的生产方法，称为 Pareto 最优或有效。$P$ 的极小元素集合称为有效生产前沿。

可以给出 Pareto 最优性的简单解释。若一个生产方法的资源向量为 $x$，另一个为 $y$，并且对所有 $i$ 有 $x_i\leq y_i$，且对某个 $i$ 有 $x_i<y_i$，则称前者优于后者。换言之，一个生产方法若对每种资源都不比另一个方法用得更多，并且至少对一种资源用得更少，则它更好。这对应于 $x\preceq y$、$x\neq y$。因此可以说：若不存在更好的生产方法，则该生产方法是 Pareto 最优或有效的。

通过在生产向量集合 $P$ 上最小化

$$
\lambda^Tx=\lambda_1x_1+\cdots+\lambda_nx_n
$$

并使用任意满足 $\lambda\succ0$ 的 $\lambda$，可以找到 Pareto 最优生产方法（即极小资源向量）。

这里向量 $\lambda$ 有简单解释：$\lambda_i$ 是第 $i$ 种资源的价格。通过在 $P$ 上最小化 $\lambda^Tx$，我们在寻找给定资源价格 $\lambda_i$ 下总成本最低的生产方法。只要价格为正，得到的生产方法保证是有效的。

这些思想如图 2.27 所示。

**图 2.27** 生产一种需要劳动力和燃料的产品的生产集 $P$。两条深色曲线表示有效生产前沿。点 $x_1,x_2,x_3$ 是有效的；点 $x_4,x_5$ 不是。点 $x_1$ 也是价格向量 $\lambda$ 下成本最低的生产方法。点 $x_2$ 是有效的，但无法通过对任何价格向量 $\lambda\succeq0$ 最小化总成本 $\lambda^Tx$ 得到。

## Bibliography

Minkowski 通常被认为最早系统研究了凸集，并引入了支撑超平面、支撑超平面定理、Minkowski 距离函数（习题 3.34）、凸集的极点等基本概念。

一些著名的早期综述包括 Bonnesen 和 Fenchel [BF48]、Eggleston [Egg58]、Klee [Kle63] 以及 Valentine [Val64]。较近的凸集几何专著包括 Lay [Lay82] 和 Webster [Web94]。Klee [Kle71]、Fenchel [Fen83]、Tikhomorov [Tik90] 和 Berger [Ber90] 对凸性及其在数学中应用的历史给出了可读性很强的概述。

线性不等式和多面集合与线性规划问题密切相关，并被广泛研究；关于线性规划的参考文献列在第 4 章末。线性不等式和线性规划历史上的一些里程碑式出版物包括 Motzkin [Mot33]、von Neumann 和 Morgenstern [vNM53]、Kantorovich [Kan60]、Koopmans [Koo51] 以及 Dantzig [Dan63]。Dantzig [Dan63, Chapter 2] 包含了截至约 1963 年的线性不等式历史综述。

广义不等式在 1960 年代被引入非线性优化中（见 Luenberger [Lue69, §8.2] 和 Isii [Isi64]），并在锥规划中被广泛使用（见第 4 章参考文献）。Bellman 和 Fan [BF63] 是关于广义线性不等式集合（相对于正半定锥）的早期论文。

关于分离超平面定理的扩展和证明，读者可参考 Rockafellar [Roc70, part III] 以及 Hiriart-Urruty 和 Lemaréchal [HUL93, volume 1, §III4]。Dantzig [Dan63, page 21] 将“择一定理”这一术语归因于 von Neumann 和 Morgenstern [vNM53, page 138]。关于择一定理的更多参考文献见第 5 章。

例 2.27 中的术语（包括 Pareto 最优性、有效生产以及 $\lambda$ 的价格解释）由 Luenberger [Lue95] 详细讨论。

凸几何在经典矩问题理论中扮演重要角色（Krein 和 Nudelman [KN77]，Karlin 和 Studden [KS66]）。一个著名例子是非负多项式锥与幂矩锥之间的对偶性；见习题 2.37。

## Exercises

### 凸性的定义

**2.1** 设 $C\subseteq\mathbf{R}^n$ 是凸集，$x_1,\ldots,x_k\in C$，并且 $\theta_1,\ldots,\theta_k\in\mathbf{R}$ 满足 $\theta_i\geq0$、$\theta_1+\cdots+\theta_k=1$。证明 $\theta_1x_1+\cdots+\theta_kx_k\in C$。（凸性的定义是 $k=2$ 时成立；你需要证明任意 $k$ 的情形。）提示：对 $k$ 使用归纳法。

**2.2** 证明：一个集合是凸的，当且仅当它与任意直线的交是凸的。证明：一个集合是仿射的，当且仅当它与任意直线的交是仿射的。

**2.3 中点凸性。** 若集合 $C$ 中任意两点 $a,b$ 的平均点或中点 $(a+b)/2$ 都在 $C$ 中，则称 $C$ 是中点凸的。显然凸集是中点凸的。可以证明，在温和条件下中点凸性蕴含凸性。作为简单情形，证明：若 $C$ 是闭的且中点凸，则 $C$ 是凸的。

**2.4** 证明集合 $S$ 的凸包是所有包含 $S$ 的凸集的交。（同样方法可用于证明 $S$ 的锥包、仿射包或线性包分别是所有包含 $S$ 的锥集、仿射集或子空间的交。）

### 例子

**2.5** 两个平行超平面 $\{x\in\mathbf{R}^n\mid a^Tx=b_1\}$ 和 $\{x\in\mathbf{R}^n\mid a^Tx=b_2\}$ 之间的距离是多少？

**2.6** 何时一个半空间包含另一个半空间？给出

$$
\{x\mid a^Tx\leq b\}\subseteq\{x\mid \tilde{a}^Tx\leq\tilde{b}\}
$$

成立的条件，其中 $a\neq0$、$\tilde{a}\neq0$。同时求出两个半空间相等的条件。

**2.7 半空间的 Voronoi 描述。** 设 $a$ 和 $b$ 是 $\mathbf{R}^n$ 中不同的点。证明所有到 $a$ 比到 $b$ 更近（按欧几里得范数）的点集合

$$
\{x\mid \|x-a\|_2\leq\|x-b\|_2\}
$$

是一个半空间。将其显式描述为 $c^Tx\leq d$ 形式的不等式，并画图。

**2.8** 下列集合 $S$ 哪些是多面体？若可能，将 $S$ 表示为 $S=\{x\mid Ax\preceq b,\ Fx=g\}$。

(a) $S=\{y_1a_1+y_2a_2\mid -1\leq y_1\leq1,\ -1\leq y_2\leq1\}$，其中 $a_1,a_2\in\mathbf{R}^n$。

(b) $S=\{x\in\mathbf{R}^n\mid x\succeq0,\ \mathbf{1}^Tx=1,\ \sum_{i=1}^n x_ia_i=b_1,\ \sum_{i=1}^n x_ia_i^2=b_2\}$，其中 $a_1,\ldots,a_n\in\mathbf{R}$，$b_1,b_2\in\mathbf{R}$。

(c) $S=\{x\in\mathbf{R}^n\mid x\succeq0,\ x^Ty\leq1\ \text{for all }y\text{ with }\|y\|_2=1\}$。

(d) $S=\{x\in\mathbf{R}^n\mid x\succeq0,\ x^Ty\leq1\ \text{for all }y\text{ with }\sum_{i=1}^n|y_i|=1\}$。

**2.9 Voronoi 集与多面体分解。** 设 $x_0,\ldots,x_K\in\mathbf{R}^n$ 互不相同。考虑到 $x_0$ 比到其他 $x_i$ 更近的点集合

$$
V=\{x\in\mathbf{R}^n\mid \|x-x_0\|_2\leq\|x-x_i\|_2,\ i=1,\ldots,K\}.
$$

$V$ 称为 $x_0$ 关于 $x_1,\ldots,x_K$ 的 Voronoi 区域。

(a) 证明 $V$ 是多面体，并将其表示为 $V=\{x\mid Ax\preceq b\}$。

(b) 反过来，给定一个内部非空的多面体 $P$，说明如何找到 $x_0,\ldots,x_K$，使得该多面体是 $x_0$ 关于 $x_1,\ldots,x_K$ 的 Voronoi 区域。

(c) 还可以考虑集合

$$
V_k=\{x\in\mathbf{R}^n\mid \|x-x_k\|_2\leq\|x-x_i\|_2,\ i\neq k\}.
$$

集合 $V_k$ 由 $\mathbf{R}^n$ 中那些在集合 $\{x_0,\ldots,x_K\}$ 中以 $x_k$ 为最近点的点组成。集合 $V_0,\ldots,V_K$ 给出 $\mathbf{R}^n$ 的一个多面体分解。更精确地，$V_k$ 是内部非空的多面体，$\bigcup_{k=0}^K V_k=\mathbf{R}^n$，且当 $i\neq j$ 时 $\operatorname{int}V_i\cap\operatorname{int}V_j=\emptyset$，即 $V_i$ 与 $V_j$ 至多沿边界相交。若 $P_1,\ldots,P_m$ 是内部非空的多面体，满足 $\bigcup_{i=1}^mP_i=\mathbf{R}^n$ 且 $\operatorname{int}P_i\cap\operatorname{int}P_j=\emptyset$（$i\neq j$），这个 $\mathbf{R}^n$ 的多面体分解是否可描述为由适当点集生成的 Voronoi 区域？

**2.10 二次不等式的解集。** 设 $C\subseteq\mathbf{R}^n$ 是二次不等式的解集

$$
C=\{x\in\mathbf{R}^n\mid x^TAx+b^Tx+c\leq0\},
$$

其中 $A\in\mathbf{S}^n$、$b\in\mathbf{R}^n$、$c\in\mathbf{R}$。

(a) 证明若 $A\succeq0$，则 $C$ 是凸的。

(b) 证明若对某个 $\lambda\in\mathbf{R}$ 有 $A+\lambda gg^T\succeq0$，则 $C$ 与超平面 $g^Tx+h=0$（其中 $g\neq0$）的交是凸的。上述命题的反命题是否成立？

**2.11 双曲集合。** 证明双曲集合 $\{x\in\mathbf{R}^2_+\mid x_1x_2\geq1\}$ 是凸的。作为推广，证明 $\{x\in\mathbf{R}^n_+\mid \prod_{i=1}^n x_i\geq1\}$ 是凸的。提示：若 $a,b\geq0$ 且 $0\leq\theta\leq1$，则 $a^\theta b^{1-\theta}\leq\theta a+(1-\theta)b$；见 §3.1.9。

**2.12** 下列集合哪些是凸的？

(a) 板，即形如 $\{x\in\mathbf{R}^n\mid \alpha\leq a^Tx\leq\beta\}$ 的集合。

(b) 矩形，即形如 $\{x\in\mathbf{R}^n\mid \alpha_i\leq x_i\leq\beta_i,\ i=1,\ldots,n\}$ 的集合。当 $n>2$ 时，矩形有时称为超矩形。

(c) 楔形，即 $\{x\in\mathbf{R}^n\mid a_1^Tx\leq b_1,\ a_2^Tx\leq b_2\}$。

(d) 到给定点比到给定集合更近的点集合：

$$
\{x\mid \|x-x_0\|_2\leq\|x-y\|_2\ \text{for all }y\in S\},
$$

其中 $S\subseteq\mathbf{R}^n$。

(e) 到一个集合比到另一个集合更近的点集合：

$$
\{x\mid \operatorname{dist}(x,S)\leq\operatorname{dist}(x,T)\},
$$

其中 $S,T\subseteq\mathbf{R}^n$，且

$$
\operatorname{dist}(x,S)=\inf\{\|x-z\|_2\mid z\in S\}.
$$

(f) [HUL93, volume 1, page 93] 集合 $\{x\mid x+S_2\subseteq S_1\}$，其中 $S_1,S_2\subseteq\mathbf{R}^n$ 且 $S_1$ 是凸的。

(g) 到 $a$ 的距离不超过到 $b$ 的距离的固定比例 $\theta$ 的点集合，即 $\{x\mid \|x-a\|_2\leq\theta\|x-b\|_2\}$。可假设 $a\neq b$ 且 $0\leq\theta\leq1$。

**2.13 外积的锥包。** 考虑秩为 $k$ 的外积集合 $\{XX^T\mid X\in\mathbf{R}^{n\times k},\ \operatorname{rank}X=k\}$。用简单术语描述它的锥包。

**2.14 扩张集与收缩集。** 设 $S\subseteq\mathbf{R}^n$，$\|\cdot\|$ 是 $\mathbf{R}^n$ 上的范数。

(a) 对 $a\geq0$，定义 $S_a=\{x\mid \operatorname{dist}(x,S)\leq a\}$，其中 $\operatorname{dist}(x,S)=\inf_{y\in S}\|x-y\|$。称 $S_a$ 为 $S$ 扩张或延伸 $a$。证明若 $S$ 是凸的，则 $S_a$ 是凸的。

(b) 对 $a\geq0$，定义 $S_{-a}=\{x\mid B(x,a)\subseteq S\}$，其中 $B(x,a)$ 是以 $x$ 为中心、半径为 $a$ 的范数球。称 $S_{-a}$ 为 $S$ 收缩或限制 $a$，因为 $S_{-a}$ 由所有距 $\mathbf{R}^n\setminus S$ 至少 $a$ 的点组成。证明若 $S$ 是凸的，则 $S_{-a}$ 是凸的。

**2.15 一些概率分布集合。** 设 $x$ 是实值随机变量，$\operatorname{prob}(x=a_i)=p_i$，$i=1,\ldots,n$，其中 $a_1<a_2<\cdots<a_n$。当然，$p\in\mathbf{R}^n$ 位于标准概率单纯形 $P=\{p\mid \mathbf{1}^Tp=1,\ p\succeq0\}$ 中。下列哪些条件关于 $p$ 是凸的？也就是说，对哪些条件，满足该条件的 $p\in P$ 的集合是凸的？

(a) $\alpha\leq \mathbf{E}f(x)\leq\beta$，其中 $\mathbf{E}f(x)=\sum_{i=1}^np_if(a_i)$。（函数 $f:\mathbf{R}\to\mathbf{R}$ 已给定。）

(b) $\operatorname{prob}(x>\alpha)\leq\beta$。

(c) $\mathbf{E}|x^3|\leq\alpha\mathbf{E}|x|$。

(d) $\mathbf{E}x^2\leq\alpha$。

(e) $\mathbf{E}x^2\geq\alpha$。

(f) $\operatorname{var}(x)\leq\alpha$，其中 $\operatorname{var}(x)=\mathbf{E}(x-\mathbf{E}x)^2$ 是 $x$ 的方差。

(g) $\operatorname{var}(x)\geq\alpha$。

(h) $\operatorname{quartile}(x)\geq\alpha$，其中 $\operatorname{quartile}(x)=\inf\{\beta\mid \operatorname{prob}(x\leq\beta)\geq0.25\}$。

(i) $\operatorname{quartile}(x)\leq\alpha$。

### 保持凸性的运算

**2.16** 证明若 $S_1$ 和 $S_2$ 是 $\mathbf{R}^{m+n}$ 中的凸集，则它们的部分和

$$
S=\{(x,y_1+y_2)\mid x\in\mathbf{R}^m,\ y_1,y_2\in\mathbf{R}^n,\ (x,y_1)\in S_1,\ (x,y_2)\in S_2\}
$$

也是凸的。

**2.17 透视函数下多面集的像。** 本题研究超平面、半空间和多面体在透视函数 $P(x,t)=x/t$ 下的像，其中 $\operatorname{dom}P=\mathbf{R}^n\times\mathbf{R}_{++}$。对下列每个集合 $C$，给出

$$
P(C)=\{v/t\mid (v,t)\in C,\ t>0\}
$$

的简单描述。

(a) 多面体 $C=\operatorname{conv}\{(v_1,t_1),\ldots,(v_K,t_K)\}$，其中 $v_i\in\mathbf{R}^n$ 且 $t_i>0$。

(b) 超平面 $C=\{(v,t)\mid f^Tv+gt=h\}$（$f$ 和 $g$ 不全为零）。

(c) 半空间 $C=\{(v,t)\mid f^Tv+gt\leq h\}$（$f$ 和 $g$ 不全为零）。

(d) 多面体 $C=\{(v,t)\mid Fv+gt\preceq h\}$。

**2.18 可逆线性分式函数。** 设 $f:\mathbf{R}^n\to\mathbf{R}^n$ 是线性分式函数

$$
f(x)=\frac{Ax+b}{c^Tx+d},\qquad
\operatorname{dom}f=\{x\mid c^Tx+d>0\}.
$$

假设矩阵

$$
Q=\begin{bmatrix}A&b\\c^T&d\end{bmatrix}
$$

非奇异。证明 $f$ 可逆，且 $f^{-1}$ 是线性分式映射。用 $A,b,c,d$ 给出 $f^{-1}$ 及其定义域的显式表达。提示：用 $Q$ 表示 $f^{-1}$ 可能更容易。

**2.19 线性分式函数与凸集。** 设 $f:\mathbf{R}^m\to\mathbf{R}^n$ 是线性分式函数

$$
f(x)=\frac{Ax+b}{c^Tx+d},\qquad
\operatorname{dom}f=\{x\mid c^Tx+d>0\}.
$$

本题研究凸集 $C$ 在 $f$ 下的逆像：

$$
f^{-1}(C)=\{x\in\operatorname{dom}f\mid f(x)\in C\}.
$$

对下列每个集合 $C\subseteq\mathbf{R}^n$，给出 $f^{-1}(C)$ 的简单描述。

(a) 半空间 $C=\{y\mid g^Ty\leq h\}$（$g\neq0$）。

(b) 多面体 $C=\{y\mid Gy\preceq h\}$。

(c) 椭球 $\{y\mid y^TP^{-1}y\leq1\}$，其中 $P\in\mathbf{S}^n_{++}$。

(d) 线性矩阵不等式的解集 $C=\{y\mid y_1A_1+\cdots+y_nA_n\preceq B\}$，其中 $A_1,\ldots,A_n,B\in\mathbf{S}^p$。

### 分离定理与支撑超平面

**2.20 线性方程的严格正解。** 设 $A\in\mathbf{R}^{m\times n}$、$b\in\mathbf{R}^m$，且 $b\in\mathcal{R}(A)$。证明存在 $x$ 满足

$$
x\succ0,\qquad Ax=b
$$

当且仅当不存在 $\lambda$ 满足

$$
A^T\lambda\succeq0,\qquad A^T\lambda\neq0,\qquad b^T\lambda\leq0.
$$

提示：先证明如下线性代数事实：对所有满足 $Ax=b$ 的 $x$，都有 $c^Tx=d$，当且仅当存在向量 $\lambda$ 使得 $c=A^T\lambda$、$d=b^T\lambda$。

**2.21 分离超平面的集合。** 设 $C$ 和 $D$ 是 $\mathbf{R}^n$ 中不相交的子集。考虑所有 $(a,b)\in\mathbf{R}^{n+1}$ 构成的集合，其中对所有 $x\in C$ 有 $a^Tx\leq b$，对所有 $x\in D$ 有 $a^Tx\geq b$。证明该集合是凸锥（若不存在分离 $C$ 和 $D$ 的超平面，则它是单点集 $\{0\}$）。

**2.22** 补充分离超平面定理在 §2.5.1 中的证明：证明任意两个不相交凸集 $C$ 和 $D$ 存在分离超平面。可以使用 §2.5.1 已证明的结果，即当两个集合中存在一对点，其距离等于两集合间距离时，存在分离超平面。

提示：若 $C$ 和 $D$ 是不相交凸集，则集合 $\{x-y\mid x\in C,\ y\in D\}$ 是凸的且不包含原点。

**2.23** 给出两个闭凸集的例子，它们不相交但不能被严格分离。

**2.24 支撑超平面。**

(a) 将闭凸集 $\{x\in\mathbf{R}^2_+\mid x_1x_2\geq1\}$ 表示为半空间的交。

(b) 设 $C=\{x\in\mathbf{R}^n\mid \|x\|_\infty\leq1\}$，即 $\mathbf{R}^n$ 中 $\ell_\infty$ 范数单位球，并设 $\hat{x}$ 是 $C$ 边界上的一点。显式识别 $C$ 在 $\hat{x}$ 处的支撑超平面。

**2.25 内外多面体近似。** 设 $C\subseteq\mathbf{R}^n$ 是闭凸集，且 $x_1,\ldots,x_K$ 位于 $C$ 的边界。假设对每个 $i$，$a_i^T(x-x_i)=0$ 定义了 $C$ 在 $x_i$ 处的支撑超平面，即 $C\subseteq\{x\mid a_i^T(x-x_i)\leq0\}$。考虑两个多面体

$$
P_{\mathrm{inner}}=\operatorname{conv}\{x_1,\ldots,x_K\},\qquad
P_{\mathrm{outer}}=\{x\mid a_i^T(x-x_i)\leq0,\ i=1,\ldots,K\}.
$$

证明 $P_{\mathrm{inner}}\subseteq C\subseteq P_{\mathrm{outer}}$。画图说明。

**2.26 支撑函数。** 集合 $C\subseteq\mathbf{R}^n$ 的支撑函数定义为

$$
S_C(y)=\sup\{y^Tx\mid x\in C\}.
$$

允许 $S_C(y)$ 取值 $+\infty$。设 $C$ 和 $D$ 是 $\mathbf{R}^n$ 中的闭凸集。证明 $C=D$ 当且仅当它们的支撑函数相等。

**2.27 反向支撑超平面定理。** 设集合 $C$ 是闭的、内部非空，并且在其边界每一点都有支撑超平面。证明 $C$ 是凸的。

### 凸锥与广义不等式

**2.28 $n=1,2,3$ 时的正半定锥。** 对 $n=1,2,3$，用矩阵系数和普通不等式显式描述正半定锥 $\mathbf{S}^n_+$。为描述 $\mathbf{S}^n$ 的一般元素，当 $n=1,2,3$ 时分别使用记号

$$
x_1,\qquad
\begin{bmatrix}x_1&x_2\\x_2&x_3\end{bmatrix},\qquad
\begin{bmatrix}x_1&x_2&x_3\\x_2&x_4&x_5\\x_3&x_5&x_6\end{bmatrix}.
$$

**2.29 $\mathbf{R}^2$ 中的锥。** 设 $K\subseteq\mathbf{R}^2$ 是闭凸锥。

(a) 用元素的极坐标给出 $K$ 的简单描述（$x=r(\cos\phi,\sin\phi)$，$r\geq0$）。

(b) 给出 $K^*$ 的简单描述，并画图说明 $K$ 与 $K^*$ 的关系。

(c) $K$ 何时是尖的？

(d) $K$ 何时是正则的（从而定义广义不等式）？画图说明当 $K$ 正则时 $x\preceq_Ky$ 的含义。

**2.30 广义不等式的性质。** 证明 §2.4.1 中列出的（非严格和严格）广义不等式性质。

**2.31 对偶锥的性质。** 设 $K^*$ 是凸锥 $K$ 的对偶锥，如 (2.19) 所定义。证明：

(a) $K^*$ 确实是凸锥。

(b) $K_1\subseteq K_2$ 蕴含 $K_2^*\subseteq K_1^*$。

(c) $K^*$ 是闭的。

(d) $K^*$ 的内部为 $\operatorname{int}K^*=\{y\mid y^Tx>0\ \text{for all }x\in\operatorname{cl}K\}$。

(e) 若 $K$ 内部非空，则 $K^*$ 是尖的。

(f) $K^{**}$ 是 $K$ 的闭包。（因此若 $K$ 是闭的，则 $K^{**}=K$。）

(g) 若 $K$ 的闭包是尖的，则 $K^*$ 内部非空。

**2.32** 求 $\{Ax\mid x\succeq0\}$ 的对偶锥，其中 $A\in\mathbf{R}^{m\times n}$。

**2.33 单调非负锥。** 定义单调非负锥为

$$
K_{m+}=\{x\in\mathbf{R}^n\mid x_1\geq x_2\geq\cdots\geq x_n\geq0\},
$$

即所有分量按非增顺序排列的非负向量。

(a) 证明 $K_{m+}$ 是正则锥。

(b) 求对偶锥 $K^*_{m+}$。提示：使用恒等式

$$
\begin{aligned}
\sum_{i=1}^n x_iy_i
&=(x_1-x_2)y_1+(x_2-x_3)(y_1+y_2)\\
&\quad +(x_3-x_4)(y_1+y_2+y_3)+\cdots\\
&\quad +(x_{n-1}-x_n)(y_1+\cdots+y_{n-1})+x_n(y_1+\cdots+y_n).
\end{aligned}
$$

**2.34 字典序锥与排序。** 字典序锥定义为

$$
K_{\mathrm{lex}}
=\{0\}\cup\{x\in\mathbf{R}^n\mid x_1=\cdots=x_k=0,\ x_{k+1}>0,\ \text{for some }k,\ 0\leq k<n\},
$$

即所有第一个非零分量（若存在）为正的向量。

(a) 验证 $K_{\mathrm{lex}}$ 是锥，但不是正则锥。

(b) 在 $\mathbf{R}^n$ 上定义字典序如下：$x\leq_{\mathrm{lex}}y$ 当且仅当 $y-x\in K_{\mathrm{lex}}$。（由于 $K_{\mathrm{lex}}$ 不是正则锥，字典序不是广义不等式。）证明字典序是线性序：对任意 $x,y\in\mathbf{R}^n$，要么 $x\leq_{\mathrm{lex}}y$，要么 $y\leq_{\mathrm{lex}}x$。因此，任意向量集合都可以相对于字典序锥排序，得到字典中熟悉的排序方式。

(c) 求 $K^*_{\mathrm{lex}}$。

**2.35 共正矩阵。** 若矩阵 $X\in\mathbf{S}^n$ 满足对所有 $z\succeq0$ 都有 $z^TXz\geq0$，则称 $X$ 是共正的。验证共正矩阵集合是正则锥。求其对偶锥。

**2.36 欧几里得距离矩阵。** 设 $x_1,\ldots,x_n\in\mathbf{R}^k$。由 $D_{ij}=\|x_i-x_j\|_2^2$ 定义的矩阵 $D\in\mathbf{S}^n$ 称为欧几里得距离矩阵。它满足一些显然性质，例如 $D_{ij}=D_{ji}$、$D_{ii}=0$、$D_{ij}\geq0$，以及由三角不等式得到的 $D_{ik}^{1/2}\leq D_{ij}^{1/2}+D_{jk}^{1/2}$。

现在提出问题：什么时候矩阵 $D\in\mathbf{S}^n$ 是某些 $\mathbf{R}^k$ 中点的欧几里得距离矩阵（对某个 $k$）？一个著名结果回答了该问题：$D\in\mathbf{S}^n$ 是欧几里得距离矩阵，当且仅当 $D_{ii}=0$，并且对所有满足 $\mathbf{1}^Tx=0$ 的 $x$ 有 $x^TDx\leq0$。（见 §8.3.3。）证明欧几里得距离矩阵集合是凸锥。

**2.37 非负多项式与 Hankel LMI。** 设 $K_{\mathrm{pol}}$ 是 $\mathbf{R}$ 上 $2k$ 次非负多项式（系数）集合：

$$
K_{\mathrm{pol}}
=\{x\in\mathbf{R}^{2k+1}\mid
x_1+x_2t+x_3t^2+\cdots+x_{2k+1}t^{2k}\geq0\ \text{for all }t\in\mathbf{R}\}.
$$

(a) 证明 $K_{\mathrm{pol}}$ 是正则锥。

(b) 一个基本结果表明，$2k$ 次多项式在 $\mathbf{R}$ 上非负，当且仅当它可以表示为两个次数不超过 $k$ 的多项式平方和。换言之，$x\in K_{\mathrm{pol}}$ 当且仅当多项式

$$
p(t)=x_1+x_2t+x_3t^2+\cdots+x_{2k+1}t^{2k}
$$

可以表示为

$$
p(t)=r(t)^2+s(t)^2,
$$

其中 $r$ 和 $s$ 是次数不超过 $k$ 的多项式。使用该结果证明

$$
K_{\mathrm{pol}}
=\left\{x\in\mathbf{R}^{2k+1}\ \middle|\ 
x_i=\sum_{m+n=i+1}Y_{mn}\ \text{for some }Y\in\mathbf{S}^{k+1}_+
\right\}.
$$

也就是说，$p(t)=x_1+x_2t+x_3t^2+\cdots+x_{2k+1}t^{2k}$ 非负，当且仅当存在矩阵 $Y\in\mathbf{S}^{k+1}_+$，使得

$$
\begin{aligned}
x_1&=Y_{11},\\
x_2&=Y_{12}+Y_{21},\\
x_3&=Y_{13}+Y_{22}+Y_{31},\\
&\ \vdots\\
x_{2k+1}&=Y_{k+1,k+1}.
\end{aligned}
$$

(c) 证明 $K_{\mathrm{pol}}^*=K_{\mathrm{han}}$，其中

$$
K_{\mathrm{han}}=\{z\in\mathbf{R}^{2k+1}\mid H(z)\succeq0\},
$$

且

$$
H(z)=
\begin{bmatrix}
z_1&z_2&z_3&\cdots&z_k&z_{k+1}\\
z_2&z_3&z_4&\cdots&z_{k+1}&z_{k+2}\\
z_3&z_4&z_5&\cdots&z_{k+2}&z_{k+3}\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots\\
z_k&z_{k+1}&z_{k+2}&\cdots&z_{2k-1}&z_{2k}\\
z_{k+1}&z_{k+2}&z_{k+3}&\cdots&z_{2k}&z_{2k+1}
\end{bmatrix}.
$$

这是由系数 $z_1,\ldots,z_{2k+1}$ 构成的 Hankel 矩阵。

(d) 设 $K_{\mathrm{mom}}$ 是所有形如 $(1,t,t^2,\ldots,t^{2k})$ 的向量的锥包，其中 $t\in\mathbf{R}$。证明 $y\in K_{\mathrm{mom}}$ 当且仅当 $y_1\geq0$ 且

$$
y=y_1(1,\mathbf{E}u,\mathbf{E}u^2,\ldots,\mathbf{E}u^{2k})
$$

对某个随机变量 $u$ 成立。换言之，$K_{\mathrm{mom}}$ 的元素是 $\mathbf{R}$ 上所有可能分布的矩向量的非负倍数。证明 $K_{\mathrm{pol}}=K_{\mathrm{mom}}^*$。

(e) 结合 (c) 和 (d) 的结果，推出 $K_{\mathrm{han}}=\operatorname{cl}K_{\mathrm{mom}}$。作为说明 $K_{\mathrm{mom}}$ 与 $K_{\mathrm{han}}$ 之间关系的例子，取 $k=2$ 和 $z=(1,0,0,0,1)$。证明 $z\in K_{\mathrm{han}}$、$z\notin K_{\mathrm{mom}}$。找出一个收敛到 $z$ 的 $K_{\mathrm{mom}}$ 中点的显式序列。

**2.38 [Roc70, pages 15, 61] 由集合构造的凸锥。**

(a) 集合 $C$ 的障碍锥定义为所有使得 $y^Tx$ 在 $x\in C$ 上有上界的向量 $y$ 的集合。换言之，非零向量 $y$ 属于障碍锥，当且仅当它是某个包含 $C$ 的半空间 $\{x\mid y^Tx\leq\alpha\}$ 的法向量。验证障碍锥是凸锥（不对 $C$ 作任何假设）。

(b) 集合 $C$ 的衰退锥（也称渐近锥）定义为所有满足如下性质的向量 $y$ 的集合：对每个 $x\in C$ 和所有 $t\geq0$，都有 $x-ty\in C$。证明凸集的衰退锥是凸锥。证明若 $C$ 非空、闭且凸，则 $C$ 的衰退锥是障碍锥的对偶。

(c) 集合 $C$ 在边界点 $x_0$ 处的法锥定义为所有满足对所有 $x\in C$ 有 $y^T(x-x_0)\leq0$ 的向量 $y$ 的集合（即所有定义 $C$ 在 $x_0$ 处支撑超平面的向量集合）。证明法锥是凸锥（不对 $C$ 作任何假设）。给出多面体 $\{x\mid Ax\preceq b\}$ 在其边界点处的法锥的简单描述。

**2.39 锥的分离。** 设 $K$ 和 $\tilde{K}$ 是两个凸锥，它们的内部非空且不相交。证明存在非零 $y$，使得 $y\in K^*$ 且 $-y\in\tilde{K}^*$。
