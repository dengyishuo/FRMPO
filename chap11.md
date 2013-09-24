## 11.5 R包概要

### 11.5.1 DEoptim和RcppDE

很明显，在数值计算的观点下，上文展示的投资组合优化问题，特别是对一个资产配置(asset allocation)的ES的等风险贡献问题(the equal-risk contributions to the ES)，寻找一个全局最优解是很复杂的。有鉴于此，同时也因为 **PortfolioAnalysis** (见11.5.3小节）对**DEoptim**的依赖性，我们需要分开展示实现差异进化算法(differential evolution)的R包。当然，实现者不只是解决了给定置信水平下(confidence level)，对ES等边际贡献的投资组合配置的优化问题，还能解决别的优化问题。

差异进化算法（DE）是有Stron和Ulrich(1997)年提出的，Price等人(2006)详细论述了它的应用。它被归类为一种遗传算法，因此是一种无导数的全局优化(derivative-free global optimizer)。从某种意义上来说， DE也是一种进化算法，其中候选解的初始种群接受变异和选择操作从而产生目标函数的全局最小值。因为DE算法不需要导数，所以对目标函数和参数只要求它是实值的，不要求目标函数可导或连续。因为初始种群是随机生成的，所以要重复实验结果，就需要给随机数生成器选定种子。对遗传算法的详细解释可以参考Mitchell(1998)。

历史上，CRAN第一个实现DE算法的R包是Mullen等人(2011)的DEoptim。这个包在CRAN的优化任务视图(Optimization Task View)中。DEoptim的早期版本中，DE算法部分完全用R实现（由David Ardia编写）。为了提高运算速度，后来的版本中，DE过程用C编写，类似于Price等人(2006)的MS Visual C++实现。这个包使用了S3类和方法，附带两个简介短文(vignettes)。在其中一篇简介短文中，Ardia等人（2010）给出了一个例子，展示了用DEoptim()解决大规模投资组合优化问题。这个例子中，如果使用了梯度优化方法，DEoptim()会失效。

DEoptim()的参数列表需要制定四个参数，此外还有省略号参数(ellipsis argument)。参数fn表示需要最小化的目标函数。这个函数必须返回一个实数值。参数的下界和上界可以通过参数lower和upper设定。目标函数中必须以惩罚项的形式指定约束(constraint)而不是箱约束(box constraint)。约束可以直接以可加项的形式出现在目标函数中。每个约束都会乘以一个正常数因子。违反限制会增大要求最小值的函数的值。函数DEoptim.control()返回的列表对象(list object)可以赋值给参数control，对算法进行修改。

DEoptim()返回一个S3类的列表对象DEoptim，summary()和plot()函数可以用于这个对象。这个列表对象由两个命名元素optim和member组成。这两个元素也是列表对象。前者包括优化解(bestmem)、目标函数在最优解上的值(bestval)、进化的代数(nfeval)和迭代的次数(iter)。后者包括指定的下界和上界(lower, upper)、每次迭代中目标函数的最小值(bestvalit)和相关的参数值(bestmenit),此外还有最后一次迭代中生成的种群(pop)和一个保存中间种群的列表(storepop)。如上面指出的一样，使用者可以通过DEoptim.control()返回的列表对象调整算法。其中，优化过程可以通过为目标函数设定一个特定的值来终止（参数VTR），这个值默认是-Inf，如果修改成更高的值，优化过程会在迭代次数达到最大值（可以通过itermax修改，默认值是200）或者目标函数值小于VTR的时候终止。此外，使用者可以设置相对收敛限度(relative convergence tolerance,参数reltol)来控制优化过程的终止。根据**R**运行的平台不同，这个值大约是1-e8。此外，使用者可以用参数strategy选择在六种种群繁衍的策略中的一种。这些策略的细节可以在包手册和简介短文中找到，它们基于Price等人(2006)的详细说明。参数NP决定了种群大小。如果不改变它的默认值NA，那么程序会生成十组随机的参数向量。另一方面，使用者可以为initpot提供一个矩阵对象(matrix object)确定初始种群的参数。这个矩阵的行对应目标函数的参数，列对应种群。参数storepopfrom和stroepopfreq用于控制迭代间处理参数的过程。前者设置从哪一代开始保存中间种群，默认设置下不保存中间种群。后者决定保存种群的频率，默认值为1，每一代都会保存。参数F的值决定每次迭代的步长。参数列表中其余的参数与继承的参数有关，有些也和选择的繁衍策略的有关。读者可以参阅用户手册获得更多关于这些控制的细节。最后，优化的过程可以通过参数trace来监控，如果这个值为TRUE（默认），每次迭代的中介结果都会输出，如果设定了一个整数值i，那么每i次迭代之后会输出结果一次。

**RcppDE**和**DEoptim**很相似（参看Eddelbuettel 2012)。这个包也被包括于CRAN的优化任务视图中。这个两个包的主要不同在于DE算法的接口。包如其名，RcppDE与C++代码接合。在简介短文与基于C实现的DEoptim的比较中，RcppDE在各种优化任务中运行得更快。与DEoptim包不同，RcppDE中的DEoptim()函数有一个额外的参数env。这个参数用于确定目标函数运行的环境。如果不指定，fn的函数会在一个新环境中运行。对类属性DEoptim有定义的函数仍然可用，DEoptim.control()也是。

###11.5.2 FRAPO包###
这个小节介绍了**FRAPO**包中与本章主题有关的函数。

函数dr()、cr()和rhow()可以用于度量给定权重向量(weighted vector)和方差-协方差矩阵(variacne-covariance matrix)的多样化程度(the degree of diversification)。这些函数返回11.2节引入的多样化比率(diversification ratio),集中度(concentraition rate)和波动性加权平均相关(volatility-weighted average correlation)。这些函数都有相同的参数列表，包括权重向量(weight)和方差-协方差矩阵(Sigma)。最多样化投资组合(the most diversified portfolio)可以用PMD()计算，返回的阵列(array)由参数Returns确定。逻辑值参数percentage确定返回的权重向量是以百分比表示（默认）还是以十进制小数表示。省略号参数会传给cov()，使用者可以控制方差-协方差矩阵的计算。结果可以通过函数cov2cor()转换成相关性矩阵。PMD返回一个常规的S4类PortSol，对于这个类，show()、weight()、solution()和update()都是可用的。

给定资产配置对风险的边际贡献(the marginal congtribution)和散度矩阵由函数mrc()返回。除了权重向量(weight)和方差-协方差矩阵(Sigma)之外，使用者还可以指定贡献是以百分比形式或是十进制小数形式返回。函数PERC()计算等风险投资组合的解，这个函数有三个参数：Sigma，方差-协方差矩阵；par，初始权重向量；省略号参数，传递给实施优化的nlminb()函数。除了目标值和初始值之外，优化过程中还设定了lower=0，upper=1。

最后，函数tdc()和PMTD()可以用于构造优化尾部依赖投资组合(optimal tail-dependent portfolio)。前者返回二元的TDC，使用者可以用method参数确定使用经验尾部关联(empirical tail copula)或是稳定尾部函数(stable tail function)来进行非参数估计。使用者必须提供返回的阵列，因为参数x值要求能转换成矩阵对象。返回下尾部依赖系数(lower tail dependent coefficient)或(upper tail dependent coefficient)由逻辑值参数lower确定，这个参数的默认值为TRUE。阈值k默认设为NULL，如果不改变这个值，程序会使用样本数的平方根。省略号参数会被传到rank()函数中，使用者可以指定是否计算域数据相关的序统计量。【The function returns a matrix object with the bivariate TDCs and 1s on the diagonal as by definition a series has a dependence of unity with itself.】对一个做多(long-only)的投资者,最小尾部依赖投资组合(the minimum tail-dependent portifolio)的解通过调用PMTD()函数得到。这个函数总体与PGMV()相似，返回一个全局最小方差投资组合(a global minimum-variance portfolio)。但它使用tdc()返回的TDC矩阵代替方差协方差矩阵来度量离差(dispersion)。除了以上参数之外，还需为参数Returns指明返回的阵列。PMTD返回正规的S4类PortSol。上面提到的函数对它仍然可用。

###PortfolioAnalytics包###
**PortfolioAnalytics**包求解投资组合优化的方法和之前提到的包有着概念性的不同。它可以让使用者获得复限制或/与复目标函数的数值解。再比如，对这章中展示的风险度量相关方法(the risk measure related method)，我们可以直接用一个返回问题中的风险度量的函数作为目标函数。数值优化解由**DEoptim**（参看 Mullen 等人 2011年的文章)包中的差异进化算法优化器得到，或者通过随机生成符合某些特殊限制的投资组合得到。前一种算法由Storn和Ulrich（1997）提出。一般的方法是依限制函数修改目标函数，限制函数会被乘上一个惩罚因子。这个包的主页在R-Forge上，附带一篇简介短文。这篇短文说明了它是如何优化一个带CVaR预算的投资组合的。

这个包中，最基础的三个函数分别是： constraints(), 定义投资组合权重上的限制； add.object()，将constraints()返回的限制的列表加到目标函数上； optimize.portfolio()， 计算出投资组合问题的数值解。 constraints()函数返回一个非正式S3类constraint的列表对象。is()函数和update()函数对于这种类型的对象（constraint对象)是可用的。optimize.portfolio()函数返回非正式S3类optimize.portfolio.DEoptim、optimize.portfolio.random或optimize.portfolio的对象。plot()和extractStats()可以用于这些对象。extractStats()函数返回优化投资组合解的统计量。给定一组权重下目标函数的值可调用constrained_objective()函数计算得。有两个实验的函数值得特别说明一下：optimize.portfolio.rebalancing()，能以特定的频率计算投资组合再平衡； optimize.portfolio.parallel()在多核计算机上多线程调用optimize.portfolio()【这句需推敲】，该函数默认设定使用4个节点。包中其他没有介绍的函数主要与作图和优化投资组合的解中查询/提取信息有关。

