本次作业用高斯消元法和LU分解进行了矩阵求逆（暂时不知道怎么用nrulib里面的子程序来做SVD）。
	具体做法是封装了一个类Matrix，里面包含了创建矩阵creat()，展示矩阵display()及各种求逆矩阵的函数，可以直接调用。（矩阵维数为10以上不显示矩阵）
	测试结果如下。在较低维度时（小于500），用时相近，但随着维度增加，LU分解所用的时间比高斯消元法所用的时间明显更多。这里有可能是我写的程序并没有最优。但我认为LU分解的时间确实比高斯消元法多。
