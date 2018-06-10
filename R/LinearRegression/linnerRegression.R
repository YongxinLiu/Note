# 加载所需的包
library(ggplot2)
library(car)
library(caret)
library(corrplot)

# 读取数据
data(mtcars)
str(mtcars)
head(mtcars)
summary(mtcars)

# 数据准备
# 确保分类变量存储为因子。 在下面的程序中，我们将变量转换为因子。
mtcars$am <- as.factor(mtcars$am)
mtcars$cyl <- as.factor(mtcars$cyl)
mtcars$vs <- as.factor(mtcars$vs)
mtcars$gear <-  as.factor(mtcars$gear)

# 识别和修正共线性
# 在这一步中，我们正在识别彼此高度相关的自变量。 由于mpg是一个因变量，我们将在下面的代码中删除它。
mtcars_a <- subset(mtcars, select = -c(mpg))
# 筛选数值变量
numericData <- mtcars_a[sapply(mtcars_a, is.numeric)]

# 统计相关变绘图
descrCor <- cor(numericData)
print(descrCor)
corrplot(descrCor, order = "FPC", method = "color",
         type = "lower", tl.cex = 0.7, tl.col = rgb(0,0,0))

# 筛选相关高于0.7的组，并删除
highlyCorrelated <- findCorrelation(descrCor, cutoff=0.7)
highlyCorCol <- colnames(numericData)[highlyCorrelated]
dat3 <- mtcars[, -which(colnames(mtcars) %in% highlyCorCol)]
dim(dat3)
# 有三个变量“hp”“disp”“wt”被发现是高度相关的。 我们已经删除它们以避免共线。 现在，我们有7个独立变量和1个因变量。


# 开发回归模型
# 在这一步，我们正在建立多元线性回归模型。
fit <- lm(mpg ~ ., data=dat3)
# 查看线性回归模型的系数和ANOVA表
summary(fit)
summary(fit)$coeff
anova(fit)
# 性回归模型测试估计等于零的零假设。 具有小于0.05的p值的独立变量意味着我们拒绝5％显着性水平的零假设。 这意味着该变量的系数不等于0.大p值意味着变量对预测目标变量没有意义。
par(mfrow=c(2,2))
plot(fit)
# 线性回归模型的关键图如下所示 :-
# 残留物与拟合值
# 正常Q-Q
# 缩放位置
# 残差与杠杆


# 计算模型性能度量
summary(fit)$r.squared
summary(fit)$adj.r.squared
AIC(fit)
BIC(fit)
# 更高的R平方和调整的R平方值，更好的模型。 然而，更低的AIC和BIC得分，更好的模型。
# AIC和BIC是拟合度的衡量标准。 他们惩罚复杂的模型。 换句话说，它会惩罚更多的估计参数。 它相信一个概念，即一个具有较少参数的模型将被一个具有更多参数的模型要好。 一般来说，BIC比AIC更为免费参数惩罚模型。 两个标准都取决于估计模型的似然函数L的最大值。

# 变量选择方法
# 有三种变量选择方法 - 向前，向后，逐步。
# 1.以单个变量开始，然后基于AIC（“前进”）一次添加一个变量
# 2.从所有变量开始，基于AIC（’后退’）迭代地去除那些重要性低的变量
# 3.双向运行（’逐步’）
library(MASS)
step <- stepAIC(fit, direction="both")
summary(step)
step <- stepAIC(fit, direction="forward")
summary(step)
n <- dim(dat3)[1]
stepBIC <- stepAIC(fit,k=log(n))
summary(stepBIC)
# 在基于BIC执行逐步选择之后，查看以上估计值。 变量已经减少，但调整后的R-Squared保持不变（稍微改善）。 AIC和BIC分数也下降，这表明一个更好的模型。
AIC(stepBIC)
BIC(stepBIC)

# 计算标准化系数
# 标准化系数有助于根据标准化估计值的绝对值排列预测值。 值越高，变量越重要。
#使用QuantPsyc包的lm.beta函数计算
library(QuantPsyc)
lm.beta(stepBIC)
#自定义函数计算
stdz.coff <- function (regmodel){
  b <- summary(regmodel)$coef[-1,1]
  sx <- sapply(regmodel$model[-1], sd)
  sy <- sapply(regmodel$model[1], sd)
  beta <- b * sx / sy
  return(beta)
}
std.Coeff <- data.frame(Standardized.Coeff = stdz.coff(stepBIC))
std.Coeff <- cbind(Variable = row.names(std.Coeff), std.Coeff)
row.names(std.Coeff) <- NULL

# 计算方差膨胀因子（VIF）
# 与独立变量高度相关的情况相比，差异膨胀因子衡量的是系数的变化幅度。 它应该小于5。
vif(stepBIC)

# 测试其它假设 Autocorrelation Test
durbinWatsonTest(stepBIC)
#Normality Of Residuals (Should be > 0.05)
res=residuals(stepBIC,type="pearson")
shapiro.test(res)

#Testing for heteroscedasticity (Should be > 0.05)
ncvTest(stepBIC)

#Outliers – Bonferonni test
outlierTest(stepBIC)

#See Residuals
resid = residuals(stepBIC)

#Relative Importance
library(relaimpo)
calc.relimp(stepBIC)

# 查看实际值和预测值 See Predicted Value
pred <- predict(stepBIC,dat3)
#See Actual vs. Predicted Value
finaldata <-  cbind(mtcars,pred)
print(head(subset(finaldata, select = c(mpg,pred))))

# 其它有用的函数
#Calculating RMSE
rmse = sqrt(mean((dat3$mpg - pred)^2))
print(rmse)

#Calculating Rsquared manually
y = dat3[,c("mpg")]
R.squared = 1 - sum((y-pred)^2)/sum((y-mean(y))^2)
print(R.squared)

# Calculating Adj. Rsquared manually
n = dim(dat3)[1]
p = dim(summary(stepBIC)$coeff)[1] - 1
adj.r.squared = 1 - (1 - R.squared) * ((n - 1)/(n-p-1))
print(adj.r.squared)

#Box Cox Transformation
library(lmSupport)
modelBoxCox(stepBIC)

# K-fold交叉验证
# 在下面的程序中，我们正在进行5倍交叉验证。 在5倍交叉验证中，数据被随机分成5个相同大小的样本。 在5个样本中，随机20％数据的单个样本保留为验证数据，其余80％用作训练数据。 然后这个过程重复5次，5个样本中的每一个都只用作验证数据一次。 稍后我们将结果平均。
library(DAAG)
kfold = cv.lm(data=dat3, stepBIC, m=5)
