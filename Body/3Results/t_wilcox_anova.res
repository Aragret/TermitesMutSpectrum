x = Ter; y = Cock
t test:
# t = 8.5145, df = 54.382, p-value = 1.392e-11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.1847356 0.2985031
# sample estimates:
#  mean of x  mean of y 
# 0.32780823 0.08618889 
wilcox:
# W = 684, p-value = 1.581e-07
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   0.1905307 0.3128353
# sample estimates:
#   difference in location 
# 0.2542703 
anova:
1) summary
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# Group        2 0.8141  0.4071   27.81 4.98e-09 ***
# Residuals   54 0.7904  0.0146                    
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
2)TukeyHSD
#                        diff        lwr        upr     p adj
# Higher Ter-Cockroach  0.27380441  0.18372103 0.36388778 0.0000000
# Lower Ter-Cockroach   0.19441459  0.09482367 0.29400551 0.0000531
# Lower Ter-Higher Ter -0.07938982 -0.17702138 0.01824174 0.1321502

