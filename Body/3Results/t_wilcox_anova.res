x = Ter; y = Cock
t test:
# t = 8.5145, df = 54.382, p-value = 2.338e-11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.1791942 0.2917100
# sample estimates:
#  mean of x  mean of y 
# 0.32780823 0.09235614 
wilcox:
# W = 684, p-value = 2.078e-07
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   0.1821497 0.3051049
# sample estimates:
#   difference in location 
# 0.2477359 
anova:
1) summary
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# Group        2 0.7759  0.3880   26.73 8.51e-09 ***
# Residuals   54 0.7836  0.0145                    
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
2)TukeyHSD
#                        diff        lwr        upr     p adj
# Higher Ter-Cockroach  0.26763715  0.17794260 0.35733170 0.0000000
# Lower Ter-Cockroach   0.18824733  0.08908627 0.28740839 0.0000829
# Lower Ter-Higher Ter -0.07938982 -0.17659998 0.01782033 0.1299455

