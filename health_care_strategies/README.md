OVERVIEW:

Using a dataset of 29 health plans, I compared the effectiveness of two strategies of controlling prescription
costs: (1) generic substitutes, and (2) restricting drugs that physicians can prescribe. The probability that
generic drugs lowered costs was 99.9%, compared to a 14.4% to 20% probability that restricting drugs lowered
costs. Holding all other variables constant, increasing restrictions had
negligable effects, while each additional percentage of generic drugs reduced costs by an average of $.03
(with a 95% confidence interval of $.015 to $.04) per day. Given that the average cost of
prescriptions to the plan was $1.23 per day, this is an average daily projected savings
of 2.19% (1.2% to 3.2%). However, our dataset was small, and had restricted ranges of member ages (25.2 to
32.4), percent female (49% to 57%), etc. Therefore, although I would recommend
generic drugs for reducing prescription drug costs, the level of savings cannot be accurately estimated. I
would also cautiously recommend against physician restrictions, with the caveat that more data is needed
to be certain of the strategy's ineffectiveness.

METHODOLOGY:

Exploratory analysis indicated that multicollinearity wasn't a problem, but the presence of
outliers was. I ran tests to identify outliers and deleted two. I then
separated the remaining 27 health plans into a test and training dataset, to best assess the predictive ability
of different regression models. One of the drawbacks of this method is that, in small datasets such as ours,
splitting into test and training does not always work. We have to be sure that the two sets cover about the
same space and have similar statistical properties.

I then compared different regression methods: OLS, g-prior, lasso, BMA, and ridge regressions. Of all these, the OLS regression seemed optimal, yielding the smallest MSE on the test dataset. I then identified which set of variables should be included in the model, without the risk of overfitting. The methods I used included minimizing the BIC, maximizing
the adjusted R2, and using Bayesian Model Averaging.
