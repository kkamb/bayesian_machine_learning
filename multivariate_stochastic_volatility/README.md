This analysis explores financial markets risk through a latent factor stochastic volatility
model, in order to explore its underlying structure and interdependencies. I look at
seven different markets – the S&P500, the US Bond index, US 10yr treasury, Japanese
Yen, Emerging Markets Index, oil, and gold. I first derive a Bayesian latent factor
model based on the weekly returns for a ten year period, and find that the seven markets
can be characterized by three underlying factors. I then use the factor loadings to derive
a volatility model to explore each market’s idiosyncratic volatilities, as well as the three
common underlying factor volatilities. I use Gibbs sampling and FFBS algorithms for
the derivations and inferences. My results indicate that, while straightforward extensions
of this model could effectively capture the underlying structure of well-functioning markets,
they are not robust at times of market failure.

To run this analysis, download all the files in one directory, and run 
<a href="https://github.com/kkamb/bayesian_machine_learning/blob/master/multivariate_stochastic_volatility/projectmultivariate.r">projectmultivariate.r</a>.
