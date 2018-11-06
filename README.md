# White-Box Clustering Algorithm Design

This is a framework for representative based component based cluster algorithm design. Representative based cluster algorithms are one of the largest class of clustering algorithms (K-Means as a typical example) which can be decomposed to RCs (Reusable Components). Identified RCs are:
* (Optional) Normalization,
* Initialization of representatives,
* Measuring distance (and assigning of data points to representatives),
* Updating representatives,
* Cluster Evaluation.

All of the RCs have different implementation as representative based clustering algorithms evolved. Every improvement tried to solve major issues of K-Means algorithm just by improving one of the RC leaving other unchanged. This R package aims to develop all of those improvemets in such manner that one can extend it by implementing new solution for specific RC and use it with already developed solutions of other RCs.

## Getting Started

Hopefully, one will find all needed instructions for running and developing improvements for _whiboclustering_ R package.

### Installing

Installment of R package should be as easy as every other package.

```
install.package("whiboclustering")
```

After that you can use it by running command:

```
library("whiboclustering")
```

### Using the package

Once you installed and loaded the library you can use it for your clustering experiments. 

Simple example:

```
data <- iris[, 1:4] #Selecting only numeric columns

model <- whiboclustering(data = data, k = 3)
```

With this step you created WhiBo Cluster model which will show information about clustering. You can access brief information about cluster model by providing:

```
print(model)
```

### Implementations of RCs

Here one can find implemented solutions for RCs. Formulas and more detailed explanation will be available on GitHub page (soon to be). If you know any other RC solution that is available, but not implemented please contact us.

Normalization:
* No normalization
* Z transformation
* $L_{2}$ normalization
* $L_{1}$ normalization
* $L_{\infty}$ normalization
* Max-Min normalization
* Mean normalization
* Logarithmic normalization
* Non-monotonic normalization
* Comprehensive normalization
* Decimal Scaling
* Sigmoid transformation
* Softmax transformation

Cluster Initialization:
* Random
* Forgy
* K-Means++
* KKZ
* PCA
* AGNES
* DIANA
* Ward
* Quantile
* CCIA

Measuring Distance and Cluster Assignment:
* Euclidean
* Squared Euclidean
* Manhattan
* Canbera
* Chebyshev
* Cosine
* Correlation
* Sorensen
* Soergel
* Kulczynski
* Lorentzian
* Gower
* Inersection
* Czekanowski
* Motika
* Ruzicka
* Tanimoto
* Inner Product
* Jaccard
* Dice
* Fidelity
* Bhattacharyya
* Hellinger
* Whittaker

Update (Recalculate) Cluster Representatives:
* Mean
* Median
* Trimmed mean
* Geometric mean
* Harmonic mean
* Quadratic mean
* Trimean
* Midhinge
* Midrange
* Online Mean
* Online Median
* Online Trimmed mean
* Online Geometric mean
* Online Harmonic mean
* Online Quadratic mean
* Online Trimean
* Online Midhinge
* Online Midrange

Cluster evaluation:
* Total Sum of Squares
* Within Sum of Squares
* Between Sum of Squares
* Ball-Hall index
* Banfeld-Raftery index
* C index
* Calinski-Harabasz index
* Davies-Bouldin score
* Det ratio
* Dunn index
* Gamma index
* G+ index
* Silhouette score
* Xie-Beni index

## Extending _whiboclustering_

This part is explained in more details on GitHub page (soon to be) or you can contact us.

## Authors

* **Sandro Radovanovic** [LinkedIn] (https://www.linkedin.com/in/sandroradovanovic/)
* **Milan Vukicevic** [LinkedIn] (https://www.linkedin.com/in/milan-vukicevic-1a60685/)

## Acknowledgements

We had a help from other members:
* Milija Suknovic
* Boris Delibasic [ResearchGate] (https://www.researchgate.net/profile/Boris_Delibasic/info)
* Milos Jovanovic [LinkedIn] (https://www.linkedin.com/in/harmonija/)
