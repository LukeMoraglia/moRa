data("iris")

col4o <- c(rep("red", 50), rep("blue", 50), rep("green", 50))
col4g <- c(setosa = 'red', versicolor = 'blue', virginica = 'green')

moRa::mora_pca(iris[,1:4], iris[,5], col4obs = col4o, col4group = col4g, want34 = TRUE)
