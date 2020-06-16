# data for testing matrix_plotter

data("iris")
data("mtcars")
data("Harman74.cor")

cor_iris <- cor(iris[,1:4])
cor_cars <- cor(mtcars)
cor_har <- Harman74.cor$cov


rand_corrmat <- function(n){
   cor(purrr::map_dfc(1:n, ~ sample(1:10)) %>%
       stats::setNames(paste0("v", 1:n)))
}

cor5 <- rand_corrmat(5)
cor10 <- rand_corrmat(10)
cor20 <- rand_corrmat(20)
cor50 <- rand_corrmat(50)
cor60 <- rand_corrmat(60)
cor75 <- rand_corrmat(75)
cor100 <- rand_corrmat(100)
cor1000 <- rand_corrmat(1000)
