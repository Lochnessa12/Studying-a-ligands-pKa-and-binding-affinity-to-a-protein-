#_______________________ Part I - pH-dependence of Neutral red __________________________

# Import Data
x_pH <- c(5.0,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12)
y_Abs <- c(0.231,0.224,0.151,0.077,0.033,0.021,0.03,0.029,0.027)
plot(x_pH, y_Abs, type="b",ylab="Absorbance", xlab="pH", main="Absorbance of Neutral Red at different pH
     (wavelength 532 nm)")

#Estimate constants: A= Abs at low pH, B= Abs at high pH, Read off pKa 
A <-max(y_Abs)
B <- min(y_Abs)
pKa <-7
  
#Non-linear regression (Nonlinear Least Sqaures, NLS) 
model <- nls(y_Abs ~ ((a*10^-x_pH + b*10^-c)/(10^-c+10^-x_pH)), start = list (a=A, b=B, c=pKa))
lines(x_pH, predict(model), col='red')

legend("topright",
       c("Non-linear fitting"), 
       col=c("Red"), lty=c(1))

summary(model)


#Linear regression
HA <- c(A - y_Abs)
A_Base <- c(y_Abs - B)

Log_y <- log10(HA/A_Base)
Log_y

plot(x_pH, Log_y,ylab="Log(Absorbance)", xlab="pH",type="b",main="Logarithmic absorbance of Neutral Red at different pH
     (wavelength 532 nm)")

#Choose the subset of the data that shows a linear tendency. Leave out Inf values. The chosen values are value 2-5.
SubLog_y <- Log_y[c(2,3,4,5)]
Subx_pH <- x_pH[c(2,3,4,5)]

plot(Subx_pH, SubLog_y,ylab="Log(Absorbance)", xlab="pH",type="b", main="Logarithmic absorbance for Neutral Red at different pH (linear tendency), w. linear fitting
     (wavelength 532 nm)")

#Perform the linear regresion (Linear models, lm)
Log_model <- lm(SubLog_y ~ Subx_pH)
lines(Subx_pH, predict(Log_model), col='green')
legend("topleft",
       c("Linear fitting"), 
       col=c("Green"), lty=c(1))

summary(Log_model)

#_______________________ Part II - Protein Binding __________________________

#Import data. Valgt bølgelængde 550
X_Conc <- c(0.0,3.07,6.08,9.03,11.92,14.76,17.55,20.28,22.96,25.60,28.18,30.72,33.21)
Y_dA <- c(0.0127,0.0165,0.0269,0.032,0.0343,0.0345,0.0371,0.035,0.0363,0.0402,0.0425,0.043,0.0369)
plot(X_Conc, Y_dA, type="b", ylab="Absorbance", xlab="Concentration (Riboflavin binding protein) [uM]", main="Absorbance for Neutral Red in different Riboflavin binding protein concentrations
     (wavelength 550 nm)")

#Define dAmax and enter an estimated value for k based on the plot
dAmax <- max(Y_dA)
k <-6.5

#Non-linear regression on Equation 20 from the lab manual
kdmodel <- nls(Y_dA ~ ((dAmax * X_Conc)/(k + X_Conc)), start = list (dAmax = dAmax, k = k))
lines(X_Conc, predict(kdmodel), col='red')
legend("topleft",
       c("Non-linear fitting"), 
       col=c("Red"), lty=c(1))
summary(kdmodel)


