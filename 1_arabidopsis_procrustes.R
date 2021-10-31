#require(devtools)
#install_version("shapes", version = "1.1-10", repos = "http://cran.us.r-project.org")

#install.packages("./shape", repos = NULL, type="source")


#Read in the 'shapes' library
library(shapes)

#SPECIFY LANDMARK NUMBER (K), LANDMARK DIMENSIONS (M), and NUMBER OF SAMPLES (N)
#########

k <- 100
m <- 2
n <- 11

#READ IN DATA WITH ASSOCIATED DIMENSION PARAMETERS
#########

data <- read.in("./interpolated_points.txt",k,m)

#PERFORM PROCRUSTES ANALYSIS, NO SCALING, ALLOW REFLECTION
########

GPA <- procGPA(data, scale=FALSE, reflect=TRUE)

#SAVE PROCRUSTES COORDINATES FOR EACH LEAF AS A SEPARATE FILE
########

leaf1 <- as.data.frame(GPA$rotated[,,1])
leaf2 <- as.data.frame(GPA$rotated[,,2])
leaf3 <- as.data.frame(GPA$rotated[,,3])
leaf4 <- as.data.frame(GPA$rotated[,,4])
leaf5 <- as.data.frame(GPA$rotated[,,5])
leaf6 <- as.data.frame(GPA$rotated[,,6])
leaf7 <- as.data.frame(GPA$rotated[,,7])
leaf8 <- as.data.frame(GPA$rotated[,,8])
leaf9 <- as.data.frame(GPA$rotated[,,9])
leaf10 <- as.data.frame(GPA$rotated[,,10])
leaf11 <- as.data.frame(GPA$rotated[,,11])

#CENTER THE BASE OF EACH LEAF TO THE ORIGIN
#THE BASE OF EACH LEAF IS THE MEAN OF THE FIRST AND LAST POINTS
########

leaf1[,1] <- leaf1[,1] - mean( leaf1[1,1], leaf1[100,1] )
leaf1[,2] <- leaf1[,2] - mean( leaf1[1,2], leaf1[100,2] )

leaf2[,1] <- leaf2[,1] - mean( leaf2[1,1], leaf2[100,1] )
leaf2[,2] <- leaf2[,2] - mean( leaf2[1,2], leaf2[100,2] )

leaf3[,1] <- leaf3[,1] - mean( leaf3[1,1], leaf3[100,1] )
leaf3[,2] <- leaf3[,2] - mean( leaf3[1,2], leaf3[100,2] )

leaf4[,1] <- leaf4[,1] - mean( leaf4[1,1], leaf4[100,1] )
leaf4[,2] <- leaf4[,2] - mean( leaf4[1,2], leaf4[100,2] )

leaf5[,1] <- leaf5[,1] - mean( leaf5[1,1], leaf5[100,1] )
leaf5[,2] <- leaf5[,2] - mean( leaf5[1,2], leaf5[100,2] )

leaf6[,1] <- leaf6[,1] - mean( leaf6[1,1], leaf6[100,1] )
leaf6[,2] <- leaf6[,2] - mean( leaf6[1,2], leaf6[100,2] )

leaf7[,1] <- leaf7[,1] - mean( leaf7[1,1], leaf7[100,1] )
leaf7[,2] <- leaf7[,2] - mean( leaf7[1,2], leaf7[100,2] )

leaf8[,1] <- leaf8[,1] - mean( leaf8[1,1], leaf8[100,1] )
leaf8[,2] <- leaf8[,2] - mean( leaf8[1,2], leaf8[100,2] )

leaf9[,1] <- leaf9[,1] - mean( leaf9[1,1], leaf9[100,1] )
leaf9[,2] <- leaf9[,2] - mean( leaf9[1,2], leaf9[100,2] )

leaf10[,1] <- leaf10[,1] - mean( leaf10[1,1], leaf10[100,1] )
leaf10[,2] <- leaf10[,2] - mean( leaf10[1,2], leaf10[100,2] )

leaf11[,1] <- leaf11[,1] - mean( leaf11[1,1], leaf11[100,1] )
leaf11[,2] <- leaf11[,2] - mean( leaf11[1,2], leaf11[100,2] )

#CHECK THAT LEAVES ARE PROCRUSTES ALIGNED AND CENTERED ON ORIGIN
########

# Important note:
# Leaves were aligned on their side
# By transposing data and changing sign of the y values
# leaves can be reoriented to be upright
# data will be saved transposed and -y values so leaves are upright oriented

library(ggplot2)

p <- ggplot(leaf1, aes(V2, -V1))
p + geom_point() +
geom_point(data=leaf2, aes(V2, -V1))+
geom_point(data=leaf3, aes(V2, -V1))+
geom_point(data=leaf4, aes(V2, -V1))+
geom_point(data=leaf5, aes(V2, -V1))+
geom_point(data=leaf6, aes(V2, -V1))+
geom_point(data=leaf7, aes(V2, -V1))+
geom_point(data=leaf8, aes(V2, -V1))+
geom_point(data=leaf9, aes(V2, -V1))+
geom_point(data=leaf10, aes(V2, -V1))+
geom_point(data=leaf11, aes(V2, -V1))+
coord_fixed()

#RECREATE DATA IN TRANSPOSED FORM, ALIGNED AT LEAF BASES
########

leaf1 <- cbind(leaf1$V2, -leaf1$V1)
leaf2 <- cbind(leaf2$V2, -leaf2$V1)
leaf3 <- cbind(leaf3$V2, -leaf3$V1)
leaf4 <- cbind(leaf4$V2, -leaf4$V1)
leaf5 <- cbind(leaf5$V2, -leaf5$V1)
leaf6 <- cbind(leaf6$V2, -leaf6$V1)
leaf7 <- cbind(leaf7$V2, -leaf7$V1)
leaf8 <- cbind(leaf8$V2, -leaf8$V1)
leaf9 <- cbind(leaf9$V2, -leaf9$V1)
leaf10 <- cbind(leaf10$V2, -leaf10$V1)
leaf11 <- cbind(leaf11$V2, -leaf11$V1))

#SAVE DATA IN FORM OF SEPARATE SETS OF XVALS AND YVALS
#WHERE COLS ARE XVAL/YVAL 1 THRU 100
#AND ROWS ARE 1-11 FOR LEAVES 1-11
########

xvals_for_modeling <- matrix(nrow=11, ncol=100)

for(i in 1:100) {

	print(i)
	
	xvals_for_modeling[1,i] <- leaf1[i,1]
	xvals_for_modeling[2,i] <- leaf2[i,1]
	xvals_for_modeling[3,i] <- leaf3[i,1]
	xvals_for_modeling[4,i] <- leaf4[i,1]
	xvals_for_modeling[5,i] <- leaf5[i,1]
	xvals_for_modeling[6,i] <- leaf6[i,1]
	xvals_for_modeling[7,i] <- leaf7[i,1]
	xvals_for_modeling[8,i] <- leaf8[i,1]
	xvals_for_modeling[9,i] <- leaf9[i,1]
	xvals_for_modeling[10,i] <- leaf10[i,1]
	xvals_for_modeling[11,i] <- leaf11[i,1]
}


yvals_for_modeling <- matrix(nrow=11, ncol=100)

for(i in 1:100) {

	print(i)
	
	yvals_for_modeling[1,i] <- leaf1[i,2]
	yvals_for_modeling[2,i] <- leaf2[i,2]
	yvals_for_modeling[3,i] <- leaf3[i,2]
	yvals_for_modeling[4,i] <- leaf4[i,2]
	yvals_for_modeling[5,i] <- leaf5[i,2]
	yvals_for_modeling[6,i] <- leaf6[i,2]
	yvals_for_modeling[7,i] <- leaf7[i,2]
	yvals_for_modeling[8,i] <- leaf8[i,2]
	yvals_for_modeling[9,i] <- leaf9[i,2]
	yvals_for_modeling[10,i] <- leaf10[i,2]
	yvals_for_modeling[11,i] <- leaf11[i,2]
}

#WRITE OUT DATA
########

write.csv(xvals_for_modeling, "xvals_for_modeling.csv")
write.csv(yvals_for_modeling, "yvals_for_modeling.csv")


