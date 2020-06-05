
# http://scikit-image.org/docs/dev/auto_examples/features_detection/plot_blob.html

library(dplyr)
library(imager)

orig_image <- load.image("data/fish.png")
# orig_image <- load.image("data/BP3G_AON_2017_dr.jpg")

# resize image and grayscale
resize_factor <- -50 # if negative interpreted as %
my_imageRGB <- resize(orig_image, resize_factor, resize_factor)
my_image <- G(my_imageRGB)

plot(my_image)
save.image(my_image, file = "output/grayscale_image.jpg")

# now most of the birds should be much darker than the rest of the image
# set a manual threshold to turn birds black and rest white
binarized_image <- my_image < .8
binarized_image <- 1 - binarized_image
plot(binarized_image)
save.image(binarized_image, "output/processed_image.jpg")

# we work with this

# function to extract center of a blob
get.centers <- function(im,thr="99%")
{
  dt <- imhessian(im) %$% { xx*yy - xy^2 } %>% threshold(thr) %>% label
  as.data.frame(dt) %>% subset(value>0) %>% dplyr::group_by(value) %>% dplyr::summarise(mx=mean(x),my=mean(y))
}

# focus in on a small patch of the image, this should smooth over noise
small_patch <- imsub(1-my_image, x %in% 0:200, y %in% 0:200)
plot(small_patch)

# blur the patch so we can use the blob detector, need to experiment
# with the sigma parameter, depends on size of small_patch and of 
# object you want to detect.
blurry_patch <- isoblur(small_patch, sigma = 4)
plot(blurry_patch)

# get the blob centers, need to experiment with threshold to get it right
centers <- get.centers(blurry_patch, thr = "97%") 

# does poorly at edges but that's ok because we'll shift the patch along
# bit by bit. So remove centers found at the margins, if they're real
# they will be picked up in other patches.
centers <- centers %>% dplyr::filter(mx > 3, mx < 197, my > 3, my < 197)

# plot the estimated centers on top of the image, and check
plot(imsub(my_imageRGB, x %in% 0:200, y %in% 0:200))
points(centers$mx,centers$my,col="red",pch=20,cex=1)

