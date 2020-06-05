
# http://scikit-image.org/docs/dev/auto_examples/features_detection/plot_blob.html

library(dplyr)
library(imager)

orig_image <- load.image("data/gentoo.jpg")

# resize image and grayscale
resize_factor <- -50 # if negative interpreted as %
my_image <- resize(orig_image, resize_factor, resize_factor)
my_image <- grayscale(my_image)

plot(my_image)
save.image(my_image, file = "output/grayscale_image.jpg")

# we want to identify the smallest polygon containing the penguin 
# colony, something like the convex hull. Probably simplest to do
# this manually, since that's not the time consuming bit. This
# works, sort of

# blur the image to get a mostly white central bit
blurry <- isoblur(my_image, 10) 
plot(blurry)

# binarize image by turning all "whitish" pixels white, just need
# to fiddle with parameter so that the (probably longest) contour 
# goes around the area you want
binarized_blurry <- blurry > .59
plot(binarized_blurry)
contours <- contours(binarized_blurry)

# find the longest contour 
longest_contour_id <- which.max(unlist(lapply(contours, function(v) length(v$x))))
lines(contours[[longest_contour_id]]$x,
      contours[[longest_contour_id]]$y,col="red")

# make a new image which is white along the contour and black outside,
# this is surprisingly painful, must be an easier way.

# pixels in contour of the area
contour_pixels <- data.frame(x = round(contours[[longest_contour_id]]$x, 0),
                           y = round(contours[[longest_contour_id]]$y, 0),
                           value = 1)

# find xmin, xmax, ymin, ymax of area of interest
xmin <- min(contour_pixels$x)
xmax <- max(contour_pixels$x)
ymin <- min(contour_pixels$y)
ymax <- max(contour_pixels$y)
aoi_dims <- data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
save(aoi_dims, file = "output/aoi_dims.RData")

# start with a black image
all_pixels <- data.frame(expand.grid(x = 1:nrow(binarized_blurry),
                             y = 1:ncol(binarized_blurry)),
                         value = 0)

# all pixels except the contour pixels 
noncontour_pixels <- dplyr::anti_join(all_pixels, contour_pixels, by = c("x","y"))

# now have a data frame with the contour on it
all_pixels <- rbind(contour_pixels, noncontour_pixels) %>% arrange(x,y)

# turn that into an image
area_of_interest <- as.cimg(all_pixels, dims = dim(binarized_blurry))

# voila
plot(area_of_interest)

# choose a pixel outside the contour and flood fill to get a pixel 
# set indicating all pixels inside the area of interest. I assume
# that the pixel (1,1) will be OUTSIDE the area of interest!
pixels_of_interest <- px.flood(area_of_interest, 1, 1, sigma=.1)
plot(pixels_of_interest)

# now go back to the original image and turn all pixels outside the
# area of interest white (= 1)
my_image_orig <- my_image
plot(my_image_orig)
my_image[pixels_of_interest == TRUE] <- 1 # 1 is white
plot(my_image)

# now most of the birds should be much darker than the rest of the image
# set a manual threshold to turn birds black and rest white
binarized_image <- my_image < .28
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
small_patch <- imsub(binarized_image, x %in% 900:1100, y %in% 600:800)
plot(small_patch)

# blur the patch so we can use the blob detector, need to experiment
# with the sigma parameter, depends on size of small_patch and of 
# object you want to detect.
blurry_patch <- isoblur(1 - small_patch, sigma = 5)
plot(blurry_patch)

# get the blob centers, need to experiment with threshold to get it right
centers <- get.centers(blurry_patch, thr = "96%") 

# does poorly at edges but that's ok because we'll shift the patch along
# bit by bit. So remove centers found at the margins, if they're real
# they will be picked up in other patches.
centers <- centers %>% dplyr::filter(mx > 3, mx < 197, my > 3, my < 197)

# plot the estimated centers on top of the image, and check
plot(small_patch)
points(centers$mx,centers$my,col="red",pch=20,cex=1)

