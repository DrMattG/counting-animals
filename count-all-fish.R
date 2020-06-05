library(imager)
library(purrr)

img <- load.image("data/fish.png")
#load("output/aoi_dims.RData")
plot(img)

# function to extract center of a blob
get.centers <- function(im, thr)
{
  dt <- imhessian(im) %$% { xx*yy - xy^2 } %>% threshold(thr) %>% label
  as.data.frame(dt) %>% subset(value>0) %>% 
    dplyr::group_by(value) %>% 
    dplyr::summarise(mx = mean(x), my = mean(y))
}


count_patch <- function(xorig, yorig, patch_size, blur_sigma, thr_centers, im){
  
  # focus in on a small patch of the image, this should smooth over noise
  small_patch <- imsub(1 - im, x %in% xorig:(xorig + patch_size), y %in% yorig:(yorig + patch_size))
  
  # blur the patch so we can use the blob detector, need to experiment
  # with the sigma parameter, depends on size of small_patch and of 
  # object you want to detect.
  blurry_patch <- isoblur(small_patch, sigma = blur_sigma)
  
  # get the blob centers, need to experiment with threshold to get it right
  centers <- get.centers(blurry_patch, thr = thr_centers) 
  
  # does poorly at edges but that's ok because we'll shift the patch along
  # bit by bit. So remove centers found at the 2% margins, if they're real
  # they will be picked up in other patches.
  centers <- centers %>% dplyr::filter(mx > trunc(0.02 * patch_size), mx < trunc(0.98 * patch_size), 
                                       my > trunc(0.02 * patch_size), my < trunc(0.98 * patch_size))
  
  # record origin x and y so that can add back later to recover the true positions
  centers <- centers %>% dplyr::mutate(xorig = xorig, yorig = yorig)
  
  return(centers)
  
}

# makes the origins of all possible patches, a grid covering the area of interest
xt <- seq(from = 0, dim(img)[1], by = 50) 
yt <- seq(from = 0, dim(img)[2], by = 50) 
patch_origins <- expand.grid(x = xt, y = yt)

# get the centers in each patch using purrr
all_centers <- map2_dfr(patch_origins$x, patch_origins$y, .f = count_patch, patch_size = 200, 
                    blur_sigma = 4, thr_centers = "97%", im = img)
xac<-all_centers
all_centers<-xac
# recover the real (x, y) positions of the centers by adding back the patch origins
all_centers <- all_centers %>% mutate(mx = round(xorig + mx, 0), 
                                      my = round(yorig + my, 0)) %>% 
  dplyr::select(-xorig, -yorig) %>% arrange(mx, my)

# add a variable that tells you in how many patches each center was found, can use this
# to identify reliable centers (should have been picked up in many patches)
all_centers <- all_centers %>% dplyr::group_by(mx, my) %>% dplyr::count() %>% dplyr::ungroup()

# remove any centers picked up less than 3 times
all_centers <- all_centers %>% dplyr::filter(n > 1)

# group together any centers that are separated by a distance of less than d_sep units
d_sep <- 10
dist_between_centers <- dist(all_centers)
all_centers$cluster <- as.factor(cutree(hclust(dist_between_centers, method = "complete"), h = d_sep))

# find the cluster centroids and use these as the new positions for the centers
all_centers_grouped <- all_centers %>% dplyr::group_by(cluster) %>% 
  dplyr::summarize(mx = weighted.mean(mx, n), my = weighted.mean(my, n), n = sum(n))

# can do another round of filtering if you want (n here but be > n above to have an effect)
all_centers_grouped <- all_centers_grouped %>% dplyr::filter(n > 3)

# plot the estimated centers on top of the image, and check
xorig_check <- 0
yorig_check <- 0
patch_size_check <- 400
small_patch <- imsub(img, x %in% (xorig_check):(xorig_check + patch_size_check), 
                     y %in% (yorig_check):(yorig_check + patch_size_check))
plot(small_patch)
points(all_centers_grouped$mx - xorig_check, all_centers_grouped$my - yorig_check, 
       col = 'red', pch = 20, cex = 1)

# plot the whole image
plot(img) 
points(all_centers_grouped$mx, all_centers_grouped$my, 
       col = 'red', pch = 20, cex = 0.5)


library(ggplot2)

# saving the result needs converting to ggplot, I think
p1 <- ggplot(my_image %>% as.data.frame(), aes(x, y)) + 
  geom_raster(aes(fill=value)) + 
  geom_point(data = all_centers_grouped, aes(x = mx, y = my), col = "red", size = 0.5) + 
  scale_fill_gradient(low="black",high="white") + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

h_to_w_ratio <- dim(my_image)[2] / dim(my_image)[1]
w <- 12
h <- w * h_to_w_ratio
ggsave("output/Counted_BP3G_AON_2017_dr.jpg", p1, width = w, height = h, dpi = 600)

# lower resolution plot
resize_factor <- 0.5

small_my_image <- resize(my_image, -100*resize_factor, -100*resize_factor)

# saving the result needs converting to ggplot, I think
p1s <- ggplot(small_my_image %>% as.data.frame(), aes(x, y)) + 
  geom_raster(aes(fill=value)) + 
  geom_point(data = all_centers_grouped, aes(x = mx * resize_factor, y = my * resize_factor), col = "red", size = 0.5) + 
  scale_fill_gradient(low="black",high="white") + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

small_my_image_df <- small_my_image %>% as.data.frame()
save(small_my_image_df, all_centers_grouped, file = "output/for_plotly.Rdata")

h_to_w_ratio <- dim(my_image)[2] / dim(my_image)[1]
w <- 12
h <- w * h_to_w_ratio
ggsave("output/smallCounted_BP3G_AON_2017_dr.jpg", p1s, width = w, height = h, dpi = 600)


# not run, just trying out stuff
xorig <- 905
yorig <- 649
centers <- map2_dfr(xorig, yorig, .f = count_patch, patch_size = 200, 
                        blur_sigma = 5, thr_centers = "96%", im = binarized_image)
small_patch <- imsub(binarized_image, x %in% xorig:(xorig + 200), y %in% yorig:(yorig + 200))
plot(small_patch)
points(centers$mx, centers$my, col = 'red', pch = 20, cex = 1)


