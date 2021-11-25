
TurkeyPlot <- function(seurat_object, 
                       group.by,
                       feather_length = 6,
                       feather_width = 1.5,
                       body_width = 5,
                       pt.size = 0.5) {
    
    #parameter setting
    origin <- c(0,0)
    plot_center_x = origin[1]
    plot_center_y = origin[2]
    feather_x_radius = feather_width
    feather_y_radius = feather_length
    body_radius <- body_width
    extend_radius <- body_radius+feather_y_radius
    
    #getting cluster info
    Idents(seurat_object) <- group.by
    num_clusters <- length(unique(Idents(seurat_object)))
    cluster_counts <- table(Idents(seurat_object))
    
    #making body and head
    head_size <- floor(cluster_counts[1]/4)
    body_size <- cluster_counts[1] - head_size
    
    body_theta <- 2*pi*runif(body_size, 0, 1)
    body_d <- body_radius * sqrt(runif(body_size, 0, 1))
    x = (plot_center_x + body_d * cos(body_theta)) 
    y = (plot_center_y + body_d * sin(body_theta)) 
    cluster_ident <- names(cluster_counts)[1]
    body_coords <- data.frame(cbind(x, y), cluster = cluster_ident)
    colnames(body_coords)[1:2] <- c("Turkey_1", "Turkey_2")
    
    head_theta = 2*pi * runif(head_size,0, 1)
    head_d = body_radius/2 * sqrt(runif(head_size, 0, 1))
    head_center_x <- -1.1*body_radius
    head_center_y <- 1.1*body_radius
    x = (head_center_x + head_d * cos(head_theta)) 
    y = (head_center_y + head_d * sin(head_theta)) 
    cluster_ident <- names(cluster_counts)[1]
    head_coords <- data.frame(cbind(x, y), cluster = cluster_ident)
    colnames(head_coords)[1:2] <- c("Turkey_1", "Turkey_2")
    
    all_coords <- rbind(body_coords, head_coords)
    
    #setting feathers
    feather_clusters <- cluster_counts[-1]
    num_feathers <- length(feather_clusters)
    thetas <- seq(pi/2, pi/8, length.out = num_feathers)
    feather_coords <- lapply(1:length(feather_clusters), function(i){
        N <- feather_clusters[i]
        theta <- thetas[i]
        t = 2*pi * runif(N,0, 1)
        d = sqrt(runif(N, 0, 1))
        x = (plot_center_x  + feather_x_radius * d * cos(t)) 
        y = (plot_center_y  + feather_y_radius * d * sin(t)) 
        coord_mat <- cbind(x,y)
        rot_theta <- theta-(pi/2)
        rot_matrix <- matrix(c(cos(rot_theta), -sin(rot_theta), sin(rot_theta), cos(rot_theta)), nrow = 2, ncol = 2)
        new_mat <- coord_mat %*% rot_matrix
        new_mat[,1] <- new_mat[,1] + (extend_radius * cos(theta))
        new_mat[,2] <- new_mat[,2] + (extend_radius * sin(theta))
        cluster_ident <- names(feather_clusters)[i]
        new_mat <- data.frame(new_mat, cluster = cluster_ident)
        colnames(new_mat)[1:2] <- c("Turkey_1", "Turkey_2")
        return(new_mat)
    }) %>% do.call(rbind, .)
    
    #making the beak
    beak_coord_1x <- min(head_coords[,1])
    beak_coord_1y <- head_center_y
    beak_coord_2x <- (beak_coord_1x + head_center_x)/2
    beak_coord_2y <- (beak_coord_1y + min(head_coords[,2]))/2.5
    beak_coord_3x <- beak_coord_1x - abs(beak_coord_1x - beak_coord_2x)
    beak_coord_3y <- beak_coord_2y
    beak_df <- data.frame(Turkey_1 = c(beak_coord_1x, beak_coord_2x, beak_coord_3x),
                          Turkey_2 = c(beak_coord_1y, beak_coord_2y, beak_coord_3y))
    
    #making the legs
    r_leg_x1 <- body_radius * cos(-3*pi/8)
    r_leg_x2 <- body_radius * cos(-3*pi/8) * 1.5
    r_leg_y1 <- body_radius * sin(-3*pi/8)
    r_leg_y2 <- body_radius * sin(-3*pi/8) * 1.5
    r_leg_df <- data.frame(Turkey_1 = c(r_leg_x1, r_leg_x2),
                           Turkey_2 = c(r_leg_y1, r_leg_y2))
    
    l_leg_x1 <- body_radius * cos(-5*pi/8)
    l_leg_x2 <- body_radius * cos(-5*pi/8) * 1.5
    l_leg_y1 <- body_radius * sin(-5*pi/8)
    l_leg_y2 <- body_radius * sin(-5*pi/8) * 1.5
    l_leg_df <- data.frame(Turkey_1 = c(l_leg_x1, l_leg_x2),
                           Turkey_2 = c(l_leg_y1, l_leg_y2))
    
    

    all_coords <- rbind(all_coords, feather_coords)
    all_coords$cluster <- factor(all_coords$cluster, levels = names(cluster_counts))
    color_pal <- c("chocolate4", pals::brewer.set1(num_clusters-1))
    ggplot(all_coords, aes(Turkey_1, Turkey_2, color = cluster))+
        geom_point(size = pt.size)+
        scale_color_manual(values = color_pal)+
        geom_polygon(data = beak_df, color = "black", fill = "gold", alpha = 0.5)+
        geom_path(data = r_leg_df, color = "black", size = 2)+
        geom_path(data = l_leg_df, color = "black", size = 2)+
        theme_classic()
    
}











