data = read.csv('/Users/haticeyazici/Downloads/DATASHARE29_Fullnd.csv')
family <- data$Family

features = data[, c('PDBx',
              'Ts.Tm', 'C26.Tet', 'C28.C29', 'PAH.RI', 'SC2D.SC2P',
              'SC3D.SC3P', 'C28.C29TT', 'C20.C23TT', 'C22.C21TT', 'C24.C23TT',
              'C26.C25TT', 'C31S.H', 'C29H.H', 'C35.C34S.hopanes',
              'BNH.H', 'OI.H', 'G.H', 'C29Ts.C29H')]
corr_features = data[, c('PDBx',
                      'Ts.Tm', 'C26.Tet', 'C28.C29', 'PAH.RI',
                       'SC3D.SC3P', 'C28.C29TT', 'C20.C23TT', 'C24.C23TT',
                       'C31S.H', 'C29H.H', 
                       'BNH.H', 'OI.H', 'C29Ts.C29H')]
corr_features2 = data[, c('PDBx',
                         'Ts.Tm', 'C26.Tet', 'C28.C29', 'PAH.RI',
                         'SC3D.SC3P', 'C28.C29TT', 'C20.C23TT', 'C24.C23TT',
                         'C31S.H', 'C29H.H', 'C35.C34S.hopanes',
                         'BNH.H', 'OI.H', 'C29Ts.C29H')]
library(kohonen)

set.seed(61)

oil.grid <- somgrid(xdim = 12, ydim = 11, topo = "hexagonal")

# make a SOM model

ads.model <- som(as.matrix(features), oil.grid, rlen = 100, radius = 1, 
                 keep.data = TRUE,
                 dist.fcts = "euclidean")

#shows the samples on their respective SOM nodes
plot(ads.model, type = "mapping", pchs = data$Sample, col=family,
     shape = "straight", 
     palette.name = rainbow)


#quantization error: mean distance between each observation and its
#Best Matching Unit (BMU).

distances1 = ads.model$distances  
quant.error = mean(distances1)

#variance explained
data_matrix1 = ads.model$data[[1]]  # Extracts input data matrix from model

total_variance1 = sum(apply(data_matrix1, 2, var))  # overall variance
variance_explained1 = 1 - (quant.error / total_variance1)
variance_explained1
# --- 1. Define Features and Color Vector ---
nfeatures <- ncol(ads.model$codes[[1]])
feature_names <- colnames(ads.model$codes[[1]])
feature_colors <- rainbow(nfeatures) # Use rainbow for a vibrant, standard palette

# --- 2. Set up a larger plotting area ---
par(mfrow = c(1, 1), mar = c(1, 1, 1, 15)) 

# --- 3. Generate the Star Plot ---
# Key Fix: Explicitly use the 'rainbow' palette.
plot(ads.model, 
     type = "codes", 
     main = "Codebook Vectors (19-Feature Star Plots)",
     codeRendering = "stars", 
     shape = "straight",
     palette.name = rainbow # <--- Using 'rainbow'
)

# --- 4. Draw the Legend ---
legend("right", 
       legend = feature_names,
       col = feature_colors,
       pch = 15,
       cex = 0.75,
       bty = "n",
       inset = c(0, 0)
)

# --- 5. Reset plot layout ---
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

#generate a new model with reduced features
ads.model2 = som(as.matrix(corr_features), oil.grid, rlen = 100, 
                              radius = 1, 
                              keep.data = TRUE,
                              dist.fcts = "euclidean")


#shows the samples on their respective SOM nodes
plot(ads.model2, type = "mapping", labels = data$Sample, col=family,
     shape = "straight", 
     palette.name = rainbow)

plot(ads.model2, type = "codes",
     shape = "straight", 
     palette.name = rainbow)

#U-matrix
plot(ads.model2, type="dist.neighbours", main = "SOM neighbour distances",
     shape='straight')

plot(ads.model2, type="property",
     shape='straight', property = getCodes(ads.model2, 1)[,1])

plot(ads.model2, type="quality", shape='straight')

#quantization error: mean distance between each observation and its
#Best Matching Unit (BMU).

distances = ads.model2$distances  
quant.error2 = mean(distances)

#variance explained
data_matrix = ads.model2$data[[1]]  # Extracts input data matrix from model

total_variance = sum(apply(data_matrix, 2, var))  # overall variance
variance_explained = 1 - (quant.error2 / total_variance)
variance_explained
#attempt to extract the elements from the nodes on the mapping plot:

#som.prediction = predict(ads.model2)
plot(ads.model2, type = "mapping", labels = data$Sample, col=family,
     shape = "straight", 
     palette.name = rainbow,
     classif =som.prediction$unit.classif
     )

# this vector should supposedly tell us which node index each row was mapped to.
node.assignment.vector = ads.model2$unit.classif
sample.names = row.names(ads.model$data[[1]])

mapping_df <- data.frame(
  ObservationName = data$Sample,
  NodeIndex = node.assignment.vector
)

# Group by NodeIndex and list all the ObservationNames in a new column.
names_per_node_list <- aggregate(
  ObservationName ~ NodeIndex,
  data = mapping_df,
  FUN = function(x) paste(x, collapse = ", ")
)

# Convert the resulting data frame into a list for easier access
names_per_node <- split(
  mapping_df$ObservationName,
  mapping_df$NodeIndex
)

write.csv(names_per_node_list, 
          file = "/Users/haticeyazici/Downloads/Mapping_info.csv", 
          row.names = FALSE)

node.features = as.data.frame(ads.model$codes[[1]])
#node.features$node = 1:nrow(node.features)

write.csv(node.features, 
          file = "/Users/haticeyazici/Downloads/Codes_info.csv", 
          row.names = TRUE)

