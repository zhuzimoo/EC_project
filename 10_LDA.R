library(topicmodels)
library(tidytext)
library('tidyr')
library(ggplot2)
library(dplyr)
library(pheatmap)

k = 15 # a num indicating the number of topics
output = "./LDA_plot/" # the folder for saving result
commu_file = "./data/lda/liana/liana_LDA_matrix_num_norm_receiver.csv" # row = communication pairs, columns = tissue, samples

k = as.integer(k)
# communication topic matrix
commu_summary = read.csv(commu_file, row.names = 1)
commu_summary = t(commu_summary) # rows are tissue, cols are m_s pairs

ap_lda <- LDA(commu_summary, k = k, control = list(seed = 1234, iter = 2000), method = 'Gibbs')
ap_topics <- tidy(ap_lda, matrix = "beta")
ap_documents <- tidy(ap_lda, matrix = "gamma")

logl <- logLik(ap_lda)
top_terms <- ap_topics %>%
  group_by(topic) %>%
  slice_max(beta, n = 15) %>%
  ungroup() %>%
  arrange(topic, -beta)

# save the result matrices
write.table(ap_topics, file = "./data/lda/liana/liana_num_norm_receiver_LDA_topics_15t.csv", sep = ",", row.names = FALSE)
write.table(ap_documents, file = "./data/lda/liana/liana_num_norm_receiver_LDA_documents_15t.csv", sep = ",", row.names = FALSE)


df = read.csv("./data/lda/liana/liana_num_norm_receiver_LDA_documents_15t.csv")
new <- df %>%
  pivot_wider(
    id_cols = c('document'),
    names_from = 'topic',
    values_from = 'gamma'
  )

df_new <- as.matrix(new[, 2:15])
rownames(df_new) <- new$document

# Fig 4A and 5A
hm <- pheatmap::pheatmap(df_new, color=colorRampPalette(c("white", "light blue", "dodgerblue", "mediumorchid4"))(100), scale="none",
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         fontsize_row=10, fontsize_col=8)
setHook("grid.newpage", NULL, "replace")
grid.text("topic", y=-0.03, x=0.45,  gp=gpar(fontsize=10))
grid.text("tissue", x=-0.05, rot=90, gp=gpar(fontsize=10))
