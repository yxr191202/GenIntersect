#' Calculate Intersections of Gene Sets
#'
#' Calculates intersections of gene sets from multiple datasets.
#' Outputs a data frame listing each combination of datasets.
#'
#' @param list_of_datasets List of data frames with gene names.
#' @param names_of_datasets Names corresponding to datasets.
#' @return Data frame of dataset combinations and intersection genes.
#' @export
#' @examples
#' df1 <- data.frame(Gene = c("Gene1", "Gene2", "Gene3", "Gene4"))
#' df2 <- data.frame(Gene = c("Gene2", "Gene3", "Gene5"))
#' df3 <- data.frame(Gene = c("Gene3", "Gene6"))
#' df4 <- data.frame(Gene = c("Gene1", "Gene3", "Gene7"))
#' list_of_datasets <- list(df1, df2, df3, df4)
#' names_of_datasets <- c("Dataset1", "Dataset2","Dataset3", "Dataset4")
#' results <- calculate_combination_intersects(list_of_datasets, names_of_datasets)
#' print(results)

calculate_combination_intersects <- function(list_of_datasets, names_of_datasets) {
  # 验证输入是否正确
  if (length(list_of_datasets) == 0) {
    stop("No datasets provided")
  }

  # 验证数据集名称是否正确
  if (length(names_of_datasets) != length(list_of_datasets)) {
    stop("Names of datasets do not match the number of datasets provided")
  }

  # 存储结果的列表
  results_list <- list()
  combinations_list <- list()

  # 对2个、3个和4个数据集的所有组合进行操作
  for (i in 2:length(list_of_datasets)) {
    # 获取当前组合数的所有组合
    combn_indices <- combn(seq_along(list_of_datasets), i, simplify = FALSE)

    # 对每种组合计算交集
    for (indices in combn_indices) {
      # 提取组合中的基因名列表
      gene_lists <- lapply(list_of_datasets[indices], function(df) df$Gene)
      # 计算交集
      intersect_genes <- Reduce(intersect, gene_lists)

      # 组合名称
      combination_name <- paste(names_of_datasets[indices], collapse = " & ")

      # 生成交集结果
      if (length(intersect_genes) > 0) {
        results_list[[combination_name]] <- intersect_genes
      } else {
        results_list[[combination_name]] <- character(0)
      }

      # 记录组合名
      combinations_list[[combination_name]] <- combination_name
    }
  }

  # 创建数据框以输出结果
  results_df <- data.frame(
    "Dataset Combination" = names(combinations_list),
    "Intersection Genes" = sapply(results_list, function(genes) paste(genes, collapse=", ")),
    "Number of Genes" = sapply(results_list, length),
    stringsAsFactors = FALSE
  )

  return(results_df)
}

