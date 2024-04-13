## Description

计算多个数据集中基因集的交集。输出列出每个数据集组合的数据框。

## Usage

calculate_combination_intersects(list_of_datasets, names_of_datasets)

## Arguments

| 列表                                | 描述                             |
| ----------------------------------- | -------------------------------- |
| list_of_datasets                    | names_of_datasets                |
| List of data frames with gene names | Names corresponding to datasets. |

## Value

Data frame of dataset combinations and intersection genes.

## Examples

```R
> df1 <- data.frame(Gene = c("Gene1", "Gene2", "Gene3", "Gene4"))
> df2 <- data.frame(Gene = c("Gene2", "Gene3", "Gene5"))
> df3 <- data.frame(Gene = c("Gene3", "Gene6"))
> df4 <- data.frame(Gene = c("Gene1", "Gene3", "Gene7"))
> list_of_datasets <- list(df1, df2, df3, df4)
> names_of_datasets <- c("Dataset1", "Dataset2","Dataset3", "Dataset4")
> results <- calculate_combination_intersects(list_of_datasets, names_of_datasets)
> print(results)
                                                                Dataset.Combination Intersection.Genes Number.of.Genes
Dataset1 & Dataset2                                             Dataset1 & Dataset2       Gene2, Gene3               2
Dataset1 & Dataset3                                             Dataset1 & Dataset3              Gene3               1
Dataset1 & Dataset4                                             Dataset1 & Dataset4       Gene1, Gene3               2
Dataset2 & Dataset3                                             Dataset2 & Dataset3              Gene3               1
Dataset2 & Dataset4                                             Dataset2 & Dataset4              Gene3               1
Dataset3 & Dataset4                                             Dataset3 & Dataset4              Gene3               1
Dataset1 & Dataset2 & Dataset3                       Dataset1 & Dataset2 & Dataset3              Gene3               1
Dataset1 & Dataset2 & Dataset4                       Dataset1 & Dataset2 & Dataset4              Gene3               1
Dataset1 & Dataset3 & Dataset4                       Dataset1 & Dataset3 & Dataset4              Gene3               1
Dataset2 & Dataset3 & Dataset4                       Dataset2 & Dataset3 & Dataset4              Gene3               1
Dataset1 & Dataset2 & Dataset3 & Dataset4 Dataset1 & Dataset2 & Dataset3 & Dataset4              Gene3  
1
```




