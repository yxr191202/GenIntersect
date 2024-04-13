Calculate Intersections of Gene Sets
Description
Calculates intersections of gene sets from multiple datasets. Outputs a data frame listing each combination of datasets.

Usage
calculate_combination_intersects(list_of_datasets, names_of_datasets)
Arguments
list_of_datasets	
List of data frames with gene names.

names_of_datasets	
Names corresponding to datasets.

Value
Data frame of dataset combinations and intersection genes.

Examples

```R
df1 <- data.frame(Gene = c("Gene1", "Gene2", "Gene3", "Gene4"))
df2 <- data.frame(Gene = c("Gene2", "Gene3", "Gene5"))
df3 <- data.frame(Gene = c("Gene3", "Gene6"))
df4 <- data.frame(Gene = c("Gene1", "Gene3", "Gene7"))
list_of_datasets <- list(df1, df2, df3, df4)
names_of_datasets <- c("Dataset1", "Dataset2","Dataset3", "Dataset4")
results <- calculate_combination_intersects(list_of_datasets, names_of_datasets)
print(results)
```

