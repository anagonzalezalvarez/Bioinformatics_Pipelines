
# Functional Characterization 

## Content

```

```

##  Workflow SST-VIP

### 1. Load LIANA outputs of final robust interactions

```
robust_all_df_p10 <- read.csv(paste0(basepath, "2_liana/results/robust_all_df_p10.csv"))
robust_all_df_p28 <- read.csv(paste0(basepath, "2_liana/results/robust_all_df_p28.csv"))
```


```
> head(robust_all_df_p10)
  source target ligand receptor aggregate_rank  mean_rank call_cellchat.pval
1  Lamp5    Sst  L1cam    Cntn1   0.0001857196   86.66667               0.43
2  Lamp5    Sst   Jam3    Itgb1   0.0037766857  153.83333               0.00
3  Lamp5    Sst  Vegfc    Itgb1   0.0041634512  157.50000               0.00
4  Lamp5    Sst  Nxph1    Nrxn1   0.0058327916 1921.00000                 NA
5  Lamp5  Pvalb  L1cam    Cntn1   0.0018684714  139.83333               0.63
6  Lamp5  Pvalb  Lamc1    Itga9   0.0026755554  141.16667               0.00
  call_cellchat.rank cellphonedb.pvalue cellphonedb.rank sca.LRscore sca.rank
1               17.0              0.000            122.0   0.8537607      121
2                6.5              0.000            122.0   0.6924550      333
3                6.5              0.000            122.0   0.6852910      344
4             5434.0              0.096            327.0   0.9590056        2
5               20.0              0.007            263.5   0.8453654      136
6                9.0              0.000            117.5   0.7186082      297
              interaction
1   Lamp5-Sst-L1cam-Cntn1
2    Lamp5-Sst-Jam3-Itgb1
3   Lamp5-Sst-Vegfc-Itgb1
4   Lamp5-Sst-Nxph1-Nrxn1
5 Lamp5-Pvalb-L1cam-Cntn1
6 Lamp5-Pvalb-Lamc1-Itga9
```



### 2. Create separe objects  for interactions of sub-types of interest, considering the direction of interactions 

```
SSTtoVIP_p10 <- subset(robust_all_df_p10, robust_all_df_p10$source == "Sst" & robust_all_df_p10$target == "Vip")
SSTtoVIP_p28 <- subset(robust_all_df_p28, robust_all_df_p28$source == "Sst" & robust_all_df_p28$target == "Vip")
VIPtoSST_p10 <- subset(robust_all_df_p10, robust_all_df_p10$source == "Vip" & robust_all_df_p10$target == "Sst")
VIPtoSST_p28 <- subset(robust_all_df_p28, robust_all_df_p28$source == "Vip" & robust_all_df_p28$target == "Sst")


```


### 3. Visualize number of unique and shared interactions across the two stages of developemnt: p10 and p28

```
#--- Sst -> Vip 
x <- list(p10 = SSTtoVIP_p10$interaction,
          p28 = SSTtoVIP_p28$interaction)
ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "maroon2") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom") + coord_flip()
  
#----  Vip -> Sst
x <- list(p10 = VIPtoSST_p10$interaction,
          p28 = VIPtoSST_p28$interaction)
ggVennDiagram(x[1:2], label = "count", label_alpha = 0, set_size = 6,label_size = 6) + 
  scale_fill_gradient(low = "#F4FAFE", high = "maroon2") +  theme(legend.title = element_text(color = "black"), legend.position = "bottom") + coord_flip()
```

![picture alt](./content/1_vendiagrams.png.png)
