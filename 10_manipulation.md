# (PART) Focus Topics {-}


# Data Manipulation {#datamanipulation}

<script>
document.addEventListener("click", function (event) {
    if (event.target.classList.contains("rebook-collapse")) {
        event.target.classList.toggle("active");
        var content = event.target.nextElementSibling;
        if (content.style.display === "block") {
            content.style.display = "none";
        } else {
            content.style.display = "block";
        }
    }
})
</script>

<style>
.rebook-collapse {
  background-color: #eee;
  color: #444;
  cursor: pointer;
  padding: 18px;
  width: 100%;
  border: none;
  text-align: left;
  outline: none;
  font-size: 15px;
}

.rebook-content {
  padding: 0 18px;
  display: none;
  overflow: hidden;
  background-color: #f1f1f1;
}
</style>


## Tidying and subsetting

### Tidy data

For several custom analysis and visualization packages, such as those from the 
`tidyverse`, the `SE` data can be converted to long data.frame format with 
`meltAssay`.    



```r
library(mia)
data(GlobalPatterns, package="mia")
tse <- GlobalPatterns
tse <- transformSamples(tse, method="relabundance")

molten_tse <- meltAssay(tse,
                        add_row_data = TRUE,
                        add_col_data = TRUE,
                        abund_values = "relabundance")
molten_tse
```

```
## # A tibble: 499,616 × 17
##    FeatureID SampleID relabundance Kingdom Phylum    Class    Order Family Genus
##    <fct>     <fct>           <dbl> <chr>   <chr>     <chr>    <chr> <chr>  <chr>
##  1 549322    CL3                 0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  2 549322    CC1                 0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  3 549322    SV1                 0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  4 549322    M31Fcsw             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  5 549322    M11Fcsw             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  6 549322    M31Plmr             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  7 549322    M11Plmr             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  8 549322    F21Plmr             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
##  9 549322    M31Tong             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
## 10 549322    M11Tong             0 Archaea Crenarch… Thermop… <NA>  <NA>   <NA> 
## # … with 499,606 more rows, and 8 more variables: Species <chr>,
## #   X.SampleID <fct>, Primer <fct>, Final_Barcode <fct>,
## #   Barcode_truncated_plus_T <fct>, Barcode_full_length <fct>,
## #   SampleType <fct>, Description <fct>
```

### Subsetting

**Subsetting** data helps to draw the focus of analysis on particular
  sets of samples and / or features. When dealing with large data
  sets, the subset of interest can be extracted and investigated
  separately. This might improve performance and reduce the
  computational load.

Load:

* mia
* dplyr
* knitr
* data `GlobalPatterns`



Let us store `GlobalPatterns` into `tse` and check its original number of features (rows) and samples (columns). **Note**: when subsetting by sample, expect the number of columns to decrease; when subsetting by feature, expect the number of rows to decrease.


```r
# store data into se and check dimensions
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns
# show dimensions (features x samples)
dim(tse) 
```

```
## [1] 19216    26
```

#### Subset by sample (column-wise)

For the sake of demonstration, here we will extract a subset containing only the samples of human origin (feces, skin or tongue), stored as `SampleType` within `colData(tse)` and also in `tse`.

First, we would like to see all the possible values that `SampleType` can take on and how frequent those are: 


```r
# inspect possible values for SampleType
unique(tse$SampleType)
```

```
## [1] Soil               Feces              Skin               Tongue            
## [5] Freshwater         Freshwater (creek) Ocean              Sediment (estuary)
## [9] Mock              
## 9 Levels: Feces Freshwater Freshwater (creek) Mock ... Tongue
```

```r
# show recurrence for each value
tse$SampleType %>% table()
```
<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> . </th>
   <th style="text-align:right;"> Freq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Feces </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Freshwater </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Freshwater (creek) </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mock </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ocean </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sediment (estuary) </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Skin </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Soil </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tongue </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
</tbody>
</table></div>

**Note**: after subsetting, expect the number of columns to equal the
  sum of the recurrences of the samples that you are interested
  in. For instance, `ncols = Feces + Skin + Tongue = 4 + 3 + 2 = 9`.

Next, we _logical index_ across the columns of `tse` (make sure to
leave the first index empty to select all rows) and filter for the
samples of human origin. For this, we use the information on the
samples from the meta data `colData(tse)`.


```r
# subset by sample
tse_subset_by_sample <- tse[ , tse$SampleType %in% c("Feces", "Skin", "Tongue")]

# show dimensions
dim(tse_subset_by_sample)
```

```
## [1] 19216     9
```

As a sanity check, the new object `tse_subset_by_sample` should have
the original number of features (rows) and a number of samples
(columns) equal to the sum of the samples of interest (in this case
9).

Several characteristics can be used to subset by sample:

* origin
* sampling time
* sequencing method
* DNA / RNA barcode
* cohort

#### Subset by feature (row-wise)

Similarly, here we will extract a subset containing only the features
that belong to the Phyla "Actinobacteria" and "Chlamydiae", stored as
`Phylum` within `rowData(tse)`. However, subsetting by feature implies
a few more obstacles, such as the presence of NA elements and the
possible need for agglomeration.

As previously, we would first like to see all the possible values that
`Phylum` can take on and how frequent those are:
  

```r
# inspect possible values for Phylum
unique(rowData(tse)$Phylum)
```

```
##  [1] "Crenarchaeota"    "Euryarchaeota"    "Actinobacteria"   "Spirochaetes"    
##  [5] "MVP-15"           "Proteobacteria"   "SBR1093"          "Fusobacteria"    
##  [9] "Tenericutes"      "ZB3"              "Cyanobacteria"    "GOUTA4"          
## [13] "TG3"              "Chlorobi"         "Bacteroidetes"    "Caldithrix"      
## [17] "KSB1"             "SAR406"           "LCP-89"           "Thermi"          
## [21] "Gemmatimonadetes" "Fibrobacteres"    "GN06"             "AC1"             
## [25] "TM6"              "OP8"              "Elusimicrobia"    "NC10"            
## [29] "SPAM"             NA                 "Acidobacteria"    "CCM11b"          
## [33] "Nitrospirae"      "NKB19"            "BRC1"             "Hyd24-12"        
## [37] "WS3"              "PAUC34f"          "GN04"             "GN12"            
## [41] "Verrucomicrobia"  "Lentisphaerae"    "LD1"              "Chlamydiae"      
## [45] "OP3"              "Planctomycetes"   "Firmicutes"       "OP9"             
## [49] "WPS-2"            "Armatimonadetes"  "SC3"              "TM7"             
## [53] "GN02"             "SM2F11"           "ABY1_OD1"         "ZB2"             
## [57] "OP11"             "Chloroflexi"      "SC4"              "WS1"             
## [61] "GAL15"            "AD3"              "WS2"              "Caldiserica"     
## [65] "Thermotogae"      "Synergistetes"    "SR1"
```

```r
# show recurrence for each value
rowData(tse)$Phylum %>% table()
```
<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> . </th>
   <th style="text-align:right;"> Freq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ABY1_OD1 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AC1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acidobacteria </td>
   <td style="text-align:right;"> 1021 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Actinobacteria </td>
   <td style="text-align:right;"> 1631 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AD3 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Armatimonadetes </td>
   <td style="text-align:right;"> 61 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bacteroidetes </td>
   <td style="text-align:right;"> 2382 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BRC1 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caldiserica </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Caldithrix </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCM11b </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chlamydiae </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chlorobi </td>
   <td style="text-align:right;"> 64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Chloroflexi </td>
   <td style="text-align:right;"> 437 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Crenarchaeota </td>
   <td style="text-align:right;"> 106 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cyanobacteria </td>
   <td style="text-align:right;"> 393 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Elusimicrobia </td>
   <td style="text-align:right;"> 31 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Euryarchaeota </td>
   <td style="text-align:right;"> 102 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fibrobacteres </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Firmicutes </td>
   <td style="text-align:right;"> 4356 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fusobacteria </td>
   <td style="text-align:right;"> 37 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GAL15 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gemmatimonadetes </td>
   <td style="text-align:right;"> 191 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GN02 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GN04 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GN06 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GN12 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GOUTA4 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hyd24-12 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KSB1 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCP-89 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LD1 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lentisphaerae </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MVP-15 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NC10 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nitrospirae </td>
   <td style="text-align:right;"> 74 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NKB19 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OP11 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OP3 </td>
   <td style="text-align:right;"> 30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OP8 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OP9 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PAUC34f </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Planctomycetes </td>
   <td style="text-align:right;"> 638 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Proteobacteria </td>
   <td style="text-align:right;"> 6416 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAR406 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SBR1093 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SC3 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SC4 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SM2F11 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SPAM </td>
   <td style="text-align:right;"> 22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Spirochaetes </td>
   <td style="text-align:right;"> 124 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SR1 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Synergistetes </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tenericutes </td>
   <td style="text-align:right;"> 143 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TG3 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermi </td>
   <td style="text-align:right;"> 46 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thermotogae </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TM6 </td>
   <td style="text-align:right;"> 27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TM7 </td>
   <td style="text-align:right;"> 32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Verrucomicrobia </td>
   <td style="text-align:right;"> 470 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WPS-2 </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WS1 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WS2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WS3 </td>
   <td style="text-align:right;"> 70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZB2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZB3 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
</tbody>
</table></div>

**Note**: after subsetting, expect the number of columns to equal the
  sum of the recurrences of the feature(s) that you are interested
  in. For instance, `nrows = Actinobacteria + Chlamydiae = 1631 + 21 =
  1652`.

Depending on your research question, you might need to or need not
agglomerate the data in the first place: if you want to find the
abundance of each and every feature that belongs to Actinobacteria and
Chlamydiae, agglomeration is not needed; if you want to find the total
abundance of all the features that belong to Actinobacteria or
Chlamydiae, agglomeration is recommended.

##### Non-agglomerated data

Next, we _logical index_ across the rows of `tse` (make sure to leave
the second index empty to select all columns) and filter for the
features that fall in either Actinobacteria or Chlamydiae. For this,
we use the information on the samples from the meta data
`rowData(tse)`.

The first term with the `%in%` operator are includes all the features
of interest, whereas the second term after the AND operator `&`
filters out all the features that present a NA in place of Phylum.


```r
# subset by feature
tse_subset_by_feature <- tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse)$Phylum), ]

# show dimensions
dim(tse_subset_by_feature)
```

```
## [1] 1652   26
```

As a sanity check, the new object `tse_subset_by_feature` should have the original number of samples (columns) and a number of features (rows) equal to the sum of the features of interest (in this case 1652).

##### Agglomerated data

When total abundances of certain Phyla are of relevance, the data is initially agglomerated by Phylum. Then, similar steps as in the case of not agglomerated data are followed.


```r
# agglomerate by Phylum
tse_phylum <- tse %>% agglomerateByRank(rank = "Phylum")

# subset by feature and get rid of NAs
tse_phylum_subset_by_feature <- tse_phylum[rowData(tse_phylum)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse_phylum)$Phylum), ]

# show dimensions
dim(tse_phylum_subset_by_feature)
```

```
## [1]  2 26
```

**Note**: as data was agglomerated, the number of rows equal the
  number of Phyla used to index (in this case, just 2)

Alternatively:


```r
# store features of interest into phyla
phyla <- c("Phylum:Actinobacteria", "Phylum:Chlamydiae")
# subset by feature
tse_phylum_subset_by_feature <- tse_phylum[phyla, ]
# show dimensions
dim(tse_subset_by_feature)
```

```
## [1] 1652   26
```

The code above returns the not agglomerated version of the data.

Fewer characteristics can be used to subset by feature:

* Taxonomic rank
* Meta-taxonomic group

For subsetting by Kingdom, agglomeration does not apply, whereas for
the other ranks it can be applied if necessary.

#### Subset by sample and feature

Finally, we can subset data by sample and feature at once. The
resulting subset contains all the samples of human origin and all the
features of Phyla "Actinobacteria" or "Chlamydiae".


```r
# subset by sample and feature and get rid of NAs
tse_subset_by_sample_feature <- tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse)$Phylum), tse$SampleType %in% c("Feces", "Skin", "Tongue")]

# show dimensions
dim(tse_subset_by_sample_feature)
```

```
## [1] 1652    9
```

**Note**: the dimensions of `tse_subset_by_sample_feature` agree with
  those of the previous subsets (9 columns filtered by sample and 1652
  rows filtered by feature).

If a study was to consider and quantify the presence of Actinobacteria
as well as Chlamydiae in different sites of the human body,
`tse_subset_by_sample_feature` might be a suitable subset to start
with.

