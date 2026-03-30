# Niche overlap calculation

`nicheoverlap.R` is used for the calculation of niche overlap.

### Environment
* **Platform:** RStudio / R

### Data Preparation
1.  **Generate data:** Use Matlab to produce the `Depth+Phyto+AOA.txt` file.
    ```matlab
    writematrix(trab, 'trab.txt')
    ```
2.  **Data arrangement:**
    
| Format | 1st Column | 2nd to X Column | X+1 to End Column |
| :--- | :--- | :--- | :--- |
| **Standard** | Depth | AOA abundance | Phyto abundance |
| **Time-series (k days)** | Depth | AOA (2 ~ k/2+1) | Phyto (k/2+2 ~ k+1) |

---

### Running Commands
In your R console, run the following:

```r
# 1. Load the function
source('path/nicheoverlap.R') 

# 2. Load the depth and abundance data (csv format)
trab <- read.table('path/trab.txt', sep=",", header=FALSE) 

# 3. Convert to matrix and calculate
trab <- as.matrix(trab)
NO <- nicheoverlap(trab)
