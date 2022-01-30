这是一个RNA-seq检测的实例

使用RNA-seq数据进行nt检测，发现第2和3个样品 大部分reads并没有map到指定物种上，并且对应物种比对率普遍低于使用hisat2比对率。

| Sample        | Highest proportion species | %Percentage | paired in sequencing passed | mapped reads | %mapped reads |
|---------------|----------------------------|-------------|-----------------------------|--------------|---------------|
| control-input | Homo sapiens               | 69.929      | 94417764                    | 91243519     | 96.6381%      |
| control-m6A-1 | Saccharum perrieri         | 20.567      | 110003114                   | 41336042     | 37.5772%      |
| control-m6A-2 | Sapajus apella             | 28.893      | 103461516                   | 94926062     | 91.7501%      |
| treat-input   | Homo sapiens               | 71.988      | 101275712                   | 97750882     | 96.5196%      |
| treat-m6A-1   | Homo sapiens               | 34.262      | 86830578                    | 75957919     | 87.4783%      |
| treat-m6A-2   | Homo sapiens               | 43.191      | 92827518                    | 76440278     | 82.3466%      |

表头说明：
- Highest proportion species : nt比对得分最高物种
- %Percentage : nt比对得分最高物种百分比
- paired in sequencing passed : clean_data reads数量，可以理解为总reads数
- mapped reads : hisat2比对上参考基因组的reads数量
- %mapped reads : hisat2比对上参考基因组的reads数量百分比

