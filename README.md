# OrthoBridge
基于序列以及blast的同源基因转换


## 基于蛋白质序列进行同源基因转换

对于亲缘关系非常远的物种，基于DNA序列进行blast不再合适，需要根据蛋白质序列进行blast，以对非模式动物的基因进行功能注释与解析。

### 1.从NCBI下载nr.gz数据库

首先使用`ascp`高速下载nr.gz数据库，优点是：速度较快；允许断点续传。

需要注意的是`-i`参数提供**链接所需的密令文件**，这个文件需要根据自己**aspera**安装的位置更换路径，其他参数不需要改变。

目前NCBI不再提供单独物种的蛋白质序列，需要将所有的序列下载后单独进行构建。

此外，下载过程中**NCBI**似乎对访问流量进行了限制，实际下载过程中，大约每下载**40G**就会中断一次，需要重新下载，所以断点续传功能非常重要。

```
ascp -v -k 1 -QT -l 300M -i /home/bio/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/FASTA/nr.gz .
```

这是NCBI对应的[nr.gz](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)所在的ftp链接位置。

### 2.从nr.gz中提取特定物种信息构建相应的比对数据库

目前

### 3.基于基因列表批量从NCBI上下载对应蛋白质序列

