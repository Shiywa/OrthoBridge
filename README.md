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

我们在*Ixodes scapularis*物种的基因进行分析时，发现无法直接通过其基因ID在**uniprot**上下载对应的序列，只有从**NCBI**上方便查阅。

NCBI提供了`esearch`等命令，方便用户在终端通过网络访问**NCBI**进行查询。

```
esearch -db gene -query "8033311 [uid] AND Ixodes scapularis [ORGN]" |   elink -target protein |   efetch -format fasta
```

这个命令可以直接通过**Entrze ID**例如**8033311**获得其对应的蛋白质序列。

如果知道的是**Symbol**而不是**Entrze ID**，只需要将**uid**换成**GENE**：

```
esearch -db gene -query "FASN1 [GENE] AND Ixodes scapularis [ORGN]" |   elink -target protein |   efetch -format fasta
```

这个命令非常方便，但是有几个问题，首先是：1.一个基因可能对应一个蛋白质，也可能对应多个蛋白质，还可能因为是非编码RNA不对应蛋白质，尽管多个蛋白质能够顺利输出，但是我们进行批量下载的时候，很难再回溯回去，将我们的**Entrze ID**
和蛋白质的**acc id**对应起来，所以需要我们进行代码的优化，添加id转换的文件；2.此外，我们在批量下载过程中发现shell脚本中的`while do; done`无法顺利的循环输出，常规的循环也会导致对**NCBI**节点瞬时大批量访问，导致请求被拒绝，最终
才用了一个很笨的办法，即对每个单独的**Entrze ID**生成一个查询并下载的脚本，在依次提交，这样会导致查询速度很慢，但是能够稳定的下载序列，查询速度大概是：16000个基因花了两天半的时间。

下面是对应的总脚本，主要功能是：1.根据**Entrze ID**列表文件生成每个基因的单独查询文件；2.将**Entrze ID**与**Acc id**对应关系输出到**gene_protein_map.txt**文件；3.最终通过**串行执行**还是**并行执行**提交生成的脚本，目前仅支持
**串行执行**，防止瞬时访问次数过多。

```
#!/bin/bash

# 检查是否传入基因列表文件
if [ "$#" -ne 1 ]; then
    echo "❌ 请提供基因列表文件作为参数"
    echo "用法: $0 gene_list.txt"
    exit 1
fi

# 输入基因列表文件
gene_list_file=$1

# 创建存放脚本的目录
mkdir -p gene_scripts
rm -f gene_scripts/run_*.sh gene_protein_map.txt  # 清空旧脚本和映射文件

# 读取基因列表文件并遍历
while IFS= read -r gene || [[ -n "$gene" ]]; do
    # 确保 gene 变量正确赋值
    if [[ -z "$gene" ]]; then
        continue
    fi

    script_file="gene_scripts/run_${gene}.sh"

    # 生成独立的运行脚本
    cat <<EOF > "$script_file"
#!/bin/bash
echo '🔍 处理基因: $gene'
protein_ids=\$(esearch -db gene -query "$gene [UID] AND Ixodes scapularis [ORGN]" | \\
  elink -target protein | \\
  efetch -format uid | grep -E '^[0-9]+\$')

if [[ -n "\$protein_ids" ]]; then
    refseq_ids=\$(for id in \$protein_ids; do efetch -db protein -id "\$id" -format acc; done | tr '\n' ',')  # 将所有acc ID合并成一行
    echo "  ✅ 找到蛋白 RefSeq IDs: \$refseq_ids"
    echo "$gene \$refseq_ids" >> gene_protein_map.txt  # 记录基因 ID 和蛋白 RefSeq ID
    for refseq_id in \$refseq_ids; do
        efetch -db protein -id "\$refseq_id" -format fasta >> Ixodes_proteins.fasta
    done
else
    echo "  ⚠️ 未找到蛋白质 RefSeq ID，跳过 $gene"
    echo "$gene NONE" >> gene_protein_map.txt  # 记录未找到的情况
fi
EOF

    chmod +x "$script_file"  # 赋予可执行权限
    echo "✅ 生成任务: $script_file"  # 让用户看到生成的脚本
done < "$gene_list_file"

echo "✅ 所有任务脚本已生成，存放于 gene_scripts/ 目录下。"

# 选择执行方式
echo "如何执行这些任务？"
echo "1) 串行执行（适用于小量任务）"
echo "2) 并行执行（适用于大量任务，需要 parallel）"
read -p "选择执行模式 (1/2): " mode

if [[ "$mode" == "1" ]]; then
    echo "🚀 串行执行所有任务..."
    for script in gene_scripts/run_*.sh; do
        bash "$script"
    done
elif [[ "$mode" == "2" ]]; then
    if command -v parallel &>/dev/null; then
        echo "🚀 并行执行所有任务..."
        parallel bash ::: gene_scripts/run_*.sh
    else
        echo "⚠️ 未安装 parallel，请使用 'sudo apt install parallel' 或 'brew install parallel' 安装。"
    fi
else
    echo "❌ 无效选择，脚本已退出。"
fi
```
执行方式为`sh generate_and_run.sh your_gene_list.txt`，`your_gene_list.txt`文件类似于：

```
8030992
8033311
121835630
8041617
115309236
```

执行过程中，命令行会打印执行的进度，类似于：

```
🔍 处理基因: 8054078
  ✅ 找到蛋白 RefSeq IDs: XP_029848343.1,
🔍 处理基因: 8054079
  ✅ 找到蛋白 RefSeq IDs: XP_029848341.1,XP_029848342.1,
🔍 处理基因: 8054081
  ✅ 找到蛋白 RefSeq IDs: XP_029848345.1,XP_029848346.1,
🔍 处理基因: 8054090
  ✅ 找到蛋白 RefSeq IDs: XP_040356995.1,XP_040356997.1,
🔍 处理基因: LGTVgp1
  ⚠️ 未找到蛋白质 RefSeq ID，跳过 LGTVgp1
```

输出的`gene_protein_map.txt`文件类似于：

```
115308162 XP_040067104.2,
115308164 XP_029842035.2,XP_029842036.2,
115308165 XP_029826729.1,
115308168 XP_029833004.2,XP_029833006.2,XP_029833007.2,XP_042145751.1,
115308169 XP_040360673.2,
```








