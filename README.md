# OrthoBridge
基于序列以及blast的同源基因转换


## 基于蛋白质序列进行同源基因转换

对于亲缘关系非常远的物种，基于DNA序列进行blast不再合适，需要根据蛋白质序列进行blast，以对非模式动物的基因进行功能注释与解析。

主要步骤包括

- [从NCBI下载nr.gz数据库](##1.从NCBI下载nr.gz数据库)
- [从nr.gz中提取特定物种蛋白质序列](#2.从nr.gz中提取特定物种蛋白质序列)
- [删除冗余header中的非人信息](#3.删除冗余header中的非人信息)
- [一步构建参考数据库](#4.一步构建参考数据库)
- [基于基因列表批量从NCBI上下载对应蛋白质序列](#5.基于基因列表批量从NCBI上下载对应蛋白质序列)
- [blastp比对获取结果](#6.blastp比对获取结果)

### 1.从NCBI下载nr.gz数据库

首先使用`ascp`高速下载nr.gz数据库，优点是：速度较快；允许断点续传。

需要注意的是`-i`参数提供**链接所需的密令文件**，这个文件需要根据自己**aspera**安装的位置更换路径，其他参数不需要改变。

目前NCBI不再提供单独物种的蛋白质序列，需要将所有的序列下载后单独进行构建。

此外，下载过程中**NCBI**似乎对访问流量进行了限制，实际下载过程中，大约每下载**40G**就会中断一次，需要重新下载，所以断点续传功能非常重要。

```
ascp -v -k 1 -QT -l 300M -i /home/bio/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/FASTA/nr.gz .
```

这是NCBI对应的[nr.gz](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)所在的ftp链接位置。

### 2.从nr.gz中提取特定物种蛋白质序列

目的是将非模式物种的蛋白和人的蛋白序列进行blast，所以要构建人的reference。

首先提取人的序列，**两步法**构建：第一步提取人的蛋白的header，

```
zgrep "^>.*Homo sapiens" nr.gz > human_headers.txt
```

第二步，根据headers，使用`seqkit`根据全名提取，注意的是seqkit识别行名的时候不能有`>`，所以要先去掉headers文件里面的`>`

```
sed 's/^>//' human_headers.txt > human_headers2.txt
seqkit grep -n -f human_headers2.txt nr.fa -o human_proteins.fa
```

**一步法**构建：

```
seqkit grep -nrp "\[Homo sapiens\]" nr.fa > human_proteins.fasta
```

```
seqkit grep -nirp "Homo sapiens" nr.fa > human_proteins_2.fasta
```

出来的结果不完全一样,目前还没有发现到底是什么问题，只知道**Homo sapiens**是否加\[\]的结果不一样。但是所差条目很小，决定用第二种**一步法**构建的结果

**nr（非冗余蛋白数据库）里面会出现header里面包含多个id，但是均代表一条序列的情况，这是因为不同的数据库条目（如 RefSeq、GenBank、PDB）可能会存储相同的蛋白质序列，但具有不同的来源或注释。**

### 3.删除冗余header中的非人信息

由于我们的目的是构建人的reference，有的蛋白质序列包含了多个物种的信息，我们希望去掉这些信息，仅仅保留一条人的注释信息，因为后续比对也仅仅可以生成一个**acc id**的比对结果。

具体的思路就是提取header，根据特定分隔符分开不同的物种注释信息，匹配**Homo sapiens**，并将匹配到的第一个输出为新的header。

**值得一提的是，nr.gz文件中header其实是有一个分隔符分开不同的注释信息的，这个分隔符就是\x01，可视化页面可能会看到^A。** 

所以我们根据这个来进行分隔，基于`awk`命令写了一个脚本：

```
#!/bin/bash

# 检查是否提供了输入文件
if [ $# -eq 0 ]; then
    echo "Usage: $0 <fasta_file>"
    exit 1
fi

# 获取输入文件名
input_file=$1

# 处理FASTA文件
awk '
/^>/ {
    # 移除 ">" 并按 "]" 分割字符串
    header = substr($0, 2)  # 去掉开头的 ">"
    split(header, parts, "\x01")

    # 遍历切分后的部分，查找包含 "Homo sapiens" 的元素
    for (i in parts) {
        if (parts[i] ~ /Homo sapiens/) {
            # 如果找到 "Homo sapiens"，直接输出该部分
            print ">" parts[i]  # 输出包含 "Homo sapiens" 的 header
            break
        }
    }
}
!/^>/ { print }  # 输出原始的序列部分
' "$input_file" > "${input_file%.fasta}_cleaned.fasta"

echo "Processed file saved as ${input_file%.fasta}_cleaned.fasta"
```

运行`sh process_fasta.sh human_proteins_2.fasta`可以生成一个`human_proteins_2_cleaned.fasta`文件。

### 4.一步构建参考数据库

新建好一个`human_protein`文件夹后，运行以下命令，可以生成对应的参考数据库。

```
makeblastdb -in human_proteins_2_cleaned.fasta -dbtype prot -title human_proteins -out human_proteins -parse_seqids
```
生成的数据库，be like”

```
$ ll
总用量 599596
drwxrwxr-x 2 bio bio     20480 3月  19 17:34 ./
drwxrwxr-x 4 bio bio      4096 3月  19 09:48 ../
-rw-rw-r-- 1 bio bio  90947584 3月  19 17:34 human_proteins.pdb
-rw-rw-r-- 1 bio bio 217405424 3月  19 17:34 human_proteins.phr
-rw-rw-r-- 1 bio bio  15047216 3月  19 17:34 human_proteins.pin
-rw-rw-r-- 1 bio bio       575 3月  19 17:34 human_proteins.pjs
-rw-rw-r-- 1 bio bio   7523588 3月  19 17:34 human_proteins.pog
-rw-rw-r-- 1 bio bio  37176871 3月  19 17:34 human_proteins.pos
-rw-rw-r-- 1 bio bio  22570676 3月  19 17:34 human_proteins.pot
-rw-rw-r-- 1 bio bio 215727409 3月  19 17:34 human_proteins.psq
-rw-rw-r-- 1 bio bio     16384 3月  19 17:34 human_proteins.ptf
-rw-rw-r-- 1 bio bio   7523560 3月  19 17:34 human_proteins.pto
```

### 5.基于基因列表批量从NCBI上下载对应蛋白质序列

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

### 6.blastp比对获取结果

```
blastp -query Ixodes_proteins.fasta -db /home/Toshiba4/blast_db/human_protein/human_proteins -out panjun_Ixodes_to_human_transfer.txt -evalue 1e-5 -outfmt 6 -num_threads 64
```
运行以上代码，直接获得表格格式的比对结果。结果类似于：

```
XP_040067104.2	NP_055569.1	57.019	463	179	5	9	451	7	469	0.0	519
XP_040067104.2	BAF85725.1	56.803	463	180	5	9	451	7	469	0.0	519
XP_040067104.2	NP_001277154.1	54.969	322	131	3	144	451	2	323	3.93e-125	369
XP_040067104.2	AAH04502.2	54.517	321	132	3	145	451	7	327	1.34e-123	365
XP_040067104.2	KAI2550776.1	54.508	244	106	2	201	439	2	245	2.65e-90	279
XP_040067104.2	NP_110410.1	34.054	370	232	5	7	367	33	399	5.53e-70	233
XP_040067104.2	NP_001316477.1	35.191	341	209	5	36	367	2	339	3.00e-67	224
XP_040067104.2	XP_047283599.1	35.191	341	209	5	36	367	34	371	5.92e-67	224
XP_040067104.2	NP_001316473.1	32.642	386	232	6	7	367	33	415	1.14e-65	222
XP_040067104.2	XP_024304467.1	33.613	357	209	6	36	367	2	355	3.56e-63	214
XP_040067104.2	NP_001316474.1	39.267	191	114	2	178	367	1	190	3.85e-43	156
XP_040067104.2	EAW91753.1	64.706	85	25	1	9	88	7	91	2.71e-27	108
XP_029842035.2	CAG5082931.1	32.530	166	100	3	41	199	30	190	1.45e-10	65.9
XP_029842036.2	CAG5082931.1	32.530	166	100	3	29	187	30	190	1.01e-10	65.9
XP_029826729.1	NP_055302.1	69.565	253	75	1	1	253	1	251	9.18e-134	385
XP_029826729.1	BAD97329.1	69.170	253	76	1	1	253	1	251	6.28e-133	383
XP_029826729.1	KAI2535999.1	72.489	229	63	0	1	229	1	229	5.43e-128	367
XP_029826729.1	NP_001278931.1	71.585	183	50	1	71	253	11	191	5.93e-95	285
XP_029826729.1	BAG57430.1	71.038	183	51	1	71	253	11	191	5.42e-94	282
XP_029826729.1	KAI2536000.1	70.588	153	45	0	1	153	1	153	2.01e-81	246
XP_029826729.1	KAI2536001.1	69.286	140	43	0	1	140	1	140	6.07e-72	221
XP_029826729.1	KAI2536002.1	60.000	50	20	0	1	50	1	50	1.11e-11	63.2
XP_029833004.2	2Y3I_A	72.330	412	113	1	4	414	5	416	0.0	608
XP_029833004.2	NP_000282.1	72.330	412	113	1	4	414	5	416	0.0	608
XP_029833004.2	3C39_A	72.330	412	113	1	4	414	8	419	0.0	608
XP_029833004.2	2WZB_A	72.330	412	113	1	4	414	4	415	0.0	608
```

每个序列会输出top match的一些结果。

1. 使用`awk '!seen[$1]++' panjun_Ixodes_to_human_transfer.txt > best_hits.txt` 可以获得top1的match结果

2. 使用`awk '$3 >= 90' results.txt > high_identity_hits.txt` 可以根据序列相似度过滤特定结果。




