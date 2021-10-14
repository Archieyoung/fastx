# Fastx - An ultrafast toolkit for manipulating Fasta/q file.

## 作者

yangqi(Email: yangqi735@berrygenomics.com)

## 依赖

- linux
- gcc 4.8 or later
- cmake 3.2 or later
- boost 1.74 or later
- external libraries
    - zlib
    - libbz2
    - liblzma
    - libcurl

## 安装

```sh
cd fastx
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/boost ..
make
```
build/bin目录下将生成可执行文件`fastx`。如果出现链接错误，请检查是否正确安装了zlib, libbz2, liblzma, libcurl。

## 使用说明

`fastx sample` fasta/fastq序列抽样。可以抽取给定碱基数目的reads或者给定比例的reads。

示例：

```sh
# 按照碱基数进行抽样
fastx sample -i t.R1.fastq.gz \
    -I t.R2.fastq.gz \
    -o t.s.R1.fastq.gz \
    -O t.s.R2.fastq.gz \
    -b 3G \        # 抽取3Gbp（3000000000bp）的碱基，支持的单位有K(k)/M(m)/G(g), 如果不提供单位，例如-b 3000就是抽取3000bp的碱基
    -p /path/to/pigz

# 按照比例进行抽样
fastx sample -i t.R1.fastq.gz \
    -I t.R2.fastq.gz \
    -o t.s.R1.fastq.gz \
    -O t.s.R2.fastq.gz \
    -f 0.2 \        # 抽取 0.2 * 总碱基数 的碱基
    -p /path/to/pigz
```

`-b`和`-f`参数不能同时使用。

## 问题和不足

- 该项目目前还在开发中，后续将加入更多的功能。
