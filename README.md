# DMRArrayMetric

This is an R package to evaluate those DMRs(differential methylation regions) obtained using 450K methylation array data, contains four main functions including preinput(), compute_diff(), compute_metric() and rank_graph().  

First, the user needs to convert the original DMR into the formatted DMR, the user can use the function preinput().  

Second, the user use the formatted DMR to compute the diff of DMR by function compute_diff(). the user should provide a beta matrix, the name of case group and control group.  

Third, the user uses the DMR obtained from the second step to compute metric DMRn.  

Finally, the package provides a function rank_graph to draw a line chart for DMRn.  

The following Running the test sections will guide the user how to use the R package.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

This R package is based on R-4.1.3

```
You need install R-4.1.3
```

### Installing
Please enter R CMD, user can input the follwing command to install our R package.

```
install.packages("devtools")
install.packages("ggplot2")
library("devtools")
install_github("xqpeng/DMRArrayMetric")
```

## Running the tests


### preinput

This function generatea DMR data that includes chromosome, starting and end position, DMR length, and probe distribution.

```
library("DMRArrayMetric")
DMR_file = system.file("extdata","raw_DMR.csv",package = 'DMRArrayMetric')
DMR <- read.csv(DMR_file)
DMR$seqnames = as.character(DMR$seqnames)
convert_DMR <- preinput(DMR, chr=3, start=4, end=5, minProbeLength=3)
```
And repeat

```
[2022-07-19 16:30:40] Generate DMRs data that can be used for metric calculations #
Obtain the probes included in the DMR...
remove the number of probe in DMRs < 3 ...
Output converted DMR...
[2022-07-19 16:30:57] Done

```
### compute_diff

This function calculate the differential methylation value of the DMR and the number of probes and CpGs included and the length of the DMR.

```
beta_file = system.file("extdata","example_beta.RDS",package = 'DMRArrayMetric')
example_beta = readRDS(beta_file)
case = colnames(example_beta)[1:5]
control = colnames(example_beta)[6:10]
res_DMR = compute_diff(DMR = convert_DMR,beta = example_beta, human_ref = 'Cancer', Case_name=case, Control_name=control)

```
And repeat

```
[2022-07-19 16:32:16] # Calculate the number of cpg and the difference methylation #
Calculate the difference methylation value for each probe...
Calculate the difference methylation value for each DMR...
[2022-07-19 16:33:56] Done

```
### compute_metric

This function calculate the metric DMRn.
```
DMRn = compute_metric(DMR = res_DMR)

```
And repeat

```
[2022-07-19 16:34:48] # Calculate the metric DMRn #
The group differences of DMR range from 0 to 1 ,so the threshold will be 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 ...
Calculate the total number of cpg, probe and DMR length for each interval...
Calculate the DMRn at the threshold from 0 to 0.9 
[2022-07-19 16:34:48] Done

```

### rank_graph
This function will draw a line chart for DMRn. 
```
DMRn_file = read.table(system.file('extdata','rank_DMRn.txt',package='DMRArrayMetric'))
for(i in 1:7){
   DMRn_file[i,1] = eval(parse(text=DMRn_file[i,1]))

 }
rank_graph(DMRn_file,threshold=seq(0,0.9,0.1),Title='GroupA-GroupB')

```
And repeat

```
[2022-07-19 17:25:08] # Draw a line chart for DMRn. #
There are 7 methods in the file...
The plot is stored in the path: /home as pdf file...
[2022-07-19 17:25:08] Done

```
## Built With

* [R](https://www.r-project.org/) - R is a free software environment for statistical computing and graphics.


## Authors
* **xq Peng** - *Thesis basic framework and conception work* - [Central South University](https://life.csu.edu.cn/jsxx.jsp?urltype=news.NewsContentUrl&wbtreeid=1815&wbnewsid=3625)
* **wx Cui** - *Thesis code writing work* [Central South University](https://cse.csu.edu.cn/)





