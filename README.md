** if you need to transfer data from lab laptop,
**use teamviewer app  316 720 782 password 3ygs75
**find directory 
**open file transfer from pull down on top bar

but you shouldnt have to,

** data is on nupacs1:/data2/bacon/sipmRuns/

root library is in directory bobj. type make in directory bobj to make library

raw data files *.txt should be all in one directory under top directory nupacs1:/data2/bacon/sipmRuns/<dir>

  run as "root readRaw.cc"

makes root output file  /data2/bacon/sipmRuns/rootData/<dir>.root

then on output "root readRun.cc"

then on output "root post.cc"

