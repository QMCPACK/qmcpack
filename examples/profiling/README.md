# Profiling QMCPACK

## LLVM-XRAY based
With CUDA:
### LLVM 4.0.1 and NVCC 9.1 -- allows simultaneous nvprof and xray instrumentation of CPU code
* Build: see [build_script](./build_llvm_nvcc_xray_instrumented.sh)
* Running: With the performance tests set up.
``` shell
[epd@oxygen ]$ export QMCPACK_APP=/path/to/your/instrumented/bin/qmcpack
[epd@oxygen ]$ cd /your_scratch/performace/NiO/sample/dmc-a4-e48-gpu
[epd@oxygen ]$ OMP_NUM_THREADS=12 XRAY_OPTIONS="patch_premain=true xray_mode=xray-basic" nvprof -o simul.nvprof $QMCPACK_APP NiO-fcc-S1-dmc.xml
[epd@oxygen ]$ llvm-xray account xray-log.qmcpack.yrxHG1 -instr_map=$QMCPACK_APP -sort=med -top=10 -sortorder=dsc

Functions with latencies: 514
   funcid      count [      min,       med,       90p,       99p,       max]       sum  function
        1          1 [12.220955, 12.220955, 12.220955, 12.220955, 12.220955] 12.220955  <invalid>:0:0: main
       55          1 [ 5.188287,  5.188287,  5.188287,  5.188287,  5.188287]  5.188287  <invalid>:0:0: qmcplusplus::QMCMain::execute()
       59          1 [ 2.152482,  2.152482,  2.152482,  2.152482,  2.152482]  2.152482  <invalid>:0:0: qmcplusplus::QMCMain::validateXML()
       64          1 [ 1.795763,  1.795763,  1.795763,  1.795763,  1.795763]  1.795763  <invalid>:0:0: qmcplusplus::WaveFunctionPool::put(_xmlNode*)
     1266          1 [ 1.795672,  1.795672,  1.795672,  1.795672,  1.795672]  1.795672  <invalid>:0:0: qmcplusplus::WaveFunctionFactory::build(_xmlNode*, bool)
     1267          1 [ 1.488744,  1.488744,  1.488744,  1.488744,  1.488744]  1.488744  <invalid>:0:0: qmcplusplus::WaveFunctionFactory::addFermionTerm(_xmlNode*)
     1767          1 [ 1.488722,  1.488722,  1.488722,  1.488722,  1.488722]  1.488722  <invalid>:0:0: qmcplusplus::SlaterDetBuilder::put(_xmlNode*)
      335          1 [ 1.198301,  1.198301,  1.198301,  1.198301,  1.198301]  1.198301  <invalid>:0:0: qmcplusplus::VMCcuda::runWithDrift()
       56          3 [ 0.793548,  1.041833,  1.198912,  1.198912,  1.198912]  3.034292  <invalid>:0:0: qmcplusplus::QMCMain::executeQMCSection(_xmlNode*, bool)
       58          3 [ 0.793514,  1.041820,  1.198901,  1.198901,  1.198901]  3.034235  <invalid>:0:0: qmcplusplus::QMCMain::runQMC(_xmlNode*)
```

If you also have a new version of llvm installed you can use the newer xray tools with the old trace. The best of these allows you to convert to a event-trace format that can then be massaged further to view in a graphical tool.

``` shell
/home/epd/opt/llvm-7/bin/llvm-xray convert -output-format=trace_event -instr_map=$QMCPACK_APP -symbolize -sort -output=dmc-a4-e48-gpu.trace xray-log.qmcpack.yrxHG1
```

Then convert into chrome viewable html doc using the [catapult](https://github.com/catapult-project/catapult) tool.

``` shell
~/codes/catapult/tracing/bin/trace2html dmc-a4-e48-gpu.trace -o dmc-a4-e48-gpu.html
```

[example -- only works with chrome](http://cdash-minimal.ornl.gov/profiling/dmc-a4-e48-gpu.html)

You can then download the nvprof trace:
``` shell
$ scp /your_scratch/performace/NiO/sample/dmc-a4-e48-gpu/simul.nvprof ./
$ nvvp
```
You'll need to import simul.nvprof to a new session in nvvp.
Then you'll be able look at the traces from the CPU and GPU point of view at the same time.

