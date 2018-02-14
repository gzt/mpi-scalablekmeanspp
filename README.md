# mpi-scalablekmeanspp
This is a port of the scalable *k*-means++ (*k*-means||) to the OpenMPI framework

See the paper by Bahmani et al: https://arxiv.org/abs/1203.6402

I put it into C OpenMPI because I didn't find any implementation in C/OpenMPI of it elsewhere. I have included a wrapper with some other guy's code to allow you to run it in the context of a problem, but I recommend disentangling from it to use it yourself. The header files I use for their matrix and vector macros may conflict with what you have for what you are doing, however you may be able to disentangle from that fairly easily. The weighted *k*-means part is based on the Hartigan and Wong algorithm and runs on a single processor, you may find that useful on its own.
